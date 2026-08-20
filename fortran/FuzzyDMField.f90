    ! Archived experimental Klein-Gordon + EFA path for ultralight axion DM.
    !
    ! Solves the KG equation for V(phi) = (1/2) m^2 phi^2 to get correct background
    ! evolution (frozen w=-1 at early times, matter-like w=0 after oscillation onset).
    ! After a_match (where m/H = N_match), switches to effective fluid approximation
    ! for both background and perturbations.
    !
    ! Active axion configurations are disabled pending abundance normalization,
    ! transition-matching fixes, and an independent reference comparison.
    !
    ! References:
    !   Hu, Barkana, Gruzinov 2000; Hlozek+ 2015, 2018
    !   Passaglia & Hu 2022 (improved EFA); Marsh 2016

    module FuzzyDMField
    use DarkEnergyInterface
    use results
    use constants
    use classes
    use Interpolation
    use config, only: GlobalError, error_unsupported_params
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    real(dl), parameter :: Tpl_FDM = sqrt(kappa*hbar/c**5)
    real(dl), parameter :: Vunits = MPC_in_sec**2 / Tpl_FDM**2

    type, extends(TDarkEnergyModel) :: TFuzzyDMField
        ! User parameters (must come first for Python ctypes mapping)
        real(dl) :: m_axion = 1e-22_dl       ! axion mass [eV]
        real(dl) :: f_axion = 0._dl          ! fraction of total CDM+axion budget
        real(dl) :: omega_axion_h2 = 0._dl   ! alternative: direct Omega_ax h^2
        integer  :: N_match = 100            ! KG->EFA transition criterion m/H
        integer  :: npoints_bg = 5000        ! baseline background integration points
        integer  :: min_steps_per_osc_bg = 10
        ! AxiCLASS improvements
        integer  :: potential_type = 1       ! 1=quadratic, 2=cosine [1-cos(phi/f)]^n
        integer  :: n_potential = 1          ! power of cosine potential
        real(dl) :: f_decay = 1._dl          ! decay constant [M_Pl] for cosine potential
        logical  :: use_improved_efa = .true. ! Passaglia-Hu improved EFA

        ! Cached quantities (set in Init, not mapped in Python)
        real(dl) :: m_planck = 0._dl
        real(dl) :: m_conf = 0._dl
        real(dl) :: frac_lambda0 = 1._dl
        real(dl) :: a_osc = 0._dl
        real(dl) :: a_match = 1._dl
        real(dl) :: grho_match = 0._dl
        real(dl) :: grho_ax_today = 0._dl
        real(dl) :: cached_grhov = 0._dl
        real(dl) :: log_astart = 0._dl

        ! Background density spline arrays
        integer :: bg_n = 0
        real(dl), dimension(:), allocatable :: bg_a, bg_rho, bg_ddrho
        ! w(a) spline for improved EFA
        real(dl), dimension(:), allocatable :: bg_w, bg_ddw

        ! State pointer for background integration
        class(CAMBdata), pointer :: bg_State => null()
    contains
    procedure :: Init => TFuzzyDMField_Init
    procedure :: Validate => TFuzzyDMField_Validate
    procedure :: BackgroundDensityAndPressure => TFuzzyDMField_BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy => TFuzzyDMField_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TFuzzyDMField_PerturbationEvolve
    procedure :: PerturbationInitial => TFuzzyDMField_PerturbationInitial
    procedure :: ReadParams => TFuzzyDMField_ReadParams
    procedure, nopass :: PythonClass => TFuzzyDMField_PythonClass
    procedure, nopass :: SelfPointer => TFuzzyDMField_SelfPointer
    procedure :: PrintFeedback => TFuzzyDMField_PrintFeedback
    procedure :: VofPhi_val
    procedure :: VofPhi_deriv1
    procedure :: grho_ax_t_at
    end type TFuzzyDMField

    ! dverk ODE solver (from subroutines.f90)
    procedure(TClassDverk) :: dverk

    public TFuzzyDMField
    contains

    subroutine TFuzzyDMField_Validate(this, OK)
    class(TFuzzyDMField), intent(in) :: this
    logical, intent(inout) :: OK

    call this%TDarkEnergyModel%Validate(OK)
    if (.not. ieee_is_finite(this%m_axion) .or. .not. ieee_is_finite(this%f_axion) .or. &
        .not. ieee_is_finite(this%omega_axion_h2) .or. .not. ieee_is_finite(this%f_decay) .or. &
        this%m_axion <= 0._dl .or. this%f_axion < 0._dl .or. this%f_axion >= 1._dl .or. &
        this%omega_axion_h2 < 0._dl .or. this%N_match <= 0 .or. &
        (this%potential_type /= 1 .and. this%potential_type /= 2) .or. &
        this%n_potential < 1 .or. this%f_decay <= 0._dl) then
        OK = .false.
        write(*,*) 'FuzzyDMField parameters are outside their finite physical ranges.'
    end if
    if (this%f_axion > 0._dl .or. this%omega_axion_h2 > 0._dl) then
        OK = .false.
        write(*,*) 'Active FuzzyDMField is disabled pending abundance and KG-to-EFA validation.'
    end if

    end subroutine TFuzzyDMField_Validate


    function VofPhi_val(this, phi) result(V)
    ! V(phi) for different potential types
    ! type 1: (Vunits/2)*m_pl^2*phi^2  (quadratic)
    ! type 2: Lambda^4 * [1-cos(phi/f)]^n  (cosine, AxiCLASS)
    class(TFuzzyDMField), intent(in) :: this
    real(dl), intent(in) :: phi
    real(dl) :: V
    real(dl) :: Lambda4, cos_term

    select case(this%potential_type)
    case(2)
        ! Cosine potential: Lambda^4 = m^2 * f^2 / n (matches quadratic at small phi)
        Lambda4 = Vunits * this%m_planck**2 * this%f_decay**2 / this%n_potential
        cos_term = 1._dl - cos(phi / this%f_decay)
        V = Lambda4 * cos_term**this%n_potential + this%frac_lambda0 * this%cached_grhov
    case default
        ! Quadratic: V = (1/2) m^2 phi^2
        V = Vunits * 0.5_dl * this%m_planck**2 * phi**2 &
            + this%frac_lambda0 * this%cached_grhov
    end select
    end function VofPhi_val

    function VofPhi_deriv1(this, phi) result(Vp)
    ! dV/dphi for different potential types
    class(TFuzzyDMField), intent(in) :: this
    real(dl), intent(in) :: phi
    real(dl) :: Vp
    real(dl) :: Lambda4, cos_term

    select case(this%potential_type)
    case(2)
        Lambda4 = Vunits * this%m_planck**2 * this%f_decay**2 / this%n_potential
        cos_term = 1._dl - cos(phi / this%f_decay)
        ! dV/dphi = Lambda^4 * n * [1-cos]^(n-1) * sin(phi/f) / f
        if (this%n_potential == 1) then
            Vp = Lambda4 * sin(phi / this%f_decay) / this%f_decay
        else
            Vp = Lambda4 * this%n_potential * cos_term**(this%n_potential - 1) &
                * sin(phi / this%f_decay) / this%f_decay
        end if
    case default
        Vp = Vunits * this%m_planck**2 * phi
    end select
    end function VofPhi_deriv1


    ! KG background evolution: d/da of [phi, a^2*phidot]
    subroutine FDM_EvolveBackground(self, num, a, y, yprime)
    class(TFuzzyDMField) :: self
    integer num
    real(dl) a, y(num), yprime(num)
    real(dl) a2, phi, phidot, grhode, tot, adot

    a2 = a**2
    phi = y(1)
    phidot = y(2) / a2
    grhode = a2 * (0.5_dl * phidot**2 + a2 * self%VofPhi_val(phi))
    tot = self%bg_State%grho_no_de(a) + grhode
    adot = sqrt(tot / 3._dl)
    yprime(1) = phidot / adot         ! dphi/da
    yprime(2) = -a2**2 * self%VofPhi_deriv1(phi) / adot  ! d(a^2*phidot)/da

    end subroutine FDM_EvolveBackground

    ! Same but in log(a) variable
    subroutine FDM_EvolveBackgroundLog(self, num, loga, y, yprime)
    class(TFuzzyDMField) :: self
    integer num
    real(dl) loga, y(num), yprime(num)
    real(dl) a

    a = exp(loga)
    call FDM_EvolveBackground(self, num, a, y, yprime)
    yprime = yprime * a  ! chain rule: dy/dloga = a * dy/da

    end subroutine FDM_EvolveBackgroundLog


    function grho_ax_t_at(this, a) result(grho_ax_t)
    ! Interpolate axion grhov_t = 8*pi*G*rho_ax*a^2 at scale factor a
    class(TFuzzyDMField), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: grho_ax_t
    real(dl) :: da, a0, b0, ho2o6
    integer :: ix, i

    if (a < exp(this%log_astart) .or. this%bg_n == 0) then
        grho_ax_t = 0._dl
        return
    end if

    if (a >= this%a_match) then
        ! Matter-like: grho_ax_t proportional to 1/a
        grho_ax_t = this%grho_match * (this%a_match / a)
        return
    end if

    ! Cubic spline interpolation in bg_a table
    ! Find interval by bisection
    ix = 1
    do i = 1, this%bg_n - 1
        if (this%bg_a(i+1) > a) then
            ix = i
            exit
        end if
        ix = i
    end do
    ix = min(ix, this%bg_n - 1)

    da = this%bg_a(ix+1) - this%bg_a(ix)
    a0 = (this%bg_a(ix+1) - a) / da
    b0 = 1._dl - a0
    ho2o6 = da**2 / 6._dl
    grho_ax_t = b0 * this%bg_rho(ix+1) + a0 * (this%bg_rho(ix) &
        - b0 * ((a0 + 1) * this%bg_ddrho(ix) &
        + (2 - a0) * this%bg_ddrho(ix+1)) * ho2o6)

    end function grho_ax_t_at


    subroutine TFuzzyDMField_Init(this, State)
    class(TFuzzyDMField), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) :: h2, grho_rad, adotoa_test, grhoa2
    real(dl) :: alo, ahi, amid
    real(dl) :: initial_phi
    real(dl) :: aend, afrom, astart
    integer, parameter :: NumEqs = 2
    real(dl) :: cc(24), w(NumEqs,9), y(NumEqs)
    integer :: ind, i, ix, npoints, tot_points
    real(dl), parameter :: splZero = 0._dl
    real(dl) :: lastsign, da_osc, last_a, a2, phi, phidot
    real(dl) :: dloga, max_a_log, da_lin, integrate_tol
    real(dl), dimension(:), allocatable :: tmp_a, tmp_rho
    integer :: npoints_log, npoints_linear

    this%is_cosmological_constant = .false.
    this%num_perturb_equations = 2
    astart = 1e-7_dl
    integrate_tol = 1e-6_dl

    if (this%f_axion > 0._dl .or. this%omega_axion_h2 > 0._dl) then
        call GlobalError('Active FuzzyDMField is disabled pending abundance and KG-to-EFA validation.', &
            error_unsupported_params)
        return
    end if

    select type(S => State)
    class is (CAMBdata)

    this%bg_State => S
    this%cached_grhov = S%grhov
    this%log_astart = log(astart)

    ! Unit conversions
    this%m_planck = this%m_axion / 2.435e27_dl
    this%m_conf = this%m_axion / 6.396e-30_dl

    ! Target axion density at a=1
    h2 = (S%CP%H0 / 100._dl)**2
    if (this%omega_axion_h2 > 0._dl) then
        this%grho_ax_today = S%grhocrit * this%omega_axion_h2 / h2
    else if (this%f_axion > 0._dl .and. this%f_axion < 1._dl) then
        this%grho_ax_today = S%grhoc * this%f_axion / (1._dl - this%f_axion)
    else if (this%f_axion >= 1._dl) then
        this%grho_ax_today = S%grhov * 0.9_dl
    else
        this%is_cosmological_constant = .true.
        this%num_perturb_equations = 0
        return
    end if

    ! Lambda fraction of grhov
    if (S%grhov > 0._dl) then
        this%frac_lambda0 = 1._dl - this%grho_ax_today / S%grhov
    else
        this%frac_lambda0 = 0._dl
    end if

    ! a_osc estimate: m_conf = 3*adotoa in radiation domination
    grho_rad = S%grhog + S%grhornomass
    this%a_osc = sqrt(3._dl * sqrt(grho_rad / 3._dl) / this%m_conf)

    ! If field never oscillates (ultra-light), treat as cosmological constant
    if (this%a_osc >= 0.999_dl) then
        this%is_cosmological_constant = .true.
        this%num_perturb_equations = 0
        return
    end if
    this%a_osc = min(this%a_osc, 0.99_dl)

    ! Find a_match where m_conf/adotoa = N_match, bisection
    alo = this%a_osc
    ahi = 1.0_dl
    do i = 1, 100
        amid = 0.5_dl * (alo + ahi)
        grhoa2 = S%grho_no_de(amid) + S%grhov * amid**4
        adotoa_test = sqrt(grhoa2 / 3._dl) / amid
        if (this%m_conf / adotoa_test > this%N_match) then
            ahi = amid
        else
            alo = amid
        end if
    end do
    this%a_match = 0.5_dl * (alo + ahi)
    this%a_match = min(this%a_match, 0.999_dl)

    ! Ensure astart is well before a_osc (adapt for heavy masses)
    astart = min(astart, this%a_osc * 0.01_dl)
    astart = max(astart, 1e-10_dl)
    this%log_astart = log(astart)

    ! Initial field value from frozen density estimate
    ! V_ax(phi0) * a_osc^3 ~ grho_ax_today
    if (this%potential_type == 2) then
        ! For cosine: V = Lambda^4 * (1-cos(phi/f))^n, use small-angle to estimate
        initial_phi = sqrt(2._dl * this%grho_ax_today / &
            (this%a_osc**3 * Vunits * this%m_planck**2))
        ! Clamp to < pi*f_decay to avoid crossing the potential peak
        initial_phi = min(initial_phi, 3._dl * this%f_decay)
    else
        initial_phi = sqrt(2._dl * this%grho_ax_today / &
            (this%a_osc**3 * Vunits * this%m_planck**2))
    end if

    ! === KG background integration ===
    dloga = (-this%log_astart) / (this%npoints_bg - 1)
    npoints = this%npoints_bg
    allocate(tmp_a(npoints), tmp_rho(npoints))

    y(1) = initial_phi
    y(2) = 0._dl  ! a^2 * phidot, initially zero (frozen)

    a2 = astart**2
    phi = y(1)
    phidot = 0._dl
    tmp_a(1) = astart
    ! Use VofPhi_val - but it adds frac_lambda0*grhov, so subtract it for axion-only part
    tmp_rho(1) = a2 * (this%VofPhi_val(phi) - this%frac_lambda0 * this%cached_grhov)

    da_osc = 1._dl
    last_a = astart
    lastsign = 0._dl

    ind = 1
    afrom = this%log_astart
    npoints_log = 1
    do i = 1, npoints - 1
        aend = this%log_astart + dloga * i
        if (exp(aend) > this%a_match) aend = log(this%a_match)
        ! Skip if target equals current position (floating point)
        if (abs(aend - afrom) < 1e-14_dl) then
            npoints_log = ix
            exit
        end if
        ix = i + 1
        tmp_a(ix) = exp(aend)
        a2 = tmp_a(ix)**2

        call dverk(this, NumEqs, FDM_EvolveBackgroundLog, afrom, y, aend, &
            integrate_tol, ind, cc, NumEqs, w)
        if (global_error_flag /= 0) then
            write(*,*) 'TFuzzyDMField: KG integration error at a=', tmp_a(ix)
            return
        end if
        call FDM_EvolveBackgroundLog(this, NumEqs, aend, y, w(:,1))

        phi = y(1)
        phidot = y(2) / a2
        ! Axion part of grhov_t (without Lambda)
        tmp_rho(ix) = 0.5_dl * phidot**2 &
            + a2 * (this%VofPhi_val(phi) - this%frac_lambda0 * this%cached_grhov)

        ! Track oscillation scale
        if (i == 1) then
            lastsign = y(2)
        elseif (y(2) * lastsign < 0) then
            da_osc = min(da_osc, tmp_a(ix) - last_a)
            last_a = tmp_a(ix)
            lastsign = y(2)
        end if

        npoints_log = ix
        if (tmp_a(ix) >= this%a_match * 0.9999_dl) exit
        if (tmp_a(ix) * (exp(dloga) - 1) * this%min_steps_per_osc_bg > da_osc) exit
    end do

    max_a_log = tmp_a(npoints_log)

    ! Linear spacing for remaining steps up to a_match
    if (max_a_log < this%a_match * 0.999_dl) then
        da_lin = min(max_a_log * (exp(dloga) - 1), &
            da_osc / this%min_steps_per_osc_bg, &
            (this%a_match - max_a_log) / max(npoints - npoints_log, 1))
        npoints_linear = int((this%a_match - max_a_log) / da_lin) + 1
        da_lin = (this%a_match - max_a_log) / npoints_linear
    else
        npoints_linear = 0
        da_lin = 0._dl
    end if

    tot_points = npoints_log + npoints_linear

    ! Allocate final arrays
    if (allocated(this%bg_a)) deallocate(this%bg_a, this%bg_rho, this%bg_ddrho)
    if (allocated(this%bg_w)) deallocate(this%bg_w, this%bg_ddw)
    allocate(this%bg_a(tot_points), this%bg_rho(tot_points), this%bg_ddrho(tot_points))
    allocate(this%bg_w(tot_points), this%bg_ddw(tot_points))
    this%bg_a(1:npoints_log) = tmp_a(1:npoints_log)
    this%bg_rho(1:npoints_log) = tmp_rho(1:npoints_log)

    ! Linear steps
    if (npoints_linear > 0) then
        ind = 1
        afrom = max_a_log
        do i = 1, npoints_linear
            ix = npoints_log + i
            aend = max_a_log + da_lin * i
            if (aend > this%a_match) aend = this%a_match
            if (abs(aend - afrom) < 1e-14_dl) exit
            a2 = aend**2
            this%bg_a(ix) = aend

            call dverk(this, NumEqs, FDM_EvolveBackground, afrom, y, aend, &
                integrate_tol, ind, cc, NumEqs, w)
            if (global_error_flag /= 0) then
                write(*,*) 'TFuzzyDMField: KG integration error (linear) at a=', aend
                return
            end if
            call FDM_EvolveBackground(this, NumEqs, aend, y, w(:,1))

            phi = y(1)
            phidot = y(2) / a2
            this%bg_rho(ix) = 0.5_dl * phidot**2 &
                + a2 * (this%VofPhi_val(phi) - this%frac_lambda0 * this%cached_grhov)
        end do
    end if

    this%bg_n = tot_points

    ! Build cubic spline for bg_rho(a)
    call spline(this%bg_a, this%bg_rho, tot_points, splZero, splZero, this%bg_ddrho)

    ! Compute w(a) from rho(a): p = -rho - a/3 * drho/da  =>  w = p/rho
    ! grho_ax_t = rho_ax * a^2.  d(grho_ax_t)/da = (2*rho + a*drho/da)*a = a^2*(2*rho/a + drho/da)
    ! For matter-like (w=0): grho_t ~ 1/a => d(grho_t)/da ~ -1/a^2
    ! w = -1 - a*d(ln grho_t)/da / 3 - 2/3  => w = -1 - (a/(3*grho_t)) * d(grho_t)/da - 2/3
    ! Actually: grho_t = 8piG rho a^2, so rho = grho_t/a^2
    ! rho' = (grho_t'/a^2 - 2*grho_t/a^3) = (grho_t' - 2*grho_t/a) / a^2
    ! p + rho = -a*rho'/3 => p = -rho - a*rho'/3
    ! w = p/rho = -1 - a*rho'/(3*rho) = -1 - a*(d ln rho/da)/3
    ! = -1 - a/(3*grho_t) * (d(grho_t)/da - 2*grho_t/a) / 1
    ! = -1 - (a*grho_t' - 2*grho_t)/(3*grho_t)
    ! = -1/3 - a*grho_t'/(3*grho_t)
    ! At a_match, w ~ 0 (matter-like) so grho_t ~ 1/a
    do i = 1, tot_points
        if (i == 1) then
            this%bg_w(i) = -1._dl  ! frozen field
        else if (i == tot_points) then
            this%bg_w(i) = 0._dl  ! matter-like at match
        else
            block
                real(dl) :: grho_plus, grho_minus, da_w, dgrho_da, a_w
                a_w = this%bg_a(i)
                da_w = this%bg_a(i+1) - this%bg_a(i-1)
                dgrho_da = (this%bg_rho(i+1) - this%bg_rho(i-1)) / da_w
                if (this%bg_rho(i) > 0._dl) then
                    this%bg_w(i) = -1._dl/3._dl - a_w * dgrho_da / (3._dl * this%bg_rho(i))
                else
                    this%bg_w(i) = 0._dl
                end if
                ! Clamp w to physical range
                this%bg_w(i) = max(min(this%bg_w(i), 1._dl), -1._dl)
            end block
        end if
    end do
    call spline(this%bg_a, this%bg_w, tot_points, splZero, splZero, this%bg_ddw)

    ! Record axion grhov_t at a_match (instantaneous, ~cycle average for large N_match)
    this%grho_match = this%bg_rho(tot_points)

    end select
    deallocate(tmp_a, tmp_rho)

    end subroutine TFuzzyDMField_Init


    subroutine TFuzzyDMField_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    class(TFuzzyDMField), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w
    real(dl) :: a2, grho_ax_t, grho_lam_t

    if (this%is_cosmological_constant) then
        grhov_t = grhov * a * a
        if (present(w)) w = -1._dl
        return
    end if

    a2 = a**2
    grho_lam_t = this%frac_lambda0 * grhov * a2
    grho_ax_t = this%grho_ax_t_at(a)
    grhov_t = grho_ax_t + grho_lam_t

    if (present(w)) then
        if (grhov_t > 0._dl) then
            ! w = (p_ax*a^2 + p_lam*a^2) / grhov_t
            ! p_ax ~ 0 (cycle-averaged), p_lam = -rho_lam
            w = -grho_lam_t / grhov_t
        else
            w = -1._dl
        end if
    end if

    end subroutine TFuzzyDMField_BackgroundDensityAndPressure


    subroutine TFuzzyDMField_PerturbationInitial(this, y, a, tau, k)
    class(TFuzzyDMField), intent(in) :: this
    real(dl), intent(out) :: y(:)
    real(dl), intent(in) :: a, tau, k

    ! Frozen field: no perturbations initially
    y = 0._dl

    end subroutine TFuzzyDMField_PerturbationInitial


    subroutine TFuzzyDMField_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TFuzzyDMField), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix
    real(dl) :: delta_ax, v_ax, grho_ax_t

    if (this%is_cosmological_constant .or. a < this%a_osc) then
        dgrhoe = 0._dl
        dgqe = 0._dl
        return
    end if

    delta_ax = ay(w_ix)
    v_ax = ay(w_ix + 1)
    grho_ax_t = this%grho_ax_t_at(a)

    dgrhoe = grho_ax_t * delta_ax
    dgqe = grho_ax_t * v_ax

    end subroutine TFuzzyDMField_PerturbedStressEnergy


    subroutine TFuzzyDMField_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TFuzzyDMField), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) :: delta_ax, v_ax, cs2, cs2_phi, Hv3_over_k
    real(dl) :: w_efa, ca2, H_over_m, one_plus_w

    if (this%is_cosmological_constant) return

    if (a < this%a_osc) then
        ayprime(w_ix) = 0._dl
        ayprime(w_ix + 1) = 0._dl
        return
    end if

    delta_ax = y(w_ix)
    v_ax = y(w_ix + 1)

    ! Quantum pressure sound speed (always present)
    cs2_phi = k**2 / (4._dl * this%m_conf**2 * a**2 + k**2)

    if (this%use_improved_efa) then
        ! Passaglia & Hu 2022 improved EFA corrections
        H_over_m = adotoa / (this%m_conf * a)

        ! Improved effective EOS: w_efa = (3/2)(H/m)^2
        w_efa = 1.5_dl * H_over_m**2

        ! Improved sound speed: cs2 = cs2_phi + (5/4)(H/m)^2
        cs2 = cs2_phi + 1.25_dl * H_over_m**2

        ! Adiabatic sound speed: ca2 = w - a*dw/da / (3*(1+w))
        ! For improved EFA, ca2 ~ w_efa (since w varies slowly)
        ca2 = w_efa
        one_plus_w = max(1._dl + w_efa, 1e-10_dl)

        ! GDM form perturbation equations
        Hv3_over_k = 3._dl * adotoa * v_ax / k
        ! delta' = -(1+w)(k*v + k*z) - 3H(cs2-w)(delta + 3H*v/k)
        ayprime(w_ix) = -one_plus_w * (k * v_ax + k * z) &
            - 3._dl * adotoa * (cs2 - w_efa) * (delta_ax + Hv3_over_k)
        ! v' = -(1-3w)H*v + cs2*k*delta/(1+w)
        ayprime(w_ix + 1) = -(1._dl - 3._dl * w_efa) * adotoa * v_ax &
            + cs2 * k * delta_ax / one_plus_w
    else
        ! Original EFA: w=0, cs2 = quantum pressure only
        cs2 = cs2_phi
        Hv3_over_k = 3._dl * adotoa * v_ax / k
        ayprime(w_ix) = -k * v_ax - k * z &
            - 3._dl * adotoa * cs2 * (delta_ax + Hv3_over_k)
        ayprime(w_ix + 1) = -adotoa * v_ax + cs2 * k * delta_ax
    end if

    end subroutine TFuzzyDMField_PerturbationEvolve


    subroutine TFuzzyDMField_ReadParams(this, Ini)
    use IniObjects
    class(TFuzzyDMField) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyModel%ReadParams(Ini)
    this%m_axion = Ini%Read_Double('FuzzyDMField_m_axion', 1e-22_dl)
    this%f_axion = Ini%Read_Double('FuzzyDMField_f_axion', 0._dl)
    this%omega_axion_h2 = Ini%Read_Double('FuzzyDMField_omega_axion_h2', 0._dl)
    this%N_match = Ini%Read_Int('FuzzyDMField_N_match', 100)
    this%potential_type = Ini%Read_Int('FuzzyDMField_potential_type', 1)
    this%n_potential = Ini%Read_Int('FuzzyDMField_n_potential', 1)
    this%f_decay = Ini%Read_Double('FuzzyDMField_f_decay', 1._dl)

    end subroutine TFuzzyDMField_ReadParams


    function TFuzzyDMField_PythonClass()
    character(LEN=:), allocatable :: TFuzzyDMField_PythonClass
    TFuzzyDMField_PythonClass = 'FuzzyDMField'
    end function TFuzzyDMField_PythonClass

    subroutine TFuzzyDMField_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TFuzzyDMField), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TFuzzyDMField_SelfPointer


    subroutine TFuzzyDMField_PrintFeedback(this, FeedbackLevel)
    class(TFuzzyDMField) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("FuzzyDMField: m_axion = ",ES12.4," eV")') this%m_axion
        write(*,'("  f_axion = ",F8.5)') this%f_axion
        write(*,'("  a_osc = ",ES10.3,"  a_match = ",ES10.3)') this%a_osc, this%a_match
        write(*,'("  frac_lambda0 = ",F10.6)') this%frac_lambda0
        write(*,'("  grho_ax_today = ",ES12.4)') this%grho_ax_today
        write(*,'("  grho_match = ",ES12.4)') this%grho_match
    end if

    end subroutine TFuzzyDMField_PrintFeedback

    end module FuzzyDMField
