    module FuzzyDM
    use DarkMatterInteraction
    use results
    use constants
    use classes
    use Interpolation
    use config, only: GlobalError, error_unsupported_params, global_error_flag
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    integer, parameter :: default_match_ratio = 50
    integer, parameter :: background_points = 6000
    integer, parameter :: fluid_points = 3000
    real(dl), parameter :: integrate_tol = 2e-8_dl

    type, extends(TDarkMatterModel) :: TFuzzyDM
        ! Public parameters: keep first for the Python ctypes mapping.
        real(dl) :: m_axion = 0._dl
        real(dl) :: omega_axion_h2 = 0._dl
        real(dl) :: f_axion = 0._dl
        integer :: match_ratio = default_match_ratio
        real(dl) :: m_conf = 0._dl
        real(dl) :: a_osc = 0._dl
        real(dl) :: a_match = 1._dl
        real(dl) :: a_exact_end = 1._dl
        real(dl) :: a_start = 1e-10_dl
        real(dl) :: f_effective = 0._dl
        real(dl) :: grho_ax_today = 0._dl
        real(dl) :: initial_phi = 0._dl
        logical :: background_ready = .false.
        class(CAMBdata), pointer, private :: bg_state => null()
        real(dl), allocatable, private :: exact_a(:), exact_phi(:), exact_phidot(:)
        real(dl), allocatable, private :: exact_ddphi(:), exact_ddphidot(:)
        real(dl), allocatable, private :: fluid_a(:), fluid_rho(:), fluid_w(:), fluid_dwdloga(:)
        real(dl), allocatable, private :: fluid_ddrho(:), fluid_ddw(:), fluid_dddwdloga(:)
    contains
    procedure :: ReadParams => TFuzzyDM_ReadParams
    procedure, nopass :: PythonClass => TFuzzyDM_PythonClass
    procedure, nopass :: SelfPointer => TFuzzyDM_SelfPointer
    procedure :: Init => TFuzzyDM_Init
    procedure :: Validate => TFuzzyDM_Validate
    procedure :: BackgroundDensityAndPressure => TFuzzyDM_BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy => TFuzzyDM_PerturbedStressEnergy
    procedure :: PerturbationInitial => TFuzzyDM_PerturbationInitial
    procedure :: PerturbationEvolve => TFuzzyDM_PerturbationEvolve
    procedure :: has_switch => TFuzzyDM_has_switch
    procedure :: Switch => TFuzzyDM_Switch
    procedure :: ScalarSwitchScaleFactor => TFuzzyDM_ScalarSwitchScaleFactor
    procedure :: PrintFeedback => TFuzzyDM_PrintFeedback
    procedure :: cs2_eff => TFuzzyDM_cs2_eff
    procedure, private :: ExactValsAt => TFuzzyDM_ExactValsAt
    procedure, private :: FluidValsAt => TFuzzyDM_FluidValsAt
    procedure, private :: FluidHVals => TFuzzyDM_FluidHVals
    procedure, private :: CalcAuxiliary => TFuzzyDM_CalcAuxiliary
    procedure, private :: EvaluateAmplitude => TFuzzyDM_EvaluateAmplitude
    procedure, private :: BuildBackground => TFuzzyDM_BuildBackground
    procedure, private :: ExactHAt => TFuzzyDM_ExactHAt
    procedure, private :: RefineMatchSurface => TFuzzyDM_RefineMatchSurface
    end type TFuzzyDM

    procedure(TClassDverk) :: dverk
    public TFuzzyDM

    contains

    subroutine TFuzzyDM_Validate(this, OK)
    class(TFuzzyDM), intent(in) :: this
    logical, intent(inout) :: OK
    logical :: active

    call this%TDarkMatterModel%Validate(OK)
    active = this%f_axion > 0._dl .or. this%omega_axion_h2 > 0._dl
    if (.not. ieee_is_finite(this%m_axion) .or. .not. ieee_is_finite(this%f_axion) .or. &
        .not. ieee_is_finite(this%omega_axion_h2) .or. this%f_axion < 0._dl .or. &
        this%omega_axion_h2 < 0._dl .or. &
        (this%f_axion > 0._dl .and. this%omega_axion_h2 > 0._dl)) then
        OK = .false.
        write(*,*) 'FuzzyDM parameters must be finite, non-negative, and specify only one abundance.'
    end if
    if (active .and. (this%m_axion < 1e-23_dl .or. this%m_axion > 1e-21_dl)) then
        OK = .false.
        write(*,*) 'Active FuzzyDM requires 1e-23 <= m_axion/eV <= 1e-21.'
    end if
    if (this%f_axion > 0._dl .and. (this%f_axion < 1e-3_dl .or. this%f_axion > 0.1_dl)) then
        OK = .false.
        write(*,*) 'FuzzyDM f_axion must be in the supported guard range 1e-3 <= f_axion <= 0.1.'
    end if
    if (active .and. (this%match_ratio < 50 .or. this%match_ratio > 75)) then
        OK = .false.
        write(*,*) 'FuzzyDM fdm_match_ratio must be in the supported range 50...75.'
    end if
    end subroutine TFuzzyDM_Validate

    subroutine TFuzzyDM_ReadParams(this, Ini)
    use IniObjects
    class(TFuzzyDM) :: this
    class(TIniFile), intent(in) :: Ini
    this%m_axion = Ini%Read_Double('m_axion', 0._dl)
    this%omega_axion_h2 = Ini%Read_Double('omega_axion_h2', 0._dl)
    this%f_axion = Ini%Read_Double('f_axion', 0._dl)
    this%match_ratio = Ini%Read_Int('fdm_match_ratio', default_match_ratio)
    end subroutine TFuzzyDM_ReadParams

    function TFuzzyDM_PythonClass()
    character(LEN=:), allocatable :: TFuzzyDM_PythonClass
    TFuzzyDM_PythonClass = 'FuzzyDM'
    end function TFuzzyDM_PythonClass

    subroutine TFuzzyDM_SelfPointer(cptr, P)
    use iso_c_binding
    type(c_ptr) :: cptr
    type(TFuzzyDM), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P
    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TFuzzyDM_SelfPointer

    function TFuzzyDM_cs2_eff(this, k, a) result(cs2)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: k, a
    real(dl) :: cs2, alpha, H, dHdt, w, dwdloga, rho

    alpha = (k / max(a*this%m_conf, 1e-60_dl))**2
    if (alpha > 1e-10_dl) then
        cs2 = (sqrt(1._dl + alpha) - 1._dl)**2 / alpha
    else
        cs2 = k**2 / (k**2 + 4._dl*this%m_conf**2*a**2)
    end if
    if (allocated(this%fluid_a)) then
        call this%FluidValsAt(max(a,this%a_match), rho, w, dwdloga)
    else
        rho = this%grho_ax_today/max(a,1e-30_dl)
    end if
    call this%FluidHVals(a, rho, H, dHdt, w, dwdloga)
    cs2 = min(1._dl, max(0._dl, cs2 + 1.25_dl*(H/this%m_conf)**2))
    end function TFuzzyDM_cs2_eff

    subroutine TFuzzyDM_Init(this, State)
    class(TFuzzyDM), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    logical :: OK
    real(dl) :: alo, ahi, amid, H, grhoa2, new_match, match_residual
    integer :: i, match_iteration

    OK = .true.
    call this%Validate(OK)
    if (.not. OK) then
        call GlobalError('FuzzyDM parameters are outside the KG--EFA validation contract.', error_unsupported_params)
        return
    end if
    this%background_ready = .false.
    this%f_effective = 0._dl
    this%is_standard_cdm = (this%f_axion == 0._dl .and. this%omega_axion_h2 == 0._dl)
    this%num_dm_equations = 0
    this%num_dr_equations = 0
    this%has_cdm_velocity = .false.
    this%has_background_pressure = .false.
    this%num_perturb_equations = 0
    if (this%is_standard_cdm) return

    select type(S => State)
    class is (CAMBdata)
        this%bg_state => S
        if (.not. S%CP%DarkEnergy%is_cosmological_constant) then
            call GlobalError('Exact FuzzyDM currently requires a cosmological-constant dark-energy background.', &
                error_unsupported_params)
            return
        end if
        if (S%CP%omch2 <= 0._dl) then
            call GlobalError('FuzzyDM requires positive total omch2.', error_unsupported_params)
            return
        end if
        if (this%omega_axion_h2 > 0._dl) then
            this%f_effective = this%omega_axion_h2/S%CP%omch2
        else
            this%f_effective = this%f_axion
        end if
        if (this%f_effective < 1e-3_dl .or. this%f_effective > 0.1_dl) then
            call GlobalError('FuzzyDM abundance resolves outside 1e-3 <= f_axion <= 0.1.', error_unsupported_params)
            return
        end if
        ! Convert m [eV] to the conformal inverse-Mpc convention used by CAMB.
        this%m_conf = this%m_axion*eV*Mpc/(hbar*c)
        this%grho_ax_today = this%f_effective*S%grhoc
        this%num_dm_equations = 4
        this%num_perturb_equations = 4
        this%has_background_pressure = .true.

        alo = 1e-12_dl
        ahi = 1._dl
        do i = 1, 100
            amid = sqrt(alo*ahi)
            grhoa2 = S%grho_no_dm_de(amid) + (1._dl-this%f_effective)*S%grhoc*amid + &
                S%grhov*amid**4 + this%grho_ax_today*amid
            H = sqrt(grhoa2/3._dl)/amid**2
            if (this%m_conf/H > 3._dl) then
                ahi = amid
            else
                alo = amid
            end if
        end do
        this%a_osc = sqrt(alo*ahi)
        alo = this%a_osc
        ahi = 1._dl
        do i = 1, 100
            amid = sqrt(alo*ahi)
            grhoa2 = S%grho_no_dm_de(amid) + (1._dl-this%f_effective)*S%grhoc*amid + &
                S%grhov*amid**4 + this%grho_ax_today*amid
            H = sqrt(grhoa2/3._dl)/amid**2
            if (this%m_conf/H > real(this%match_ratio,dl)) then
                ahi = amid
            else
                alo = amid
            end if
        end do
        this%a_match = sqrt(alo*ahi)
        this%a_exact_end = min(1._dl,3._dl*this%a_match)
        this%a_start = max(1e-12_dl,min(1e-8_dl,this%a_osc/100._dl))
        ! The first estimate uses a dust-like axion density. Rebuild the
        ! shooting solution until the common surface is self-consistent with
        ! the exact KG contribution to H(a).
        do match_iteration = 1, 6
            call this%BuildBackground()
            if (global_error_flag /= 0) return
            call this%RefineMatchSurface(new_match,match_residual)
            if (global_error_flag /= 0) return
            if (match_residual <= 1e-8_dl) exit
            if (match_iteration == 6) then
                call GlobalError('FuzzyDM exact m/H matching surface did not converge.',error_unsupported_params)
                return
            end if
            this%a_match = new_match
            this%a_exact_end = min(1._dl,3._dl*this%a_match)
        end do
        this%background_ready = .true.
    end select
    end subroutine TFuzzyDM_Init

    subroutine FDM_EvolveBackgroundLog(self, num, loga, y, yprime)
    class(TFuzzyDM) :: self
    integer, intent(in) :: num
    real(dl), intent(in) :: loga, y(num)
    real(dl), intent(out) :: yprime(num)
    real(dl) :: a, a2, phidot, rho_ax_t, grhoa2, adot
    a = exp(loga)
    yprime = 0._dl
    if (.not. ieee_is_finite(a) .or. a <= 0._dl .or. .not. all(ieee_is_finite(y))) then
        call GlobalError('FuzzyDM KG background integration produced a non-finite state.',error_unsupported_params)
        return
    end if
    a2 = a*a
    phidot = y(2)/a2
    rho_ax_t = 0.5_dl*phidot**2 + 0.5_dl*a2*self%m_conf**2*y(1)**2
    grhoa2 = self%bg_state%grho_no_dm_de(a) + &
        (1._dl-self%f_effective)*self%bg_state%grhoc*a + self%bg_state%grhov*a**4 + a2*rho_ax_t
    if (.not. ieee_is_finite(rho_ax_t) .or. rho_ax_t <= 0._dl .or. &
        .not. ieee_is_finite(grhoa2) .or. grhoa2 <= 0._dl) then
        call GlobalError('FuzzyDM KG background integration encountered a non-positive density.', &
            error_unsupported_params)
        return
    end if
    adot = sqrt(grhoa2/3._dl)
    yprime(1) = a*phidot/adot
    yprime(2) = -a**5*self%m_conf**2*y(1)/adot
    end subroutine FDM_EvolveBackgroundLog

    subroutine TFuzzyDM_EvaluateAmplitude(this, phi0, final_rho, save_tables)
    class(TFuzzyDM), intent(inout) :: this
    real(dl), intent(in) :: phi0
    real(dl), intent(out) :: final_rho
    logical, intent(in) :: save_tables
    integer, parameter :: neq = 2
    real(dl) :: c(24), work(neq,9), y(neq), afrom, ato, a
    real(dl) :: rho_match, p_match, phic, phis, dphicdt, dphisdt, D, H
    integer :: ind, i

    final_rho = 0._dl
    if (.not. ieee_is_finite(phi0) .or. phi0 <= 0._dl) then
        call GlobalError('FuzzyDM abundance shooting received an invalid field amplitude.',error_unsupported_params)
        return
    end if
    y = [phi0,0._dl]
    afrom = log(this%a_start)
    ato = log(this%a_match)
    ind = 1
    call dverk(this,neq,FDM_EvolveBackgroundLog,afrom,y,ato,integrate_tol,ind,c,neq,work)
    if (global_error_flag /= 0) return
    if (.not. all(ieee_is_finite(y))) then
        call GlobalError('FuzzyDM abundance shooting returned a non-finite KG state.',error_unsupported_params)
        return
    end if
    call this%CalcAuxiliary(this%a_match,y(1),y(2)/this%a_match**2,rho_match,p_match, &
        phic,phis,dphicdt,dphisdt,D,H)
    if (global_error_flag /= 0) return
    if (.not. all(ieee_is_finite([rho_match,p_match,phic,phis,dphicdt,dphisdt,D,H])) .or. &
        rho_match <= 0._dl) then
        call GlobalError('FuzzyDM matching map returned an invalid density or phase envelope.', &
            error_unsupported_params)
        return
    end if
    call EvolveFluidToToday(this,rho_match,final_rho)
    if (global_error_flag /= 0) return
    if (.not. ieee_is_finite(final_rho) .or. final_rho <= 0._dl) then
        call GlobalError('FuzzyDM abundance shooting returned an invalid final density.',error_unsupported_params)
        return
    end if
    if (.not. save_tables) return

    if (allocated(this%exact_a)) &
        deallocate(this%exact_a,this%exact_phi,this%exact_phidot,this%exact_ddphi,this%exact_ddphidot)
    allocate(this%exact_a(background_points),this%exact_phi(background_points), &
        this%exact_phidot(background_points),this%exact_ddphi(background_points), &
        this%exact_ddphidot(background_points))
    y = [phi0,0._dl]
    afrom = log(this%a_start)
    this%exact_a(1) = this%a_start
    this%exact_phi(1) = y(1)
    this%exact_phidot(1) = 0._dl
    ind = 1
    do i = 2, background_points
        ato = log(this%a_start)+real(i-1,dl)/real(background_points-1,dl)*log(this%a_exact_end/this%a_start)
        call dverk(this,neq,FDM_EvolveBackgroundLog,afrom,y,ato,integrate_tol,ind,c,neq,work)
        if (global_error_flag /= 0) return
        a = exp(ato)
        if (.not. ieee_is_finite(a) .or. .not. all(ieee_is_finite(y)) .or. a <= 0._dl) then
            call GlobalError('FuzzyDM exact background table contains a non-finite state.',error_unsupported_params)
            return
        end if
        this%exact_a(i) = a
        this%exact_phi(i) = y(1)
        this%exact_phidot(i) = y(2)/a**2
    end do
    call spline(this%exact_a,this%exact_phi,background_points,0._dl,0._dl,this%exact_ddphi)
    call spline(this%exact_a,this%exact_phidot,background_points,0._dl,0._dl,this%exact_ddphidot)
    call BuildFluidTable(this,rho_match)
    end subroutine TFuzzyDM_EvaluateAmplitude

    subroutine TFuzzyDM_BuildBackground(this)
    class(TFuzzyDM), intent(inout) :: this
    real(dl) :: phi,rho_final,loglo,loghi,logmid,reslo,reshi,resmid
    integer :: i
    phi = sqrt(2._dl*this%grho_ax_today/(max(this%a_osc,1e-30_dl)**3*this%m_conf**2))
    loglo = log(phi)-log(4._dl)
    loghi = log(phi)+log(4._dl)
    call this%EvaluateAmplitude(exp(loglo),rho_final,.false.)
    if (global_error_flag /= 0) return
    if (.not. ieee_is_finite(rho_final) .or. rho_final <= 0._dl) then
        call GlobalError('FuzzyDM lower shooting bracket has an invalid density.',error_unsupported_params)
        return
    end if
    reslo = log(rho_final/this%grho_ax_today)
    call this%EvaluateAmplitude(exp(loghi),rho_final,.false.)
    if (global_error_flag /= 0) return
    if (.not. ieee_is_finite(rho_final) .or. rho_final <= 0._dl) then
        call GlobalError('FuzzyDM upper shooting bracket has an invalid density.',error_unsupported_params)
        return
    end if
    reshi = log(rho_final/this%grho_ax_today)
    do i = 1, 8
        if (reslo < 0._dl .and. reshi > 0._dl) exit
        if (reslo >= 0._dl) then
            loglo = loglo-log(4._dl)
            call this%EvaluateAmplitude(exp(loglo),rho_final,.false.)
            if (global_error_flag /= 0) return
            if (.not. ieee_is_finite(rho_final) .or. rho_final <= 0._dl) then
                call GlobalError('FuzzyDM expanded lower bracket has an invalid density.',error_unsupported_params)
                return
            end if
            reslo = log(rho_final/this%grho_ax_today)
        end if
        if (reshi <= 0._dl) then
            loghi = loghi+log(4._dl)
            call this%EvaluateAmplitude(exp(loghi),rho_final,.false.)
            if (global_error_flag /= 0) return
            if (.not. ieee_is_finite(rho_final) .or. rho_final <= 0._dl) then
                call GlobalError('FuzzyDM expanded upper bracket has an invalid density.',error_unsupported_params)
                return
            end if
            reshi = log(rho_final/this%grho_ax_today)
        end if
    end do
    if (.not. ieee_is_finite(reslo) .or. .not. ieee_is_finite(reshi) .or. &
        reslo >= 0._dl .or. reshi <= 0._dl) then
        call GlobalError('FuzzyDM could not bracket the log-abundance shooting residual.',error_unsupported_params)
        return
    end if
    do i = 1, 60
        logmid = 0.5_dl*(loglo+loghi)
        call this%EvaluateAmplitude(exp(logmid),rho_final,.false.)
        if (global_error_flag /= 0) return
        if (.not. ieee_is_finite(rho_final) .or. rho_final <= 0._dl) then
            call GlobalError('FuzzyDM bisection returned an invalid density.',error_unsupported_params)
            return
        end if
        resmid = log(rho_final/this%grho_ax_today)
        if (abs(resmid) < 2e-8_dl) exit
        if (resmid > 0._dl) then
            loghi = logmid
        else
            loglo = logmid
        end if
    end do
    this%initial_phi = exp(logmid)
    call this%EvaluateAmplitude(this%initial_phi,rho_final,.true.)
    if (global_error_flag /= 0) return
    if (.not. ieee_is_finite(rho_final) .or. rho_final <= 0._dl .or. &
        abs(rho_final/this%grho_ax_today-1._dl) > 1e-6_dl) &
        call GlobalError('FuzzyDM abundance shooting failed its 1e-6 closure test.',error_unsupported_params)
    end subroutine TFuzzyDM_BuildBackground

    subroutine EvolveFluidToToday(this,rho_start,rho_final)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: rho_start
    real(dl), intent(out) :: rho_final
    real(dl) :: x,dx,rho,k1,k2,k3,k4
    integer :: i
    rho_final = 0._dl
    if (.not. ieee_is_finite(rho_start) .or. rho_start <= 0._dl) then
        call GlobalError('FuzzyDM EFA received an invalid matching density.',error_unsupported_params)
        return
    end if
    x = log(this%a_match)
    dx = -x/real(fluid_points,dl)
    rho = rho_start
    do i = 1, fluid_points
        k1 = FluidDerivative(this,x,rho)
        k2 = FluidDerivative(this,x+dx/2,rho+dx*k1/2)
        k3 = FluidDerivative(this,x+dx/2,rho+dx*k2/2)
        k4 = FluidDerivative(this,x+dx,rho+dx*k3)
        rho = rho+dx*(k1+2*k2+2*k3+k4)/6
        if (.not. ieee_is_finite(rho) .or. rho <= 0._dl) then
            call GlobalError('FuzzyDM EFA integration produced an invalid density.',error_unsupported_params)
            return
        end if
        x = x+dx
    end do
    rho_final = rho
    end subroutine EvolveFluidToToday

    subroutine TFuzzyDM_ExactHAt(this,a,H)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: H
    real(dl) :: phi,phidot,rho,grhoa2

    H = 0._dl
    if (.not. ieee_is_finite(a) .or. a <= 0._dl) then
        call GlobalError('FuzzyDM exact H(a) received an invalid scale factor.',error_unsupported_params)
        return
    end if
    call this%ExactValsAt(a,phi,phidot)
    rho = 0.5_dl*phidot**2+0.5_dl*a*a*this%m_conf**2*phi**2
    grhoa2 = this%bg_state%grho_no_dm_de(a)+(1._dl-this%f_effective)*this%bg_state%grhoc*a + &
        this%bg_state%grhov*a**4+a*a*rho
    if (.not. all(ieee_is_finite([phi,phidot,rho,grhoa2])) .or. rho <= 0._dl .or. grhoa2 <= 0._dl) then
        call GlobalError('FuzzyDM exact H(a) encountered an invalid density.',error_unsupported_params)
        return
    end if
    H = sqrt(grhoa2/3._dl)/a**2
    end subroutine TFuzzyDM_ExactHAt

    subroutine TFuzzyDM_RefineMatchSurface(this,new_match,residual)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(out) :: new_match,residual
    real(dl) :: alo,ahi,amid,H,flo,fhi,fmid,target
    integer :: i

    target = real(this%match_ratio,dl)
    call this%ExactHAt(this%a_match,H)
    if (global_error_flag /= 0) return
    residual = abs(this%m_conf/H/target-1._dl)
    new_match = this%a_match
    if (residual <= 1e-8_dl) return

    alo = max(this%a_start*(1._dl+1e-8_dl),this%a_osc)
    ahi = this%a_exact_end
    call this%ExactHAt(alo,H)
    if (global_error_flag /= 0) return
    flo = this%m_conf/H-target
    call this%ExactHAt(ahi,H)
    if (global_error_flag /= 0) return
    fhi = this%m_conf/H-target
    if (.not. all(ieee_is_finite([flo,fhi])) .or. flo > 0._dl .or. fhi < 0._dl) then
        call GlobalError('FuzzyDM could not bracket the exact m/H matching surface.',error_unsupported_params)
        return
    end if
    do i = 1, 80
        amid = sqrt(alo*ahi)
        call this%ExactHAt(amid,H)
        if (global_error_flag /= 0) return
        fmid = this%m_conf/H-target
        if (.not. ieee_is_finite(fmid)) then
            call GlobalError('FuzzyDM exact m/H root returned a non-finite residual.',error_unsupported_params)
            return
        end if
        if (fmid > 0._dl) then
            ahi = amid
        else
            alo = amid
        end if
    end do
    new_match = sqrt(alo*ahi)
    end subroutine TFuzzyDM_RefineMatchSurface

    function FluidDerivative(this,loga,rho) result(drhodloga)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: loga,rho
    real(dl) :: drhodloga,H,dHdt,w,dwdloga
    call this%FluidHVals(exp(loga),rho,H,dHdt,w,dwdloga)
    drhodloga = -(1._dl+3._dl*w)*rho
    end function FluidDerivative

    subroutine BuildFluidTable(this,rho_start)
    class(TFuzzyDM), intent(inout) :: this
    real(dl), intent(in) :: rho_start
    real(dl) :: x,dx,rho,k1,k2,k3,k4,H,dHdt
    integer :: i
    if (.not. ieee_is_finite(rho_start) .or. rho_start <= 0._dl) then
        call GlobalError('FuzzyDM fluid table received an invalid matching density.',error_unsupported_params)
        return
    end if
    if (allocated(this%fluid_a)) deallocate(this%fluid_a,this%fluid_rho,this%fluid_w, &
        this%fluid_dwdloga,this%fluid_ddrho,this%fluid_ddw,this%fluid_dddwdloga)
    allocate(this%fluid_a(fluid_points+1),this%fluid_rho(fluid_points+1), &
        this%fluid_w(fluid_points+1),this%fluid_dwdloga(fluid_points+1), &
        this%fluid_ddrho(fluid_points+1),this%fluid_ddw(fluid_points+1), &
        this%fluid_dddwdloga(fluid_points+1))
    x = log(this%a_match)
    dx = -x/real(fluid_points,dl)
    rho = rho_start
    do i = 1, fluid_points+1
        this%fluid_a(i) = exp(x)
        this%fluid_rho(i) = rho
        call this%FluidHVals(exp(x),rho,H,dHdt,this%fluid_w(i),this%fluid_dwdloga(i))
        if (i <= fluid_points) then
            k1 = FluidDerivative(this,x,rho)
            k2 = FluidDerivative(this,x+dx/2,rho+dx*k1/2)
            k3 = FluidDerivative(this,x+dx/2,rho+dx*k2/2)
            k4 = FluidDerivative(this,x+dx,rho+dx*k3)
            rho = rho+dx*(k1+2*k2+2*k3+k4)/6
            if (.not. ieee_is_finite(rho) .or. rho <= 0._dl) then
                call GlobalError('FuzzyDM fluid table integration produced an invalid density.',error_unsupported_params)
                return
            end if
            x = x+dx
        end if
    end do
    call spline(this%fluid_a,this%fluid_rho,fluid_points+1,0._dl,0._dl,this%fluid_ddrho)
    call spline(this%fluid_a,this%fluid_w,fluid_points+1,0._dl,0._dl,this%fluid_ddw)
    call spline(this%fluid_a,this%fluid_dwdloga,fluid_points+1,0._dl,0._dl,this%fluid_dddwdloga)
    end subroutine BuildFluidTable

    subroutine TFuzzyDM_FluidHVals(this,a,rho,H,dHdt,w,dwdloga)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: a,rho
    real(dl), intent(out) :: H,dHdt,w,dwdloga
    real(dl) :: grho_tot,gpres_tot
    grho_tot = this%bg_state%grho_no_dm_de(a) + &
        (1._dl-this%f_effective)*this%bg_state%grhoc*a + this%bg_state%grhov*a**4+a*a*rho
    H = sqrt(max(grho_tot,1e-60_dl)/3._dl)/a**2
    w = 1.5_dl*(H/this%m_conf)**2
    gpres_tot = this%bg_state%gpres_no_dm_de(a)-this%bg_state%grhov*a**4+a*a*w*rho
    dHdt = -0.5_dl*(grho_tot/3._dl+gpres_tot)/a**4-H**2
    dwdloga = 3._dl*dHdt/this%m_conf**2
    end subroutine TFuzzyDM_FluidHVals

    subroutine TFuzzyDM_CalcAuxiliary(this,a,phi,phidot,rho,p,phic,phis,dphicdt,dphisdt,D,H)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: a,phi,phidot
    real(dl), intent(out) :: rho,p,phic,phis,dphicdt,dphisdt,D,H
    real(dl) :: dphidt,dHdt,grho_tot,gpres_tot,phase_denom,phase_scale
    integer :: i
    dphidt = phidot/a
    rho = 0.5_dl*phidot**2+0.5_dl*a*a*this%m_conf**2*phi**2
    p = 0.5_dl*phidot**2-0.5_dl*a*a*this%m_conf**2*phi**2
    do i = 1, 3
        grho_tot = this%bg_state%grho_no_dm_de(a) + &
            (1._dl-this%f_effective)*this%bg_state%grhoc*a + this%bg_state%grhov*a**4+a*a*rho
        gpres_tot = this%bg_state%gpres_no_dm_de(a)-this%bg_state%grhov*a**4+a*a*p
        if (.not. all(ieee_is_finite([grho_tot,gpres_tot])) .or. grho_tot <= 0._dl) then
            call GlobalError('FuzzyDM phase map encountered an invalid total density.',error_unsupported_params)
            return
        end if
        H = sqrt(grho_tot/3._dl)/a**2
        dHdt = -0.5_dl*(grho_tot/3._dl+gpres_tot)/a**4-H**2
        D = -H/2._dl*(3._dl-2._dl*dHdt/H**2)
        phase_denom = D**2+3*D*H+4*this%m_conf**2
        phase_scale = max(abs(D**2),abs(3*D*H),4*this%m_conf**2,tiny(1._dl))
        if (.not. all(ieee_is_finite([H,dHdt,D,phase_denom,phase_scale])) .or. H <= 0._dl .or. &
            abs(phase_denom) <= 100._dl*epsilon(1._dl)*phase_scale) then
            call GlobalError('FuzzyDM phase map encountered a singular envelope denominator.', &
                error_unsupported_params)
            return
        end if
        phic = phi
        phis = (dphidt*((D+3*H)**2+4*this%m_conf**2)+6*H*this%m_conf**2*phi) / &
            (this%m_conf*phase_denom)
        dphisdt = 3*H*this%m_conf*(-2*dphidt+D*phi)/phase_denom
        dphicdt = -3*H*(D*dphidt+3*dphidt*H+2*this%m_conf**2*phi)/phase_denom
        rho = 0.5_dl*a*a*(0.5_dl*(dphisdt**2+dphicdt**2)+this%m_conf*(-phic*dphisdt+phis*dphicdt) + &
            this%m_conf**2*(phic**2+phis**2))
        p = 0.5_dl*a*a*(0.5_dl*(dphisdt**2+dphicdt**2)+this%m_conf*(-phic*dphisdt+phis*dphicdt))
        if (.not. all(ieee_is_finite([rho,p,phic,phis,dphicdt,dphisdt])) .or. rho <= 0._dl) then
            call GlobalError('FuzzyDM phase map produced an invalid envelope state.',error_unsupported_params)
            return
        end if
    end do
    end subroutine TFuzzyDM_CalcAuxiliary

    subroutine TFuzzyDM_ExactValsAt(this,a,phi,phidot)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: phi,phidot
    if (a <= this%a_start .or. .not. allocated(this%exact_a)) then
        phi = this%initial_phi
        phidot = 0._dl
    else
        phi = SplineValue(a,this%exact_a,this%exact_phi,this%exact_ddphi)
        phidot = SplineValue(a,this%exact_a,this%exact_phidot,this%exact_ddphidot)
    end if
    end subroutine TFuzzyDM_ExactValsAt

    subroutine TFuzzyDM_FluidValsAt(this,a,rho,w,dwdloga)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: rho,w,dwdloga
    rho = SplineValue(a,this%fluid_a,this%fluid_rho,this%fluid_ddrho)
    w = SplineValue(a,this%fluid_a,this%fluid_w,this%fluid_ddw)
    dwdloga = SplineValue(a,this%fluid_a,this%fluid_dwdloga,this%fluid_dddwdloga)
    end subroutine TFuzzyDM_FluidValsAt

    function SplineValue(x,xa,ya,y2a) result(value)
    real(dl), intent(in) :: x,xa(:),ya(:),y2a(:)
    real(dl) :: value,h,aa,bb
    integer :: lo,hi,mid,n
    n = size(xa)
    if (x <= xa(1)) then
        value = ya(1)
        return
    else if (x >= xa(n)) then
        value = ya(n)
        return
    end if
    lo = 1
    hi = n
    do while (hi-lo > 1)
        mid = (hi+lo)/2
        if (xa(mid) > x) then
            hi = mid
        else
            lo = mid
        end if
    end do
    h = xa(hi)-xa(lo)
    aa = (xa(hi)-x)/h
    bb = 1._dl-aa
    value = aa*ya(lo)+bb*ya(hi)+((aa**3-aa)*y2a(lo)+(bb**3-bb)*y2a(hi))*h*h/6._dl
    end function SplineValue

    function ExactBlendWeight(this,a,dweight_dloga) result(weight)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out), optional :: dweight_dloga
    real(dl) :: weight,u,width
    if (a <= this%a_match) then
        weight = 1._dl
        if (present(dweight_dloga)) dweight_dloga = 0._dl
    else if (a >= this%a_exact_end) then
        weight = 0._dl
        if (present(dweight_dloga)) dweight_dloga = 0._dl
    else
        width = log(this%a_exact_end/this%a_match)
        u = log(a/this%a_match)/width
        weight = 1._dl-3._dl*u**2+2._dl*u**3
        if (present(dweight_dloga)) dweight_dloga = (-6._dl*u+6._dl*u**2)/width
    end if
    end function ExactBlendWeight

    subroutine AxionBackgroundAt(this,a,rho,p)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: rho,p
    real(dl) :: phi,phidot,rho_fluid,w_fluid,dwdloga,weight,dweight,aeval
    real(dl) :: rho_exact,p_exact,drho_exact,drho_fluid,drho_mix
    aeval = min(a,this%a_exact_end)
    call this%ExactValsAt(aeval,phi,phidot)
    rho_exact = 0.5_dl*phidot**2+0.5_dl*aeval**2*this%m_conf**2*phi**2
    p_exact = 0.5_dl*phidot**2-0.5_dl*aeval**2*this%m_conf**2*phi**2
    if (.not. all(ieee_is_finite([rho_exact,p_exact])) .or. rho_exact <= 0._dl) then
        rho = 0._dl
        p = 0._dl
        call GlobalError('FuzzyDM background interpolation returned an invalid exact density.', &
            error_unsupported_params)
        return
    end if
    rho = rho_exact
    p = p_exact
    if (a >= this%a_match) then
        call this%FluidValsAt(a,rho_fluid,w_fluid,dwdloga)
        if (.not. all(ieee_is_finite([rho_fluid,w_fluid,dwdloga])) .or. rho_fluid <= 0._dl) then
            rho = 0._dl
            p = 0._dl
            call GlobalError('FuzzyDM background interpolation returned an invalid EFA density.', &
                error_unsupported_params)
            return
        end if
        weight = ExactBlendWeight(this,a,dweight)
        drho_exact = -(1._dl+3._dl*p_exact/rho_exact)*rho_exact
        drho_fluid = -(1._dl+3._dl*w_fluid)*rho_fluid
        rho = weight*rho_exact+(1._dl-weight)*rho_fluid
        drho_mix = weight*drho_exact+(1._dl-weight)*drho_fluid+dweight*(rho_exact-rho_fluid)
        ! Derive pressure from continuity so the blend, including W', is
        ! conservative rather than independently interpolating rho and p.
        p = (-drho_mix-rho)/3._dl
    end if
    end subroutine AxionBackgroundAt

    subroutine TFuzzyDM_BackgroundDensityAndPressure(this,grhodm,a,grhodm_t,gpres_dm)
    class(TFuzzyDM), intent(inout) :: this
    real(dl), intent(in) :: grhodm,a
    real(dl), intent(out) :: grhodm_t
    real(dl), optional, intent(out) :: gpres_dm
    real(dl) :: rho_ax,p_ax
    if (a <= 0._dl) then
        grhodm_t = 0._dl
        if (present(gpres_dm)) gpres_dm = 0._dl
    else if (this%is_standard_cdm) then
        grhodm_t = grhodm/a
        if (present(gpres_dm)) gpres_dm = 0._dl
    else if (.not. this%background_ready) then
        grhodm_t = (1._dl-this%f_effective)*grhodm/a
        if (present(gpres_dm)) gpres_dm = 0._dl
    else
        call AxionBackgroundAt(this,a,rho_ax,p_ax)
        grhodm_t = (1._dl-this%f_effective)*grhodm/a+rho_ax
        if (present(gpres_dm)) gpres_dm = p_ax
    end if
    end subroutine TFuzzyDM_BackgroundDensityAndPressure

    function TFuzzyDM_has_switch(this,a) result(has_switch)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: a
    logical :: has_switch
    has_switch = .not. this%is_standard_cdm .and. this%background_ready .and. a >= this%a_match
    end function TFuzzyDM_has_switch

    function TFuzzyDM_ScalarSwitchScaleFactor(this) result(a_switch)
    class(TFuzzyDM), intent(in) :: this
    real(dl) :: a_switch
    if (this%background_ready) then
        a_switch = this%a_match
    else
        a_switch = 0._dl
    end if
    end function TFuzzyDM_ScalarSwitchScaleFactor

    subroutine TFuzzyDM_Switch(this,dm_ix,a,k,z,y)
    class(TFuzzyDM), intent(inout) :: this
    integer, intent(in) :: dm_ix
    real(dl), intent(in) :: a,k,z
    real(dl), intent(inout) :: y(:)
    real(dl) :: phi,phidot,rho,p,phic,phis,dphicdt,dphisdt,D,H
    real(dl) :: delphi,ddelphidt,dhsdt,delphic,delphis,ddelphicdt,ddelphisdt
    real(dl) :: dgrhoe,dgqe,denom,a2,m,q,qscale,delta_efa,heat_efa
    if (dm_ix <= 0) return
    y(dm_ix+2:dm_ix+3) = 0._dl
    if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(k) .or. .not. ieee_is_finite(z) .or. &
        a <= 0._dl .or. .not. all(ieee_is_finite(y(dm_ix:dm_ix+1)))) then
        call GlobalError('FuzzyDM KG-to-EFA switch received a non-finite state.',error_unsupported_params)
        return
    end if
    a2 = a*a
    m = this%m_conf
    delphi = y(dm_ix)
    ddelphidt = y(dm_ix+1)/a
    dhsdt = 2._dl*k*z/a
    call this%ExactValsAt(a,phi,phidot)
    call this%CalcAuxiliary(a,phi,phidot,rho,p,phic,phis,dphicdt,dphisdt,D,H)
    if (global_error_flag /= 0) return
    if (.not. all(ieee_is_finite([phi,phidot,rho,p,phic,phis,dphicdt,dphisdt,D,H])) .or. &
        rho <= 0._dl .or. m <= 0._dl) then
        call GlobalError('FuzzyDM KG-to-EFA switch has an invalid background state.',error_unsupported_params)
        return
    end if
    delphic = delphi
    q = 2*k**2+a2*(D**2+3*D*H+4*m**2)
    qscale = max(2*k**2,a2*(abs(D**2)+abs(3*D*H)+4*m**2),tiny(1._dl))
    denom = 2*m*q
    if (.not. all(ieee_is_finite([q,qscale,denom])) .or. &
        abs(q) <= 100._dl*epsilon(1._dl)*qscale .or. abs(denom) <= tiny(1._dl)) then
        call GlobalError('FuzzyDM KG-to-EFA switch encountered a singular phase-map denominator.', &
            error_unsupported_params)
        return
    end if
    delphis = (2*delphi*(D+3*H)*k**2+a2*(2*D**2*ddelphidt+12*D*ddelphidt*H+18*ddelphidt*H**2 + &
        4*(2*ddelphidt+3*delphi*H)*m**2+2*dhsdt*m*(-dphisdt+m*phic) + &
        D*dhsdt*(dphicdt+m*phis)+3*dhsdt*H*(dphicdt+m*phis)))/denom
    ddelphicdt = (2*(2*ddelphidt-delphi*(D+3*H))*k**2-a2*(6*D*ddelphidt*H+3*dhsdt*dphicdt*H + &
        6*H*(3*ddelphidt*H+2*delphi*m**2)+dhsdt*m*(-2*dphisdt+2*m*phic+3*H*phis) + &
        D*dhsdt*(dphicdt+m*phis)))/(4*k**2+2*a2*(D**2+3*D*H+4*m**2))
    ddelphisdt = -0.5_dl*(2*delphi*k**4+a2*k**2*(2*D*ddelphidt+6*ddelphidt*H+4*delphi*m**2 + &
        dhsdt*(dphicdt+m*phis))+a2**2*m*(12*ddelphidt*H*m-6*D*delphi*H*m + &
        D*dhsdt*(dphisdt-m*phic)+2*dhsdt*m*(dphicdt+m*phis))) / &
        (a2*m*(2*k**2+a2*(D**2+3*D*H+4*m**2)))
    dgrhoe = 0.5_dl*a2*(dphicdt*ddelphicdt+dphisdt*ddelphisdt+m*(phis*ddelphicdt-phic*ddelphisdt) + &
        m*(delphis*dphicdt-delphic*dphisdt)+2*m**2*(phis*delphis+phic*delphic))
    dgqe = 0.5_dl*k*a*(m*(delphic*phis-delphis*phic)+delphic*dphicdt+delphis*dphisdt)
    delta_efa = dgrhoe/rho
    heat_efa = dgqe/rho
    if (.not. all(ieee_is_finite([delphis,ddelphicdt,ddelphisdt,dgrhoe,dgqe,delta_efa,heat_efa]))) then
        call GlobalError('FuzzyDM KG-to-EFA switch produced a non-finite mapped state.',error_unsupported_params)
        return
    end if
    y(dm_ix+2) = delta_efa
    y(dm_ix+3) = heat_efa
    end subroutine TFuzzyDM_Switch

    subroutine TFuzzyDM_PerturbationInitial(this,y,a,tau,k,dm_ix,dr_ix,vc_ix)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    real(dl), intent(in) :: a,tau,k
    integer, intent(in) :: dm_ix,dr_ix,vc_ix
    if (.not. this%is_standard_cdm .and. dm_ix > 0) y(dm_ix:dm_ix+3) = 0._dl
    end subroutine TFuzzyDM_PerturbationInitial

    subroutine TFuzzyDM_PerturbedStressEnergy(this,dgrhoe,dgqe,a,dgq,dgrho,grho,grhoc_t,adotoa,k, &
        ay,ayprime,dm_ix,dr_ix,vc_ix)
    class(TFuzzyDM), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe,dgqe
    real(dl), intent(in) :: a,dgq,dgrho,grho,grhoc_t,adotoa,k
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: dm_ix,dr_ix,vc_ix
    real(dl) :: phi,phidot,rho_ax,p_ax,rho_fluid,w,dwdloga,weight,dgrho_exact,dgq_exact
    dgrhoe = 0._dl
    dgqe = 0._dl
    if (this%is_standard_cdm .or. dm_ix <= 0) return
    call AxionBackgroundAt(this,a,rho_ax,p_ax)
    call this%ExactValsAt(min(a,this%a_exact_end),phi,phidot)
    dgrho_exact = phidot*ay(dm_ix+1)+a*a*this%m_conf**2*phi*ay(dm_ix)
    dgq_exact = k*phidot*ay(dm_ix)
    dgrhoe = dgrho_exact
    dgqe = dgq_exact
    if (a >= this%a_match) then
        call this%FluidValsAt(a,rho_fluid,w,dwdloga)
        weight = ExactBlendWeight(this,a)
        dgrhoe = weight*dgrho_exact+(1._dl-weight)*rho_fluid*ay(dm_ix+2)
        dgqe = weight*dgq_exact+(1._dl-weight)*rho_fluid*ay(dm_ix+3)
    end if
    dgrhoe = dgrhoe-rho_ax*ay(2)
    end subroutine TFuzzyDM_PerturbedStressEnergy

    subroutine TFuzzyDM_PerturbationEvolve(this,ayprime,a,adotoa,k,z,y,dm_ix,dr_ix,vc_ix,clxc,vb, &
        grhoc_t,grhob_t,sigma,high_ktau_dr)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a,adotoa,k,z,y(:)
    integer, intent(in) :: dm_ix,dr_ix,vc_ix
    real(dl), intent(in) :: clxc,vb,grhoc_t,grhob_t,sigma
    logical, intent(in), optional :: high_ktau_dr
    real(dl) :: phi,phidot,rho,w,dwdloga,cs2,deriv,hv3_over_k
    if (this%is_standard_cdm .or. dm_ix <= 0) return
    call this%ExactValsAt(min(a,this%a_exact_end),phi,phidot)
    if (a <= this%a_exact_end) then
        ayprime(dm_ix) = y(dm_ix+1)
        ayprime(dm_ix+1) = -2*adotoa*y(dm_ix+1)-k*z*phidot-(k**2+a*a*this%m_conf**2)*y(dm_ix)
    else
        ayprime(dm_ix:dm_ix+1) = 0._dl
    end if
    if (a >= this%a_match) then
        call this%FluidValsAt(a,rho,w,dwdloga)
        deriv = dwdloga/(1._dl+w)
        cs2 = this%cs2_eff(k,a)
        hv3_over_k = 3*adotoa*y(dm_ix+3)/k
        ayprime(dm_ix+2) = -3*adotoa*(cs2-w)*(y(dm_ix+2)+hv3_over_k)-k*y(dm_ix+3) - &
            (1._dl+w)*k*z-adotoa*deriv*hv3_over_k
        ayprime(dm_ix+3) = -adotoa*(1._dl-3*cs2-deriv)*y(dm_ix+3)+k*cs2*y(dm_ix+2)
    else
        ayprime(dm_ix+2:dm_ix+3) = 0._dl
    end if
    end subroutine TFuzzyDM_PerturbationEvolve

    subroutine TFuzzyDM_PrintFeedback(this,FeedbackLevel)
    class(TFuzzyDM) :: this
    integer, intent(in) :: FeedbackLevel
    if (FeedbackLevel > 0 .and. .not. this%is_standard_cdm) then
        write(*,'("Fuzzy DM exact KG -> PH-EFA: m_axion = ",ES12.4," eV")') this%m_axion
        write(*,'("          resolved f_axion = ",F8.5)') this%f_effective
        write(*,'("          a_osc, a_match, m/H_match = ",2ES12.4,I5)') &
            this%a_osc,this%a_match,this%match_ratio
    end if
    end subroutine TFuzzyDM_PrintFeedback

    end module FuzzyDM
