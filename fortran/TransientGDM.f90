    module TransientGDM
    use DarkMatterInteraction
    use results
    use classes
    implicit none
    private

    type, extends(TDarkMatterModel) :: TTransientGDM
        ! Public parameters. Keep this order synchronized with camb/dark_matter.py.
        real(dl) :: f_X = 0._dl
        real(dl) :: w_p = 0._dl
        real(dl) :: kappa0_hmpc = 6.4_dl
        real(dl) :: a_i = 1.46e-7_dl
        real(dl) :: a_f = 1.46e-5_dl

        ! Derived values, initialized by Init.
        real(dl) :: a_p = 1.46e-6_dl
        real(dl) :: sigma_w = 0.7675283643313485_dl
        real(dl) :: hubble_h = 0.7_dl
    contains
    procedure :: ReadParams => TTransientGDM_ReadParams
    procedure, nopass :: PythonClass => TTransientGDM_PythonClass
    procedure, nopass :: SelfPointer => TTransientGDM_SelfPointer
    procedure :: Init => TTransientGDM_Init
    procedure :: BackgroundDensityAndPressure => TTransientGDM_BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy => TTransientGDM_PerturbedStressEnergy
    procedure :: PerturbationInitial => TTransientGDM_PerturbationInitial
    procedure :: PerturbationEvolve => TTransientGDM_PerturbationEvolve
    procedure :: PrintFeedback => TTransientGDM_PrintFeedback
    procedure :: w_of_a => TTransientGDM_w_of_a
    procedure :: ca2_of_a => TTransientGDM_ca2_of_a
    procedure :: cs2_of_a_k => TTransientGDM_cs2_of_a_k
    procedure :: density_ratio => TTransientGDM_density_ratio
    procedure :: fraction_at_a => TTransientGDM_fraction_at_a
    procedure :: NewtonianVelocityTransfer => TTransientGDM_NewtonianVelocityTransfer
    procedure :: ComponentNewtonianVelocityTransfer => TTransientGDM_ComponentNewtonianVelocityTransfer
    procedure :: NonNuDensityTransfer => TTransientGDM_NonNuDensityTransfer
    end type TTransientGDM

    public TTransientGDM

    contains

    subroutine TTransientGDM_ReadParams(this, Ini)
    use IniObjects
    class(TTransientGDM) :: this
    class(TIniFile), intent(in) :: Ini

    this%f_X = Ini%Read_Double('TransientGDM_f_X', 0._dl)
    this%w_p = Ini%Read_Double('TransientGDM_w_p', 0._dl)
    this%kappa0_hmpc = Ini%Read_Double('TransientGDM_kappa0_hmpc', 6.4_dl)
    this%a_i = Ini%Read_Double('TransientGDM_a_i', 1.46e-7_dl)
    this%a_f = Ini%Read_Double('TransientGDM_a_f', 1.46e-5_dl)

    end subroutine TTransientGDM_ReadParams

    function TTransientGDM_PythonClass()
    character(LEN=:), allocatable :: TTransientGDM_PythonClass
    TTransientGDM_PythonClass = 'TransientGDM'
    end function TTransientGDM_PythonClass

    subroutine TTransientGDM_SelfPointer(cptr, P)
    use iso_c_binding
    type(c_ptr) :: cptr
    type(TTransientGDM), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TTransientGDM_SelfPointer

    subroutine TTransientGDM_Init(this, State)
    class(TTransientGDM), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    if (this%f_X < 0._dl .or. this%f_X >= 1._dl) &
        error stop 'TransientGDM requires 0 <= f_X < 1'
    if (this%w_p <= -1._dl .or. this%w_p > 0._dl) &
        error stop 'TransientGDM requires -1 < w_p <= 0'
    if (this%kappa0_hmpc <= 0._dl) &
        error stop 'TransientGDM requires kappa0_hmpc > 0'
    if (this%a_i <= 0._dl .or. this%a_f <= this%a_i .or. this%a_f >= 1._dl) &
        error stop 'TransientGDM requires 0 < a_i < a_f < 1'

    this%a_p = sqrt(this%a_i * this%a_f)
    this%sigma_w = log(this%a_f / this%a_i) / 6._dl

    select type(S => State)
    class is (CAMBdata)
        this%hubble_h = S%CP%H0 / 100._dl
    end select

    this%is_standard_cdm = (this%f_X == 0._dl .or. this%w_p == 0._dl)
    if (this%is_standard_cdm) then
        this%has_cdm_velocity = .false.
        this%num_dr_equations = 0
        this%num_perturb_equations = 0
    else
        ! vc_ix stores v_X (theta_X=k v_X); dm_ix stores delta_X.
        ! The current dark-matter interface reserves one inert dr_ix slot when
        ! num_dr_equations > 0. This also prevents the generic all-CDM velocity
        ! source from being added to the metric.
        this%has_cdm_velocity = .true.
        this%num_dr_equations = 1
        this%num_perturb_equations = 2
    end if

    end subroutine TTransientGDM_Init

    function TTransientGDM_w_of_a(this, a) result(w)
    class(TTransientGDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: w, x

    if (a <= 0._dl .or. this%w_p == 0._dl) then
        w = 0._dl
    else
        x = log(a / this%a_p) / this%sigma_w
        if (abs(x) > 38._dl) then
            w = 0._dl
        else
            w = this%w_p * exp(-0.5_dl * x*x)
        end if
    end if
    end function TTransientGDM_w_of_a

    function TTransientGDM_ca2_of_a(this, a) result(ca2)
    class(TTransientGDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: ca2, w, log_ratio

    w = this%w_of_a(a)
    if (w == 0._dl) then
        ca2 = 0._dl
    else
        log_ratio = log(a / this%a_p)
        ca2 = w + w * log_ratio / (3._dl * this%sigma_w**2 * (1._dl + w))
    end if
    end function TTransientGDM_ca2_of_a

    function TTransientGDM_cs2_of_a_k(this, a, k) result(cs2)
    class(TTransientGDM), intent(in) :: this
    real(dl), intent(in) :: a, k
    real(dl) :: cs2, kappa_mpc

    ! k is in Mpc^-1. kappa0 is supplied in h/Mpc.
    kappa_mpc = this%kappa0_hmpc * this%hubble_h / this%a_f
    cs2 = this%ca2_of_a(a) / (1._dl + (k / (a * kappa_mpc))**2)
    end function TTransientGDM_cs2_of_a_k

    function TTransientGDM_density_ratio(this, a) result(ratio)
    ! rho_X(a) / [rho_X0 a^-3], normalized to one today.
    class(TTransientGDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: ratio, integral_w, root2, prefactor, upper, lower

    if (a <= 0._dl) then
        ratio = 0._dl
        return
    end if
    root2 = sqrt(2._dl)
    prefactor = this%w_p * this%sigma_w * sqrt(acos(-1._dl) / 2._dl)
    upper = erf(log(1._dl / this%a_p) / (root2 * this%sigma_w))
    lower = erf(log(a / this%a_p) / (root2 * this%sigma_w))
    integral_w = prefactor * (upper - lower)
    ratio = exp(3._dl * integral_w)
    end function TTransientGDM_density_ratio

    function TTransientGDM_fraction_at_a(this, a) result(fraction)
    class(TTransientGDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: fraction, ratio, denominator

    ratio = this%density_ratio(a)
    denominator = (1._dl - this%f_X) + this%f_X * ratio
    fraction = this%f_X * ratio / denominator
    end function TTransientGDM_fraction_at_a

    subroutine TTransientGDM_BackgroundDensityAndPressure(this, grhodm, a, grhodm_t, gpres_dm)
    class(TTransientGDM), intent(inout) :: this
    real(dl), intent(in) :: grhodm, a
    real(dl), intent(out) :: grhodm_t
    real(dl), optional, intent(out) :: gpres_dm
    real(dl) :: ratio, rho_x_t

    if (a <= 0._dl) then
        grhodm_t = 0._dl
        if (present(gpres_dm)) gpres_dm = 0._dl
        return
    end if
    ratio = this%density_ratio(a)
    grhodm_t = grhodm / a * ((1._dl - this%f_X) + this%f_X * ratio)
    if (present(gpres_dm)) then
        rho_x_t = grhodm / a * this%f_X * ratio
        gpres_dm = rho_x_t * this%w_of_a(a)
    end if
    end subroutine TTransientGDM_BackgroundDensityAndPressure

    subroutine TTransientGDM_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhoc_t, adotoa, k, ay, ayprime, dm_ix, dr_ix, vc_ix)
    class(TTransientGDM), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhoc_t, adotoa, k
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl) :: fraction, w

    dgrhoe = 0._dl
    dgqe = 0._dl
    if (this%is_standard_cdm .or. dm_ix <= 0 .or. vc_ix <= 0) return

    fraction = this%fraction_at_a(a)
    w = this%w_of_a(a)
    ! The baseline source already treats all grhoc_t as standard CDM. Replace
    ! the X fraction by its own density and momentum perturbations.
    dgrhoe = grhoc_t * fraction * (ay(dm_ix) - ay(2))
    dgqe = grhoc_t * fraction * (1._dl + w) * ay(vc_ix)
    end subroutine TTransientGDM_PerturbedStressEnergy

    subroutine TTransientGDM_PerturbationInitial(this, y, a, tau, k, dm_ix, dr_ix, vc_ix)
    class(TTransientGDM), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer, intent(in) :: dm_ix, dr_ix, vc_ix

    if (this%is_standard_cdm) return
    if (dm_ix > 0) y(dm_ix) = y(2)
    if (vc_ix > 0) y(vc_ix) = 0._dl
    if (dr_ix > 0) y(dr_ix) = 0._dl
    end subroutine TTransientGDM_PerturbationInitial

    subroutine TTransientGDM_PerturbationEvolve(this, ayprime, a, adotoa, k, z, y, &
        dm_ix, dr_ix, vc_ix, clxc, vb, grhoc_t, grhob_t, sigma, high_ktau_dr)
    class(TTransientGDM), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, y(:)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl), intent(in) :: clxc, vb, grhoc_t, grhob_t, sigma
    logical, intent(in), optional :: high_ktau_dr
    real(dl) :: delta_x, v_x, w, ca2, cs2, one_plus_w

    if (this%is_standard_cdm .or. dm_ix <= 0 .or. vc_ix <= 0) return
    if (dr_ix > 0) ayprime(dr_ix) = 0._dl

    delta_x = y(dm_ix)
    v_x = y(vc_ix)
    w = this%w_of_a(a)
    ca2 = this%ca2_of_a(a)
    cs2 = this%cs2_of_a_k(a, k)
    one_plus_w = 1._dl + w

    ayprime(dm_ix) = -one_plus_w * k * (v_x + z) &
        - 3._dl * adotoa * (cs2 - w) * delta_x &
        - 9._dl * one_plus_w * (cs2 - ca2) * adotoa**2 * v_x / k
    ayprime(vc_ix) = -adotoa * (1._dl - 3._dl * cs2) * v_x &
        + k * cs2 * delta_x / one_plus_w
    end subroutine TTransientGDM_PerturbationEvolve

    function TTransientGDM_NewtonianVelocityTransfer(this, a, k, adotoa, sigma, ay, vc_ix) result(Tv)
    class(TTransientGDM), intent(in) :: this
    real(dl), intent(in) :: a, k, adotoa, sigma
    real(dl), intent(in) :: ay(*)
    integer, intent(in) :: vc_ix
    real(dl) :: Tv, x_momentum_fraction

    x_momentum_fraction = 0._dl
    if (.not. this%is_standard_cdm .and. vc_ix > 0) then
        x_momentum_fraction = this%fraction_at_a(a) * (1._dl + this%w_of_a(a))
    end if
    Tv = -k * (sigma + x_momentum_fraction * ay(vc_ix)) / adotoa
    end function TTransientGDM_NewtonianVelocityTransfer

    function TTransientGDM_ComponentNewtonianVelocityTransfer(this, a, k, adotoa, sigma, ay, vc_ix) result(Tv)
    class(TTransientGDM), intent(in) :: this
    real(dl), intent(in) :: a, k, adotoa, sigma
    real(dl), intent(in) :: ay(*)
    integer, intent(in) :: vc_ix
    real(dl) :: Tv

    if (this%is_standard_cdm .or. vc_ix <= 0) then
        Tv = 0._dl
    else
        Tv = -k * (sigma + ay(vc_ix)) / adotoa
    end if
    end function TTransientGDM_ComponentNewtonianVelocityTransfer

    function TTransientGDM_NonNuDensityTransfer(this, a, grhob_t, grhodm_t, clxb, clxc, ay, dm_ix) result(Tdelta)
    class(TTransientGDM), intent(in) :: this
    real(dl), intent(in) :: a, grhob_t, grhodm_t, clxb, clxc
    real(dl), intent(in) :: ay(*)
    integer, intent(in) :: dm_ix
    real(dl) :: Tdelta, fraction

    if (this%is_standard_cdm .or. dm_ix <= 0) then
        Tdelta = (grhob_t * clxb + grhodm_t * clxc) / (grhob_t + grhodm_t)
    else
        fraction = this%fraction_at_a(a)
        Tdelta = (grhob_t * clxb + grhodm_t * &
            ((1._dl - fraction) * clxc + fraction * ay(dm_ix))) / &
            (grhob_t + grhodm_t)
    end if
    end function TTransientGDM_NonNuDensityTransfer

    subroutine TTransientGDM_PrintFeedback(this, FeedbackLevel)
    class(TTransientGDM) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0 .and. .not. this%is_standard_cdm) then
        write(*,'("Transient GDM: f_X = ",ES12.4,"  w_p = ",F8.4)') this%f_X, this%w_p
        write(*,'("               kappa0 = ",F8.3," h/Mpc")') this%kappa0_hmpc
        write(*,'("               a_p = ",ES12.4,"  sigma_w = ",F8.5)') this%a_p, this%sigma_w
    end if
    end subroutine TTransientGDM_PrintFeedback

    end module TransientGDM
