    module DarkEnergyChaplygin
    use DarkEnergyInterface
    use DarkEnergyFluid
    use results
    use constants
    use classes
    use config, only: GlobalError, error_unsupported_params
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    ! Generalized Chaplygin gas (GCG) as the dark-energy component.
    !
    ! Equation of state   p = -A / rho^alpha   (alpha >= 0).
    ! Writing As = A / rho_de0^(1+alpha) in (0,1) the closed-form background is
    !
    !   rho_de(a) = rho_de0 * [ As + (1-As) a^(-3(1+alpha)) ]^(1/(1+alpha))
    !   w(a)      = -As / [ As + (1-As) a^(-3(1+alpha)) ]
    !   cs2(a)    = dp/drho = -alpha * w(a)   (>= 0 for w<0, alpha>=0)
    !
    ! cs2 is the adiabatic (rest-frame) sound speed fed to the fluid perturbation
    ! equations. No shooting is needed: both w(a) and cs2(a) are analytic, so we
    ! evaluate them directly on a log-a grid, fill the base-class tabulated-w
    ! machinery (equation_of_state + logdensity splines) and a companion cs2(a)
    ! spline read by the Fang-Hu-Lewis fluid perturbation equations.
    !
    ! Limits:
    !   As -> 1                   -> w -> -1 (pure Lambda) for any alpha.
    !   early (a -> 0)            -> w -> 0^- (dust-like), cs2 -> 0^+.
    !   today  (a = 1)            -> w = -As, cs2 = alpha*As.
    !
    ! For alpha -> 0 the gas reduces to LCDM plus a tiny early dust-like phase.
    ! For larger alpha the late-time sound speed cs2 = -alpha*w = O(alpha) grows,
    ! seeding the well-known small-scale oscillations/blow-up of unified GCG
    ! (Sandvik et al. 2004, astro-ph/0212114) — a property of the model, hence
    ! realistic unified-DM/DE usage needs alpha << 1.
    !
    ! (As, alpha) share the background convention lagRamses will reuse in its
    ! kessence_tabulate machinery (scalar_de_commons.f90); no Chaplygin section
    ! exists there yet, so the closed forms above are the matching contract.

    type, extends(TDarkEnergyFluid) :: TChaplygin
        real(dl) :: As    = 0.75_dl    ! A/rho_de0^(1+alpha), must be in (0,1)
        real(dl) :: alpha = 0.01_dl    ! GCG exponent, must be >= 0
        integer  :: n_tab = 4096       ! background table size (log-a grid; matches lagRamses)
        real(dl) :: cs2_0 = 0._dl      ! cs2 today (diagnostic, set in Init)
        Type(TCubicSpline) :: cs2ofa   ! rest-frame cs2 as function of log(a)
    contains
    procedure :: ReadParams => TChaplygin_ReadParams
    procedure :: Validate => TChaplygin_Validate
    procedure, nopass :: PythonClass => TChaplygin_PythonClass
    procedure, nopass :: SelfPointer => TChaplygin_SelfPointer
    procedure :: Init => TChaplygin_Init
    procedure :: PerturbationEvolve => TChaplygin_PerturbationEvolve
    procedure :: SetChaplyginParams => TChaplygin_SetChaplyginParams
    end type TChaplygin

    contains

    subroutine TChaplygin_Validate(this, OK)
    class(TChaplygin), intent(in) :: this
    logical, intent(inout) :: OK

    call this%TDarkEnergyFluid%Validate(OK)
    if (.not. ieee_is_finite(this%As) .or. .not. ieee_is_finite(this%alpha) .or. &
        this%As <= 0._dl .or. this%As >= 1._dl .or. this%alpha < 0._dl .or. &
        this%alpha*this%As > 1._dl) then
        OK = .false.
        write(*,*) 'Chaplygin requires finite 0<As<1, alpha>=0, and alpha*As<=1.'
    end if

    end subroutine TChaplygin_Validate

    subroutine TChaplygin_SetChaplyginParams(this, Asin, alphain)
    class(TChaplygin), intent(inout) :: this
    real(dl), intent(in) :: Asin, alphain

    this%As = Asin
    this%alpha = alphain

    end subroutine TChaplygin_SetChaplyginParams


    subroutine TChaplygin_ReadParams(this, Ini)
    use IniObjects
    class(TChaplygin) :: this
    class(TIniFile), intent(in) :: Ini

    ! Chaplygin builds its own w(a)/cs2(a) tables in Init; only As, alpha are read.
    this%As    = Ini%Read_Double('chaplygin_As', 0.75d0)
    this%alpha = Ini%Read_Double('chaplygin_alpha', 0.01d0)

    end subroutine TChaplygin_ReadParams


    function TChaplygin_PythonClass()
    character(LEN=:), allocatable :: TChaplygin_PythonClass

    TChaplygin_PythonClass = 'Chaplygin'

    end function TChaplygin_PythonClass


    subroutine TChaplygin_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TChaplygin), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TChaplygin_SelfPointer


    subroutine TChaplygin_Init(this, State)
    use classes
    class(TChaplygin), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl), allocatable :: a(:), w(:), cs2(:)
    real(dl) :: opa, denom, lna_min, dlna
    integer  :: n, i

    if (.not. ieee_is_finite(this%As) .or. .not. ieee_is_finite(this%alpha) .or. &
        this%As <= 0._dl .or. this%As >= 1._dl .or. this%alpha < 0._dl .or. &
        this%alpha*this%As > 1._dl) then
        call GlobalError('Chaplygin requires finite 0<As<1, alpha>=0, and alpha*As<=1.', &
            error_unsupported_params)
        return
    end if

    n = this%n_tab
    allocate(a(n), w(n), cs2(n))

    ! log-a grid, amin = 1e-9, forcing a(n)=1 exactly (SetwTable needs a end=1)
    lna_min = log(1.d-9)
    dlna = -lna_min/dble(n-1)
    do i = 1, n
        a(i) = exp(lna_min + dble(i-1)*dlna)
    end do
    a(n) = 1.d0

    ! Closed-form GCG background: no shooting, evaluate w(a) and cs2(a) directly.
    opa = 1.d0 + this%alpha
    do i = 1, n
        denom  = this%As + (1.d0 - this%As)*a(i)**(-3.d0*opa)
        w(i)   = -this%As/denom
        cs2(i) = -this%alpha*w(i)          ! adiabatic dp/drho, >= 0
    end do

    ! Fill base-class tabulated-w machinery (equation_of_state + logdensity
    ! splines, sets use_tabulated_w=.true., w_lam/wa to today's values).
    call this%SetwTable(a, w, n)

    ! Companion rest-frame cs2(a) spline (same log-a abscissa as w table).
    call this%cs2ofa%Init(log(a), cs2)
    this%cs2_0   = cs2(n)
    this%cs2_lam = cs2(n)   ! today's cs2 (for any halofit/diagnostic use)

    ! Parent fluid Init: checks (no) w-crossing, sets num_perturb_equations.
    call this%TDarkEnergyFluid%Init(State)

    deallocate(a, w, cs2)

    end subroutine TChaplygin_Init


    subroutine TChaplygin_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TChaplygin), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) :: Hv3_over_k, loga, cs2

    ! time-varying rest-frame sound speed from the Chaplygin table
    loga = log(a)
    if (loga <= this%cs2ofa%Xmin_interp) then
        cs2 = this%cs2ofa%F(1)
    else if (loga >= this%cs2ofa%Xmax_interp) then
        cs2 = this%cs2ofa%F(this%cs2ofa%n)
    else
        cs2 = this%cs2ofa%Value(loga)
    end if

    Hv3_over_k = 3*adotoa* y(w_ix + 1) / k
    !density perturbation (Fang-Hu-Lewis rest-frame fluid, cs2 -> cs2(a))
    ayprime(w_ix) = -3 * adotoa * (cs2 - w) * (y(w_ix) + (1 + w) * Hv3_over_k) &
        - (1 + w) * k * y(w_ix + 1) - (1 + w) * k * z
    !account for derivatives of w (use_tabulated_w is always true here)
    if (loga > this%equation_of_state%Xmin_interp .and. &
        loga < this%equation_of_state%Xmax_interp) then
        ayprime(w_ix) = ayprime(w_ix) - adotoa*this%equation_of_state%Derivative(loga)*Hv3_over_k
    end if
    !velocity
    if (abs(w+1) > 1e-6) then
        ayprime(w_ix + 1) = -adotoa * (1 - 3 * cs2) * y(w_ix + 1) + &
            k * cs2 * y(w_ix) / (1 + w)
    else
        ayprime(w_ix + 1) = 0
    end if

    end subroutine TChaplygin_PerturbationEvolve

    end module DarkEnergyChaplygin
