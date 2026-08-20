    module DarkEnergyRunningVacuum
    use DarkEnergyInterface
    use results
    use constants
    use classes
    use config, only: GlobalError, error_unsupported_params
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    ! Running Vacuum Model (RVM), Sola group.
    ! Vacuum energy density scales with the Hubble rate:
    !     Lambda(H^2) = c0 + nu*H^2,
    !     rho_L(H) = (3/(8 pi G)) (c0 + nu H^2),
    ! and the vacuum exchanges energy with COLD DARK MATTER ONLY (baryons and
    ! radiation stay uncoupled, as observationally required). Local conservation
    !     rho_c' + 3H rho_c = -rho_L'
    ! together with rho_L(H) closes into the exact background solution
    !     rho_c(a) = rho_c0 * a^{-3(1-nu)}
    !     rho_L(a) = rho_L0 + [nu/(1-nu)] * rho_c0 * (a^{-3(1-nu)} - 1)
    ! (check: d rho_c/dln a = -3(1-nu) rho_c, d rho_L/dln a = -3 nu rho_c, so
    !  d(rho_c+rho_L)/dln a = -3 rho_c and conservation closes). The equation of
    ! state of the vacuum is w = -1 identically; the dark-energy component is a
    ! smooth vacuum with a time-varying density. At early times
    ! rho_L(a->0) ~ [nu/(1-nu)] rho_c(a), a small constant fraction nu of the CDM.
    !
    ! Perturbations: the vacuum is SMOOTH (delta rho_L = 0 in the CDM rest frame),
    ! the standard Gomez-Valent/Sola/Basilakos linear treatment. The energy
    ! transfer Q = -rho_L' feeds unclustered particles into the CDM, so the
    ! perturbed CDM continuity gains a dilution term (added in equations.f90).
    !
    ! References: Sola 2013 (arXiv:1306.1527); Gomez-Valent & Sola 2015
    !             (arXiv:1409.7048); Basilakos, Sola & Gomez-Valent.

    type, extends(TDarkEnergyModel) :: TRunningVacuum
        real(dl) :: nu = 0._dl        ! Lambda(H^2) = c0 + nu*H^2 running coefficient
        real(dl) :: grhoc_rvm = 0._dl ! cached 8 pi G rho_c0 = State%grhoc (for rho_L(a))
    contains
    procedure :: ReadParams => TRunningVacuum_ReadParams
    procedure :: Validate => TRunningVacuum_Validate
    procedure, nopass :: PythonClass => TRunningVacuum_PythonClass
    procedure, nopass :: SelfPointer => TRunningVacuum_SelfPointer
    procedure :: Init => TRunningVacuum_Init
    procedure :: BackgroundDensityAndPressure => TRunningVacuum_BackgroundDensityAndPressure
    procedure :: CDM_BackgroundCorrection => TRunningVacuum_CDM_BackgroundCorrection
    procedure :: PrintFeedback => TRunningVacuum_PrintFeedback
    procedure :: Effective_w_wa => TRunningVacuum_Effective_w_wa
    end type TRunningVacuum

    contains

    subroutine TRunningVacuum_Validate(this, OK)
    class(TRunningVacuum), intent(in) :: this
    logical, intent(inout) :: OK

    call this%TDarkEnergyModel%Validate(OK)
    if (.not. ieee_is_finite(this%nu) .or. abs(this%nu) > 0.01_dl) then
        OK = .false.
        write(*,*) 'RunningVacuum nu must be finite and satisfy |nu|<=0.01.'
    end if

    end subroutine TRunningVacuum_Validate

    subroutine TRunningVacuum_ReadParams(this, Ini)
    use IniObjects
    class(TRunningVacuum) :: this
    class(TIniFile), intent(in) :: Ini

    this%nu = Ini%Read_Double('rvm_nu', 0._dl)

    end subroutine TRunningVacuum_ReadParams

    function TRunningVacuum_PythonClass()
    character(LEN=:), allocatable :: TRunningVacuum_PythonClass
    TRunningVacuum_PythonClass = 'RunningVacuum'
    end function TRunningVacuum_PythonClass

    subroutine TRunningVacuum_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TRunningVacuum), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TRunningVacuum_SelfPointer

    subroutine TRunningVacuum_Init(this, State)
    use classes
    class(TRunningVacuum), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    if (.not. ieee_is_finite(this%nu) .or. abs(this%nu) > 0.01_dl) then
        call GlobalError('RunningVacuum nu must be finite and satisfy |nu|<=0.01.', error_unsupported_params)
        return
    end if

    ! Cache 8 pi G rho_c0 (= State%grhoc, set before DarkEnergy%Init) so the
    ! smooth-vacuum density rho_L(a) can add its CDM-tracking excess. grhoc is
    ! not passed to BackgroundDensityAndPressure, unlike CDM_BackgroundCorrection.
    select type(S => State)
    class is (CAMBdata)
        this%grhoc_rvm = S%grhoc
    end select

    ! nu = 0 is bit-identical LCDM: treat it as a cosmological constant so the
    ! background and perturbation code take the untouched Lambda paths.
    this%is_cosmological_constant = (this%nu == 0._dl)
    ! Smooth vacuum: no dark-energy perturbation variables (delta rho_L = 0).
    this%num_perturb_equations = 0

    end subroutine TRunningVacuum_Init

    subroutine TRunningVacuum_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    ! grhov_t = 8 pi G rho_L(a) a^2 for the smooth running vacuum (w = -1).
    !   rho_L(a) = rho_L0 + [nu/(1-nu)] rho_c0 (a^{-3(1-nu)} - 1)
    ! In grho units (grhov = 8 pi G rho_L0, grhoc_rvm = 8 pi G rho_c0):
    !   grhov_t = grhov*a^2 + [nu/(1-nu)] grhoc_rvm (a^{-3(1-nu)} - 1) a^2.
    ! nu = 0 reduces exactly to grhov*a^2 (cosmological constant).
    class(TRunningVacuum), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w
    real(dl) :: a2

    if (present(w)) w = -1._dl

    a2 = a * a
    if (this%nu == 0._dl .or. a <= 0._dl) then
        grhov_t = grhov * a2
        return
    end if

    grhov_t = grhov * a2 + (this%nu / (1._dl - this%nu)) * this%grhoc_rvm * &
        (a**(-3._dl*(1._dl - this%nu)) - 1._dl) * a2

    end subroutine TRunningVacuum_BackgroundDensityAndPressure

    function TRunningVacuum_CDM_BackgroundCorrection(this, grhoc, grhov, a) result(dgrhoc_t)
    ! CDM excess over the naive a^{-3} scaling from the vacuum-CDM exchange.
    !   rho_c(a) = rho_c0 a^{-3(1-nu)}, so
    !   dgrhoc_t = 8 pi G a^2 [rho_c(a) - rho_c0 a^{-3}]
    !            = grhoc * (a^{-3(1-nu)} - a^{-3}) * a^2.
    ! Vanishes at a = 1 (rho_c0 fixed to today's CDM) and for nu = 0.
    class(TRunningVacuum) :: this
    real(dl), intent(in) :: grhoc, grhov, a
    real(dl) :: dgrhoc_t

    dgrhoc_t = 0._dl
    if (this%nu == 0._dl .or. a <= 0._dl) return

    dgrhoc_t = grhoc * (a**(-3._dl*(1._dl - this%nu)) - a**(-3._dl)) * a * a

    end function TRunningVacuum_CDM_BackgroundCorrection

    subroutine TRunningVacuum_Effective_w_wa(this, w, wa)
    class(TRunningVacuum), intent(inout) :: this
    real(dl), intent(out) :: w, wa

    ! Smooth vacuum, exactly w = -1 (used only for non-linear/halofit fitting).
    w = -1._dl
    wa = 0._dl

    end subroutine TRunningVacuum_Effective_w_wa

    subroutine TRunningVacuum_PrintFeedback(this, FeedbackLevel)
    class(TRunningVacuum) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) write(*,'("Running Vacuum: Lambda(H^2)=c0+nu H^2, nu = ",ES12.4)') this%nu

    end subroutine TRunningVacuum_PrintFeedback

    end module DarkEnergyRunningVacuum
