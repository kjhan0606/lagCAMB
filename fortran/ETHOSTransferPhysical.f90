    module ETHOSTransferPhysical
    use DarkMatterInteraction
    use results
    use classes
    implicit none

    ! Physics-parameterized surrogate for ETHOS DM-DR small-scale suppression.
    !
    ! Maps physical ETHOS parameters (a_dark_n, xi_dr, n_dark) to the Murgia
    ! 2017 transfer function T(k) = [1 + (alpha k)^beta]^gamma via a
    ! power-law fit:
    !
    !   alpha [Mpc/h] = A0 * (a_dark_n / a_ref)^p_a * (xi_dr / xi_ref)^p_xi
    !                    * (omega_dmdr_h2 / 0.12)^p_om
    !
    ! The fit constants (A0, p_a, p_xi, p_om) are user-exposed so each science
    ! application can recalibrate against the full DMDR_ETHOS Boltzmann
    ! solution for its parameter ranges. Defaults give an order-of-magnitude
    ! correct mapping for ETHOS n_dark=2 (dark Thomson-like opacity).
    !
    ! This is a fast surrogate intended for parameter scans; it does NOT
    ! reproduce dark acoustic oscillations or CMB effects. Use DMDR_ETHOS for
    ! the full self-consistent calculation.

    type, extends(TDarkMatterModel) :: TETHOSTransferPhysical
        ! Physical inputs
        real(dl) :: a_dark_n    = 0._dl    ! ETHOS opacity coefficient [Mpc^-1]
        real(dl) :: xi_dr       = 0.5_dl   ! T_DR / T_gamma at recomb
        real(dl) :: omdmdrh2    = 0.12_dl  ! interacting DM density Omega_dmdr h^2
        integer  :: n_dark      = 2        ! opacity index (advisory only)
        ! Surrogate fit constants (override to recalibrate)
        real(dl) :: A0_fit      = 0.5_dl   ! alpha [Mpc/h] at reference point
        real(dl) :: a_ref       = 1.e-4_dl ! reference a_dark_n [Mpc^-1]
        real(dl) :: xi_ref      = 0.5_dl   ! reference xi_dr
        real(dl) :: p_a         = 0.25_dl  ! power on a_dark_n
        real(dl) :: p_xi        = 1._dl    ! power on xi_dr
        real(dl) :: p_om        = -0.1_dl  ! power on omega_dmdr_h2
        ! Murgia shape parameters (defaults from Murgia+2017)
        real(dl) :: beta_mur    = 2.24_dl
        real(dl) :: gamma_mur   = -4.46_dl
        ! Derived
        real(dl) :: alpha_mpch  = 0._dl
    contains
    procedure, nopass :: PythonClass => TETHOSTransferPhysical_PythonClass
    procedure, nopass :: SelfPointer => TETHOSTransferPhysical_SelfPointer
    procedure :: Init => TETHOSTransferPhysical_Init
    procedure :: PrintFeedback => TETHOSTransferPhysical_PrintFeedback
    procedure :: TransferFunction => TETHOSTransferPhysical_TransferFunction
    end type TETHOSTransferPhysical

    contains

    function TETHOSTransferPhysical_PythonClass()
    character(LEN=:), allocatable :: TETHOSTransferPhysical_PythonClass
    TETHOSTransferPhysical_PythonClass = 'ETHOSTransferPhysical'
    end function TETHOSTransferPhysical_PythonClass

    subroutine TETHOSTransferPhysical_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TETHOSTransferPhysical), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TETHOSTransferPhysical_SelfPointer

    subroutine TETHOSTransferPhysical_Init(this, State)
    !alpha_mpch is computed in Python set_params(); Fortran Init only sets
    !the flag fields based on the already-stored alpha.
    class(TETHOSTransferPhysical), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    this%is_standard_cdm = (this%alpha_mpch <= 0._dl)
    this%has_cdm_velocity = .false.
    this%num_perturb_equations = 0
    this%num_dr_equations = 0
    end subroutine TETHOSTransferPhysical_Init

    function TETHOSTransferPhysical_TransferFunction(this, k_h) result(Tk)
    class(TETHOSTransferPhysical), intent(in) :: this
    real(dl), intent(in) :: k_h
    real(dl) :: Tk, x

    if (this%is_standard_cdm .or. this%alpha_mpch <= 0._dl) then
        Tk = 1._dl
        return
    end if

    x = this%alpha_mpch * k_h
    Tk = (1._dl + x**this%beta_mur)**this%gamma_mur
    end function TETHOSTransferPhysical_TransferFunction

    subroutine TETHOSTransferPhysical_PrintFeedback(this, FeedbackLevel)
    class(TETHOSTransferPhysical) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("ETHOS transfer (physical surrogate):")')
        write(*,'("   a_dark_n  = ",ES10.3," Mpc^-1   (n_dark=",I2,")")') &
            this%a_dark_n, this%n_dark
        write(*,'("   xi_dr     = ",F8.4)') this%xi_dr
        write(*,'("   -> alpha  = ",ES10.3," Mpc/h")') this%alpha_mpch
        write(*,'("      beta   = ",F8.3,"   gamma = ",F8.3)') &
            this%beta_mur, this%gamma_mur
    end if
    end subroutine TETHOSTransferPhysical_PrintFeedback

    end module ETHOSTransferPhysical
