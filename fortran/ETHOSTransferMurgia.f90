    module ETHOSTransferMurgia
    use DarkMatterInteraction
    use results
    use classes
    implicit none

    ! Murgia et al. 2017 (arXiv:1704.07838) phenomenological transfer
    ! function used to describe ETHOS-like small-scale suppression:
    !
    !   T(k) = [1 + (alpha * k)^beta]^gamma
    !
    ! with k in h/Mpc, alpha in Mpc/h. The matter power spectrum is
    ! multiplied by T(k)^2. This is a pure post-processing model: no
    ! background or perturbation modification, no CMB effect. For self
    ! consistent CMB + matter behavior including dark acoustic
    ! oscillations use the full DMDR_ETHOS model instead.

    type, extends(TDarkMatterModel) :: TETHOSTransferMurgia
        real(dl) :: alpha_mpch = 0._dl    ! free-streaming scale [Mpc/h]
        real(dl) :: beta_mur  = 2.24_dl   ! shape parameter
        real(dl) :: gamma_mur = -4.46_dl  ! cutoff sharpness
    contains
    procedure, nopass :: PythonClass => TETHOSTransferMurgia_PythonClass
    procedure, nopass :: SelfPointer => TETHOSTransferMurgia_SelfPointer
    procedure :: Init => TETHOSTransferMurgia_Init
    procedure :: PrintFeedback => TETHOSTransferMurgia_PrintFeedback
    procedure :: TransferFunction => TETHOSTransferMurgia_TransferFunction
    end type TETHOSTransferMurgia

    contains

    function TETHOSTransferMurgia_PythonClass()
    character(LEN=:), allocatable :: TETHOSTransferMurgia_PythonClass
    TETHOSTransferMurgia_PythonClass = 'ETHOSTransferMurgia'
    end function TETHOSTransferMurgia_PythonClass

    subroutine TETHOSTransferMurgia_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TETHOSTransferMurgia), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TETHOSTransferMurgia_SelfPointer

    subroutine TETHOSTransferMurgia_Init(this, State)
    class(TETHOSTransferMurgia), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    this%is_standard_cdm = (this%alpha_mpch <= 0._dl)
    this%has_cdm_velocity = .false.
    this%num_perturb_equations = 0
    this%num_dr_equations = 0
    end subroutine TETHOSTransferMurgia_Init

    function TETHOSTransferMurgia_TransferFunction(this, k_h) result(Tk)
    class(TETHOSTransferMurgia), intent(in) :: this
    real(dl), intent(in) :: k_h
    real(dl) :: Tk, x

    if (this%is_standard_cdm .or. this%alpha_mpch <= 0._dl) then
        Tk = 1._dl
        return
    end if

    x = this%alpha_mpch * k_h
    Tk = (1._dl + x**this%beta_mur)**this%gamma_mur
    end function TETHOSTransferMurgia_TransferFunction

    subroutine TETHOSTransferMurgia_PrintFeedback(this, FeedbackLevel)
    class(TETHOSTransferMurgia) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("ETHOS transfer (Murgia 2017): alpha = ",ES10.3," Mpc/h")') &
            this%alpha_mpch
        write(*,'("                              beta  = ",F8.3)') this%beta_mur
        write(*,'("                              gamma = ",F8.3)') this%gamma_mur
    end if
    end subroutine TETHOSTransferMurgia_PrintFeedback

    end module ETHOSTransferMurgia
