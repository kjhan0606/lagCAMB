    module WarmDM
    use DarkMatterInteraction
    use results
    use constants
    use classes
    implicit none

    ! Warm Dark Matter (WDM) model
    ! Thermal relic WDM with free-streaming that suppresses small-scale structure.
    !
    ! Two approaches:
    ! 1) Transfer function fitting formula (Viel+ 2005, Bode+ 2001):
    !    T(k) = [1 + (alpha*k)^(2*nu)]^(-5/nu)
    !    where alpha depends on m_WDM and Omega_WDM
    !    This is applied as a post-processing step to the matter power spectrum.
    !
    ! 2) Full Boltzmann hierarchy (like massive neutrinos with different mass/temperature):
    !    WDM as an extra massive species with T_WDM and m_WDM
    !    This approach uses the existing massive neutrino machinery.
    !
    ! We implement approach 1 (transfer function) as a DM model that modifies P(k),
    ! plus store parameters needed for approach 2 (Boltzmann) if user wants full treatment.
    !
    ! References:
    !   Bode+ 2001, Viel+ 2005, Viel+ 2013 (Ly-alpha constraints)
    !
    ! Parameters:
    !   m_wdm [keV]: WDM particle mass
    !   T_wdm_ratio: temperature ratio T_WDM/T_nu (default 1 for standard thermal relic)
    !   Omega_wdm_h2: WDM density (defaults to total CDM if 0)

    type, extends(TDarkMatterModel) :: TWarmDM
        real(dl) :: m_wdm = 0._dl        ! WDM mass in keV
        real(dl) :: T_wdm_ratio = 1._dl  ! T_WDM / T_nu temperature ratio
        real(dl) :: Omega_wdm_h2 = 0._dl ! WDM density parameter (0 = use omch2)
        ! Cached transfer function parameters
        real(dl) :: alpha_wdm = 0._dl    ! free-streaming scale [Mpc/h]
        real(dl) :: nu_wdm = 1.12_dl     ! shape parameter
    contains
    procedure :: ReadParams => TWarmDM_ReadParams
    procedure, nopass :: PythonClass => TWarmDM_PythonClass
    procedure, nopass :: SelfPointer => TWarmDM_SelfPointer
    procedure :: Init => TWarmDM_Init
    procedure :: PrintFeedback => TWarmDM_PrintFeedback
    procedure :: TransferFunction => TWarmDM_TransferFunction
    end type TWarmDM

    contains

    subroutine TWarmDM_ReadParams(this, Ini)
    use IniObjects
    class(TWarmDM) :: this
    class(TIniFile), intent(in) :: Ini

    this%m_wdm = Ini%Read_Double('m_wdm', 0._dl)
    this%T_wdm_ratio = Ini%Read_Double('T_wdm_ratio', 1._dl)
    this%Omega_wdm_h2 = Ini%Read_Double('Omega_wdm_h2', 0._dl)

    end subroutine TWarmDM_ReadParams

    function TWarmDM_PythonClass()
    character(LEN=:), allocatable :: TWarmDM_PythonClass
    TWarmDM_PythonClass = 'WarmDM'
    end function TWarmDM_PythonClass

    subroutine TWarmDM_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TWarmDM), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TWarmDM_SelfPointer

    subroutine TWarmDM_Init(this, State)
    class(TWarmDM), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) :: h, Omega_wdm, g_wdm

    this%is_standard_cdm = (this%m_wdm <= 0._dl)

    if (.not. this%is_standard_cdm) then
        ! No extra perturbation equations needed for transfer function approach
        this%has_cdm_velocity = .false.
        this%num_perturb_equations = 0
        this%num_dr_equations = 0

        ! Compute the free-streaming scale alpha (Viel+ 2005 Eq. 7)
        ! alpha = 0.049 * (m_WDM/keV)^(-1.11) * (Omega_WDM/0.25)^(0.11) * (h/0.7)^(1.22) [Mpc/h]
        select type(S => State)
        class is (CAMBdata)
            h = S%CP%H0 / 100._dl
            if (this%Omega_wdm_h2 > 0._dl) then
                Omega_wdm = this%Omega_wdm_h2 / h**2
            else
                Omega_wdm = S%CP%omch2 / h**2
            end if

            ! Degrees of freedom correction for non-standard T_WDM
            ! g_wdm effectively shifts the mass constraint
            g_wdm = this%T_wdm_ratio**3  ! ~ (T_wdm/T_nu)^3

            this%alpha_wdm = 0.049_dl * (this%m_wdm)**(-1.11_dl) * &
                (Omega_wdm / 0.25_dl)**(0.11_dl) * (h / 0.7_dl)**(1.22_dl) * &
                g_wdm**(0.15_dl)  ! g_wdm correction from Viel+ 2005

            this%nu_wdm = 1.12_dl  ! Viel+ 2005 shape parameter
        end select
    end if

    end subroutine TWarmDM_Init

    function TWarmDM_TransferFunction(this, k_h) result(Tk)
    ! WDM transfer function T(k) relative to CDM
    ! T(k) = [1 + (alpha*k)^(2*nu)]^(-5/nu)
    ! k_h is in units of h/Mpc
    class(TWarmDM), intent(in) :: this
    real(dl), intent(in) :: k_h
    real(dl) :: Tk
    real(dl) :: x

    if (this%is_standard_cdm .or. this%alpha_wdm <= 0._dl) then
        Tk = 1._dl
        return
    end if

    x = this%alpha_wdm * k_h
    Tk = (1._dl + x**(2._dl * this%nu_wdm))**(-5._dl / this%nu_wdm)

    end function TWarmDM_TransferFunction

    subroutine TWarmDM_PrintFeedback(this, FeedbackLevel)
    class(TWarmDM) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("Warm DM: m_WDM = ",F10.4," keV")') this%m_wdm
        write(*,'("         T_WDM/T_nu = ",F8.4)') this%T_wdm_ratio
        if (this%alpha_wdm > 0) then
            write(*,'("         alpha = ",ES10.3," Mpc/h (free-streaming scale)")') this%alpha_wdm
            write(*,'("         k_half = ",F8.3," h/Mpc")') &
                (1._dl / this%alpha_wdm) * (2._dl**(this%nu_wdm/5._dl) - 1._dl)**(0.5_dl/this%nu_wdm)
        end if
    end if

    end subroutine TWarmDM_PrintFeedback

    end module WarmDM
