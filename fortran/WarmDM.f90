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
    !    Applied as post-processing to the matter power spectrum.
    !
    ! 2) Full Boltzmann hierarchy (boltzmann_mode=.true.):
    !    WDM registered as an extra massive neutrino species with T_WDM and m_WDM.
    !    Uses existing massive neutrino machinery (momentum integration, hierarchy).
    !    Exact for T_wdm_ratio=1 (standard thermal relic); approximate for T_ratio/=1
    !    because the q-sampling weights assume standard Fermi-Dirac shape.
    !
    ! References:
    !   Bode+ 2001, Viel+ 2005, Viel+ 2013, Lesgourgues & Tram 2011

    type, extends(TDarkMatterModel) :: TWarmDM
        real(dl) :: m_wdm = 0._dl        ! WDM mass in keV
        real(dl) :: T_wdm_ratio = 1._dl  ! T_WDM / T_nu temperature ratio
        real(dl) :: Omega_wdm_h2 = 0._dl ! WDM density parameter (0 = use omch2)
        logical  :: boltzmann_mode = .false. ! True: full Boltzmann; False: transfer function
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
    procedure :: ModifyNeutrinoParams => TWarmDM_ModifyNeutrinoParams
    end type TWarmDM

    contains

    subroutine TWarmDM_ReadParams(this, Ini)
    use IniObjects
    class(TWarmDM) :: this
    class(TIniFile), intent(in) :: Ini

    this%m_wdm = Ini%Read_Double('m_wdm', 0._dl)
    this%T_wdm_ratio = Ini%Read_Double('T_wdm_ratio', 1._dl)
    this%Omega_wdm_h2 = Ini%Read_Double('Omega_wdm_h2', 0._dl)
    this%boltzmann_mode = Ini%Read_Logical('boltzmann_mode', .false.)

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

    subroutine TWarmDM_ModifyNeutrinoParams(this, max_nu_val, &
        Nu_mass_eigenstates, Nu_mass_degeneracies, Nu_mass_fractions, &
        Nu_mass_numbers, Num_Nu_Massive, Num_Nu_Massless, omnuh2, omch2, H0)
    !Register WDM as an extra massive neutrino eigenstate (Boltzmann mode).
    !Called from CAMBdata_SetParams before neutrino validation.
    !
    !Key physics: for a thermal relic WDM with mass m_wdm and density Omega_wdm h^2,
    !the effective temperature ratio is T_ratio = (Omega_wdm h^2 * 94.07 / m_wdm_eV)^(1/3).
    !This ensures find_nu_mass_for_rho gives the PHYSICAL mass (not an effective lighter mass),
    !so that the free-streaming scale is correct.
    class(TWarmDM), intent(inout) :: this
    integer, intent(in) :: max_nu_val
    integer, intent(inout) :: Nu_mass_eigenstates, Num_Nu_Massive
    real(dl), intent(inout) :: Nu_mass_degeneracies(:), Nu_mass_fractions(:)
    integer, intent(inout) :: Nu_mass_numbers(:)
    real(dl), intent(inout) :: Num_Nu_Massless, omnuh2, omch2, H0
    real(dl) :: omega_wdm_h2, old_omnuh2, T_ratio_eff, m_wdm_eV
    integer :: i
    ! neutrino_mass_fac: Omega_nu h^2 = sum(m_nu) / 94.07 eV (per species with deg=1)
    real(dl), parameter :: neutrino_mass_density_fac = 94.07_dl

    if (.not. this%boltzmann_mode .or. this%m_wdm <= 0._dl) return

    ! Determine WDM density
    if (this%Omega_wdm_h2 > 0._dl) then
        omega_wdm_h2 = this%Omega_wdm_h2
    else
        omega_wdm_h2 = omch2  ! all CDM is WDM
    end if

    if (omega_wdm_h2 <= 0._dl) return
    if (omega_wdm_h2 > omch2) then
        write(*,*) 'WarmDM WARNING: omega_wdm_h2 > omch2, capping to omch2'
        omega_wdm_h2 = omch2
    end if

    ! Check array bounds
    if (Nu_mass_eigenstates + 1 > max_nu_val) then
        call GlobalError('WarmDM: max_Nu exceeded when adding WDM species', error_unsupported_params)
        return
    end if

    m_wdm_eV = this%m_wdm * 1e3_dl  ! keV -> eV

    ! Compute the effective temperature ratio from mass and density.
    ! For a Fermi-Dirac thermal relic: Omega h^2 = T_ratio^3 * m / 94.07 eV
    ! This ensures find_nu_mass_for_rho gives nu_mass = m_phys/(k_B*T_wdm),
    ! so the free-streaming scale matches the physical WDM mass.
    T_ratio_eff = (omega_wdm_h2 * neutrino_mass_density_fac / m_wdm_eV) ** (1._dl/3._dl)

    ! Save old omnuh2 for fraction renormalization
    old_omnuh2 = omnuh2

    ! Move density from CDM to massive neutrino sector
    omch2 = omch2 - omega_wdm_h2
    omnuh2 = omnuh2 + omega_wdm_h2

    ! Add new eigenstate
    Nu_mass_eigenstates = Nu_mass_eigenstates + 1

    ! Effective degeneracy = T_ratio^4 (determines radiation density contribution)
    Nu_mass_degeneracies(Nu_mass_eigenstates) = T_ratio_eff**4

    ! Physical number of species
    Nu_mass_numbers(Nu_mass_eigenstates) = 1
    Num_Nu_Massive = Num_Nu_Massive + 1

    ! Mass fractions must sum to 1 over all eigenstates
    if (old_omnuh2 > 0._dl) then
        ! Renormalize existing fractions
        do i = 1, Nu_mass_eigenstates - 1
            Nu_mass_fractions(i) = Nu_mass_fractions(i) * old_omnuh2 / omnuh2
        end do
        Nu_mass_fractions(Nu_mass_eigenstates) = omega_wdm_h2 / omnuh2
    else
        ! Only WDM contributes to omnuh2
        Nu_mass_fractions(Nu_mass_eigenstates) = 1._dl
    end if

    ! In Boltzmann mode, WDM perturbations are handled by neutrino machinery
    this%is_standard_cdm = .true.

    end subroutine TWarmDM_ModifyNeutrinoParams

    subroutine TWarmDM_Init(this, State)
    class(TWarmDM), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) :: h, Omega_wdm, g_wdm

    if (this%boltzmann_mode) then
        ! In Boltzmann mode, WDM is handled as a neutrino species
        ! is_standard_cdm was already set in ModifyNeutrinoParams
        this%has_cdm_velocity = .false.
        this%num_perturb_equations = 0
        this%num_dr_equations = 0
        return
    end if

    this%is_standard_cdm = (this%m_wdm <= 0._dl)

    if (.not. this%is_standard_cdm) then
        this%has_cdm_velocity = .false.
        this%num_perturb_equations = 0
        this%num_dr_equations = 0

        ! Compute the free-streaming scale alpha (Viel+ 2005 Eq. 7)
        select type(S => State)
        class is (CAMBdata)
            h = S%CP%H0 / 100._dl
            if (this%Omega_wdm_h2 > 0._dl) then
                Omega_wdm = this%Omega_wdm_h2 / h**2
            else
                Omega_wdm = S%CP%omch2 / h**2
            end if

            g_wdm = this%T_wdm_ratio**3

            this%alpha_wdm = 0.049_dl * (this%m_wdm)**(-1.11_dl) * &
                (Omega_wdm / 0.25_dl)**(0.11_dl) * (h / 0.7_dl)**(1.22_dl) * &
                g_wdm**(0.15_dl)

            this%nu_wdm = 1.12_dl
        end select
    end if

    end subroutine TWarmDM_Init

    function TWarmDM_TransferFunction(this, k_h) result(Tk)
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
        if (this%boltzmann_mode) then
            write(*,'("         Mode: Full Boltzmann (massive neutrino machinery)")')
        else
            write(*,'("         Mode: Transfer function")')
            if (this%alpha_wdm > 0) then
                write(*,'("         alpha = ",ES10.3," Mpc/h (free-streaming scale)")') this%alpha_wdm
                write(*,'("         k_half = ",F8.3," h/Mpc")') &
                    (1._dl / this%alpha_wdm) * (2._dl**(this%nu_wdm/5._dl) - 1._dl)**(0.5_dl/this%nu_wdm)
            end if
        end if
    end if

    end subroutine TWarmDM_PrintFeedback

    end module WarmDM
