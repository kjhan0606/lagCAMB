    module DMNeutrino
    use DarkMatterInteraction
    use results
    use constants
    use classes
    implicit none

    ! DM-Neutrino scattering model
    ! DM interacts with neutrinos via scattering: sigma ~ sigma_0 * (E_nu/E_0)^n_dmnu
    ! This adds:
    !   - Drag term in CDM velocity equation: -R_dmnu * (v_c - v_nu)
    !   - Collision term in neutrino Boltzmann hierarchy (damping of higher multipoles)
    !
    ! The drag rate on DM from neutrinos:
    !   R_DM = (4/3) * rho_nu / rho_c * n_nu * sigma * c / a
    ! And the back-reaction on neutrinos:
    !   R_nu = rho_c / rho_nu * R_DM * (rho_nu/rho_c) adjustment
    !
    ! For massive neutrinos, the interaction rate depends on the neutrino energy
    ! and the cross section energy dependence.
    !
    ! References:
    !   Mangano+ 2006, Serra+ 2010, Wilkinson+ 2014,
    !   Mosbech & Hannestad 2021, He+ 2023 (Nature Astronomy 2025)
    !
    ! Parameters:
    !   sigma_dmnu [cm^2]: cross section at reference energy E_0 = 1 MeV
    !   n_dmnu: energy power-law index (typically 0 or 2)
    !   m_dm [GeV]: DM particle mass

    type, extends(TDarkMatterModel) :: TDMNeutrinoScattering
        real(dl) :: sigma_dmnu = 0._dl    ! DM-nu cross section [cm^2] at E_0=1MeV
        integer  :: n_dmnu = 0            ! energy dependence power-law index
        real(dl) :: m_dm = 100._dl        ! DM particle mass [GeV]
        ! Cached
        real(dl) :: sigma_0_conf = 0._dl  ! sigma_dmnu in conformal units
        real(dl) :: grhornomass_stored = 0._dl  ! cached massless neutrino density
        real(dl) :: grhormass_total = 0._dl     ! cached massive neutrino density
    contains
    procedure :: ReadParams => TDMNeutrinoScattering_ReadParams
    procedure, nopass :: PythonClass => TDMNeutrinoScattering_PythonClass
    procedure, nopass :: SelfPointer => TDMNeutrinoScattering_SelfPointer
    procedure :: Init => TDMNeutrinoScattering_Init
    procedure :: PerturbationEvolve => TDMNeutrinoScattering_PerturbationEvolve
    procedure :: PerturbationInitial => TDMNeutrinoScattering_PerturbationInitial
    procedure :: PrintFeedback => TDMNeutrinoScattering_PrintFeedback
    procedure :: DragRate_DM_nu => TDMNeutrinoScattering_DragRate
    procedure :: NuCollisionRate => TDMNeutrinoScattering_NuCollisionRate
    end type TDMNeutrinoScattering

    contains

    subroutine TDMNeutrinoScattering_ReadParams(this, Ini)
    use IniObjects
    class(TDMNeutrinoScattering) :: this
    class(TIniFile), intent(in) :: Ini

    this%sigma_dmnu = Ini%Read_Double('sigma_dmnu', 0._dl)
    this%n_dmnu = Ini%Read_Int('n_dmnu', 0)
    this%m_dm = Ini%Read_Double('m_dm_nu', 100._dl)

    end subroutine TDMNeutrinoScattering_ReadParams

    function TDMNeutrinoScattering_PythonClass()
    character(LEN=:), allocatable :: TDMNeutrinoScattering_PythonClass
    TDMNeutrinoScattering_PythonClass = 'DMNeutrinoScattering'
    end function TDMNeutrinoScattering_PythonClass

    subroutine TDMNeutrinoScattering_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TDMNeutrinoScattering), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDMNeutrinoScattering_SelfPointer

    subroutine TDMNeutrinoScattering_Init(this, State)
    class(TDMNeutrinoScattering), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    this%is_standard_cdm = (this%sigma_dmnu == 0._dl)

    if (.not. this%is_standard_cdm) then
        ! CDM gets a velocity perturbation from neutrino drag
        this%has_cdm_velocity = .true.
        this%num_perturb_equations = 1  ! CDM velocity v_c
        this%num_dr_equations = 0       ! no separate DR hierarchy

        this%sigma_0_conf = this%sigma_dmnu

        ! Cache neutrino densities from State
        select type(S => State)
        class is (CAMBdata)
            this%grhornomass_stored = S%grhornomass
            this%grhormass_total = sum(S%grhormass(1:S%CP%Nu_mass_eigenstates))
        end select
    else
        this%has_cdm_velocity = .false.
        this%num_perturb_equations = 0
    end if

    end subroutine TDMNeutrinoScattering_Init

    function TDMNeutrinoScattering_DragRate(this, a, adotoa, grhoc_t, grhonu_t) result(R_drag)
    ! Drag rate on DM from neutrino scattering [Mpc^-1]
    ! R_DM = (4/3) * (rho_nu/rho_c) * n_nu * <sigma v> / a
    ! For relativistic neutrinos: <sigma v> ~ sigma_0 * c * (T_nu/E_0)^n
    ! n_nu ~ (3*zeta(3)/(4*pi^2)) * T_nu^3 per species
    class(TDMNeutrinoScattering), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t, grhonu_t
    real(dl) :: R_drag
    real(dl) :: T_nu, E_ref, rho_ratio, sigma_eff
    real(dl), parameter :: Mpc_in_cm = 3.0856776e24_dl
    real(dl), parameter :: kB_eV = 8.617333e-5_dl  ! Boltzmann constant [eV/K]
    real(dl), parameter :: E_0_eV = 1.0e6_dl       ! reference energy 1 MeV in eV
    real(dl), parameter :: T_nu0_K = 1.9517_dl     ! neutrino temperature today [K]

    if (this%sigma_dmnu == 0._dl .or. grhoc_t == 0._dl) then
        R_drag = 0._dl
        return
    end if

    ! Neutrino temperature at scale factor a
    T_nu = T_nu0_K / a  ! [K]

    ! Energy ratio (T_nu in eV / E_0 in eV)^n
    sigma_eff = this%sigma_dmnu  ! [cm^2]
    if (this%n_dmnu /= 0) then
        E_ref = kB_eV * T_nu / E_0_eV  ! T_nu/1MeV (dimensionless)
        sigma_eff = sigma_eff * E_ref**this%n_dmnu
    end if

    ! Density ratio (using CAMB grho quantities which are 8*pi*G*rho*a^2)
    rho_ratio = grhonu_t / max(grhoc_t, 1e-30_dl)

    ! R_DM [Mpc^-1] ~ (4/3) * rho_nu/rho_c * sigma_eff * n_nu * Mpc_cm / a^2
    ! n_nu is proportional to T_nu^3 ~ a^{-3}
    ! For the neutrino number density: n_nu,0 ~ 112 per species [cm^-3]
    ! n_nu(a) = n_nu,0 / a^3 [cm^-3]
    ! R_DM = (4/3) * rho_nu/rho_c * sigma_eff [cm^2] * n_nu(a) * c * Mpc [cm] / a
    ! But in CAMB units we can simplify using energy densities directly:
    ! R_DM ~ sigma_eff * (rho_nu/(m_dm*c^2)) * c / a (with appropriate unit conversions)

    ! Number density of DM: n_DM = rho_c / (m_dm * c^2)
    ! Interaction rate per DM particle: Gamma_int = n_nu * sigma * c
    ! n_nu(a) ~ 112/cm^3 / a^3 (per species, 3 species active)
    ! Total: n_nu_total ~ 336/cm^3 / a^3
    R_drag = 4._dl/3._dl * rho_ratio * sigma_eff * 336._dl / a**3 * Mpc_in_cm / a

    ! Cap for numerical stability
    R_drag = min(R_drag, 1e4_dl * adotoa)

    end function TDMNeutrinoScattering_DragRate

    function TDMNeutrinoScattering_NuCollisionRate(this, a, adotoa, grhoc_t, grhonu_t) result(R_nu)
    ! Collision rate for neutrino hierarchy damping [Mpc^-1]
    ! Back-reaction: R_nu = (rho_c / rho_nu) * R_DM * (3/4)
    ! This damps neutrino multipoles l >= 2
    class(TDMNeutrinoScattering), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t, grhonu_t
    real(dl) :: R_nu, R_dm

    R_dm = this%DragRate_DM_nu(a, adotoa, grhoc_t, grhonu_t)

    if (grhonu_t > 0._dl) then
        ! Momentum conservation: rho_c * R_DM * (v_c - v_nu) = rho_nu * R_nu * (v_nu - v_c)
        ! => R_nu = (rho_c / rho_nu) * R_DM
        R_nu = R_dm * grhoc_t / grhonu_t
    else
        R_nu = 0._dl
    end if

    ! Cap for stability
    R_nu = min(R_nu, 1e4_dl * adotoa)

    end function TDMNeutrinoScattering_NuCollisionRate

    subroutine TDMNeutrinoScattering_PerturbationInitial(this, y, a, tau, k, dm_ix, dr_ix, vc_ix)
    class(TDMNeutrinoScattering), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer, intent(in) :: dm_ix, dr_ix, vc_ix

    ! CDM velocity starts at zero (adiabatic IC)
    if (this%has_cdm_velocity .and. vc_ix > 0) then
        y(vc_ix) = 0._dl
    end if

    end subroutine TDMNeutrinoScattering_PerturbationInitial

    subroutine TDMNeutrinoScattering_PerturbationEvolve(this, ayprime, a, adotoa, k, z, y, &
        dm_ix, dr_ix, vc_ix, clxc, vb, grhoc_t, grhob_t, sigma, high_ktau_dr)
    ! Evolve CDM velocity with neutrino drag
    ! The neutrino hierarchy modification is done in equations.f90 via NuCollisionRate
    class(TDMNeutrinoScattering), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, y(:)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl), intent(in) :: clxc, vb, grhoc_t, grhob_t, sigma
    logical, intent(in), optional :: high_ktau_dr
    real(dl) :: vc, R_dm, qr_approx, grhonu_t

    if (.not. this%has_cdm_velocity .or. vc_ix <= 0) return

    vc = y(vc_ix)

    ! Neutrino energy density using cached values from Init
    grhonu_t = this%grhornomass_stored / a**2 + this%grhormass_total / a**2

    R_dm = this%DragRate_DM_nu(a, adotoa, grhoc_t, grhonu_t)

    ! CDM velocity equation: v_c' = -aH*v_c - R_DM*(v_c - v_nu)
    ! v_nu is the bulk neutrino velocity qr/4 (for massless) or summed over q bins
    ! We approximate v_nu ~ 0 at leading order (valid for k*tau >> 1)
    ! The full implementation would access qr from the neutrino hierarchy
    ! For now, the drag term drives v_c -> 0 which is the main effect
    ayprime(vc_ix) = -adotoa * vc - R_dm * vc

    end subroutine TDMNeutrinoScattering_PerturbationEvolve

    subroutine TDMNeutrinoScattering_PrintFeedback(this, FeedbackLevel)
    class(TDMNeutrinoScattering) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("DM-Neutrino scattering: sigma_0 = ",ES12.4," cm^2, n = ",I3)') &
            this%sigma_dmnu, this%n_dmnu
        write(*,'("                        m_DM = ",F8.2," GeV")') this%m_dm
    end if

    end subroutine TDMNeutrinoScattering_PrintFeedback

    end module DMNeutrino
