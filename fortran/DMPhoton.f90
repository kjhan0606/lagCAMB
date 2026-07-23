    module DMPhoton
    use DarkMatterInteraction
    use results
    use constants
    use classes
    implicit none

    ! DM-Photon scattering model (Boddy & Gluscevic 2018)
    ! sigma_DM-gamma = sigma_Th * u_idm_g * (m_DM / 100 GeV)
    ! Photon opacity from DM: mu_dot = sigma_0_conf / a^(2+n)
    ! DM drag: R_DM = (4/3)(rho_gamma/rho_c) * mu_dot
    type, extends(TDarkMatterModel) :: TDMPhotonScattering
        real(dl) :: u_idm_g = 0._dl        ! dimensionless coupling constant
        integer  :: n_idm_g = 0            ! temperature power index
        real(dl) :: m_dm = 100._dl         ! DM mass [GeV]
        ! Cached quantities (set in Init)
        real(dl) :: sigma_0_conf = 0._dl   ! conformal Mpc^-1 cross section
        real(dl) :: grhog_stored = 0._dl   ! State%grhog
    contains
    procedure :: ReadParams => TDMPhotonScattering_ReadParams
    procedure, nopass :: PythonClass => TDMPhotonScattering_PythonClass
    procedure, nopass :: SelfPointer => TDMPhotonScattering_SelfPointer
    procedure :: Init => TDMPhotonScattering_Init
    procedure :: PerturbationEvolve => TDMPhotonScattering_PerturbationEvolve
    procedure :: PerturbationInitial => TDMPhotonScattering_PerturbationInitial
    procedure :: PrintFeedback => TDMPhotonScattering_PrintFeedback
    procedure :: DragRate_Photon => TDMPhotonScattering_DragRate_Photon
    end type TDMPhotonScattering

    contains

    subroutine TDMPhotonScattering_ReadParams(this, Ini)
    use IniObjects
    class(TDMPhotonScattering) :: this
    class(TIniFile), intent(in) :: Ini

    this%u_idm_g = Ini%Read_Double('u_idm_g', 0._dl)
    this%n_idm_g = Ini%Read_Int('n_idm_g', 0)
    this%m_dm = Ini%Read_Double('m_dm_photon', 100._dl)

    end subroutine TDMPhotonScattering_ReadParams

    function TDMPhotonScattering_PythonClass()
    character(LEN=:), allocatable :: TDMPhotonScattering_PythonClass
    TDMPhotonScattering_PythonClass = 'DMPhotonScattering'
    end function TDMPhotonScattering_PythonClass

    subroutine TDMPhotonScattering_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TDMPhotonScattering), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDMPhotonScattering_SelfPointer

    subroutine TDMPhotonScattering_Init(this, State)
    class(TDMPhotonScattering), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl), parameter :: sigma_Th_cm2 = 6.6524587e-25_dl  ! Thomson cross section [cm^2]
    real(dl), parameter :: Mpc_in_cm = 3.0856776e24_dl

    this%is_standard_cdm = (this%u_idm_g == 0._dl)
    if (.not. this%is_standard_cdm) then
        this%has_cdm_velocity = .true.
        this%num_perturb_equations = 1  ! CDM velocity v_c
        this%num_dr_equations = 0

        ! sigma_0 [cm^2] = sigma_Th * u_idm_g * (m_dm / 100 GeV)
        ! Convert to conformal Mpc^-1 units: sigma_0_conf = sigma_0 * n_gamma,0 * Mpc_cm
        ! n_gamma,0 = 2*zeta(3)/pi^2 * T_CMB^3 = 410.7 cm^-3
        ! But the opacity is mu_dot = sigma_0_conf / a^(2+n)
        ! Actually follow CAMB convention like akthom = sigma_T * Nnow * Mpc:
        ! mu_dot_idm_g = sigma_0 * n_gamma(a) * c / a * Mpc (in conformal Mpc^-1)
        ! n_gamma(a) = n_gamma,0 / a^3
        ! mu_dot = sigma_0_cm2 * n_gamma0 * Mpc_cm / a^3 / a = sigma_0_cm2 * n_gamma0 * Mpc_cm / a^(2+n_idm_g+2)
        ! Wait, let me think more carefully.
        !
        ! The DM-photon opacity contribution to photon hierarchy:
        !   mu_dot = n_DM * sigma * c  (in physical time, same structure as Thomson opacity)
        ! In conformal time (d/dtau = a * d/dt): mu_dot_conf = n_DM * sigma * c * a
        ! But CAMB opacity = n_e * sigma_T * c * a (conformal) = akthom * x_e / a
        ! where akthom = sigma_T * Nnow * Mpc is in Mpc^-1 conformal units
        !
        ! For DM-photon:
        ! n_DM(a) = n_DM,0 / a^3 = (rho_c,0 / m_DM) / a^3
        ! sigma(T) = sigma_0 * (T/T_0)^n = sigma_0 / a^n  (for CMB photon temperature T ~ 1/a)
        ! mu_dot_conf = n_DM,0 * sigma_0 * c / a^3 / a^n * a = n_DM,0 * sigma_0 * c / a^(2+n)
        !
        ! In CAMB units (Mpc^-1):
        ! mu_dot_conf = sigma_0_cm2 * (rho_c,0 / m_DM_g) / a^(2+n) * c_cm_s * Mpc_s
        ! where Mpc_s = Mpc_cm / c_cm_s, so * c * Mpc_s = Mpc_cm
        ! mu_dot_conf = sigma_0_cm2 * n_DM,0_cm3 * Mpc_cm / a^(2+n)

        select type(S => State)
        class is (CAMBdata)
            this%grhog_stored = S%grhog
            block
                real(dl) :: sigma_0_cm2, n_dm_0, m_dm_g
                real(dl), parameter :: GeV_to_g = 1.78266192e-24_dl
                real(dl), parameter :: c_cm_s = 2.99792458e10_dl

                sigma_0_cm2 = sigma_Th_cm2 * this%u_idm_g * (this%m_dm / 100._dl)
                m_dm_g = this%m_dm * GeV_to_g
                ! n_DM,0 = rho_c,0 / m_DM (but we need rho_c,0 in g/cm^3)
                ! From CAMB: grhoc = 8*pi*G * rho_c * (Mpc/c)^2 = 3*H0^2*Omega_c
                ! rho_c_phys = grhoc / (8*pi*G) * c^2 / Mpc^2 [g/cm^3 * cm^2/s^2]
                ! This is complex, so use simpler: mu_dot scales with sigma * akthom_analog
                ! akthom_DM = sigma_0_cm2 * n_DM,0 [cm^-3] * Mpc_cm
                ! But n_DM,0 depends on rho_c,0 which depends on cosmological parameters.
                !
                ! Alternative: express via grhoc = 8*pi*G*rho_c*a^4 / a^3
                ! mu_dot = (4/3)*(grhog_t/grhoc_t) * R_photon
                ! Actually let's use a simpler approach:
                ! In units where the opacity is kappa = sigma * n / a,
                ! we can express the DM-photon scattering opacity as:
                ! mu_dot = sigma_0 * n_DM,0 * Mpc_cm / a^(2+n)
                ! Then grhog_stored is used for DM drag rate calculation

                ! n_DM,0 [cm^-3]: approximate from Omega_c * rho_crit / m_DM
                ! rho_crit = 3 H0^2 / (8 pi G) ~ 1.878e-29 h^2 g/cm^3
                ! Omega_c ~ grhoc / (3 * (H0*Mpc_s/c)^2) ... complex
                ! Simplify: store sigma_0_conf = sigma_0_cm2 * Mpc_cm, then compute
                ! mu_dot at runtime using CAMB density variables
                this%sigma_0_conf = sigma_0_cm2 * Mpc_in_cm
            end block
        end select
    else
        this%has_cdm_velocity = .false.
        this%num_perturb_equations = 0
    end if

    end subroutine TDMPhotonScattering_Init

    function TDMPhotonScattering_DragRate_Photon(this, a, adotoa, grhoc_t, grhog_t) result(mu_dot)
    ! DM-photon opacity [Mpc^-1]: additional opacity for photon Boltzmann hierarchy
    ! mu_dot = n_DM * sigma(T) * c (in conformal Mpc^-1 units)
    ! = sigma_0_conf * (rho_c / m_dm) / a^(2+n)  [via akthom analog]
    !
    ! We use CAMB grho variables:
    ! grhoc_t = 8*pi*G*rho_c*a^2 => rho_c = grhoc_t / (8*pi*G*a^2)
    ! n_DM = rho_c / m_DM, and mu_dot = n_DM * sigma(T) * c
    ! sigma(T) = sigma_0 * (T/T0)^n = sigma_0 / a^n
    ! mu_dot_conf (per Mpc) = n_DM(phys) * sigma / a (extra 1/a from conformal->physical)
    !   Actually in conformal time, the opacity picks up no extra factor:
    !   dtau_opacity/dtau_conformal = n_scatter * sigma * c
    !   (because both are in conformal time)
    !
    ! So: mu_dot = sigma_0_cm2 * n_DM,0 / a^3 / a^n * Mpc_cm
    !           = sigma_0_conf * n_DM,0 / a^(3+n)
    ! But we stored sigma_0_conf = sigma_0_cm2 * Mpc_cm, and need n_DM,0.
    !
    ! Use the fact that grhoc_t = grhoc / a, and grhoc = 8*pi*G * rho_c,0 * Mpc^2 / c^2
    ! n_DM,0 = rho_c,0 / m_DM_GeV / GeV_to_g
    ! This requires extracting rho_c from grhoc, which involves G and unit conversions.
    !
    ! SIMPLER APPROACH: Follow CLASS convention where
    ! mu_dot ~ u_idm_g * (m_dm/100) * sigma_Th * n_DM,0 * Mpc_cm / a^(3+n)
    ! Pre-compute the a-independent part in Init and just divide by a^(3+n) here.
    !
    ! But we don't know n_DM,0 without the cosmology. So we compute it from grhoc_t:
    ! grhoc_t/a = grhoc/a^2 = 8*pi*G * rho_c,0 * Mpc^2 / (c^2 * a^2)
    ! -> rho_c,0 [g/cm^3] = grhoc * c^2 / (8*pi*G * Mpc^2)
    !
    ! Actually, let me use a direct physical formula:
    ! Gamma_DM-gamma [s^-1] = (4/3) * (rho_gamma / rho_DM) * n_DM * sigma * c
    ! But the way it enters the photon Boltzmann:
    ! dF_l/dtau += -mu_dot * F_l (same form as Thomson opacity)
    ! where mu_dot [Mpc^-1] = n_DM * sigma * c_light * (Mpc/c_light)
    !
    ! n_DM * sigma = (rho_c / m_DM) * sigma_0 / a^(3+n)
    ! mu_dot = rho_c / m_DM * sigma_0 / a^(3+n) * Mpc
    !
    ! From CAMB: rho_c = grhoc * c^2 / (8*pi*G * Mpc^2)
    ! mu_dot = grhoc * c^2 / (8*pi*G * Mpc^2) / m_DM * sigma_0 / a^(3+n) * Mpc
    !        = grhoc / (8*pi*G * Mpc) * c^2 / m_DM * sigma_0 / a^(3+n)
    !
    ! This is getting complex. Let me use a simpler parameterization:
    ! Following Boddy & Gluscevic 2018, define:
    ! mu_dot = u_idm_g * akthom_dm / a^(1+n)
    ! where akthom_dm = sigma_Th * (m_dm/100GeV) * Nnow * Mpc
    ! Wait, that doesn't work either since we need DM number density, not baryon.
    !
    ! Let's just compute everything numerically:
    class(TDMPhotonScattering), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t, grhog_t
    real(dl) :: mu_dot
    real(dl) :: n_dm_phys
    real(dl), parameter :: GeV_to_g = 1.78266192e-24_dl
    real(dl), parameter :: Mpc_in_cm = 3.0856776e24_dl
    real(dl), parameter :: c_cm_s = 2.99792458e10_dl
    real(dl), parameter :: G_cgs = 6.67430e-8_dl
    real(dl), parameter :: sigma_Th_cm2 = 6.6524587e-25_dl
    real(dl) :: sigma_cm2, rho_c_cgs, m_dm_g

    if (this%u_idm_g == 0._dl) then
        mu_dot = 0._dl
        return
    end if

    ! DM cross section [cm^2]
    sigma_cm2 = sigma_Th_cm2 * this%u_idm_g * (this%m_dm / 100._dl)

    ! Temperature scaling: sigma(a) = sigma_0 / a^n
    sigma_cm2 = sigma_cm2 / a**this%n_idm_g

    ! Physical DM number density: n_DM = rho_c / m_DM
    ! From CAMB: grhoc_t = 8*pi*G*rho_c*a^2 (in Mpc^-2 units)
    ! grhoc_t is actually in units where c=1, Mpc=1, so
    ! grhoc_t [Mpc^-2] = 8*pi*G/c^2 * rho_c [g/cm^3] * (Mpc_cm)^2 * a^2
    ! => rho_c = grhoc_t * c^2 / (8*pi*G * Mpc^2 * a^2)
    m_dm_g = this%m_dm * GeV_to_g

    ! mu_dot [Mpc^-1] = n_DM [cm^-3] * sigma [cm^2] * Mpc_cm
    ! n_DM = rho_c / m_DM = grhoc_t * c^2 / (8*pi*G * Mpc^2 * a^2) / m_DM
    ! mu_dot = grhoc_t * c^2 / (8*pi*G * Mpc^2 * a^2 * m_DM) * sigma * Mpc_cm
    !        = grhoc_t * c^2 * sigma / (8*pi*G * Mpc * a^2 * m_DM)
    rho_c_cgs = grhoc_t * c_cm_s**2 / (8._dl * const_pi * G_cgs * Mpc_in_cm**2 * a**2)
    n_dm_phys = rho_c_cgs / m_dm_g  ! [cm^-3]

    mu_dot = n_dm_phys * sigma_cm2 * Mpc_in_cm

    ! Cap for numerical stability
    mu_dot = min(mu_dot, 1e4_dl * adotoa)

    end function TDMPhotonScattering_DragRate_Photon

    subroutine TDMPhotonScattering_PerturbationInitial(this, y, a, tau, k, dm_ix, dr_ix, vc_ix)
    class(TDMPhotonScattering), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer, intent(in) :: dm_ix, dr_ix, vc_ix

    if (this%has_cdm_velocity .and. vc_ix > 0) then
        y(vc_ix) = 0._dl
    end if

    end subroutine TDMPhotonScattering_PerturbationInitial

    subroutine TDMPhotonScattering_PerturbationEvolve(this, ayprime, a, adotoa, k, z, y, &
        dm_ix, dr_ix, vc_ix, clxc, vb, grhoc_t, grhob_t, sigma, high_ktau_dr)
    ! CDM velocity: v_c' = -aH*v_c - R_DM_gamma * v_c
    ! The +R_DM_gamma * qg/4 correction is added in equations.f90
    class(TDMPhotonScattering), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, y(:)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl), intent(in) :: clxc, vb, grhoc_t, grhob_t, sigma
    logical, intent(in), optional :: high_ktau_dr
    real(dl) :: vc, R_dm_gamma, mu_dot_ph, grhog_t

    if (.not. this%has_cdm_velocity .or. vc_ix <= 0) return

    vc = y(vc_ix)

    ! DM drag from photons: R_DM = (4/3) * (rho_gamma/rho_c) * mu_dot
    grhog_t = this%grhog_stored / a**2
    mu_dot_ph = this%DragRate_Photon(a, adotoa, grhoc_t, grhog_t)
    R_dm_gamma = 4._dl / 3._dl * grhog_t / max(grhoc_t, 1e-30_dl) * mu_dot_ph

    ! v_c' = -aH*v_c - R_DM*(v_c)  [qg/4 correction in equations.f90]
    ayprime(vc_ix) = -adotoa * vc - R_dm_gamma * vc

    end subroutine TDMPhotonScattering_PerturbationEvolve

    subroutine TDMPhotonScattering_PrintFeedback(this, FeedbackLevel)
    class(TDMPhotonScattering) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("DM-Photon scattering: u_idm_g = ",ES12.4,", n = ",I3)') &
            this%u_idm_g, this%n_idm_g
        write(*,'("                      m_DM = ",F8.2," GeV")') this%m_dm
    end if

    end subroutine TDMPhotonScattering_PrintFeedback

    end module DMPhoton
