    module DMBaryon
    use DarkMatterInteraction
    use results
    use constants
    use classes
    implicit none

    ! DM-Baryon scattering model: sigma_MT = sigma_0 * (v/c)^n_dmb
    ! Following Boddy & Gluscevic 2018
    type, extends(TDarkMatterModel) :: TDMBaryonScattering
        real(dl) :: sigma_dmb = 0._dl     ! momentum-transfer cross section sigma_0 [cm^2]
        integer  :: n_dmb = 0             ! velocity power-law index (-4,-2,0,2,4)
        real(dl) :: m_dm = 100._dl        ! DM particle mass in GeV
        real(dl) :: Nnow_stored = 0._dl   ! baryon number density today [Mpc^-3], set in Init
    contains
    procedure :: ReadParams => TDMBaryonScattering_ReadParams
    procedure, nopass :: PythonClass => TDMBaryonScattering_PythonClass
    procedure, nopass :: SelfPointer => TDMBaryonScattering_SelfPointer
    procedure :: Init => TDMBaryonScattering_Init
    procedure :: PerturbationEvolve => TDMBaryonScattering_PerturbationEvolve
    procedure :: PerturbationInitial => TDMBaryonScattering_PerturbationInitial
    procedure :: PrintFeedback => TDMBaryonScattering_PrintFeedback
    procedure :: DragRate_DM => TDMBaryonScattering_DragRate_DM
    procedure :: DragRate_Baryon => TDMBaryonScattering_DragRate_Baryon
    procedure :: VelFactor => TDMBaryonScattering_VelFactor
    end type TDMBaryonScattering

    contains

    subroutine TDMBaryonScattering_ReadParams(this, Ini)
    use IniObjects
    class(TDMBaryonScattering) :: this
    class(TIniFile), intent(in) :: Ini

    this%sigma_dmb = Ini%Read_Double('sigma_dmb', 0._dl)
    this%n_dmb = Ini%Read_Int('n_dmb', 0)
    this%m_dm = Ini%Read_Double('m_dm', 100._dl)

    end subroutine TDMBaryonScattering_ReadParams

    function TDMBaryonScattering_PythonClass()
    character(LEN=:), allocatable :: TDMBaryonScattering_PythonClass
    TDMBaryonScattering_PythonClass = 'DMBaryonScattering'
    end function TDMBaryonScattering_PythonClass

    subroutine TDMBaryonScattering_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TDMBaryonScattering), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDMBaryonScattering_SelfPointer

    subroutine TDMBaryonScattering_Init(this, State)
    class(TDMBaryonScattering), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    this%is_standard_cdm = (this%sigma_dmb == 0._dl)
    if (.not. this%is_standard_cdm) then
        this%has_cdm_velocity = .true.
        this%num_perturb_equations = 1  ! CDM velocity v_c
    else
        this%has_cdm_velocity = .false.
        this%num_perturb_equations = 0
    end if
    this%num_dr_equations = 0

    ! Store Nnow from State for use in DragRate_DM
    select type(State)
    class is (CAMBdata)
        this%Nnow_stored = State%Nnow
    end select

    end subroutine TDMBaryonScattering_Init

    function TDMBaryonScattering_VelFactor(this, a) result(Fvel)
    ! Velocity-dependent thermal average factor F(n, T)
    ! Includes the velocity-averaged momentum transfer rate coefficient
    class(TDMBaryonScattering), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: Fvel
    real(dl) :: Tcmb_now, T_b, mu_GeV, thermal_v2
    real(dl), parameter :: m_proton_GeV = 0.938272_dl
    real(dl), parameter :: kB_over_GeV = 8.617333e-14_dl  ! Boltzmann constant [GeV/K]
    real(dl), parameter :: c_light = 2.99792458e10_dl     ! cm/s

    ! Baryon temperature (approx T_CMB * (1+z)^2 / (1+z) = T_CMB/a after decoupling)
    ! Before decoupling T_b ~ T_CMB/a, after decoupling T_b ~ T_CMB/a^2
    ! Use simple approximation T_b = T_CMB / a (valid when coupled to photons)
    Tcmb_now = 2.7255_dl  ! K (COBE value)
    T_b = Tcmb_now / a

    ! Reduced mass in GeV
    mu_GeV = this%m_dm * m_proton_GeV / (this%m_dm + m_proton_GeV)

    ! Thermal velocity squared: v_th^2 = k_B T / mu (in natural units then convert)
    thermal_v2 = kB_over_GeV * T_b / mu_GeV  ! dimensionless (v/c)^2

    select case(this%n_dmb)
    case(-4)
        ! Coulomb-like: F ~ T^(-3/2)
        if (thermal_v2 > 0) then
            Fvel = (4._dl * const_pi / 3._dl) * thermal_v2**(-1.5_dl)
        else
            Fvel = 0._dl
        end if
    case(-2)
        if (thermal_v2 > 0) then
            Fvel = (2._dl / sqrt(const_pi)) * thermal_v2**(-0.5_dl)
        else
            Fvel = 0._dl
        end if
    case(0)
        ! Constant cross section
        Fvel = 1._dl
    case(2)
        Fvel = 1.5_dl * sqrt(thermal_v2)
    case(4)
        Fvel = (15._dl / 4._dl) * thermal_v2
    case default
        Fvel = 1._dl
    end select

    end function TDMBaryonScattering_VelFactor

    function TDMBaryonScattering_DragRate_DM(this, a, adotoa, grhoc_t, grhob_t) result(R_drag)
    ! Drag rate on DM from baryons [Mpc^-1] (conformal time units)
    ! Following CAMB convention: akthom = sigma_T * Nnow * Mpc, dotmu = xe * akthom / a^2
    ! Here: R_DM = sigma_dmb [cm^2] * Nnow * Mpc [cm] / a^2 * Fvel * rho_ratio
    !            = sigma_dmb * Nnow * Mpc / a^2 * Fvel * rho_b/(rho_c+rho_b)
    class(TDMBaryonScattering), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t, grhob_t
    real(dl) :: R_drag
    real(dl) :: Fvel, rho_ratio
    real(dl), parameter :: Mpc_in_cm = 3.0856776e24_dl

    if (this%sigma_dmb == 0._dl) then
        R_drag = 0._dl
        return
    end if

    Fvel = this%VelFactor(a)

    ! rho_b / (rho_c + rho_b) using grhox_t = 8*pi*G*rho*a^2
    rho_ratio = grhob_t / (grhoc_t + grhob_t)

    ! R_DM [Mpc^-1] = sigma_dmb [cm^2] * Nnow * Mpc [cm] / a^2 * Fvel * rho_ratio
    ! Same unit convention as Thomson: akthom = sigma_T * Nnow * Mpc
    R_drag = this%sigma_dmb * this%Nnow_stored * Mpc_in_cm / a**2 * Fvel * rho_ratio

    ! Cap drag rate to prevent ODE stiffness. When R >> adotoa, DM and baryons are
    ! tightly coupled (v_c ≈ v_b). The physical momentum transfer is independent of
    ! R in this limit, so capping preserves the correct physics.
    R_drag = min(R_drag, 1e4_dl * adotoa)

    end function TDMBaryonScattering_DragRate_DM

    function TDMBaryonScattering_DragRate_Baryon(this, a, adotoa, grhoc_t, grhob_t) result(R_drag)
    ! Back-reaction drag on baryons from DM: R_b_chi
    ! By Newton's third law: rho_c * R_chi_b = rho_b * R_b_chi (momentum conservation)
    ! So R_b_chi = (rho_c/rho_b) * R_chi_b
    class(TDMBaryonScattering), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t, grhob_t
    real(dl) :: R_drag, R_dm

    R_dm = this%DragRate_DM(a, adotoa, grhoc_t, grhob_t)
    ! Momentum conservation: rho_c * R_DM * (v_c - v_b) + rho_b * R_baryon * (v_b - v_c) = 0
    ! => R_baryon = (rho_c / rho_b) * R_DM
    ! But R_DM already has rho_b/(rho_c+rho_b) factor
    ! So R_baryon = R_DM * rho_c/rho_b = sigma * n_b * F * rho_c/(rho_c+rho_b) * (1/a^2)
    ! Actually simpler: R_baryon = sigma * n_DM * F / a^2 * rho_c/(rho_c+rho_b)
    ! By symmetry with R_DM but swapping b<->c in the density ratio
    if (grhob_t > 0) then
        R_drag = R_dm * (grhoc_t / grhob_t)
    else
        R_drag = 0._dl
    end if

    end function TDMBaryonScattering_DragRate_Baryon

    subroutine TDMBaryonScattering_PerturbationInitial(this, y, a, tau, k, dm_ix, dr_ix, vc_ix)
    class(TDMBaryonScattering), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer, intent(in) :: dm_ix, dr_ix, vc_ix

    ! CDM velocity starts at zero in synchronous gauge (adiabatic IC)
    if (this%has_cdm_velocity .and. vc_ix > 0) then
        y(vc_ix) = 0._dl
    end if

    end subroutine TDMBaryonScattering_PerturbationInitial

    subroutine TDMBaryonScattering_PerturbationEvolve(this, ayprime, a, adotoa, k, z, y, &
        dm_ix, dr_ix, vc_ix, clxc, vb, grhoc_t, grhob_t, sigma)
    class(TDMBaryonScattering), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, y(:)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl), intent(in) :: clxc, vb, grhoc_t, grhob_t, sigma
    real(dl) :: vc, R_dm

    if (.not. this%has_cdm_velocity .or. vc_ix <= 0) return

    vc = y(vc_ix)
    R_dm = this%DragRate_DM(a, adotoa, grhoc_t, grhob_t)

    ! CDM velocity equation: v_c' = -aH v_c - R_DM (v_c - v_b)
    ! R_dm is already capped in DragRate_DM to prevent stiffness
    ayprime(vc_ix) = -adotoa * vc - R_dm * (vc - vb)

    end subroutine TDMBaryonScattering_PerturbationEvolve

    subroutine TDMBaryonScattering_PrintFeedback(this, FeedbackLevel)
    class(TDMBaryonScattering) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("DM-Baryon scattering: sigma_0 = ",ES12.4," cm^2, n = ",I3)') &
            this%sigma_dmb, this%n_dmb
        write(*,'("                      m_DM = ",F8.2," GeV")') this%m_dm
    end if

    end subroutine TDMBaryonScattering_PrintFeedback

    end module DMBaryon
