    module MultiInteractingDM
    use DarkMatterInteraction
    use results
    use constants
    use classes
    implicit none

    integer, parameter :: max_ethos_coeffs_multi = 6

    ! Multi-Interacting DM: single DM fluid with up to 3 simultaneous interaction channels
    ! 1. DR channel (ETHOS-style): DM-Dark Radiation Boltzmann hierarchy
    ! 2. Baryon channel: DM-Baryon momentum transfer
    ! 3. Photon channel: DM-Photon opacity addition to photon hierarchy
    !
    ! References: Becker+ 2021 (arXiv:2010.04074)
    type, extends(TDarkMatterModel) :: TMultiInteractingDM
        ! DR channel (ETHOS-based)
        real(dl) :: omdmdrh2 = 0._dl
        real(dl) :: N_dark = 0._dl
        real(dl) :: a_dark(max_ethos_coeffs_multi) = 0._dl
        real(dl) :: alpha_l_dark = 1._dl
        integer  :: lmax_dr = 15
        ! Baryon channel
        real(dl) :: sigma_dmb = 0._dl      ! [cm^2]
        integer  :: n_dmb = 0
        ! Photon channel
        real(dl) :: u_idm_g = 0._dl
        integer  :: n_idm_g = 0
        ! Common
        real(dl) :: m_dm = 100._dl         ! [GeV]
        real(dl) :: cs2_dm = 0._dl
        ! Channel flags (set in Init)
        logical :: has_dr_channel = .false.
        logical :: has_baryon_channel = .false.
        logical :: has_photon_channel = .false.
        ! Cached background densities
        real(dl) :: grhodmdr_bg = 0._dl
        real(dl) :: grhodr_bg = 0._dl
        real(dl) :: Nnow_stored = 0._dl
        real(dl) :: grhog_stored = 0._dl
    contains
    procedure :: ReadParams => TMultiInteractingDM_ReadParams
    procedure, nopass :: PythonClass => TMultiInteractingDM_PythonClass
    procedure, nopass :: SelfPointer => TMultiInteractingDM_SelfPointer
    procedure :: Init => TMultiInteractingDM_Init
    procedure :: PerturbedStressEnergy => TMultiInteractingDM_PerturbedStressEnergy
    procedure :: PerturbationInitial => TMultiInteractingDM_PerturbationInitial
    procedure :: PerturbationEvolve => TMultiInteractingDM_PerturbationEvolve
    procedure :: PrintFeedback => TMultiInteractingDM_PrintFeedback
    procedure :: DragRate_Baryon => TMultiInteractingDM_DragRate_Baryon
    procedure :: DragRate_Photon => TMultiInteractingDM_DragRate_Photon
    procedure :: SetBackgroundDensities => TMultiInteractingDM_SetBackgroundDensities
    procedure :: Opacity_DMDR => TMultiInteractingDM_Opacity
    procedure :: BaryonDragRate_DM => TMultiInteractingDM_BaryonDragRate_DM
    procedure :: PhotonMuDot => TMultiInteractingDM_PhotonMuDot
    end type TMultiInteractingDM

    contains

    subroutine TMultiInteractingDM_ReadParams(this, Ini)
    use IniObjects
    class(TMultiInteractingDM) :: this
    class(TIniFile), intent(in) :: Ini
    integer :: i

    this%omdmdrh2 = Ini%Read_Double('omdmdrh2', 0._dl)
    this%N_dark = Ini%Read_Double('N_dark', 0._dl)
    this%lmax_dr = Ini%Read_Int('lmax_dr', 15)
    this%alpha_l_dark = Ini%Read_Double('alpha_l_dark', 1._dl)
    this%cs2_dm = Ini%Read_Double('cs2_dm', 0._dl)
    do i = 1, max_ethos_coeffs_multi
        this%a_dark(i) = 0._dl
    end do
    this%sigma_dmb = Ini%Read_Double('sigma_dmb', 0._dl)
    this%n_dmb = Ini%Read_Int('n_dmb', 0)
    this%u_idm_g = Ini%Read_Double('u_idm_g', 0._dl)
    this%n_idm_g = Ini%Read_Int('n_idm_g', 0)
    this%m_dm = Ini%Read_Double('m_dm', 100._dl)

    end subroutine TMultiInteractingDM_ReadParams

    function TMultiInteractingDM_PythonClass()
    character(LEN=:), allocatable :: TMultiInteractingDM_PythonClass
    TMultiInteractingDM_PythonClass = 'MultiInteractingDM'
    end function TMultiInteractingDM_PythonClass

    subroutine TMultiInteractingDM_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TMultiInteractingDM), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TMultiInteractingDM_SelfPointer

    subroutine TMultiInteractingDM_Init(this, State)
    class(TMultiInteractingDM), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    ! Determine active channels
    this%has_dr_channel = (this%omdmdrh2 > 0._dl .or. this%N_dark > 0._dl)
    this%has_baryon_channel = (this%sigma_dmb > 0._dl)
    this%has_photon_channel = (this%u_idm_g > 0._dl)

    this%is_standard_cdm = (.not. this%has_dr_channel .and. .not. this%has_baryon_channel &
        .and. .not. this%has_photon_channel)

    if (.not. this%is_standard_cdm) then
        this%has_cdm_velocity = .true.

        if (this%has_dr_channel) then
            this%num_dr_equations = this%lmax_dr + 1
            ! v_dm + delta_dm (for DR channel) + DR hierarchy
            this%num_perturb_equations = 1 + 1 + this%num_dr_equations
        else
            this%num_dr_equations = 0
            ! Just v_dm (using standard clxc for density)
            this%num_perturb_equations = 1
        end if

        select type(S => State)
        class is (CAMBdata)
            this%grhodmdr_bg = S%grhodmdr
            this%grhodr_bg = S%grhodr
            this%Nnow_stored = S%Nnow
            this%grhog_stored = S%grhog
        end select
    else
        this%has_cdm_velocity = .false.
        this%num_dr_equations = 0
        this%num_perturb_equations = 0
    end if

    end subroutine TMultiInteractingDM_Init

    function TMultiInteractingDM_Opacity(this, a) result(kappa)
    class(TMultiInteractingDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: kappa
    real(dl) :: T_ratio
    integer :: i

    T_ratio = 1._dl / a
    kappa = 0._dl
    do i = 1, max_ethos_coeffs_multi
        if (this%a_dark(i) /= 0._dl) then
            kappa = kappa + this%a_dark(i) * T_ratio**(i+1)
        end if
    end do

    end function TMultiInteractingDM_Opacity

    function TMultiInteractingDM_BaryonDragRate_DM(this, a, adotoa, grhoc_t, grhob_t) result(R_drag)
    ! Drag rate on DM from baryon scattering (same physics as DMBaryonScattering)
    class(TMultiInteractingDM), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t, grhob_t
    real(dl) :: R_drag
    real(dl) :: Fvel, rho_ratio
    real(dl), parameter :: Mpc_in_cm = 3.0856776e24_dl
    real(dl), parameter :: m_proton_GeV = 0.938272_dl
    real(dl), parameter :: kB_over_GeV = 8.617333e-14_dl
    real(dl) :: Tcmb_now, T_b, mu_GeV, thermal_v2

    if (.not. this%has_baryon_channel) then
        R_drag = 0._dl
        return
    end if

    Tcmb_now = 2.7255_dl
    T_b = Tcmb_now / a
    mu_GeV = this%m_dm * m_proton_GeV / (this%m_dm + m_proton_GeV)
    thermal_v2 = kB_over_GeV * T_b / mu_GeV

    select case(this%n_dmb)
    case(-4)
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
        Fvel = 1._dl
    case(2)
        Fvel = 1.5_dl * sqrt(thermal_v2)
    case(4)
        Fvel = (15._dl / 4._dl) * thermal_v2
    case default
        Fvel = 1._dl
    end select

    rho_ratio = grhob_t / (grhoc_t + grhob_t)
    R_drag = this%sigma_dmb * this%Nnow_stored * Mpc_in_cm / a**2 * Fvel * rho_ratio
    R_drag = min(R_drag, 1e4_dl * adotoa)

    end function TMultiInteractingDM_BaryonDragRate_DM

    function TMultiInteractingDM_DragRate_Baryon(this, a, adotoa, grhoc_t, grhob_t) result(R_drag)
    ! Back-reaction on baryons from DM
    class(TMultiInteractingDM), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t, grhob_t
    real(dl) :: R_drag, R_dm

    if (.not. this%has_baryon_channel) then
        R_drag = 0._dl
        return
    end if

    R_dm = this%BaryonDragRate_DM(a, adotoa, grhoc_t, grhob_t)
    if (grhob_t > 0) then
        R_drag = R_dm * (grhoc_t / grhob_t)
    else
        R_drag = 0._dl
    end if

    end function TMultiInteractingDM_DragRate_Baryon

    function TMultiInteractingDM_PhotonMuDot(this, a, adotoa, grhoc_t) result(mu_dot)
    ! DM-photon opacity (same physics as DMPhotonScattering)
    class(TMultiInteractingDM), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t
    real(dl) :: mu_dot
    real(dl), parameter :: sigma_Th_cm2 = 6.6524587e-25_dl
    real(dl), parameter :: Mpc_in_cm = 3.0856776e24_dl
    real(dl), parameter :: GeV_to_g = 1.78266192e-24_dl
    real(dl), parameter :: c_cm_s = 2.99792458e10_dl
    real(dl), parameter :: G_cgs = 6.67430e-8_dl
    real(dl) :: sigma_cm2, rho_c_cgs, n_dm_phys, m_dm_g

    if (.not. this%has_photon_channel) then
        mu_dot = 0._dl
        return
    end if

    sigma_cm2 = sigma_Th_cm2 * this%u_idm_g * (this%m_dm / 100._dl)
    sigma_cm2 = sigma_cm2 / a**this%n_idm_g

    m_dm_g = this%m_dm * GeV_to_g
    rho_c_cgs = grhoc_t * c_cm_s**2 / (8._dl * const_pi * G_cgs * Mpc_in_cm**2 * a**2)
    n_dm_phys = rho_c_cgs / m_dm_g

    mu_dot = n_dm_phys * sigma_cm2 * Mpc_in_cm
    mu_dot = min(mu_dot, 1e4_dl * adotoa)

    end function TMultiInteractingDM_PhotonMuDot

    function TMultiInteractingDM_DragRate_Photon(this, a, adotoa, grhoc_t, grhog_t) result(mu_dot)
    class(TMultiInteractingDM), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t, grhog_t
    real(dl) :: mu_dot

    mu_dot = this%PhotonMuDot(a, adotoa, grhoc_t)

    end function TMultiInteractingDM_DragRate_Photon

    subroutine TMultiInteractingDM_SetBackgroundDensities(this, grhocrit, grhor, h2, grhodmdr, grhodr)
    class(TMultiInteractingDM), intent(in) :: this
    real(dl), intent(in) :: grhocrit, grhor, h2
    real(dl), intent(out) :: grhodmdr, grhodr

    if (this%has_dr_channel) then
        grhodmdr = grhocrit * this%omdmdrh2 / h2
        grhodr = grhor * this%N_dark
    else
        grhodmdr = 0._dl
        grhodr = 0._dl
    end if

    end subroutine TMultiInteractingDM_SetBackgroundDensities

    subroutine TMultiInteractingDM_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhoc_t, adotoa, k, ay, ayprime, dm_ix, dr_ix, vc_ix)
    class(TMultiInteractingDM), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhoc_t, adotoa, k
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl) :: grhodmdr_t, grhodr_t, v_dm, F0_dr

    if (this%is_standard_cdm) then
        dgrhoe = 0
        dgqe = 0
        return
    end if

    dgrhoe = 0._dl
    dgqe = 0._dl
    v_dm = 0._dl
    if (vc_ix > 0) v_dm = ay(vc_ix)

    if (this%has_dr_channel .and. dm_ix > 0 .and. dr_ix > 0) then
        grhodmdr_t = this%grhodmdr_bg / a
        grhodr_t = this%grhodr_bg / (a*a)

        dgrhoe = grhodmdr_t * ay(dm_ix) + grhodr_t * ay(dr_ix)
        ! F_1 follows CAMB qr convention (see DMDR_ETHOS); no factor of 3.
        dgqe = grhodmdr_t * v_dm + grhodr_t * ay(dr_ix + 1)
    end if

    end subroutine TMultiInteractingDM_PerturbedStressEnergy

    subroutine TMultiInteractingDM_PerturbationInitial(this, y, a, tau, k, dm_ix, dr_ix, vc_ix)
    class(TMultiInteractingDM), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl) :: x2, x3

    if (this%is_standard_cdm) return

    x2 = (k*tau)**2
    x3 = (k*tau)**3

    if (vc_ix > 0) y(vc_ix) = 0._dl

    if (this%has_dr_channel) then
        if (dm_ix > 0) y(dm_ix) = -x2/2._dl * 0.75_dl
        if (dr_ix > 0) then
            ! Match CAMB massless-neutrino IC (synchronous gauge, leading order)
            y(dr_ix) = -x2/3._dl
            y(dr_ix + 1) = -x3/27._dl
        end if
    end if

    end subroutine TMultiInteractingDM_PerturbationInitial

    subroutine TMultiInteractingDM_PerturbationEvolve(this, ayprime, a, adotoa, k, z, y, &
        dm_ix, dr_ix, vc_ix, clxc, vb, grhoc_t, grhob_t, sigma, high_ktau_dr)
    class(TMultiInteractingDM), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, y(:)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl), intent(in) :: clxc, vb, grhoc_t, grhob_t, sigma
    logical, intent(in), optional :: high_ktau_dr
    real(dl) :: v_dm, clx_dm, kappa, Gamma_DR, Gamma_DR_tilde
    real(dl) :: Gamma_b, Gamma_g, mu_dot_ph, grhog_t
    real(dl) :: grhodmdr_t, grhodr_t
    real(dl) :: denl_l, denl2_l, cothxor, a2
    real(dl) :: R4, Vdot
    integer :: l, ix
    logical :: use_ufa, use_tca

    if (this%is_standard_cdm) return
    if (vc_ix <= 0) return

    a2 = a*a
    cothxor = adotoa  ! ~ 1/tau in radiation era; matches DMDR_ETHOS convention
    v_dm = y(vc_ix)

    ! Get DM density perturbation
    if (this%has_dr_channel .and. dm_ix > 0) then
        clx_dm = y(dm_ix)
    else
        clx_dm = clxc
    end if

    ! ============ DM velocity equation ============
    ayprime(vc_ix) = -adotoa * v_dm + this%cs2_dm * k * clx_dm

    ! --- DR channel ---
    if (this%has_dr_channel .and. dm_ix > 0 .and. dr_ix > 0) then
        grhodmdr_t = this%grhodmdr_bg / a
        grhodr_t = this%grhodr_bg / a2

        kappa = this%Opacity_DMDR(a)
        if (grhodmdr_t > 0 .and. grhodr_t > 0) then
            Gamma_DR = a * kappa * grhodr_t / grhodmdr_t
            Gamma_DR_tilde = a * kappa * grhodmdr_t / grhodr_t * (4._dl/3._dl)
        else
            Gamma_DR = 0._dl
            Gamma_DR_tilde = 0._dl
        end if

        ! DR Boltzmann hierarchy
        use_ufa = .false.
        if (present(high_ktau_dr)) use_ufa = high_ktau_dr
        ! TCA: when collision rates dominate, v_dm and F_1/4 are locked
        R4 = (4._dl/3._dl) * grhodr_t / max(grhodmdr_t, 1e-30_dl)
        use_tca = ((Gamma_DR + Gamma_DR_tilde/4._dl) > 50._dl * k)

        if (use_tca) then
            ! ============ TIGHT COUPLING ============
            ayprime(dm_ix) = -k * v_dm - k * z
            Vdot = (-adotoa * v_dm + this%cs2_dm * k * clx_dm &
                + R4 * k * y(dr_ix) / 12._dl) / (1._dl + R4)
            ayprime(vc_ix) = Vdot
            ayprime(dr_ix) = -(4._dl/3._dl) * k * z - 4._dl * k * v_dm
            ayprime(dr_ix + 1) = 4._dl * Vdot
            if (use_ufa) then
                ayprime(dr_ix + 2) = 0._dl
            else
                do l = 2, this%lmax_dr
                    ayprime(dr_ix + l) = 0._dl
                end do
            end if
        else if (use_ufa) then
            ! UFA: only evolve F_0, F_1, F_2
            ayprime(dm_ix) = -k * v_dm - k * z
            ayprime(vc_ix) = ayprime(vc_ix) - Gamma_DR * (v_dm - y(dr_ix + 1) / 4._dl)
            ayprime(dr_ix) = -(4._dl/3._dl)*k*z - k*y(dr_ix+1)
            ayprime(dr_ix+1) = k/3._dl*(y(dr_ix) - 2._dl*y(dr_ix+2)) &
                + Gamma_DR_tilde*(v_dm - y(dr_ix+1)/4._dl)
            ayprime(dr_ix+2) = -3._dl*cothxor*y(dr_ix+2) - ayprime(dr_ix) &
                + 8._dl/15._dl*k*sigma &
                - this%alpha_l_dark * Gamma_DR_tilde * y(dr_ix+2)
        else
            ayprime(dm_ix) = -k * v_dm - k * z
            ayprime(vc_ix) = ayprime(vc_ix) - Gamma_DR * (v_dm - y(dr_ix + 1) / 4._dl)
            ayprime(dr_ix) = -(4._dl/3._dl)*k*z - k*y(dr_ix+1)

            if (this%lmax_dr >= 2) then
                ayprime(dr_ix + 1) = k/3._dl * (y(dr_ix) - 2._dl * y(dr_ix + 2)) &
                    + Gamma_DR_tilde * (v_dm - y(dr_ix + 1) / 4._dl)
            else
                ayprime(dr_ix + 1) = k/3._dl * y(dr_ix) &
                    + Gamma_DR_tilde * (v_dm - y(dr_ix + 1) / 4._dl)
            end if

            do l = 2, this%lmax_dr - 1
                ix = dr_ix + l
                denl_l = k * real(l, dl) / real(2*l+1, dl)
                denl2_l = k * real(l+1, dl) / real(2*l+1, dl)
                ayprime(ix) = denl_l * y(ix-1) - denl2_l * y(ix+1) &
                    - this%alpha_l_dark * Gamma_DR_tilde * y(ix)
                if (l == 2) ayprime(ix) = ayprime(ix) + 8._dl/15._dl * k * sigma
            end do

            if (this%lmax_dr >= 2) then
                ix = dr_ix + this%lmax_dr
                ayprime(ix) = k * y(ix-1) - real(this%lmax_dr + 1, dl) * cothxor * y(ix) &
                    - this%alpha_l_dark * Gamma_DR_tilde * y(ix)
            end if
        end if
    end if

    ! --- Baryon channel ---
    if (this%has_baryon_channel) then
        Gamma_b = this%BaryonDragRate_DM(a, adotoa, grhoc_t, grhob_t)
        ayprime(vc_ix) = ayprime(vc_ix) - Gamma_b * (v_dm - vb)
    end if

    ! --- Photon channel ---
    ! v_dm += -R_gamma * v_dm (qg/4 correction added in equations.f90)
    if (this%has_photon_channel) then
        grhog_t = this%grhog_stored / a2
        mu_dot_ph = this%PhotonMuDot(a, adotoa, grhoc_t)
        Gamma_g = 4._dl/3._dl * grhog_t / max(grhoc_t, 1e-30_dl) * mu_dot_ph
        ayprime(vc_ix) = ayprime(vc_ix) - Gamma_g * v_dm
    end if

    end subroutine TMultiInteractingDM_PerturbationEvolve

    subroutine TMultiInteractingDM_PrintFeedback(this, FeedbackLevel)
    class(TMultiInteractingDM) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("Multi-Interacting DM: m_DM = ",F8.2," GeV")') this%m_dm
        if (this%has_dr_channel) then
            write(*,'("  DR channel: omdmdrh2 = ",ES12.4,", N_dark = ",F8.4)') &
                this%omdmdrh2, this%N_dark
        end if
        if (this%has_baryon_channel) then
            write(*,'("  Baryon channel: sigma_dmb = ",ES12.4," cm^2, n = ",I3)') &
                this%sigma_dmb, this%n_dmb
        end if
        if (this%has_photon_channel) then
            write(*,'("  Photon channel: u_idm_g = ",ES12.4,", n = ",I3)') &
                this%u_idm_g, this%n_idm_g
        end if
    end if

    end subroutine TMultiInteractingDM_PrintFeedback

    end module MultiInteractingDM
