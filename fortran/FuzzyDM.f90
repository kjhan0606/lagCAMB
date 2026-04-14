    module FuzzyDM
    use DarkMatterInteraction
    use results
    use constants
    use classes
    implicit none

    ! Fuzzy / Ultralight Axion Dark Matter (ULDM)
    ! Ultra-light scalar field DM with m ~ 10^{-22} eV.
    ! The de Broglie wavelength is astrophysically large, leading to a quantum
    ! pressure (Jeans scale) that suppresses structure below the Jeans scale.
    !
    ! Background:
    !   At early times (m << H): field frozen, behaves as cosmological constant
    !   At late times (m >> H): field oscillates, <rho> ~ a^{-3}, behaves as CDM
    !   Transition at a_osc where m ~ 3H(a_osc)
    !
    ! Perturbations (effective fluid approach, Hlozek+ 2015):
    !   The axion field can be described as an effective fluid with
    !   scale-dependent sound speed:
    !     cs2_eff(k,a) = k^2 / (4 * m_a^2 * a^2)   [in natural units]
    !   This produces a Jeans scale:
    !     k_J(a) = (2 * m_a * a * H)^{1/2}
    !   Modes with k > k_J oscillate and are suppressed.
    !
    ! Implementation:
    !   We treat ULDM as a modification to CDM with extra perturbation equations
    !   (density and velocity) that include the effective sound speed.
    !   The background is standard CDM (since ULDM behaves as CDM at late times
    !   when m >> H), with corrections at early times.
    !
    ! References:
    !   Hu+ 2000, Hlozek+ 2015, 2018, axionCAMB (Grin+ 2022), AxiECAMB
    !
    ! Parameters:
    !   m_axion: axion mass [eV] (typically 10^{-22} to 10^{-18})
    !   omega_axion_h2: axion density parameter Omega_a h^2
    !   f_axion: fraction of CDM that is ULDM (alternative to omega_axion_h2)

    type, extends(TDarkMatterModel) :: TFuzzyDM
        real(dl) :: m_axion = 0._dl        ! axion mass [eV]
        real(dl) :: omega_axion_h2 = 0._dl ! axion density Omega_a h^2 (0 = use f_axion*omch2)
        real(dl) :: f_axion = 0._dl        ! fraction of CDM that is axion (0 to 1)
        ! Cached quantities
        real(dl) :: m_conf = 0._dl         ! m_axion in CAMB units [Mpc^-1]
        real(dl) :: a_osc = 0._dl          ! scale factor at oscillation onset
        real(dl) :: grho_ax0 = 0._dl       ! axion background density (8*pi*G*rho*a^4 at a=1)
    contains
    procedure :: ReadParams => TFuzzyDM_ReadParams
    procedure, nopass :: PythonClass => TFuzzyDM_PythonClass
    procedure, nopass :: SelfPointer => TFuzzyDM_SelfPointer
    procedure :: Init => TFuzzyDM_Init
    procedure :: PerturbedStressEnergy => TFuzzyDM_PerturbedStressEnergy
    procedure :: PerturbationInitial => TFuzzyDM_PerturbationInitial
    procedure :: PerturbationEvolve => TFuzzyDM_PerturbationEvolve
    procedure :: PrintFeedback => TFuzzyDM_PrintFeedback
    procedure :: cs2_eff => TFuzzyDM_cs2_eff
    procedure :: JeansScale => TFuzzyDM_JeansScale
    end type TFuzzyDM

    contains

    subroutine TFuzzyDM_ReadParams(this, Ini)
    use IniObjects
    class(TFuzzyDM) :: this
    class(TIniFile), intent(in) :: Ini

    this%m_axion = Ini%Read_Double('m_axion', 0._dl)
    this%omega_axion_h2 = Ini%Read_Double('omega_axion_h2', 0._dl)
    this%f_axion = Ini%Read_Double('f_axion', 0._dl)

    end subroutine TFuzzyDM_ReadParams

    function TFuzzyDM_PythonClass()
    character(LEN=:), allocatable :: TFuzzyDM_PythonClass
    TFuzzyDM_PythonClass = 'FuzzyDM'
    end function TFuzzyDM_PythonClass

    subroutine TFuzzyDM_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TFuzzyDM), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TFuzzyDM_SelfPointer

    function TFuzzyDM_cs2_eff(this, k, a) result(cs2)
    ! Effective sound speed squared: cs2 = k^2 / (4 * m^2 * a^2)
    ! In natural units (hbar=c=1): m in eV, k in eV (but we need CAMB units)
    ! In CAMB units: k [Mpc^-1], m_conf [Mpc^-1]
    ! cs2 = (k / (2*m_conf*a))^2 = k^2 / (4 * m_conf^2 * a^2)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: k, a
    real(dl) :: cs2

    if (this%m_conf > 0._dl .and. a > 0._dl) then
        cs2 = k**2 / (4._dl * this%m_conf**2 * a**2)
        ! Cap cs2 at 1 (cannot exceed speed of light)
        cs2 = min(cs2, 1._dl)
    else
        cs2 = 0._dl
    end if

    end function TFuzzyDM_cs2_eff

    function TFuzzyDM_JeansScale(this, a, adotoa) result(k_J)
    ! Jeans wavenumber: k_J = (2*m*a*H)^{1/2} * a
    ! In CAMB conformal units: k_J = sqrt(2 * m_conf * a * adotoa) * a
    ! Actually: k_J^2 = 4*m^2*a^2 * H*a (from cs2*k_J^2 = ... gravitational instability)
    ! More precisely: k_J = (4*m_conf*a^2*adotoa)^{1/2}  [Mpc^-1]
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(in) :: a, adotoa
    real(dl) :: k_J

    if (this%m_conf > 0._dl) then
        k_J = sqrt(4._dl * this%m_conf * a**2 * adotoa)
    else
        k_J = 0._dl
    end if

    end function TFuzzyDM_JeansScale

    subroutine TFuzzyDM_Init(this, State)
    class(TFuzzyDM), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) :: h2, hbar_eV_s, eV_per_Mpc

    this%is_standard_cdm = (this%m_axion <= 0._dl .and. this%f_axion <= 0._dl)

    if (.not. this%is_standard_cdm) then
        ! Convert axion mass from eV to CAMB units [Mpc^-1]
        ! m [eV] -> m [s^-1] = m[eV] / hbar[eV*s]
        ! m [Mpc^-1] = m[s^-1] / c[Mpc/s] = m[eV] / (hbar*c) [eV*Mpc]
        ! hbar*c = 1.9733e-7 eV*m = 1.9733e-7 * 3.0857e22 eV*Mpc = 6.5822e-16 eV*s * 3e8 m/s
        ! Actually: hbar*c = 1.9733e-7 eV*m, and 1 Mpc = 3.0857e22 m
        ! So hbar*c [eV*Mpc] = 1.9733e-7 * 3.0857e22 = 6.0857e15 eV*Mpc? No...
        ! hbar*c = 197.3 MeV*fm = 197.3e-15 m * 1e6 eV = 1.973e-7 eV*m
        ! 1 Mpc = 3.0857e22 m
        ! So: eV_per_Mpc_inv = hbar*c / Mpc = 1.973e-7 / 3.0857e22 = 6.4e-30 eV
        ! Therefore: m [Mpc^-1] = m [eV] / 6.4e-30 eV = m[eV] * 1.563e29
        ! This is a huge number. For m = 1e-22 eV: m_conf = 1.563e7 Mpc^-1
        eV_per_Mpc = 1.973e-7_dl / (3.0857e22_dl) ! hbar*c in eV*Mpc -> nope
        ! Let me be more careful:
        ! hbar = 6.582e-16 eV*s, c = 2.998e8 m/s = 2.998e8 / 3.0857e22 Mpc/s = 9.716e-15 Mpc/s
        ! m [Mpc^-1] = m [eV] / (hbar [eV*s] * c [Mpc/s])
        ! hbar*c = 6.582e-16 * 9.716e-15 = 6.396e-30 eV*Mpc
        ! m_conf = m_axion / (hbar*c [eV*Mpc]) = m_axion / 6.396e-30
        this%m_conf = this%m_axion / 6.396e-30_dl

        ! CDM velocity + axion density perturbation
        this%has_cdm_velocity = .true.
        ! Equations: v_axion (1) + delta_axion (1) = 2
        ! We use dm_ix for delta_axion, vc_ix for v_axion
        this%num_dr_equations = 1  ! delta_axion stored at dm_ix
        this%num_perturb_equations = 1 + 1  ! v_axion + delta_axion

        ! Scale factor at oscillation onset: m ~ 3*H(a_osc)
        ! H(a_osc) ~ H_0 * sqrt(Omega_r/a_osc^4 + Omega_m/a_osc^3)
        ! For radiation domination: a_osc ~ (m / (3*H_eq))^{1/2} * a_eq
        ! Approximate: a_osc from m = 3*H_0*sqrt(Omega_r)/a_osc^2
        select type(S => State)
        class is (CAMBdata)
            h2 = (S%CP%H0/100._dl)**2
            if (this%omega_axion_h2 > 0._dl) then
                this%grho_ax0 = S%grhocrit * this%omega_axion_h2 / h2
            else
                this%grho_ax0 = S%grhocrit * S%CP%omch2 / h2 * this%f_axion
            end if
            ! a_osc: m = 3*H(a_osc), approximate in radiation domination
            ! H^2 = (8*pi*G/3) * rho_rad / a^4 => H = sqrt(grhog+grhor)/(sqrt(3)*a^2)
            ! m_conf = 3 * sqrt(grhog+grhornomass+sum(grhormass)) / (sqrt(3)*a_osc^2)
            ! a_osc^2 = 3*sqrt((grhog+grhornomass)/3) / m_conf
            this%a_osc = (3._dl * sqrt((S%grhog + S%grhornomass) / 3._dl) / this%m_conf)**0.5_dl
            this%a_osc = min(this%a_osc, 1._dl)
        end select
    else
        this%has_cdm_velocity = .false.
        this%num_perturb_equations = 0
        this%num_dr_equations = 0
    end if

    end subroutine TFuzzyDM_Init

    subroutine TFuzzyDM_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhoc_t, adotoa, k, ay, ayprime, dm_ix, dr_ix, vc_ix)
    class(TFuzzyDM), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhoc_t, adotoa, k
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl) :: grho_ax_t, f_ax

    if (this%is_standard_cdm) then
        dgrhoe = 0
        dgqe = 0
        return
    end if

    ! Axion fraction of total CDM
    f_ax = this%f_axion
    if (this%omega_axion_h2 > 0._dl .and. grhoc_t > 0._dl) then
        grho_ax_t = this%grho_ax0 / a  ! axion density * a^2
        f_ax = grho_ax_t / (grhoc_t + grho_ax_t)
    end if

    dgrhoe = 0._dl
    dgqe = 0._dl

    ! Axion replaces f_ax fraction of CDM. Standard dgrho already has grhoc_t*clxc
    ! which includes the axion part. Correction: f*(delta_ax - clxc)
    ! ix_clxc = 2 is a fixed index in equations.f90
    if (dm_ix > 0) then
        dgrhoe = grhoc_t * f_ax * (ay(dm_ix) - ay(2))
    end if

    ! Axion velocity (CDM velocity = 0 in synchronous gauge, so full v_ax is correction)
    if (vc_ix > 0) then
        dgqe = grhoc_t * f_ax * ay(vc_ix)
    end if

    end subroutine TFuzzyDM_PerturbedStressEnergy

    subroutine TFuzzyDM_PerturbationInitial(this, y, a, tau, k, dm_ix, dr_ix, vc_ix)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer, intent(in) :: dm_ix, dr_ix, vc_ix

    if (this%is_standard_cdm) return

    ! Adiabatic IC: axion starts frozen (like cosmological constant) if a < a_osc
    ! For a > a_osc, axion IC same as CDM
    if (dm_ix > 0) then
        if (a < this%a_osc) then
            y(dm_ix) = 0._dl  ! frozen field: no density perturbation
        else
            y(dm_ix) = -0.5_dl * (k*tau)**2  ! same as CDM adiabatic
        end if
    end if

    if (vc_ix > 0) y(vc_ix) = 0._dl

    ! dr_ix is allocated but unused by FuzzyDM; initialize to 0
    if (dr_ix > 0) y(dr_ix) = 0._dl

    end subroutine TFuzzyDM_PerturbationInitial

    subroutine TFuzzyDM_PerturbationEvolve(this, ayprime, a, adotoa, k, z, y, &
        dm_ix, dr_ix, vc_ix, clxc, vb, grhoc_t, grhob_t, sigma)
    class(TFuzzyDM), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, y(:)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl), intent(in) :: clxc, vb, grhoc_t, grhob_t, sigma
    real(dl) :: cs2, delta_ax, v_ax, k_J

    if (this%is_standard_cdm) return
    if (dm_ix <= 0 .or. vc_ix <= 0) return

    ! dr_ix allocated but unused; set derivative to 0
    if (dr_ix > 0) ayprime(dr_ix) = 0._dl

    delta_ax = y(dm_ix)
    v_ax = y(vc_ix)

    ! Effective sound speed (quantum pressure)
    cs2 = this%cs2_eff(k, a)

    ! Before oscillation onset (a < a_osc), the field is frozen
    ! Perturbations don't grow; delta_ax ~ 0, v_ax ~ 0
    if (a < this%a_osc) then
        ayprime(dm_ix) = 0._dl
        ayprime(vc_ix) = 0._dl
        return
    end if

    ! After oscillation onset: behaves as fluid with effective sound speed
    ! delta_ax' = -k*v_ax - k*z (synchronous gauge)
    ayprime(dm_ix) = -k * v_ax - k * z

    ! v_ax' = -aH*v_ax + cs2_eff * k * delta_ax
    ! The cs2 term provides the quantum Jeans pressure
    ayprime(vc_ix) = -adotoa * v_ax + cs2 * k * delta_ax

    end subroutine TFuzzyDM_PerturbationEvolve

    subroutine TFuzzyDM_PrintFeedback(this, FeedbackLevel)
    class(TFuzzyDM) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("Fuzzy DM: m_axion = ",ES12.4," eV")') this%m_axion
        if (this%omega_axion_h2 > 0) then
            write(*,'("          Omega_a h^2 = ",ES12.4)') this%omega_axion_h2
        else
            write(*,'("          f_axion = ",F8.5," (fraction of CDM)")') this%f_axion
        end if
        write(*,'("          a_osc = ",ES10.3)') this%a_osc
    end if

    end subroutine TFuzzyDM_PrintFeedback

    end module FuzzyDM
