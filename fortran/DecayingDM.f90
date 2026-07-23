    module DecayingDM
    use DarkMatterInteraction
    use results
    use constants
    use classes
    implicit none

    ! Decaying Dark Matter (DCDM) with WARM-DAUGHTER Boltzmann fluid.
    !
    ! Parent CDM -> daughter (massive) + DR (massless):
    !   m_parent -> m_daughter = epsilon * m_parent  + E_DR = (1-eps^2)/2 * m_p c^2
    !   daughter kick: v_kick = (1-eps^2)/(2*eps)  (non-rel limit, c=1)
    !
    ! Background: rho_parent ~ a^-3 * exp(-Gamma*t).
    !   rho_daughter_rest = eps * (parent_lost)
    !   rho_DR ~ a^-4 (sourced over decay history; matter-dom approx gam_inc(5/3,Gamma t))
    !
    ! Perturbations (synchronous gauge, fluid moments F_0=delta, F_1=v):
    !  - Daughter fluid: warm with time-dependent w(a), cs2(a) (Hu 1998 GDM).
    !    Free-streaming pressure suppresses small-scale clustering.
    !  - Daughter w_eff(a) interpolates between v_kick^2/3 at creation and
    !    v_kick^2/3*(a_star/a)^2 at late times (momentum redshift).
    !  - DR fluid: F_0, F_1 hierarchy with decay source from parent perturbation.
    !
    ! State variables (allocated by equations.f90 SetupScalarArrayIndices):
    !   vc_ix       -> v_d  (daughter velocity)
    !   dm_ix       -> delta_d  (daughter density perturbation)
    !   dr_ix       -> F_0  (DR monopole)
    !   dr_ix+1     -> F_1  (DR dipole = q_r = (4/3) v_r)
    !
    ! References: Blackadder & Koushiappas 2014, Audren+ 2014 (CLASS dcdm),
    !             Poulin+ 2016, Hu 1998 (GDM), Lesgourgues & Tram 2011.

    type, extends(TDarkMatterModel) :: TDecayingDM
        real(dl) :: Gamma_dcdm = 0._dl    ! decay rate [km/s/Mpc]
        real(dl) :: epsilon_dcdm = 1._dl  ! mass ratio m_daughter/m_parent
        real(dl) :: f_dcdm = 1._dl        ! fraction of CDM that decays
        ! Cached (set in Init)
        real(dl) :: Gamma_conf = 0._dl    ! Gamma in conformal units [Mpc^-1]
        real(dl) :: v_kick = 0._dl        ! (1-eps^2)/(2*eps)
        real(dl) :: w_kick = 0._dl        ! v_kick^2 / 3 (capped at 1/3)
        real(dl) :: a_star = 1._dl        ! decay scale factor: (3*H0/(2*Gamma))^(2/3)
        real(dl) :: H0_stored = 67._dl
    contains
    procedure :: ReadParams => TDecayingDM_ReadParams
    procedure, nopass :: PythonClass => TDecayingDM_PythonClass
    procedure, nopass :: SelfPointer => TDecayingDM_SelfPointer
    procedure :: Init => TDecayingDM_Init
    procedure :: BackgroundDensityAndPressure => TDecayingDM_BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy => TDecayingDM_PerturbedStressEnergy
    procedure :: PerturbationInitial => TDecayingDM_PerturbationInitial
    procedure :: PerturbationEvolve => TDecayingDM_PerturbationEvolve
    procedure :: PrintFeedback => TDecayingDM_PrintFeedback
    procedure :: SetBackgroundDensities => TDecayingDM_SetBackgroundDensities
    procedure :: DecayAux => TDecayingDM_DecayAux
    end type TDecayingDM

    contains

    subroutine TDecayingDM_ReadParams(this, Ini)
    use IniObjects
    class(TDecayingDM) :: this
    class(TIniFile), intent(in) :: Ini

    this%Gamma_dcdm = Ini%Read_Double('Gamma_dcdm', 0._dl)
    this%epsilon_dcdm = Ini%Read_Double('epsilon_dcdm', 1._dl)
    this%f_dcdm = Ini%Read_Double('f_dcdm', 1._dl)

    end subroutine TDecayingDM_ReadParams

    function TDecayingDM_PythonClass()
    character(LEN=:), allocatable :: TDecayingDM_PythonClass
    TDecayingDM_PythonClass = 'DecayingDM'
    end function TDecayingDM_PythonClass

    subroutine TDecayingDM_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TDecayingDM), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDecayingDM_SelfPointer

    subroutine TDecayingDM_DecayAux(this, a, decay_exp, frac_decayed, x, &
        w_eff, rho_DR_factor)
    ! Helper: compute auxiliary quantities at scale factor a.
    !   decay_exp    = exp(-Gamma*t(a)) (matter-dom approx)
    !   frac_decayed = 1 - decay_exp
    !   x            = Gamma*t(a) (dimensionless)
    !   w_eff        = daughter effective EOS (capped at 1/3)
    !   rho_DR_factor= dimensionless DR density factor (rho_DR*a^4 / grhodm)
    class(TDecayingDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: decay_exp, frac_decayed, x, w_eff, rho_DR_factor
    real(dl) :: a_dec_eff_sq, gam_inc
    real(dl), parameter :: G53 = 0.9027452929509336_dl
    real(dl), parameter :: x0 = 1.281_dl

    if (a <= 0._dl .or. this%Gamma_dcdm <= 0._dl) then
        decay_exp = 1._dl; frac_decayed = 0._dl; x = 0._dl
        w_eff = 0._dl; rho_DR_factor = 0._dl
        return
    end if

    ! Gamma*t(a) in matter-dom approx: t = (2/3)*H0^-1 * a^(3/2)
    x = this%Gamma_dcdm / this%H0_stored * (2._dl/3._dl) * a**1.5_dl
    decay_exp = exp(-x)
    frac_decayed = 1._dl - decay_exp

    ! Daughter effective EOS: w ~ <v^2>/3 ~ v_kick^2/3 * <(a_c/a)^2>
    ! For mostly-recent decays (a < a_star): a_c ~ a, so <(a_c/a)^2> -> 1
    ! For Gamma*t >> 1 (a > a_star): decays happened around a_star, <(a_c/a)^2> -> (a_star/a)^2
    ! Smooth interpolation:
    a_dec_eff_sq = min(1._dl, (this%a_star / a)**2)
    w_eff = min(1._dl/3._dl, this%w_kick * a_dec_eff_sq)

    ! DR density: rho_DR*a^4 = (1-eps^2)/2 * f * grhodm0 * a_star * gam_inc(5/3, Gamma t)
    ! (matter-dom approx; 1-eps^2)/2 is energy fraction lost to DR per decay)
    if (x > 1.e-30_dl) then
        gam_inc = G53 * (1._dl - exp(-(x/x0)**(5._dl/3._dl)))
        ! Dimensionless factor: rho_DR*a^2 / grhodm = (1-eps^2)/2 * f * a_star * gam_inc / a^2
        rho_DR_factor = (1._dl - this%epsilon_dcdm**2) * 0.5_dl * this%f_dcdm * &
                        this%a_star * gam_inc
    else
        rho_DR_factor = 0._dl
    end if

    end subroutine TDecayingDM_DecayAux

    subroutine TDecayingDM_Init(this, State)
    class(TDecayingDM), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    this%is_standard_cdm = (this%Gamma_dcdm == 0._dl .or. this%f_dcdm == 0._dl)

    if (.not. this%is_standard_cdm) then
        this%Gamma_conf = this%Gamma_dcdm * 1000._dl / c
        this%v_kick = (1._dl - this%epsilon_dcdm**2) / &
                      (2._dl * max(this%epsilon_dcdm, 1.e-6_dl))
        this%w_kick = this%v_kick**2 / 3._dl  ! NR limit; capped at 1/3 inside DecayAux

        this%has_cdm_velocity = .true.
        ! State allocation: vc_ix(1) + dm_ix(1) + dr_ix(2) = 4 slots
        this%num_dr_equations = 2  ! DR F_0, F_1 only (no F_2 = no anisotropic stress)
        this%num_perturb_equations = 1 + 1 + this%num_dr_equations

        select type(S => State)
        class is (CAMBdata)
            this%H0_stored = S%CP%H0
            ! Decay scale factor (3 H0 / 2 Gamma)^(2/3); the a at which Gamma*t ~ 1
            this%a_star = (1.5_dl * S%CP%H0 / max(this%Gamma_dcdm, 1.e-30_dl))**(2._dl/3._dl)
        end select
    else
        this%has_cdm_velocity = .false.
        this%num_dr_equations = 0
        this%num_perturb_equations = 0
    end if

    end subroutine TDecayingDM_Init

    subroutine TDecayingDM_BackgroundDensityAndPressure(this, grhodm, a, grhodm_t, gpres_dm)
    ! Returns total dark-matter-sector density in 8*pi*G*rho*a^2 units:
    ! (1) Stable CDM + (2) surviving parent + (3) daughter rest mass + (4) DR
    ! All wrapped into grhodm_t so that the Hubble rate stays correct.
    class(TDecayingDM), intent(inout) :: this
    real(dl), intent(in) :: grhodm, a
    real(dl), intent(out) :: grhodm_t
    real(dl), optional, intent(out) :: gpres_dm
    real(dl) :: decay_exp, frac_decayed, x, w_eff, rho_DR_factor
    real(dl) :: rho_stable_total, grhoDR_t

    if (a <= 0._dl) then
        grhodm_t = 0._dl
        if (present(gpres_dm)) gpres_dm = 0._dl
        return
    end if

    if (this%is_standard_cdm) then
        grhodm_t = grhodm / a
        if (present(gpres_dm)) gpres_dm = 0._dl
        return
    end if

    call this%DecayAux(a, decay_exp, frac_decayed, x, w_eff, rho_DR_factor)

    ! Total rest-mass-like density / a (matter-like):
    !   (1-f)*grhodm + f*(decay_exp + eps*(1-decay_exp))*grhodm
    ! = grhodm * [1 - f*(1-eps)*(1-decay_exp)]
    rho_stable_total = grhodm / a * &
        (1._dl - this%f_dcdm * (1._dl - this%epsilon_dcdm) * frac_decayed)

    ! DR contribution (radiation-like, scales as 1/a^2 in 8piG units):
    grhoDR_t = grhodm * rho_DR_factor / (a*a)

    grhodm_t = rho_stable_total + grhoDR_t

    if (present(gpres_dm)) then
        ! Pressure from daughter (warm) + DR (relativistic).
        ! Daughter pressure = w_eff * rho_d_rest
        gpres_dm = w_eff * grhodm / a * this%f_dcdm * this%epsilon_dcdm * frac_decayed &
                 + grhoDR_t / 3._dl
    end if

    end subroutine TDecayingDM_BackgroundDensityAndPressure

    subroutine TDecayingDM_SetBackgroundDensities(this, grhocrit, grhor, h2, grhodmdr, grhodr)
    class(TDecayingDM), intent(in) :: this
    real(dl), intent(in) :: grhocrit, grhor, h2
    real(dl), intent(out) :: grhodmdr, grhodr

    ! No separate density pool; daughter/DR tracked via BackgroundDensityAndPressure.
    grhodmdr = 0._dl
    grhodr = 0._dl

    end subroutine TDecayingDM_SetBackgroundDensities

    subroutine TDecayingDM_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhoc_t, adotoa, k, ay, ayprime, dm_ix, dr_ix, vc_ix)
    ! Add daughter and DR perturbations on top of the implicit grhoc_t*clxc weighting.
    ! Since grhoc_t (from BackgroundDensityAndPressure) includes daughter+DR rest-mass
    ! density, the framework's dgrho_matter = grhoc_t*clxc effectively assumes
    ! delta_d = delta_DR = clxc. We correct this by adding:
    !   dgrho_e = -(rho_d + rho_DR)*clxc + rho_d*delta_d + rho_DR*F0_DR
    !   dgq_e   = rho_d*v_d + rho_DR*F1_DR
    class(TDecayingDM), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhoc_t, adotoa, k
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl) :: decay_exp, frac_decayed, x, w_eff, rho_DR_factor
    real(dl) :: grho_d_t, grho_DR_t, grhodm_total, clxc_local
    real(dl) :: delta_d, v_d, F0_DR, F1_DR
    real(dl) :: grhodm_orig

    if (this%is_standard_cdm) then
        dgrhoe = 0; dgqe = 0; return
    end if
    if (dm_ix <= 0 .or. dr_ix <= 0 .or. vc_ix <= 0) then
        dgrhoe = 0; dgqe = 0; return
    end if

    call this%DecayAux(a, decay_exp, frac_decayed, x, w_eff, rho_DR_factor)

    ! Need the original grhodm (constant) -- recover by undoing the (1 - ...) factor:
    !   grhoc_t = grhodm/a * (1 - f*(1-eps)*frac) + grhodm * rho_DR_factor / a^2
    ! Solve numerically: grhodm = (grhoc_t - rho_DR_t) * a / (1 - f*(1-eps)*frac)
    grho_DR_t = 0._dl
    if (rho_DR_factor > 0._dl) then
        ! Two-step: compute grhodm_orig assuming grho_DR_t = 0 first, then refine.
        grhodm_orig = grhoc_t * a / max(1._dl - this%f_dcdm*(1._dl - this%epsilon_dcdm)*frac_decayed, 1.e-30_dl)
        grho_DR_t = grhodm_orig * rho_DR_factor / (a*a)
        ! One iteration to refine
        grhodm_orig = (grhoc_t - grho_DR_t) * a / &
            max(1._dl - this%f_dcdm*(1._dl - this%epsilon_dcdm)*frac_decayed, 1.e-30_dl)
        grho_DR_t = grhodm_orig * rho_DR_factor / (a*a)
    else
        grhodm_orig = grhoc_t * a / &
            max(1._dl - this%f_dcdm*(1._dl - this%epsilon_dcdm)*frac_decayed, 1.e-30_dl)
    end if

    ! Daughter rest-mass density (matter-like, 1/a):
    grho_d_t = grhodm_orig / a * this%f_dcdm * this%epsilon_dcdm * frac_decayed

    ! Read perturbation variables
    delta_d = ay(dm_ix)
    v_d     = ay(vc_ix)
    F0_DR   = ay(dr_ix)
    F1_DR   = ay(dr_ix + 1)

    ! clxc is the cold-CDM (parent + stable) perturbation; we pull it from the
    ! global state via the standard ix_clxc=2 slot (CAMB convention).
    clxc_local = ay(2)

    ! Correction: dgrho framework uses grhoc_t*clxc, which implicitly weights
    ! the daughter+DR by clxc too. Subtract that and add proper weightings.
    dgrhoe = -(grho_d_t + grho_DR_t) * clxc_local &
           + grho_d_t * (1._dl + w_eff) * delta_d &
           + grho_DR_t * F0_DR

    ! dgq: cold parent has v_c=0 in synch gauge (not added implicitly here);
    ! daughter contributes ~ rho_d*(1+w)*v_d; DR contributes rho_DR*F1.
    dgqe = grho_d_t * (1._dl + w_eff) * v_d + grho_DR_t * F1_DR

    end subroutine TDecayingDM_PerturbedStressEnergy

    subroutine TDecayingDM_PerturbationInitial(this, y, a, tau, k, dm_ix, dr_ix, vc_ix)
    ! Adiabatic IC: daughter inherits parent's perturbation; DR starts with
    ! standard radiation IC. v_d = 0 initially.
    class(TDecayingDM), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl) :: x2, x3

    if (this%is_standard_cdm) return

    ! Daughter perturbations: inherit parent's adiabatic IC.
    x2 = (k*tau)**2
    if (vc_ix > 0) y(vc_ix) = 0._dl          ! v_d = 0 (parent at rest in synch)
    if (dm_ix > 0) y(dm_ix) = -x2 / 2._dl * 0.75_dl  ! delta_d = clxc leading order

    ! DR radiation IC (synchronous gauge, leading order in x = k*tau)
    if (dr_ix > 0) then
        x3 = (k*tau)**3
        y(dr_ix) = -x2 / 3._dl       ! F_0
        y(dr_ix + 1) = -x3 / 27._dl  ! F_1 = qr
    end if

    end subroutine TDecayingDM_PerturbationInitial

    subroutine TDecayingDM_PerturbationEvolve(this, ayprime, a, adotoa, k, z, y, &
        dm_ix, dr_ix, vc_ix, clxc, vb, grhoc_t, grhob_t, sigma, high_ktau_dr)
    ! Evolve daughter fluid (delta_d, v_d) and DR hierarchy (F_0, F_1).
    ! Daughter: Hu 1998 GDM with time-dependent w(a), cs2(a) and decay source.
    ! DR: radiation hierarchy with decay source proportional to parent perturbation.
    class(TDecayingDM), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, y(:)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl), intent(in) :: clxc, vb, grhoc_t, grhob_t, sigma
    logical, intent(in), optional :: high_ktau_dr
    real(dl) :: decay_exp, frac_decayed, x, w_eff, rho_DR_factor
    real(dl) :: Gamma_a, Gamma_drag_d, Gamma_drag_DR
    real(dl) :: delta_d, v_d, F0_DR, F1_DR
    real(dl) :: cs2_d, cs2_minus_w
    real(dl) :: rho_p_over_d, rho_p_over_DR
    real(dl) :: regul

    if (this%is_standard_cdm) return
    if (dm_ix <= 0 .or. dr_ix <= 0 .or. vc_ix <= 0) return

    call this%DecayAux(a, decay_exp, frac_decayed, x, w_eff, rho_DR_factor)

    ! Conformal-time decay rate: Gamma_phys * dt/dtau = Gamma_conf * a
    Gamma_a = this%Gamma_conf * a

    ! Drag rates: source = (Gamma_a * rho_parent / rho_X) * (delta_p - delta_X)
    ! For daughter (rest mass): rho_p/rho_d_rest = decay_exp / (eps * frac_decayed)
    ! For DR: rho_p/rho_DR (massless) is much larger but the rate is bounded by
    ! Gamma_a (we cap drag at 50*adotoa for numerical stability, same threshold as CAMB TCA).
    regul = max(frac_decayed, 1.e-6_dl)
    rho_p_over_d  = decay_exp / (max(this%epsilon_dcdm, 1.e-6_dl) * regul)
    Gamma_drag_d  = Gamma_a * rho_p_over_d
    Gamma_drag_d  = min(Gamma_drag_d, 50._dl * max(adotoa, 1.e-30_dl))

    if (rho_DR_factor > 1.e-30_dl) then
        ! For DR: source = Gamma_a*(1-eps^2)/2 * rho_p / rho_DR; using DR conv F0 = delta_r.
        ! Approx normalization: source ~ Gamma_a * decay_exp / max(rho_DR_factor/a^2 ... ).
        ! Simpler: rate-cap directly.
        Gamma_drag_DR = Gamma_a * decay_exp / max(rho_DR_factor, 1.e-30_dl) * a*a * 0.5_dl
        Gamma_drag_DR = min(Gamma_drag_DR, 50._dl * max(adotoa, 1.e-30_dl))
    else
        Gamma_drag_DR = 0._dl
    end if

    delta_d = y(dm_ix)
    v_d     = y(vc_ix)
    F0_DR   = y(dr_ix)
    F1_DR   = y(dr_ix + 1)

    cs2_d = w_eff                 ! adiabatic-ish; cs2 = w for simple warm fluid
    cs2_minus_w = 0._dl           ! adiabatic limit

    ! ============ Daughter fluid (Hu 1998 GDM, F_2 -> 0 closure) ============
    ! delta_d' = -(1+w)*(k*v_d + k*z) - 3H*(cs2-w)*delta_d - Gamma_drag*(delta_d - clxc)
    ayprime(dm_ix) = -(1._dl + w_eff) * (k * v_d + k * z) &
                    - 3._dl * adotoa * cs2_minus_w * delta_d &
                    - Gamma_drag_d * (delta_d - clxc)

    ! v_d' = -H*(1-3*cs2)*v_d + cs2*k*delta_d/(1+w) - Gamma_drag*v_d
    !       (parent v_p = 0 in synch gauge -> drag is one-sided)
    ayprime(vc_ix) = -adotoa * (1._dl - 3._dl * cs2_d) * v_d &
                   + cs2_d * k * delta_d / max(1._dl + w_eff, 1.e-30_dl) &
                   - Gamma_drag_d * v_d

    ! ============ DR hierarchy (massless, decay-sourced) ============
    ! F_0' = -(4/3)*k*z - k*F_1 + Gamma_drag_DR*(clxc - F_0)
    ayprime(dr_ix) = -(4._dl/3._dl) * k * z - k * F1_DR &
                   + Gamma_drag_DR * (clxc - F0_DR)

    ! F_1' = k/3 * F_0 - Gamma_drag_DR * F_1
    !   (closure: F_2 = 0; full hierarchy would have +8/15*k*sigma source via F_2)
    ayprime(dr_ix + 1) = k / 3._dl * F0_DR - Gamma_drag_DR * F1_DR

    end subroutine TDecayingDM_PerturbationEvolve

    subroutine TDecayingDM_PrintFeedback(this, FeedbackLevel)
    class(TDecayingDM) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("Decaying DM (warm daughter): Gamma = ",ES12.4," km/s/Mpc")') this%Gamma_dcdm
        write(*,'("    epsilon = ",F8.5,", f_dcdm = ",F8.5)') &
            this%epsilon_dcdm, this%f_dcdm
        write(*,'("    v_kick = ",ES10.3,"  w_kick = ",ES10.3, "  a_star = ",ES10.3)') &
            this%v_kick, this%w_kick, this%a_star
    end if

    end subroutine TDecayingDM_PrintFeedback

    end module DecayingDM
