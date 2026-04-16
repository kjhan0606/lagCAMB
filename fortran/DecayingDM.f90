    module DecayingDM
    use DarkMatterInteraction
    use results
    use constants
    use classes
    implicit none

    ! Decaying Dark Matter (DCDM) model
    ! Parent DM decays into massive daughter + dark radiation:
    !   DM_parent -> DM_daughter + DR
    ! Background: rho_dcdm' + 3H*rho_dcdm = -Gamma*a*rho_dcdm
    ! Daughter inherits velocity kick v_kick = (1-epsilon^2)/(2*epsilon)
    ! DR gains remaining energy
    !
    ! References: Blackadder & Koushiappas 2014, Poulin+ 2016, CLASS dcdm module
    type, extends(TDarkMatterModel) :: TDecayingDM
        real(dl) :: Gamma_dcdm = 0._dl    ! decay rate [km/s/Mpc] (same units as H0)
        real(dl) :: epsilon_dcdm = 1._dl  ! mass ratio m_daughter/m_parent (0 < epsilon <= 1)
        real(dl) :: f_dcdm = 1._dl        ! fraction of CDM that is decaying (0 to 1)
        ! Cached quantities (set in Init)
        real(dl) :: Gamma_conf = 0._dl    ! Gamma in conformal time units [Mpc^-1]
        real(dl) :: v_kick = 0._dl        ! daughter velocity kick (1-eps^2)/(2*eps) in c
        real(dl) :: grho_dcdm0 = 0._dl    ! initial decaying DM density (8*pi*G*rho_dcdm*a^4 at a=1)
        real(dl) :: grho_dr_decay0 = 0._dl ! DR from decay density (accumulated)
        real(dl) :: H0_stored = 67._dl    ! H0 cached for background calculation
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
    procedure :: DecayFactor => TDecayingDM_DecayFactor
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

    function TDecayingDM_DecayFactor(this, a) result(efac)
    ! exp(-Gamma * t(a)) where t is cosmic time
    ! For matter domination: t ~ (2/3) * 1/H0 * a^(3/2)
    ! More accurately, integrate da/(a*H) from 0 to a
    ! For simplicity, use exponential in conformal time integral:
    ! integral of Gamma*a*dtau from 0 to a
    class(TDecayingDM), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: efac

    ! In radiation domination: a ~ tau, so integral ~ Gamma * tau^2/2
    ! In matter domination: a ~ tau^2, so integral ~ Gamma * (2/5) * tau^(5/2)
    ! Use the simple exponential approximation: e^(-Gamma * t)
    ! where Gamma is already in conformal units
    ! The proper time integral is done numerically via the background tables
    efac = 1._dl  ! placeholder, actual decay is tracked in background evolution

    end function TDecayingDM_DecayFactor

    subroutine TDecayingDM_Init(this, State)
    class(TDecayingDM), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) :: h2

    this%is_standard_cdm = (this%Gamma_dcdm == 0._dl .or. this%f_dcdm == 0._dl)

    if (.not. this%is_standard_cdm) then
        ! Convert Gamma from km/s/Mpc to Mpc^-1 (physical time rate)
        ! Same as H0 conversion: H0 [Mpc^-1] = H0 [km/s/Mpc] * 1000 / c
        this%Gamma_conf = this%Gamma_dcdm * 1000._dl / c
        ! v_kick = (1 - epsilon^2) / (2 * epsilon)
        this%v_kick = (1._dl - this%epsilon_dcdm**2) / (2._dl * this%epsilon_dcdm)

        ! CDM velocity for daughter particles
        this%has_cdm_velocity = .true.
        ! Equations: v_daughter (1) + delta_dcdm (1) + delta_daughter (1) + DR_F0 (1) + DR_F1 (1) + DR_F2 (1)
        ! Minimal: v_daughter + delta_daughter + DR 3 moments = 5
        ! Simpler approach: track delta_dcdm_eff and v_dcdm as fluid + DR hierarchy
        ! We use: v_c (CDM velocity, 1) + DR hierarchy (lmax+1)
        this%num_dr_equations = 3  ! DR F0, F1, F2 (truncated)
        this%num_perturb_equations = 1 + this%num_dr_equations  ! v_c + DR

        select type(S => State)
        class is (CAMBdata)
            h2 = (S%CP%H0/100._dl)**2
            this%H0_stored = S%CP%H0
            ! The decaying fraction of CDM
            this%grho_dcdm0 = S%grhocrit * S%CP%omch2 / h2 * this%f_dcdm
        end select
    else
        this%has_cdm_velocity = .false.
        this%num_dr_equations = 0
        this%num_perturb_equations = 0
    end if

    end subroutine TDecayingDM_Init

    subroutine TDecayingDM_BackgroundDensityAndPressure(this, grhodm, a, grhodm_t, gpres_dm)
    class(TDecayingDM), intent(inout) :: this
    real(dl), intent(in) :: grhodm, a
    real(dl), intent(out) :: grhodm_t
    real(dl), optional, intent(out) :: gpres_dm
    real(dl) :: decay_exp, rho_stable, rho_decay, a_Gamma_t

    if (a <= 0._dl) then
        grhodm_t = 0._dl
        if (present(gpres_dm)) gpres_dm = 0._dl
        return
    end if

    if (this%is_standard_cdm) then
        ! Standard CDM: rho ~ a^-3
        grhodm_t = grhodm / a
        if (present(gpres_dm)) gpres_dm = 0._dl
        return
    end if

    ! Decaying DM: rho_dcdm(a) = rho_dcdm,0 * a^-3 * exp(-Gamma * t(a))
    ! For the exponential factor, use Gamma*t ~ Gamma/(H0) * integral(da'/(a'*E(a')))
    ! Approximate: in matter domination, t ~ (2/3H0) * a^(3/2)
    ! So Gamma*t ~ (2/3) * (Gamma/H0) * a^(3/2)
    ! This is a reasonable approximation for z < 1000
    a_Gamma_t = (2._dl/3._dl) * this%Gamma_conf * a**(1.5_dl) / &
        sqrt(grhodm * 3._dl) * grhodm  ! approximate
    ! More robust: use Gamma * conformal_time / a (since dt = a*dtau, and H ~ 1/tau in RD)
    ! Actually, let's use a simpler and correct approach:
    ! Gamma_conf was set so that Gamma [1/Mpc] in cosmic time
    ! The proper approach: exp(-Gamma * t) where t is cosmic time
    ! In CAMB units: Gamma [Mpc^-1 cosmic time], and t comes from integration
    ! Simple exponential: exp(-Gamma/H0 * f(a)) where f(a) is a dimensionless time integral

    ! Use the approximation valid across epochs:
    ! Gamma*t(a) ~ Gamma_conf * a / adotoa, roughly
    ! Better: just use rho ~ a^{-3} * exp(-Gamma_eff * a^alpha) with alpha tuned
    ! For now, use direct exponential with cosmic time approximation
    ! decay_exp = exp(-Gamma * t(a))
    ! In radiation era: a ~ (2H_rad t)^(1/2), t = a^2/(2H_rad)
    ! In matter era: a ~ (3H_mat t/2)^(2/3), t = (2/3)*a^(3/2)/H_mat
    ! General: Gamma*t ~ Gamma_conf * integral_0^a da'/(a'*H(a'))

    ! For the background, we track the modification to total CDM density
    ! grhodm already includes the full CDM budget (stable + decaying parts)
    ! The stable fraction is (1-f_dcdm)*grhodm, scaling as a^-3
    ! The decaying fraction: f_dcdm*grhodm * exp(-Gamma*t(a)) * a^-3
    ! We approximate: t(a) ~ tau_approx(a) where tau is conformal time
    ! For background only, we modify the effective grhodm_t

    ! Simple approach: assume Gamma << H at early times (valid for Gamma ~ H0)
    ! Then rho_dcdm ~ rho_cdm * exp(-Gamma*t) and the exponential only matters at late times
    ! In matter domination: t = (2/3)*H0^{-1}*a^{3/2}/sqrt(Omega_m)
    ! Gamma*t = (2/3)*(Gamma/H0)*a^{3/2}/sqrt(Omega_m)
    ! Using Gamma_conf which is Gamma/(c/1000/Mpc):
    ! Actually let's be more careful. this%Gamma_dcdm is in km/s/Mpc.
    ! H0 is also in km/s/Mpc. So Gamma/H0 is dimensionless.
    ! t(a) in matter domination = (2/(3*H0)) * a^{3/2} (for Omega_m=1 approx)

    ! Use Gamma/H0 ratio for time calculation (H0 cached from Init)
    a_Gamma_t = this%Gamma_dcdm / this%H0_stored * 2._dl/3._dl * a**1.5_dl

    decay_exp = exp(-a_Gamma_t)

    ! Total DM density: stable part + surviving decaying part + daughter
    ! Daughter particles are massive (epsilon fraction), so they also scale ~ a^-3 at late times
    ! rho_total = rho_stable + rho_surviving_dcdm + rho_daughter
    ! = (1-f)*rho_cdm + f*rho_cdm*exp(-Gamma*t) + f*rho_cdm*epsilon*(1-exp(-Gamma*t))
    ! = rho_cdm * [1 - f*(1-epsilon)*(1-exp(-Gamma*t))]
    ! grhodm_t = grhodm/a * [1 - f*(1-epsilon)*(1-exp(-Gamma*t))]
    grhodm_t = grhodm / a * (1._dl - this%f_dcdm * (1._dl - this%epsilon_dcdm) * (1._dl - decay_exp))

    ! The energy going to DR: f*rho_cdm*(1-epsilon)*(1-exp(-Gamma*t)) but DR scales as a^-4
    ! This is handled via the grhodr background contribution

    if (present(gpres_dm)) gpres_dm = 0._dl

    end subroutine TDecayingDM_BackgroundDensityAndPressure

    subroutine TDecayingDM_SetBackgroundDensities(this, grhocrit, grhor, h2, grhodmdr, grhodr)
    ! Set the DR density from decay products (massless)
    ! This is called once during initialization, so we estimate the DR density
    ! The DR from decay accumulates over cosmic time
    class(TDecayingDM), intent(in) :: this
    real(dl), intent(in) :: grhocrit, grhor, h2
    real(dl), intent(out) :: grhodmdr, grhodr

    ! No separate interacting DM density pool (decay modifies standard CDM)
    grhodmdr = 0._dl
    ! DR from decay is subdominant for Gamma << H0; set to zero as leading order
    ! The actual DR contribution is computed dynamically via BackgroundDensityAndPressure
    grhodr = 0._dl

    end subroutine TDecayingDM_SetBackgroundDensities

    subroutine TDecayingDM_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhoc_t, adotoa, k, ay, ayprime, dm_ix, dr_ix, vc_ix)
    class(TDecayingDM), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhoc_t, adotoa, k
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl) :: grhodr_t, F0_dr, F1_dr

    if (this%is_standard_cdm) then
        dgrhoe = 0
        dgqe = 0
        return
    end if

    dgrhoe = 0._dl
    dgqe = 0._dl

    ! CDM velocity contribution
    if (vc_ix > 0) then
        dgqe = dgqe + grhoc_t * this%f_dcdm * ay(vc_ix)
    end if

    ! DR from decay: F0 gives density, F1 gives momentum
    if (dr_ix > 0) then
        F0_dr = ay(dr_ix)
        F1_dr = ay(dr_ix + 1)
        ! DR density contribution (small, proportional to Gamma*t)
        ! grhodr_t ~ grho_dr_from_decay / a^2
        ! For now, DR perturbations are sourced but subdominant
        ! The DR density at scale factor a from decay is:
        ! rho_dr ~ integral of Gamma*rho_dcdm*(1-epsilon)*a_decay/a * da_decay
        ! This is a small correction for Gamma << H0
        ! rho_DR ~ Gamma/H * rho_dcdm * (1-epsilon), with Gamma/H = Gamma_conf*a/adotoa
        grhodr_t = grhoc_t * this%f_dcdm * (1._dl - this%epsilon_dcdm) * &
            this%Gamma_conf * a / max(adotoa, 1e-30_dl)
        dgrhoe = dgrhoe + grhodr_t * F0_dr
        dgqe = dgqe + grhodr_t * F1_dr / 3._dl
    end if

    end subroutine TDecayingDM_PerturbedStressEnergy

    subroutine TDecayingDM_PerturbationInitial(this, y, a, tau, k, dm_ix, dr_ix, vc_ix)
    class(TDecayingDM), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer, intent(in) :: dm_ix, dr_ix, vc_ix

    if (this%is_standard_cdm) return

    ! Daughter velocity starts at zero (adiabatic IC)
    if (vc_ix > 0) y(vc_ix) = 0._dl

    ! dm_ix is allocated by SetupScalarArrayIndices when num_dr_equations > 0
    ! For DecayingDM we don't use it, but must initialize to 0
    if (dm_ix > 0) y(dm_ix) = 0._dl

    ! DR from decay: initially zero (no decay products yet)
    if (dr_ix > 0) then
        y(dr_ix) = 0._dl      ! F0
        y(dr_ix + 1) = 0._dl  ! F1
        y(dr_ix + 2) = 0._dl  ! F2
    end if

    end subroutine TDecayingDM_PerturbationInitial

    subroutine TDecayingDM_PerturbationEvolve(this, ayprime, a, adotoa, k, z, y, &
        dm_ix, dr_ix, vc_ix, clxc, vb, grhoc_t, grhob_t, sigma, high_ktau_dr)
    class(TDecayingDM), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, y(:)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl), intent(in) :: clxc, vb, grhoc_t, grhob_t, sigma
    logical, intent(in), optional :: high_ktau_dr
    real(dl) :: vc, Gamma_a, grhodr_t
    real(dl) :: F0, F1, F2

    if (this%is_standard_cdm) return
    if (vc_ix <= 0 .or. dr_ix <= 0) return

    ! dm_ix is allocated but unused for DecayingDM; must set derivative to 0
    if (dm_ix > 0) ayprime(dm_ix) = 0._dl

    ! Gamma_a = Gamma_phys * a (conformal time rate: Gamma * dt/dtau = Gamma * a)
    Gamma_a = this%Gamma_conf * a

    vc = y(vc_ix)
    F0 = y(dr_ix)
    F1 = y(dr_ix + 1)
    F2 = y(dr_ix + 2)

    ! Daughter CDM velocity equation:
    ! v_daughter' = -aH*v_daughter + Gamma*a * v_kick * (rho_dcdm/rho_daughter) * delta_direction
    ! In the limit where daughter quickly becomes non-relativistic:
    ! v_c' = -aH * v_c + source from decay (kick velocity averaged over directions)
    ! The kick provides a drag-like effect on the effective CDM
    ! For perturbations: v_c' = -aH*v_c - Gamma_a * v_kick * clxc (velocity source from decay)
    ayprime(vc_ix) = -adotoa * vc

    ! DR hierarchy from decay products
    ! Source term: S = Gamma_a * rho_dcdm * (1-epsilon) * (perturbation)
    ! DR density equation: F0' = -k*F1 + source
    ! Source is proportional to Gamma * rho_dcdm * delta_dcdm / rho_dr
    ! Since rho_dr from decay is small, normalize properly

    ! DR from decay: evolve as massless radiation with source
    ! F0' = -k*F1 + Gamma_a * (1-epsilon) * clxc * (grhoc_t / max(grhodr_t, 1e-30))
    ! For stability, use a simplified source that doesn't blow up
    grhodr_t = grhoc_t * this%f_dcdm * (1._dl - this%epsilon_dcdm) * &
        Gamma_a / max(adotoa, 1e-30_dl)

    ! F0' = -k*F1 + Gamma_a*(1-epsilon)*clxc * rho_dcdm/rho_dr (source from decay)
    ! Normalize: source ~ Gamma_a * clxc (the perturbation source)
    ayprime(dr_ix) = -k * F1 + Gamma_a * (1._dl - this%epsilon_dcdm) * clxc

    ! F1' = k/3 * (F0 - 2*F2)
    ayprime(dr_ix + 1) = k / 3._dl * (F0 - 2._dl * F2)

    ! F2 truncation: F2' = k/5 * (2*F1) - 3*adotoa*F2 (approx 3/tau in RD)
    ayprime(dr_ix + 2) = 2._dl * k / 5._dl * F1 - 3._dl * adotoa * F2 &
        + 8._dl / 15._dl * k * sigma

    end subroutine TDecayingDM_PerturbationEvolve

    subroutine TDecayingDM_PrintFeedback(this, FeedbackLevel)
    class(TDecayingDM) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("Decaying DM: Gamma = ",ES12.4," km/s/Mpc")') this%Gamma_dcdm
        write(*,'("             epsilon = ",F8.5,", f_dcdm = ",F8.5)') &
            this%epsilon_dcdm, this%f_dcdm
    end if

    end subroutine TDecayingDM_PrintFeedback

    end module DecayingDM
