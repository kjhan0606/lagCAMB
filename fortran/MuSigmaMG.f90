    module MuSigmaMG
    use precision
    use classes
    use Interpolation, only : TCubicSpline
    implicit none

    private

    ! Phenomenological Modified Gravity: mu-Sigma parameterization
    !
    ! Modifies the Poisson and lensing potential equations:
    !   k^2 Psi = -4*pi*G * mu(a,k) * a^2 * rho * Delta    (Poisson)
    !   k^2 (Phi + Psi) = -8*pi*G * Sigma(a,k) * a^2 * rho * Delta  (Lensing)
    !
    ! In GR: mu = Sigma = 1
    !
    ! Parameterization (following MGCAMB, Zhao+ 2009, Planck 2018 MG):
    !   mu(a,k) = 1 + mu_0 * Omega_DE(a) * [1 + c1*(lambda_mu*k/a)^2 / (1 + (lambda_mu*k/a)^2)]
    !   Sigma(a,k) = 1 + Sigma_0 * Omega_DE(a) * [1 + c2*(lambda_Sigma*k/a)^2 / (1 + (lambda_Sigma*k/a)^2)]
    !
    ! For scale-independent case (lambda -> 0):
    !   mu(a) = 1 + mu_0 * Omega_DE(a)
    !   Sigma(a) = 1 + Sigma_0 * Omega_DE(a)
    !
    ! TABULATED mode (this file):
    !   In addition to the analytic parameterization above, a user may supply
    !   tabulated mu(a) and/or Sigma(a) via MuSigmaMG_SetMuTable / MuSigmaMG_SetSigmaTable
    !   (called from python MuSigmaMGParams.set_mu_a_table / set_sigma_a_table).
    !   When a table is present the corresponding mu(a) or Sigma(a) is obtained by
    !   cubic-spline interpolation in ln(a) (scale independent). Outside the table
    !   range the value is CLAMPED to the nearest table edge: below the smallest a
    !   it takes the a->0 (first-row) value; above a=1 (or the largest tabulated a)
    !   it takes the last-row value. A single-row table is treated as a constant.
    !   The two tables are independent: if only a mu table is set, Sigma still
    !   follows its analytic (1 + Sigma_0*Omega_DE) path, and vice versa.
    !
    !   The tables live as module-global state (a POD CAMB_Structure field cannot
    !   hold allocatables), mirroring how DarkEnergyEqnOfState%SetwTable builds a
    !   spline. As a consequence the tables are process-global: setting a table
    !   affects every subsequent get_results() until MuSigmaMG_ClearTables is
    !   called (set_mu_a_table/set_sigma_a_table with empty arrays, or clear_tables
    !   on the python side).
    !
    ! This is NOT a dark energy model per se; it modifies the gravitational sector.
    ! It is implemented as a separate module that equations.f90 calls to rescale potentials.
    !
    ! References: Zhao+ 2009, Simpson+ 2013, Planck 2015/2018 XIV, MGCAMB

    type :: TMuSigmaMG
        logical :: is_active = .false.
        real(dl) :: mu_0 = 0._dl         ! mu deviation amplitude
        real(dl) :: sigma_0 = 0._dl      ! Sigma deviation amplitude
        real(dl) :: lambda_mu = 0._dl    ! mu scale [Mpc] (0 = scale-independent)
        real(dl) :: lambda_sigma = 0._dl ! Sigma scale [Mpc] (0 = scale-independent)
        real(dl) :: c1 = 1._dl           ! mu k-dependence coefficient
        real(dl) :: c2 = 1._dl           ! Sigma k-dependence coefficient
        ! ---- Theory-specific QSA models (MGCAMB-style) ----
        ! All assume the LCDM background (standard for viable parameter space):
        !   1 = f(R) Hu-Sawicki:  mu = (1+4Q/3)/(1+Q), Q = k^2/(a^2 M^2(a)),
        !       M^2 = 1/(3 f_RR) = R^{n+2}/(3(n+1)|fR0| R0^{n+1}); Sigma = 1.
        !   2 = nDGP normal branch: mu = 1 + 1/(3 beta(a)),
        !       beta = 1 + 2 H rc (1 + Hdot/(3H^2)); Sigma = 1 (scale-independent).
        !   3 = Symmetron: mu = 1 for a<=a_ssb, else
        !       mu = 1 + 2 beta*(a)^2 k^2/(k^2 + a^2 m(a)^2),
        !       beta*(a) = beta_sym sqrt(1-(a_ssb/a)^3), m(a) = sqrt(1-(a_ssb/a)^3)/L_sym;
        !       Sigma = 1.
        ! Sigma=1 with the induced slip eta = 2/mu - 1 reproduces each theory's QSA
        ! slip exactly (f(R): (1+2Q/3)/(1+4Q/3); nDGP: (3beta-1)/(3beta+1)).
        integer :: model = 0             ! 0=phenomenological mu-Sigma above
        real(dl) :: fR0 = 0._dl          ! f(R): |f_R0| today (e.g. 1e-5 for F5)
        real(dl) :: fR_n = 1._dl         ! f(R): Hu-Sawicki index n
        real(dl) :: Omega_rc = 0._dl     ! nDGP: Omega_rc = 1/(4 rc^2 H0^2)
        real(dl) :: beta_sym = 0._dl     ! Symmetron: coupling beta_*
        real(dl) :: a_ssb = 0.5_dl       ! Symmetron: symmetry-breaking scale factor
        real(dl) :: L_sym = 1._dl        ! Symmetron: Compton range today [Mpc]
        ! Derived LCDM background constants, filled by Init (from results.f90)
        real(dl) :: bg_H0 = 0._dl        ! H0 [Mpc^-1]
        real(dl) :: bg_omm = 0._dl       ! Omega_m (baryons + CDM + massive nu)
        real(dl) :: bg_omv = 0._dl       ! Omega_Lambda (dark energy)
        real(dl) :: bg_omr = 0._dl       ! Omega_r (photons + massless nu)
    contains
    procedure :: Init => TMuSigmaMG_Init
    procedure :: mu_of_a_k => TMuSigmaMG_mu
    procedure :: Sigma_of_a_k => TMuSigmaMG_Sigma
    procedure :: eta_of_a_k => TMuSigmaMG_eta  ! eta = Phi/Psi = 2*Sigma/mu - 1
    procedure :: PrintFeedback => TMuSigmaMG_PrintFeedback
    end type TMuSigmaMG

    ! Global instance (set from CAMBparams)
    type(TMuSigmaMG), save, target :: MG_params

    ! ---- Tabulated-mode module-global state (see header notes) ----
    logical, save :: mg_mu_table_set = .false.
    logical, save :: mg_mu_is_const  = .false.
    real(dl), save :: mg_mu_const = 1._dl
    real(dl), save :: mg_mu_lna_min = 0._dl, mg_mu_lna_max = 0._dl
    real(dl), save :: mg_mu_edge_lo = 1._dl, mg_mu_edge_hi = 1._dl
    type(TCubicSpline), save :: mg_mu_spline

    logical, save :: mg_sigma_table_set = .false.
    logical, save :: mg_sigma_is_const  = .false.
    real(dl), save :: mg_sigma_const = 1._dl
    real(dl), save :: mg_sigma_lna_min = 0._dl, mg_sigma_lna_max = 0._dl
    real(dl), save :: mg_sigma_edge_lo = 1._dl, mg_sigma_edge_hi = 1._dl
    type(TCubicSpline), save :: mg_sigma_spline

    public TMuSigmaMG, MG_params
    public MuSigmaMG_SetMuTable, MuSigmaMG_SetSigmaTable, MuSigmaMG_ClearTables

    contains

    subroutine TMuSigmaMG_Init(this, H0_Mpc, omm, omv, omr)
    class(TMuSigmaMG), intent(inout) :: this
    real(dl), intent(in) :: H0_Mpc, omm, omv, omr

    this%bg_H0 = H0_Mpc
    this%bg_omm = omm
    this%bg_omv = omv
    this%bg_omr = omr

    select case (this%model)
    case (1)  ! f(R) Hu-Sawicki
        this%is_active = this%fR0 /= 0._dl
    case (2)  ! nDGP
        this%is_active = this%Omega_rc > 0._dl
    case (3)  ! Symmetron
        this%is_active = this%beta_sym /= 0._dl .and. this%L_sym > 0._dl
    case default
        ! Active if the analytic amplitudes are non-zero OR a table has been supplied.
        this%is_active = (this%mu_0 /= 0._dl .or. this%sigma_0 /= 0._dl &
            .or. mg_mu_table_set .or. mg_sigma_table_set)
    end select

    end subroutine TMuSigmaMG_Init

    function TMuSigmaMG_mu(this, a, k, Omega_DE_a) result(mu)
    ! mu(a,k): modification to Poisson equation
    ! Psi = mu * Psi_GR
    class(TMuSigmaMG), intent(in) :: this
    real(dl), intent(in) :: a, k, Omega_DE_a
    real(dl) :: mu
    real(dl) :: koa2, scale_factor
    real(dl) :: Q, lnQ, R_a, R0, E2, beta_dgp, s

    mu = 1._dl
    if (.not. this%is_active) return

    select case (this%model)
    case (1)
        ! f(R) Hu-Sawicki QSA on LCDM background:
        ! mu = (1+4Q/3)/(1+Q), Q = k^2/(a^2 M^2(a)),
        ! M^2(a) = 1/(3 f_RR) = R^{n+2}/(3(n+1)|fR0| R0^{n+1}),
        ! R(a) = 3 H0^2 (Omega_m/a^3 + 4 Omega_L) (radiation is traceless).
        ! Q spans many decades, so work in log space to avoid overflow.
        if (a <= 0._dl .or. k <= 0._dl) return
        R_a = 3*this%bg_H0**2 * (this%bg_omm/a**3 + 4*this%bg_omv)
        R0  = 3*this%bg_H0**2 * (this%bg_omm + 4*this%bg_omv)
        lnQ = 2*log(k/a) - log(R_a/(3*(this%fR_n+1)*abs(this%fR0))) &
            - (this%fR_n+1)*log(R_a/R0)
        if (lnQ > 30._dl) then
            mu = 4._dl/3._dl
        else if (lnQ > -30._dl) then
            Q = exp(lnQ)
            mu = (1 + 4*Q/3)/(1 + Q)
        end if
    case (2)
        ! nDGP normal branch on LCDM background: mu = 1 + 1/(3 beta(a)),
        ! beta = 1 + 2 H rc (1 + Hdot/(3H^2)), 2 H rc = sqrt(E2/Omega_rc).
        ! In (1 + Hdot/3H^2): radiation era -> 1/3, matter era -> 1/2, Lambda -> 1,
        ! so beta > 1 always and mu -> 1 at early times.
        if (a <= 0._dl) return
        E2 = this%bg_omm/a**3 + this%bg_omr/a**4 + this%bg_omv
        beta_dgp = 1 + sqrt(E2/this%Omega_rc) * &
            (1 - (1.5_dl*this%bg_omm/a**3 + 2*this%bg_omr/a**4)/(3*E2))
        mu = 1 + 1/(3*beta_dgp)
    case (3)
        ! Symmetron: symmetric (fully screened) before a_ssb, then
        ! mu = 1 + 2 beta*(a)^2 k^2/(k^2 + a^2 m(a)^2) with
        ! beta*(a) = beta_sym sqrt(s), m(a) = sqrt(s)/L_sym, s = 1-(a_ssb/a)^3.
        if (a <= this%a_ssb) return
        s = 1 - (this%a_ssb/a)**3
        mu = 1 + 2*this%beta_sym**2*s * k**2/(k**2 + a**2*s/this%L_sym**2)
    case default
        ! Tabulated mode takes precedence when a mu table is present (scale independent).
        if (mg_mu_table_set) then
            mu = mg_eval_table(a, mg_mu_is_const, mg_mu_const, mg_mu_spline, &
                mg_mu_lna_min, mg_mu_lna_max, mg_mu_edge_lo, mg_mu_edge_hi)
            return
        end if

        mu = 1._dl + this%mu_0 * Omega_DE_a

        ! Scale dependence
        if (this%lambda_mu > 0._dl .and. a > 0._dl) then
            koa2 = (this%lambda_mu * k / a)**2
            scale_factor = this%c1 * koa2 / (1._dl + koa2)
            mu = 1._dl + this%mu_0 * Omega_DE_a * (1._dl + scale_factor)
        end if
    end select

    end function TMuSigmaMG_mu

    function TMuSigmaMG_Sigma(this, a, k, Omega_DE_a) result(Sig)
    ! Sigma(a,k): modification to lensing potential
    ! (Phi + Psi) = Sigma * (Phi + Psi)_GR
    class(TMuSigmaMG), intent(in) :: this
    real(dl), intent(in) :: a, k, Omega_DE_a
    real(dl) :: Sig
    real(dl) :: koa2, scale_factor

    if (.not. this%is_active) then
        Sig = 1._dl
        return
    end if

    ! f(R)/nDGP/Symmetron: the Weyl potential sourced by matter is unmodified
    ! in QSA (Sigma = 1 exactly); the slip eta = 2/mu - 1 then matches each
    ! theory's QSA gamma. Lensing still changes indirectly via modified growth.
    if (this%model /= 0) then
        Sig = 1._dl
        return
    end if

    ! Tabulated mode takes precedence when a Sigma table is present.
    if (mg_sigma_table_set) then
        Sig = mg_eval_table(a, mg_sigma_is_const, mg_sigma_const, mg_sigma_spline, &
            mg_sigma_lna_min, mg_sigma_lna_max, mg_sigma_edge_lo, mg_sigma_edge_hi)
        return
    end if

    Sig = 1._dl + this%sigma_0 * Omega_DE_a

    ! Scale dependence
    if (this%lambda_sigma > 0._dl .and. a > 0._dl) then
        koa2 = (this%lambda_sigma * k / a)**2
        scale_factor = this%c2 * koa2 / (1._dl + koa2)
        Sig = 1._dl + this%sigma_0 * Omega_DE_a * (1._dl + scale_factor)
    end if

    end function TMuSigmaMG_Sigma

    function TMuSigmaMG_eta(this, a, k, Omega_DE_a) result(eta)
    ! eta = Phi/Psi = 2*Sigma/mu - 1
    ! In GR: eta = 1 (no anisotropic stress from matter)
    class(TMuSigmaMG), intent(in) :: this
    real(dl), intent(in) :: a, k, Omega_DE_a
    real(dl) :: eta, mu_val, sig_val

    mu_val = this%mu_of_a_k(a, k, Omega_DE_a)
    sig_val = this%Sigma_of_a_k(a, k, Omega_DE_a)

    if (mu_val > 0._dl) then
        eta = 2._dl * sig_val / mu_val - 1._dl
    else
        eta = 1._dl
    end if

    end function TMuSigmaMG_eta

    ! ---------- Tabulated-mode helpers ----------

    function mg_eval_table(a, is_const, const_val, spline, lna_min, lna_max, &
        edge_lo, edge_hi) result(val)
    ! Interpolate a tabulated mu/Sigma at scale factor a. Cubic spline in ln(a);
    ! clamp to the nearest edge value outside [lna_min, lna_max]. a<=0 -> a->0 edge.
    real(dl), intent(in) :: a, const_val, lna_min, lna_max, edge_lo, edge_hi
    logical, intent(in) :: is_const
    type(TCubicSpline) :: spline   ! no intent: %Value may update cached spline state
    real(dl) :: val, lna

    if (is_const) then
        val = const_val
        return
    end if

    if (a <= 0._dl) then
        val = edge_lo
        return
    end if

    lna = log(a)
    if (lna <= lna_min) then
        val = edge_lo
    else if (lna >= lna_max) then
        val = edge_hi
    else
        val = spline%Value(lna)
    end if

    end function mg_eval_table

    subroutine MuSigmaMG_SetMuTable(a, mu, n)
    ! Install a tabulated mu(a). n==0 clears the table; n==1 is a constant mu.
    ! a must be strictly increasing and positive. Called from python via ctypes.
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n), mu(n)

    if (mg_mu_table_set .and. .not. mg_mu_is_const) call mg_mu_spline%Clear()
    mg_mu_table_set = .false.
    mg_mu_is_const = .false.

    if (n <= 0) return

    if (n == 1) then
        mg_mu_is_const = .true.
        mg_mu_const = mu(1)
    else
        call mg_mu_spline%Init(log(a), mu)
        mg_mu_lna_min = log(a(1))
        mg_mu_lna_max = log(a(n))
        mg_mu_edge_lo = mu(1)
        mg_mu_edge_hi = mu(n)
    end if
    mg_mu_table_set = .true.

    end subroutine MuSigmaMG_SetMuTable

    subroutine MuSigmaMG_SetSigmaTable(a, sig, n)
    ! Install a tabulated Sigma(a). n==0 clears; n==1 is constant.
    integer, intent(in) :: n
    real(dl), intent(in) :: a(n), sig(n)

    if (mg_sigma_table_set .and. .not. mg_sigma_is_const) call mg_sigma_spline%Clear()
    mg_sigma_table_set = .false.
    mg_sigma_is_const = .false.

    if (n <= 0) return

    if (n == 1) then
        mg_sigma_is_const = .true.
        mg_sigma_const = sig(1)
    else
        call mg_sigma_spline%Init(log(a), sig)
        mg_sigma_lna_min = log(a(1))
        mg_sigma_lna_max = log(a(n))
        mg_sigma_edge_lo = sig(1)
        mg_sigma_edge_hi = sig(n)
    end if
    mg_sigma_table_set = .true.

    end subroutine MuSigmaMG_SetSigmaTable

    subroutine MuSigmaMG_ClearTables()
    ! Drop both tables and return to the analytic path.
    if (mg_mu_table_set .and. .not. mg_mu_is_const) call mg_mu_spline%Clear()
    if (mg_sigma_table_set .and. .not. mg_sigma_is_const) call mg_sigma_spline%Clear()
    mg_mu_table_set = .false.
    mg_mu_is_const = .false.
    mg_sigma_table_set = .false.
    mg_sigma_is_const = .false.

    end subroutine MuSigmaMG_ClearTables

    subroutine TMuSigmaMG_PrintFeedback(this, FeedbackLevel)
    class(TMuSigmaMG) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0 .and. this%is_active) then
        select case (this%model)
        case (1)
            write(*,'("Modified Gravity: f(R) Hu-Sawicki (QSA)")')
            write(*,'("  |fR0| = ",ES10.3,"  n = ",F6.2)') abs(this%fR0), this%fR_n
            return
        case (2)
            write(*,'("Modified Gravity: nDGP normal branch (QSA)")')
            write(*,'("  Omega_rc = ",ES10.3)') this%Omega_rc
            return
        case (3)
            write(*,'("Modified Gravity: Symmetron (QSA)")')
            write(*,'("  a_ssb = ",F8.4,"  beta = ",F8.4,"  L = ",ES10.3," Mpc")') &
                this%a_ssb, this%beta_sym, this%L_sym
            return
        end select
        write(*,'("Modified Gravity (mu-Sigma):")')
        if (mg_mu_table_set) then
            write(*,'("  mu(a): tabulated")')
        else
            write(*,'("  mu_0 = ",F10.5)') this%mu_0
        end if
        if (mg_sigma_table_set) then
            write(*,'("  Sigma(a): tabulated")')
        else
            write(*,'("  Sigma_0 = ",F10.5)') this%sigma_0
        end if
        if (this%lambda_mu > 0) write(*,'("  lambda_mu = ",ES10.3," Mpc")') this%lambda_mu
        if (this%lambda_sigma > 0) write(*,'("  lambda_sigma = ",ES10.3," Mpc")') this%lambda_sigma
    end if

    end subroutine TMuSigmaMG_PrintFeedback

    end module MuSigmaMG
