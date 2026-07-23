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

    subroutine TMuSigmaMG_Init(this)
    class(TMuSigmaMG), intent(inout) :: this

    ! Active if the analytic amplitudes are non-zero OR a table has been supplied.
    this%is_active = (this%mu_0 /= 0._dl .or. this%sigma_0 /= 0._dl &
        .or. mg_mu_table_set .or. mg_sigma_table_set)

    end subroutine TMuSigmaMG_Init

    function TMuSigmaMG_mu(this, a, k, Omega_DE_a) result(mu)
    ! mu(a,k): modification to Poisson equation
    ! Psi = mu * Psi_GR
    class(TMuSigmaMG), intent(in) :: this
    real(dl), intent(in) :: a, k, Omega_DE_a
    real(dl) :: mu
    real(dl) :: koa2, scale_factor

    if (.not. this%is_active) then
        mu = 1._dl
        return
    end if

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
