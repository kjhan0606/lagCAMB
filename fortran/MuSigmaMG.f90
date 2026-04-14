    module MuSigmaMG
    use precision
    use classes
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

    public TMuSigmaMG, MG_params

    contains

    subroutine TMuSigmaMG_Init(this)
    class(TMuSigmaMG), intent(inout) :: this

    this%is_active = (this%mu_0 /= 0._dl .or. this%sigma_0 /= 0._dl)

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

    subroutine TMuSigmaMG_PrintFeedback(this, FeedbackLevel)
    class(TMuSigmaMG) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0 .and. this%is_active) then
        write(*,'("Modified Gravity (mu-Sigma):")')
        write(*,'("  mu_0 = ",F10.5,", Sigma_0 = ",F10.5)') this%mu_0, this%sigma_0
        if (this%lambda_mu > 0) write(*,'("  lambda_mu = ",ES10.3," Mpc")') this%lambda_mu
        if (this%lambda_sigma > 0) write(*,'("  lambda_sigma = ",ES10.3," Mpc")') this%lambda_sigma
    end if

    end subroutine TMuSigmaMG_PrintFeedback

    end module MuSigmaMG
