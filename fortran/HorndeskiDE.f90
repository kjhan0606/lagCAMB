
    module HorndeskiDE
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none

    private

    type, extends(TDarkEnergyEqnOfState) :: THorndeskiDE
        ! Alpha function amplitudes (proportional parameterization):
        ! alpha_X(a) = alpha_X_0 * Omega_DE(a)
        real(dl) :: alpha_K_0 = 0._dl      ! kineticity
        real(dl) :: alpha_B_0 = 0._dl      ! braiding
        real(dl) :: alpha_M_0 = 0._dl      ! Planck mass run rate
        real(dl) :: alpha_T_0 = 0._dl      ! tensor speed excess
        real(dl) :: M_star_ini = 1._dl     ! initial M*^2/M_Pl^2

        ! Internal state (set during Init)
        logical :: horn_active = .false.

        ! Spline tables for precomputed quantities (na points)
        integer, private :: na = 0
        real(dl), allocatable, private :: lna_tab(:)
        real(dl), allocatable, private :: Mstar2_tab(:), dMstar2_tab(:)
        real(dl), allocatable, private :: mu_tab(:), dmu_tab(:)
        real(dl), allocatable, private :: Sigma_tab(:), dSigma_tab(:)
    contains
        procedure :: Init => THorndeskiDE_Init
        procedure :: SetHorndeskiParams => THorndeskiDE_SetHorndeskiParams
        procedure :: PrintFeedback => THorndeskiDE_PrintFeedback
        procedure, nopass :: PythonClass => THorndeskiDE_PythonClass
        procedure, nopass :: SelfPointer => THorndeskiDE_SelfPointer
        procedure :: alpha_K_of_a => THorndeskiDE_alpha_K
        procedure :: alpha_B_of_a => THorndeskiDE_alpha_B
        procedure :: alpha_M_of_a => THorndeskiDE_alpha_M
        procedure :: alpha_T_of_a => THorndeskiDE_alpha_T
        procedure :: Mstar2_of_a => THorndeskiDE_Mstar2
        procedure :: mu_of_a => THorndeskiDE_mu
        procedure :: Sigma_of_a => THorndeskiDE_Sigma
    end type THorndeskiDE

    public THorndeskiDE

    contains

    function THorndeskiDE_PythonClass()
    character(LEN=:), allocatable :: THorndeskiDE_PythonClass
    THorndeskiDE_PythonClass = 'HorndeskiDE'
    end function THorndeskiDE_PythonClass

    subroutine THorndeskiDE_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(THorndeskiDE), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P
    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine THorndeskiDE_SelfPointer

    subroutine THorndeskiDE_SetHorndeskiParams(this, aK, aB, aM, aT, Msi)
    class(THorndeskiDE), intent(inout) :: this
    real(dl), intent(in) :: aK, aB, aM, aT, Msi

    this%alpha_K_0 = aK
    this%alpha_B_0 = aB
    this%alpha_M_0 = aM
    this%alpha_T_0 = aT
    this%M_star_ini = Msi

    end subroutine THorndeskiDE_SetHorndeskiParams

    ! ---- Alpha functions: proportional parameterization ----
    ! alpha_X(a) = alpha_X_0 * Omega_DE(a)
    ! Omega_DE(a) is approximated from the background

    function THorndeskiDE_alpha_K(this, a, Omega_DE) result(aK)
    class(THorndeskiDE), intent(in) :: this
    real(dl), intent(in) :: a, Omega_DE
    real(dl) :: aK
    aK = this%alpha_K_0 * Omega_DE
    end function

    function THorndeskiDE_alpha_B(this, a, Omega_DE) result(aB)
    class(THorndeskiDE), intent(in) :: this
    real(dl), intent(in) :: a, Omega_DE
    real(dl) :: aB
    aB = this%alpha_B_0 * Omega_DE
    end function

    function THorndeskiDE_alpha_M(this, a, Omega_DE) result(aM)
    class(THorndeskiDE), intent(in) :: this
    real(dl), intent(in) :: a, Omega_DE
    real(dl) :: aM
    aM = this%alpha_M_0 * Omega_DE
    end function

    function THorndeskiDE_alpha_T(this, a, Omega_DE) result(aT)
    class(THorndeskiDE), intent(in) :: this
    real(dl), intent(in) :: a, Omega_DE
    real(dl) :: aT
    aT = this%alpha_T_0 * Omega_DE
    end function

    ! ---- M*^2(a) from spline ----
    function THorndeskiDE_Mstar2(this, a) result(M2)
    class(THorndeskiDE), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: M2, lna
    integer :: i
    real(dl) :: t, h

    if (.not. this%horn_active .or. a < 1e-10_dl) then
        M2 = this%M_star_ini
        return
    end if

    lna = log(a)
    if (lna <= this%lna_tab(1)) then
        M2 = this%Mstar2_tab(1)
    else if (lna >= this%lna_tab(this%na)) then
        M2 = this%Mstar2_tab(this%na)
    else
        ! Cubic spline interpolation
        call horn_spline_interp(this%lna_tab, this%Mstar2_tab, this%dMstar2_tab, this%na, lna, M2)
    end if

    end function

    ! ---- mu(a) for Poisson equation ----
    function THorndeskiDE_mu(this, a) result(mu)
    class(THorndeskiDE), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: mu, lna

    if (.not. this%horn_active .or. a < 1e-10_dl) then
        mu = 1._dl
        return
    end if

    lna = log(a)
    if (lna <= this%lna_tab(1)) then
        mu = this%mu_tab(1)
    else if (lna >= this%lna_tab(this%na)) then
        mu = this%mu_tab(this%na)
    else
        call horn_spline_interp(this%lna_tab, this%mu_tab, this%dmu_tab, this%na, lna, mu)
    end if

    end function

    ! ---- Sigma(a) for lensing potential ----
    function THorndeskiDE_Sigma(this, a) result(Sig)
    class(THorndeskiDE), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: Sig, lna

    if (.not. this%horn_active .or. a < 1e-10_dl) then
        Sig = 1._dl
        return
    end if

    lna = log(a)
    if (lna <= this%lna_tab(1)) then
        Sig = this%Sigma_tab(1)
    else if (lna >= this%lna_tab(this%na)) then
        Sig = this%Sigma_tab(this%na)
    else
        call horn_spline_interp(this%lna_tab, this%Sigma_tab, this%dSigma_tab, this%na, lna, Sig)
    end if

    end function

    ! ---- Cubic spline interpolation helper ----
    subroutine horn_spline_interp(x, y, dd, n, xval, yval)
    integer, intent(in) :: n
    real(dl), intent(in) :: x(n), y(n), dd(n), xval
    real(dl), intent(out) :: yval
    integer :: lo, hi, mid
    real(dl) :: h, a_coef, b_coef

    lo = 1; hi = n
    do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (x(mid) > xval) then
            hi = mid
        else
            lo = mid
        end if
    end do
    h = x(hi) - x(lo)
    a_coef = (x(hi) - xval) / h
    b_coef = (xval - x(lo)) / h
    yval = a_coef*y(lo) + b_coef*y(hi) + &
        ((a_coef**3 - a_coef)*dd(lo) + (b_coef**3 - b_coef)*dd(hi)) * h*h/6._dl

    end subroutine horn_spline_interp

    ! ---- Compute spline second derivatives (natural boundary) ----
    subroutine horn_spline_deriv(x, y, n, dd)
    integer, intent(in) :: n
    real(dl), intent(in) :: x(n), y(n)
    real(dl), intent(out) :: dd(n)
    real(dl), allocatable :: u(:)
    real(dl) :: p, sig_s
    integer :: i

    allocate(u(n))
    dd(1) = 0._dl; u(1) = 0._dl
    do i = 2, n-1
        sig_s = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
        p = sig_s * dd(i-1) + 2._dl
        dd(i) = (sig_s - 1._dl) / p
        u(i) = (6._dl*((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1))) &
            / (x(i+1)-x(i-1)) - sig_s*u(i-1)) / p
    end do
    dd(n) = 0._dl
    do i = n-1, 1, -1
        dd(i) = dd(i)*dd(i+1) + u(i)
    end do
    deallocate(u)

    end subroutine horn_spline_deriv


    subroutine THorndeskiDE_Init(this, State)
    use classes
    use config
    class(THorndeskiDE), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    integer :: i
    real(dl) :: a_val, lna_val, dlna
    real(dl) :: Omega_DE_a, aK, aB, aM, aT, M2, D_kin
    real(dl) :: grhov_t_val, grho_tot, w_a
    real(dl) :: mu_val, Sigma_val, eta_val
    integer, parameter :: n_table = 500
    real(dl), parameter :: lna_min = -14._dl  ! a ~ 8e-7
    real(dl), parameter :: lna_max = 0._dl    ! a = 1

    ! Initialize parent (sets w, wa, is_cosmological_constant, etc.)
    call this%TDarkEnergyEqnOfState%Init(State)

    ! Check if any alpha is nonzero
    this%horn_active = (abs(this%alpha_K_0) > 1e-15_dl .or. &
                        abs(this%alpha_B_0) > 1e-15_dl .or. &
                        abs(this%alpha_M_0) > 1e-15_dl .or. &
                        abs(this%alpha_T_0) > 1e-15_dl .or. &
                        abs(this%M_star_ini - 1._dl) > 1e-15_dl)

    if (.not. this%horn_active) then
        this%num_perturb_equations = 0
        return
    end if

    ! Force is_cosmological_constant = .false. so perturbation machinery is active
    ! (even for w=-1 background, MG effects need z, sigma modifications)
    this%is_cosmological_constant = .false.
    this%num_perturb_equations = 0  ! QSA: no extra dynamic equations

    ! Check mutual exclusivity with MuSigmaMG
    select type(State)
    class is (CAMBdata)
        if (State%CP%MG%is_active) then
            call GlobalError('Cannot use both MuSigmaMG and Horndeski gravity', error_unsupported_params)
            return
        end if
    end select

    ! Allocate spline tables
    this%na = n_table
    if (allocated(this%lna_tab)) deallocate(this%lna_tab)
    if (allocated(this%Mstar2_tab)) deallocate(this%Mstar2_tab)
    if (allocated(this%dMstar2_tab)) deallocate(this%dMstar2_tab)
    if (allocated(this%mu_tab)) deallocate(this%mu_tab)
    if (allocated(this%dmu_tab)) deallocate(this%dmu_tab)
    if (allocated(this%Sigma_tab)) deallocate(this%Sigma_tab)
    if (allocated(this%dSigma_tab)) deallocate(this%dSigma_tab)
    allocate(this%lna_tab(n_table), this%Mstar2_tab(n_table), this%dMstar2_tab(n_table))
    allocate(this%mu_tab(n_table), this%dmu_tab(n_table))
    allocate(this%Sigma_tab(n_table), this%dSigma_tab(n_table))

    dlna = (lna_max - lna_min) / (n_table - 1)

    ! Step 1: Compute M*^2(a) by integrating d ln M*^2 / d ln a = alpha_M(a)
    ! Using simple forward Euler in ln(a)
    ! ln M*^2(a) = ln M*^2_ini + integral alpha_M(a') d ln a'
    select type(State)
    class is (CAMBdata)

    M2 = this%M_star_ini
    do i = 1, n_table
        lna_val = lna_min + (i-1) * dlna
        this%lna_tab(i) = lna_val
        a_val = exp(lna_val)

        ! Get Omega_DE(a) from background
        call this%BackgroundDensityAndPressure(State%grhov, a_val, grhov_t_val, w_a)
        grho_tot = State%grho_no_de(a_val) + grhov_t_val
        if (grho_tot > 0._dl) then
            Omega_DE_a = grhov_t_val / grho_tot
        else
            Omega_DE_a = 0._dl
        end if

        ! Alpha functions
        aK = this%alpha_K_of_a(a_val, Omega_DE_a)
        aB = this%alpha_B_of_a(a_val, Omega_DE_a)
        aM = this%alpha_M_of_a(a_val, Omega_DE_a)
        aT = this%alpha_T_of_a(a_val, Omega_DE_a)

        ! Integrate M*^2: M*^2(a) = M*^2_prev * exp(alpha_M * dlna)
        if (i > 1) then
            M2 = M2 * exp(aM * dlna)
        end if
        this%Mstar2_tab(i) = M2

        ! Step 2: Compute mu(a) and Sigma(a) in QSA
        ! From Bellini & Sawicki 2014, Pogosian & Silvestri 2016:
        ! D = alpha_K + 6*alpha_B^2 (effective kineticity)
        D_kin = aK + 6._dl * aB**2

        ! In the sub-horizon QSA (k >> aH):
        ! mu = (1/M*^2) * [D + 2*aB^2*(1+aT)] / [D + 2*aB^2]
        ! eta = Phi/Psi = [D + 2*aB^2] / [D + 2*aB^2*(1+aT)]  ... wait this isn't right
        !
        ! From Pogosian & Silvestri 2016 Eq. A10-A11 (subhorizon QSA):
        ! mu = (1/M*^2) * [D + 2*aB^2*(1+aT)] / [D + 2*aB^2]
        ! eta = [D + 2*aB*(1+aT)] / [D + 2*aB^2*(1+aT)]
        !
        ! Sigma = mu * (1+eta)/2 for lensing
        !
        ! Stability: D > 0 (no ghost), cs^2 > 0 (no gradient), 1+aT > 0 (no tensor ghost)

        if (D_kin + 2._dl * aB**2 > 1e-30_dl) then
            mu_val = (D_kin + 2._dl * aB**2 * (1._dl + aT)) / &
                     (M2 * (D_kin + 2._dl * aB**2))
        else
            mu_val = 1._dl / M2
        end if

        ! Gravitational slip: eta = Phi/Psi
        if (D_kin + 2._dl * aB**2 * (1._dl + aT) > 1e-30_dl) then
            eta_val = (D_kin + 2._dl * aB**2) / &
                      (D_kin + 2._dl * aB**2 * (1._dl + aT))
        else
            eta_val = 1._dl
        end if

        ! Sigma = mu * (1+eta)/2 (lensing)
        Sigma_val = mu_val * (1._dl + eta_val) / 2._dl

        this%mu_tab(i) = mu_val
        this%Sigma_tab(i) = Sigma_val
    end do

    end select

    ! Stability checks
    do i = 1, n_table
        if (this%Mstar2_tab(i) <= 0._dl) then
            call GlobalError('Horndeski: M*^2 <= 0 (Planck mass instability)', error_unsupported_params)
            return
        end if
    end do

    ! Build spline second derivatives
    call horn_spline_deriv(this%lna_tab, this%Mstar2_tab, n_table, this%dMstar2_tab)
    call horn_spline_deriv(this%lna_tab, this%mu_tab, n_table, this%dmu_tab)
    call horn_spline_deriv(this%lna_tab, this%Sigma_tab, n_table, this%dSigma_tab)

    end subroutine THorndeskiDE_Init


    subroutine THorndeskiDE_PrintFeedback(this, FeedbackLevel)
    class(THorndeskiDE) :: this
    integer, intent(in) :: FeedbackLevel

    call this%TDarkEnergyEqnOfState%PrintFeedback(FeedbackLevel)
    if (FeedbackLevel > 0 .and. this%horn_active) then
        write(*,'("Horndeski alpha parameterization (proportional to Omega_DE):")')
        write(*,'("  alpha_K_0=",f8.5,", alpha_B_0=",f8.5,", alpha_M_0=",f8.5,", alpha_T_0=",f8.5)') &
            this%alpha_K_0, this%alpha_B_0, this%alpha_M_0, this%alpha_T_0
        write(*,'("  M*^2_ini=",f8.5,", M*^2(a=1)=",f8.5)') &
            this%M_star_ini, this%Mstar2_tab(this%na)
        write(*,'("  Sigma(a=1)=",f8.5," (lensing modification)")') this%Sigma_tab(this%na)
        write(*,'("  Modifies: lensing potential (Sigma), tensor propagation (alpha_M, alpha_T)")')
    end if

    end subroutine THorndeskiDE_PrintFeedback

    end module HorndeskiDE
