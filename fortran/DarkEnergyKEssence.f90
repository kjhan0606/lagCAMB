    module DarkEnergyKEssence
    use DarkEnergyInterface
    use DarkEnergyFluid
    use results
    use constants
    use classes
    use config, only: GlobalError, error_unsupported_params
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    ! Purely kinetic k-essence dark energy.
    !
    ! Lagrangian   P(X) = M^4 (-Xt + Xt^2),   Xt = X/M^4  (dimensionless)
    ! (overall amplitude M^4 factored out; it only fixes the DE budget).
    !
    ! rho = 2 X P_X - P  =>  rho = M^4 (3 Xt^2 - Xt)
    ! w(a)   = (Xt - 1)/(3 Xt - 1)
    ! cs2(a) = P_X/(P_X + 2 X P_XX) = (2 Xt - 1)/(6 Xt - 1)   (rest frame)
    !
    ! Shift-symmetry (pure kinetic) EoM integrates algebraically:
    !   (2 Xt - 1) sqrt(Xt) = C0 a^-3 ,   C0 = (2 x0 - 1) sqrt(x0)
    ! with x0 = Xt(a=1) the single free shape parameter (must be > 1/2 so
    ! that rho>0, cs2>0, w in (-1, 1/3]).  The amplitude M^4 cancels in
    ! f_de = rho(X(a))/rho(x0), so no shooting is needed: solve the cubic
    ! for Xt(a) on a log-a grid, fill the tabulated-w machinery (base class)
    ! and a companion cs2(a) spline used by the fluid perturbation equations.
    !
    ! Limits:  early Xt ~ a^-2 -> rho ~ a^-4, w -> 1/3, cs2 -> 1/3
    !          late  Xt -> 1/2+          -> w -> -1,  cs2 -> 0+
    !
    ! x0 <=> kes_x0 in the N-body solver scalar_de_commons.f90 (identical
    ! background convention and algebraic solve).

    type, extends(TDarkEnergyFluid) :: TKEssence
        real(dl) :: x0 = 0.500001_dl   ! viable near-Lambda default; Xt today must be > 1/2
        integer  :: n_tab = 4096       ! background table size (log-a grid; matches lagRamses)
        real(dl) :: cs2_0 = 0._dl      ! cs2 today (diagnostic, set in Init)
        Type(TCubicSpline) :: cs2ofa   ! rest-frame cs2 as function of log(a)
    contains
    procedure :: ReadParams => TKEssence_ReadParams
    procedure :: Validate => TKEssence_Validate
    procedure, nopass :: PythonClass => TKEssence_PythonClass
    procedure, nopass :: SelfPointer => TKEssence_SelfPointer
    procedure :: Init => TKEssence_Init
    procedure :: PerturbationEvolve => TKEssence_PerturbationEvolve
    procedure :: SetKEssenceParams => TKEssence_SetKEssenceParams
    end type TKEssence

    contains

    subroutine TKEssence_Validate(this, OK)
    class(TKEssence), intent(in) :: this
    logical, intent(inout) :: OK

    call this%TDarkEnergyFluid%Validate(OK)
    if (.not. ieee_is_finite(this%x0) .or. this%x0 <= 0.5_dl .or. this%x0 > 0.50001_dl) then
        OK = .false.
        write(*,*) 'KEssence x0 must satisfy 0.5 < x0 <= 0.50001.'
    end if

    end subroutine TKEssence_Validate

    subroutine TKEssence_SetKEssenceParams(this, x0in)
    class(TKEssence), intent(inout) :: this
    real(dl), intent(in) :: x0in

    this%x0 = x0in

    end subroutine TKEssence_SetKEssenceParams


    subroutine TKEssence_ReadParams(this, Ini)
    use IniObjects
    class(TKEssence) :: this
    class(TIniFile), intent(in) :: Ini

    ! k-essence builds its own w(a)/cs2(a) tables in Init; only x0 is read.
    this%x0 = Ini%Read_Double('kessence_x0', 0.500001d0)

    end subroutine TKEssence_ReadParams


    function TKEssence_PythonClass()
    character(LEN=:), allocatable :: TKEssence_PythonClass

    TKEssence_PythonClass = 'KEssence'

    end function TKEssence_PythonClass


    subroutine TKEssence_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TKEssence), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TKEssence_SelfPointer


    subroutine TKEssence_Init(this, State)
    use classes
    class(TKEssence), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl), allocatable :: a(:), w(:), cs2(:)
    real(dl) :: c0, x, targ, g, gp, dxs, lna_min, dlna, lna
    integer  :: n, i, it

    if (.not. ieee_is_finite(this%x0) .or. this%x0 <= 0.5_dl .or. this%x0 > 0.50001_dl) then
        call GlobalError('KEssence x0 must satisfy 0.5 < x0 <= 0.50001.', error_unsupported_params)
        return
    end if

    n = this%n_tab
    allocate(a(n), w(n), cs2(n))

    ! log-a grid, amin = 1e-9, forcing a(n)=1 exactly (SetwTable needs a end=1)
    lna_min = log(1.d-9)
    dlna = -lna_min/dble(n-1)
    do i = 1, n
        a(i) = exp(lna_min + dble(i-1)*dlna)
    end do
    a(n) = 1.d0

    ! Algebraic background: (2X-1)sqrt(X) = C0 a^-3, C0 = (2 x0-1)sqrt(x0).
    ! March from a=1 (X=x0) to early times; X grows monotonically as a drops.
    ! Newton with damped step (mirrors scalar_de_commons.f90 kessence_tabulate).
    c0 = (2.d0*this%x0 - 1.d0)*sqrt(this%x0)
    x = this%x0
    do i = n, 1, -1
        lna = log(a(i))
        targ = c0*exp(-3.d0*lna)
        x = max(x, 0.5d0 + 0.25d0*(targ/2.d0)**(2.d0/3.d0))
        if (targ > 10.d0) x = max(x, (targ/2.d0)**(2.d0/3.d0))
        do it = 1, 100
            g  = (2.d0*x - 1.d0)*sqrt(x) - targ
            gp = 2.d0*sqrt(x) + (2.d0*x - 1.d0)/(2.d0*sqrt(x))
            dxs = g/gp
            if (x - dxs <= 0.5d0) dxs = 0.5d0*(x - 0.5d0)  ! damped step
            x = x - dxs
            if (abs(dxs) < 1.d-14*max(x, 1.d0)) exit
        end do
        w(i)   = (x - 1.d0)/(3.d0*x - 1.d0)
        cs2(i) = (2.d0*x - 1.d0)/(6.d0*x - 1.d0)
    end do

    ! Fill base-class tabulated-w machinery (equation_of_state + logdensity
    ! splines, sets use_tabulated_w=.true., w_lam/wa to today's values).
    call this%SetwTable(a, w, n)

    ! Companion rest-frame cs2(a) spline (same log-a abscissa as w table).
    call this%cs2ofa%Init(log(a), cs2)
    this%cs2_0   = cs2(n)
    this%cs2_lam = cs2(n)   ! today's cs2 (for any halofit/diagnostic use)

    ! Parent fluid Init: checks (no) w-crossing, sets num_perturb_equations.
    call this%TDarkEnergyFluid%Init(State)

    deallocate(a, w, cs2)

    end subroutine TKEssence_Init


    subroutine TKEssence_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TKEssence), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) :: Hv3_over_k, loga, cs2

    ! time-varying rest-frame sound speed from the k-essence table
    loga = log(a)
    if (loga <= this%cs2ofa%Xmin_interp) then
        cs2 = this%cs2ofa%F(1)
    else if (loga >= this%cs2ofa%Xmax_interp) then
        cs2 = this%cs2ofa%F(this%cs2ofa%n)
    else
        cs2 = this%cs2ofa%Value(loga)
    end if

    Hv3_over_k = 3*adotoa* y(w_ix + 1) / k
    !density perturbation (Fang-Hu-Lewis rest-frame fluid, cs2 -> cs2(a))
    ayprime(w_ix) = -3 * adotoa * (cs2 - w) * (y(w_ix) + (1 + w) * Hv3_over_k) &
        - (1 + w) * k * y(w_ix + 1) - (1 + w) * k * z
    !account for derivatives of w (use_tabulated_w is always true here)
    if (loga > this%equation_of_state%Xmin_interp .and. &
        loga < this%equation_of_state%Xmax_interp) then
        ayprime(w_ix) = ayprime(w_ix) - adotoa*this%equation_of_state%Derivative(loga)*Hv3_over_k
    end if
    !velocity
    if (abs(w+1) > 1e-6) then
        ayprime(w_ix + 1) = -adotoa * (1 - 3 * cs2) * y(w_ix + 1) + &
            k * cs2 * y(w_ix) / (1 + w)
    else
        ayprime(w_ix + 1) = 0
    end if

    end subroutine TKEssence_PerturbationEvolve

    end module DarkEnergyKEssence
