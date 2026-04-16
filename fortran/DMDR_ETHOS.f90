    module DMDR_ETHOS
    use DarkMatterInteraction
    use results
    use constants
    use classes
    implicit none

    integer, parameter :: max_ethos_coeffs = 6

    ! DM-Dark Radiation (ETHOS) model
    ! DM coupled to a dark radiation sector via scattering
    ! DM: fluid equations (delta, v)
    ! DR: Boltzmann hierarchy (F_0, F_1, ..., F_lmax)
    type, extends(TDarkMatterModel) :: TDMDR_ETHOS
        real(dl) :: omdmdrh2 = 0._dl      ! interacting DM density Omega_dmdr h^2
        real(dl) :: N_dark = 0._dl         ! dark radiation effective number Delta_N_eff
        real(dl) :: a_dark(max_ethos_coeffs) = 0._dl  ! ETHOS opacity coefficients
        real(dl) :: alpha_l_dark = 1._dl   ! angular damping coefficient for l>=2
        integer  :: lmax_dr = 15           ! DR Boltzmann hierarchy truncation
        real(dl) :: cs2_dm = 0._dl         ! DM effective sound speed squared (usually ~0)
        ! Cached background densities (set in Init)
        real(dl) :: grhodmdr_bg = 0._dl    ! 8*pi*G*rho_dmdr,0 * a0^4 (CAMB convention)
        real(dl) :: grhodr_bg = 0._dl      ! 8*pi*G*rho_dr,0 * a0^4
    contains
    procedure :: ReadParams => TDMDR_ETHOS_ReadParams
    procedure, nopass :: PythonClass => TDMDR_ETHOS_PythonClass
    procedure, nopass :: SelfPointer => TDMDR_ETHOS_SelfPointer
    procedure :: Init => TDMDR_ETHOS_Init
    procedure :: PerturbedStressEnergy => TDMDR_ETHOS_PerturbedStressEnergy
    procedure :: PerturbationInitial => TDMDR_ETHOS_PerturbationInitial
    procedure :: PerturbationEvolve => TDMDR_ETHOS_PerturbationEvolve
    procedure :: PrintFeedback => TDMDR_ETHOS_PrintFeedback
    procedure :: Opacity_DMDR => TDMDR_ETHOS_Opacity
    procedure :: SetBackgroundDensities => TDMDR_ETHOS_SetBackgroundDensities
    end type TDMDR_ETHOS

    contains

    subroutine TDMDR_ETHOS_ReadParams(this, Ini)
    use IniObjects
    class(TDMDR_ETHOS) :: this
    class(TIniFile), intent(in) :: Ini
    integer :: i

    this%omdmdrh2 = Ini%Read_Double('omdmdrh2', 0._dl)
    this%N_dark = Ini%Read_Double('N_dark', 0._dl)
    this%lmax_dr = Ini%Read_Int('lmax_dr', 15)
    this%alpha_l_dark = Ini%Read_Double('alpha_l_dark', 1._dl)
    this%cs2_dm = Ini%Read_Double('cs2_dm', 0._dl)
    do i = 1, max_ethos_coeffs
        this%a_dark(i) = 0._dl
    end do
    ! Read individual coefficients if present
    call Ini%Read('a_dark_2', this%a_dark(1))
    call Ini%Read('a_dark_3', this%a_dark(2))
    call Ini%Read('a_dark_4', this%a_dark(3))
    call Ini%Read('a_dark_5', this%a_dark(4))
    call Ini%Read('a_dark_6', this%a_dark(5))

    end subroutine TDMDR_ETHOS_ReadParams

    function TDMDR_ETHOS_PythonClass()
    character(LEN=:), allocatable :: TDMDR_ETHOS_PythonClass
    TDMDR_ETHOS_PythonClass = 'DMDR_ETHOS'
    end function TDMDR_ETHOS_PythonClass

    subroutine TDMDR_ETHOS_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TDMDR_ETHOS), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TDMDR_ETHOS_SelfPointer

    subroutine TDMDR_ETHOS_Init(this, State)
    class(TDMDR_ETHOS), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) :: h2

    this%is_standard_cdm = (this%omdmdrh2 == 0._dl .and. this%N_dark == 0._dl)

    if (.not. this%is_standard_cdm) then
        ! DM velocity (1) + DR hierarchy (lmax_dr + 1 moments: F_0 to F_lmax)
        this%has_cdm_velocity = .true.
        this%num_dr_equations = this%lmax_dr + 1
        this%num_perturb_equations = 1 + this%num_dr_equations  ! v_dm + DR hierarchy

        ! Cache background densities (already set in results.f90 SetParams)
        select type(S => State)
        class is (CAMBdata)
            this%grhodmdr_bg = S%grhodmdr
            this%grhodr_bg = S%grhodr
        end select
    else
        this%has_cdm_velocity = .false.
        this%num_dr_equations = 0
        this%num_perturb_equations = 0
    end if

    end subroutine TDMDR_ETHOS_Init

    function TDMDR_ETHOS_Opacity(this, a) result(kappa)
    ! DM-DR opacity: kappa_DM-DR(T_DR)
    ! kappa = Sum_n a_n * (T_DR/T_0)^n
    ! where T_DR/T_0 = 1/a (assuming DR temperature scales as a^-1)
    class(TDMDR_ETHOS), intent(in) :: this
    real(dl), intent(in) :: a
    real(dl) :: kappa
    real(dl) :: T_ratio
    integer :: i

    T_ratio = 1._dl / a  ! T_DR / T_0

    kappa = 0._dl
    do i = 1, max_ethos_coeffs
        if (this%a_dark(i) /= 0._dl) then
            kappa = kappa + this%a_dark(i) * T_ratio**(i+1)  ! a_2 * T^2, a_3 * T^3, etc.
        end if
    end do

    end function TDMDR_ETHOS_Opacity

    subroutine TDMDR_ETHOS_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhoc_t, adotoa, k, ay, ayprime, dm_ix, dr_ix, vc_ix)
    ! Returns density and velocity perturbation contributions from the DM-DR sector
    class(TDMDR_ETHOS), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhoc_t, adotoa, k
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl) :: grhodmdr_t, grhodr_t, clx_dmdr, v_dm, F0_dr, F1_dr
    real(dl) :: a2, h2

    if (this%is_standard_cdm) then
        dgrhoe = 0
        dgqe = 0
        return
    end if

    a2 = a*a

    ! Background densities at this scale factor
    grhodmdr_t = this%grhodmdr_bg / a  ! 8*pi*G*rho_dmdr * a^2
    grhodr_t = this%grhodr_bg / a2     ! 8*pi*G*rho_dr * a^2

    ! DM density perturbation
    if (dm_ix > 0) then
        clx_dmdr = ay(dm_ix)
    else
        clx_dmdr = 0
    end if

    ! DR monopole and dipole
    if (dr_ix > 0) then
        F0_dr = ay(dr_ix)
        F1_dr = ay(dr_ix + 1)
    else
        F0_dr = 0
        F1_dr = 0
    end if

    dgrhoe = grhodmdr_t * clx_dmdr + grhodr_t * F0_dr

    v_dm = 0._dl
    if (vc_ix > 0) v_dm = ay(vc_ix)
    ! (rho+p)*v for DR: (4/3)*rho_dr * (F1/4)
    dgqe = grhodmdr_t * v_dm + grhodr_t * F1_dr / 3._dl

    end subroutine TDMDR_ETHOS_PerturbedStressEnergy

    subroutine TDMDR_ETHOS_PerturbationInitial(this, y, a, tau, k, dm_ix, dr_ix, vc_ix)
    class(TDMDR_ETHOS), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl) :: x2, x3

    if (this%is_standard_cdm) return

    x2 = (k*tau)**2
    x3 = (k*tau)**3

    ! Adiabatic initial conditions
    ! DM density: same as standard CDM
    if (dm_ix > 0) then
        y(dm_ix) = -x2/2._dl * 0.75_dl  ! delta_dmdr = (3/4)*delta_g = -(3/4)*(1/3)*k^2*tau^2
    end if

    ! DM velocity: starts at zero
    if (vc_ix > 0) then
        y(vc_ix) = 0._dl
    end if

    ! DR hierarchy: adiabatic, like photons/massless neutrinos
    if (dr_ix > 0) then
        y(dr_ix) = -2._dl/3._dl * x2     ! F_0 = -(2/3) k^2 tau^2 (like clxr)
        y(dr_ix + 1) = 2._dl/9._dl * x3  ! F_1 = (2/9) k^3 tau^3
        ! F_l = 0 for l >= 2 (already zero)
    end if

    end subroutine TDMDR_ETHOS_PerturbationInitial

    subroutine TDMDR_ETHOS_PerturbationEvolve(this, ayprime, a, adotoa, k, z, y, &
        dm_ix, dr_ix, vc_ix, clxc, vb, grhoc_t, grhob_t, sigma, high_ktau_dr)
    ! Evolve the coupled DM-DR system
    class(TDMDR_ETHOS), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, y(:)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl), intent(in) :: clxc, vb, grhoc_t, grhob_t, sigma
    logical, intent(in), optional :: high_ktau_dr
    real(dl) :: clx_dm, v_dm, kappa, Gamma_DR, Gamma_DR_tilde
    real(dl) :: grhodmdr_t, grhodr_t
    real(dl) :: denl_l, denl2_l
    integer :: l, ix
    real(dl) :: a2, cothxor
    logical :: use_ufa

    if (this%is_standard_cdm) return
    if (dm_ix <= 0 .or. dr_ix <= 0 .or. vc_ix <= 0) return

    a2 = a*a
    cothxor = adotoa / a  ! 1/tau approximation for flat

    grhodmdr_t = this%grhodmdr_bg / a
    grhodr_t = this%grhodr_bg / a2

    ! Get opacity
    kappa = this%Opacity_DMDR(a)

    ! Drag rates
    ! Gamma_DR = a * kappa * rho_DR / rho_DM (drag on DM from DR)
    if (grhodmdr_t > 0 .and. grhodr_t > 0) then
        Gamma_DR = a * kappa * grhodr_t / grhodmdr_t
        ! Gamma_tilde_DR = a * kappa * rho_DM / rho_DR * (4/3) (drag on DR from DM)
        Gamma_DR_tilde = a * kappa * grhodmdr_t / grhodr_t * (4._dl/3._dl)
    else
        Gamma_DR = 0._dl
        Gamma_DR_tilde = 0._dl
    end if

    ! DM perturbation variables
    clx_dm = y(dm_ix)
    v_dm = y(vc_ix)

    ! ============ DM equations (fluid) ============
    ! delta_dm' = -k*v_dm - h'/2
    ! In synchronous gauge, h'/2 = k*z (where z = (h'+6*eta')/(2k))
    ! But actually h' contribution is through -k*z for density
    ayprime(dm_ix) = -k * v_dm - k * z

    ! v_dm' = -aH * v_dm + cs2_dm * k * delta_dm - Gamma_DR * (v_dm - F_1/4)
    ! F_1/4 is the DR velocity (v_DR)
    ayprime(vc_ix) = -adotoa * v_dm + this%cs2_dm * k * clx_dm &
        - Gamma_DR * (v_dm - y(dr_ix + 1) / 4._dl)

    ! ============ DR Boltzmann hierarchy ============
    use_ufa = .false.
    if (present(high_ktau_dr)) use_ufa = high_ktau_dr

    if (use_ufa) then
        ! UFA: only evolve F_0, F_1, F_2 (mirrors neutrino UFA, arXiv:1104.2933)
        ayprime(dr_ix) = -(4._dl/3._dl)*k*z - k*y(dr_ix+1)
        ayprime(dr_ix+1) = k/3._dl*(y(dr_ix) - 2._dl*y(dr_ix+2)) &
            + Gamma_DR_tilde*(v_dm - y(dr_ix+1)/4._dl)
        ! UFA closure for pi: pi' = -3*cothxor*pi - F_0' + shear source
        ayprime(dr_ix+2) = -3._dl*cothxor*y(dr_ix+2) - ayprime(dr_ix) &
            + 8._dl/15._dl*k*sigma &
            - this%alpha_l_dark * Gamma_DR_tilde * y(dr_ix+2)
    else
        ! Full hierarchy
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
            if (l == 2) then
                ayprime(ix) = ayprime(ix) + 8._dl/15._dl * k * sigma
            end if
        end do

        if (this%lmax_dr >= 2) then
            ix = dr_ix + this%lmax_dr
            ayprime(ix) = k * y(ix-1) - real(this%lmax_dr + 1, dl) * cothxor * y(ix) &
                - this%alpha_l_dark * Gamma_DR_tilde * y(ix)
        end if
    end if

    end subroutine TDMDR_ETHOS_PerturbationEvolve

    subroutine TDMDR_ETHOS_SetBackgroundDensities(this, grhocrit, grhor, h2, grhodmdr, grhodr)
    class(TDMDR_ETHOS), intent(in) :: this
    real(dl), intent(in) :: grhocrit, grhor, h2
    real(dl), intent(out) :: grhodmdr, grhodr

    grhodmdr = grhocrit * this%omdmdrh2 / h2
    grhodr = grhor * this%N_dark

    end subroutine TDMDR_ETHOS_SetBackgroundDensities

    subroutine TDMDR_ETHOS_PrintFeedback(this, FeedbackLevel)
    class(TDMDR_ETHOS) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("DM-DR (ETHOS): omdmdrh2 = ",ES12.4)') this%omdmdrh2
        write(*,'("               N_dark = ",F8.4)') this%N_dark
        write(*,'("               lmax_dr = ",I4)') this%lmax_dr
    end if

    end subroutine TDMDR_ETHOS_PrintFeedback

    end module DMDR_ETHOS
