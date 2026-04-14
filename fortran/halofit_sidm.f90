    module NonLinear_SIDM
    use NonLinear
    use results
    use classes
    use constants
    implicit none
    private

    ! Modified HMcode for Self-Interacting Dark Matter (SIDM)
    ! Modifies halo concentration and density profile to account for
    ! core formation due to DM self-interactions.
    !
    ! Key effects:
    ! 1. Reduced halo concentration: c_SIDM = c_CDM * f_conc
    ! 2. Core formation: NFW cusp -> isothermal core transition
    ! Both effects suppress small-scale power (k > 1 h/Mpc)
    !
    ! Parameters:
    !   sigma_over_m: self-interaction cross section / DM mass [cm^2/g]
    !   v_dep_power:  velocity dependence power (sigma ~ v^alpha)
    !   v_dep_v0:     reference velocity [km/s]

    type, extends(THalofit) :: TSIDMHalofit
        real(dl) :: sigma_over_m = 0._dl     ! sigma/m [cm^2/g]
        real(dl) :: v_dep_power = 0._dl      ! sigma(v) = sigma_0 * (v/v0)^alpha
        real(dl) :: v_dep_v0 = 100._dl       ! reference velocity [km/s]
    contains
    procedure :: ReadParams => TSIDMHalofit_ReadParams
    procedure :: GetNonLinRatios => TSIDMHalofit_GetNonLinRatios
    procedure, nopass :: PythonClass => TSIDMHalofit_PythonClass
    procedure, nopass :: SelfPointer => TSIDMHalofit_SelfPointer
    end type TSIDMHalofit

    public TSIDMHalofit

    contains

    subroutine TSIDMHalofit_ReadParams(this, Ini)
    use IniObjects
    class(TSIDMHalofit) :: this
    class(TIniFile), intent(in) :: Ini

    call this%THalofit%ReadParams(Ini)
    this%sigma_over_m = Ini%Read_Double('sigma_over_m', 0._dl)
    this%v_dep_power = Ini%Read_Double('v_dep_power', 0._dl)
    this%v_dep_v0 = Ini%Read_Double('v_dep_v0', 100._dl)

    end subroutine TSIDMHalofit_ReadParams

    function TSIDMHalofit_PythonClass()
    character(LEN=:), allocatable :: TSIDMHalofit_PythonClass
    TSIDMHalofit_PythonClass = 'SIDMHalofit'
    end function TSIDMHalofit_PythonClass

    subroutine TSIDMHalofit_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TSIDMHalofit), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TSIDMHalofit_SelfPointer

    subroutine TSIDMHalofit_GetNonLinRatios(this, State, CAMB_Pk)
    ! Apply SIDM corrections on top of standard halofit/HMcode
    class(TSIDMHalofit) :: this
    class(TCAMBdata) :: State
    type(MatterPowerData), target :: CAMB_Pk
    real(dl) :: sigma_eff, v_rms, conc_factor, r_core_factor
    real(dl) :: kh, suppression, k_core
    integer :: ik, iz
    real(dl), parameter :: sigma_ref = 1._dl  ! reference sigma/m = 1 cm^2/g
    real(dl), parameter :: beta_core = 0.5_dl ! core radius scaling exponent

    ! First compute standard nonlinear ratios
    call this%THalofit%GetNonLinRatios(State, CAMB_Pk)

    ! If no SIDM, return standard result
    if (this%sigma_over_m == 0._dl) return

    ! Apply SIDM suppression to nonlinear power spectrum
    ! The suppression is a function of k and z
    do iz = 1, CAMB_Pk%num_z
        do ik = 1, CAMB_Pk%num_k
            kh = exp(CAMB_Pk%log_kh(ik))

            ! Effective cross section at typical halo velocity
            v_rms = 200._dl  ! typical cluster velocity in km/s
            if (this%v_dep_power /= 0._dl) then
                sigma_eff = this%sigma_over_m * (v_rms / this%v_dep_v0)**this%v_dep_power
            else
                sigma_eff = this%sigma_over_m
            end if

            ! Concentration reduction factor
            ! f_conc = 1 / (1 + sigma_eff/sigma_ref)^0.3
            ! This approximately captures the effect of core formation
            ! on the halo concentration-mass relation
            conc_factor = 1._dl / (1._dl + (sigma_eff/sigma_ref)**0.3_dl)

            ! Core formation scale in k-space
            ! k_core ~ 1/r_core, where r_core grows with sigma/m
            r_core_factor = (sigma_eff / sigma_ref)**beta_core
            k_core = 10._dl / r_core_factor  ! h/Mpc, typical core scale

            ! Suppression of 1-halo term at k > k_core
            ! Model: S(k) = 1 / (1 + (k/k_core)^2)^(f_c/2) where f_c depends on conc_factor
            if (kh > 0.5_dl) then
                suppression = 1._dl / (1._dl + ((kh/k_core)**2) * (1._dl - conc_factor))**0.5_dl
            else
                suppression = 1._dl
            end if

            ! Apply suppression to the nonlinear ratio
            ! nonlin_ratio = sqrt(P_NL/P_L), so we modify it directly
            CAMB_Pk%nonlin_ratio(ik, iz) = CAMB_Pk%nonlin_ratio(ik, iz) * sqrt(suppression)
        end do
    end do

    end subroutine TSIDMHalofit_GetNonLinRatios

    end module NonLinear_SIDM
