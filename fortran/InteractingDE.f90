    module InteractingDE
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none

    ! Interacting Dark Energy (IDE) model
    ! DM-DE energy-momentum exchange: Q_mu
    ! Background: rho_c' + 3H*rho_c = Q, rho_de' + 3(1+w)H*rho_de = -Q
    ! Interaction kernels:
    !   type 1: Q = xi * H * rho_de
    !   type 2: Q = xi * H * rho_c
    !   type 3: Q = xi * H * (rho_c + rho_de)
    !
    ! Perturbations follow the covariant approach of Valiviita+ 2008
    ! In the CDM rest frame, Q^mu = (Q, 0) gives:
    !   delta_c' = -k*v_c - h'/2 + Q/(rho_c) * (delta_de - delta_c)
    !   v_c' = -H*v_c  (CDM remains pressureless in rest frame)
    !   delta_de' = -(1+w)(k*v_de + h'/2) - 3H(cs2-w)*delta_de + Q/rho_de*(delta_c-delta_de) ...
    !   v_de' = ... (additional momentum transfer terms)
    !
    ! References: Valiviita+ 2008, He+ 2009, Gavela+ 2009, Costa+ 2017
    !             IDECAMB (arXiv:2306.01593)

    integer, parameter :: ide_Q_H_rho_de = 1   ! Q = xi*H*rho_de
    integer, parameter :: ide_Q_H_rho_c = 2    ! Q = xi*H*rho_c
    integer, parameter :: ide_Q_H_rho_tot = 3  ! Q = xi*H*(rho_c+rho_de)

    type, extends(TDarkEnergyEqnOfState) :: TInteractingDE
        real(dl) :: xi_ide = 0._dl        ! coupling strength (dimensionless)
        integer  :: interaction_type = 1   ! which Q kernel (1,2,3)
        real(dl) :: cs2_ide = 1._dl       ! DE rest-frame sound speed squared
        ! Cached quantities
        real(dl) :: grhoc_init = 0._dl    ! initial CDM density for modified evolution
    contains
    procedure :: ReadParams => TInteractingDE_ReadParams
    procedure, nopass :: PythonClass => TInteractingDE_PythonClass
    procedure, nopass :: SelfPointer => TInteractingDE_SelfPointer
    procedure :: Init => TInteractingDE_Init
    procedure :: BackgroundDensityAndPressure => TInteractingDE_BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy => TInteractingDE_PerturbedStressEnergy
    procedure :: PerturbationEvolve => TInteractingDE_PerturbationEvolve
    procedure :: PerturbationInitial => TInteractingDE_PerturbationInitial
    procedure :: PrintFeedback => TInteractingDE_PrintFeedback
    procedure :: SetIDEParams => TInteractingDE_SetIDEParams
    end type TInteractingDE

    contains

    subroutine TInteractingDE_SetIDEParams(this, xi, itype, cs2)
    class(TInteractingDE), intent(inout) :: this
    real(dl), intent(in) :: xi
    integer, intent(in) :: itype
    real(dl), intent(in) :: cs2

    this%xi_ide = xi
    this%interaction_type = itype
    this%cs2_ide = cs2

    end subroutine TInteractingDE_SetIDEParams

    subroutine TInteractingDE_ReadParams(this, Ini)
    use IniObjects
    class(TInteractingDE) :: this
    class(TIniFile), intent(in) :: Ini

    call this%TDarkEnergyEqnOfState%ReadParams(Ini)
    this%xi_ide = Ini%Read_Double('xi_ide', 0._dl)
    this%interaction_type = Ini%Read_Int('interaction_type', 1)
    this%cs2_ide = Ini%Read_Double('cs2_ide', 1._dl)

    end subroutine TInteractingDE_ReadParams

    function TInteractingDE_PythonClass()
    character(LEN=:), allocatable :: TInteractingDE_PythonClass
    TInteractingDE_PythonClass = 'InteractingDE'
    end function TInteractingDE_PythonClass

    subroutine TInteractingDE_SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type(TInteractingDE), pointer :: PType
    class(TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType

    end subroutine TInteractingDE_SelfPointer

    subroutine TInteractingDE_Init(this, State)
    use classes
    class(TInteractingDE), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State
    real(dl) :: h2

    call this%TDarkEnergyEqnOfState%Init(State)

    this%is_cosmological_constant = .false.  ! always evolve perturbations

    if (this%xi_ide /= 0._dl) then
        ! Need 2 perturbation equations: delta_de and (1+w)*v_de
        this%num_perturb_equations = 2

        select type(S => State)
        class is (CAMBdata)
            h2 = (S%CP%H0/100._dl)**2
            this%grhoc_init = S%grhocrit * S%CP%omch2 / h2
        end select
    else
        if (this%TDarkEnergyEqnOfState%is_cosmological_constant) then
            this%num_perturb_equations = 0
        else
            this%num_perturb_equations = 2
        end if
    end if

    end subroutine TInteractingDE_Init

    subroutine TInteractingDE_BackgroundDensityAndPressure(this, grhov, a, grhov_t, w)
    ! Modified background: the interaction changes the standard DE evolution
    ! For Q = xi*H*rho_de (type 1):
    !   rho_de ~ a^{-3(1+w)-xi}  instead of a^{-3(1+w)}
    !   rho_c gets compensating modification
    ! For Q = xi*H*rho_c (type 2):
    !   rho_c ~ a^{-3+xi} and rho_de modified accordingly
    class(TInteractingDE), intent(inout) :: this
    real(dl), intent(in) :: grhov, a
    real(dl), intent(out) :: grhov_t
    real(dl), optional, intent(out) :: w
    real(dl) :: w_val, grho_de_standard

    if (a <= 1e-10_dl) then
        grhov_t = 0._dl
        if (present(w)) w = this%w_de(a)
        return
    end if

    w_val = this%w_de(a)
    if (present(w)) w = w_val

    if (this%xi_ide == 0._dl) then
        ! No interaction, use standard evolution
        call this%TDarkEnergyEqnOfState%BackgroundDensityAndPressure(grhov, a, grhov_t, w)
        return
    end if

    ! Modified scaling due to interaction
    ! For type 1 (Q = xi*H*rho_de): rho_de propto a^{-3(1+w)-xi}
    ! grhov_t = grhov * a^2 * (a/a0)^{-3(1+w)-xi} = grhov * a^{-1-3w-xi}
    ! Recall grhov = 8*pi*G*rho_de,0 and grhov_t = 8*pi*G*rho_de*a^2

    select case(this%interaction_type)
    case(ide_Q_H_rho_de)
        ! rho_de(a)/rho_de(1) = a^{-3(1+w)-xi}
        ! grhov_t = grhov * a^{4} * a^{-3(1+w)-xi} / a^2 = grhov * a^{1-3w-xi} (wrong)
        ! Let's be careful: grhov_t = 8*pi*G*rho_de(a)*a^2
        ! rho_de(a) = rho_de(1) * a^{-3(1+w_lam)-xi}   [for constant w]
        ! grhov_t = grhov * a^{2-3(1+w_lam)-xi}
        ! For w=-1: grhov_t = grhov * a^{2-xi}
        if (.not. this%use_tabulated_w) then
            grhov_t = grhov * a**(1._dl - 3._dl*this%w_lam - 3._dl*this%wa - this%xi_ide)
            if (this%wa /= 0._dl) then
                grhov_t = grhov_t * exp(-3._dl*this%wa*(1._dl - a))
            end if
        else
            ! For tabulated w, modify the standard result by the xi correction
            call this%TDarkEnergyEqnOfState%BackgroundDensityAndPressure(grhov, a, grho_de_standard)
            grhov_t = grho_de_standard * a**(-this%xi_ide)
        end if

    case(ide_Q_H_rho_c)
        ! Q = xi*H*rho_c modifies CDM scaling, DE gets indirect modification
        ! For leading order, DE evolution is approximately standard
        ! The CDM modification is handled in equations.f90 via the source term
        call this%TDarkEnergyEqnOfState%BackgroundDensityAndPressure(grhov, a, grhov_t)

    case(ide_Q_H_rho_tot)
        ! Mixed case: approximate as superposition
        if (.not. this%use_tabulated_w) then
            grhov_t = grhov * a**(1._dl - 3._dl*this%w_lam - 3._dl*this%wa - this%xi_ide*0.5_dl)
            if (this%wa /= 0._dl) then
                grhov_t = grhov_t * exp(-3._dl*this%wa*(1._dl - a))
            end if
        else
            call this%TDarkEnergyEqnOfState%BackgroundDensityAndPressure(grhov, a, grho_de_standard)
            grhov_t = grho_de_standard * a**(-this%xi_ide*0.5_dl)
        end if

    case default
        call this%TDarkEnergyEqnOfState%BackgroundDensityAndPressure(grhov, a, grhov_t)
    end select

    end subroutine TInteractingDE_BackgroundDensityAndPressure

    subroutine TInteractingDE_PerturbationInitial(this, y, a, tau, k)
    class(TInteractingDE), intent(in) :: this
    real(dl), intent(out) :: y(:)
    real(dl), intent(in) :: a, tau, k

    ! Adiabatic initial conditions
    y = 0

    end subroutine TInteractingDE_PerturbationInitial

    subroutine TInteractingDE_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1, ay, ayprime, w_ix)
    class(TInteractingDE), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhov_t, w, gpres_noDE, etak, adotoa, k, kf1
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: w_ix

    if (this%num_perturb_equations == 0) then
        dgrhoe = 0
        dgqe = 0
        return
    end if

    ! delta_de stored in ay(w_ix), (1+w)*v_de in ay(w_ix+1)
    dgrhoe = ay(w_ix) * grhov_t
    dgqe = ay(w_ix + 1) * grhov_t

    end subroutine TInteractingDE_PerturbedStressEnergy

    subroutine TInteractingDE_PerturbationEvolve(this, ayprime, w, w_ix, &
        a, adotoa, k, z, y)
    class(TInteractingDE), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, w, k, z, y(:)
    integer, intent(in) :: w_ix
    real(dl) :: Hv3_over_k, Q_over_rho_de, delta_de, v_de_1pw
    real(dl) :: cs2, loga

    if (this%num_perturb_equations == 0) return

    delta_de = y(w_ix)
    v_de_1pw = y(w_ix + 1)  ! (1+w)*v_de
    cs2 = this%cs2_ide

    Hv3_over_k = 3._dl * adotoa * v_de_1pw / k / max(abs(1._dl + w), 1e-6_dl)

    ! Compute Q/(rho_de * H) = xi for type 1
    ! The perturbation equations follow Valiviita+ 2008, Gavela+ 2009
    select case(this%interaction_type)
    case(ide_Q_H_rho_de)
        Q_over_rho_de = this%xi_ide * adotoa

    case(ide_Q_H_rho_c)
        ! Q = xi*H*rho_c => Q/rho_de = xi*H*(rho_c/rho_de)
        ! We don't have direct access to rho_c here, approximate using grhoc
        Q_over_rho_de = 0._dl  ! simplified: interaction enters CDM equation primarily

    case(ide_Q_H_rho_tot)
        Q_over_rho_de = this%xi_ide * adotoa * 0.5_dl

    case default
        Q_over_rho_de = 0._dl
    end select

    ! DE density perturbation equation:
    ! delta_de' = -3H(cs2-w)*delta_de - (1+w)*k*v_de - (1+w)*k*z
    !             - 3H(cs2-w)*(1+w)*3H*v_de/k  [the Hv3_over_k term]
    !             + Q_interaction_perturbation_terms
    ! For type 1 (Q=xi*H*rho_de):
    !   Additional: +xi*H*(delta_c - delta_de) where delta_c comes from constraint
    !   At leading order, the xi terms modify the effective w evolution
    ayprime(w_ix) = -3._dl * adotoa * (cs2 - w) * (delta_de + (1._dl + w) * Hv3_over_k) &
        - (1._dl + w) * k * v_de_1pw / max(abs(1._dl + w), 1e-6_dl) &
        - (1._dl + w) * k * z

    ! Add w derivative term for wa parametrization
    if (.not. this%use_tabulated_w) then
        if (this%wa /= 0._dl) then
            ayprime(w_ix) = ayprime(w_ix) + Hv3_over_k * this%wa * adotoa * a
        end if
    else
        loga = log(a)
        if (loga > this%equation_of_state%Xmin_interp .and. &
            loga < this%equation_of_state%Xmax_interp) then
            ayprime(w_ix) = ayprime(w_ix) - adotoa * this%equation_of_state%Derivative(loga) * Hv3_over_k
        end if
    end if

    ! Interaction source terms for delta_de:
    ! Type 1: -xi*H*delta_de (self-damping; the +xi*H*delta_c is added in equations.f90)
    ! Type 2 and 3: handled entirely in equations.f90 where rho_c/rho_de is available
    if (this%xi_ide /= 0._dl .and. this%interaction_type == ide_Q_H_rho_de) then
        ayprime(w_ix) = ayprime(w_ix) - this%xi_ide * adotoa * delta_de
    end if

    ! DE velocity equation:
    ! (1+w)*v_de' = -H*(1-3cs2)*(1+w)*v_de + k*cs2*delta_de
    !             + interaction terms
    if (abs(1._dl + w) > 1e-6_dl) then
        ayprime(w_ix + 1) = -adotoa * (1._dl - 3._dl * cs2) * v_de_1pw + &
            k * cs2 * delta_de
        ! Interaction term in velocity
        if (this%xi_ide /= 0._dl .and. this%interaction_type == ide_Q_H_rho_de) then
            ayprime(w_ix + 1) = ayprime(w_ix + 1) - this%xi_ide * adotoa * v_de_1pw
        end if
    else
        ayprime(w_ix + 1) = 0._dl
    end if

    end subroutine TInteractingDE_PerturbationEvolve

    subroutine TInteractingDE_PrintFeedback(this, FeedbackLevel)
    class(TInteractingDE) :: this
    integer, intent(in) :: FeedbackLevel

    if (FeedbackLevel > 0) then
        write(*,'("Interacting DE: xi = ",ES12.4,", type = ",I2)') &
            this%xi_ide, this%interaction_type
        write(*,'("                (w0, wa) = (",F8.5,", ",F8.5,")")') &
            this%w_lam, this%wa
        write(*,'("                cs2_de = ",F8.5)') this%cs2_ide
    end if

    end subroutine TInteractingDE_PrintFeedback

    end module InteractingDE
