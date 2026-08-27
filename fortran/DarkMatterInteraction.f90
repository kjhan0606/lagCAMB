    module DarkMatterInteraction
    use precision
    use classes
    implicit none

    private

    type, extends(TCambComponent) :: TDarkMatterModel
        logical :: is_standard_cdm = .true.
        integer :: num_perturb_equations = 0
        logical :: has_cdm_velocity = .false.  ! whether CDM velocity is evolved (DM-baryon)
        integer :: num_dr_equations = 0        ! number of DR hierarchy equations (DM-DR)
    contains
    procedure :: Init => TDarkMatterModel_Init
    procedure :: BackgroundDensityAndPressure
    procedure :: PerturbedStressEnergy => TDarkMatterModel_PerturbedStressEnergy
    procedure :: PerturbationInitial => TDarkMatterModel_PerturbationInitial
    procedure :: PerturbationEvolve => TDarkMatterModel_PerturbationEvolve
    procedure :: PrintFeedback => TDarkMatterModel_PrintFeedback
    procedure :: DragRate_DM => TDarkMatterModel_DragRate_DM
    procedure :: DragRate_Baryon => TDarkMatterModel_DragRate_Baryon
    procedure :: DragRate_Photon => TDarkMatterModel_DragRate_Photon
    procedure :: SetBackgroundDensities => TDarkMatterModel_SetBackgroundDensities
    procedure :: ModifyNeutrinoParams => TDarkMatterModel_ModifyNeutrinoParams
    procedure :: TransferFunction => TDarkMatterModel_TransferFunction
    procedure :: NewtonianVelocityTransfer => TDarkMatterModel_NewtonianVelocityTransfer
    procedure :: ComponentNewtonianVelocityTransfer => TDarkMatterModel_ComponentNewtonianVelocityTransfer
    procedure :: NonNuDensityTransfer => TDarkMatterModel_NonNuDensityTransfer
    end type TDarkMatterModel

    public TDarkMatterModel

    contains

    subroutine TDarkMatterModel_Init(this, State)
    class(TDarkMatterModel), intent(inout) :: this
    class(TCAMBdata), intent(in), target :: State

    end subroutine TDarkMatterModel_Init

    subroutine BackgroundDensityAndPressure(this, grhodm, a, grhodm_t, gpres_dm)
    class(TDarkMatterModel), intent(inout) :: this
    real(dl), intent(in) :: grhodm, a
    real(dl), intent(out) :: grhodm_t
    real(dl), optional, intent(out) :: gpres_dm

    ! Standard CDM: pressureless dust, rho ~ a^-3
    if (a > 0) then
        grhodm_t = grhodm / a
    else
        grhodm_t = 0._dl
    end if
    if (present(gpres_dm)) gpres_dm = 0._dl

    end subroutine BackgroundDensityAndPressure

    subroutine TDarkMatterModel_PerturbedStressEnergy(this, dgrhoe, dgqe, &
        a, dgq, dgrho, grho, grhoc_t, adotoa, k, ay, ayprime, dm_ix, dr_ix, vc_ix)
    class(TDarkMatterModel), intent(inout) :: this
    real(dl), intent(out) :: dgrhoe, dgqe
    real(dl), intent(in) :: a, dgq, dgrho, grho, grhoc_t, adotoa, k
    real(dl), intent(in) :: ay(*)
    real(dl), intent(inout) :: ayprime(*)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix

    dgrhoe = 0
    dgqe = 0

    end subroutine TDarkMatterModel_PerturbedStressEnergy

    subroutine TDarkMatterModel_PerturbationInitial(this, y, a, tau, k, dm_ix, dr_ix, vc_ix)
    class(TDarkMatterModel), intent(in) :: this
    real(dl), intent(inout) :: y(:)
    real(dl), intent(in) :: a, tau, k
    integer, intent(in) :: dm_ix, dr_ix, vc_ix

    end subroutine TDarkMatterModel_PerturbationInitial

    subroutine TDarkMatterModel_PerturbationEvolve(this, ayprime, a, adotoa, k, z, y, &
        dm_ix, dr_ix, vc_ix, clxc, vb, grhoc_t, grhob_t, sigma, high_ktau_dr)
    class(TDarkMatterModel), intent(in) :: this
    real(dl), intent(inout) :: ayprime(:)
    real(dl), intent(in) :: a, adotoa, k, z, y(:)
    integer, intent(in) :: dm_ix, dr_ix, vc_ix
    real(dl), intent(in) :: clxc, vb, grhoc_t, grhob_t, sigma
    logical, intent(in), optional :: high_ktau_dr

    end subroutine TDarkMatterModel_PerturbationEvolve

    subroutine TDarkMatterModel_PrintFeedback(this, FeedbackLevel)
    class(TDarkMatterModel) :: this
    integer, intent(in) :: FeedbackLevel

    end subroutine TDarkMatterModel_PrintFeedback

    function TDarkMatterModel_DragRate_DM(this, a, adotoa, grhoc_t, grhob_t) result(R_drag)
    class(TDarkMatterModel), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t, grhob_t
    real(dl) :: R_drag
    R_drag = 0._dl
    end function TDarkMatterModel_DragRate_DM

    function TDarkMatterModel_DragRate_Baryon(this, a, adotoa, grhoc_t, grhob_t) result(R_drag)
    class(TDarkMatterModel), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t, grhob_t
    real(dl) :: R_drag
    R_drag = 0._dl
    end function TDarkMatterModel_DragRate_Baryon

    function TDarkMatterModel_DragRate_Photon(this, a, adotoa, grhoc_t, grhog_t) result(mu_dot)
    class(TDarkMatterModel), intent(in) :: this
    real(dl), intent(in) :: a, adotoa, grhoc_t, grhog_t
    real(dl) :: mu_dot
    mu_dot = 0._dl
    end function TDarkMatterModel_DragRate_Photon

    subroutine TDarkMatterModel_SetBackgroundDensities(this, grhocrit, grhor, h2, grhodmdr, grhodr)
    class(TDarkMatterModel), intent(in) :: this
    real(dl), intent(in) :: grhocrit, grhor, h2
    real(dl), intent(out) :: grhodmdr, grhodr
    grhodmdr = 0._dl
    grhodr = 0._dl
    end subroutine TDarkMatterModel_SetBackgroundDensities

    subroutine TDarkMatterModel_ModifyNeutrinoParams(this, max_nu_val, &
        Nu_mass_eigenstates, Nu_mass_degeneracies, Nu_mass_fractions, &
        Nu_mass_numbers, Num_Nu_Massive, Num_Nu_Massless, omnuh2, omch2, H0)
    class(TDarkMatterModel), intent(inout) :: this
    integer, intent(in) :: max_nu_val
    integer, intent(inout) :: Nu_mass_eigenstates, Num_Nu_Massive
    real(dl), intent(inout) :: Nu_mass_degeneracies(:), Nu_mass_fractions(:)
    integer, intent(inout) :: Nu_mass_numbers(:)
    real(dl), intent(inout) :: Num_Nu_Massless, omnuh2, omch2, H0
    ! Default: do nothing
    end subroutine TDarkMatterModel_ModifyNeutrinoParams

    function TDarkMatterModel_TransferFunction(this, k_h) result(Tk)
    !Optional small-scale matter power suppression T(k); P(k) is multiplied by T(k)^2.
    !Default returns 1 (no suppression). Override in models that act as a
    !phenomenological transfer-function modifier on top of CDM.
    class(TDarkMatterModel), intent(in) :: this
    real(dl), intent(in) :: k_h
    real(dl) :: Tk
    Tk = 1._dl
    end function TDarkMatterModel_TransferFunction

    function TDarkMatterModel_NewtonianVelocityTransfer(this, a, k, adotoa, sigma, ay, vc_ix) result(Tv)
    ! Newtonian-gauge velocity transfer written to the standard CAMB v_cdm
    ! column. Models with a pressureless subcomponent velocity may override
    ! this with the mass-weighted total-DM velocity needed by DMO ICs.
    class(TDarkMatterModel), intent(in) :: this
    real(dl), intent(in) :: a, k, adotoa, sigma
    real(dl), intent(in) :: ay(*)
    integer, intent(in) :: vc_ix
    real(dl) :: Tv

    Tv = -k * sigma / adotoa
    end function TDarkMatterModel_NewtonianVelocityTransfer

    function TDarkMatterModel_ComponentNewtonianVelocityTransfer(this, a, k, adotoa, sigma, ay, vc_ix) result(Tv)
    ! Optional Newtonian velocity of a separately evolved DM subcomponent.
    class(TDarkMatterModel), intent(in) :: this
    real(dl), intent(in) :: a, k, adotoa, sigma
    real(dl), intent(in) :: ay(*)
    integer, intent(in) :: vc_ix
    real(dl) :: Tv

    Tv = 0._dl
    end function TDarkMatterModel_ComponentNewtonianVelocityTransfer

    function TDarkMatterModel_NonNuDensityTransfer(this, a, grhob_t, grhodm_t, clxb, clxc, ay, dm_ix) result(Tdelta)
    ! Density transfer of the particle species (baryons plus all DM), excluding
    ! massive neutrinos. Models with a separately evolved DM fraction override
    ! this to include that component in the mass-weighted density.
    class(TDarkMatterModel), intent(in) :: this
    real(dl), intent(in) :: a, grhob_t, grhodm_t, clxb, clxc
    real(dl), intent(in) :: ay(*)
    integer, intent(in) :: dm_ix
    real(dl) :: Tdelta

    Tdelta = (grhob_t * clxb + grhodm_t * clxc) / (grhob_t + grhodm_t)
    end function TDarkMatterModel_NonNuDensityTransfer

    end module DarkMatterInteraction
