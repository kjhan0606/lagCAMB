    module MassiveNu
    use constants
    implicit none
    private

    real(dl), parameter  :: fermi_dirac_const  = 7._dl/120*const_pi**4 ! 5.68219698_dl
    !fermi_dirac_const = int q^3 F(q) dq = 7/120*pi^4
    real(dl), parameter  :: const2 = 5._dl/7._dl/const_pi**2   !0.072372274_dl

    !Steps for spline interpolation (use series outside this range)
    integer, parameter  :: nrhopn=400
    real(dl), parameter :: am_min = 0.3_dl
    !smallest a*m_nu to integrate distribution function rather than using series
    real(dl), parameter :: am_max = 70._dl
    !max a*m_nu to integrate

    !Actual range for using series (to avoid inaccurate ends of spline)
    real(dl), parameter :: am_minp=am_min + am_max/(nrhopn-1)*1.01_dl
    real(dl), parameter :: am_maxp=am_max*0.9_dl
    !Optimized 8-point background quadrature (derived via minimax/least-squares fit)
    ! set legacy toggle true to use (less accurate and slower) original 100-step grid.
    logical, parameter :: use_legacy_nu_background_grid = .false.
    real(dl), parameter :: nu_background_q(8) = (/0.2937822_dl, 0.73583979_dl, 1.49222507_dl, 2.68795368_dl, &
        4.30678084_dl, 4.63078102_dl, 7.37122449_dl, 11.91683009_dl/)
    real(dl), parameter :: nu_background_rho_weights(8) = (/0.000640376236953801_dl, 0.01312614_dl, 0.10233804_dl, &
        0.31935253_dl, 0.12193422_dl, 0.27616318_dl, 0.15455734_dl, 0.01189232_dl/)
    real(dl), parameter :: nu_background_pressure_weights(8) = (/0.0002028435952467642_dl, 0.00440170_dl, 0.03405480_dl, &
        0.10658874_dl, 0.03987692_dl, 0.09276772_dl, 0.05146992_dl, 0.00396969_dl/)

    Type TNuPerturbations
        !Sample for massive neutrino momentum
        !Default settings appear to be OK for P_k accuate at 1e-3 level
        integer nqmax !actual number of q modes evolves
        real(dl), allocatable ::  nu_q(:), nu_int_kernel(:)
    contains
    procedure :: init => TNuPerturbations_init
    end type TNuPerturbations

    Type TThermalNuBackground
        !Quantities for the neutrino background momentum distribution assuming thermal
        real(dl) dam !step in a*m
        real(dl), dimension(:), allocatable ::  r1,p1,dr1,dp1
        real(dl), dimension(:), allocatable ::  pp1, dpp1 !p_pseudo spline for CLASS UFA
        real(dl), dimension(:), allocatable ::  iv21, div21 !I_v2 spline for CLASS UFA
        real(dl), private :: target_rho
    contains
    procedure :: init => ThermalNuBackground_init
    procedure :: rho_P => ThermalNuBackground_rho_P
    procedure :: rho => ThermalNuBackground_rho
    procedure :: drho => ThermalNuBackground_drho
    procedure :: find_nu_mass_for_rho => ThermalNuBackground_find_nu_mass_for_rho
    procedure :: ppseudo => ThermalNuBackground_ppseudo
    procedure :: Iv2 => ThermalNuBackground_Iv2
    end type TThermalNuBackground

    Type TCustomNuPSD
        logical :: active = .false.
        integer :: nq_tab = 0
        real(dl), allocatable :: q_tab(:), f0_tab(:)
        real(dl), allocatable :: nu_int_kernel(:)
        real(dl) :: dam_bg = 0._dl
        real(dl), dimension(:), allocatable :: r1, p1, dr1, dp1
        real(dl), dimension(:), allocatable :: pp1, dpp1, iv21, div21
        real(dl), private :: target_rho = 0._dl
    contains
    procedure :: Init => TCustomNuPSD_Init
    procedure :: rho_P => TCustomNuPSD_rho_P
    procedure :: rho => TCustomNuPSD_rho
    procedure :: drho => TCustomNuPSD_drho
    procedure :: ppseudo => TCustomNuPSD_ppseudo
    procedure :: Iv2 => TCustomNuPSD_Iv2
    procedure :: find_nu_mass_for_rho => TCustomNuPSD_find_nu_mass
    end type TCustomNuPSD

    integer, parameter :: max_nu_custom = 10

    Type(TThermalNuBackground), target :: ThermalNuBackground
    class(TThermalNuBackground), pointer :: ThermalNuBack !ifort workaround
    Type(TCustomNuPSD), save, target :: CustomNuPSD(max_nu_custom)

    public fermi_dirac_const,  sum_mnu_for_m1, neutrino_mass_fac, TNuPerturbations, &
        ThermalNuBackground, ThermalNuBack, TCustomNuPSD, CustomNuPSD, &
        SetCustomNuPSD, ClearCustomNuPSD
    contains

    subroutine sum_mnu_for_m1(summnu,dsummnu, m1, targ, sgn)
    use constants
    real(dl), intent(in) :: m1, targ, sgn
    real(dl), intent(out) :: summnu, dsummnu
    real(dl) :: m2,m3

    m2 = sqrt(m1**2 + delta_mnu21)
    m3 = sqrt(m1**2 + sgn*delta_mnu31)
    summnu = m1 + m2 + m3 - targ
    dsummnu = m1/m2+m1/m3 + 1

    end subroutine sum_mnu_for_m1

    subroutine TNuPerturbations_init(this,Accuracy)
    !Set up which momenta to integrate the neutrino perturbations, depending on accuracy
    !Using three optimized momenta works very well in most cases
    class(TNuPerturbations) :: this
    real(dl), intent(in) :: Accuracy
    real(dl) :: dq,dlfdlq, q
    integer i

    this%nqmax=3
    if (Accuracy>1) this%nqmax=4
    if (Accuracy>2) this%nqmax=5
    if (Accuracy>3) this%nqmax=nint(Accuracy*10)
    !note this may well be worse than the 5 optimized points

    !We evolve evolve 4F_l/dlfdlq(i), so kernel includes dlfdlnq factor
    !Integration scheme gets (Fermi-Dirac thing)*q^n exact,for n=-4, -2..2
    !see CAMB notes and https://camb.info/maple/nu_integration_kernels.py
    if (allocated(this%nu_q)) deallocate(this%nu_q, this%nu_int_kernel)
    allocate(this%nu_q(this%nqmax))
    allocate(this%nu_int_kernel(this%nqmax))

    if (this%nqmax==3) then
        !Accurate at 2e-4 level
        this%nu_q(1:3) = (/0.913201, 3.37517, 7.79184/)
        this%nu_int_kernel(1:3) = (/0.0687359, 3.31435, 2.29911/)
    else if (this%nqmax==4) then
        !Free-node least-squares fit for n=-4,-2..2 and v(am/q), 1/v(am/q)
        !Original rule kept here for reference:
        !this%nu_q(1:4) = (/0.7, 2.62814, 5.90428, 12.0/)
        !this%nu_int_kernel(1:4) = (/0.0200251, 1.84539, 3.52736, 0.289427/)
        this%nu_q(1:4) = (/0.5802007037903776_dl, 2.2150938570691223_dl, 4.948032138986023_dl, 9.65253759848097_dl/)
        this%nu_int_kernel(1:4) = (/0.0082119845039711_dl, 1.1143258498419168_dl, &
            3.6819104154615907_dl, 0.8777790167504481_dl/)
    else if (this%nqmax==5) then
        !Exact for n=-4,-2..2 with remaining freedom fit to v(am/q), 1/v(am/q)
        !Original rule kept here for reference:
        !this%nu_q(1:5) = (/0.583165, 2.0, 4.0, 7.26582, 13.0/)
        !this%nu_int_kernel(1:5) = (/0.0081201, 0.689407, 2.8063, 2.05156, 0.126817/)
        this%nu_q(1:5) = (/0.4620995950854295_dl, 1.7331898360630928_dl, 3.7956972681313816_dl, &
            7.2113928588584990_dl, 13.2665914595911080_dl/)
        this%nu_int_kernel(1:5) = (/0.0026946402277849193_dl, 0.46041071394952310_dl, 2.9207114780286405_dl, &
            2.1821643017186352_dl, 0.11621584305889110_dl/)
    else
        dq = (12 + this%nqmax/5)/real(this%nqmax)
        do i=1,this%nqmax
            q=(i-0.5d0)*dq
            this%nu_q(i) = q
            dlfdlq=-q/(1._dl+exp(-q))
            this%nu_int_kernel(i)=dq*q**3/(exp(q)+1._dl) * (-0.25_dl*dlfdlq) !now evolve 4F_l/dlfdlq(i)
        end do
    end if
    this%nu_int_kernel=this%nu_int_kernel/fermi_dirac_const

    end subroutine TNuPerturbations_init

    subroutine ThermalNuBackground_init(this)
    use splines
    class(TThermalNuBackground) :: this
    !  Initialize interpolation tables for massive neutrino background.
    integer i
    real(dl) am, rhonu,pnu, drhonu_dam, dpnu_dam
    real(dl) ppnu, iv2nu, dppnu_dam, div2nu_dam
    real(dl) spline_data(nrhopn)

    if (allocated(this%r1)) return
    ThermalNuBack => ThermalNuBackground !ifort bug workaround

    allocate(this%r1(nrhopn),this%p1(nrhopn),this%dr1(nrhopn),this%dp1(nrhopn))
    allocate(this%pp1(nrhopn),this%dpp1(nrhopn),this%iv21(nrhopn),this%div21(nrhopn))
    this%dam=(am_max-am_min)/(nrhopn-1)

    if (use_legacy_nu_background_grid) then
        !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), &
        !$OMP& PRIVATE(am,rhonu,pnu,ppnu,iv2nu,dppnu_dam,div2nu_dam)
        do i=1,nrhopn
            am=am_min + (i-1)*this%dam
            call nuRhoPres(am,rhonu,pnu)
            this%r1(i)=rhonu
            this%p1(i)=pnu
            call nuPseudoPres_Iv2(am,ppnu,iv2nu,dppnu_dam,div2nu_dam)
            this%pp1(i)=ppnu
            this%iv21(i)=iv2nu
            this%dpp1(i)=dppnu_dam*this%dam
            this%div21(i)=div2nu_dam*this%dam
        end do
        !$OMP END PARALLEL DO

        call splini(spline_data,nrhopn)
        call splder(this%r1,this%dr1,nrhopn,spline_data)
        call splder(this%p1,this%dp1,nrhopn,spline_data)
    else
        !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), &
        !$OMP& PRIVATE(am,rhonu,pnu,drhonu_dam,dpnu_dam,ppnu,iv2nu,dppnu_dam,div2nu_dam)
        do i=1,nrhopn
            am=am_min + (i-1)*this%dam
            call nuRhoPres_8point(am,rhonu,pnu)
            call nuRhoPres_8point_derivs(am,drhonu_dam,dpnu_dam)
            this%r1(i)=rhonu
            this%p1(i)=pnu
            this%dr1(i)=drhonu_dam*this%dam
            this%dp1(i)=dpnu_dam*this%dam
            call nuPseudoPres_Iv2(am,ppnu,iv2nu,dppnu_dam,div2nu_dam)
            this%pp1(i)=ppnu
            this%iv21(i)=iv2nu
            this%dpp1(i)=dppnu_dam*this%dam
            this%div21(i)=div2nu_dam*this%dam
        end do
        !$OMP END PARALLEL DO
    end if

    end subroutine ThermalNuBackground_init

    subroutine nuRhoPres(am,rhonu,pnu)
    !  Compute the density and pressure of one eigenstate of massive neutrinos,
    !  in units of the mean density of one flavor of massless neutrinos.
    use splines
    real(dl),  parameter :: qmax=30._dl
    integer, parameter :: nq=100
    real(dl) dum1(nq+1),dum2(nq+1)
    real(dl), intent(in) :: am
    real(dl), intent(out) ::  rhonu,pnu
    integer i
    real(dl) q,aq,v,aqdn,adq

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    !  Integrate up to qmax and then use asymptotic expansion for remainder.
    adq=qmax/nq
    dum1(1)=0._dl
    dum2(1)=0._dl
    do  i=1,nq
        q=i*adq
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
        aqdn=adq*q*q*q/(exp(q)+1._dl)
        dum1(i+1)=aqdn/v
        dum2(i+1)=aqdn*v
    end do
    call splint(dum1,rhonu,nq+1)
    call splint(dum2,pnu,nq+1)
    !  Apply asymptotic corrrection for q>qmax and normalize by relativistic
    !  energy density.
    rhonu=(rhonu+dum1(nq+1)/adq)/fermi_dirac_const
    pnu=(pnu+dum2(nq+1)/adq)/fermi_dirac_const/3._dl

    end subroutine nuRhoPres

    subroutine nuRhoPres_8point(am,rhonu,pnu)
    !  Optimized 8-point shared-node quadrature for ThermalNuBackground table generation.
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu, pnu
    real(dl) inv_v, v
    integer i

    rhonu = 0._dl
    pnu = 0._dl
    do i=1,size(nu_background_q)
        inv_v = sqrt(1._dl + (am/nu_background_q(i))**2)
        v = 1._dl/inv_v
        rhonu = rhonu + nu_background_rho_weights(i)*inv_v
        pnu = pnu + nu_background_pressure_weights(i)*v
    end do

    end subroutine nuRhoPres_8point

    subroutine nuRhoPres_8point_derivs(am,drhonu_dam,dpnu_dam)
    !  Exact a*m derivatives of the optimized 8-point background quadrature.
    real(dl), intent(in) :: am
    real(dl), intent(out) :: drhonu_dam, dpnu_dam
    real(dl) :: inv_v, v, q2
    integer i

    drhonu_dam = 0._dl
    dpnu_dam = 0._dl
    do i=1,size(nu_background_q)
        q2 = nu_background_q(i)**2
        inv_v = sqrt(1._dl + (am/nu_background_q(i))**2)
        v = 1._dl/inv_v
        drhonu_dam = drhonu_dam + nu_background_rho_weights(i)*am*v/q2
        dpnu_dam = dpnu_dam - nu_background_pressure_weights(i)*am*v**3/q2
    end do

    end subroutine nuRhoPres_8point_derivs

    subroutine nuPseudoPres_Iv2(am, ppnu, iv2nu, dppnu_dam, div2nu_dam)
    !  Compute p_pseudo = (1/3)*int q^3*f0*v^3 dq / fermi_dirac_const
    !  and I_v2 = int q^3*f0*v^2 dq / fermi_dirac_const
    !  plus their derivatives w.r.t. am, using 100-point quadrature.
    use splines
    real(dl), intent(in) :: am
    real(dl), intent(out) :: ppnu, iv2nu, dppnu_dam, div2nu_dam
    real(dl), parameter :: qmax = 30._dl
    integer, parameter :: nq = 100
    real(dl) :: dum_pp(nq+1), dum_iv2(nq+1), dum_dpp(nq+1), dum_div2(nq+1)
    real(dl) :: adq, q, aq, v, v2, aqdn
    integer :: i

    adq = qmax/nq
    dum_pp(1) = 0._dl
    dum_iv2(1) = 0._dl
    dum_dpp(1) = 0._dl
    dum_div2(1) = 0._dl
    do i = 1, nq
        q = i*adq
        aq = am/q
        v2 = 1._dl/(1._dl + aq*aq)
        v = sqrt(v2)
        aqdn = adq*q*q*q/(exp(q) + 1._dl)
        dum_pp(i+1) = aqdn*v*v2          ! q^3 f0 v^3
        dum_iv2(i+1) = aqdn*v2            ! q^3 f0 v^2
        ! d(q^3 f0 v^3)/d(am) = q^3 f0 * 3v^2 * (-am/q^2) * v^3 = -3*am*q*f0*v^5
        dum_dpp(i+1) = -3._dl*am*adq*q/(exp(q) + 1._dl)*v2*v2*v
        ! d(q^3 f0 v^2)/d(am) = q^3 f0 * 2v * (-am/q^2) * v^3 = -2*am*q*f0*v^4
        dum_div2(i+1) = -2._dl*am*adq*q/(exp(q) + 1._dl)*v2*v2
    end do
    call splint(dum_pp, ppnu, nq+1)
    call splint(dum_iv2, iv2nu, nq+1)
    call splint(dum_dpp, dppnu_dam, nq+1)
    call splint(dum_div2, div2nu_dam, nq+1)
    ppnu = (ppnu + dum_pp(nq+1)/adq)/fermi_dirac_const/3._dl
    iv2nu = (iv2nu + dum_iv2(nq+1)/adq)/fermi_dirac_const
    dppnu_dam = (dppnu_dam + dum_dpp(nq+1)/adq)/fermi_dirac_const/3._dl
    div2nu_dam = (div2nu_dam + dum_div2(nq+1)/adq)/fermi_dirac_const

    end subroutine nuPseudoPres_Iv2

    subroutine ThermalNuBackground_rho_P(this,am,rhonu,pnu)
    class(TThermalNuBackground) :: this
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu, pnu
    real(dl) d, logam, am2
    integer i
    !  Compute massive neutrino density and pressure in units of the mean
    !  density of one eigenstate of massless neutrinos.  Use cubic splines to
    !  interpolate from a table. Accuracy generally better than 1e-5.

    if (am <= am_minp) then
        if (am< 0.01_dl) then
            rhonu=1._dl + const2*am**2
            pnu=(2-rhonu)/3._dl
        else
            !Higher order expansion result less obvious, Appendix A of arXiv:0911.2714
            am2=am**2
            logam = log(am)
            rhonu = 1+am2*(const2+am2*(.1099926669d-1*logam-.3492416767d-2-.5866275571d-2*am))
            pnu = (1+am2*(-const2+am2*(-.3299780009d-1*logam-.5219952794d-3+.2346510229d-1*am)))/3
        end if
        return
    else if (am >= am_maxp) then
        !Simple series solution (expanded in 1/(a*m))
        rhonu = 3/(2*fermi_dirac_const)*(zeta3*am + ((15*zeta5)/2 - 945._dl/16*zeta7/am**2)/am)
        pnu = 900._dl/120._dl/fermi_dirac_const*(zeta5-63._dl/4*Zeta7/am**2)/am
        return
    end if

    d=(am-am_min)/this%dam+1._dl
    i=int(d)
    d=d-i

    !  Cubic spline interpolation.
    rhonu=this%r1(i)+d*(this%dr1(i)+d*(3._dl*(this%r1(i+1)-this%r1(i))-2._dl*this%dr1(i) &
        -this%dr1(i+1)+d*(this%dr1(i)+this%dr1(i+1)+2._dl*(this%r1(i)-this%r1(i+1)))))
    pnu=this%p1(i)+d*(this%dp1(i)+d*(3._dl*(this%p1(i+1)-this%p1(i))-2._dl*this%dp1(i) &
        -this%dp1(i+1)+d*(this%dp1(i)+this%dp1(i+1)+2._dl*(this%p1(i)-this%p1(i+1)))))

    end subroutine ThermalNuBackground_rho_P

    subroutine ThermalNuBackground_rho(this,am,rhonu)
    class(TThermalNuBackground) :: this
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu
    real(dl) d, am2
    integer i

    !  Compute massive neutrino density in units of the mean
    !  density of one eigenstate of massless neutrinos.  Use series solutions or
    !  cubic splines to interpolate from a table.

    if (am <= am_minp) then
        if (am < 0.01_dl) then
            rhonu=1._dl + const2*am**2
        else
            am2=am**2
            rhonu = 1+am2*(const2+am2*(.1099926669d-1*log(am)-.3492416767d-2-.5866275571d-2*am))
        end if
        return
    else if (am >= am_maxp) then
        rhonu = 3/(2*fermi_dirac_const)*(zeta3*am + ((15*zeta5)/2 - 945._dl/16*zeta7/am**2)/am)
        return
    end if

    d=(am-am_min)/this%dam+1._dl
    i=int(d)
    d=d-i

    !  Cubic spline interpolation.
    rhonu=this%r1(i)+d*(this%dr1(i)+d*(3._dl*(this%r1(i+1)-this%r1(i))-2._dl*this%dr1(i) &
        -this%dr1(i+1)+d*(this%dr1(i)+this%dr1(i+1)+2._dl*(this%r1(i)-this%r1(i+1)))))

    end subroutine ThermalNuBackground_rho

    function rho_err(this, nu_mass)
    class(TThermalNuBackground) :: this
    real(dl) rho_err, nu_mass, rhonu

    call this%rho(nu_mass, rhonu)
    rho_err = rhonu - this%target_rho

    end function rho_err

    function ThermalNuBackground_find_nu_mass_for_rho(this,rho) result(nu_mass)
    !  Get eigenstate mass given input density (rho is neutrino density in units of one massless)
    !  nu_mass=m_n*c**2/(k_B*T_nu0).
    !  Get number density n of neutrinos from
    !  rho_massless/n = int q^3/(1+e^q) / int q^2/(1+e^q)=7/180 pi^4/Zeta(3)
    !  then m = Omega_nu/N_nu rho_crit /n if non-relativistic
    use MathUtils
    use config
    class(TThermalNuBackground) :: this
    real(dl), intent(in) :: rho
    real(dl) nu_mass, rhonu, rhonu1, delta
    real(dl) fzero
    integer iflag

    if (rho <= 1.001_dl) then
        !energy density all accounted for by massless result
        nu_mass=0
    else
        !Get mass assuming fully non-relativistic
        nu_mass=fermi_dirac_const/(1.5d0*zeta3)*rho

        if (nu_mass>4) then
            !  perturbative correction for velocity when nearly non-relativistic
            !  Error due to velocity < 1e-5 for mnu~0.06 but can easily correct (assuming non-relativistic today)
            !  Note that python does not propagate mnu to omnuh2 consistently to the same accuracy; but this makes
            !  fortran more internally consistent between input and computed Omega_nu h^2

            !Make perturbative correction for the tiny error due to the neutrino velocity
            call this%rho(nu_mass, rhonu)
            call this%rho(nu_mass*0.9, rhonu1)
            delta = rhonu - rho
            nu_mass = nu_mass*(1 + delta/((rhonu1 - rhonu)/0.1) )
        else
            !Directly solve to avoid issues with perturbative result when no longer very relativistic
            this%target_rho = rho
            call brentq(this,rho_err,0._dl,nu_mass,0.01_dl,nu_mass,fzero,iflag)
            if (iflag/=0) call GlobalError('find_nu_mass_for_rho failed to find neutrino mass')
        end if
    end if

    end function ThermalNuBackground_find_nu_mass_for_rho


    function ThermalNuBackground_ppseudo(this,am) result(ppnu)
    !  Compute pseudo-pressure p_pseudo = (1/3)*int q^3*f0*v^3 dq / fermi_dirac_const
    !  Used for CLASS-style UFA closure: c_g^2 = w*(5-p_pseudo/p)/(3*(1+w))
    class(TThermalNuBackground) :: this
    real(dl), intent(in) :: am
    real(dl) :: ppnu
    real(dl) :: d, am2, ppnu_dum, iv2_dum, dpp_dum, div2_dum
    integer :: i
    !  NR leading coefficient: 4050*zeta7/pi^4
    real(dl), parameter :: pp_nr_coeff = 4050._dl*zeta7/const_pi**4

    if (am <= am_minp) then
        if (am < 0.01_dl) then
            ppnu = (1._dl - 3._dl*const2*am**2)/3._dl
        else
            call nuPseudoPres_Iv2(am,ppnu,iv2_dum,dpp_dum,div2_dum)
        end if
        return
    else if (am >= am_maxp) then
        ppnu = pp_nr_coeff/am**3
        return
    end if

    d = (am - am_min)/this%dam + 1._dl
    i = int(d)
    d = d - i
    ppnu = this%pp1(i) + d*(this%dpp1(i) + d*(3._dl*(this%pp1(i+1) - this%pp1(i)) - 2._dl*this%dpp1(i) &
        - this%dpp1(i+1) + d*(this%dpp1(i) + this%dpp1(i+1) + 2._dl*(this%pp1(i) - this%pp1(i+1)))))

    end function ThermalNuBackground_ppseudo

    function ThermalNuBackground_Iv2(this,am) result(iv2nu)
    !  Compute I_v2 = int q^3*f0*v^2 dq / fermi_dirac_const
    !  Used for CLASS-style UFA: dynamic G11_t = q * I_v2
    class(TThermalNuBackground) :: this
    real(dl), intent(in) :: am
    real(dl) :: iv2nu
    real(dl) :: d, ppnu_dum, iv2_dum, dpp_dum, div2_dum
    integer :: i
    !  NR leading coefficient: 310*pi^2/147
    real(dl), parameter :: iv2_nr_coeff = 310._dl*const_pi**2/147._dl

    if (am <= am_minp) then
        if (am < 0.01_dl) then
            iv2nu = 1._dl - 2._dl*const2*am**2
        else
            call nuPseudoPres_Iv2(am,ppnu_dum,iv2nu,dpp_dum,div2_dum)
        end if
        return
    else if (am >= am_maxp) then
        iv2nu = iv2_nr_coeff/am**2
        return
    end if

    d = (am - am_min)/this%dam + 1._dl
    i = int(d)
    d = d - i
    iv2nu = this%iv21(i) + d*(this%div21(i) + d*(3._dl*(this%iv21(i+1) - this%iv21(i)) - 2._dl*this%div21(i) &
        - this%div21(i+1) + d*(this%div21(i) + this%div21(i+1) + 2._dl*(this%iv21(i) - this%iv21(i+1)))))

    end function ThermalNuBackground_Iv2

    function ThermalNuBackground_drho(this,am,adotoa) result (rhonudot)
    !  Compute the time derivative of the mean density in massive neutrinos
    class(TThermalNuBackground) :: this
    real(dl) adotoa,rhonudot
    real(dl) am2, rhonu, pnu
    real(dl), intent(IN) :: am

    if (am< am_minp) then
        !rhonudot = 2*const2*am**2*adotoa
        am2 = am**2
        rhonudot = am2 * (2 * const2 + am2 * (.4399706676d-1 * log(am) &
            - .2970400378d-2 - .29331377855d-1 * am)) * adotoa
    else if (am>am_maxp) then
        rhonudot = 3/(2*fermi_dirac_const)*(zeta3*am +( -(15*zeta5)/2 + 2835._dl/16*zeta7/am**2)/am)*adotoa
    else
        call this%rho_P(am,rhonu,pnu)
        ! am * (d rho_nu / d am) analytically simplifies exactly to (rho_nu - 3 P_nu)
        rhonudot = (rhonu - 3._dl*pnu)*adotoa
    end if

    end function ThermalNuBackground_drho

    !=====================================================================
    ! Custom neutrino PSD: module-level setter/clearer
    !=====================================================================

    subroutine SetCustomNuPSD(nu_i, nq, q_tab, f0_tab)
    integer, intent(in) :: nu_i, nq
    real(dl), intent(in) :: q_tab(nq), f0_tab(nq)

    if (nu_i < 1 .or. nu_i > max_nu_custom) then
        write(*,*) 'SetCustomNuPSD: invalid nu_i =', nu_i
        return
    end if
    CustomNuPSD(nu_i)%active = .true.
    CustomNuPSD(nu_i)%nq_tab = nq
    if (allocated(CustomNuPSD(nu_i)%q_tab)) deallocate(CustomNuPSD(nu_i)%q_tab)
    if (allocated(CustomNuPSD(nu_i)%f0_tab)) deallocate(CustomNuPSD(nu_i)%f0_tab)
    allocate(CustomNuPSD(nu_i)%q_tab(nq), CustomNuPSD(nu_i)%f0_tab(nq))
    CustomNuPSD(nu_i)%q_tab = q_tab
    CustomNuPSD(nu_i)%f0_tab = f0_tab

    end subroutine SetCustomNuPSD

    subroutine ClearCustomNuPSD()
    integer :: i
    do i = 1, max_nu_custom
        CustomNuPSD(i)%active = .false.
        CustomNuPSD(i)%nq_tab = 0
        if (allocated(CustomNuPSD(i)%q_tab)) deallocate(CustomNuPSD(i)%q_tab)
        if (allocated(CustomNuPSD(i)%f0_tab)) deallocate(CustomNuPSD(i)%f0_tab)
        if (allocated(CustomNuPSD(i)%nu_int_kernel)) deallocate(CustomNuPSD(i)%nu_int_kernel)
        if (allocated(CustomNuPSD(i)%r1)) deallocate(CustomNuPSD(i)%r1, CustomNuPSD(i)%p1, &
            CustomNuPSD(i)%dr1, CustomNuPSD(i)%dp1)
        if (allocated(CustomNuPSD(i)%pp1)) deallocate(CustomNuPSD(i)%pp1, CustomNuPSD(i)%dpp1, &
            CustomNuPSD(i)%iv21, CustomNuPSD(i)%div21)
    end do
    end subroutine ClearCustomNuPSD

    !=====================================================================
    ! Custom PSD helpers: interpolation
    !=====================================================================

    function custom_f0_at_q(psd, q) result(f0)
    type(TCustomNuPSD), intent(in) :: psd
    real(dl), intent(in) :: q
    real(dl) :: f0
    real(dl) :: t
    integer :: i

    if (q <= psd%q_tab(1)) then
        f0 = psd%f0_tab(1)
        return
    end if
    if (q >= psd%q_tab(psd%nq_tab)) then
        f0 = 0._dl
        return
    end if
    do i = 1, psd%nq_tab - 1
        if (q < psd%q_tab(i+1)) exit
    end do
    t = (q - psd%q_tab(i)) / (psd%q_tab(i+1) - psd%q_tab(i))
    f0 = psd%f0_tab(i) * (1._dl - t) + psd%f0_tab(i+1) * t

    end function custom_f0_at_q

    function custom_dlnf0dlnq(psd, q) result(dlnf)
    type(TCustomNuPSD), intent(in) :: psd
    real(dl), intent(in) :: q
    real(dl) :: dlnf
    real(dl) :: dq, f0p, f0m, f0c

    f0c = custom_f0_at_q(psd, q)
    if (f0c <= 0._dl) then
        dlnf = -q
        return
    end if
    dq = q * 0.01_dl
    if (dq < 1e-8_dl) dq = 1e-8_dl
    f0p = custom_f0_at_q(psd, q + dq)
    f0m = custom_f0_at_q(psd, q - dq)
    if (f0m <= 0._dl) f0m = f0c
    if (f0p <= 0._dl) f0p = f0c * exp(-2._dl * dq / q)
    dlnf = log(f0p / f0m) / (2._dl * dq) * q

    end function custom_dlnf0dlnq

    !=====================================================================
    ! Custom PSD: numerical integration for background
    !=====================================================================

    subroutine custom_rhopres(psd, am, rhonu, pnu, drhonu_dam, dpnu_dam)
    type(TCustomNuPSD), intent(in) :: psd
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu, pnu, drhonu_dam, dpnu_dam
    real(dl) :: q, dq_loc, aq, inv_v, v, f0
    real(dl) :: rho_prev, p_prev, drho_prev, dp_prev
    real(dl) :: rho_cur, p_cur, drho_cur, dp_cur
    integer :: i

    rhonu = 0._dl
    pnu = 0._dl
    drhonu_dam = 0._dl
    dpnu_dam = 0._dl

    q = psd%q_tab(1)
    f0 = psd%f0_tab(1)
    aq = am / q
    inv_v = sqrt(1._dl + aq*aq)
    v = 1._dl / inv_v
    rho_prev = q**3 * f0 * inv_v
    p_prev = q**3 * f0 * v
    drho_prev = q * f0 * am * v
    dp_prev = -q * f0 * am * v**3

    do i = 2, psd%nq_tab
        q = psd%q_tab(i)
        f0 = psd%f0_tab(i)
        aq = am / q
        inv_v = sqrt(1._dl + aq*aq)
        v = 1._dl / inv_v
        rho_cur = q**3 * f0 * inv_v
        p_cur = q**3 * f0 * v
        drho_cur = q * f0 * am * v
        dp_cur = -q * f0 * am * v**3

        dq_loc = psd%q_tab(i) - psd%q_tab(i-1)
        rhonu = rhonu + 0.5_dl * dq_loc * (rho_prev + rho_cur)
        pnu = pnu + 0.5_dl * dq_loc * (p_prev + p_cur)
        drhonu_dam = drhonu_dam + 0.5_dl * dq_loc * (drho_prev + drho_cur)
        dpnu_dam = dpnu_dam + 0.5_dl * dq_loc * (dp_prev + dp_cur)

        rho_prev = rho_cur
        p_prev = p_cur
        drho_prev = drho_cur
        dp_prev = dp_cur
    end do

    rhonu = rhonu / fermi_dirac_const
    pnu = pnu / fermi_dirac_const / 3._dl
    drhonu_dam = drhonu_dam / fermi_dirac_const
    dpnu_dam = dpnu_dam / fermi_dirac_const / 3._dl

    end subroutine custom_rhopres

    subroutine custom_ppiv2(psd, am, ppnu, iv2nu, dppnu_dam, div2nu_dam)
    type(TCustomNuPSD), intent(in) :: psd
    real(dl), intent(in) :: am
    real(dl), intent(out) :: ppnu, iv2nu, dppnu_dam, div2nu_dam
    real(dl) :: q, dq_loc, aq, v, v2, f0
    real(dl) :: pp_prev, iv2_prev, dpp_prev, div2_prev
    real(dl) :: pp_cur, iv2_cur, dpp_cur, div2_cur
    integer :: i

    ppnu = 0._dl
    iv2nu = 0._dl
    dppnu_dam = 0._dl
    div2nu_dam = 0._dl

    q = psd%q_tab(1)
    f0 = psd%f0_tab(1)
    aq = am / q
    v2 = 1._dl / (1._dl + aq*aq)
    v = sqrt(v2)
    pp_prev = q**3 * f0 * v * v2
    iv2_prev = q**3 * f0 * v2
    dpp_prev = -3._dl * am * q * f0 * v2*v2*v
    div2_prev = -2._dl * am * q * f0 * v2*v2

    do i = 2, psd%nq_tab
        q = psd%q_tab(i)
        f0 = psd%f0_tab(i)
        aq = am / q
        v2 = 1._dl / (1._dl + aq*aq)
        v = sqrt(v2)
        pp_cur = q**3 * f0 * v * v2
        iv2_cur = q**3 * f0 * v2
        dpp_cur = -3._dl * am * q * f0 * v2*v2*v
        div2_cur = -2._dl * am * q * f0 * v2*v2

        dq_loc = psd%q_tab(i) - psd%q_tab(i-1)
        ppnu = ppnu + 0.5_dl * dq_loc * (pp_prev + pp_cur)
        iv2nu = iv2nu + 0.5_dl * dq_loc * (iv2_prev + iv2_cur)
        dppnu_dam = dppnu_dam + 0.5_dl * dq_loc * (dpp_prev + dpp_cur)
        div2nu_dam = div2nu_dam + 0.5_dl * dq_loc * (div2_prev + div2_cur)

        pp_prev = pp_cur
        iv2_prev = iv2_cur
        dpp_prev = dpp_cur
        div2_prev = div2_cur
    end do

    ppnu = ppnu / fermi_dirac_const / 3._dl
    iv2nu = iv2nu / fermi_dirac_const
    dppnu_dam = dppnu_dam / fermi_dirac_const / 3._dl
    div2nu_dam = div2nu_dam / fermi_dirac_const

    end subroutine custom_ppiv2

    !=====================================================================
    ! TCustomNuPSD type-bound methods
    !=====================================================================

    subroutine TCustomNuPSD_Init(this, nu_q, nu_int_kernel_std, nqmax)
    class(TCustomNuPSD), intent(inout) :: this
    real(dl), intent(in) :: nu_q(:), nu_int_kernel_std(:)
    integer, intent(in) :: nqmax
    real(dl) :: norm, rescale, q, f0, dlnf0, f0_fd, dlnf0_fd, R_ratio
    real(dl) :: am, rhonu, pnu, drhonu_dam, dpnu_dam
    real(dl) :: ppnu, iv2nu, dppnu_dam, div2nu_dam
    real(dl) :: dq_loc
    integer :: i, j

    if (.not. this%active .or. this%nq_tab < 2) return

    ! 1. Normalize f0 so that int q^3 f0 dq = fermi_dirac_const
    norm = 0._dl
    do i = 1, this%nq_tab - 1
        dq_loc = this%q_tab(i+1) - this%q_tab(i)
        norm = norm + 0.5_dl * dq_loc * &
            (this%q_tab(i)**3 * this%f0_tab(i) + this%q_tab(i+1)**3 * this%f0_tab(i+1))
    end do
    if (norm > 0._dl) then
        rescale = fermi_dirac_const / norm
        this%f0_tab = this%f0_tab * rescale
    end if

    ! 2. Compute kernel weights at NuPerturbations q-points
    if (allocated(this%nu_int_kernel)) deallocate(this%nu_int_kernel)
    allocate(this%nu_int_kernel(nqmax))

    do i = 1, nqmax
        q = nu_q(i)
        f0 = custom_f0_at_q(this, q)
        dlnf0 = custom_dlnf0dlnq(this, q)
        f0_fd = 1._dl / (exp(q) + 1._dl)
        dlnf0_fd = -q / (1._dl + exp(-q))
        if (abs(f0_fd * dlnf0_fd) > 1e-30_dl .and. abs(f0 * dlnf0) > 1e-30_dl) then
            R_ratio = (f0 * abs(dlnf0)) / (f0_fd * abs(dlnf0_fd))
        else
            R_ratio = 0._dl
        end if
        this%nu_int_kernel(i) = nu_int_kernel_std(i) * R_ratio
    end do

    ! 3. Compute background spline tables
    if (allocated(this%r1)) deallocate(this%r1, this%p1, this%dr1, this%dp1)
    if (allocated(this%pp1)) deallocate(this%pp1, this%dpp1, this%iv21, this%div21)
    allocate(this%r1(nrhopn), this%p1(nrhopn), this%dr1(nrhopn), this%dp1(nrhopn))
    allocate(this%pp1(nrhopn), this%dpp1(nrhopn), this%iv21(nrhopn), this%div21(nrhopn))
    this%dam_bg = (am_max - am_min) / (nrhopn - 1)

    do j = 1, nrhopn
        am = am_min + (j - 1) * this%dam_bg
        call custom_rhopres(this, am, rhonu, pnu, drhonu_dam, dpnu_dam)
        call custom_ppiv2(this, am, ppnu, iv2nu, dppnu_dam, div2nu_dam)
        this%r1(j) = rhonu
        this%p1(j) = pnu
        this%dr1(j) = drhonu_dam * this%dam_bg
        this%dp1(j) = dpnu_dam * this%dam_bg
        this%pp1(j) = ppnu
        this%iv21(j) = iv2nu
        this%dpp1(j) = dppnu_dam * this%dam_bg
        this%div21(j) = div2nu_dam * this%dam_bg
    end do

    end subroutine TCustomNuPSD_Init

    subroutine TCustomNuPSD_rho_P(this, am, rhonu, pnu)
    class(TCustomNuPSD) :: this
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu, pnu
    real(dl) :: d, dum1, dum2
    integer :: i

    if (am <= am_minp .or. am >= am_maxp) then
        call custom_rhopres(this, am, rhonu, pnu, dum1, dum2)
        return
    end if

    d = (am - am_min) / this%dam_bg + 1._dl
    i = int(d)
    d = d - i
    rhonu = this%r1(i) + d*(this%dr1(i) + d*(3._dl*(this%r1(i+1)-this%r1(i))-2._dl*this%dr1(i) &
        - this%dr1(i+1) + d*(this%dr1(i)+this%dr1(i+1)+2._dl*(this%r1(i)-this%r1(i+1)))))
    pnu = this%p1(i) + d*(this%dp1(i) + d*(3._dl*(this%p1(i+1)-this%p1(i))-2._dl*this%dp1(i) &
        - this%dp1(i+1) + d*(this%dp1(i)+this%dp1(i+1)+2._dl*(this%p1(i)-this%p1(i+1)))))

    end subroutine TCustomNuPSD_rho_P

    subroutine TCustomNuPSD_rho(this, am, rhonu)
    class(TCustomNuPSD) :: this
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu
    real(dl) :: pnu
    call this%rho_P(am, rhonu, pnu)
    end subroutine TCustomNuPSD_rho

    function TCustomNuPSD_drho(this, am, adotoa) result(rhonudot)
    class(TCustomNuPSD) :: this
    real(dl), intent(in) :: am, adotoa
    real(dl) :: rhonudot, rhonu, pnu
    call this%rho_P(am, rhonu, pnu)
    rhonudot = (rhonu - 3._dl*pnu)*adotoa
    end function TCustomNuPSD_drho

    function TCustomNuPSD_ppseudo(this, am) result(ppnu)
    class(TCustomNuPSD) :: this
    real(dl), intent(in) :: am
    real(dl) :: ppnu
    real(dl) :: d, dum1, dum2, dum3
    integer :: i

    if (am <= am_minp .or. am >= am_maxp) then
        call custom_ppiv2(this, am, ppnu, dum1, dum2, dum3)
        return
    end if

    d = (am - am_min) / this%dam_bg + 1._dl
    i = int(d)
    d = d - i
    ppnu = this%pp1(i) + d*(this%dpp1(i) + d*(3._dl*(this%pp1(i+1)-this%pp1(i))-2._dl*this%dpp1(i) &
        - this%dpp1(i+1) + d*(this%dpp1(i)+this%dpp1(i+1)+2._dl*(this%pp1(i)-this%pp1(i+1)))))

    end function TCustomNuPSD_ppseudo

    function TCustomNuPSD_Iv2(this, am) result(iv2nu)
    class(TCustomNuPSD) :: this
    real(dl), intent(in) :: am
    real(dl) :: iv2nu
    real(dl) :: d, dum1, dum2, dum3
    integer :: i

    if (am <= am_minp .or. am >= am_maxp) then
        call custom_ppiv2(this, am, dum1, iv2nu, dum2, dum3)
        return
    end if

    d = (am - am_min) / this%dam_bg + 1._dl
    i = int(d)
    d = d - i
    iv2nu = this%iv21(i) + d*(this%div21(i) + d*(3._dl*(this%iv21(i+1)-this%iv21(i))-2._dl*this%div21(i) &
        - this%div21(i+1) + d*(this%div21(i)+this%div21(i+1)+2._dl*(this%iv21(i)-this%iv21(i+1)))))

    end function TCustomNuPSD_Iv2

    function custom_rho_err(this, nu_mass)
    class(TCustomNuPSD) :: this
    real(dl) :: custom_rho_err, nu_mass, rhonu
    call this%rho(nu_mass, rhonu)
    custom_rho_err = rhonu - this%target_rho
    end function custom_rho_err

    function TCustomNuPSD_find_nu_mass(this, rho) result(nu_mass)
    use MathUtils
    use config
    class(TCustomNuPSD) :: this
    real(dl), intent(in) :: rho
    real(dl) :: nu_mass, fzero
    integer :: iflag

    if (rho <= 1.001_dl) then
        nu_mass = 0
        return
    end if
    nu_mass = fermi_dirac_const / (1.5d0*zeta3) * rho
    this%target_rho = rho
    call brentq(this, custom_rho_err, 0._dl, nu_mass*2, 0.01_dl, nu_mass, fzero, iflag)
    if (iflag /= 0) call GlobalError('TCustomNuPSD: find_nu_mass_for_rho failed')

    end function TCustomNuPSD_find_nu_mass

    end module MassiveNu
