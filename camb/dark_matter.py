"""
Dark matter interaction models for CAMB.

Implements:
- DM-Baryon scattering (Boddy & Gluscevic 2018)
- DM-Dark Radiation (ETHOS framework, Cyr-Racine et al. 2016)
- Decaying Dark Matter (Blackadder & Koushiappas 2014, CLASS dcdm)
- DM-Neutrino scattering (Mangano+ 2006, He+ 2023)
- Warm Dark Matter (Viel+ 2005, transfer function approach)
- Fuzzy/Ultralight Axion Dark Matter (Hu+ 2000, Hlozek+ 2015)
"""

from ctypes import c_bool, c_double, c_int

from .baseconfig import F2003Class, fortran_class

max_ethos_coeffs = 6


class DarkMatterModel(F2003Class):
    """
    Abstract base class for dark matter interaction model implementations.
    """

    _fields_ = (
        ("__is_standard_cdm", c_bool),
        ("__num_perturb_equations", c_int),
        ("__has_cdm_velocity", c_bool),
        ("__num_dr_equations", c_int),
    )

    def validate_params(self) -> None:
        pass


@fortran_class
class DMBaryonScattering(DarkMatterModel):
    """
    DM-Baryon scattering model with momentum-transfer cross section
    sigma_MT = sigma_0 * (v/c)^n.

    Following Boddy & Gluscevic (2018).

    Parameters:
        sigma_dmb: momentum-transfer cross section sigma_0 [cm^2]
        n_dmb: velocity power-law index (-4, -2, 0, 2, 4)
        m_dm: DM particle mass [GeV]
    """

    _fields_ = (
        ("sigma_dmb", c_double, "Momentum-transfer cross section sigma_0 [cm^2]"),
        ("n_dmb", c_int, "Velocity power-law index (-4,-2,0,2,4)"),
        ("m_dm", c_double, "DM particle mass [GeV]"),
    )

    _fortran_class_module_ = "DMBaryon"
    _fortran_class_name_ = "TDMBaryonScattering"

    def set_params(self, sigma_dmb=0.0, n_dmb=0, m_dm=100.0):
        """
        Set DM-baryon scattering parameters.

        :param sigma_dmb: momentum-transfer cross section [cm^2]
        :param n_dmb: velocity dependence index (-4,-2,0,2,4)
        :param m_dm: DM particle mass [GeV]
        """
        self.sigma_dmb = sigma_dmb
        self.n_dmb = n_dmb
        self.m_dm = m_dm
        self.validate_params()

    def validate_params(self):
        if self.n_dmb not in (-4, -2, 0, 2, 4):
            from .baseconfig import CAMBError

            raise CAMBError(f"n_dmb must be -4, -2, 0, 2, or 4, got {self.n_dmb}")
        if self.sigma_dmb < 0:
            from .baseconfig import CAMBError

            raise CAMBError("sigma_dmb must be non-negative")
        if self.m_dm <= 0:
            from .baseconfig import CAMBError

            raise CAMBError("m_dm must be positive")


@fortran_class
class DMDR_ETHOS(DarkMatterModel):
    """
    DM-Dark Radiation interaction model in the ETHOS framework.

    DM is coupled to a dark radiation sector via scattering.
    DM: fluid equations (delta, v)
    DR: Boltzmann hierarchy (F_0, F_1, ..., F_lmax)

    Following Cyr-Racine et al. (2016).

    Parameters:
        omdmdrh2: interacting DM density Omega_dmdr h^2
        N_dark: dark radiation effective number (Delta N_eff in dark sector)
        a_dark: ETHOS opacity coefficients [a_2, a_3, a_4, a_5, a_6, reserved]
        alpha_l_dark: angular damping coefficient for l>=2
        lmax_dr: DR Boltzmann hierarchy truncation
        cs2_dm: DM effective sound speed squared
    """

    _fields_ = (
        ("omdmdrh2", c_double, "Interacting DM density Omega_dmdr h^2"),
        ("N_dark", c_double, "Dark radiation effective number"),
        ("a_dark", c_double * max_ethos_coeffs, "ETHOS opacity coefficients"),
        ("alpha_l_dark", c_double, "Angular damping coefficient for l>=2"),
        ("lmax_dr", c_int, "DR Boltzmann hierarchy truncation"),
        ("cs2_dm", c_double, "DM effective sound speed squared"),
    )

    _fortran_class_module_ = "DMDR_ETHOS"
    _fortran_class_name_ = "TDMDR_ETHOS"

    def set_params(
        self,
        omdmdrh2=0.0,
        N_dark=0.0,
        a_dark=None,
        alpha_l_dark=1.0,
        lmax_dr=15,
        cs2_dm=0.0,
    ):
        """
        Set DM-DR ETHOS parameters.

        :param omdmdrh2: interacting DM density Omega_dmdr h^2
        :param N_dark: dark radiation effective number
        :param a_dark: list/array of ETHOS opacity coefficients [a_2, ..., a_6]
        :param alpha_l_dark: angular damping coefficient for l>=2
        :param lmax_dr: DR hierarchy truncation (default 15)
        :param cs2_dm: DM effective sound speed squared
        """
        self.omdmdrh2 = omdmdrh2
        self.N_dark = N_dark
        if a_dark is not None:
            for i, val in enumerate(a_dark):
                if i < max_ethos_coeffs:
                    self.a_dark[i] = val
        self.alpha_l_dark = alpha_l_dark
        self.lmax_dr = lmax_dr
        self.cs2_dm = cs2_dm
        self.validate_params()

    def validate_params(self):
        if self.omdmdrh2 < 0:
            from .baseconfig import CAMBError

            raise CAMBError("omdmdrh2 must be non-negative")
        if self.N_dark < 0:
            from .baseconfig import CAMBError

            raise CAMBError("N_dark must be non-negative")
        if self.lmax_dr < 2:
            from .baseconfig import CAMBError

            raise CAMBError("lmax_dr must be >= 2")


@fortran_class
class DecayingDM(DarkMatterModel):
    """
    Decaying Dark Matter model.

    Parent DM decays into massive daughter + dark radiation:
    DM_parent -> DM_daughter (mass fraction epsilon) + DR

    Background: rho_dcdm(a) = rho_0 * a^{-3} * exp(-Gamma*t(a))

    References: Blackadder & Koushiappas 2014, Poulin+ 2016, CLASS dcdm module.

    Parameters:
        Gamma_dcdm: decay rate [km/s/Mpc]
        epsilon_dcdm: mass ratio m_daughter/m_parent (0 < epsilon <= 1)
        f_dcdm: fraction of CDM that is decaying (0 to 1)
    """

    _fields_ = (
        ("Gamma_dcdm", c_double, "Decay rate [km/s/Mpc]"),
        ("epsilon_dcdm", c_double, "Mass ratio m_daughter/m_parent"),
        ("f_dcdm", c_double, "Fraction of CDM that is decaying"),
    )

    _fortran_class_module_ = "DecayingDM"
    _fortran_class_name_ = "TDecayingDM"

    def set_params(self, Gamma_dcdm=0.0, epsilon_dcdm=1.0, f_dcdm=1.0):
        """
        Set decaying DM parameters.

        :param Gamma_dcdm: decay rate [km/s/Mpc]
        :param epsilon_dcdm: daughter/parent mass ratio (0 < eps <= 1)
        :param f_dcdm: fraction of CDM that decays
        """
        self.Gamma_dcdm = Gamma_dcdm
        self.epsilon_dcdm = epsilon_dcdm
        self.f_dcdm = f_dcdm
        self.validate_params()

    def validate_params(self):
        if self.Gamma_dcdm < 0:
            from .baseconfig import CAMBError

            raise CAMBError("Gamma_dcdm must be non-negative")
        if not (0 < self.epsilon_dcdm <= 1):
            from .baseconfig import CAMBError

            raise CAMBError("epsilon_dcdm must be in (0, 1]")
        if not (0 <= self.f_dcdm <= 1):
            from .baseconfig import CAMBError

            raise CAMBError("f_dcdm must be in [0, 1]")


@fortran_class
class DMNeutrinoScattering(DarkMatterModel):
    """
    DM-Neutrino scattering model.

    DM interacts with neutrinos: sigma ~ sigma_0 * (E_nu/1MeV)^n.
    Adds drag in CDM velocity equation and collision damping in neutrino hierarchy.

    References: Mangano+ 2006, Wilkinson+ 2014, He+ 2023 (Nature Astronomy 2025).

    Parameters:
        sigma_dmnu: cross section at E_0=1 MeV [cm^2]
        n_dmnu: energy power-law index (typically 0 or 2)
        m_dm: DM particle mass [GeV]
    """

    _fields_ = (
        ("sigma_dmnu", c_double, "DM-nu cross section at 1 MeV [cm^2]"),
        ("n_dmnu", c_int, "Energy dependence index"),
        ("m_dm", c_double, "DM particle mass [GeV]"),
    )

    _fortran_class_module_ = "DMNeutrino"
    _fortran_class_name_ = "TDMNeutrinoScattering"

    def set_params(self, sigma_dmnu=0.0, n_dmnu=0, m_dm=100.0):
        """
        Set DM-neutrino scattering parameters.

        :param sigma_dmnu: cross section at E_0=1 MeV [cm^2]
        :param n_dmnu: energy dependence index (0 or 2)
        :param m_dm: DM mass [GeV]
        """
        self.sigma_dmnu = sigma_dmnu
        self.n_dmnu = n_dmnu
        self.m_dm = m_dm
        self.validate_params()

    def validate_params(self):
        if self.sigma_dmnu < 0:
            from .baseconfig import CAMBError

            raise CAMBError("sigma_dmnu must be non-negative")
        if self.m_dm <= 0:
            from .baseconfig import CAMBError

            raise CAMBError("m_dm must be positive")


@fortran_class
class WarmDM(DarkMatterModel):
    """
    Warm Dark Matter model.

    Two modes:
    1. Transfer function (default): T(k) = [1 + (alpha*k)^(2*nu)]^(-5/nu)
    2. Full Boltzmann (boltzmann_mode=True): WDM registered as extra massive
       neutrino species, using CAMB's neutrino Boltzmann hierarchy.
       Exact for T_wdm_ratio=1; approximate for other values.

    References: Bode+ 2001, Viel+ 2005, Lesgourgues & Tram 2011.

    Parameters:
        m_wdm: WDM particle mass [keV]
        T_wdm_ratio: temperature ratio T_WDM/T_nu (default 1)
        Omega_wdm_h2: WDM density (0 = use all of omch2)
        boltzmann_mode: use full Boltzmann instead of transfer function
    """

    _fields_ = (
        ("m_wdm", c_double, "WDM mass [keV]"),
        ("T_wdm_ratio", c_double, "T_WDM/T_nu ratio"),
        ("Omega_wdm_h2", c_double, "WDM density Omega_wdm h^2"),
        ("boltzmann_mode", c_bool, "Full Boltzmann mode"),
    )

    _fortran_class_module_ = "WarmDM"
    _fortran_class_name_ = "TWarmDM"

    def set_params(self, m_wdm=3.0, T_wdm_ratio=1.0, Omega_wdm_h2=0.0, boltzmann_mode=False):
        """
        Set warm DM parameters.

        :param m_wdm: WDM mass in keV
        :param T_wdm_ratio: temperature ratio T_WDM/T_nu
        :param Omega_wdm_h2: WDM density (0 uses all omch2)
        :param boltzmann_mode: True for full Boltzmann, False for transfer function
        """
        self.m_wdm = m_wdm
        self.T_wdm_ratio = T_wdm_ratio
        self.Omega_wdm_h2 = Omega_wdm_h2
        self.boltzmann_mode = boltzmann_mode
        self.validate_params()

    def validate_params(self):
        if self.m_wdm <= 0:
            from .baseconfig import CAMBError

            raise CAMBError("m_wdm must be positive")
        if self.T_wdm_ratio <= 0:
            from .baseconfig import CAMBError

            raise CAMBError("T_wdm_ratio must be positive")


@fortran_class
class FuzzyDM(DarkMatterModel):
    """
    Fuzzy/Ultralight Axion Dark Matter (ULDM).

    Ultra-light boson DM (m ~ 10^-22 eV) with astrophysically large
    de Broglie wavelength. Effective fluid with scale-dependent sound speed
    cs2(k,a) = k^2 / (4*m^2*a^2) producing a quantum Jeans scale.

    References: Hu+ 2000, Hlozek+ 2015, axionCAMB.

    Parameters:
        m_axion: axion mass [eV]
        omega_axion_h2: axion density Omega_a h^2 (0 = use f_axion)
        f_axion: fraction of CDM that is ULDM
    """

    _fields_ = (
        ("m_axion", c_double, "Axion mass [eV]"),
        ("omega_axion_h2", c_double, "Axion density Omega_a h^2"),
        ("f_axion", c_double, "Fraction of CDM that is axion"),
    )

    _fortran_class_module_ = "FuzzyDM"
    _fortran_class_name_ = "TFuzzyDM"

    def set_params(self, m_axion=1e-22, omega_axion_h2=0.0, f_axion=1.0):
        """
        Set fuzzy DM parameters.

        :param m_axion: axion mass [eV]
        :param omega_axion_h2: axion density (0 uses f_axion * omch2)
        :param f_axion: fraction of CDM that is ULDM
        """
        self.m_axion = m_axion
        self.omega_axion_h2 = omega_axion_h2
        self.f_axion = f_axion
        self.validate_params()

    def validate_params(self):
        if self.m_axion <= 0:
            from .baseconfig import CAMBError

            raise CAMBError("m_axion must be positive")
        if self.f_axion < 0 or self.f_axion > 1:
            from .baseconfig import CAMBError

            raise CAMBError("f_axion must be in [0, 1]")


@fortran_class
class TransientGDM(DarkMatterModel):
    """Transient early-dark-energy-like generalized dark matter.

    A fraction ``f_X`` of the total density supplied through ``omch2`` follows
    the Gaussian equation-of-state dip and scale-dependent rest-frame sound
    speed of Bhattacharjee et al. (2026), arXiv:2608.20288. The component
    returns to pressureless matter after the transient while retaining its own
    density and velocity perturbations.
    """

    _fields_ = (
        ("f_X", c_double, "Present fraction of total omch2 assigned to X"),
        ("w_p", c_double, "Minimum transient equation of state"),
        ("kappa0_hmpc", c_double, "Comoving regulator scale at a_f [h/Mpc]"),
        ("a_i", c_double, "Start of the +/-3-sigma transient interval"),
        ("a_f", c_double, "End of the +/-3-sigma transient interval"),
    )

    _fortran_class_module_ = "TransientGDM"
    _fortran_class_name_ = "TTransientGDM"

    def set_params(
        self,
        f_X=0.00155,
        w_p=-0.45,
        kappa0_hmpc=6.4,
        a_i=1.46e-7,
        a_f=1.46e-5,
    ):
        self.f_X = f_X
        self.w_p = w_p
        self.kappa0_hmpc = kappa0_hmpc
        self.a_i = a_i
        self.a_f = a_f
        self.validate_params()
        return self

    def validate_params(self):
        from .baseconfig import CAMBError

        if not 0 <= self.f_X < 1:
            raise CAMBError("TransientGDM f_X must be in [0, 1)")
        if not -1 < self.w_p <= 0:
            raise CAMBError("TransientGDM w_p must be in (-1, 0]")
        if self.kappa0_hmpc <= 0:
            raise CAMBError("TransientGDM kappa0_hmpc must be positive")
        if not 0 < self.a_i < self.a_f < 1:
            raise CAMBError("TransientGDM requires 0 < a_i < a_f < 1")


@fortran_class
class MultiInteractingDM(DarkMatterModel):
    """
    Multi-Interacting DM: single DM fluid with up to 3 simultaneous channels.
    - DR (ETHOS): DM-Dark Radiation Boltzmann hierarchy
    - Baryon: DM-Baryon momentum transfer
    - Photon: DM-Photon opacity (Boddy & Gluscevic 2018)

    Reference: Becker+ 2021 (arXiv:2010.04074)
    """

    _fields_ = [
        ("omdmdrh2", c_double),
        ("N_dark", c_double),
        ("a_dark", c_double * 6),
        ("alpha_l_dark", c_double),
        ("lmax_dr", c_int),
        ("sigma_dmb", c_double),
        ("n_dmb", c_int),
        ("u_idm_g", c_double),
        ("n_idm_g", c_int),
        ("m_dm", c_double),
        ("cs2_dm", c_double),
    ]
    _fortran_class_module_ = "MultiInteractingDM"
    _fortran_class_name_ = "TMultiInteractingDM"

    def set_params(
        self,
        omdmdrh2=0.0,
        N_dark=0.0,
        a_dark=None,
        alpha_l_dark=1.0,
        lmax_dr=15,
        sigma_dmb=0.0,
        n_dmb=0,
        u_idm_g=0.0,
        n_idm_g=0,
        m_dm=100.0,
        cs2_dm=0.0,
    ):
        self.omdmdrh2 = omdmdrh2
        self.N_dark = N_dark
        if a_dark is not None:
            for i, v in enumerate(a_dark[:6]):
                self.a_dark[i] = v
        self.alpha_l_dark = alpha_l_dark
        self.lmax_dr = lmax_dr
        self.sigma_dmb = sigma_dmb
        self.n_dmb = n_dmb
        self.u_idm_g = u_idm_g
        self.n_idm_g = n_idm_g
        self.m_dm = m_dm
        self.cs2_dm = cs2_dm

    def validate_params(self):
        if self.m_dm <= 0:
            from .baseconfig import CAMBError

            raise CAMBError("m_dm must be positive")


@fortran_class
class ETHOSTransferMurgia(DarkMatterModel):
    """
    Murgia et al. 2017 phenomenological ETHOS transfer function applied as
    post-processing on the matter power spectrum:

        T(k) = [1 + (alpha*k)^beta]^gamma

    P(k) is multiplied by T(k)^2. CMB and background are unchanged. For full
    self-consistent CMB + matter + dark acoustic oscillations use
    :class:`DMDR_ETHOS` instead.

    Reference: Murgia, Merle, Viel, Totzauer & Schneider (arXiv:1704.07838).

    Parameters:
        alpha_mpch: free-streaming scale alpha [Mpc/h]
        beta_mur: shape parameter beta (default 2.24)
        gamma_mur: cutoff slope gamma (default -4.46)
    """

    _fields_ = [
        ("alpha_mpch", c_double, "Free-streaming scale alpha [Mpc/h]"),
        ("beta_mur", c_double, "Shape parameter beta"),
        ("gamma_mur", c_double, "Cutoff slope gamma"),
    ]
    _fortran_class_module_ = "ETHOSTransferMurgia"
    _fortran_class_name_ = "TETHOSTransferMurgia"

    def set_params(self, alpha_mpch=0.0, beta_mur=2.24, gamma_mur=-4.46):
        """
        Set Murgia ETHOS transfer function parameters.

        :param alpha_mpch: free-streaming scale [Mpc/h] (0 = LCDM)
        :param beta_mur: shape parameter
        :param gamma_mur: cutoff sharpness
        """
        self.alpha_mpch = alpha_mpch
        self.beta_mur = beta_mur
        self.gamma_mur = gamma_mur
        self.validate_params()

    def validate_params(self):
        if self.alpha_mpch < 0:
            from .baseconfig import CAMBError

            raise CAMBError("alpha_mpch must be non-negative")
        if self.beta_mur <= 0:
            from .baseconfig import CAMBError

            raise CAMBError("beta_mur must be positive")


@fortran_class
class ETHOSTransferPhysical(DarkMatterModel):
    """
    Physics-parameterized ETHOS transfer-function surrogate.

    Maps physical ETHOS-like parameters (a_dark_n, xi_dr, omega_dmdr_h2)
    to the Murgia 2017 transfer function via a power-law fit:

        alpha [Mpc/h] = A0 * (a_dark_n / a_ref)^p_a * (xi_dr / xi_ref)^p_xi
                          * (omdmdrh2 / 0.12)^p_om

    then P(k) is multiplied by [1 + (alpha k)^beta]^(2 gamma).

    The fit constants (A0_fit, a_ref, xi_ref, p_a, p_xi, p_om) are exposed
    so each application can recalibrate the surrogate against the full
    :class:`DMDR_ETHOS` Boltzmann solution for its parameter ranges.
    Defaults give an order-of-magnitude correct mapping for ETHOS
    n_dark=2.

    This is a fast surrogate for parameter scans. It does not reproduce
    dark acoustic oscillations or CMB effects -- use :class:`DMDR_ETHOS`
    for the full self-consistent calculation.

    Parameters:
        a_dark_n: ETHOS opacity coefficient [Mpc^-1]
        xi_dr: dark-photon temperature ratio T_DR / T_gamma
        omdmdrh2: interacting DM density Omega_dmdr h^2
        n_dark: opacity index (advisory only; 2 = dark-Thomson)
        A0_fit: alpha [Mpc/h] at the reference point
        a_ref, xi_ref: reference parameter values for the fit
        p_a, p_xi, p_om: power-law exponents
        beta_mur, gamma_mur: Murgia shape parameters
    """

    _fields_ = [
        ("a_dark_n", c_double, "ETHOS opacity coefficient [Mpc^-1]"),
        ("xi_dr", c_double, "T_DR / T_gamma at recomb"),
        ("omdmdrh2", c_double, "Omega_dmdr h^2"),
        ("n_dark", c_int, "Opacity power-law index"),
        ("A0_fit", c_double, "alpha [Mpc/h] at reference point"),
        ("a_ref", c_double, "Reference a_dark_n [Mpc^-1]"),
        ("xi_ref", c_double, "Reference xi_dr"),
        ("p_a", c_double, "Power on a_dark_n"),
        ("p_xi", c_double, "Power on xi_dr"),
        ("p_om", c_double, "Power on omega_dmdr_h2"),
        ("beta_mur", c_double, "Murgia beta"),
        ("gamma_mur", c_double, "Murgia gamma"),
        ("alpha_mpch", c_double, "Derived alpha [Mpc/h] (output)"),
    ]
    _fortran_class_module_ = "ETHOSTransferPhysical"
    _fortran_class_name_ = "TETHOSTransferPhysical"

    def set_params(
        self,
        a_dark_n=0.0,
        xi_dr=0.5,
        omdmdrh2=0.12,
        n_dark=2,
        A0_fit=0.5,
        a_ref=1e-4,
        xi_ref=0.5,
        p_a=0.25,
        p_xi=1.0,
        p_om=-0.1,
        beta_mur=2.24,
        gamma_mur=-4.46,
    ):
        """
        Set ETHOS physical-surrogate parameters. The derived ``alpha_mpch``
        is computed here and stored, so it can be read back from the
        instance after :meth:`set_params`.

        :param a_dark_n: ETHOS opacity coefficient [Mpc^-1] (0 = LCDM)
        :param xi_dr: dark-photon temperature ratio
        :param omdmdrh2: Omega_dmdr h^2
        :param n_dark: opacity power-law index
        :param A0_fit, a_ref, xi_ref, p_a, p_xi, p_om: surrogate fit constants
        :param beta_mur, gamma_mur: Murgia shape parameters
        """
        self.a_dark_n = a_dark_n
        self.xi_dr = xi_dr
        self.omdmdrh2 = omdmdrh2
        self.n_dark = n_dark
        self.A0_fit = A0_fit
        self.a_ref = a_ref
        self.xi_ref = xi_ref
        self.p_a = p_a
        self.p_xi = p_xi
        self.p_om = p_om
        self.beta_mur = beta_mur
        self.gamma_mur = gamma_mur
        if a_dark_n > 0 and xi_dr > 0:
            self.alpha_mpch = A0_fit * (a_dark_n / a_ref) ** p_a * (xi_dr / xi_ref) ** p_xi * (omdmdrh2 / 0.12) ** p_om
        else:
            self.alpha_mpch = 0.0
        self.validate_params()

    def validate_params(self):
        from .baseconfig import CAMBError

        if self.a_dark_n < 0:
            raise CAMBError("a_dark_n must be non-negative")
        if self.xi_dr <= 0:
            raise CAMBError("xi_dr must be positive")
        if self.beta_mur <= 0:
            raise CAMBError("beta_mur must be positive")


@fortran_class
class DMPhotonScattering(DarkMatterModel):
    """
    DM-Photon scattering model (Boddy & Gluscevic 2018).
    sigma = sigma_Th * u_idm_g * (m_DM / 100 GeV) * (T/T_0)^n_idm_g
    """

    _fields_ = [
        ("u_idm_g", c_double),
        ("n_idm_g", c_int),
        ("m_dm", c_double),
    ]
    _fortran_class_module_ = "DMPhoton"
    _fortran_class_name_ = "TDMPhotonScattering"

    def set_params(self, u_idm_g=0.0, n_idm_g=0, m_dm=100.0):
        self.u_idm_g = u_idm_g
        self.n_idm_g = n_idm_g
        self.m_dm = m_dm
        self.validate_params()

    def validate_params(self):
        if self.u_idm_g < 0:
            from .baseconfig import CAMBError

            raise CAMBError("u_idm_g must be non-negative")
        if self.m_dm <= 0:
            from .baseconfig import CAMBError

            raise CAMBError("m_dm must be positive")


# Register short names
F2003Class._class_names.update(
    {
        "dm_baryon": DMBaryonScattering,
        "dm_dr_ethos": DMDR_ETHOS,
        "decaying_dm": DecayingDM,
        "dm_neutrino": DMNeutrinoScattering,
        "warm_dm": WarmDM,
        "fuzzy_dm": FuzzyDM,
        "transient_gdm": TransientGDM,
        "dm_photon": DMPhotonScattering,
        "multi_idm": MultiInteractingDM,
        "ethos_transfer_murgia": ETHOSTransferMurgia,
        "ethos_transfer_physical": ETHOSTransferPhysical,
    }
)
