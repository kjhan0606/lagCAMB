from ctypes import POINTER, byref, c_bool, c_double, c_int

from .baseconfig import AllocatableArrayDouble, CAMBError, F2003Class, f_pointer, fortran_class, np, numpy_1d


class DarkEnergyModel(F2003Class):
    """
    Abstract base class for dark energy model implementations.
    """

    _fields_ = (("__is_cosmological_constant", c_bool), ("__num_perturb_equations", c_int))

    def validate_params(self) -> None:
        pass


class DarkEnergyEqnOfState(DarkEnergyModel):
    """
    Abstract base class for models using w and wa parameterization with use w(a) = w + (1-a)*wa parameterization,
    or call set_w_a_table to set another tabulated w(a). If tabulated w(a) is used, w and wa are set
    to approximate values at z=0.

    """

    _fortran_class_module_ = "DarkEnergyInterface"
    _fortran_class_name_ = "TDarkEnergyEqnOfState"

    _fields_ = (
        ("w", c_double, "w(0)"),
        ("wa", c_double, "-dw/da(0)"),
        ("cs2", c_double, "fluid rest-frame sound speed squared"),
        ("use_tabulated_w", c_bool, "using an interpolated tabulated w(a) rather than w, wa above"),
        ("__no_perturbations", c_bool, "turn off perturbations (unphysical, so hidden in Python)"),
    )

    _methods_ = (("SetWTable", [numpy_1d, numpy_1d, POINTER(c_int)]),)

    def set_params(self, w=-1.0, wa=0, cs2=1.0, use_tabulated_w=False, wde_a_array=None, wde_w_array=None):
        """
         Set the parameters so that P(a)/rho(a) = w(a) = w + (1-a)*wa

        :param w: w(0)
        :param wa: -dw/da(0)
        :param cs2: fluid rest-frame sound speed squared
        :param use_tabulated_w: whether use interpolated w
        :param wde_a_array: array of scale factors
        :param wde_w_array: array of w(a)
        """
        self.use_tabulated_w = use_tabulated_w
        if self.use_tabulated_w:
            if w != -1.0 or wa != 0:
                raise ValueError("cannot use w, wa as well as use_tabulated_w")
            self.set_w_a_table(wde_a_array, wde_w_array)
        else:
            self.w = w
            self.wa = wa
        self.cs2 = cs2
        self.validate_params()

    def validate_params(self):
        if not self.use_tabulated_w and self.wa + self.w > 0:
            raise CAMBError("dark energy model has w + wa > 0, giving w>0 at high redshift")

    def set_w_a_table(self, a, w) -> "DarkEnergyEqnOfState":
        """
        Set w(a) from numerical values (used as cubic spline). Note this is quite slow.

        :param a: array of scale factors
        :param w: array of w(a)
        :return: self
        """
        a = np.ascontiguousarray(a, dtype=np.float64)
        w = np.ascontiguousarray(w, dtype=np.float64)

        if len(a) != len(w):
            raise ValueError("Dark energy w(a) table non-equal sized arrays")
        if not np.isclose(a[-1], 1):
            raise ValueError("Dark energy w(a) arrays must end at a=1")
        if np.any(a <= 0):
            raise ValueError("Dark energy w(a) table cannot be set for a<=0")

        self.f_SetWTable(a, w, byref(c_int(len(a))))
        return self

    def __getstate__(self):
        if self.use_tabulated_w:
            raise TypeError("Cannot save class with splines")
        return super().__getstate__()


@fortran_class
class DarkEnergyFluid(DarkEnergyEqnOfState):
    """
    Class implementing the w, wa or splined w(a) parameterization using the constant sound-speed single fluid model
    (as for single-field quintessence).

    """

    _fortran_class_module_ = "DarkEnergyFluid"
    _fortran_class_name_ = "TDarkEnergyFluid"

    def validate_params(self) -> None:
        super().validate_params()
        if not self.use_tabulated_w and self.wa and (self.w < -1 - 1e-6 or 1 + self.w + self.wa < -1e-6):
            raise CAMBError("fluid dark energy model does not support w crossing -1")

    def set_w_a_table(self, a, w) -> "DarkEnergyEqnOfState":
        # check w array has elements that do not cross -1
        if np.sign(1 + np.max(w)) - np.sign(1 + np.min(w)) == 2:
            raise CAMBError("fluid dark energy model does not support w crossing -1")
        return super().set_w_a_table(a, w)


@fortran_class
class DarkEnergyPPF(DarkEnergyEqnOfState):
    """
    Class implementing the w, wa or splined w(a) parameterization in the PPF perturbation approximation
    (`arXiv:0808.3125 <https://arxiv.org/abs/0808.3125>`_)
    Use inherited methods to set parameters or interpolation table.

    Note PPF is not a physical model and just designed to allow crossing -1 in an ad hoc smooth way. For models
    with w>-1 but far from cosmological constant, it can give quite different answers to the fluid model with c_s^2=1.

    """

    # cannot declare c_Gamma_ppf directly here as have not defined all fields in DarkEnergyEqnOfState (TCubicSpline)
    _fortran_class_module_ = "DarkEnergyPPF"
    _fortran_class_name_ = "TDarkEnergyPPF"


@fortran_class
class AxionEffectiveFluid(DarkEnergyModel):
    """
    Example implementation of a specific (early) dark energy fluid model
    (`arXiv:1806.10608 <https://arxiv.org/abs/1806.10608>`_).
    Not well tested, but should serve to demonstrate how to make your own custom classes.
    """

    _fields_ = (
        ("w_n", c_double, "effective equation of state parameter"),
        ("fde_zc", c_double, "energy density fraction at z=zc"),
        ("zc", c_double, "decay transition redshift (not same as peak of energy density fraction)"),
        ("theta_i", c_double, "initial condition field value"),
    )

    _fortran_class_name_ = "TAxionEffectiveFluid"
    _fortran_class_module_ = "DarkEnergyFluid"

    def set_params(self, w_n, fde_zc, zc, theta_i=None):
        self.w_n = w_n
        self.fde_zc = fde_zc
        self.zc = zc
        if theta_i is not None:
            self.theta_i = theta_i


# base class for scalar field quintessence models
class Quintessence(DarkEnergyModel):
    r"""
    Abstract base class for single scalar field quintessence models.

    For each model the field value and derivative are stored and splined at sampled scale factor values.

    To implement a new model, need to define a new derived class in Fortran,
    defining Vofphi and setting up initial conditions and interpolation tables (see TEarlyQuintessence as example).

    """

    _fields_ = (
        ("DebugLevel", c_int),
        ("astart", c_double),
        ("integrate_tol", c_double),
        ("sampled_a", AllocatableArrayDouble),
        ("phi_a", AllocatableArrayDouble),
        ("phidot_a", AllocatableArrayDouble),
        ("__npoints_linear", c_int),
        ("__npoints_log", c_int),
        ("__dloga", c_double),
        ("__da", c_double),
        ("__log_astart", c_double),
        ("__max_a_log", c_double),
        ("__ddphi_a", AllocatableArrayDouble),
        ("__ddphidot_a", AllocatableArrayDouble),
        ("__state", f_pointer),
    )
    _fortran_class_module_ = "Quintessence"

    def __getstate__(self):
        raise TypeError("Cannot save class with splines")


@fortran_class
class EarlyQuintessence(Quintessence):
    r"""
    Example early quintessence (axion-like, as `arXiv:1908.06995 <https://arxiv.org/abs/1908.06995>`_) with potential

     V(\phi) = m^2f^2 (1 - cos(\phi/f))^n + \Lambda_{cosmological constant}

    """

    _fields_ = (
        ("n", c_double, "power index for potential"),
        ("f", c_double, r"f/Mpl (sqrt(8\piG)f); only used for initial search value when use_zc is True"),
        (
            "m",
            c_double,
            "mass parameter in reduced Planck mass units; only used for initial search value when use_zc is True",
        ),
        ("theta_i", c_double, "phi/f initial field value"),
        ("frac_lambda0", c_double, "fraction of dark energy in cosmological constant today (approximated as 1)"),
        ("use_zc", c_bool, "solve for f, m to get specific critical redshift zc and fde_zc"),
        ("zc", c_double, "redshift of peak fractional early dark energy density"),
        ("fde_zc", c_double, "fraction of early dark energy density to total at peak"),
        ("npoints", c_int, "number of points for background integration spacing"),
        ("min_steps_per_osc", c_int, "minimum number of steps per background oscillation scale"),
        (
            "fde",
            AllocatableArrayDouble,
            "after initialized, the calculated background early dark energy fractions at sampled_a",
        ),
        ("__ddfde", AllocatableArrayDouble),
    )
    _fortran_class_name_ = "TEarlyQuintessence"

    def set_params(self, n, f=0.05, m=5e-54, theta_i=0.0, use_zc=True, zc=None, fde_zc=None):
        self.n = n
        self.f = f
        self.m = m
        self.theta_i = theta_i
        self.use_zc = use_zc
        if use_zc:
            if zc is None or fde_zc is None:
                raise ValueError("must set zc and fde_zc if using 'use_zc'")
            self.zc = zc
            self.fde_zc = fde_zc


@fortran_class
class TrackerQuintessence(Quintessence):
    r"""
    Tracker / thawing quintessence with a single scalar field and no cosmological constant.

    Two potentials, selected by ``pot_type``:

      * ``pot_type=1`` (Ratra-Peebles): :math:`V(\phi) = A\,\phi^{-\alpha}`
      * ``pot_type=2`` (exponential):   :math:`V(\phi) = A\,e^{-\lambda\phi}`

    with :math:`\phi` in reduced Planck units. The amplitude :math:`A` is fixed by shooting
    (bisection on :math:`\ln A`) so that :math:`\Omega_\phi(a=1)` equals the dark-energy budget.
    ``ic_mode=0`` starts the field frozen (:math:`d\phi/dN=0`) at ``phi_ini`` at
    :math:`a=10^{-6}`. ``ic_mode=1`` selects the matter-era Ratra--Peebles tracker
    and computes the initial field from the trial amplitude at every shooting step;
    this is the physical :math:`\phi`CDM path.
    """

    _fields_ = (
        ("pot_type", c_int, "1=Ratra-Peebles (phi^-alpha), 2=exponential (exp(-lam*phi))"),
        ("ic_mode", c_int, "0=frozen phi_ini, 1=matter-era Ratra-Peebles tracker (phiCDM)"),
        ("alpha", c_double, "Ratra-Peebles exponent (V ~ phi^-alpha)"),
        ("lam", c_double, "exponential slope (V ~ exp(-lam*phi))"),
        ("phi_ini", c_double, "initial field value at a=astart (reduced Planck units)"),
        ("lnA", c_double, "log amplitude in grhocrit units, set by shooting (read-only output)"),
        ("omega_solved", c_double, "Omega_phi(a=1) achieved by the shooting (diagnostic)"),
        ("npoints", c_int, "number of points for background integration spacing"),
        (
            "fde",
            AllocatableArrayDouble,
            "after initialized, the calculated background dark energy fractions at sampled_a",
        ),
        ("__ddfde", AllocatableArrayDouble),
    )
    _fortran_class_name_ = "TTrackerQuintessence"

    def set_params(self, pot_type=1, ic_mode=0, alpha=1.0, lam=1.0, phi_ini=0.01):
        """
        Set tracker quintessence parameters.

        :param pot_type: 1=Ratra-Peebles (V ~ phi^-alpha), 2=exponential (V ~ exp(-lam*phi))
        :param ic_mode: 0=frozen phi_ini; 1=matter-era Ratra-Peebles tracker (physical phiCDM)
        :param alpha: Ratra-Peebles exponent
        :param lam: exponential slope
        :param phi_ini: initial (frozen) field value at a=1e-6 in reduced Planck units
        """
        if ic_mode not in (0, 1):
            raise ValueError("ic_mode must be 0 (frozen) or 1 (Ratra-Peebles tracker)")
        if ic_mode == 1 and pot_type != 1:
            raise ValueError("ic_mode=1 is defined only for the Ratra-Peebles potential")
        if alpha <= 0 or lam <= 0 or (ic_mode == 0 and phi_ini <= 0):
            raise ValueError("active tracker-quintessence parameters must be positive")
        self.pot_type = pot_type
        self.ic_mode = ic_mode
        self.alpha = alpha
        self.lam = lam
        self.phi_ini = phi_ini


@fortran_class
class CoupledQuintessence(TrackerQuintessence):
    r"""
    Conformal (Amendola-type) coupled quintessence: the CDM mass runs with the field,

    .. math:: m_{\rm dm}(\phi) = m_0\, e^{-\beta\phi},

    with :math:`\phi` in reduced Planck units (:math:`\kappa=1`), so the CDM density dilutes as
    :math:`\rho_c = \rho_{c0}\,a^{-3} e^{-\beta(\phi-\phi_0)}` (:math:`\phi_0=\phi(a{=}1)`) and the
    Klein-Gordon equation gains a source :math:`+\beta\rho_c`. The tracker potential
    (``pot_type`` 1=Ratra-Peebles, 2=exponential) and ``ln A`` shooting are inherited from
    :class:`TrackerQuintessence`; the coupled background is solved by a fixed-point outer loop and
    the excess CDM density flows through the polymorphic ``CDM_BackgroundCorrection`` hook so
    :math:`H(a)` is exact. Perturbations add the CDM fifth force (:math:`G_{\rm eff}=G(1+2\beta^2)`)
    and the field density source in synchronous gauge.

    The ``beta`` sign convention matches the N-body solver ``scalar_de_commons.f90``
    (``use_coupled_de``, ``beta_cde``): :math:`\rho_c\propto a^{-3}e^{-\beta\phi}` with KG source
    :math:`+\beta\rho_c`. ``beta=0`` reproduces :class:`TrackerQuintessence` exactly.

    Usage::

        pars.DarkEnergy = CoupledQuintessence()
        pars.DarkEnergy.set_params(pot_type=1, alpha=1.0, beta=0.05)

    References: Amendola PRD 62 043511 (2000); Amendola astro-ph/0311175.
    """

    _fields_ = (("beta", c_double, "conformal coupling; m_dm ~ exp(-beta*phi) (lagRamses beta_cde sign)"),)
    _fortran_class_name_ = "TCoupledQuintessence"

    def set_params(self, pot_type=1, alpha=1.0, lam=1.0, phi_ini=0.01, beta=0.0):
        """
        Set coupled tracker quintessence parameters.

        :param pot_type: 1=Ratra-Peebles (V ~ phi^-alpha), 2=exponential (V ~ exp(-lam*phi))
        :param alpha: Ratra-Peebles exponent
        :param lam: exponential slope
        :param phi_ini: initial (frozen) field value at a=1e-6 in reduced Planck units
        :param beta: conformal coupling of the running CDM mass m_dm ~ exp(-beta*phi)
        """
        super().set_params(pot_type=pot_type, alpha=alpha, lam=lam, phi_ini=phi_ini)
        self.beta = beta


@fortran_class
class InteractingDE(DarkEnergyEqnOfState):
    """
    Interacting Dark Energy model with DM-DE energy-momentum exchange.

    Dark energy and dark matter exchange energy: Q_mu.
    Background: rho_c' + 3H*rho_c = Q, rho_de' + 3(1+w)H*rho_de = -Q

    Interaction types:
        1: Q = xi * H * rho_de
        2: Q = xi * H * rho_c
        3: Q = xi * H * (rho_c + rho_de)

    References: Valiviita+ 2008, Costa+ 2017, IDECAMB (arXiv:2306.01593).
    """

    # Cannot add _fields_ here because DarkEnergyEqnOfState has unmapped
    # TCubicSpline fields in Fortran (same issue as DarkEnergyPPF).
    # Parameters are pushed to Fortran via SetIDEParams method instead.
    # IMPORTANT: Use set_params() to set IDE parameters (not individual attributes),
    # because pars.DarkEnergy returns a new Python wrapper each time.

    _fortran_class_module_ = "InteractingDE"
    _fortran_class_name_ = "TInteractingDE"

    _methods_ = (("SetIDEParams", [POINTER(c_double), POINTER(c_int), POINTER(c_double)]),)

    def set_params(self, w=-1.0, wa=0, xi_ide=0.0, interaction_type=1, cs2_ide=1.0, **kwargs):
        """
        Set interacting DE parameters. Must be called as a single method
        (not via individual attribute assignment).

        :param w: w(0) equation of state
        :param wa: -dw/da(0)
        :param xi_ide: DM-DE coupling strength
        :param interaction_type: Q kernel type (1=xi*H*rho_de, 2=xi*H*rho_c, 3=xi*H*(rho_c+rho_de))
        :param cs2_ide: DE rest-frame sound speed squared
        """
        self.w = w
        self.wa = wa
        if interaction_type not in (1, 2, 3):
            raise CAMBError("interaction_type must be 1, 2, or 3")
        self.f_SetIDEParams(
            byref(c_double(float(xi_ide))),
            byref(c_int(int(interaction_type))),
            byref(c_double(float(cs2_ide))),
        )


@fortran_class
class RunningVacuum(DarkEnergyModel):
    r"""
    Running Vacuum Model (RVM, Sola group) with :math:`\Lambda(H^2) = c_0 + \nu H^2`,
    i.e. :math:`\rho_\Lambda(H) = (3/8\pi G)(c_0 + \nu H^2)`. The vacuum exchanges
    energy with cold dark matter only (baryons and radiation stay uncoupled), giving
    the exact background

    .. math::
        \rho_c(a) &= \rho_{c0}\,a^{-3(1-\nu)}, \\
        \rho_\Lambda(a) &= \rho_{\Lambda 0} + \frac{\nu}{1-\nu}\,\rho_{c0}\,
                           \left(a^{-3(1-\nu)} - 1\right).

    The vacuum has :math:`w=-1` identically and is smooth
    (:math:`\delta\rho_\Lambda=0` in the CDM frame; the standard
    Gomez-Valent/Sola/Basilakos linear treatment). The energy transfer feeds
    unclustered particles into the CDM, adding a dilution term
    :math:`-3\nu\mathcal{H}\,\delta_c` to the perturbed CDM continuity. For
    :math:`\nu=0` the model is bit-identical to :math:`\Lambda`CDM.

    Usage::

        pars.DarkEnergy = RunningVacuum()
        pars.DarkEnergy.set_params(nu=0.001)

    References: Sola 2013 (arXiv:1306.1527), Gomez-Valent & Sola 2015
    (arXiv:1409.7048).
    """

    _fields_ = (
        ("nu", c_double, "running coefficient nu in Lambda(H^2)=c0+nu*H^2 (0 = LCDM)"),
        ("__grhoc_rvm", c_double, "cached 8 pi G rho_c0 (internal)"),
    )

    _fortran_class_module_ = "DarkEnergyRunningVacuum"
    _fortran_class_name_ = "TRunningVacuum"

    def set_params(self, nu=0.0):
        """
        Set the running-vacuum parameter.

        :param nu: dimensionless running coefficient of Lambda(H^2)=c0+nu*H^2
                   (typical |nu| <= few x 1e-3; nu=0 gives LCDM).
        """
        self.nu = nu


@fortran_class
class KEssence(DarkEnergyEqnOfState):
    r"""
    Purely kinetic k-essence dark energy with Lagrangian

    .. math::
        P(X) = M^4\,(-\tilde X + \tilde X^2), \qquad \tilde X = X/M^4,

    (the overall amplitude :math:`M^4` only fixes the DE budget). The background is
    algebraic (no shooting): the shift-symmetric equation of motion integrates to
    :math:`(2\tilde X - 1)\sqrt{\tilde X} = C_0\,a^{-3}` with
    :math:`C_0 = (2x_0-1)\sqrt{x_0}`, where :math:`x_0 = \tilde X(a{=}1)` is the single
    shape parameter (must be > 1/2). This gives

    .. math::
        w(a) = \frac{\tilde X - 1}{3\tilde X - 1}, \qquad
        c_s^2(a) = \frac{2\tilde X - 1}{6\tilde X - 1},

    with the early limit :math:`w\to 1/3,\ c_s^2\to 1/3` (dark-radiation-like) and the
    late limit :math:`w\to -1^+,\ c_s^2\to 0^+`. The time-varying rest-frame sound
    speed makes the DE cluster below the horizon down to its (small) Jeans scale,
    unlike a :math:`c_s^2=1` fluid.

    Matches the N-body solver convention ``kes_x0`` in scalar_de_commons.f90.

    Usage::

        pars.DarkEnergy = KEssence()
        pars.DarkEnergy.set_params(x0=0.6)
    """

    # Cannot declare extra _fields_ here: DarkEnergyEqnOfState ends with unmapped
    # TCubicSpline fields in Fortran (same constraint as DarkEnergyPPF/InteractingDE).
    # x0 is pushed to Fortran via the SetKEssenceParams method instead.
    _fortran_class_module_ = "DarkEnergyKEssence"
    _fortran_class_name_ = "TKEssence"

    _methods_ = (("SetKEssenceParams", [POINTER(c_double)]),)

    def set_params(self, x0=0.6):
        """
        Set purely-kinetic k-essence parameters.

        :param x0: dimensionless kinetic term Xt = X/M^4 today; must be > 1/2
                   (x0 -> 1/2+ approaches LCDM w=-1; larger x0 => w further from -1).
        """
        if x0 <= 0.5:
            raise CAMBError("KEssence x0 (= Xt today) must be > 1/2 for rho>0, cs2>0")
        self.f_SetKEssenceParams(byref(c_double(float(x0))))


@fortran_class
class Chaplygin(DarkEnergyEqnOfState):
    r"""
    Generalized Chaplygin gas (GCG) dark energy with equation of state

    .. math::
        p = -A\,\rho^{-\alpha}, \qquad \alpha \ge 0.

    Writing :math:`A_s = A/\rho_{de,0}^{1+\alpha}\in(0,1)` the background is fully
    closed-form (no shooting):

    .. math::
        \rho_{de}(a) &= \rho_{de,0}\,[A_s + (1-A_s)\,a^{-3(1+\alpha)}]^{1/(1+\alpha)},\\
        w(a) &= \frac{-A_s}{A_s + (1-A_s)\,a^{-3(1+\alpha)}},\\
        c_s^2(a) &= \frac{dp}{d\rho} = -\alpha\,w(a).

    The adiabatic sound speed :math:`c_s^2=-\alpha w` is used as the rest-frame
    sound speed of the fluid perturbations. Limits: :math:`A_s\to 1` gives
    :math:`w=-1` (pure :math:`\Lambda`) for any :math:`\alpha`; the gas is
    dust-like (:math:`w\to 0^-`) at early times and :math:`w=-A_s` today.

    As :math:`\alpha` grows the late-time sound speed :math:`\sim\alpha A_s`
    drives the small-scale oscillations/blow-up of unified GCG
    (`Sandvik et al. 2004 <https://arxiv.org/abs/astro-ph/0212114>`_), so
    unified-DM/DE usage requires :math:`\alpha\ll 1`.

    Matches the tabulated background convention lagRamses reuses in
    scalar_de_commons.f90.

    Usage::

        pars.DarkEnergy = Chaplygin()
        pars.DarkEnergy.set_params(As=0.75, alpha=0.01)
    """

    # Cannot declare extra _fields_ here: DarkEnergyEqnOfState ends with unmapped
    # TCubicSpline fields in Fortran (same constraint as DarkEnergyPPF/KEssence).
    # As, alpha are pushed to Fortran via the SetChaplyginParams method instead.
    _fortran_class_module_ = "DarkEnergyChaplygin"
    _fortran_class_name_ = "TChaplygin"

    _methods_ = (("SetChaplyginParams", [POINTER(c_double), POINTER(c_double)]),)

    def set_params(self, As=0.75, alpha=0.01):
        """
        Set generalized Chaplygin gas parameters.

        :param As: A/rho_de0^(1+alpha), must be in (0,1) (As -> 1 approaches LCDM)
        :param alpha: GCG exponent, must be >= 0 (alpha -> 0 approaches LCDM;
                      larger alpha raises the late-time sound speed cs2 = -alpha*w)
        """
        if not (0.0 < As < 1.0):
            raise CAMBError("Chaplygin As (= A/rho_de0^(1+alpha)) must be in (0,1)")
        if alpha < 0.0:
            raise CAMBError("Chaplygin alpha must be >= 0")
        self.f_SetChaplyginParams(byref(c_double(float(As))), byref(c_double(float(alpha))))


@fortran_class
class FuzzyDMField(DarkEnergyModel):
    """
    Fuzzy/Ultralight Axion Dark Matter via Klein-Gordon background + EFA perturbations.

    Solves the KG equation for V(phi) = (1/2)m^2 phi^2 to get correct background
    evolution (frozen w=-1 at early times, matter-like w=0 after oscillation onset).
    After a_match (where m/H = N_match), switches to effective fluid approximation
    with Passaglia-Hu sound speed cs2 = k^2/(4m^2 a^2 + k^2).

    Placed in DarkEnergy slot. Usage::

        pars.set_cosmology(omch2=0.12*(1-f_axion), ...)
        pars.DarkEnergy = FuzzyDMField()
        pars.DarkEnergy.set_params(m_axion=1e-22, f_axion=0.05)

    References: Hu+ 2000, Hlozek+ 2015, Passaglia & Hu 2022, Marsh 2016.
    """

    _fields_ = (
        ("m_axion", c_double, "Axion mass [eV]"),
        ("f_axion", c_double, "Fraction of CDM+axion budget that is axion"),
        ("omega_axion_h2", c_double, "Axion density Omega_a h^2 (0 = use f_axion)"),
        ("N_match", c_int, "KG->EFA transition criterion m/H"),
        ("npoints_bg", c_int, "Background integration points"),
        ("min_steps_per_osc_bg", c_int, "Min steps per oscillation"),
        ("potential_type", c_int, "1=quadratic, 2=cosine [1-cos(phi/f)]^n"),
        ("n_potential", c_int, "Power of cosine potential"),
        ("f_decay", c_double, "Decay constant [M_Pl] for cosine potential"),
        ("use_improved_efa", c_bool, "Use Passaglia-Hu improved EFA"),
    )

    _fortran_class_module_ = "FuzzyDMField"
    _fortran_class_name_ = "TFuzzyDMField"

    def set_params(
        self,
        m_axion=1e-22,
        f_axion=0.0,
        omega_axion_h2=0.0,
        N_match=100,
        potential_type=1,
        n_potential=1,
        f_decay=1.0,
        use_improved_efa=True,
    ):
        """
        Set fuzzy DM field parameters.

        :param m_axion: axion mass [eV]
        :param f_axion: fraction of CDM+axion budget that is axion
        :param omega_axion_h2: axion density (0 uses f_axion * omch2)
        :param N_match: KG->EFA transition criterion m/H (default 100)
        :param potential_type: 1=quadratic, 2=cosine [1-cos(phi/f)]^n
        :param n_potential: power of cosine potential
        :param f_decay: decay constant [M_Pl] for cosine potential
        :param use_improved_efa: use Passaglia-Hu improved EFA corrections
        """
        self.m_axion = m_axion
        self.f_axion = f_axion
        self.omega_axion_h2 = omega_axion_h2
        self.N_match = N_match
        self.potential_type = potential_type
        self.n_potential = n_potential
        self.f_decay = f_decay
        self.use_improved_efa = use_improved_efa

    def __getstate__(self):
        raise TypeError("Cannot save class with splines")


@fortran_class
class HorndeskiDE(DarkEnergyEqnOfState):
    """
    Horndeski scalar-tensor gravity with alpha parameterization (QSA).

    Implements the Bellini & Sawicki (2014) alpha parameterization:
    alpha_X(a) = alpha_X_0 * Omega_DE(a) (proportional parameterization).

    Modifications:
        - Modified Poisson equation via mu(a) from alpha_K, alpha_B, alpha_T
        - Modified lensing (Weyl) potential via Sigma(a)
        - Modified tensor propagation: (2+alpha_M)*H friction, (1+alpha_T)*k^2 speed
        - Running Planck mass M*^2(a) from integrating alpha_M

    Mutually exclusive with MuSigmaMG.

    References: Bellini & Sawicki 2014 (arXiv:1404.3713),
    Pogosian & Silvestri 2016 (arXiv:1606.05339),
    Zumalacárregui+ 2017 (hi_class, arXiv:1605.06102).
    """

    _fortran_class_module_ = "HorndeskiDE"
    _fortran_class_name_ = "THorndeskiDE"

    _methods_ = (
        (
            "SetHorndeskiParams",
            [
                POINTER(c_double),
                POINTER(c_double),
                POINTER(c_double),
                POINTER(c_double),
                POINTER(c_double),
            ],
        ),
    )

    def set_params(self, w=-1.0, wa=0, alpha_K=0.0, alpha_B=0.0, alpha_M=0.0, alpha_T=0.0, M_star_ini=1.0, **kwargs):
        """
        Set Horndeski gravity parameters.

        :param w: w(0) DE equation of state
        :param wa: -dw/da(0)
        :param alpha_K: kineticity amplitude (proportional to Omega_DE)
        :param alpha_B: braiding amplitude
        :param alpha_M: Planck mass run rate amplitude
        :param alpha_T: tensor speed excess amplitude
        :param M_star_ini: initial M*^2/M_Pl^2
        """
        self.w = w
        self.wa = wa
        self.f_SetHorndeskiParams(
            byref(c_double(float(alpha_K))),
            byref(c_double(float(alpha_B))),
            byref(c_double(float(alpha_M))),
            byref(c_double(float(alpha_T))),
            byref(c_double(float(M_star_ini))),
        )

    def __getstate__(self):
        raise TypeError("Cannot save class with splines")


# short names for models that support w/wa
F2003Class._class_names.update(
    {
        "fluid": DarkEnergyFluid,
        "ppf": DarkEnergyPPF,
        "kessence": KEssence,
        "chaplygin": Chaplygin,
        "interacting_de": InteractingDE,
        "running_vacuum": RunningVacuum,
        "fuzzy_dm_field": FuzzyDMField,
        "horndeski": HorndeskiDE,
    }
)
