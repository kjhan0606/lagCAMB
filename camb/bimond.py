"""Experimental homogeneous BIMOND background, independent of CAMB's Boltzmann API.

The fixed model is alpha=-1, beta=1, f=1, flat FLRW, pressureless ordinary
baryons plus radiation, and empty twin matter. The signed negative-Z extension
is a research specification, not a derived cosmological completion of QUMOND.
There are no perturbations, transfer functions, or CAMBparams integration here.

Coordinates are q=(ln a, ln b, ln(X/N)); velocities are dq/dtau, tau=H_ref*t,
with N=1. Density coefficients are 8*pi*G*rho_0/H_ref**2. No constraint
projection, density adjustment, or regularization of singular branches is used.
"""

from ctypes import CDLL, POINTER, byref, c_double, c_int
from dataclasses import dataclass

import numpy as np
from numpy.ctypeslib import ndpointer
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .baseconfig import camblib

_vector = ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
_matrix = ndpointer(dtype=np.float64, ndim=2, flags="F_CONTIGUOUS")
# BIND(C) symbols must bypass the loader's compiler-specific module mangling.
_function = CDLL.__getitem__(camblib, "camb_bimond_function")
_function.argtypes = [POINTER(c_double)] * 3 + [_vector, POINTER(c_int)]
_function.restype = None
_background = CDLL.__getitem__(camblib, "camb_bimond_background")
_background.argtypes = [
    _vector,
    _vector,
    _vector,
    POINTER(c_double),
    _vector,
    _vector,
    _matrix,
    _matrix,
    POINTER(c_int),
]
_background.restype = None


class BIMONDDomainError(ValueError):
    """The specified state is outside the supported regular background domain."""


@dataclass(frozen=True)
class BIMONDParameters:
    """Fixed model coefficients; a0bar=a0_phys/(c*H_ref), M(0)=m0."""

    a0bar: float
    baryons: float
    radiation: float
    m0: float = 0.0

    def __post_init__(self):
        if not np.all(np.isfinite(self.as_array())):
            raise ValueError("BIMOND parameters must be finite")
        if self.a0bar <= 0 or self.baryons < 0 or self.radiation < 0:
            raise ValueError("Require a0bar > 0 and nonnegative baryon/radiation coefficients")

    def as_array(self):
        return np.array([self.a0bar, self.baryons, self.radiation, self.m0], dtype=np.float64)


@dataclass(frozen=True)
class BackgroundEvaluation:
    lagrangian: float
    constraint: float
    constraint_scale: float
    z: float
    gradient: np.ndarray
    mixed: np.ndarray
    hessian: np.ndarray
    hessian_condition: float

    @property
    def normalized_constraint(self):
        return self.constraint / max(self.constraint_scale, np.finfo(float).tiny)

    def acceleration(self, velocity):
        return np.linalg.solve(self.hessian, self.gradient - self.mixed @ velocity)


@dataclass(frozen=True)
class BackgroundSolution:
    """Sampled q/velocities and diagnostics, plus accepted-step constraint extrema.

    States have shape (number of samples, 6). The homogeneous velocity Hessian
    condition is a numerical regularity diagnostic, not a perturbative ghost test.
    """

    tau: np.ndarray
    states: np.ndarray
    constraint: np.ndarray
    normalized_constraint: np.ndarray
    z: np.ndarray
    hessian_condition: np.ndarray
    max_abs_normalized_constraint: float
    max_hessian_condition: float
    nfev: int


@dataclass(frozen=True)
class TensorPrincipalMatrices:
    """Tensor kinetic/gradient matrices for the selected action only."""

    kinetic: np.ndarray
    gradient: np.ndarray
    m: float
    s: float
    t: float
    kinetic_eigenvalues: np.ndarray


def _positive(value, name):
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def _state_vector(value, size, name):
    array = np.ascontiguousarray(value, dtype=np.float64)
    if array.shape != (size,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain {size} finite values")
    return array


def _check_status(status):
    if status.value == 2:
        raise BIMONDDomainError("The state reaches the excluded |Z| <= zmin surface")
    if status.value:
        raise BIMONDDomainError("Invalid or nonfinite state in the BIMOND Fortran evaluator")


def interaction_function(z, m0=0.0, *, zmin=1e-12):
    """Return (M, M', M'') for the signed primitive of simple-nu minus one.

    M'(Z)=nu(sqrt(abs(Z)))-1 > 0; M'' changes sign across Z=0. zmin is
    a rejection threshold, never an interpolation or smoothing parameter.
    """
    _positive(zmin, "zmin")
    values = np.empty(3, dtype=np.float64)
    status = c_int()
    _function(byref(c_double(z)), byref(c_double(m0)), byref(c_double(zmin)), values, byref(status))
    _check_status(status)
    return values


def tensor_principal_matrices(q, z, parameters, *, zmin=1e-12):
    """Return the two-metric tensor principal matrices, without evolving modes.

    For one polarization normalized by e_ij*e_ij=2, the quadratic action has
    prefactor a**3/(32*pi*G), kinetic matrix K, and spatial term -G/a**2.
    Here r=b/a, X=exp(w), s=r**2/X**2, t=r**3/X, and
    m=sqrt(X)*r**1.5*M'(Z). These formulae apply only to the fixed
    alpha=-1, beta=1, f=1 action used in this module.

    Since the specified interaction has M'>0, this action cannot have both
    positive definite K and a stable tensor gradient matrix: K>0 requires
    m<1, for which det(G)=-(t/s)*(1-m)-m<0. Other cases have a negative or
    degenerate kinetic direction. Real propagation speeds alone would not
    demonstrate absence of ghosts. This tensor test is independent of the
    homogeneous velocity Hessian; no regular-background test is imposed here.
    """
    q = _state_vector(q, 3, "q")
    if np.max(np.abs(q)) > 100:
        raise BIMONDDomainError("Coordinates exceed the supported exponential range")
    u, v, w = q
    mp = interaction_function(z, parameters.m0, zmin=zmin)[1]
    s = np.exp(2 * (v - u - w))
    t = np.exp(3 * (v - u) - w)
    m = np.exp(0.5 * w + 1.5 * (v - u)) * mp
    with np.errstate(over="ignore", invalid="ignore"):
        kinetic = np.array([[1 - m, m * s], [m * s, -t - m * (2 * s - 1)]])
        gradient = np.array([[1 - m, m], [m, -np.exp(v - u + w) - m]])
    if not np.all(np.isfinite(kinetic)) or not np.all(np.isfinite(gradient)):
        raise BIMONDDomainError("Nonfinite tensor principal matrix")
    return TensorPrincipalMatrices(kinetic, gradient, float(m), float(s), float(t), np.linalg.eigvalsh(kinetic))


def evaluate_background(q, velocity, parameters, *, zmin=1e-12, max_condition=1e12):
    """Evaluate the Fortran action and derivatives, rejecting degenerate Hessians."""
    _positive(zmin, "zmin")
    _positive(max_condition, "max_condition")
    q = _state_vector(q, 3, "q")
    velocity = _state_vector(velocity, 3, "velocity")
    scalars, gradient = np.empty(4), np.empty(3)
    mixed, hessian = np.empty((3, 3), order="F"), np.empty((3, 3), order="F")
    status = c_int()
    _background(
        q, velocity, parameters.as_array(), byref(c_double(zmin)), scalars, gradient, mixed, hessian, byref(status)
    )
    _check_status(status)
    condition = float(np.linalg.cond(hessian))
    if not np.isfinite(condition) or condition > max_condition:
        raise BIMONDDomainError(f"Singular or ill-conditioned homogeneous Hessian: condition={condition:.6g}")
    return BackgroundEvaluation(*map(float, scalars), gradient, mixed, hessian, condition)


def constraint_initial_velocity(
    q,
    ordinary_rate,
    twin_rate,
    parameters,
    *,
    branch,
    z_bracket,
    zmin=1e-12,
    max_condition=1e12,
    constraint_tolerance=1e-10,
):
    """Solve E=0 for W at fixed q,U,V without changing model coefficients.

    The caller supplies a sign (+1 or -1) and a finite, ordered Z bracket on
    that sign branch. A missing root or a singular (s-1)*V is an error. The
    monotonic equation is 2*Z*M'(Z)-M(Z)=-(T+rho)/C.
    """
    q = _state_vector(q, 3, "q")
    rates = _state_vector([ordinary_rate, twin_rate], 2, "rates")
    _positive(zmin, "zmin")
    _positive(constraint_tolerance, "constraint_tolerance")
    if branch not in (-1, 1):
        raise ValueError("branch must be +1 or -1")
    bracket = _state_vector(z_bracket, 2, "z_bracket")
    if bracket[0] >= bracket[1] or np.any(branch * bracket <= zmin):
        raise ValueError("z_bracket must be ordered and strictly within the requested nonzero branch")
    if np.max(np.abs(q)) > 100:
        raise BIMONDDomainError("Coordinates exceed the supported exponential range")
    u, v, w = q
    ordinary_rate, twin_rate = rates
    s_minus_one = np.expm1(2 * (v - u - w))
    denominator = s_minus_one * twin_rate
    if abs(denominator) <= 32 * np.finfo(float).eps * max(1, abs(twin_rate), abs(s_minus_one)):
        raise BIMONDDomainError("Cannot solve the constraint where (s-1)*V is singular")
    s = s_minus_one + 1
    kinetic = -3 * np.exp(3 * u) * ordinary_rate**2 + 3 * np.exp(3 * v - w) * twin_rate**2
    amplitude = parameters.a0bar**2 * np.exp(1.5 * (u + v) + 0.5 * w)
    density = parameters.baryons + parameters.radiation * np.exp(-u)
    target = -(kinetic + density) / amplitude

    def residual(z):
        m, mp, _ = interaction_function(z, parameters.m0, zmin=zmin)
        return 2 * z * mp - m - target

    try:
        z = brentq(residual, *bracket, xtol=np.nextafter(0.0, 1.0), rtol=8 * np.finfo(float).eps)
    except ValueError as error:
        raise BIMONDDomainError("The specified Z bracket does not contain a regular constraint root") from error
    k = 1.5 / parameters.a0bar**2
    lapse_rate = (
        z / k - 2 * ordinary_rate**2 + (3 + s) * ordinary_rate * twin_rate - (1 + s) * twin_rate**2
    ) / denominator
    velocity = np.array([ordinary_rate, twin_rate, lapse_rate])
    evaluation = evaluate_background(q, velocity, parameters, zmin=zmin, max_condition=max_condition)
    if branch * evaluation.z <= zmin or abs(evaluation.normalized_constraint) > constraint_tolerance:
        raise BIMONDDomainError("Constraint root lost accuracy when reconstructed as an initial velocity")
    return velocity


def integrate_background(
    parameters,
    initial_state,
    tau_span,
    *,
    sample_times,
    rtol,
    atol,
    zmin=1e-12,
    max_condition=1e12,
    constraint_tolerance=1e-8,
    max_step=np.inf,
):
    """Integrate unprojected Euler-Lagrange equations with DOP853.

    Explicit tolerances and sample times are required. Reject initial or
    accepted-step constraint violations, branch crossing (including internal
    RHS stages), |Z|<=zmin, and ill-conditioned Hessians. Constraint checking
    does not change the trajectory. This is a local background solver, not
    evidence that an observationally viable cosmological history exists.
    """
    initial = _state_vector(initial_state, 6, "initial_state")
    span = _state_vector(tau_span, 2, "tau_span")
    if span[1] <= span[0]:
        raise ValueError("tau_span must be increasing")
    samples = np.asarray(sample_times, dtype=np.float64)
    if (
        samples.ndim != 1
        or len(samples) == 0
        or not np.all(np.isfinite(samples))
        or np.any(np.diff(samples) <= 0)
        or samples[0] < span[0]
        or samples[-1] > span[1]
    ):
        raise ValueError("sample_times must be finite, strictly increasing, and within tau_span")
    _positive(rtol, "rtol")
    _positive(constraint_tolerance, "constraint_tolerance")
    absolute_tolerance = np.asarray(atol, dtype=np.float64)
    if (
        absolute_tolerance.shape not in ((), (6,))
        or np.any(~np.isfinite(absolute_tolerance))
        or np.any(absolute_tolerance <= 0)
    ):
        raise ValueError("atol must be positive and finite, either scalar or length six")
    if np.isnan(max_step) or max_step <= 0:
        raise ValueError("max_step must be positive")
    initial_evaluation = evaluate_background(
        initial[:3], initial[3:], parameters, zmin=zmin, max_condition=max_condition
    )
    branch = np.sign(initial_evaluation.z)
    if abs(initial_evaluation.normalized_constraint) > constraint_tolerance:
        raise BIMONDDomainError("Initial state violates the lapse constraint E=0")
    max_seen_condition = initial_evaluation.hessian_condition

    def evaluate(state):
        nonlocal max_seen_condition
        result = evaluate_background(state[:3], state[3:], parameters, zmin=zmin, max_condition=max_condition)
        if branch * result.z <= zmin:
            raise BIMONDDomainError("Integration attempted to cross the Z=0 branch surface")
        max_seen_condition = max(max_seen_condition, result.hessian_condition)
        return result

    def rhs(tau, state):
        return np.concatenate((state[3:], evaluate(state).acceleration(state[3:])))

    solution = solve_ivp(
        rhs,
        span,
        initial,
        method="DOP853",
        rtol=rtol,
        atol=absolute_tolerance,
        max_step=max_step,
        dense_output=True,
    )
    if not solution.success:
        raise BIMONDDomainError(f"Background integration failed: {solution.message}")
    accepted = [evaluate(state) for state in solution.y.T]
    states = solution.sol(samples).T
    sampled = [evaluate(state) for state in states]
    maximum_constraint = max(abs(item.normalized_constraint) for item in accepted + sampled)
    if maximum_constraint > constraint_tolerance:
        raise BIMONDDomainError(f"Unprojected constraint drift exceeded tolerance: {maximum_constraint:.6g}")
    return BackgroundSolution(
        samples.copy(),
        states,
        np.array([item.constraint for item in sampled]),
        np.array([item.normalized_constraint for item in sampled]),
        np.array([item.z for item in sampled]),
        np.array([item.hessian_condition for item in sampled]),
        maximum_constraint,
        max_seen_condition,
        solution.nfev,
    )
