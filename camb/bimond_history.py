"""Arbitrary-precision background and scalar fundamental-matrix diagnostics.

This is a finite-interval evolution of the positive-Einstein--Hilbert research
candidate, not a cosmological transfer function. Identity initial data define
basis responses, not physical primordial modes. No ghost directions are
removed. The complete local operator, matter inventory, gauge, and admission
gates are inherited from ``bimond_scalar``.

Background and scalar evolution are one-way coupled at linear order. A dense
high-precision background supplies every scalar collocation point. Background
constraints are monitored without projection; matter continuity is integrated
analytically. All expensive work occurs on explicit construction/call.
"""

from dataclasses import dataclass
from typing import Any

import mpmath as mp
import sympy as sp


class HistoryDomainError(ValueError):
    """The requested diagnostic history left its admitted numerical domain."""


def exact_number(value, name):
    """Accept exact symbolic numbers or declared decimals, never binary64 seeds."""
    if isinstance(value, (bool, float, sp.Float)):
        raise TypeError(f"{name} must be an exact number or decimal string, not a floating-point seed")
    if isinstance(value, sp.Expr):
        if value.free_symbols or value.has(sp.Float):
            raise ValueError(f"{name} must be an exact closed symbolic expression")
        result = mp.mpf(str(value.evalf(mp.mp.dps + 5)))
    elif isinstance(value, (str, int)):
        result = mp.mpf(value)
    else:
        raise TypeError(f"{name} must be an integer, decimal string, or exact SymPy expression")
    if not mp.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(frozen=True)
class SignedRationalInteraction:
    """The audited research function M'=2/(4+sqrt(|Z|)), with explicit signed continuation."""

    m0: Any = "0"

    def jet(self, z):
        if z == 0:
            raise HistoryDomainError("Interaction derivatives are singular at Z=0")
        y, sign = mp.sqrt(abs(z)), mp.sign(z)
        denominator = 4 + y
        return (
            exact_number(self.m0, "M0") + sign * (4 * y - 16 * mp.log1p(y / 4)),
            2 / denominator,
            -sign / (y * denominator**2),
            1 / (2 * y**3 * denominator**2) + 1 / (y**2 * denominator**3),
            -sign
            * (3 / (4 * y**5 * denominator**2) + 3 / (2 * y**4 * denominator**3) + 3 / (2 * y**3 * denominator**4)),
        )


class BackgroundHistory:
    """Dense unprojected six-variable background in tau=H_ref*t.

    Supply exact/decimal a,b,X,U,V,W,dust_density,radiation_rate,a0bar and
    radiation_kappa. alpha=beta=1 are explicit defaults for this diagnostic;
    other coefficients are rejected here rather than silently change models.
    Dust density*a**3 and radiation_rate*a are conserved exactly by their
    continuity integrals. This does not solve or project the lapse constraint.
    """

    def __init__(
        self,
        system,
        initial_values,
        *,
        k,
        end_time,
        dps=80,
        background_tolerance="1e-65",
        on_shell_tolerance="1e-45",
        max_condition="1e30",
        zmin="1e-20",
        degree=30,
        interaction=None,
    ):
        if dps < 50:
            raise ValueError("At least 50 working decimal digits are required")
        if degree < 8:
            raise ValueError("Background Taylor degree must be at least eight")
        self.system = system
        self.dps = dps
        self.interaction = interaction or SignedRationalInteraction()
        self.initial_exact = {str(key): value for key, value in initial_values.items()}
        self.background_evaluations = 0
        self.operator_evaluations = 0
        self.maximum_constraint = mp.mpf(0)
        self.maximum_tadpole = mp.mpf(0)
        self.maximum_conditions = {}
        self._operator_cache = {}
        self._state_cache = {}
        self.degree = degree
        with mp.workdps(dps):
            self.end_time = exact_number(end_time, "end_time")
            self.k = exact_number(k, "k")
            self.background_tolerance = exact_number(background_tolerance, "background_tolerance")
            self.on_shell_tolerance = exact_number(on_shell_tolerance, "on_shell_tolerance")
            self.max_condition = exact_number(max_condition, "max_condition")
            self.zmin = exact_number(zmin, "zmin")
            if self.end_time <= 0 or self.k <= 0 or self.zmin <= 0:
                raise ValueError("end_time, k, and zmin must be positive")
            if not 0 < self.background_tolerance < self.on_shell_tolerance * mp.mpf("1e-10"):
                raise ValueError("Background tolerance must be at least ten orders below on-shell tolerance")
            if self.background_tolerance < mp.mpf(10) ** (-dps + 10):
                raise ValueError("Working precision needs at least ten guard digits above the background tolerance")
            required = ("a", "b", "X", "U", "V", "W", "dust_density", "radiation_rate", "a0bar", "radiation_kappa")
            initial = {name: exact_number(self.initial_exact[name], name) for name in required}
            for name in ("alpha", "beta"):
                if exact_number(self.initial_exact.get(name, 1), name) != 1:
                    raise ValueError("This history implements only the audited alpha=beta=1 candidate")
                initial[name] = mp.mpf(1)
            for name in ("a", "b", "X", "dust_density", "radiation_rate", "a0bar", "radiation_kappa"):
                if initial[name] <= 0:
                    raise ValueError(f"Positive {name} is required")
            self.constants = {name: initial[name] for name in ("alpha", "beta", "a0bar", "radiation_kappa")}
            self.dust_constant = initial["dust_density"] * initial["a"] ** 3
            self.radiation_constant = initial["radiation_rate"] * initial["a"]
            self.initial_state = tuple(mp.log(initial[name]) for name in ("a", "b", "X")) + tuple(
                initial[name] for name in ("U", "V", "W")
            )
            self.branch = None
            _, z = self.values_from_state(self.initial_state)
            self.branch = int(mp.sign(z))
            # The initial state is checked before construction of the dense ODE.
            self.operator_from_state(self.initial_state)
            self._dense = mp.odefun(self._rhs, 0, self.initial_state, tol=self.background_tolerance, degree=degree)

    def values_from_state(self, state):
        """Evaluate model jets using the active precision, including temporary guard precision."""
        if len(state) != 6 or any(not mp.isfinite(value) for value in state):
            raise HistoryDomainError("Invalid six-component background state")
        a, b, X = (mp.exp(value) for value in state[:3])
        values = {self.system.symbols[name]: value for name, value in self.constants.items()}
        values.update(
            zip((self.system.symbols[name] for name in ("a", "b", "X", "U", "V", "W")), (a, b, X, *state[3:]))
        )
        values[self.system.symbols["dust_density"]] = self.dust_constant / a**3
        values[self.system.symbols["radiation_rate"]] = self.radiation_constant / a
        values[self.system.symbols["k"]] = self.k
        z = self.system._evaluate("history_Z", self.system.z, values)
        if abs(z) <= self.zmin or (self.branch is not None and mp.sign(z) != self.branch):
            raise HistoryDomainError("Background reached the excluded Z=0 surface or changed branch")
        values.update(
            zip(
                (self.system.symbols[name] for name in ("Mbar", "Mprime", "Msecond", "Mthird", "Mfourth")),
                self.interaction.jet(z),
            )
        )
        return values, z

    def _rhs(self, time, state):
        # mpmath's Taylor generator temporarily raises precision substantially.
        # Do not wrap this callback in workdps(self.dps), which would destroy
        # its high-order derivative calculation.
        self.background_evaluations += 1
        values, _ = self.values_from_state(state)
        H = self.system._evaluate("background_H", self.system.background_hessian, values)
        rhs = self.system._evaluate("background_rhs", self.system.background_rhs, values)
        acceleration = mp.lu_solve(H, rhs)
        return list(state[3:]) + list(acceleration)

    def state(self, time):
        with mp.workdps(self.dps):
            time = mp.mpf(time)
            if not 0 <= time <= self.end_time:
                raise HistoryDomainError("Time outside the explicitly requested finite interval")
            if time not in self._state_cache:
                self._state_cache[time] = tuple(self._dense(time))
            return self._state_cache[time]

    def operator_from_state(self, state):
        """Apply all existing on-shell/conditioning gates at the active precision."""
        self.operator_evaluations += 1
        values, _ = self.values_from_state(state)
        operator = self.system.evaluate(
            values, dps=max(self.dps, mp.mp.dps), tolerance=self.on_shell_tolerance, max_condition=self.max_condition
        )
        self.maximum_constraint = max(self.maximum_constraint, operator.diagnostics["energy_residual"])
        self.maximum_tadpole = max(self.maximum_tadpole, operator.diagnostics["tadpole_residual"])
        for name, value in operator.diagnostics.items():
            if name.endswith("_condition"):
                self.maximum_conditions[name] = max(value, self.maximum_conditions.get(name, 0))
        return operator

    def operator(self, time):
        with mp.workdps(self.dps):
            time = mp.mpf(time)
            if time not in self._operator_cache:
                self._operator_cache[time] = self.operator_from_state(self.state(time))
            return self._operator_cache[time]

    def flow(self, state, operator=None):
        operator = operator or self.operator_from_state(state)
        return mp.matrix(
            [*state[3:], *(operator.background[self.system.symbols[f"{name}_acc"]] for name in ("U", "V", "W"))]
        )


def first_order_matrix(operator):
    """State order (x,xdot), x=(zeta,phihat,zetahat,ellhat,chi)."""
    result = mp.matrix(10, 10)
    damping = -operator.kinetic_inverse * operator.damping
    stiffness = -operator.kinetic_inverse * operator.stiffness
    for i in range(5):
        result[i, 5 + i] = 1
        for j in range(5):
            result[5 + i, j] = stiffness[i, j]
            result[5 + i, 5 + j] = damping[i, j]
    return result


def gauss_coefficients(stages):
    """Arbitrary-precision Gauss collocation coefficients from Lagrange polynomials."""
    if stages < 2 or stages > 10:
        raise ValueError("Use between two and ten Gauss stages")
    nodes, _ = mp.gauss_quadrature(stages, "legendre")
    nodes = [(item + 1) / 2 for item in nodes]
    polynomials = []
    for j in range(stages):
        polynomial = [mp.mpf(1)]
        denominator = mp.mpf(1)
        for other in range(stages):
            if other == j:
                continue
            updated = [mp.mpf(0)] * (len(polynomial) + 1)
            for power, coefficient in enumerate(polynomial):
                updated[power] -= nodes[other] * coefficient
                updated[power + 1] += coefficient
            polynomial = updated
            denominator *= nodes[j] - nodes[other]
        polynomials.append([coefficient / denominator for coefficient in polynomial])
    integrated = [
        [mp.mpf(0)] + [coefficient / (power + 1) for power, coefficient in enumerate(polynomial)]
        for polynomial in polynomials
    ]
    A = mp.matrix([[mp.polyval(polynomial[::-1], node) for polynomial in integrated] for node in nodes])
    weights = [sum(polynomial) for polynomial in integrated]
    return nodes, weights, A, polynomials, integrated


def solve_multiple(matrix, rhs):
    """One public LU decomposition followed by triangular solves for all columns."""
    P, L, U = mp.lu(matrix)
    result = P * rhs
    for i in range(matrix.rows):
        for column in range(rhs.cols):
            result[i, column] = (result[i, column] - sum(L[i, j] * result[j, column] for j in range(i))) / L[i, i]
    for i in range(matrix.rows - 1, -1, -1):
        for column in range(rhs.cols):
            result[i, column] = (
                result[i, column] - sum(U[i, j] * result[j, column] for j in range(i + 1, matrix.rows))
            ) / U[i, i]
    return result


@dataclass
class CollocationStep:
    start: Any
    width: Any
    initial: Any
    derivatives: list
    polynomials: list
    integrated: list

    def evaluate(self, time):
        theta = (time - self.start) / self.width
        value = self.initial.copy()
        derivative = mp.matrix(self.initial.rows, self.initial.cols)
        for stage, polynomial, integrated in zip(self.derivatives, self.polynomials, self.integrated):
            value += self.width * mp.polyval(integrated[::-1], theta) * stage
            derivative += mp.polyval(polynomial[::-1], theta) * stage
        return value, derivative


class ScalarHistory:
    """Gauss-Legendre history of the full ten-dimensional fundamental matrix.

    The initial identity is a mathematical basis-response normalization only.
    The last five coordinates are derivatives with respect to tau=H_ref*t;
    k is the constant dimensionless comoving wavenumber of cos(kz)/sin(kz).
    Endpoint order is twice the number of stages. The dense polynomial has
    lower order and its differential defect is measured separately.
    """

    def __init__(self, background, *, steps, stages=6):
        if not isinstance(steps, int) or steps <= 0:
            raise ValueError("steps must be a positive integer")
        self.background = background
        self.dps = background.dps
        self.steps = []
        self.stage_count = stages
        self.maximum_linear_solve_residual = mp.mpf(0)
        self.maximum_dense_defect_absolute = mp.mpf(0)
        self.maximum_dense_defect_relative = mp.mpf(0)
        self.endpoint_times = []
        self.endpoint_matrices = []
        with mp.workdps(self.dps):
            nodes, weights, A, polynomials, integrated = gauss_coefficients(stages)
            width = background.end_time / steps
            fundamental = mp.eye(10)
            self.endpoint_times.append(mp.mpf(0))
            self.endpoint_matrices.append(fundamental.copy())
            for index in range(steps):
                start = index * width
                operators = [first_order_matrix(background.operator(start + node * width)) for node in nodes]
                block = mp.eye(10 * stages)
                rhs = mp.matrix(10 * stages, 10)
                for i in range(stages):
                    for row in range(10):
                        for column in range(10):
                            rhs[10 * i + row, column] = fundamental[row, column]
                            for j in range(stages):
                                block[10 * i + row, 10 * j + column] -= width * A[i, j] * operators[j][row, column]
                stage_values = solve_multiple(block, rhs)
                residual = mp.norm(block * stage_values - rhs) / max(1, mp.norm(rhs))
                self.maximum_linear_solve_residual = max(self.maximum_linear_solve_residual, residual)
                if residual > background.on_shell_tolerance:
                    raise HistoryDomainError("Collocation stage linear solve failed its residual tolerance")
                derivatives = [operators[i] * stage_values[10 * i : 10 * (i + 1), :] for i in range(stages)]
                step = CollocationStep(start, width, fundamental.copy(), derivatives, polynomials, integrated)
                self.steps.append(step)
                for weight, derivative in zip(weights, derivatives):
                    fundamental += width * weight * derivative
                time = (index + 1) * width
                # Explicitly check the accepted endpoint background as well.
                background.operator(time)
                self.endpoint_times.append(time)
                self.endpoint_matrices.append(fundamental.copy())
                probe_time = start + width / 2
                value, derivative = step.evaluate(probe_time)
                expected = first_order_matrix(background.operator(probe_time)) * value
                absolute = mp.norm(derivative - expected)
                relative = absolute / max(1, mp.norm(derivative), mp.norm(expected))
                self.maximum_dense_defect_absolute = max(self.maximum_dense_defect_absolute, absolute)
                self.maximum_dense_defect_relative = max(self.maximum_dense_defect_relative, relative)

    def evaluate(self, time):
        with mp.workdps(self.dps):
            time = mp.mpf(time)
            if not 0 <= time <= self.background.end_time:
                raise HistoryDomainError("Time outside the scalar history")
            width = self.steps[0].width
            index = min(len(self.steps) - 1, int(mp.floor(time / width)))
            return self.steps[index].evaluate(time)


def shift_and_rate(operator, fundamental):
    """All ten basis responses of the algebraic shifts and their exact ODE rates."""
    x, xd = fundamental[:5, :], fundamental[5:, :]
    xdd = -operator.kinetic_inverse * (operator.damping * xd + operator.stiffness * x)
    C, E = operator.blocks["C"], operator.blocks["E"]
    shift = -operator.auxiliary_inverse * (C.T * xd + E.T * x)
    rate = -operator.auxiliary_inverse * (
        operator.derivatives["C"].T * xd
        + C.T * xdd
        + operator.derivatives["E"].T * x
        + E.T * xd
        + operator.derivatives["Ay"] * shift
    )
    return shift, rate


def shift_second_derivative(history, time, *, difference_step="1e-28", extra_dps=40):
    """Differentiate the exact shift-rate function along the coupled ODE tangent.

    A symmetric fourth-order directional difference is evaluated with guard
    precision at h and h/2 and Richardson extrapolated. Probes use the actual
    background/scalar tangent; their O(h**2) background-constraint departure
    must still pass the existing on-shell gate. They are derivative probes,
    not replacement background trajectories. Step and precision convergence
    must be checked by the caller/report; neither shift acceleration is zeroed.
    """
    background = history.background
    with mp.workdps(history.dps + extra_dps):
        step = exact_number(difference_step, "difference_step")
        if step <= 0:
            raise ValueError("difference_step must be positive")
        state = mp.matrix(background.state(time))
        fundamental, _ = history.evaluate(time)
        operator = background.operator_from_state(state)
        background_flow = background.flow(state, operator)
        scalar_flow = first_order_matrix(operator) * fundamental

        def value(offset):
            displaced = state + offset * background_flow
            perturbed = background.operator_from_state(displaced)
            return shift_and_rate(perturbed, fundamental + offset * scalar_flow)[1]

        def estimate(h):
            return (value(-2 * h) - 8 * value(-h) + 8 * value(h) - value(2 * h)) / (12 * h)

        coarse, fine = estimate(step), estimate(step / 2)
        result = (16 * fine - coarse) / 15
        absolute = mp.norm(fine - coarse)
        relative = absolute / max(1, mp.norm(result))
        return result, {
            "step_difference_absolute": absolute,
            "step_difference_relative": relative,
            "difference_step": step,
            "working_dps": history.dps + extra_dps,
        }


def reconstruction_residuals(history, time, shift_acceleration):
    """Original retained/lapse residuals for every basis column, with explicit shift acceleration."""
    with mp.workdps(history.dps):
        fundamental, _ = history.evaluate(time)
        operator = history.background.operator(time)
        system = operator.system
        maximum_retained = mp.mpf(0)
        maximum_lapse = mp.mpf(0)
        density = mp.matrix(1, 10)
        for column in range(10):
            reconstruction = operator.reconstruct(
                fundamental[:5, column], fundamental[5:, column], shift_acceleration=shift_acceleration[:, column]
            )
            terms = (
                operator.raw["K"] * reconstruction.raw_acceleration,
                (operator.raw["Kdot"] + operator.raw["B"] - operator.raw["B"].T) * reconstruction.raw_velocity,
                (operator.raw["Bdot"] - operator.raw["D"]) * reconstruction.raw_position,
            )
            retained = mp.norm(operator.retained_residual(reconstruction)) / max(1, *(mp.norm(term) for term in terms))
            local = dict(operator.background)
            local.update(zip(system.fields, reconstruction.raw_position))
            local.update(zip(system.velocities, reconstruction.raw_velocity))
            local.update(zip(system.accelerations, reconstruction.raw_acceleration))
            local[system.dust_symbol] = reconstruction.dust_density_perturbation
            lapse = system._evaluate("history_lapse", system.lapse_equation, local)
            lapse_terms = sp.Matrix(sp.Add.make_args(sp.expand(system.lapse_equation)))
            scale = max(1, sum(abs(item) for item in system._evaluate("history_lapse_terms", lapse_terms, local)))
            maximum_retained = max(maximum_retained, retained)
            maximum_lapse = max(maximum_lapse, abs(lapse) / scale)
            density[0, column] = reconstruction.dust_density_perturbation
        return {"retained_relative": maximum_retained, "lapse_relative": maximum_lapse}, density
