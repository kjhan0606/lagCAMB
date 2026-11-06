"""Exact local constrained scalar operators for the BIMOND perfect-fluid action.

This module reduces the complete quadratic action in dust-clock and ordinary
spatial-shear gauge, without dropping finite-wavenumber or background-derivative
terms. It exports symbolic coefficient blocks and evaluates their Schur
complement at an on-shell local background jet. Time derivatives follow the
background equations, including accelerations and jerks; they never freeze
U,V,W. There is no photon/neutrino hierarchy, transfer function, initial-mode
prescription, or cosmological-viability claim. Previously audited local
backgrounds have scalar ghosts; this operator does not remove them.

Call symbolic construction and numerical evaluation only on the designated
scientific-computation scheduler. Importing this module performs no derivation.
"""

from dataclasses import dataclass, field
from itertools import pairwise
from typing import Any

import mpmath as mp
import sympy as sp


class ScalarReductionError(ValueError):
    """The background or auxiliary system does not pass the reduction gates."""


@dataclass
class ScalarSystem:
    """Exact coefficient blocks before any potentially large symbolic inverse.

    Raw field order is (zeta,S,phihat,T,zetahat,ellhat,chi). Dynamical order is
    (zeta,phihat,zetahat,ellhat,chi); algebraic order is (S,T). Fourier scalars
    use cos(kz), shifts sin(kz), and the action is spatially averaged, giving
    the usual factor one half. k is constant comoving wavenumber.
    """

    symbols: dict
    fields: tuple
    velocities: tuple
    accelerations: tuple
    dynamic_indices: tuple
    auxiliary_indices: tuple
    lagrangian: sp.Expr
    raw_blocks: dict
    blocks: dict
    block_derivatives: dict
    background_hessian: sp.Matrix
    background_rhs: sp.Matrix
    background_hessian_dot: sp.Matrix
    background_rhs_dot: sp.Matrix
    energy: sp.Expr
    z: sp.Expr
    background_flow: dict
    auxiliary_gates: sp.Matrix
    tadpole_left: sp.Matrix
    tadpole_right: sp.Matrix
    dust_reconstruction: sp.Expr
    lapse_equation: sp.Expr
    dust_symbol: sp.Symbol
    _compiled: dict = field(default_factory=dict, repr=False)

    @property
    def dynamic_fields(self):
        return tuple(self.fields[index] for index in self.dynamic_indices)

    @property
    def auxiliary_fields(self):
        return tuple(self.fields[index] for index in self.auxiliary_indices)

    def dt(self, expression):
        """Exact local background chain derivative; M jets are needed through Mfourth."""
        if expression.has(self.symbols["Mfourth"]):
            raise ScalarReductionError("Another derivative would require M fifth, beyond this local operator")
        return sp.expand(sum(sp.diff(expression, variable) * speed for variable, speed in self.background_flow.items()))

    def _evaluate(self, key, expression, values):
        """Lambdify once per block, using arbitrary-precision real arithmetic."""
        if key not in self._compiled:
            variables = tuple(sorted(expression.free_symbols, key=str))
            self._compiled[key] = (variables, sp.lambdify(variables, expression, modules="mpmath", cse=True))
        variables, function = self._compiled[key]
        missing = [str(variable) for variable in variables if variable not in values]
        if missing:
            raise ScalarReductionError(f"Missing local-jet values: {', '.join(missing)}")
        return function(*(values[variable] for variable in variables))

    def coefficient_export(self):
        """JSON-compatible exact block/operator export, with no expanded matrix inverse.

        Expressions are SymPy strings, not untrusted executable input. The
        blocks plus the explicit Schur/derivative rules fully reconstruct the
        operator. Background accelerations and jerks must be solved as stated.
        """

        def rows(matrix):
            return [[sp.sstr(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)]

        return {
            "schema": "bimond-perfect-fluid-scalar-operator-v1",
            "scope": "local constrained scalar operator; not a transfer-function or stability verdict",
            "raw_fields": list(map(str, self.fields)),
            "dynamic_fields": list(map(str, self.dynamic_fields)),
            "auxiliary_fields": list(map(str, self.auxiliary_fields)),
            "symbols": {name: str(symbol) for name, symbol in self.symbols.items()},
            "raw_blocks": {name: rows(matrix) for name, matrix in self.raw_blocks.items()},
            "blocks": {name: rows(matrix) for name, matrix in self.blocks.items()},
            "block_time_derivatives": {name: rows(matrix) for name, matrix in self.block_derivatives.items()},
            "background_flow": {str(variable): sp.sstr(speed) for variable, speed in self.background_flow.items()},
            "background_hessian": rows(self.background_hessian),
            "background_rhs": rows(self.background_rhs),
            "background_hessian_dot": rows(self.background_hessian_dot),
            "background_rhs_dot": rows(self.background_rhs_dot),
            "background_energy_constraint": sp.sstr(self.energy),
            "Z": sp.sstr(self.z),
            "background_closure": ["H*acc=rhs", "H*jerk=rhs_dot-H_dot*acc"],
            "operator": {
                "Keff": "Kd-C*inv(Ay)*C.T",
                "Beff": "Bd-C*inv(Ay)*E.T",
                "Deff": "Dd-E*inv(Ay)*E.T",
                "inverse_derivative": "dt(inv(Ay))=-inv(Ay)*dt(Ay)*inv(Ay)",
                "equation": "Keff*xddot+(dt(Keff)+Beff-Beff.T)*xdot+(dt(Beff)-Deff)*x=0",
                "shift": "y=-inv(Ay)*(C.T*xdot+E.T*x)",
                "shift_dot": "-inv(Ay)*(dt(C).T*xdot+C.T*xddot+dt(E).T*x+E.T*xdot+dt(Ay)*y)",
            },
            "dust_density_reconstruction": sp.sstr(self.dust_reconstruction),
            "dust_lapse_equation": sp.sstr(self.lapse_equation),
            "auxiliary_gates": rows(self.auxiliary_gates),
            "tadpole_left": rows(self.tadpole_left),
            "tadpole_right": rows(self.tadpole_right),
            "assumptions": [
                "tau=0 fixes dust clock; ell=zeta fixes ordinary spatial shear",
                "dust multiplier imposes phi=0; phi equation reconstructs dlambda",
                "k is nonzero and constant; shifts have no velocity Hessian or antisymmetric symplectic pair",
                "background energy, tadpole, and invertibility gates are mandatory",
                "Mbar through Mfourth specify a consistent local interaction-function jet",
                "dust and P(Y) radiation only; no twin matter or Boltzmann hierarchy",
            ],
        }

    def evaluate(self, values, *, dps=80, tolerance="1e-45", max_condition="1e40"):
        """Evaluate the full operator on a supplied on-shell local background.

        Supply background symbols, k, and interaction derivatives Mbar through
        Mfourth as exact SymPy numbers or decimal strings. Accelerations and
        jerks are solved internally, never accepted as independently frozen
        inputs. No inverse is regularized. All matrices use mpmath arithmetic.
        """
        if dps < 50:
            raise ValueError("Use at least 50 digits for the cancellation-sensitive scalar reduction")
        with mp.workdps(dps):
            tol, condition_limit = mp.mpf(tolerance), mp.mpf(max_condition)
            if not 0 < tol < 1 or condition_limit <= 1:
                raise ValueError("Require 0<tolerance<1 and max_condition>1")
            supplied = {str(key): value for key, value in values.items()}
            local = {}
            for name, variable in self.symbols.items():
                if name in supplied:
                    raw = supplied[name]
                    local[variable] = mp.mpf(str(raw.evalf(dps))) if isinstance(raw, sp.Expr) else mp.mpf(raw)
                    if not mp.isfinite(local[variable]):
                        raise ScalarReductionError(f"Nonfinite background value: {name}")
            for name in ("U_acc", "V_acc", "W_acc", "U_jerk", "V_jerk", "W_jerk"):
                if name in supplied:
                    raise ScalarReductionError("Background accelerations and jerks must be solved, not supplied")
            for name in ("a", "b", "X", "a0bar", "dust_density", "radiation_rate", "radiation_kappa", "k"):
                if self.symbols[name] not in local or local[self.symbols[name]] <= 0:
                    raise ScalarReductionError(f"Require a positive supplied {name}")
            diagnostics = {}

            def inverse(matrix, label):
                singular = mp.svd(matrix, compute_uv=False)
                high, low = max(singular), min(singular)
                if low <= 0 or high / low > condition_limit:
                    raise ScalarReductionError(f"Singular or unresolved {label}; no regularization")
                diagnostics[f"{label}_condition"] = high / low
                return matrix**-1

            H = self._evaluate("background_H", self.background_hessian, local)
            rhs = self._evaluate("background_rhs", self.background_rhs, local)
            Hi = inverse(H, "background_hessian")
            acceleration = Hi * rhs
            local.update({self.symbols[f"{name}_acc"]: acceleration[i] for i, name in enumerate(("U", "V", "W"))})
            Hdot = self._evaluate("background_Hdot", self.background_hessian_dot, local)
            rhsdot = self._evaluate("background_rhsdot", self.background_rhs_dot, local)
            jerk = Hi * (rhsdot - Hdot * acceleration)
            local.update({self.symbols[f"{name}_jerk"]: jerk[i] for i, name in enumerate(("U", "V", "W"))})
            diagnostics["background_acceleration_residual"] = mp.norm(H * acceleration - rhs) / max(1, mp.norm(rhs))
            diagnostics["background_jerk_residual"] = mp.norm(H * jerk + Hdot * acceleration - rhsdot) / max(
                1, mp.norm(rhsdot)
            )
            energy = self._evaluate("energy", self.energy, local)
            energy_terms = sp.Matrix(sp.Add.make_args(sp.expand(self.energy)))
            energy_scale = max(1, sum(abs(term) for term in self._evaluate("energy_terms", energy_terms, local)))
            diagnostics["energy_residual"] = abs(energy) / energy_scale
            left = self._evaluate("tadpole_left", self.tadpole_left, local)
            right = self._evaluate("tadpole_right", self.tadpole_right, local)
            diagnostics["tadpole_residual"] = max(abs(l - r) / max(1, abs(l), abs(r)) for l, r in zip(left, right))
            gates = self._evaluate("auxiliary_gates", self.auxiliary_gates, local)
            diagnostics["auxiliary_gate_residual"] = max(abs(item) for item in gates)
            if any(
                diagnostics[name] > tol for name in ("energy_residual", "tadpole_residual", "auxiliary_gate_residual")
            ):
                raise ScalarReductionError("Background energy/tadpole or auxiliary Hessian/symplectic gate failed")
            blocks = {name: self._evaluate(f"block_{name}", matrix, local) for name, matrix in self.blocks.items()}
            derivatives = {
                name: self._evaluate(f"dot_{name}", matrix, local) for name, matrix in self.block_derivatives.items()
            }
            for name in ("Kd", "Dd", "Ay"):
                if mp.norm(blocks[name] - blocks[name].T) > tol * max(1, mp.norm(blocks[name])):
                    raise ScalarReductionError(f"Nonsymmetric {name}")
            Ai = inverse(blocks["Ay"], "shift_block")
            Kd, Bd, Dd, C, E = (blocks[name] for name in ("Kd", "Bd", "Dd", "C", "E"))
            K = Kd - C * Ai * C.T
            B = Bd - C * Ai * E.T
            D = Dd - E * Ai * E.T
            Aidot = -Ai * derivatives["Ay"] * Ai
            Kdot = derivatives["Kd"] - derivatives["C"] * Ai * C.T - C * Aidot * C.T - C * Ai * derivatives["C"].T
            Bdot = derivatives["Bd"] - derivatives["C"] * Ai * E.T - C * Aidot * E.T - C * Ai * derivatives["E"].T
            Ki = inverse(K, "reduced_kinetic")
            raw = {name: self._evaluate(f"raw_{name}", matrix, local) for name, matrix in self.raw_blocks.items()}
            return ScalarOperator(self, dps, local, blocks, derivatives, raw, K, B, D, Kdot, Bdot, Ki, Ai, diagnostics)


@dataclass
class ScalarReconstruction:
    shift: Any
    shift_dot: Any
    dust_density_perturbation: Any
    raw_position: Any
    raw_velocity: Any
    raw_acceleration: Any


@dataclass
class ScalarOperator:
    """Evaluated full local operator; coefficients evolve along the solved background."""

    system: ScalarSystem
    dps: int
    background: dict
    blocks: dict
    derivatives: dict
    raw: dict
    kinetic: Any
    mixed: Any
    potential: Any
    kinetic_dot: Any
    mixed_dot: Any
    kinetic_inverse: Any
    auxiliary_inverse: Any
    diagnostics: dict

    @property
    def damping(self):
        with mp.workdps(self.dps):
            return self.kinetic_dot + self.mixed - self.mixed.T

    @property
    def stiffness(self):
        with mp.workdps(self.dps):
            return self.mixed_dot - self.potential

    def acceleration(self, position, velocity):
        with mp.workdps(self.dps):
            return -self.kinetic_inverse * (self.damping * mp.matrix(velocity) + self.stiffness * mp.matrix(position))

    def reconstruct(self, position, velocity, acceleration=None, *, shift_acceleration=None):
        """Reconstruct S,T, their rates, and dlambda from the original lapse equation.

        Dynamic accelerations default to the reduced equation. Shift second
        derivatives are immaterial to the retained equations (their Hessian
        columns vanish); if the lapse reconstruction needs them, the caller
        must supply them explicitly rather than silently freeze a shift.
        """
        with mp.workdps(self.dps):
            x, xd = mp.matrix(position), mp.matrix(velocity)
            xdd = self.acceleration(x, xd) if acceleration is None else mp.matrix(acceleration)
            C, E = self.blocks["C"], self.blocks["E"]
            y = -self.auxiliary_inverse * (C.T * xd + E.T * x)
            yd = -self.auxiliary_inverse * (
                self.derivatives["C"].T * xd
                + C.T * xdd
                + self.derivatives["E"].T * x
                + E.T * xd
                + self.derivatives["Ay"] * y
            )
            system = self.system
            q, dq, ddq = (
                mp.matrix(len(system.fields), 1),
                mp.matrix(len(system.fields), 1),
                mp.matrix(len(system.fields), 1),
            )
            for i, raw_index in enumerate(system.dynamic_indices):
                q[raw_index], dq[raw_index], ddq[raw_index] = x[i], xd[i], xdd[i]
            required = any(system.dust_reconstruction.has(system.accelerations[i]) for i in system.auxiliary_indices)
            if required and shift_acceleration is None:
                raise ScalarReductionError(
                    "Dust reconstruction requires shift second derivatives; supply them explicitly"
                )
            for i, raw_index in enumerate(system.auxiliary_indices):
                q[raw_index], dq[raw_index] = y[i], yd[i]
                ddq[raw_index] = 0 if shift_acceleration is None else shift_acceleration[i]
            local = dict(self.background)
            local.update(zip(system.fields, q))
            local.update(zip(system.velocities, dq))
            local.update(zip(system.accelerations, ddq))
            density = system._evaluate("dust_reconstruction", system.dust_reconstruction, local)
            return ScalarReconstruction(y, yd, density, q, dq, ddq)

    def retained_residual(self, reconstruction):
        """Original seven retained EL equations evaluated on reconstructed fields."""
        with mp.workdps(self.dps):
            return (
                self.raw["K"] * reconstruction.raw_acceleration
                + (self.raw["Kdot"] + self.raw["B"] - self.raw["B"].T) * reconstruction.raw_velocity
                + (self.raw["Bdot"] - self.raw["D"]) * reconstruction.raw_position
            )


def constrained_scalar_system(action) -> ScalarSystem:
    """Reduce a complete ``quadratic_action('scalar', include_matter=True)`` result.

    No alpha/beta or interaction-function choice is imposed here. Generate the
    full action separately so scientific work is explicit and schedulable.
    Only dust-clock and ordinary shear gauges are fixed. The dust multiplier
    imposes phi=0; the omitted phi equation reconstructs its density multiplier.
    """
    expected = ("phi", "S", "zeta", "ell", "phihat", "T", "zetahat", "ellhat", "tau", "dlambda", "chi")
    if action.sector != "scalar" or not action.include_matter or tuple(map(str, action.fields)) != expected:
        raise ValueError("A complete scalar action with the agreed dust/radiation fields is required")
    sy = dict(action.symbols)
    sy.update(
        {
            name: sp.Symbol(name, real=True)
            for name in ("Mthird", "Mfourth", "U_acc", "V_acc", "W_acc", "U_jerk", "V_jerk", "W_jerk")
        }
    )
    sy["k"] = sp.Symbol("k", positive=True)
    a, b, X, U, V, W = (sy[name] for name in ("a", "b", "X", "U", "V", "W"))
    rho, rate, kap = (sy[name] for name in ("dust_density", "radiation_rate", "radiation_kappa"))
    z = action.Z0
    positions, rates = (a, b, X), (U, V, W)
    acc = tuple(sy[f"{name}_acc"] for name in ("U", "V", "W"))
    jerk = tuple(sy[f"{name}_jerk"] for name in ("U", "V", "W"))
    flow = dict(zip(positions + rates + acc, (a * U, b * V, X * W) + acc + jerk))
    flow.update({rho: -3 * U * rho, rate: -U * rate})
    zdot = sp.expand(sum(sp.diff(z, variable) * speed for variable, speed in flow.items()))
    interaction_names = ("Mbar", "Mprime", "Msecond", "Mthird", "Mfourth")
    flow.update({sy[left]: sy[right] * zdot for left, right in pairwise(interaction_names)})

    def dt(expression):
        return sp.expand(sum(sp.diff(expression, variable) * speed for variable, speed in flow.items()))

    def derivative(expression, variable):
        return sp.diff(expression, variable) + sp.diff(z, variable) * sum(
            sp.diff(expression, sy[left]) * sy[right] for left, right in pairwise(interaction_names)
        )

    Lg = sp.expand(action.lagrangian0 - a**3 * kap * rate**4 / 4)
    momenta = [derivative(Lg, velocity) for velocity in rates]
    H = sp.Matrix([[derivative(momentum, velocity) for velocity in rates] for momentum in momenta])
    rhs = sp.Matrix(
        [
            positions[i] * derivative(Lg, positions[i])
            - sum(positions[j] * derivative(momenta[i], positions[j]) * rates[j] for j in range(3))
            + (3 * a**3 * kap * rate**4 / 4 if i == 0 else 0)
            for i in range(3)
        ]
    )
    energy = sp.expand(
        sum(velocity * momentum for velocity, momentum in zip(rates, momenta))
        - Lg
        + a**3 * (rho + 3 * kap * rate**4 / 4)
    )
    field_by_name = {str(item): item for item in action.fields}
    dots = dict(zip(action.fields, action.time_derivatives))
    spatial = dict(zip(action.fields, action.spatial_derivatives))
    tadpole_left = sp.Matrix([sp.diff(action.lagrangian1, item) for item in action.fields])
    tadpole_right = sp.Matrix([dt(sp.diff(action.lagrangian1, dots[item])) for item in action.fields])
    gauge = {}
    tau = field_by_name["tau"]
    gauge.update({tau: 0, dots[tau]: 0, spatial[tau]: 0})
    ell, zeta = field_by_name["ell"], field_by_name["zeta"]
    gauge.update({ell: zeta, dots[ell]: dots[zeta], spatial[ell]: spatial[zeta]})
    multiplier = sp.expand(sp.diff(action.lagrangian2, field_by_name["dlambda"]).subs(gauge, simultaneous=True))
    if sp.expand(multiplier + a**3 * field_by_name["phi"]) != 0:
        raise ScalarReductionError("Dust-time multiplier does not impose phi=0")
    full = action.lagrangian2.subs(gauge, simultaneous=True)
    full_names = ("phi", "zeta", "S", "phihat", "T", "zetahat", "ellhat", "dlambda", "chi")
    cosine, sine = sp.symbols("mode_cosine mode_sine", real=True)
    fourier = {}
    for name in full_names:
        item = field_by_name[name]
        odd = name in ("S", "T")
        fourier.update(
            {
                item: item * (sine if odd else cosine),
                dots[item]: dots[item] * (sine if odd else cosine),
                spatial[item]: sy["k"] * item * cosine if odd else -sy["k"] * item * sine,
            }
        )
    full = sp.expand(full.subs(fourier, simultaneous=True))
    full = sp.expand((full.subs({cosine: 1, sine: 0}) + full.subs({cosine: 0, sine: 1})) / 2)
    full_q = tuple(field_by_name[name] for name in full_names)
    full_dq = tuple(dots[item] for item in full_q)
    full_ddq = tuple(sp.Symbol(f"{item}_ddot", real=True) for item in full_q)

    def total_dt(expression):
        return dt(expression) + sum(
            sp.diff(expression, item) * speed for item, speed in zip(full_q + full_dq, full_dq + full_ddq)
        )

    phi, dlambda = field_by_name["phi"], field_by_name["dlambda"]
    lapse_equation = sp.expand(total_dt(sp.diff(full, dots[phi])) - sp.diff(full, phi))
    density_coefficient = sp.diff(lapse_equation, dlambda)
    if density_coefficient == 0:
        raise ScalarReductionError("Lapse equation cannot reconstruct the dust density multiplier")
    eliminate = {item: 0 for item in (phi, dlambda, dots[phi], dots[dlambda], full_ddq[0], full_ddq[7])}
    density_expression = sp.expand(-lapse_equation.subs(eliminate, simultaneous=True) / density_coefficient)
    reconstructed_lapse = lapse_equation.subs({item: value for item, value in eliminate.items() if item != dlambda})
    L = sp.expand(full.subs(eliminate, simultaneous=True))
    names = ("zeta", "S", "phihat", "T", "zetahat", "ellhat", "chi")
    q = tuple(field_by_name[name] for name in names)
    dq = tuple(dots[item] for item in q)
    ddq = tuple(sp.Symbol(f"{item}_ddot", real=True) for item in q)
    K = sp.hessian(L, dq)
    B = sp.Matrix([[sp.diff(L, speed, item) for item in q] for speed in dq])
    D = sp.hessian(L, q)
    dyn, aux = (0, 2, 4, 5, 6), (1, 3)
    Baa = B.extract(aux, aux)
    gates = sp.Matrix(
        [K[i, j] for i in aux for j in range(7)] + [K[j, i] for i in aux for j in range(7)] + list(Baa - Baa.T)
    )
    blocks = {
        "Kd": K.extract(dyn, dyn),
        "Bd": B.extract(dyn, dyn),
        "Dd": D.extract(dyn, dyn),
        "C": B.extract(dyn, aux) - B.extract(aux, dyn).T,
        "E": D.extract(dyn, aux) - B.extract(aux, dyn).T.applyfunc(dt),
        "Ay": D.extract(aux, aux) - Baa.applyfunc(dt),
    }
    raw = {"K": K, "B": B, "D": D, "Kdot": K.applyfunc(dt), "Bdot": B.applyfunc(dt)}
    return ScalarSystem(
        sy,
        q,
        dq,
        ddq,
        dyn,
        aux,
        L,
        raw,
        blocks,
        {name: matrix.applyfunc(dt) for name, matrix in blocks.items()},
        H,
        rhs,
        H.applyfunc(dt),
        rhs.applyfunc(dt),
        energy,
        z,
        flow,
        gates,
        tadpole_left,
        tadpole_right,
        density_expression,
        reconstructed_lapse,
        dlambda,
    )
