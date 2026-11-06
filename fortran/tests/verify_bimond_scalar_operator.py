"""Scheduled validation/export of the complete local perfect-fluid scalar operator.

Uses only the two already-audited alpha=beta=1 local backgrounds, with k=10,
100,1000. This is not a background survey or a transfer-function calculation.
Load both source modules without importing CAMB's compiled extension. Outputs
are created exclusively; old reports are never overwritten.
"""

import argparse
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

import mpmath as mp
import sympy as sp


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def established_fixture(system, sign):
    sy = system.symbols
    magnitude = sp.Symbol("fixture_abs_Z", positive=True)
    y = 4 * (sp.sqrt(2) - 1)
    primitive = sign * (4 * sp.sqrt(magnitude) - 16 * sp.log(1 + sp.sqrt(magnitude) / 4))
    values = {
        sy["alpha"]: 1,
        sy["beta"]: 1,
        sy["a0bar"]: 1,
        sy["a"]: 1,
        sy["b"]: 1 / sp.sqrt(2),
        sy["X"]: sp.sqrt(2),
        sy["U"]: 1,
        sy["V"]: sp.Rational(1, 2),
        sy["radiation_rate"]: 1,
        sy["radiation_kappa"]: sp.Rational(2, 75),
    }
    for order, name in enumerate(("Mbar", "Mprime", "Msecond", "Mthird", "Mfourth")):
        values[sy[name]] = sign**order * sp.diff(primitive, magnitude, order).subs(magnitude, y**2)
    w = sp.solve(system.z.subs(values) - sign * y**2, sy["W"])
    if len(w) != 1:
        raise RuntimeError("Fixture no longer uniquely specifies W")
    values[sy["W"]] = w[0]
    density = sp.solve(system.energy.subs(values), sy["dust_density"])
    if len(density) != 1:
        raise RuntimeError("Fixture no longer uniquely specifies the dust constraint")
    values[sy["dust_density"]] = density[0]
    return values


def assert_small(value, tolerance, label):
    if not mp.isfinite(value) or value > tolerance:
        raise AssertionError(f"{label}: {mp.nstr(value, 12)} > {tolerance}")


def direct_euler_lagrange(system):
    """Differentiate the raw scalar action, independently of the Schur blocks."""
    residuals = []
    for coordinate, velocity in zip(system.fields, system.velocities):
        momentum = sp.diff(system.lagrangian, velocity)
        total = system.dt(momentum)
        total += sum(
            sp.diff(momentum, variable) * speed
            for variable, speed in zip(system.fields + system.velocities, system.velocities + system.accelerations)
        )
        residuals.append(sp.expand(total - sp.diff(system.lagrangian, coordinate)))
    return sp.Matrix(residuals)


def derivative_check(system, operator):
    """Central difference along the evolving EOM tangent, not a frozen-rate test.

    Nearby coefficient evaluations are algebraic directional probes, not new
    on-shell backgrounds or numerical evolution. acc_dot is the solved jerk;
    all M jets, dust density, and radiation rate use their complete chain rules.
    """
    local = operator.background
    direction = {
        variable: system._evaluate(f"flow_{variable}", speed, local)
        for variable, speed in system.background_flow.items()
    }
    step = mp.mpf("1e-18")
    results = []
    for sign in (-1, 1):
        nearby = dict(local)
        for variable, speed in direction.items():
            nearby[variable] += sign * step * speed
        blocks = {
            name: system._evaluate(f"block_{name}", expression, nearby) for name, expression in system.blocks.items()
        }
        Ai = blocks["Ay"] ** -1
        C, E = blocks["C"], blocks["E"]
        results.append((blocks["Kd"] - C * Ai * C.T, blocks["Bd"] - C * Ai * E.T))
    differences = [(results[1][i] - results[0][i]) / (2 * step) for i in range(2)]
    exact = (operator.kinetic_dot, operator.mixed_dot)
    return [mp.norm(numeric - analytic) / max(1, mp.norm(analytic)) for numeric, analytic in zip(differences, exact)]


def audit(system, dps):
    q, dq, ddq = map(sp.Matrix, (system.fields, system.velocities, system.accelerations))
    raw = system.raw_blocks
    reconstructed_action = (dq.T * raw["K"] * dq)[0] / 2 + (dq.T * raw["B"] * q)[0] + (q.T * raw["D"] * q)[0] / 2
    if sp.expand(reconstructed_action - system.lagrangian) != 0:
        raise AssertionError("K/B/D do not reproduce the full finite-k action")
    direct = direct_euler_lagrange(system)
    matrix_residual = raw["K"] * ddq + (raw["Kdot"] + raw["B"] - raw["B"].T) * dq + (raw["Bdot"] - raw["D"]) * q
    if any(sp.expand(item) != 0 for item in direct - matrix_residual):
        raise AssertionError("Raw matrix equations differ from direct Euler-Lagrange variation")
    records = []
    tol = mp.mpf("1e-40")
    for sign in (1, -1):
        fixture = established_fixture(system, sign)
        for wave in (10, 100, 1000):
            values = dict(fixture)
            values[system.symbols["k"]] = wave
            operator = system.evaluate(values, dps=dps)
            # Nonzero deterministic probes exercise position, velocity, and
            # acceleration terms simultaneously; they are not initial modes.
            x = mp.matrix([mp.mpf(i + 1) / 7 for i in range(5)])
            xd = mp.matrix([(-1) ** i * mp.mpf(i + 2) / 11 for i in range(5)])
            xdd = operator.acceleration(x, xd)
            reconstruction = operator.reconstruct(x, xd, xdd)
            local = dict(operator.background)
            local.update(zip(system.fields, reconstruction.raw_position))
            local.update(zip(system.velocities, reconstruction.raw_velocity))
            local.update(zip(system.accelerations, reconstruction.raw_acceleration))
            local[system.dust_symbol] = reconstruction.dust_density_perturbation
            full_residual = system._evaluate("direct_el_validation", direct, local)
            retained = operator.retained_residual(reconstruction)
            raw_terms = (
                operator.raw["K"] * reconstruction.raw_acceleration,
                (operator.raw["Kdot"] + operator.raw["B"] - operator.raw["B"].T) * reconstruction.raw_velocity,
                (operator.raw["Bdot"] - operator.raw["D"]) * reconstruction.raw_position,
            )
            scale = max(1, *(mp.norm(term) for term in raw_terms))
            retained_error = mp.norm(full_residual) / scale
            matrix_error = mp.norm(full_residual - retained) / scale
            assert_small(retained_error, tol, "Original retained equations after reconstruction")
            assert_small(matrix_error, tol, "Direct versus matrix retained residual")
            lapse = system._evaluate("lapse_validation", system.lapse_equation, local)
            lapse_terms = sp.Matrix(sp.Add.make_args(sp.expand(system.lapse_equation)))
            lapse_scale = max(
                1, sum(abs(item) for item in system._evaluate("lapse_terms_validation", lapse_terms, local))
            )
            lapse_error = abs(lapse) / lapse_scale
            assert_small(lapse_error, tol, "Original lapse/dust-density reconstruction")
            blocks, y = operator.blocks, reconstruction.shift
            original_reduced_action = (
                (xd.T * blocks["Kd"] * xd)[0] / 2
                + (xd.T * blocks["Bd"] * x)[0]
                + (x.T * blocks["Dd"] * x)[0] / 2
                + (xd.T * blocks["C"] * y)[0]
                + (x.T * blocks["E"] * y)[0]
                + (y.T * blocks["Ay"] * y)[0] / 2
            )
            physical_action = (
                (xd.T * operator.kinetic * xd)[0] / 2
                + (xd.T * operator.mixed * x)[0]
                + (x.T * operator.potential * x)[0] / 2
            )
            action_error = abs(original_reduced_action - physical_action) / max(
                1, abs(original_reduced_action), abs(physical_action)
            )
            assert_small(action_error, tol, "Full Schur action identity")
            derivative_errors = derivative_check(system, operator)
            for error in derivative_errors:
                assert_small(error, mp.mpf("1e-26"), "Evolving-background coefficient derivative")
            eigenvalues = mp.eigsy(operator.kinetic, eigvals_only=True)
            if sum(value < 0 for value in eigenvalues) != 2:
                raise AssertionError("Established local fixture no longer has its two scalar ghost directions")
            # The public gate must reject an off-shell supplied matter density.
            invalid = dict(values)
            invalid[system.symbols["dust_density"]] *= sp.Rational(101, 100)
            try:
                system.evaluate(invalid, dps=dps)
            except ValueError:
                pass
            else:
                raise AssertionError("Off-shell local background was silently accepted")

            def rows(matrix):
                return [[mp.nstr(matrix[i, j], dps) for j in range(matrix.cols)] for i in range(matrix.rows)]

            records.append(
                {
                    "Z_sign": sign,
                    "k": wave,
                    "background": {
                        str(variable): mp.nstr(value, dps) for variable, value in operator.background.items()
                    },
                    "diagnostics": {name: mp.nstr(value, dps) for name, value in operator.diagnostics.items()},
                    "Keff": rows(operator.kinetic),
                    "Beff": rows(operator.mixed),
                    "Deff": rows(operator.potential),
                    "Keff_dot": rows(operator.kinetic_dot),
                    "Beff_dot": rows(operator.mixed_dot),
                    "damping": rows(operator.damping),
                    "stiffness": rows(operator.stiffness),
                    "negative_kinetic_directions": sum(value < 0 for value in eigenvalues),
                    "retained_EL_relative_residual": mp.nstr(retained_error, dps),
                    "matrix_EL_relative_residual": mp.nstr(matrix_error, dps),
                    "lapse_reconstruction_relative_residual": mp.nstr(lapse_error, dps),
                    "Schur_action_relative_residual": mp.nstr(action_error, dps),
                    "coefficient_derivative_relative_errors": [mp.nstr(value, dps) for value in derivative_errors],
                    "probe_shift": [mp.nstr(value, dps) for value in reconstruction.shift],
                    "probe_shift_dot": [mp.nstr(value, dps) for value in reconstruction.shift_dot],
                    "probe_dust_density_perturbation": mp.nstr(reconstruction.dust_density_perturbation, dps),
                }
            )
    return {
        "raw_action_identity": True,
        "direct_EL_matrix_identity": True,
        "dynamic_field_order": list(map(str, system.dynamic_fields)),
        "auxiliary_field_order": list(map(str, system.auxiliary_fields)),
        "points": records,
        "limitations": [
            "Two established local backgrounds only, both with scalar ghosts",
            "No time integration, initial-mode prescription, or global history",
            "No photon/neutrino hierarchy or cosmological transfer functions",
            "Finite differences are local directional checks, not new on-shell background samples",
            "No inverse was regularized; background and auxiliary gates were enforced",
        ],
    }


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--coefficients", type=Path)
    parser.add_argument("--generator", type=Path, default=root / "camb" / "bimond_symbolic.py")
    parser.add_argument("--operator", type=Path, default=root / "camb" / "bimond_scalar.py")
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()
    if args.dps < 70:
        parser.error("Use at least 70 digits for the independent cancellation-sensitive audit")
    report = {"scope": "complete constrained local scalar operator; not a production transfer model", "dps": args.dps}
    with args.output.open("x") as stream:
        try:
            report["source_sha256"] = {
                str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                for path in (args.generator, args.operator, Path(__file__))
            }
            generator = load_module("standalone_bimond_symbolic", args.generator)
            implementation = load_module("standalone_bimond_scalar", args.operator)
            with mp.workdps(args.dps):
                system = implementation.constrained_scalar_system(
                    generator.quadratic_action("scalar", include_matter=True)
                )
                if args.coefficients:
                    with args.coefficients.open("x") as coefficients:
                        json.dump(system.coefficient_export(), coefficients, indent=2)
                        coefficients.write("\n")
                report.update(audit(system, args.dps))
        except Exception:
            report["fatal_error"] = traceback.format_exc()
            json.dump(report, stream, indent=2)
            stream.write("\n")
            raise
        json.dump(report, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
