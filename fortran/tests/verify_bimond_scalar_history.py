"""Finite-interval, high-precision BIMOND history verification on allocated CPUs.

Two established ghost-bearing backgrounds only, k=10, tau in [0,0.001].
Fundamental identity data are mathematical basis responses, not primordial
initial modes. Source modules are loaded without the CAMB compiled extension.
"""

import argparse
import hashlib
import importlib.util
import json
import sys
import time
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


def norms(left, right):
    absolute = mp.norm(left - right)
    return {"absolute": absolute, "relative": absolute / max(mp.mpf("1e-100"), mp.norm(right))}


def require(value, tolerance, name):
    if not mp.isfinite(value) or value > mp.mpf(tolerance):
        raise AssertionError(f"{name}: {mp.nstr(value, 12)} exceeds {tolerance}")


def rows(matrix, dps):
    return [[mp.nstr(matrix[i, j], dps) for j in range(matrix.cols)] for i in range(matrix.rows)]


def serialized(value, dps):
    if isinstance(value, dict):
        return {str(key): serialized(item, dps) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialized(item, dps) for item in value]
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    return mp.nstr(value, dps)


def trajectory_shift_second(module, history, when, delta):
    """Independent time difference of shift rates on the computed history."""

    def value(time):
        fundamental, _ = history.evaluate(time)
        return module.shift_and_rate(history.background.operator(time), fundamental)[1]

    return (value(when - 2 * delta) - 8 * value(when - delta) + 8 * value(when + delta) - value(when + 2 * delta)) / (
        12 * delta
    )


def run_case(system, implementation, fixture, sign, dps):
    with mp.workdps(dps + 20):
        exact = {str(key): value for key, value in fixture.items() if not str(key).startswith("M")}
        base = implementation.BackgroundHistory(
            system, exact, k=10, end_time="0.001", dps=dps, background_tolerance=f"1e-{dps - 15}"
        )
        counts = (2, 4, 8)
        histories = []
        for count in counts:
            print(f"branch={sign:+d} dps={dps} scalar_steps={count}", flush=True)
            histories.append(implementation.ScalarHistory(base, steps=count, stages=6))
        high = implementation.BackgroundHistory(
            system, exact, k=10, end_time="0.001", dps=dps + 20, background_tolerance=f"1e-{dps + 5}"
        )
        print(f"branch={sign:+d} dps={dps + 20} scalar_steps=8", flush=True)
        high_history = implementation.ScalarHistory(high, steps=8, stages=6)
        end = base.end_time
        endpoint_errors = [
            norms(history.endpoint_matrices[-1], histories[-1].endpoint_matrices[-1]) for history in histories[:-1]
        ]
        precision_error = norms(histories[-1].endpoint_matrices[-1], high_history.endpoint_matrices[-1])
        require(endpoint_errors[-1]["relative"], "1e-18", "Scalar endpoint step convergence")
        if endpoint_errors[-1]["absolute"] > endpoint_errors[0]["absolute"] / 32:
            raise AssertionError("Endpoint error did not decrease under step refinement")
        require(precision_error["relative"], "1e-55", "Scalar endpoint precision convergence")
        background_comparisons = []
        for numerator in range(9):
            when_base = base.end_time * numerator / 8
            when_high = high.end_time * numerator / 8
            comparison = norms(mp.matrix(base.state(when_base)), mp.matrix(high.state(when_high)))
            background_comparisons.append(comparison)
            require(comparison["absolute"], "1e-55", "Independent Taylor precision/tolerance refinement")
        # Compare Taylor backgrounds at every actual scalar collocation point,
        # not only the requested output grid. The lower-precision evaluations
        # have already passed all on-shell/branch/conditioning gates.
        collocation_background_error = mp.mpf(0)
        for when in list(base._operator_cache):
            when_high = high.end_time if when == base.end_time else high.end_time * (when / base.end_time)
            comparison = norms(mp.matrix(base.state(when)), mp.matrix(high.state(when_high)))
            collocation_background_error = max(collocation_background_error, comparison["absolute"])
            require(comparison["absolute"], "1e-55", "Taylor refinement at scalar collocation/query point")
            high.operator(when_high)
        defects = [history.maximum_dense_defect_relative for history in histories]
        require(defects[-1], "1e-12", "Dense scalar polynomial differential defect")
        if defects[-1] > defects[0] / 16:
            raise AssertionError("Dense differential defect failed to converge")
        when = end / 2
        shift_second = []
        directional_errors = []
        time_comparisons = []
        for history in histories:
            second, error = implementation.shift_second_derivative(history, when, difference_step="1e-28", extra_dps=40)
            refined, refined_error = implementation.shift_second_derivative(
                history, when, difference_step="5e-29", extra_dps=40
            )
            require(norms(second, refined)["relative"], "1e-60", "Directional shift derivative step refinement")
            shift_second.append(refined)
            directional_errors.append(
                {"h_vs_half_h": norms(second, refined), "coarse_internal": error, "fine_internal": refined_error}
            )
            time_difference = trajectory_shift_second(implementation, history, when, mp.mpf("1e-8"))
            time_refined = trajectory_shift_second(implementation, history, when, mp.mpf("5e-9"))
            time_comparisons.append(
                {
                    "trajectory_vs_direction": norms(time_refined, refined),
                    "trajectory_difference_step_check": norms(time_difference, time_refined),
                }
            )
        high_second, high_derivative_error = implementation.shift_second_derivative(
            high_history, high.end_time / 2, difference_step="1e-32", extra_dps=40
        )
        shift_precision = norms(shift_second[-1], high_second)
        require(shift_precision["relative"], "1e-50", "Directional shift derivative precision refinement")
        require(
            time_comparisons[-1]["trajectory_vs_direction"]["relative"],
            "1e-11",
            "Trajectory check of shift second derivative",
        )
        if (
            time_comparisons[-1]["trajectory_vs_direction"]["absolute"]
            > time_comparisons[0]["trajectory_vs_direction"]["absolute"] / 8
        ):
            raise AssertionError("Trajectory shift-second-derivative error failed to converge")
        samples = []
        for numerator in (0, 4, 8):
            time = end * numerator / 8
            second, derivative_error = implementation.shift_second_derivative(
                histories[-1], time, difference_step="5e-29", extra_dps=40
            )
            residuals, density = implementation.reconstruction_residuals(histories[-1], time, second)
            require(residuals["retained_relative"], "1e-55", "Original retained equations along history")
            require(residuals["lapse_relative"], "1e-55", "Original lapse reconstruction along history")
            operator = base.operator(time)
            eigenvalues = mp.eigsy(operator.kinetic, eigvals_only=True)
            negatives = sum(value < 0 for value in eigenvalues)
            if negatives != 2:
                raise AssertionError("The established diagnostic's two ghost directions changed")
            fundamental, _ = histories[-1].evaluate(time)
            shift, rate = implementation.shift_and_rate(operator, fundamental)
            samples.append(
                {
                    "tau": mp.nstr(time, dps),
                    "background_log_state": [mp.nstr(value, dps) for value in base.state(time)],
                    "Z": mp.nstr(system._evaluate("history_Z", system.z, operator.background), dps),
                    "fundamental_matrix": rows(fundamental, dps),
                    "shift": rows(shift, dps),
                    "shift_dot": rows(rate, dps),
                    "shift_ddot": rows(second, dps),
                    "dust_density_basis_response": rows(density, dps),
                    "reconstruction_residuals": residuals,
                    "shift_derivative_diagnostics": derivative_error,
                    "negative_kinetic_directions": negatives,
                    "kinetic_eigenvalues": list(eigenvalues),
                    "energy_constraint": system._evaluate("energy", system.energy, operator.background),
                    "operator_diagnostics": operator.diagnostics,
                    "dust_density": operator.background[system.symbols["dust_density"]],
                    "radiation_rate": operator.background[system.symbols["radiation_rate"]],
                }
            )
        # Real guard failures, not smoothing, are required for invalid requests.
        wrong = dict(exact)
        wrong["U"] = 1.0
        try:
            implementation.BackgroundHistory(system, wrong, k=10, end_time="0.001", dps=dps)
        except TypeError:
            pass
        else:
            raise AssertionError("Binary64 initial condition was accepted")
        bad = dict(exact)
        bad["dust_density"] *= sp.Rational(101, 100)
        try:
            implementation.BackgroundHistory(system, bad, k=10, end_time="0.001", dps=dps)
        except ValueError:
            pass
        else:
            raise AssertionError("Off-shell initial background was accepted")
        try:
            implementation.BackgroundHistory(system, exact, k=10, end_time="0.001", dps=dps, max_condition="1.01")
        except ValueError:
            pass
        else:
            raise AssertionError("Ill-conditioned background/operator was accepted")
        # Z is affine in W at the established state; choose the other sign's W
        # algebraically, without an off-shell evolution or background survey.
        initial_values, initial_z = base.values_from_state(base.initial_state)
        w = system.symbols["W"]
        coefficient = system._evaluate("test_dZ_dW", sp.diff(system.z, w), initial_values)
        displaced = list(base.initial_state)
        displaced[5] -= 2 * initial_z / coefficient
        try:
            base.values_from_state(displaced)
        except ValueError:
            pass
        else:
            raise AssertionError("A branch-changing background state was accepted")
        return {
            "Z_sign": sign,
            "k": 10,
            "tau_span": ["0", "0.001"],
            "dps": dps,
            "comparison_dps": dps + 20,
            "step_counts": counts,
            "Gauss_stages": 6,
            "endpoint_order": 12,
            "scalar_endpoint_step_errors": endpoint_errors,
            "scalar_endpoint_precision_error": precision_error,
            "background_output_grid_refinement": background_comparisons,
            "maximum_background_collocation_refinement_absolute": collocation_background_error,
            "background_taylor_tolerance": base.background_tolerance,
            "comparison_background_taylor_tolerance": high.background_tolerance,
            "maximum_background_constraint_residual": base.maximum_constraint,
            "maximum_background_tadpole_residual": base.maximum_tadpole,
            "maximum_matrix_conditions": base.maximum_conditions,
            "maximum_collocation_linear_solve_relative_residual": max(
                history.maximum_linear_solve_residual for history in histories
            ),
            "dense_scalar_defect_absolute": [history.maximum_dense_defect_absolute for history in histories],
            "dense_scalar_defect_relative": defects,
            "directional_shift_second_step_checks": directional_errors,
            "directional_shift_second_precision_check": shift_precision,
            "high_precision_shift_second_diagnostics": high_derivative_error,
            "trajectory_shift_second_checks": time_comparisons,
            "background_rhs_evaluations": base.background_evaluations,
            "operator_evaluations": base.operator_evaluations,
            "gated_collocation_and_query_points": len(base._operator_cache),
            "samples": samples,
            "negative_tests": ["binary64_seed", "off_shell_dust_density", "conditioning", "Z_branch_change"],
        }


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()
    if args.dps < 80:
        parser.error("This bounded precision-convergence validation starts at 80 digits")
    paths = {
        "generator": root / "camb/bimond_symbolic.py",
        "operator": root / "camb/bimond_scalar.py",
        "history": root / "camb/bimond_history.py",
        "fixtures": root / "fortran/tests/verify_bimond_scalar_operator.py",
        "verification": Path(__file__).resolve(),
    }
    started = time.monotonic()
    report = {
        "scope": "finite-interval fundamental-matrix diagnostic on known ghost-bearing backgrounds",
        "physical_transfer_functions": False,
        "primordial_initial_modes": False,
        "ghost_projection": False,
        "variable_order": [
            "zeta",
            "phihat",
            "zetahat",
            "ellhat",
            "chi",
            "zeta_dot",
            "phihat_dot",
            "zetahat_dot",
            "ellhat_dot",
            "chi_dot",
        ],
        "time_unit": "tau=H_ref*t",
        "wavenumber": "dimensionless constant comoving k in cos(k*z)/sin(k*z)",
        "initial_matrix": "identity basis normalization, not a physical perturbation realization",
        "source_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths.values()},
    }
    with args.output.open("x") as stream:
        try:
            generator = load_module("history_symbolic", paths["generator"])
            operator = load_module("history_operator", paths["operator"])
            history = load_module("history_implementation", paths["history"])
            fixtures = load_module("history_fixtures", paths["fixtures"])
            system = operator.constrained_scalar_system(generator.quadratic_action("scalar", include_matter=True))
            report["cases"] = []
            for sign in (1, -1):
                report["cases"].append(
                    run_case(system, history, fixtures.established_fixture(system, sign), sign, args.dps)
                )
            report["status"] = "passed"
        except Exception:
            report["status"] = "failed"
            report["fatal_error"] = traceback.format_exc()
            report["elapsed_seconds"] = time.monotonic() - started
            json.dump(serialized(report, args.dps), stream, indent=2)
            stream.write("\n")
            raise
        report["elapsed_seconds"] = time.monotonic() - started
        json.dump(serialized(report, args.dps), stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
