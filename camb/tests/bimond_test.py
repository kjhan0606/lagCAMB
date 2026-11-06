"""Independent action, quadrature, and unprojected background checks.

Run with ``python -m unittest camb.tests.bimond_test`` on the CPU scheduler.
Set BIMOND_TEST_REPORT to a new JSON path for measured numerical diagnostics.
These local background tests do not validate perturbations or cosmological ICs.
"""

import json
import os
import unittest
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad

from camb.bimond import (
    BIMONDDomainError,
    BIMONDParameters,
    constraint_initial_velocity,
    evaluate_background,
    integrate_background,
    interaction_function,
    tensor_principal_matrices,
)

_REPORT = {}


def tearDownModule():
    destination = os.environ.get("BIMOND_TEST_REPORT")
    if destination:
        # Preserve old validation records rather than silently overwrite them.
        with Path(destination).open("x") as stream:
            json.dump(_REPORT, stream, indent=2, sort_keys=True)
            stream.write("\n")


def independent_function(z, m0):
    """Quadrature after x=t**4 removes the endpoint singularity of nu(sqrt(x))."""
    tmax = abs(z) ** 0.25
    primitive, error = quad(
        lambda r: 8 * r**2 / (np.sqrt(tmax**2 * r**2 + 4) + tmax * r),
        0,
        1,
        epsabs=2e-13,
        epsrel=2e-13,
    )
    primitive *= tmax**3
    root = np.sqrt(abs(z))
    nu_minus_one = 1 / (root * (np.sqrt(0.25 + 1 / root) + 0.5))
    second = -np.sign(z) / (4 * abs(z) ** 1.5 * np.sqrt(0.25 + 1 / root))
    return np.array([m0 + np.sign(z) * primitive, nu_minus_one, second]), error * tmax**3


def symbolic_action_oracle():
    """Differentiate the general-lapse action before setting n=ln(N)=0.

    The oracle uses an abstract function F and independent numerical values of
    F,F',F''. It does not duplicate the manually differentiated Fortran arrays.
    """
    u, v, w, U, V, W, n, a0, B, R = sp.symbols("u v w U V W n a0 B R", real=True)
    q, velocity = (u, v, w), (U, V, W)
    s = sp.exp(2 * (v - u - w))
    kinetic = -3 * sp.exp(3 * u) * U**2 + 3 * sp.exp(3 * v - w) * V**2
    amplitude = a0**2 * sp.exp(sp.Rational(3, 2) * (u + v) + w / 2)
    z = sp.Rational(3, 2) / a0**2 * (2 * U**2 - (3 + s) * U * V + (1 + s) * V**2 + (s - 1) * V * W)
    F = sp.Function("F")
    action = sp.exp(-n) * kinetic + sp.exp(n) * amplitude * F(sp.exp(-2 * n) * z) - sp.exp(n) * (B + R * sp.exp(-u))
    expressions = [action, -sp.diff(action, n), z]
    expressions += [sp.diff(action, coordinate) for coordinate in q]
    expressions += [sp.diff(action, rate, coordinate) for rate in velocity for coordinate in q]
    expressions += [sp.diff(action, rate1, rate2) for rate1 in velocity for rate2 in velocity]
    m, mp, mpp = sp.symbols("m mp mpp")
    output = []
    for expression in expressions:
        expression = expression.subs(n, 0)
        replacements = {}
        for derivative in expression.atoms(sp.Subs):
            if isinstance(derivative.expr, sp.Derivative) and derivative.expr.expr.func == F:
                order = sum(count for _, count in derivative.expr.variable_count)
                replacements[derivative] = {1: mp, 2: mpp}[order]
        expression = expression.xreplace(replacements).replace(lambda node: node.func == F, lambda node: m)
        output.append(expression)
    return sp.lambdify((*q, *velocity, a0, B, R, m, mp, mpp), output, modules="numpy", cse=True)


def regular_initial_state(branch):
    parameters = BIMONDParameters(a0bar=1.0, baryons=0.3, radiation=0.02, m0=0.0)
    q = np.array([0.0, 0.2, -0.1])
    twin_rate = 0.65 if branch == 1 else 1.0
    bracket = (1e-8, 1e6) if branch == 1 else (-1e6, -1e-8)
    velocity = constraint_initial_velocity(q, 1.0, twin_rate, parameters, branch=branch, z_bracket=bracket)
    return parameters, np.concatenate((q, velocity))


class InteractionTests(unittest.TestCase):
    def test_signed_primitive_quadrature_and_derivatives(self):
        errors = []
        derivative_errors = []
        for magnitude in np.logspace(-10, 10, 25):
            for sign in (-1, 1):
                z = sign * magnitude
                actual = interaction_function(z, zmin=1e-14)
                expected, _ = independent_function(z, 0)
                np.testing.assert_allclose(actual, expected, rtol=4e-11, atol=0)
                self.assertGreater(actual[1], 0)
                self.assertLess(sign * actual[2], 0)
                errors.append(abs(actual[0] / expected[0] - 1))
                derivative_errors.append(float(np.max(np.abs(actual[1:] / expected[1:] - 1))))
        # The additive constant does not alter either derivative on either branch.
        for z in (-2.0, 2.0):
            np.testing.assert_allclose(interaction_function(z, 0.7) - interaction_function(z), [0.7, 0, 0], atol=1e-15)
        _REPORT["interaction"] = {
            "max_relative_primitive_error": max(errors),
            "max_relative_derivative_error": max(derivative_errors),
            "abs_z_range": [1e-10, 1e10],
        }

    def test_excluded_surface(self):
        for z in (0, 1e-12, -1e-12):
            with self.assertRaises(BIMONDDomainError):
                interaction_function(z)


class ActionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.oracle = staticmethod(symbolic_action_oracle())

    def test_independent_lapse_constraint_and_action_derivatives(self):
        errors = []
        parameters = BIMONDParameters(0.8, 0.3, 0.02, 0.17)
        # W changes the sign of Z while keeping both metrics nondegenerate.
        for q, velocity in (
            ([0.1, 0.3, -0.2], [0.7, 0.4, 0.3]),
            ([0.1, 0.3, -0.2], [0.7, 0.4, -2.0]),
            ([-0.2, 0.1, 0.15], [1.1, -0.3, 0.2]),
        ):
            result = evaluate_background(q, velocity, parameters)
            functions, _ = independent_function(result.z, parameters.m0)
            expected = np.array(
                self.oracle(*q, *velocity, parameters.a0bar, parameters.baryons, parameters.radiation, *functions)
            )
            actual = np.concatenate(
                (
                    [result.lagrangian, result.constraint, result.z],
                    result.gradient,
                    result.mixed.ravel(),
                    result.hessian.ravel(),
                )
            )
            np.testing.assert_allclose(actual, expected, rtol=3e-11, atol=3e-12)
            np.testing.assert_allclose(result.hessian, result.hessian.T, rtol=1e-14, atol=1e-14)
            errors.append(float(np.max(abs(actual - expected) / np.maximum(1, abs(expected)))))
        _REPORT["action"] = {"max_scaled_symbolic_error": max(errors), "general_lapse_differentiated": True}

    def test_singular_hessian_at_s_one_even_with_nonzero_z(self):
        with self.assertRaisesRegex(BIMONDDomainError, "Hessian"):
            evaluate_background([0, 0, 0], [1, 0.3, 0.2], BIMONDParameters(1, 0.3, 0.02))


class IntegrationTests(unittest.TestCase):
    def test_both_constrained_branches_and_tolerance_refinement(self):
        for branch in (1, -1):
            parameters, initial = regular_initial_state(branch)
            samples = np.linspace(0, 0.01, 21)
            runs = []
            for tolerance in (1e-5, 1e-8, 1e-11):
                runs.append(
                    integrate_background(
                        parameters,
                        initial,
                        (samples[0], samples[-1]),
                        sample_times=samples,
                        rtol=tolerance,
                        atol=tolerance * 1e-2,
                        constraint_tolerance=1e-5,
                        max_step=0.01,
                    )
                )
            loose, medium, tight = runs
            self.assertTrue(np.all(branch * tight.z > 1e-12))
            self.assertLess(tight.max_abs_normalized_constraint, 1e-9)
            self.assertLess(tight.max_hessian_condition, 1e12)
            loose_error = float(np.max(abs(loose.states - tight.states)))
            medium_error = float(np.max(abs(medium.states - tight.states)))
            self.assertLessEqual(medium_error, loose_error + 2e-11)
            self.assertLess(medium_error, 1e-7)
            _REPORT[f"integration_branch_{branch:+d}"] = {
                "initial_state": initial.tolist(),
                "tau_span": [0, 0.01],
                "initial_z": float(tight.z[0]),
                "final_z": float(tight.z[-1]),
                "loose_vs_tight_max_state_error": loose_error,
                "medium_vs_tight_max_state_error": medium_error,
                "max_abs_normalized_constraint": [run.max_abs_normalized_constraint for run in runs],
                "max_hessian_condition": tight.max_hessian_condition,
                "nfev": [run.nfev for run in runs],
            }

    def test_reject_nonconstrained_initial_state(self):
        parameters, initial = regular_initial_state(1)
        initial[3] *= 1.01
        with self.assertRaisesRegex(BIMONDDomainError, "Initial state"):
            integrate_background(parameters, initial, (0, 0.01), sample_times=[0, 0.01], rtol=1e-9, atol=1e-11)

    def test_reject_trajectory_approaching_hessian_degeneracy(self):
        # Grammar diagnostic 1274015 located the old V=.4 fixture's s=1
        # endpoint near tau=.00320474627. Preserve that failure as a separate
        # rejection test; the regular-branch test retains the full .01 interval.
        parameters, initial = regular_initial_state(1)
        initial[3:] = constraint_initial_velocity(initial[:3], 1.0, 0.4, parameters, branch=1, z_bracket=(1e-8, 1e6))
        with self.assertRaisesRegex(BIMONDDomainError, "Hessian") as rejected:
            integrate_background(
                parameters,
                initial,
                (0, 0.01),
                sample_times=[0, 0.01],
                rtol=1e-10,
                atol=1e-12,
                max_step=0.0001,
            )
        _REPORT["degenerate_trajectory"] = {
            "initial_state": initial.tolist(),
            "requested_tau_span": [0, 0.01],
            "max_step": 0.0001,
            "rejection": str(rejected.exception),
        }

    def test_reject_singular_initial_solve_and_wrong_branch_bracket(self):
        parameters = BIMONDParameters(1, 0.3, 0.02)
        with self.assertRaisesRegex(BIMONDDomainError, "singular"):
            constraint_initial_velocity([0, 0, 0], 1, 0.5, parameters, branch=1, z_bracket=(1e-8, 1e5))
        with self.assertRaises(ValueError):
            constraint_initial_velocity([0, 0.2, -0.1], 1, 0.5, parameters, branch=-1, z_bracket=(1e-8, 1e5))
        with self.assertRaisesRegex(BIMONDDomainError, "root"):
            constraint_initial_velocity([0, 0.2, -0.1], 1, 0.4, parameters, branch=-1, z_bracket=(-1e6, -1e-8))


class TensorPrincipalTests(unittest.TestCase):
    def test_determinant_identities_and_regular_fixtures(self):
        measured = []
        for branch in (1, -1):
            parameters, initial = regular_initial_state(branch)
            background = evaluate_background(initial[:3], initial[3:], parameters)
            tensor = tensor_principal_matrices(initial[:3], background.z, parameters)
            m, s, t = tensor.m, tensor.s, tensor.t
            expected_k = -t * (1 - m) - m * (2 * s - 1) - m**2 * (s - 1) ** 2
            expected_g = -(t / s) * (1 - m) - m
            np.testing.assert_allclose(np.linalg.det(tensor.kinetic), expected_k, rtol=2e-13, atol=2e-13)
            np.testing.assert_allclose(np.linalg.det(tensor.gradient), expected_g, rtol=2e-13, atol=2e-13)
            np.testing.assert_array_equal(tensor.kinetic, tensor.kinetic.T)
            np.testing.assert_array_equal(tensor.gradient, tensor.gradient.T)
            # Positive definite K would require an indefinite gradient matrix.
            if tensor.kinetic_eigenvalues[0] > 0:
                self.assertLess(m, 1)
                self.assertLess(expected_g, 0)
            measured.append(
                {
                    "branch": branch,
                    "m": m,
                    "s": s,
                    "t": t,
                    "kinetic_eigenvalues": tensor.kinetic_eigenvalues.tolist(),
                    "kinetic_determinant": float(np.linalg.det(tensor.kinetic)),
                    "gradient_determinant": float(np.linalg.det(tensor.gradient)),
                }
            )
        _REPORT["tensor_principal_matrices"] = measured

    def test_formal_equal_metric_ratio_limit(self):
        # Treat m as a finite constant at r=X=1. This is an algebraic check,
        # not a background solution at the singular equal-metric state.
        parameters = BIMONDParameters(1, 0.3, 0.02)
        for z in (-100.0, -1.0, 0.01, 1.0, 100.0):
            tensor = tensor_principal_matrices([0, 0, 0], z, parameters)
            self.assertAlmostEqual(np.linalg.det(tensor.kinetic), -1, places=12)
            np.testing.assert_array_equal(tensor.kinetic, tensor.gradient)
            self.assertLess(tensor.kinetic_eigenvalues[0], 0)
            self.assertGreater(tensor.kinetic_eigenvalues[1], 0)


if __name__ == "__main__":
    unittest.main()
