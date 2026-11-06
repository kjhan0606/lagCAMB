"""Verify unreduced BIMOND scalar/vector actions on the designated CPU scheduler."""

import importlib.util
import json
import sys
import time
from pathlib import Path

import sympy as sp


def load_generator(path):
    spec = importlib.util.spec_from_file_location("bimond_symbolic", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def zero(expression):
    return sp.cancel(sp.expand(expression)) == 0


def store_action(action, directory):
    payload = {
        "sector": action.sector,
        "fields": list(map(str, action.fields)),
        "time_derivatives": list(map(str, action.time_derivatives)),
        "spatial_derivatives": list(map(str, action.spatial_derivatives)),
        "symbols": {
            name: {"positive": symbol.is_positive, "real": symbol.is_real} for name, symbol in action.symbols.items()
        },
        "expressions": {
            name: str(getattr(action, name)) for name in ("lagrangian0", "lagrangian1", "lagrangian2", "Z0", "Z1", "Z2")
        },
        "normalization": "Lagrangian coefficients multiply epsilon^n and equal 8*pi*G times the action density.",
        "scope": "Unreduced gravitational action, no matter perturbations or background equations imposed.",
    }
    with (directory / f"{action.sector}_action.json").open("x") as stream:
        json.dump(payload, stream, indent=2)
        stream.write("\n")


def check_background(action):
    a, b, X, a0, U, V, W, alpha, beta, M = (
        action.symbols[name] for name in ("a", "b", "X", "a0bar", "U", "V", "W", "alpha", "beta", "Mbar")
    )
    s = b**2 / (a**2 * X**2)
    expected_z = sp.Rational(3, 2) / a0**2 * (2 * U**2 - (3 + s) * U * V + (1 + s) * V**2 + (s - 1) * V * W)
    expected_l = (
        -3 * beta * a**3 * U**2 - 3 * alpha * b**3 * V**2 / X + a0**2 * sp.sqrt(X) * (a * b) ** sp.Rational(3, 2) * M
    )
    return {"background_Z": zero(action.Z0 - expected_z), "background_L": zero(action.lagrangian0 - expected_l)}


def check_vector(module, action):
    checks = check_background(action)
    a, b, X, Mp, alpha, beta = (action.symbols[name] for name in ("a", "b", "X", "Mprime", "alpha", "beta"))
    _P, S, _Q, T = action.fields
    Pd, Sd, Qd, Td = action.time_derivatives
    Pz, Sz, Qz, Tz = action.spatial_derivatives
    s = b**2 / (a**2 * X**2)
    t = b**3 / (a**3 * X)
    m = sp.sqrt(X) * (b / a) ** sp.Rational(3, 2) * Mp
    checks["linear_vector_Z_zero"] = zero(action.Z1)
    checks["linear_vector_L_zero"] = zero(action.lagrangian1)
    principal = module.principal_part(action).subs({alpha: 1, beta: 1})
    # At frozen principal coefficients, Sd*Pz and Pd*Sz differ by a boundary.
    for shift_dot, shift_z in ((Sd, Sz), (Td, Tz)):
        for field_dot, field_z in ((Pd, Pz), (Qd, Qz)):
            coefficient = sp.diff(principal, shift_dot, field_z)
            principal = sp.expand(principal - coefficient * shift_dot * field_z + coefficient * field_dot * shift_z)
    u, v, d = Pd - Sz, Qd - Tz, Pd - Qd
    expected = a**3 / 4 * ((1 - m) * u**2 + 2 * m * u * v + (t - m) * v**2 - 2 * m * (1 - s) * d * v)
    checks["vector_principal_action"] = zero(principal - expected)
    checks["no_full_shift_velocity_hessian"] = all(
        zero(sp.diff(action.lagrangian2, shift_dot, derivative))
        for shift_dot in (Sd, Td)
        for derivative in action.time_derivatives
    )
    checks["no_shift_symplectic_pair"] = zero(sp.diff(action.lagrangian2, Sd, T) - sp.diff(action.lagrangian2, Td, S))
    # Derive the principal shift elimination without assuming the answer.
    ms, ss, ts, us, vs, ds = sp.symbols("m s t u v d", real=True)
    lag = (1 - ms) * us**2 + 2 * ms * us * vs + (ts - ms) * vs**2 - 2 * ms * (1 - ss) * ds * vs
    constraints = sp.solve([sp.diff(lag, us), sp.diff(lag, vs)], (us, vs))
    reduced = sp.factor(lag.subs(constraints))
    coefficient = -(ms**2) * (1 - ss) ** 2 * (1 - ms) / (ts - ms * (1 + ts))
    checks["reduced_vector_kinetic_coefficient"] = zero(reduced - coefficient * ds**2)
    checks["no_principal_vector_spatial_gradient"] = not principal.has(Pz, Qz)
    return checks, {
        "reduced_vector_coefficient": str(coefficient),
        "shift_solution": {str(k): str(v) for k, v in constraints.items()},
    }


def check_scalar_homogeneous(action):
    s = action.symbols
    a, b, X, a0, U, V, W = (s[name] for name in ("a", "b", "X", "a0bar", "U", "V", "W"))
    phi, S, zeta, ell, phihat, T, zetahat, ellhat = action.fields
    phid, Sd, zetad, elld, phihatd, Td, zetahatd, ellhatd = action.time_derivatives
    epsilon, du, dv, dn, dw, dU, dV, dN, dW = sp.symbols("epsilon du dv dn dw dU dV dN dW", real=True)
    substitutions = dict.fromkeys(action.spatial_derivatives, 0)
    substitutions.update({phi: dn, phihat: dn + dw, zeta: du, ell: du, zetahat: dv, ellhat: dv, S: 0, T: 0})
    substitutions.update({phid: dN, phihatd: dN + dW, zetad: dU, elld: dU, zetahatd: dV, ellhatd: dV, Sd: 0, Td: 0})
    snew = b**2 / (a**2 * X**2) * sp.exp(2 * epsilon * (dv - du - dw))
    unew, vnew, wnew = U + epsilon * dU, V + epsilon * dV, W + epsilon * dW
    znew = (
        sp.Rational(3, 2)
        / a0**2
        * sp.exp(-2 * epsilon * dn)
        * (2 * unew**2 - (3 + snew) * unew * vnew + (1 + snew) * vnew**2 + (snew - 1) * vnew * wnew)
    )
    z1 = sp.diff(znew, epsilon).subs(epsilon, 0)
    z2 = sp.diff(znew, epsilon, 2).subs(epsilon, 0) / 2
    kinetic = -3 * s["beta"] * a**3 * sp.exp(epsilon * (3 * du - dn)) * unew**2
    kinetic -= 3 * s["alpha"] * b**3 / X * sp.exp(epsilon * (3 * dv - dw - dn)) * vnew**2
    interaction = (
        a0**2
        * sp.sqrt(X)
        * (a * b) ** sp.Rational(3, 2)
        * sp.exp(epsilon * (dn + sp.Rational(3, 2) * (du + dv) + dw / 2))
    )
    interaction *= s["Mbar"] + epsilon * s["Mprime"] * z1 + epsilon**2 * (s["Mprime"] * z2 + s["Msecond"] * z1**2 / 2)
    l1 = sp.diff(kinetic + interaction, epsilon).subs(epsilon, 0)
    l2 = sp.diff(kinetic + interaction, epsilon, 2).subs(epsilon, 0) / 2
    return {
        "homogeneous_Z1": zero(action.Z1.subs(substitutions) - z1),
        "homogeneous_Z2": zero(action.Z2.subs(substitutions) - z2),
        "homogeneous_L1": zero(action.lagrangian1.subs(substitutions) - l1),
        "homogeneous_L2": zero(action.lagrangian2.subs(substitutions) - l2),
        "homogeneous_common_lapse_derivative_cancels": zero(sp.diff(action.lagrangian2.subs(substitutions), dN)),
    }


def check_scalar_symbol(action, directory):
    rates, gradients = action.time_derivatives, action.spatial_derivatives
    K = sp.hessian(action.lagrangian2, rates)
    B = sp.Matrix([[sp.diff(action.lagrangian2, rate, gradient) for gradient in gradients] for rate in rates])
    D = sp.hessian(action.lagrangian2, gradients)
    omega, k = sp.symbols("omega k", real=True)
    symbol = -(omega**2) * K + omega * k * (B + B.T) - k**2 * D
    a, b, X = (action.symbols[name] for name in ("a", "b", "X"))
    time_gauge = sp.Matrix([omega, k / a**2, 0, 0, omega, k * X**2 / b**2, 0, 0])
    spatial_gauge = sp.Matrix([0, omega, 0, -k, 0, omega, 0, -k])
    checks = {
        "principal_common_time_gauge": all(zero(item) for item in symbol * time_gauge),
        "principal_common_spatial_gauge": all(zero(item) for item in symbol * spatial_gauge),
        "scalar_time_hessian_symmetric": K == K.T,
    }
    payload = {
        "field_order": list(map(str, action.fields)),
        "K_time": [[str(value) for value in row] for row in K.tolist()],
        "B_time_space": [[str(value) for value in row] for row in B.tolist()],
        "D_space": [[str(value) for value in row] for row in D.tolist()],
        "symbol": "-omega**2*K + omega*k*(B+B.T) - k**2*D",
        "scope": "Unreduced principal scalar symbol, no physical stability conclusion from unreduced eigenvalues.",
    }
    with (directory / "scalar_principal_matrices.json").open("x") as stream:
        json.dump(payload, stream, indent=2)
        stream.write("\n")
    return checks


def main():
    module = load_generator(sys.argv[1])
    directory = Path(sys.argv[2])
    directory.mkdir()
    report = {
        "sympy_version": sp.__version__,
        "scope": "Gravitational scalar/vector quadratic action, alpha/beta arbitrary, f=1",
    }
    for sector in ("vector", "scalar"):
        started = time.monotonic()
        print(f"Generating {sector} action", flush=True)
        action = module.quadratic_action(sector)
        store_action(action, directory)
        print(f"Checking {sector} action", flush=True)
        if sector == "vector":
            checks, details = check_vector(module, action)
        else:
            checks = check_background(action)
            checks.update(check_scalar_homogeneous(action))
            checks.update(check_scalar_symbol(action, directory))
            details = {}
        report[sector] = {
            "checks": checks,
            "all_pass": all(checks.values()),
            "wall_seconds": time.monotonic() - started,
            **details,
        }
        print(json.dumps(report[sector]), flush=True)
        with (directory / f"{sector}_checks.json").open("x") as stream:
            json.dump(report[sector], stream, indent=2)
            stream.write("\n")
    report["all_pass"] = report["vector"]["all_pass"] and report["scalar"]["all_pass"]
    with (directory / "summary.json").open("x") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    if not report["all_pass"]:
        raise AssertionError(report)


if __name__ == "__main__":
    main()
