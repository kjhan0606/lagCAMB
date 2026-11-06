"""Scheduled, standalone local scalar constraint audit; not a Boltzmann solver.

Run only on an allocated scientific-computation node. Imports the generator
without importing CAMB's compiled extension. Densities use 8*pi*G=1. The output
is a new JSON file (exclusive creation); no existing artifact is overwritten.
Only two local on-shell jets are tested, never the whole BIMOND action.
"""

import argparse
import importlib.util
import json
import sys
import traceback
from pathlib import Path

import mpmath as mp
import sympy as sp


def generator(path):
    spec = importlib.util.spec_from_file_location("standalone_bimond_symbolic", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--generator", type=Path, default=Path(__file__).resolve().parents[2] / "camb" / "bimond_symbolic.py"
    )
    parser.add_argument("--dps", type=int, default=80)
    args = parser.parse_args()
    if args.dps < 70:
        parser.error("Use at least 70 digits for this cancellation-sensitive audit")
    # Reserve before expensive work; exceptions also leave an explicit report.
    with args.output.open("x") as stream:
        report = {"scope": "two local jets, not a whole-action exclusion", "dps": args.dps}
        try:
            mp.mp.dps = args.dps
            report.update(audit(args.dps, args.generator))
        except Exception:
            report["fatal_error"] = traceback.format_exc()
            json.dump(report, stream, indent=2)
            stream.write("\n")
            raise
        json.dump(report, stream, indent=2)
        stream.write("\n")
    if any(point.get("error") for point in report["points"]):
        raise SystemExit(1)


def audit(dps, generator_path):
    action = generator(generator_path).quadratic_action("scalar", include_matter=True)
    sy = dict(action.symbols)
    required = ("dust_density", "radiation_rate", "radiation_kappa")
    if any(name not in sy for name in required):
        raise RuntimeError("Generator does not expose the agreed matter API")
    a, b, X, a0, U, V, W = (sy[name] for name in ("a", "b", "X", "a0bar", "U", "V", "W"))
    M, Mp, Mpp = (sy[name] for name in ("Mbar", "Mprime", "Msecond"))
    rho, rate, kap = (sy[name] for name in required)
    M3 = sp.Symbol("Mthird", real=True)
    acc = sp.symbols("U_acc V_acc W_acc", real=True)
    k = sp.Symbol("k", positive=True)
    z = action.Z0
    positions, velocities = (a, b, X), (U, V, W)
    constant = {sy["alpha"]: 1, sy["beta"]: 1, a0: 1}
    volume = sp.sqrt(X) * (a * b) ** sp.Rational(3, 2)
    kin = -3 * a**3 * U**2 - 3 * b**3 * V**2 / X
    Lg = kin + a0**2 * volume * M
    generated_gravity0 = action.lagrangian0 - a**3 * kap * rate**4 / 4
    if sp.simplify((generated_gravity0 - Lg).subs(constant)) != 0:
        raise RuntimeError("Generator homogeneous gravitational action differs from selected minisuperspace action")
    energy = kin + a0**2 * volume * (2 * z * Mp - M) + a**3 * (rho + 3 * kap * rate**4 / 4)

    def derivative(expr, variable):
        """Partial derivative on the specified free-function jet."""
        return sp.diff(expr, variable) + sp.diff(z, variable) * (
            sp.diff(expr, M) * Mp + sp.diff(expr, Mp) * Mpp + sp.diff(expr, Mpp) * M3
        )

    def qderivative(expr, index):
        return positions[index] * derivative(expr, positions[index])

    momenta = [derivative(Lg, vel) for vel in velocities]
    H = sp.Matrix([[derivative(p, vel) for vel in velocities] for p in momenta])
    rhs = sp.Matrix(
        [
            qderivative(Lg, i)
            - sum(qderivative(momenta[i], j) * velocities[j] for j in range(3))
            + (a**3 * 3 * kap * rate**4 / 4 if i == 0 else 0)
            for i in range(3)
        ]
    )

    def dt(expr):
        return (
            sum(
                derivative(expr, var) * speed for var, speed in zip(positions + velocities, (a * U, b * V, X * W) + acc)
            )
            + sp.diff(expr, rho) * (-3 * U * rho)
            + sp.diff(expr, rate) * (-U * rate)
        )

    field = {str(f): f for f in action.fields}
    dots = dict(zip(action.fields, action.time_derivatives))
    spatial = dict(zip(action.fields, action.spatial_derivatives))
    expected = ("phi", "S", "zeta", "ell", "phihat", "T", "zetahat", "ellhat", "tau", "dlambda", "chi")
    if tuple(map(str, action.fields)) != expected:
        raise RuntimeError("Unexpected matter field order")
    tadpoles = {
        str(f): sp.diff(action.lagrangian1, f) - dt(sp.diff(action.lagrangian1, dots[f])) for f in action.fields
    }
    # Spatial derivative terms in L1 integrate to boundaries: all their
    # coefficients are homogeneous background functions.
    gauge = {}
    for name in ("tau", "phi", "dlambda"):
        f = field[name]
        gauge.update({f: 0, dots[f]: 0, spatial[f]: 0})
    gauge.update(
        {
            field["ell"]: field["zeta"],
            dots[field["ell"]]: dots[field["zeta"]],
            spatial[field["ell"]]: spatial[field["zeta"]],
        }
    )
    # Verify multiplier variation really imposes phi=0 after dust-time gauge.
    dust_time = {field["tau"]: 0, dots[field["tau"]]: 0, spatial[field["tau"]]: 0}
    multiplier_equation = sp.expand(sp.diff(action.lagrangian2, field["dlambda"]).subs(dust_time))
    multiplier_check = sp.simplify(multiplier_equation + a**3 * field["phi"])
    if multiplier_check != 0:
        raise RuntimeError(f"Unexpected dust multiplier constraint: {multiplier_equation}")
    L = action.lagrangian2.subs(gauge, simultaneous=True).subs(constant)
    names = ("zeta", "S", "phihat", "T", "zetahat", "ellhat", "chi")
    q = [field[name] for name in names]
    dq = [dots[f] for f in q]
    cosine, sine = sp.symbols("cosine sine", real=True)
    fourier = {}
    for name, f, fd in zip(names, q, dq):
        odd = name in ("S", "T")
        fourier.update(
            {
                f: f * (sine if odd else cosine),
                fd: fd * (sine if odd else cosine),
                spatial[f]: k * f * cosine if odd else -k * f * sine,
            }
        )
    L = sp.expand(L.subs(fourier, simultaneous=True))
    L = sp.expand((L.subs({cosine: 1, sine: 0}) + L.subs({cosine: 0, sine: 1})) / 2)
    K = sp.hessian(L, dq)
    A = sp.Matrix([[sp.diff(L, vi, qi) for qi in q] for vi in dq])
    D = sp.hessian(L, q)
    aux, dyn = [1, 3], [0, 2, 4, 5, 6]
    Aaa = A.extract(aux, aux)
    # Integration by parts of both auxiliary-velocity blocks, with complete
    # background derivatives. These are exact finite-k operations.
    Ay = D.extract(aux, aux) - Aaa.applyfunc(dt)
    C = A.extract(dyn, aux) - A.extract(aux, dyn).T
    Kd = K.extract(dyn, dyn)
    gates = [K[i, j] for i in aux for j in range(7)] + [K[j, i] for i in aux for j in range(7)]
    gates += list(Aaa - Aaa.T)
    zero_phi = {field["phi"]: 0, dots[field["phi"]]: 0, spatial[field["phi"]]: 0}
    # A full phi EL reconstructs dlambda, not an additional metric constraint.
    # Its nonzero algebraic multiplier coefficient suffices to establish this.
    reconstruction_coefficient = sp.diff(action.lagrangian2, field["phi"], field["dlambda"])
    assert reconstruction_coefficient.subs(zero_phi) != 0

    def number(expr, values):
        result = expr.subs(values, simultaneous=True).evalf(dps)
        if result.free_symbols:
            raise RuntimeError(f"Unsubstituted coefficient symbols: {result.free_symbols}")
        return mp.mpf(str(result))

    def matrix(expr, values):
        return mp.matrix([[number(expr[i, j], values) for j in range(expr.cols)] for i in range(expr.rows)])

    def symfloat(value):
        return sp.Float(mp.nstr(value, dps), dps)

    def txt(value):
        return mp.nstr(value, dps)

    def rows(value):
        return [[txt(value[i, j]) for j in range(value.cols)] for i in range(value.rows)]

    tolerance = mp.mpf("1e-50")

    def inverse_checked(value, label):
        singular_values = mp.svd(value, compute_uv=False)
        largest, smallest = max(singular_values), min(singular_values)
        if largest == 0 or smallest <= tolerance * largest:
            raise RuntimeError(f"{label} singular or unresolved at audit precision; no regularization")
        return value**-1, {
            "largest_singular_value": txt(largest),
            "smallest_singular_value": txt(smallest),
            "condition_number": txt(largest / smallest),
        }

    points = []
    for sign in (1, -1):
        point = {"Z_sign": sign, "k_results": []}
        points.append(point)
        try:
            y = 4 * (sp.sqrt(2) - 1)
            target_z = sign * y**2
            values = dict(constant)
            values.update(
                {
                    a: 1,
                    b: 1 / sp.sqrt(2),
                    X: sp.sqrt(2),
                    U: 1,
                    V: sp.Rational(1, 2),
                    rate: 1,
                    kap: sp.Rational(2, 75),
                    M: sign * (4 * y - 16 * sp.log(1 + y / 4)),
                    Mp: 2 / (4 + y),
                    Mpp: -sign / (y * (4 + y) ** 2),
                    M3: 1 / (2 * y**3 * (4 + y) ** 2) + 1 / (y**2 * (4 + y) ** 3),
                }
            )
            z_partial = z.subs(values)
            w_solutions = sp.solve(z_partial - target_z, W)
            if len(w_solutions) != 1:
                raise RuntimeError("Z constraint does not uniquely determine W")
            values[W] = w_solutions[0]
            dust_solutions = sp.solve(energy.subs(values), rho)
            if len(dust_solutions) != 1:
                raise RuntimeError("Energy constraint does not uniquely determine dust density")
            values[rho] = dust_solutions[0]
            if number(rho, values) <= 0:
                raise RuntimeError("Nonpositive dust density")
            hm = matrix(H, values)
            hi, hinfo = inverse_checked(hm, "homogeneous velocity Hessian")
            accelerations = hi * matrix(rhs, values)
            values.update({var: symfloat(accelerations[i]) for i, var in enumerate(acc)})
            point["background"] = {str(var): txt(number(var, values)) for var in values}
            point["background"].update(
                {
                    "Z": txt(number(z, values)),
                    "radiation_density": "0.02",
                    "energy_residual": txt(number(energy, values)),
                    "Z_dot": txt(number(dt(z), values)),
                }
            )
            point["homogeneous_hessian"] = hinfo
            point["tadpole_residuals"] = {
                name: txt(number(expr.subs(constant), values)) for name, expr in tadpoles.items()
            }
            # Normalize each tadpole by its separate terms to distinguish true
            # cancellation from a dimensionful absolute tolerance.
            normalized = []
            for f in action.fields:
                left = number(sp.diff(action.lagrangian1, f).subs(constant), values)
                right = number(dt(sp.diff(action.lagrangian1, dots[f])).subs(constant), values)
                normalized.append(abs(left - right) / max(mp.mpf(1), abs(left), abs(right)))
            point["maximum_relative_tadpole_residual"] = txt(max(normalized))
            if max(normalized) > tolerance:
                raise RuntimeError("Matter-completed tadpoles do not vanish on this local jet")
            for wave in (10, 100, 1000):
                kv = dict(values)
                kv[k] = wave
                gate_residual = max(abs(number(expr, kv)) for expr in gates)
                if gate_residual > tolerance:
                    raise RuntimeError("Shift Hessian/symplectic auxiliary gate failed")
                am = matrix(Ay, kv)
                if mp.norm(am - am.T) > tolerance * max(1, mp.norm(am)):
                    raise RuntimeError("Integrated auxiliary algebraic block is not symmetric")
                ai, ainfo = inverse_checked(am, "finite-k shift algebraic block")
                cm = matrix(C, kv)
                physical = matrix(Kd, kv) - cm * ai * cm.T
                if mp.norm(physical - physical.T) > tolerance * max(1, mp.norm(physical)):
                    raise RuntimeError("Physical kinetic matrix is not symmetric")
                physical = (physical + physical.T) / 2
                eigenvalues, eigenvectors = mp.eigsy(physical)
                scale = max(abs(e) for e in eigenvalues)
                threshold = tolerance * max(mp.mpf(1), scale)
                rank = sum(abs(e) > threshold for e in eigenvalues)
                negatives = [i for i, e in enumerate(eigenvalues) if e < -threshold]
                record = {
                    "k": wave,
                    "auxiliary_gate_maximum": txt(gate_residual),
                    "auxiliary_block": ainfo,
                    "Kphys": rows(physical),
                    "eigenvalues": [txt(e) for e in eigenvalues],
                    "eigenvalues_relative_to_maximum": [txt(e / scale) if scale else "0" for e in eigenvalues],
                    "numerical_rank": rank,
                    "negative_eigenvalues": len(negatives),
                    "condition_number": txt(scale / min(abs(e) for e in eigenvalues)) if rank == 5 else None,
                }
                if rank < 5:
                    record["verdict"] = "Remaining null/first-order sector: full Dirac reduction required"
                elif negatives:
                    record["verdict"] = "Nondegenerate reduced scalar kinetic ghost at this local point only"
                    index = negatives[0]
                    certificate = eigenvectors[:, index]
                    record["negative_direction"] = [txt(certificate[i]) for i in range(5)]
                    record["negative_quadratic_form"] = txt((certificate.T * physical * certificate)[0])
                else:
                    record["verdict"] = "Local scalar kinetic test passes; gradient and global stability untested"
                point["k_results"].append(record)
        except Exception:  # noqa: BLE001 -- Record each point's failure and still evaluate the other sign.
            point["error"] = traceback.format_exc()
    return {
        "dynamic_field_order": [names[i] for i in dyn],
        "points": points,
        "limitations": [
            "No whole-action exclusion from these two points",
            "No characteristic-speed audit",
            "Dust and P(Y) radiation are not a photon/neutrino Boltzmann hierarchy",
            "No cosmological-history viability claim",
            "No inverse was regularized",
        ],
    }


if __name__ == "__main__":
    main()
