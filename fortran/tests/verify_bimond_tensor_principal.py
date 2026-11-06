"""Independently contract the connections for the principal TT quadratic action.

Run on a compute node. The diagonal polarization is exp(2h), exp(-2h), 1.
The result concerns the specified alpha=-1, beta=1, f=1 action only.
"""

import json
import sys
from pathlib import Path

import sympy as sp

a, b, X, a0 = sp.symbols("a b X a0", positive=True)
U, V, W = sp.symbols("U V W", real=True)
p, phat, dz, dzhat = sp.symbols("p phat dz dzhat", real=True)
Mp = sp.symbols("Mp", positive=True)
perturbations = (p, phat, dz, dzhat)
zero = dict.fromkeys(perturbations, 0)


def connection(scale, lapse, rate, lapse_rate, tensor_rate, tensor_gradient):
    inverse = [-1 / lapse**2, 1 / scale**2, 1 / scale**2, 1 / scale**2]
    derivatives = sp.MutableDenseNDimArray.zeros(4, 4, 4)
    derivatives[0, 0, 0] = -2 * lapse**2 * lapse_rate
    for i, polarization in ((1, 1), (2, -1), (3, 0)):
        derivatives[0, i, i] = 2 * scale**2 * (rate + polarization * tensor_rate)
        derivatives[3, i, i] = 2 * scale**2 * polarization * tensor_gradient
    gamma = sp.MutableDenseNDimArray.zeros(4, 4, 4)
    for alpha in range(4):
        for mu in range(4):
            for nu in range(4):
                gamma[alpha, mu, nu] = (
                    inverse[alpha]
                    * (derivatives[mu, alpha, nu] + derivatives[nu, alpha, mu] - derivatives[alpha, mu, nu])
                    / 2
                )
    return gamma


ordinary = connection(a, sp.S.One, U, sp.S.Zero, p, dz)
twin = connection(b, X, V, W, phat, dzhat)
C = ordinary - twin
inverse = [-1, 1 / a**2, 1 / a**2, 1 / a**2]
upsilon = sp.S.Zero
for mu in range(4):
    for alpha in range(4):
        for lam in range(4):
            upsilon += inverse[mu] * (C[alpha, mu, lam] * C[lam, mu, alpha] - C[alpha, mu, mu] * C[lam, lam, alpha])
Z = sp.expand(-upsilon / (2 * a0**2))
Z0 = Z.subs(zero)
checks = {}
checks["no_linear_tensor_Z"] = all(sp.simplify(sp.diff(Z, item).subs(zero)) == 0 for item in perturbations)
s = b**2 / (a**2 * X**2)
expected_delta = (-(p**2) + 2 * s * p * phat + (1 - 2 * s) * phat**2 + (dz - dzhat) ** 2 / a**2) / a0**2
checks["connection_contraction"] = sp.simplify(Z - Z0 - expected_delta) == 0

# Eight pi G times the quadratic action in exponential polarization variables.
# The change h_linear=2*h_exponential cancels from both normalized matrices.
lagrangian = (
    a**3 * (p**2 - dz**2 / a**2)
    - b**3 / X * (phat**2 - X**2 * dzhat**2 / b**2)
    + a0**2 * sp.sqrt(X) * (a * b) ** sp.Rational(3, 2) * Mp * (Z - Z0)
)
kinetic = sp.hessian(lagrangian, (p, phat)) / (2 * a**3)
gradient = -sp.hessian(lagrangian, (dz, dzhat)) / (2 * a)
t = b**3 / (a**3 * X)
m = sp.sqrt(X) * (b / a) ** sp.Rational(3, 2) * Mp
expected_kinetic = sp.Matrix([[1 - m, m * s], [m * s, -t - m * (2 * s - 1)]])
expected_gradient = sp.Matrix([[1 - m, m], [m, -t / s - m]])
checks["tensor_kinetic_matrix"] = all(sp.simplify(item) == 0 for item in kinetic - expected_kinetic)
checks["tensor_gradient_matrix"] = all(sp.simplify(item) == 0 for item in gradient - expected_gradient)

ms, ss, ts = sp.symbols("m s t", positive=True)
K = sp.Matrix([[1 - ms, ms * ss], [ms * ss, -ts - ms * (2 * ss - 1)]])
G = sp.Matrix([[1 - ms, ms], [ms, -ts / ss - ms]])
checks["kinetic_determinant"] = sp.simplify(K.det() - (-ts + ms * (ts + 1 - 2 * ss) - ms**2 * (ss - 1) ** 2)) == 0
checks["gradient_determinant"] = sp.simplify(G.det() + ts / ss * (1 - ms) + ms) == 0
checks["formal_equal_ratio_determinant"] = sp.simplify(K.subs({ss: 1, ts: 1}).det() + 1) == 0
checks["high_acceleration_kinetic"] = K.subs(ms, 0) == sp.diag(1, -ts)
checks["high_acceleration_gradient"] = G.subs(ms, 0) == sp.diag(1, -ts / ss)

report = {
    "checks": checks,
    "all_pass": all(checks.values()),
    "scope": "Principal TT sector of alpha=-1,beta=1,f=1 with M'(Z)>0 at Z!=0",
    "implication": "K positive definite requires m<1, implying det(G)<0 for s,t,m>0. Otherwise K is negative or degenerate in at least one tensor direction.",
    "not_claimed": "No conclusion about other BIMOND actions or a constrained theory on a degenerate tensor surface.",
    "sympy_version": sp.__version__,
}
with Path(sys.argv[1]).open("x") as stream:
    json.dump(report, stream, indent=2)
    stream.write("\n")
if not report["all_pass"]:
    raise AssertionError(report)
