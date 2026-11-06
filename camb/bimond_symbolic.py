"""Symbolic scalar/vector gravitational actions for the f=1 BIMOND model.

This module generates the complete coefficients through second order in a
formal perturbation parameter about arbitrary flat FLRW. Both Einstein-Hilbert
terms use their first-derivative connection-squared form, modulo a boundary
term, in the original Ricci-sign convention. No perturbative field equation,
gauge choice, integration by parts, or background equation is imposed.

The default result is a gravitational action. Optional ordinary irrotational
dust and a radiation scalar provide a perfect-fluid kinetic/constraint audit;
they do not model photons, neutrinos, free streaming, or recombination.
Constraint reduction is still needed for a physical scalar stability verdict.
No Boltzmann equations or transfer functions are supplied. Generation and
simplification happen only when functions are called, never at import.
"""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

import sympy as sp

Jet = tuple[sp.Expr, sp.Expr, sp.Expr]
SparseMetric = dict[tuple[int, int], Jet]
SparseConnection = dict[tuple[int, int, int], Jet]
_ZERO: Jet = (sp.S.Zero, sp.S.Zero, sp.S.Zero)


@dataclass(frozen=True)
class QuadraticAction:
    """Coefficients of 8*pi*G times the Lagrangian density.

    ``lagrangian2`` is the coefficient of epsilon**2, with no extra factorial.
    The interaction is a0bar**2*(-g)**(1/4)*(-ghat)**(1/4)*M(Z), where
    Z=-g**(mu,nu)*Upsilon(C)_(mu,nu)/(2*a0bar**2). The three interaction
    coefficients Mbar, Mprime, Msecond are independent background symbols.

    Fields depend on t and z; indices (0,1,2,3) denote (t,x,y,z). Scalar fields
    are (phi,S,zeta,ell,phihat,T,zetahat,ellhat), with exponential ADM lapses
    and diagonal spatial metrics. Vector fields are (P,S,Q,T), one transverse
    x-polarization; P=partial_z(F_x) and Q=partial_z(Fhat_x) when Fourier
    reducing. S,T are contravariant shifts in both sectors. Time and spatial
    derivatives are independent symbols in the corresponding field order.

    With include_matter=True, three scalar fields follow the original eight:
    tau (dust clock), dlambda (dust density perturbation), and chi (radiation
    perturbation). Matter couples only to the ordinary metric.
    dust_density, radiation_rate, and radiation_kappa are positive background
    symbols; the density and kappa include the action's 8*pi*G normalization.
    """

    sector: Literal["scalar", "vector"]
    symbols: Mapping[str, sp.Symbol]
    fields: tuple[sp.Symbol, ...]
    time_derivatives: tuple[sp.Symbol, ...]
    spatial_derivatives: tuple[sp.Symbol, ...]
    lagrangian0: sp.Expr
    lagrangian1: sp.Expr
    lagrangian2: sp.Expr
    Z0: sp.Expr
    Z1: sp.Expr
    Z2: sp.Expr
    include_matter: bool = False


def _nonzero(jet: Jet) -> bool:
    return any(coefficient != 0 for coefficient in jet)


def _sum(*jets: Jet) -> Jet:
    return tuple(sp.expand(sp.Add(*(jet[order] for jet in jets))) for order in range(3))


def _scale(jet: Jet, factor: sp.Expr) -> Jet:
    return tuple(sp.expand(factor * coefficient) for coefficient in jet)


def _product(left: Jet, right: Jet) -> Jet:
    """Truncated convolution; all inputs/outputs are coefficients, not derivatives."""
    return tuple(
        sp.expand(sp.Add(*(left[index] * right[order - index] for index in range(order + 1)))) for order in range(3)
    )


def _symmetric_entry(metric: SparseMetric, row: int, column: int, jet: Jet) -> None:
    if _nonzero(jet):
        metric[row, column] = jet
        metric[column, row] = jet


def _scalar_metric(scale: sp.Symbol, lapse: sp.Expr, fields: tuple[sp.Symbol, ...]) -> SparseMetric:
    phi, shift, zeta, ell = fields
    metric: SparseMetric = {}
    _symmetric_entry(
        metric,
        0,
        0,
        (-(lapse**2), -2 * lapse**2 * phi, -2 * lapse**2 * phi**2 + scale**2 * shift**2),
    )
    _symmetric_entry(metric, 0, 3, (sp.S.Zero, scale**2 * shift, 2 * scale**2 * ell * shift))
    for direction, perturbation in ((1, zeta), (2, zeta), (3, ell)):
        _symmetric_entry(
            metric, direction, direction, (scale**2, 2 * scale**2 * perturbation, 2 * scale**2 * perturbation**2)
        )
    return metric


def _vector_metric(scale: sp.Symbol, lapse: sp.Expr, fields: tuple[sp.Symbol, ...]) -> SparseMetric:
    perturbation, shift = fields
    metric: SparseMetric = {}
    _symmetric_entry(metric, 0, 0, (-(lapse**2), sp.S.Zero, scale**2 * shift**2))
    _symmetric_entry(metric, 0, 1, (sp.S.Zero, scale**2 * shift, sp.S.Zero))
    _symmetric_entry(metric, 0, 3, (sp.S.Zero, sp.S.Zero, scale**2 * perturbation * shift))
    for direction in (1, 2, 3):
        _symmetric_entry(metric, direction, direction, (scale**2, sp.S.Zero, sp.S.Zero))
    _symmetric_entry(metric, 1, 3, (sp.S.Zero, scale**2 * perturbation, sp.S.Zero))
    return metric


def _inverse(metric: SparseMetric) -> SparseMetric:
    """Full inverse jets: D, -D*g1*D, D*g1*D*g1*D-D*g2*D.

    The FLRW zeroth-order matrix is diagonal. No exact perturbed matrix inverse
    or generic expression series is formed; off-diagonal second-order inverse
    terms are retained by the recurrence.
    """
    diagonal = tuple(1 / metric[index, index][0] for index in range(4))
    inverse: SparseMetric = {}
    for row in range(4):
        for column in range(row, 4):
            entry = metric.get((row, column), _ZERO)
            quadratic = sp.Add(
                *(
                    metric.get((row, middle), _ZERO)[1] * diagonal[middle] * metric.get((middle, column), _ZERO)[1]
                    for middle in range(4)
                )
            )
            jet = (
                diagonal[row] if row == column else sp.S.Zero,
                sp.expand(-diagonal[row] * entry[1] * diagonal[column]),
                sp.expand(diagonal[row] * (quadratic - entry[2]) * diagonal[column]),
            )
            _symmetric_entry(inverse, row, column, jet)
    return inverse


def _metric_derivatives(
    metric: SparseMetric,
    symbols: Mapping[str, sp.Symbol],
    fields: tuple[sp.Symbol, ...],
    time_derivatives: tuple[sp.Symbol, ...],
    spatial_derivatives: tuple[sp.Symbol, ...],
) -> dict[int, SparseMetric]:
    a, b, X = (symbols[name] for name in ("a", "b", "X"))
    U, V, W = (symbols[name] for name in ("U", "V", "W"))
    temporal = ((a, a * U), (b, b * V), (X, X * W), *zip(fields, time_derivatives))
    spatial = tuple(zip(fields, spatial_derivatives))
    result: dict[int, SparseMetric] = {0: {}, 3: {}}
    for coordinate, chain_rule in ((0, temporal), (3, spatial)):
        for indices, jet in metric.items():
            differentiated = tuple(
                sp.expand(
                    sp.Add(
                        *(
                            sp.diff(coefficient, variable) * derivative
                            for variable, derivative in chain_rule
                            if variable in coefficient.free_symbols
                        )
                    )
                )
                for coefficient in jet
            )
            if _nonzero(differentiated):
                result[coordinate][indices] = differentiated
    return result


def _connection(inverse: SparseMetric, derivatives: dict[int, SparseMetric]) -> SparseConnection:
    connection: SparseConnection = {}
    for upper in range(4):
        inverse_row = [(column, jet) for (row, column), jet in inverse.items() if row == upper]
        for mu in range(4):
            for nu in range(mu, 4):
                pieces = []
                for contracted, inverse_jet in inverse_row:
                    numerator = _sum(
                        derivatives.get(mu, {}).get((contracted, nu), _ZERO),
                        derivatives.get(nu, {}).get((contracted, mu), _ZERO),
                        _scale(derivatives.get(contracted, {}).get((mu, nu), _ZERO), -sp.S.One),
                    )
                    if _nonzero(numerator):
                        pieces.append(_product(inverse_jet, numerator))
                jet = _scale(_sum(*pieces), sp.Rational(1, 2))
                if _nonzero(jet):
                    connection[upper, mu, nu] = jet
                    connection[upper, nu, mu] = jet
    return connection


def _contract_upsilon(inverse: SparseMetric, connection: SparseConnection) -> Jet:
    """g^(mu,nu) [A^a_(mu,l) A^l_(nu,a) - A^a_(mu,nu) A^l_(l,a)]."""
    traces = {upper: _sum(*(connection.get((index, index, upper), _ZERO) for index in range(4))) for upper in range(4)}
    pieces = []
    for (mu, nu), inverse_jet in inverse.items():
        products = []
        for upper in range(4):
            for middle in range(4):
                left = connection.get((upper, mu, middle))
                right = connection.get((middle, nu, upper))
                if left is not None and right is not None:
                    products.append(_product(left, right))
            left = connection.get((upper, mu, nu))
            if left is not None and _nonzero(traces[upper]):
                products.append(_scale(_product(left, traces[upper]), -sp.S.One))
        if products:
            pieces.append(_product(inverse_jet, _sum(*products)))
    return _sum(*pieces)


def _contract_covector(inverse: SparseMetric, covector: dict[int, Jet]) -> Jet:
    """Full inverse-metric contraction of a covector with itself through order two."""
    return _sum(
        *(
            _product(inverse_jet, _product(covector[mu], covector[nu]))
            for (mu, nu), inverse_jet in inverse.items()
            if mu in covector and nu in covector
        )
    )


def _ordinary_matter(
    inverse: SparseMetric,
    volume: Jet,
    symbols: Mapping[str, sp.Symbol],
    fields: tuple[sp.Symbol, ...],
    time_derivatives: tuple[sp.Symbol, ...],
    spatial_derivatives: tuple[sp.Symbol, ...],
) -> Jet:
    """Irrotational dust multiplier action plus the healthy kappa*Y**2 scalar."""
    clock_gradient = {
        0: (sp.S.One, time_derivatives[8], sp.S.Zero),
        3: (sp.S.Zero, spatial_derivatives[8], sp.S.Zero),
    }
    dust_normalization = _sum(_contract_covector(inverse, clock_gradient), (sp.S.One, sp.S.Zero, sp.S.Zero))
    dust_multiplier = (symbols["dust_density"], fields[9], sp.S.Zero)
    dust = _scale(_product(volume, _product(dust_multiplier, dust_normalization)), -sp.Rational(1, 2))
    radiation_gradient = {
        0: (symbols["radiation_rate"], time_derivatives[10], sp.S.Zero),
        3: (sp.S.Zero, spatial_derivatives[10], sp.S.Zero),
    }
    kinetic_scalar = _scale(_contract_covector(inverse, radiation_gradient), -sp.Rational(1, 2))
    radiation = _scale(_product(volume, _product(kinetic_scalar, kinetic_scalar)), symbols["radiation_kappa"])
    return _sum(dust, radiation)


def quadratic_action(sector: Literal["scalar", "vector"], *, include_matter: bool = False) -> QuadraticAction:
    """Generate the complete order-two action, optionally with scalar matter.

    All background functions are arbitrary: adot=a*U, bdot=b*V, Xdot=X*W.
    a,b,X,a0bar are positive symbols; alpha,beta and U,V,W are real. Ordinary
    background lapse is one and the twin background lapse is X. Neither lapse
    perturbation is fixed in the scalar ansatz. f=1, with arbitrary alpha,beta.

    Only expansion of sparse degree-two jets is used. The returned coefficients
    retain field-only and single-derivative terms as well as principal terms.
    No background/tadpole equations or cancellation by integration by parts is
    assumed. In particular, Mbar and background-rate terms must not be dropped
    before constraint analysis. Call this potentially expensive generator only
    in the designated scientific-computation environment.

    include_matter=True is supported only for the scalar sector. It adds
    -sqrt(-g)*lambda*(g^(mu,nu)*d_mu(T)*d_nu(T)+1)/2, with
    T=t+epsilon*tau and lambda=dust_density+epsilon*dlambda, and the ordinary
    radiation action sqrt(-g)*radiation_kappa*Y**2, with
    Y=-g^(mu,nu)*d_mu(chi_total)*d_nu(chi_total)/2 and
    d_t(chi_total)=radiation_rate+epsilon*chi_dot. The background radiation
    density is 3*radiation_kappa* radiation_rate**4/4 and p=rho/3. On-shell
    background conservation gives radiation_rate_dot=-U*radiation_rate and
    dust_density_dot=-3*U*dust_density; neither relation is imposed here.
    Both matter couplings have the same 8*pi*G normalization as the returned
    action. This perfect-fluid audit does not implement a photon/neutrino
    hierarchy, free streaming, recombination, or twin matter.
    """
    if sector not in ("scalar", "vector"):
        raise ValueError("sector must be 'scalar' or 'vector'")
    if include_matter and sector != "scalar":
        raise ValueError("include_matter is supported only for the scalar sector")
    symbols = {name: sp.Symbol(name, positive=True) for name in ("a", "b", "X", "a0bar")}
    symbols.update(
        {
            name: sp.Symbol(name, real=True)
            for name in (
                "U",
                "V",
                "W",
                "alpha",
                "beta",
                "Mbar",
                "Mprime",
                "Msecond",
            )
        }
    )
    a, b, X, a0 = (symbols[name] for name in ("a", "b", "X", "a0bar"))
    names = (
        ("phi", "S", "zeta", "ell", "phihat", "T", "zetahat", "ellhat") if sector == "scalar" else ("P", "S", "Q", "T")
    )
    if include_matter:
        symbols.update(
            {
                name: sp.Symbol(name, positive=True)
                for name in (
                    "dust_density",
                    "radiation_rate",
                    "radiation_kappa",
                )
            }
        )
        names += ("tau", "dlambda", "chi")
    fields = tuple(sp.Symbol(name, real=True) for name in names)
    time_derivatives = tuple(sp.Symbol(f"{name}_dot", real=True) for name in names)
    spatial_derivatives = tuple(sp.Symbol(f"{name}_z", real=True) for name in names)
    if sector == "scalar":
        metric = _scalar_metric(a, sp.S.One, fields[:4])
        twin_metric = _scalar_metric(b, X, fields[4:8])
        ordinary_exponent = fields[0] + 2 * fields[2] + fields[3]
        twin_exponent = fields[4] + 2 * fields[6] + fields[7]
        mixed_exponent = (ordinary_exponent + twin_exponent) / 2
        volume = _scale((sp.S.One, ordinary_exponent, ordinary_exponent**2 / 2), a**3)
        twin_volume = _scale((sp.S.One, twin_exponent, twin_exponent**2 / 2), X * b**3)
        mixed_volume = _scale(
            (sp.S.One, mixed_exponent, mixed_exponent**2 / 2), sp.sqrt(X) * (a * b) ** sp.Rational(3, 2)
        )
    else:
        metric = _vector_metric(a, sp.S.One, fields[:2])
        twin_metric = _vector_metric(b, X, fields[2:])
        volume = (a**3, sp.S.Zero, -(a**3) * fields[0] ** 2 / 2)
        twin_volume = (X * b**3, sp.S.Zero, -X * b**3 * fields[2] ** 2 / 2)
        mixed_volume = _scale(
            (sp.S.One, sp.S.Zero, -(fields[0] ** 2 + fields[2] ** 2) / 4),
            sp.sqrt(X) * (a * b) ** sp.Rational(3, 2),
        )
    inverse, twin_inverse = _inverse(metric), _inverse(twin_metric)
    connection = _connection(
        inverse, _metric_derivatives(metric, symbols, fields, time_derivatives, spatial_derivatives)
    )
    twin_connection = _connection(
        twin_inverse, _metric_derivatives(twin_metric, symbols, fields, time_derivatives, spatial_derivatives)
    )
    difference = {}
    for indices in sorted(connection.keys() | twin_connection.keys()):
        jet = _sum(connection.get(indices, _ZERO), _scale(twin_connection.get(indices, _ZERO), -sp.S.One))
        if _nonzero(jet):
            difference[indices] = jet
    z = _scale(_contract_upsilon(inverse, difference), -1 / (2 * a0**2))
    interaction_function = (
        symbols["Mbar"],
        symbols["Mprime"] * z[1],
        symbols["Mprime"] * z[2] + symbols["Msecond"] * z[1] ** 2 / 2,
    )
    lagrangian = _sum(
        _scale(_product(volume, _contract_upsilon(inverse, connection)), symbols["beta"] / 2),
        _scale(_product(twin_volume, _contract_upsilon(twin_inverse, twin_connection)), symbols["alpha"] / 2),
        _scale(_product(mixed_volume, interaction_function), a0**2),
    )
    if include_matter:
        lagrangian = _sum(
            lagrangian, _ordinary_matter(inverse, volume, symbols, fields, time_derivatives, spatial_derivatives)
        )
    return QuadraticAction(
        sector, symbols, fields, time_derivatives, spatial_derivatives, *lagrangian, *z, include_matter=include_matter
    )


def principal_part(action: QuadraticAction) -> sp.Expr:
    """Select total degree two in perturbation derivatives from lagrangian2.

    U,V,W are background rates and do not count toward this degree. This helper
    neither fixes a gauge nor integrates constraints or performs integrations by
    parts. Lower-derivative terms in the full action can remain essential to
    those subsequent operations.
    """
    derivative_symbols = action.time_derivatives + action.spatial_derivatives
    polynomial = sp.Poly(action.lagrangian2, *derivative_symbols)
    return sp.Add(
        *(
            coefficient * sp.Mul(*(symbol**power for symbol, power in zip(derivative_symbols, powers)))
            for powers, coefficient in polynomial.terms()
            if sum(powers) == 2
        )
    )
