# InteractingDE background-fix validation

Date 2026-07-23. Branch `interacting-de-fix` (off `kessence-px`). camb 1.6.7 (lagCAMB clone),
gfortran 13.2.0, quint venv. Two background bugs in the Valiviita-type interacting dark
energy model (`camb.dark_energy.InteractingDE`, Fortran `TInteractingDE`, `interaction_type=1`,
`Q = xi*H*rho_de`) are fixed and validated. All runs are lagCAMB-internal (each model vs its
own LCDM). Driver `results/ide_validation.py`; artifacts `results/ide_*`.

Common cosmology: H0=67.36, ombh2=0.02237, omch2=0.12, As=2.1e-9, ns=0.9649, tau=0.0544,
YHe=0.24585548197552395 (pinned float), massless nu (mnu=0, nnu=3.044, num_massive_neutrinos=0).
Linear P(k), kmax=20/Mpc. P via `get_matter_power_interpolator(nonlinear=False,
hubble_units=False, k_hunit=False)`, band k in [1e-3, 5]/Mpc. z in {0,1,9,49,99}. OMP_NUM_THREADS=8.

## The two bugs (as found by the limit-check campaign)

1. **`grhov_t` missing the `/a^2`.** `InteractingDE.f90` set `grhov_t = grhov*a^(1-3w-3wa-xi)`
   for type 1, but `grhov_t` must be `8πG rho_de(a) a^2`; the standard path divides `grho_de`
   by `a^2` (`DarkEnergyInterface.f90:96`, exponent `-1-3w-3wa`). The line was a factor `a^2`
   too large (its own comment at :156 labelled it "(wrong)"). Effect: `rho_de ~ a^{2-xi}`,
   spuriously suppressed by z>=1.
2. **CDM continuity source `Q` never applied.** The type-1 background left CDM at exactly
   `a^{-3}` (comment at :124 promised a "compensating modification" that was absent). The
   defining energy exchange `rho_c' + 3H rho_c = Q` was missing at background level.

Because both are even/monotonic in the buggy exponent rather than linear in xi, the old model
gave **both signs of xi INCREASing** sigma8 (+5.955% at +0.03, +8.576% at -0.03) — an artifact.

## What changed (files / lines)

| file | lines | change |
|------|-------|--------|
| `fortran/DarkEnergyInterface.f90` | 23, 70-80 | New polymorphic base method `CDM_BackgroundCorrection(grhoc,grhov,a)` returning 0. Lets `results.f90` (compiled before `InteractingDE`) apply the CDM excess by dynamic dispatch — no circular module dep, no `select type`. |
| `fortran/InteractingDE.f90` | 190-199 | **BUG 1 FIX.** Type-1 `grhov_t` now delegates to the base `BackgroundDensityAndPressure` and multiplies by `a^{-xi}` — exact for any w(a) (`rho_de = rho_de_std * a^{-xi}`), and correct `/a^2`. Also handles CPL and tabulated w with no duplicated algebra. |
| `fortran/InteractingDE.f90` | 205-211 | Same `/a^2` fix applied to the type-3 approximate branch. |
| `fortran/InteractingDE.f90` | 37-40 | Type fields for the CPL/tabulated CDM-excess integral table (`lna_ia`, `ia_tab`, `use_ia_table`, `n_ia`; fixed-size 4096, no heap/copy hazards). |
| `fortran/InteractingDE.f90` | 118-143 | `Init`: tabulate `I(a)=∫_a^1 (rho_de/rho_de0) a'^3 dln a'` on a 4096-pt log-a grid, **only** for type-1 CPL/tabulated w (constant w uses the closed form). |
| `fortran/InteractingDE.f90` | 221-265 | **BUG 2 FIX.** Override `CDM_BackgroundCorrection` = `-xi*grhov*I(a)/a` = `8πG a^2 [rho_c(a)-rho_c0 a^{-3}]`. Constant w: closed form `I=(1-a^{3-p})/(3-p)`, `p=3(1+w)+xi`. CPL/tabulated: table lookup. Returns 0 for xi=0 or interaction_type≠1. |
| `fortran/results.f90` | 1274-1279 | `grho_no_de` (background H(a)/distances) adds the CDM excess (`*a^2` -> ×a^4 units). |
| `fortran/results.f90` | 997-1000 | `GetBackgroundDensities` `densities(3)` (the Python `get_background_densities` CDM column) adds the excess. |
| `fortran/equations.f90` | 2357-2359 | `derivs`: perturbation-sector `grhoc_t` (feeds adotoa=H(a), Poisson `dgrho_matter`, and the `Transfer_tot`/`Transfer_nonu` weighting) adds the excess. |
| `fortran/equations.f90` | 3268-3270 | `output` (CMB-source H(a)) adds the excess for consistency. |

### Chosen hook and why

The honest chokepoint for CDM does not exist as a single function in CAMB — CDM background is
inlined as `grhoc/a` (×a^2) / `grhoc*a` (×a^4) in several places, and the `DarkMatter%Background`
override (used by DecayingDM) is not present for a DE-slot model. `results.f90` cannot
`select type(TInteractingDE)` because `InteractingDE` `use`s `results` (circular). So I added a
**polymorphic base method** `TDarkEnergyModel%CDM_BackgroundCorrection` (returns 0 by default,
overridden by `TInteractingDE`) — the same idiom as `BackgroundDensityAndPressure` /
`PerturbationEvolve`. It is applied at exactly the sites where CDM background density enters
H(a) and the matter perturbation weighting. This keeps the DE perturbation weighting untouched
(the CDM excess is **not** smuggled into `grhov_t`, which multiplies `delta_de` in
`PerturbedStressEnergy`); `grhov_t` stays pure `8πG rho_de a^2`, and `grhoc_t` carries the CDM
excess. `Init` caches nothing that depends on call order (grhov is passed in as an argument).
21cm-only spots (`results.f90:2012,2512`, guarded by `Do21cm`) are left untouched — not on the
H(a) or P(k) path for this validation (caveat below).

### DE perturbations (PerturbationEvolve)

Already correct from commit `5c106d9` and consistent with the corrected background: the type-1
Valiviita sources are `delta_c' += xi H (rho_de/rho_c)(delta_de-delta_c)` (`equations.f90:2552`)
and `delta_de' += xi H (delta_c-delta_de)` (`-xi H delta_de` in PerturbationEvolve:290 plus
`+xi H delta_c` in `equations.f90:2555`). Both use `grhov_t`/`grhoc_t`, now corrected, so the
ratios are physically consistent. No change made — verified empirically (no NaN, P(k)>0).

## Validation battery

| # | test | result | verdict |
|---|------|--------|---------|
| 1 | Regression xi=0 (w=-1) vs LCDM | sigma8 identical (0.8228158 vs 0.8228158, \|diff\|=0.0); max\|P/P_LCDM-1\|=0.0 at every z | **PASS** bit-exact |
| 2 | Background exactness (w=-1) | see below | **PASS** (machine) |
| 3 | Physical response (sigma8, P/P_LCDM) | sign-dependent, near-antisymmetric, few % | **PASS** |
| 4 | CPL w=-0.9, wa=0.1, xi=0.03 | runs clean, table path exercised, sigma8=0.7539548 | **PASS** |
| 5 | Perturbation sanity | P(k)>0 and finite at z in {0,1,9,49,99}, all 4 cases | **PASS** no NaN |

### Test 2 — background exactness (from `get_background_densities`, a in [0.01,1], 400 pts)

| quantity | xi=+0.03 | xi=-0.03 | target | verdict |
|----------|----------|----------|--------|---------|
| log-slope `d ln rho_de/d ln a` | -0.03000000 | +0.03000000 | -xi | err 1.9e-15 / 2.1e-17 |
| `rho_de` slope abs err | 1.88e-15 | 2.08e-17 | <1e-6 | **PASS** |
| `rho_c(a)` max rel err vs analytic | 2.22e-16 | 2.22e-16 | <1e-6 | **PASS** |
| H(a=1) [km/s/Mpc] | 67.36000 | 67.36000 | 67.36 | budget closes |

Analytic used: `rho_c(a) = (rho_c0 - C) a^{-3} + C a^{-p}`, `p = 3(1+w)+xi = xi` (w=-1),
`C = xi rho_de0/(3-p)`. Agreement is at machine epsilon. CSVs `results/ide_bg_xi{+,-}0.03.csv`.

### Test 3 — physical response

| config | sigma8 | shift vs LCDM |
|--------|--------|---------------|
| LCDM | 0.8228158 | -- |
| xi=+0.03 (w=-1, type1) | 0.7915718 | **-3.797%** |
| xi=-0.03 (w=-1, type1) | 0.8538709 | **+3.774%** |

Shifts are now **opposite in sign and near-antisymmetric** (contrast the buggy +5.955% / +8.576%
both-up). P(k) ratio-1 bands over k in [1e-3,5]/Mpc (`results/ide_pkratio_z{0,49}.csv`):

| z | xi=+0.03 (min..max) | xi=-0.03 (min..max) |
|---|---------------------|---------------------|
| 0  | -9.88% .. -1.09% (all negative) | +1.10% .. +10.49% (all positive) |
| 49 | -4.81% .. +4.38% | -4.03% .. +4.80% |

**Direction and the reasoning.** The implemented Q sign matches the class's documented
convention `rho_c' + 3H rho_c = +Q`, `Q = xi H rho_de` — so xi>0 transfers energy DE->CDM.
CAMB pins `omch2` to the **present** CDM density; integrating `d(rho_c a^3)/dln a = xi rho_de a^3`
back from rho_c0 then makes the growth-epoch CDM **lower** for xi>0 (rho_c(a=0.3)=0.9745 x LCDM)
and **higher** for xi<0 (1.0250 x LCDM). Less matter during structure formation -> less growth
-> **xi>0 lowers sigma8, xi<0 raises it**, which is exactly what the numbers show (and the
near-perfect antisymmetry confirms a genuine linear-in-xi coupling). Note this is the *opposite*
of the brief's naive "energy into CDM raises growth" — that intuition holds only under
early-time normalization; under CAMB's fixed-present-CDM convention the sign flips, and the
present result is the physically correct one for this code.

### Test 4 — CPL coupling

w=-0.9, wa=0.1, xi=0.03: runs clean, `use_ia_table=.true.` path (4096-pt log-a quadrature of the
CDM-excess integral) exercised, sigma8=0.7539548, all densities finite. **PASS.**

### Test 5 — perturbation sanity

P(k) finite and > 0 at z in {0,1,9,49,99} for xi=0, xi=+0.03, xi=-0.03, and CPL. No NaN.
No Valiviita early-time instability triggered (w=-1 has (1+w)=0 so the DE-perturbation sector is
inert; the CPL w=-0.9 case with xi=+0.03 also stayed finite/positive). **PASS.**

## Caveats

### Science-use quarantine added on 2026-08-20

The public Python, low-level Fortran, and INI entry points now accept only
`interaction_type=1`. Types 2 and 3 fail before a background or perturbation
result can be generated because their background continuity equations are
incomplete. Previously generated type-2 and type-3 artifacts are nonphysical
diagnostics and must not be used for inference.

The IDE perturbation system has no PPF prescription for phantom evolution or a
crossing of `w=-1`. The exact `w=-1, wa=0` case remains supported to preserve the
validated type-1 and LCDM limits. Every other CPL or tabulated input must satisfy
`1+w(a)>1e-6` over the full spline domain. The stored velocity is
`(1+w)v_de`, and the validated sound-speed configuration is `cs2_ide=1`.

The 2026-08-20 quarantine regression passed the Python setter, direct
constructor, direct-field, low-level IDE setter, and low-level w-table bypass
cases. It also covered spline extrema and endpoint rejection, table-to-CPL state
reset, `CAMBparams.validate()`, invalid type-2/type-3 and phantom/crossing INI
inputs, and a valid type-1 INI background run. The complete
`camb.tests.camb_test` suite passed 19 tests in one process after these additions.
Four pre-patch and post-patch type-1 arrays were bitwise identical. Their SHA-256 digests were
`1de22c8c...39d`, `75b556ec...33c`, `8c57afac...fc5`, and
`05017761...9ee` for the LCDM limit, positive coupling, negative coupling, and
safe CPL cases, respectively.

- **Sigma8 direction** contradicts the brief's naive expectation (xi>0 lowers, not raises, sigma8)
  — this is correct given the documented Q sign + CAMB present-day normalization; see Test 3.
- Type 2 (`Q=xi H rho_c`) and type 3 (`Q=xi H (rho_c+rho_de)`) remain in the
  source as unreachable implementation scaffolding. All supported entry points
  reject them, including the `xi=0` case.
- 21cm-window CDM densities (`results.f90:2012,2512`, `Do21cm` only) are not corrected — outside
  this validation's scope (linear matter power, no 21cm windows).
- CPL/tabulated CDM-excess integral is a 4096-pt trapezoid (log-a); constant-w uses the exact
  closed form, so Test 2's machine-precision check is on the closed form, not the table.
