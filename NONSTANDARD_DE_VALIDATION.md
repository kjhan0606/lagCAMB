# Nonstandard dark-energy physical-validation status

Date: 2026-08-20

This record separates a missing validation record from a failed physical
validation. Historical batteries were recovered from
`/home/kjhan/BACKUP/lagCAMB_validation/`; their validated model commits are in
the ancestry of the current lagCAMB tree. A compact set of decisive invariants
was rerun against the current build after adding parameter-domain guards.

## Science status

| Model/path | Status | Evidence and supported domain |
|---|---|---|
| TrackerQuintessence | Validated at the archived benchmark points | Historical commit `dfec4c6`: exponential model agrees with CLASS to 0.168% in linear P(k), RP tracker asymptote and lagRamses background agree, and shooting closes at about 1e-10. Current parameter guards reject invalid potential/IC combinations. |
| CoupledQuintessence | Validated for the tested potentials and `0<=beta<=0.1` | Historical commit `67bd5ed`: 9/9 tests pass, including the independently derived Euler sign, CDM mass correction, budget closure, and positive beta-squared growth. Current rerun gives a bit-exact beta=0 limit and growth-enhancement ratio 3.810 for beta=0.06 versus 0.03 (expected about 4). Negative or larger beta is outside the validated contract. |
| RunningVacuum | Analytic background and internal-response validation for `abs(nu)<=0.01` | Historical commit `a9f8d1e`. Current rerun: nu=0 is bit-exact LambdaCDM; maximum CDM log-slope error is 1.0e-13, vacuum-density error is 1.8e-15, and sigma8 decreases monotonically as nu increases. The perturbation response has internal consistency tests but no independent Boltzmann-code comparison. |
| KEssence | Validated for `0.5<x0<=0.50001` | Historical commit `c0341c8`: 6/6 tests pass. Current rerun at x0=0.50001 gives an analytic-w error below 1.8e-11 and a finite positive spectrum. x0=0.5001 is a robustness-only clean run, not a cosmologically viable benchmark; x0=0.6 can break recombination. The enforced default is 0.500001. |
| Chaplygin | Analytic background plus limited perturbation benchmarks | Historical commit `ed35899`: 5/5 tests pass and reproduce the known Jeans pathology. Current analytic-background errors are 8.72e-9 in density and 1.96e-11 in w; As=0.9999 and 0.99999 approach LambdaCDM with maximum P(k) differences 1.10e-3 and 1.34e-4. This is not a precision validation of the full allowed parameter volume. |
| AxionEffectiveFluid / EDE | Conditionally validated | Historical AxiCLASS comparison: the vanishing-EDE limit agrees at 18.6 ppm. Peak-matched nominal models differ at the few-percent level because the lagCAMB transition and AxiCLASS peak parameterizations are not identical. Claims must identify this as a fluid approximation. |
| InteractingDE type 1 | Validated only under the enforced contract | `Q=xi H rho_de`; exact w=-1 or `1+w(a)>1e-6`; validated sound speed is cs2=1. Types 2 and 3 are disabled because their continuity equations are incomplete. See `INTERACTINGDE_FIX_VALIDATION.md`. |
| FuzzyDM exact KG to PH-EFA | Conditionally validated at `z=0`, `f_axion=0.05` | The zero-abundance path is bit-exact standard CDM and parameter guards fail closed. The active path shoots the exact quadratic KG background to the requested density, evolves exact field perturbations, maps every mode on one self-consistent `m/H` surface, and joins to PH-EFA with a conservative C1 blend. The AxiCLASS model-to-null double-ratio residuals are 1.148%, 1.106%, and 0.870% over `k<=10,30,50 h/Mpc` for `m=1e-23,1e-22,1e-21 eV`. Match-ratio 50-to-75 changes are below 5e-5. The wider abundance interval is a guard, not a precision domain; CMB and lensing lack an independent exact comparison. |
| FuzzyDMField active axion path | Failed; disabled | For `m=1e-22 eV`, `f_axion=0.05`, only 0.578 of the target axion density is recovered and the total present-day DE-slot budget residual is -7.98e-3. The code gives `a_osc=6.43e-7`, then its incorrect `m/(aH)` matching condition collapses `a_match` to the same epoch; the intended `m/H=100` epoch is about `3.73e-6`. Density before the integration start is zero and returned pre-match pressure omits the KG pressure. Only the zero-axion null path is accepted. |
| HorndeskiDE active MG path | Failed; disabled | The scalar perturbation and stability system is incomplete, the tested braiding response is not traceable in the matter spectrum, and no independent Horndeski-solver comparison has passed. A null scalar-spectrum response from alpha_K or alpha_T alone is not diagnostic. Only exact LambdaCDM is accepted; separately tested model-specific `pars.MG` paths are distinct. |

The superseded all-time-EFA FuzzyDM baseline retained a standard-CDM
background. Its corresponding three AxiCLASS residuals were 1.145%, 1.098%,
and 0.640%. Those historical values describe neither the current exact-KG
background nor its common PH-EFA matching surface and are retained only for
implementation provenance.

### FuzzyDM benchmark provenance

The independent reference is AxiCLASS commit
`ba4ede7b1d735aa6312ab5f4355d26b5e617e70c`. The z=0 benchmark uses
`H0=67.36`, `omega_b=0.02237`, total `omega_cdm+omega_axion=0.12`,
`A_s=2.1e-9`, `n_s=0.965`, `tau_reio=0.0544`, `N_ur=3.044`, no ncdm,
and `f_axion=0.05`. AxiCLASS uses `scf_potential=axionquad`,
`scf_evolve_as_fluid=yes`, `scf_evolve_like_axionCAMB=yes`,
`threshold_scf_fluid_m_over_H=3`, shooting enabled, perturbations enabled,
and the scalar included in `delta_m`. The reported statistic is

`abs((P_lag,FDM/P_lag,LCDM)/(P_Axi,FDM/P_Axi,LCDM)-1)`

evaluated on a common 500-point log-k lagCAMB grid after log-k interpolation
of the AxiCLASS ratio and restricted
to the stated maximum k. The three mass points test only `f_axion=0.05`; the
enforced abundance interval is a supported guard contract, not a cross-code
precision validation of every abundance.

## Recovered primary records

- Tracker: `quint/VALIDATION_D.md`
- K-essence: `quint/KESSENCE_VALIDATION.md`
- Coupled quintessence: `quint/COUPLEDQUINT_VALIDATION.md`
- Running vacuum: `ext/runvac/RVM_VALIDATION.md`
- Generalized Chaplygin gas: `ext/chaplygin/GCG_VALIDATION.md`
- EDE/AxiCLASS: `ede/RESULTS_PHASE_C.md`
- Cross-model limits: `limits/RESULTS_LIMITS.md`

All paths above are relative to `/home/kjhan/BACKUP/lagCAMB_validation/`.
Historical numerical artifacts are preserved in that archive and were not
overwritten by the current rerun.

An older user-owned script, `tests/validate_models.py`, contains a section named
`TEST 12: FuzzyDMField Background & EFA`, but it only prints sigma8 and transfer
ratios and has no PASS/FAIL assertion. The formal `tests/validation/report.json`
PASS entries named `FuzzyDM` test `camb.dark_matter.FuzzyDM`, not
`FuzzyDMField`. These records explain the earlier recollection of an Opus PASS
claim but do not validate the active Klein--Gordon/EFA field path.

## Current-head rerun acceptance criteria

The current rerun checks null limits, exact backgrounds where closed forms
exist, physical response directions, positivity/finite spectra, and known
coupling scaling. Passing a numerical smoke test is not by itself sufficient:
an active model is disabled when a required physical equation or independent
reference comparison is absent. Python and Fortran validation both enforce the
same fail-closed contract so direct constructors, INI input, or low-level field
assignment cannot silently bypass it.
