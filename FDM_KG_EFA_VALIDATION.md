# FDM exact KG to PH-EFA conditional validation

## Scope

The active `camb.dark_matter.FuzzyDM` path now replaces a fraction of the
total `omch2` budget with a quadratic scalar. It evolves the exact
Klein--Gordon background and perturbations, maps to Passaglia--Hu EFA
variables on one common `m/H` surface, and blends stress-energy with a
conservative C1 window. The older DarkEnergy-slot `FuzzyDMField` remains
disabled.

The guarded domain is:

- `1e-23 <= m_axion/eV <= 1e-21`
- `1e-3 <= f_axion <= 0.1`
- `50 <= match_ratio <= 75`
- cosmological-constant dark energy

## Numerical checks

The benchmark cosmology and exact AxiCLASS input files are in
`validation/fdm/`. The AxiCLASS reference commit is
`ba4ede7b1d735aa6312ab5f4355d26b5e617e70c`.

At `f_axion=0.05` and `match_ratio=50`, the maximum absolute difference of
the lagCAMB and AxiCLASS model-to-null ratios is:

| mass [eV] | k range [h/Mpc] | maximum double-ratio residual |
|---:|---:|---:|
| 1e-23 | <= 10 | 1.1481% |
| 1e-22 | <= 30 | 1.1063% |
| 1e-21 | <= 50 | 0.8699% |

Changing `match_ratio` from 50 to 75 changes the lagCAMB spectra by at most
`4.77e-5`, `3.83e-5`, and `3.06e-5` over those respective ranges. A delayed
switch at 100 is not enabled: the `1e-21 eV` exact perturbation becomes stiff
in the current CAMB integrator before that surface.

The public background-density output is finite over sampled frozen,
transition, and EFA epochs. At `a=1`, its component sum closes exactly in the
reported arithmetic, and the reported total DM density agrees with
`get_Omega("cdm")`. The common numerical surface is iterated against the exact
KG contribution to `H(a)`; the regression requires `abs((m/H)/match_ratio-1)
<= 1e-8` and currently closes near `3e-14` at the central benchmark.

This is a conditional validation of the reported `z=0` linear-power
benchmarks at `f_axion=0.05`. The wider abundance interval is a fail-closed
guard, not a cross-code precision claim. CMB and lensing spectra do not yet
have an independent exact-solver comparison.

## Regression commands

```console
python setup.py clean
python setup.py make
python -m unittest camb.tests.camb_test
python validation/fdm/run_lagcamb.py --output /tmp/lag-null.npz
python validation/fdm/run_lagcamb.py --output /tmp/lag-fdm.npz --mass 1e-23
python validation/fdm/compare_model_to_null.py \
  --lag-model /tmp/lag-fdm.npz --lag-null /tmp/lag-null.npz \
  --axi-model /tmp/lagcamb_fdm_axi_m1e23_00_z1_pk.dat \
  --axi-null /tmp/lagcamb_fdm_lcdm_00_z1_pk.dat --kmax 10 --tolerance 0.012
```

The full CAMB unit suite completed 23 tests successfully. Sphinx HTML and the
40-page PDF manual also build; the Sphinx run retains unrelated pre-existing
docstring/reference warnings.
