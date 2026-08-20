# FuzzyDM cross-code benchmark

Reference solver: AxiCLASS commit
`ba4ede7b1d735aa6312ab5f4355d26b5e617e70c`. The `*.ini.reference` files
are the exact z=0 benchmark inputs for the LambdaCDM null and the three
`f_axion=0.05` mass points. AxiCLASS accepts these filenames directly.

The lagCAMB side uses the same cosmology, `omch2=0.12`, and the exact
KG-to-PH-EFA path
`FuzzyDM.set_params(m_axion=<mass>, f_axion=0.05, match_ratio=50)`. Generate
z=0 linear spectra containing `k` in h/Mpc and `pk` in (Mpc/h)^3, then run the
comparison script with the four model and null paths:

```console
python validation/fdm/run_lagcamb.py --output /tmp/lag-null.npz
python validation/fdm/run_lagcamb.py --output /tmp/lag-fdm.npz --mass 1e-23
python validation/fdm/compare_model_to_null.py --lag-model /tmp/lag-fdm.npz --lag-null /tmp/lag-null.npz --axi-model /tmp/axi-fdm_z1_pk.dat --axi-null /tmp/axi-null_z1_pk.dat --kmax 10 --tolerance 0.012
```

The statistic is the absolute difference of the two model-to-null ratios,
not the raw power-spectrum difference. The recorded all-time-EFA maxima were:

| mass [eV] | maximum k [h/Mpc] | maximum error |
|---:|---:|---:|
| 1e-23 | 10 | 1.145% |
| 1e-22 | 30 | 1.098% |
| 1e-21 | 50 | 0.640% |

The sampled-grid sub-percent prefix limits are approximately `k=4.6`, `26`,
and at least `50 h/Mpc`, respectively. These numbers validate only
`f_axion=0.05` and must not be
extrapolated to the full abundance guard interval.

The exact implementation gives maximum double-ratio residuals of 1.148%,
1.106%, and 0.870% over the same three k ranges. Varying `match_ratio` from
50 to 75 changes the lagCAMB spectra by at most `4.77e-5`, `3.83e-5`, and
`3.06e-5`, respectively. `match_ratio=100` is intentionally outside the
active contract because the `1e-21 eV` exact perturbation becomes numerically
stiff before the delayed switch.

Both comparison scripts validate positivity, finiteness, and grid ordering,
print SHA-256 provenance, and exit nonzero when the tolerance is exceeded.
After generating the six match-ratio spectra, the three-point convergence
regression is:

```console
python validation/fdm/compare_match_convergence.py \
  --pair 1e-23 /tmp/lag-m1e23-n50.npz /tmp/lag-m1e23-n75.npz 10 \
  --pair 1e-22 /tmp/lag-m1e22-n50.npz /tmp/lag-m1e22-n75.npz 30 \
  --pair 1e-21 /tmp/lag-m1e21-n50.npz /tmp/lag-m1e21-n75.npz 50
```

The DarkMatter-slot path is conditionally validated only for these reported
`z=0`, `f_axion=0.05` linear-power benchmarks. The wider abundance interval
is a fail-closed guard, not a precision domain. CMB and lensing outputs have
not yet received an independent exact-solver comparison.
