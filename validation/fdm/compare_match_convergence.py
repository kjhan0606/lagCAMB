#!/usr/bin/env python
"""Fail-closed comparison of FDM spectra at two KG-to-EFA match ratios."""

import argparse
import hashlib
from pathlib import Path

import numpy as np


def load_spectrum(path):
    data = np.load(path)
    k = np.asarray(data["k"])
    pk = np.asarray(data["pk"])
    if k.ndim != 1 or pk.ndim != 1 or k.size != pk.size or k.size < 2:
        raise ValueError(f"Invalid spectrum arrays in {path}")
    if not np.all(np.isfinite(k)) or not np.all(np.isfinite(pk)):
        raise ValueError(f"Non-finite spectrum in {path}")
    if np.any(k <= 0) or np.any(pk <= 0) or np.any(np.diff(k) <= 0):
        raise ValueError(f"Spectrum in {path} is not positive on an increasing k grid")
    return k, pk


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pair",
        nargs=4,
        action="append",
        metavar=("MASS_EV", "MATCH50_NPZ", "MATCH75_NPZ", "KMAX"),
        required=True,
        help="repeat once per mass point",
    )
    parser.add_argument("--tolerance", type=float, default=6e-5)
    args = parser.parse_args()

    failed = False
    for mass, path50, path75, kmax_text in args.pair:
        kmax = float(kmax_text)
        k50, pk50 = load_spectrum(path50)
        k75, pk75 = load_spectrum(path75)
        lower = max(k50.min(), k75.min())
        available_upper = min(k50.max(), k75.max())
        if available_upper < kmax * (1 - 1e-12):
            raise ValueError(
                f"Requested kmax={kmax:g} h/Mpc for m={mass} eV is not covered; "
                f"common data end at {available_upper:g} h/Mpc"
            )
        upper = kmax
        grid = k50[(k50 >= lower) & (k50 <= upper)]
        if grid.size == 0:
            raise ValueError(f"No common samples for m={mass} eV")
        log_pk75 = np.interp(np.log(grid), np.log(k75), np.log(pk75))
        error = np.abs(np.exp(log_pk75 - np.log(pk50[(k50 >= lower) & (k50 <= upper)])) - 1)
        maximum = error.max()
        print(f"mass_eV={mass} samples={grid.size} max_abs_ratio={maximum:.8e}")
        print(f"match50_sha256={digest(path50)}")
        print(f"match75_sha256={digest(path75)}")
        failed = failed or maximum > args.tolerance

    if failed:
        raise SystemExit(f"FAIL: match-ratio convergence exceeds {args.tolerance:.8e}")
    print(f"PASS: every match-ratio comparison <= {args.tolerance:.8e}")


if __name__ == "__main__":
    main()
