#!/usr/bin/env python
"""Compare lagCAMB and AxiCLASS FDM model-to-null power-spectrum ratios."""

import argparse
import hashlib
from pathlib import Path

import numpy as np


def load_lag(path):
    data = np.load(path)
    return np.asarray(data["k"]), np.asarray(data["pk"])


def load_class(path):
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def ratio_on_grid(model_k, model_pk, null_k, null_pk, grid):
    model_log_pk = np.interp(np.log(grid), np.log(model_k), np.log(model_pk))
    null_log_pk = np.interp(np.log(grid), np.log(null_k), np.log(null_pk))
    return np.exp(model_log_pk - null_log_pk)


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_spectrum(name, k, pk):
    if k.ndim != 1 or pk.ndim != 1 or k.size != pk.size or k.size < 2:
        raise ValueError(f"{name} must contain equal one-dimensional k and P(k) arrays")
    if not np.all(np.isfinite(k)) or not np.all(np.isfinite(pk)):
        raise ValueError(f"{name} contains non-finite values")
    if np.any(k <= 0) or np.any(pk <= 0) or np.any(np.diff(k) <= 0):
        raise ValueError(f"{name} requires positive P(k) on a strictly increasing positive k grid")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lag-model", required=True, help="NPZ with k and pk arrays")
    parser.add_argument("--lag-null", required=True, help="NPZ with k and pk arrays")
    parser.add_argument("--axi-model", required=True, help="AxiCLASS z=0 *_z1_pk.dat")
    parser.add_argument("--axi-null", required=True, help="AxiCLASS z=0 *_z1_pk.dat")
    parser.add_argument("--kmax", required=True, type=float, help="maximum k in h/Mpc")
    parser.add_argument("--tolerance", type=float, default=0.012)
    args = parser.parse_args()

    lag_model_k, lag_model_pk = load_lag(args.lag_model)
    lag_null_k, lag_null_pk = load_lag(args.lag_null)
    axi_model_k, axi_model_pk = load_class(args.axi_model)
    axi_null_k, axi_null_pk = load_class(args.axi_null)
    for name, k, pk in (
        ("lag model", lag_model_k, lag_model_pk),
        ("lag null", lag_null_k, lag_null_pk),
        ("AxiCLASS model", axi_model_k, axi_model_pk),
        ("AxiCLASS null", axi_null_k, axi_null_pk),
    ):
        validate_spectrum(name, k, pk)

    lower = max(lag_model_k.min(), lag_null_k.min(), axi_model_k.min(), axi_null_k.min())
    available_upper = min(lag_model_k.max(), lag_null_k.max(), axi_model_k.max(), axi_null_k.max())
    if available_upper < args.kmax * (1 - 1e-12):
        raise ValueError(
            f"Requested kmax={args.kmax:g} h/Mpc is not covered; common data end at "
            f"{available_upper:g} h/Mpc"
        )
    upper = args.kmax
    eps = 1e-12
    grid = lag_model_k[(lag_model_k >= lower * (1 - eps)) & (lag_model_k <= upper * (1 + eps))]
    if grid.size == 0:
        raise ValueError("No common k samples in the requested interval")

    lag_ratio = ratio_on_grid(lag_model_k, lag_model_pk, lag_null_k, lag_null_pk, grid)
    axi_ratio = ratio_on_grid(axi_model_k, axi_model_pk, axi_null_k, axi_null_pk, grid)
    error = np.abs(lag_ratio / axi_ratio - 1)
    cumulative = np.maximum.accumulate(error)
    accepted = grid[cumulative <= args.tolerance]

    print(f"samples={grid.size}")
    print(f"max_abs_double_ratio={error.max():.8e}")
    print(f"k_at_max_error={grid[np.argmax(error)]:.8g} h/Mpc")
    if accepted.size:
        print(f"largest_prefix_k_at_tolerance={accepted[-1]:.8g} h/Mpc")
    else:
        print("largest_prefix_k_at_tolerance=none")
    for label, path in (
        ("lag_model_sha256", args.lag_model),
        ("lag_null_sha256", args.lag_null),
        ("axi_model_sha256", args.axi_model),
        ("axi_null_sha256", args.axi_null),
    ):
        print(f"{label}={file_digest(path)}")
    if error.max() > args.tolerance:
        raise SystemExit(
            f"FAIL: maximum double-ratio error {error.max():.8e} exceeds tolerance {args.tolerance:.8e}"
        )
    print(f"PASS: maximum double-ratio error <= {args.tolerance:.8e}")


if __name__ == "__main__":
    main()
