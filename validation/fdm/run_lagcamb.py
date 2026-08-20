#!/usr/bin/env python
"""Generate a lagCAMB z=0 linear spectrum for the FDM benchmark."""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import camb
from camb import dark_matter


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, help="output NPZ path")
    parser.add_argument("--mass", type=float, help="axion mass in eV; omit for the LCDM null")
    parser.add_argument("--fraction", type=float, default=0.05)
    parser.add_argument("--match-ratio", type=int, default=50)
    parser.add_argument("--kmax", type=float, default=50.0)
    parser.add_argument("--samples", type=int, default=500)
    args = parser.parse_args()

    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=67.36,
        ombh2=0.02237,
        omch2=0.12,
        mnu=0,
        nnu=3.044,
        tau=0.0544,
    )
    pars.InitPower.set_params(As=2.1e-9, ns=0.965)
    pars.WantTransfer = True
    pars.set_matter_power(redshifts=[0], kmax=max(1.2 * args.kmax, args.kmax + 1), silent=True)
    pars.NonLinear = camb.model.NonLinear_none
    if args.mass is not None:
        pars.DarkMatter = dark_matter.FuzzyDM()
        pars.DarkMatter.set_params(
            m_axion=args.mass,
            f_axion=args.fraction,
            match_ratio=args.match_ratio,
        )

    results = camb.get_results(pars)
    k, _, pk = results.get_matter_power_spectrum(
        minkh=1e-4,
        maxkh=args.kmax,
        npoints=args.samples,
    )
    np.savez(args.output, k=k, pk=pk[0])


if __name__ == "__main__":
    main()
