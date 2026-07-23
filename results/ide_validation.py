#!/usr/bin/env python
"""InteractingDE background-fix validation battery.
Writes CSVs (prefix ide_) into this directory and prints a machine-readable summary.
Run from repo root with the quint venv active:  OMP_NUM_THREADS=8 python results/ide_validation.py
"""
import os, json, sys
import numpy as np
import camb
from camb.dark_energy import InteractingDE

HERE = os.path.dirname(os.path.abspath(__file__))
YHE = 0.24585548197552395  # pinned BBN float (numpy-2.x ctypes workaround)

COMMON = dict(H0=67.36, ombh2=0.02237, omch2=0.12, As=2.1e-9, ns=0.9649,
              tau=0.0544, YHe=YHE, mnu=0.0, nnu=3.044, num_massive_neutrinos=0)
ZS = [0.0, 1.0, 9.0, 49.0, 99.0]
KMAX = 20.0

def base_pars():
    p = camb.set_params(**COMMON)
    return p

def lcdm_pars():
    return base_pars()  # default DarkEnergy = cosmological constant (w=-1)

def ide_pars(xi, w=-1.0, wa=0.0):
    p = base_pars()
    de = InteractingDE()
    de.set_params(w=w, wa=wa, xi_ide=xi, interaction_type=1, cs2_ide=1.0)
    p.DarkEnergy = de
    return p

def pk_interp(pars, zs):
    pars.set_matter_power(redshifts=list(zs), kmax=KMAX)
    pars.NonLinear = camb.model.NonLinear_none
    res = camb.get_results(pars)
    PK = res.get_matter_power_interpolator(nonlinear=False, hubble_units=False,
                                           k_hunit=False, var1='delta_tot', var2='delta_tot')
    return res, PK

def sigma8_of(pars):
    pars.set_matter_power(redshifts=[0.0], kmax=KMAX)
    pars.NonLinear = camb.model.NonLinear_none
    res = camb.get_results(pars)
    return float(res.get_sigma8_0()), res

summary = {}

# ---------------------------------------------------------------- Test 1: regression
kgrid = np.logspace(-3, np.log10(5.0), 300)  # 1/Mpc, band [1e-3, 5]
res_l, PK_l = pk_interp(lcdm_pars(), ZS)
res_0, PK_0 = pk_interp(ide_pars(0.0), ZS)
s8_l = float(res_l.get_sigma8_0()); s8_0 = float(res_0.get_sigma8_0())
maxdev = {}
for z in ZS:
    pl = PK_l.P(z, kgrid); p0 = PK_0.P(z, kgrid)
    maxdev[z] = float(np.max(np.abs(p0/pl - 1.0)))
reg_exact = (abs(s8_0 - s8_l) == 0.0) and all(v == 0.0 for v in maxdev.values())
summary['test1_regression'] = dict(sigma8_lcdm=s8_l, sigma8_xi0=s8_0,
    sigma8_absdiff=abs(s8_0 - s8_l), maxdev_by_z=maxdev, bit_exact=bool(reg_exact))
print("[T1] LCDM sigma8=%.7f  xi0 sigma8=%.7f  |diff|=%.2e  bit_exact=%s"
      % (s8_l, s8_0, abs(s8_0 - s8_l), reg_exact))

# ---------------------------------------------------------------- Test 2: background exactness
a_arr = np.logspace(np.log10(0.01), 0.0, 400)
t2 = {}
for xi in (0.03, -0.03):
    res = camb.get_background(ide_pars(xi))
    d = res.get_background_densities(a_arr, vars=['cdm', 'de'])
    de = np.asarray(d['de']); cdm = np.asarray(d['cdm'])   # 8piG rho a^4
    rho_de = de / a_arr**4
    # log-slope of rho_de over a in [0.01,1]
    sl, _ = np.polyfit(np.log(a_arr), np.log(rho_de), 1)
    slope_err = abs(sl - (-xi))
    # analytic rho_c in "x a^4" units.  grhoc=cdm(a=1), grhov=de(a=1)
    grhoc = cdm[-1]; grhov = de[-1]
    p = xi  # w=-1 => p = 3(1+w)+xi = xi
    Ccode = xi * grhov / (3.0 - p)
    cdm_analytic = (grhoc - Ccode) * a_arr + Ccode * a_arr**(4.0 - p)
    rc_relerr = float(np.max(np.abs(cdm / cdm_analytic - 1.0)))
    # H(a=1) budget: total 8piG rho a^4 at a=1 vs grhocrit (flat) -> H0
    H0 = float(res.hubble_parameter(0.0))
    t2['xi=%+.2f' % xi] = dict(logslope_rho_de=float(sl), expected=-xi,
        slope_abserr=float(slope_err), rho_c_max_relerr=rc_relerr, H0=H0)
    # dump csv
    np.savetxt(os.path.join(HERE, 'ide_bg_xi%+.2f.csv' % xi),
        np.column_stack([a_arr, rho_de, cdm, cdm_analytic]),
        header='a, rho_de(8piGrho), cdm_code(8piGrho_a4), cdm_analytic', delimiter=',')
    print("[T2] xi=%+.2f  rho_de logslope=%.8f (exp %.8f, err %.2e)  rho_c max relerr=%.2e  H0=%.5f"
          % (xi, sl, -xi, slope_err, rc_relerr, H0))
summary['test2_background'] = t2

# ---------------------------------------------------------------- Test 3: physical response
s8p, resp = sigma8_of(ide_pars(0.03))
s8m, resm = sigma8_of(ide_pars(-0.03))
_, PK_p = pk_interp(ide_pars(0.03), ZS)
_, PK_m = pk_interp(ide_pars(-0.03), ZS)
ratios = {}
for z in (0.0, 49.0):
    rp = PK_p.P(z, kgrid) / PK_l.P(z, kgrid)
    rm = PK_m.P(z, kgrid) / PK_l.P(z, kgrid)
    ratios['z=%g' % z] = dict(xi_plus_max=float(np.max(rp-1)), xi_plus_min=float(np.min(rp-1)),
                              xi_minus_max=float(np.max(rm-1)), xi_minus_min=float(np.min(rm-1)))
    np.savetxt(os.path.join(HERE, 'ide_pkratio_z%g.csv' % z),
        np.column_stack([kgrid, rp, rm]),
        header='k[1/Mpc], P(xi+0.03)/P_LCDM, P(xi-0.03)/P_LCDM', delimiter=',')
summary['test3_response'] = dict(sigma8_lcdm=s8_l, sigma8_xi_plus=s8p, sigma8_xi_minus=s8m,
    shift_plus_pct=100*(s8p/s8_l-1), shift_minus_pct=100*(s8m/s8_l-1), pk_ratio=ratios)
print("[T3] sigma8: LCDM=%.7f  xi=+0.03 -> %.7f (%+.3f%%)  xi=-0.03 -> %.7f (%+.3f%%)"
      % (s8_l, s8p, 100*(s8p/s8_l-1), s8m, 100*(s8m/s8_l-1)))

# ---------------------------------------------------------------- Test 4: CPL coupling
s8_cpl, res_cpl = sigma8_of(ide_pars(0.03, w=-0.9, wa=0.1))
# background integral path exercised: check finite & rho_de slope near expected
res_cplb = camb.get_background(ide_pars(0.03, w=-0.9, wa=0.1))
dcpl = res_cplb.get_background_densities(a_arr, vars=['cdm', 'de'])
cpl_finite = bool(np.all(np.isfinite(dcpl['cdm'])) and np.all(np.isfinite(dcpl['de'])) and np.isfinite(s8_cpl))
summary['test4_cpl'] = dict(sigma8=s8_cpl, w=-0.9, wa=0.1, xi=0.03, finite=cpl_finite)
print("[T4] CPL w=-0.9 wa=0.1 xi=0.03  sigma8=%.7f  finite=%s" % (s8_cpl, cpl_finite))

# ---------------------------------------------------------------- Test 5: perturbation sanity
cases = {'xi=0': ide_pars(0.0), 'xi=+0.03': ide_pars(0.03), 'xi=-0.03': ide_pars(-0.03),
         'cpl_xi=0.03': ide_pars(0.03, w=-0.9, wa=0.1)}
sanity = {}
for name, pp in cases.items():
    _, PK = pk_interp(pp, ZS)
    ok = True; details = {}
    for z in ZS:
        pv = PK.P(z, kgrid)
        pos = bool(np.all(pv > 0)); fin = bool(np.all(np.isfinite(pv)))
        details['z=%g' % z] = dict(all_positive=pos, all_finite=fin, min=float(np.min(pv)))
        ok = ok and pos and fin
    sanity[name] = dict(pass_=ok, detail=details)
    print("[T5] %-12s  P(k)>0 & finite at all z: %s" % (name, ok))
summary['test5_sanity'] = sanity

with open(os.path.join(HERE, 'ide_summary.json'), 'w') as f:
    json.dump(summary, f, indent=2)
print("\nWrote", os.path.join(HERE, 'ide_summary.json'))
