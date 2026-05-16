#!/usr/bin/env python
"""Generate comparison figures for lagCAMB manual.
Each model vs LCDM: CMB TT power spectrum + matter P(k).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import camb

# ---------- style ----------
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Georgia', 'DejaVu Serif'],
    'font.size': 11,
    'axes.labelsize': 13,
    'axes.titlesize': 14,
    'legend.fontsize': 9,
    'figure.dpi': 200,
    'savefig.dpi': 200,
    'axes.grid': True,
    'grid.alpha': 0.3,
})

COLORS = {
    'lcdm': '#555555',
    'model': '#3898EC',
    'model2': '#D97757',
    'model3': '#4EC9B0',
    'ratio_line': '#FFC107',
}
FIGDIR = os.path.join(os.path.dirname(__file__), 'figures')

def get_lcdm():
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
    pars.InitPower.set_params(As=2.1e-9, ns=0.965)
    pars.set_for_lmax(2500, lens_potential_accuracy=1)
    pars.WantTransfer = True
    pars.set_matter_power(redshifts=[0], kmax=10.0)
    results = camb.get_results(pars)
    cls = results.get_cmb_power_spectra(CMB_unit='muK')['total']
    ell = np.arange(cls.shape[0])
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=10, npoints=300)
    return ell, cls[:, 0], kh, pk[0], results

LCDM = None

def ensure_lcdm():
    global LCDM
    if LCDM is None:
        print("  Computing LCDM baseline...")
        LCDM = get_lcdm()
    return LCDM

def plot_comparison(model_name, ell_m, cl_m, kh_m, pk_m,
                    labels=None, extra_cls=None, extra_pks=None):
    """Plot CMB TT and P(k) comparison. Returns fig."""
    ell_l, cl_l, kh_l, pk_l, _ = ensure_lcdm()

    fig, axes = plt.subplots(2, 2, figsize=(11, 8),
                              gridspec_kw={'height_ratios': [3, 1]})
    fig.suptitle(f'{model_name}', fontsize=16, fontweight='bold', y=0.98)

    # -- CMB TT --
    ax = axes[0, 0]
    lmax = min(len(cl_l), len(cl_m), 2501)
    ell_plot = np.arange(2, lmax)
    fac = ell_plot * (ell_plot + 1) / (2 * np.pi)
    ax.plot(ell_plot, fac * cl_l[2:lmax], color=COLORS['lcdm'],
            label=r'$\Lambda$CDM', lw=1.2)
    ax.plot(ell_plot, fac * cl_m[2:lmax], color=COLORS['model'],
            label=labels[0] if labels else model_name, lw=1.5, ls='--')
    if extra_cls:
        cols = [COLORS['model2'], COLORS['model3']]
        for i, (ecl, elab) in enumerate(extra_cls):
            lmx = min(len(ecl), lmax)
            ax.plot(np.arange(2, lmx),
                    np.arange(2, lmx) * (np.arange(2, lmx) + 1) / (2*np.pi) * ecl[2:lmx],
                    color=cols[i % 2], label=elab, lw=1.3, ls='-.')
    ax.set_xlabel(r'$\ell$')
    ax.set_ylabel(r'$\ell(\ell+1)C_\ell^{TT}/2\pi$ [$\mu$K$^2$]')
    ax.set_xlim(2, 2500)
    ax.legend(loc='upper right', framealpha=0.8)
    ax.set_title('CMB Temperature Power Spectrum')

    # -- CMB TT ratio --
    ax = axes[1, 0]
    cl_l_safe = np.where(np.abs(cl_l[2:lmax]) > 1e-30, cl_l[2:lmax], 1e-30)
    ratio = cl_m[2:lmax] / cl_l_safe
    ax.plot(ell_plot, ratio, color=COLORS['model'], lw=1.2)
    if extra_cls:
        for i, (ecl, elab) in enumerate(extra_cls):
            lmx = min(len(ecl), lmax)
            r = ecl[2:lmx] / cl_l_safe[:lmx-2]
            ax.plot(np.arange(2, lmx), r, color=cols[i % 2], lw=1.0, ls='-.')
    ax.axhline(1, color='k', lw=0.5, ls=':')
    ax.set_xlabel(r'$\ell$')
    ax.set_ylabel('Ratio to ' + r'$\Lambda$CDM')
    ax.set_xlim(2, 2500)
    ymin, ymax = ratio[np.isfinite(ratio)].min(), ratio[np.isfinite(ratio)].max()
    margin = max(0.01, (ymax - ymin) * 0.3)
    ax.set_ylim(max(0.8, ymin - margin), min(1.2, ymax + margin))

    # -- P(k) --
    ax = axes[0, 1]
    ax.loglog(kh_l, pk_l, color=COLORS['lcdm'], label=r'$\Lambda$CDM', lw=1.2)
    ax.loglog(kh_m, pk_m, color=COLORS['model'],
              label=labels[0] if labels else model_name, lw=1.5, ls='--')
    if extra_pks:
        for i, (ekh, epk, elab) in enumerate(extra_pks):
            ax.loglog(ekh, epk, color=cols[i % 2], label=elab, lw=1.3, ls='-.')
    ax.set_xlabel(r'$k$ [h/Mpc]')
    ax.set_ylabel(r'$P(k)$ [Mpc/h]$^3$')
    ax.legend(loc='lower left', framealpha=0.8)
    ax.set_title('Matter Power Spectrum (z=0)')

    # -- P(k) ratio --
    ax = axes[1, 1]
    # interpolate to common k
    from scipy.interpolate import interp1d
    pk_l_interp = interp1d(np.log(kh_l), np.log(pk_l), fill_value='extrapolate')
    ratio_pk = np.exp(np.log(pk_m) - pk_l_interp(np.log(kh_m)))
    ax.semilogx(kh_m, ratio_pk, color=COLORS['model'], lw=1.2)
    if extra_pks:
        for i, (ekh, epk, elab) in enumerate(extra_pks):
            r = np.exp(np.log(epk) - pk_l_interp(np.log(ekh)))
            ax.semilogx(ekh, r, color=cols[i % 2], lw=1.0, ls='-.')
    ax.axhline(1, color='k', lw=0.5, ls=':')
    ax.set_xlabel(r'$k$ [h/Mpc]')
    ax.set_ylabel('Ratio to ' + r'$\Lambda$CDM')
    ymin, ymax = ratio_pk[np.isfinite(ratio_pk)].min(), ratio_pk[np.isfinite(ratio_pk)].max()
    margin = max(0.01, (ymax - ymin) * 0.3)
    ax.set_ylim(max(0.5, ymin - margin), min(1.5, ymax + margin))

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    return fig


def run_model(setup_fn, model_name, labels=None, multi=False):
    """Run a single model and save figure."""
    print(f"  [{model_name}] computing...")
    try:
        if multi:
            results_list = setup_fn()
            figs_data = []
            for pars, results, lab in results_list:
                cls = results.get_cmb_power_spectra(CMB_unit='muK')['total']
                kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=10, npoints=300)
                figs_data.append((np.arange(cls.shape[0]), cls[:, 0], kh, pk[0], lab))

            # first is primary
            extra_cls = [(d[1], d[4]) for d in figs_data[1:]]
            extra_pks = [(d[2], d[3], d[4]) for d in figs_data[1:]]
            fig = plot_comparison(model_name, figs_data[0][0], figs_data[0][1],
                                  figs_data[0][2], figs_data[0][3],
                                  labels=[figs_data[0][4]],
                                  extra_cls=extra_cls, extra_pks=extra_pks)
        else:
            pars = setup_fn()
            results = camb.get_results(pars)
            cls = results.get_cmb_power_spectra(CMB_unit='muK')['total']
            ell = np.arange(cls.shape[0])
            kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=10, npoints=300)
            fig = plot_comparison(model_name, ell, cls[:, 0], kh, pk[0], labels=labels)

        slug = model_name.replace(" ", "_").replace("-","_").lower()
        out_pdf = os.path.join(FIGDIR, f'{slug}.pdf')
        out_png = os.path.join(FIGDIR, f'{slug}.png')
        fig.savefig(out_pdf, bbox_inches='tight')
        fig.savefig(out_png, bbox_inches='tight', dpi=150)
        plt.close(fig)
        print(f"    -> {out_pdf} + .png")
        return True
    except Exception as e:
        print(f"    [FAILED] {model_name}: {e}")
        import traceback; traceback.print_exc()
        return False


# ============================================================
# MODEL SETUPS
# ============================================================

def setup_dmbaryon():
    """DM-Baryon Scattering"""
    def _multi():
        out = []
        for sigma, lab in [(1e-41, r'$\sigma_0=10^{-41}$ cm$^2$'),
                           (1e-40, r'$\sigma_0=10^{-40}$ cm$^2$')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_matter import DMBaryonScattering
            pars.DarkMatter = DMBaryonScattering()
            pars.DarkMatter.set_params(sigma_dmb=sigma, n_dmb=0, m_dm=1.0)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_dmdr_ethos():
    """DM-Dark Radiation (ETHOS)"""
    def _multi():
        out = []
        for adark, lab in [([1e3], r'$a_{\rm dark}=10^3$'),
                           ([1e5], r'$a_{\rm dark}=10^5$')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_matter import DMDR_ETHOS
            pars.DarkMatter = DMDR_ETHOS()
            pars.DarkMatter.set_params(omdmdrh2=0.005, N_dark=1.0, a_dark=adark)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_decayingdm():
    """Decaying Dark Matter"""
    def _multi():
        out = []
        for gamma, eps, lab in [(100, 0.9, r'$\Gamma=100$, $\epsilon=0.9$'),
                                 (500, 0.8, r'$\Gamma=500$, $\epsilon=0.8$')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_matter import DecayingDM
            pars.DarkMatter = DecayingDM()
            pars.DarkMatter.set_params(Gamma_dcdm=gamma, epsilon_dcdm=eps, f_dcdm=0.1)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_dmneutrino():
    """DM-Neutrino Scattering"""
    def _multi():
        out = []
        for sigma, lab in [(1e-33, r'$\sigma_0=10^{-33}$ cm$^2$'),
                           (1e-31, r'$\sigma_0=10^{-31}$ cm$^2$')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_matter import DMNeutrinoScattering
            pars.DarkMatter = DMNeutrinoScattering()
            pars.DarkMatter.set_params(sigma_dmnu=sigma, n_dmnu=0, m_dm=0.1)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_fuzzydm():
    """Fuzzy Dark Matter"""
    def _multi():
        out = []
        for mass, fax, lab in [(1e-24, 0.3, r'$m=10^{-24}$ eV, $f=0.3$'),
                                (1e-25, 0.3, r'$m=10^{-25}$ eV, $f=0.3$')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_matter import FuzzyDM
            pars.DarkMatter = FuzzyDM()
            pars.DarkMatter.set_params(m_axion=mass, f_axion=fax)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_ethos_murgia():
    """ETHOS Transfer (Murgia 2017 alpha-beta-gamma fit)"""
    def _multi():
        out = []
        for alpha, lab in [(0.02, r'$\alpha=0.02$ Mpc/h'),
                            (0.05, r'$\alpha=0.05$ Mpc/h'),
                            (0.10, r'$\alpha=0.10$ Mpc/h')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_matter import ETHOSTransferMurgia
            pars.DarkMatter = ETHOSTransferMurgia()
            pars.DarkMatter.set_params(alpha_mpch=alpha, beta_mur=2.24, gamma_mur=-4.46)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_ethos_physical():
    """ETHOS Transfer (physical surrogate)"""
    def _multi():
        out = []
        for a_dark, xi, lab in [(1e-4, 0.5, r'$a_{\rm dark}=10^{-4}$, $\xi=0.5$'),
                                 (4e-4, 0.5, r'$a_{\rm dark}=4\times10^{-4}$, $\xi=0.5$'),
                                 (1e-4, 1.0, r'$a_{\rm dark}=10^{-4}$, $\xi=1.0$')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_matter import ETHOSTransferPhysical
            pars.DarkMatter = ETHOSTransferPhysical()
            pars.DarkMatter.set_params(a_dark_n=a_dark, xi_dr=xi, omdmdrh2=0.12, n_dark=2)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_warmdm():
    """Warm Dark Matter"""
    def _multi():
        out = []
        for mwdm, lab in [(3.0, r'$m_{\rm WDM}=3$ keV'),
                           (1.0, r'$m_{\rm WDM}=1$ keV')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_matter import WarmDM
            pars.DarkMatter = WarmDM()
            pars.DarkMatter.set_params(m_wdm=mwdm)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_dmphoton():
    """DM-Photon Scattering"""
    def _multi():
        out = []
        for u, lab in [(1e-4, r'$u=10^{-4}$'),
                       (1e-3, r'$u=10^{-3}$')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_matter import DMPhotonScattering
            pars.DarkMatter = DMPhotonScattering()
            pars.DarkMatter.set_params(u_idm_g=u, n_idm_g=0, m_dm=100.)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_interactingde():
    """Interacting Dark Energy"""
    def _multi():
        out = []
        for xi, itype, lab in [(0.05, 1, r'$\xi=0.05$, type 1'),
                                (0.1, 1, r'$\xi=0.1$, type 1')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_energy import InteractingDE
            pars.DarkEnergy = InteractingDE()
            pars.DarkEnergy.set_params(w=-0.9, xi_ide=xi, interaction_type=itype)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_horndeski():
    """Horndeski Gravity"""
    def _multi():
        out = []
        for aM, lab in [(0.1, r'$\alpha_M=0.1$'),
                         (0.5, r'$\alpha_M=0.5$')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_energy import HorndeskiDE
            pars.DarkEnergy = HorndeskiDE()
            pars.DarkEnergy.set_params(alpha_M=aM, M_star_ini=1.0)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_multi_idm():
    """Multi-Channel Interacting DM"""
    def _multi():
        out = []
        configs = [
            (dict(omdmdrh2=0.005, N_dark=1.0, a_dark=[1e4], m_dm=10.0),
             r'DR only, $a_{\rm dark}=10^4$'),
            (dict(omdmdrh2=0.005, N_dark=1.0, a_dark=[1e4],
                  sigma_dmb=1e-42, n_dmb=0, m_dm=10.0),
             r'DR + DM--baryon'),
        ]
        for kwargs, lab in configs:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_matter import MultiInteractingDM
            pars.DarkMatter = MultiInteractingDM()
            pars.DarkMatter.set_params(**kwargs)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_fuzzydm_field():
    """Fuzzy DM Field (Klein-Gordon)"""
    def _multi():
        out = []
        for mass, fax, lab in [(1e-22, 0.05, r'$m=10^{-22}$ eV, $f=0.05$'),
                                (1e-23, 0.05, r'$m=10^{-23}$ eV, $f=0.05$')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122 * (1.0 - fax))
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            from camb.dark_energy import FuzzyDMField
            pars.DarkEnergy = FuzzyDMField()
            pars.DarkEnergy.set_params(m_axion=mass, f_axion=fax)
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi

def setup_musigma():
    """mu-Sigma Modified Gravity"""
    def _multi():
        out = []
        for mu0, sig0, lab in [(0.3, 0.2, r'$\mu_0=0.3$, $\Sigma_0=0.2$'),
                                (0.5, 0.3, r'$\mu_0=0.5$, $\Sigma_0=0.3$')]:
            pars = camb.CAMBparams()
            pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
            pars.InitPower.set_params(As=2.1e-9, ns=0.965)
            pars.set_for_lmax(2500, lens_potential_accuracy=1)
            pars.WantTransfer = True
            pars.set_matter_power(redshifts=[0], kmax=10.0)
            pars.MG.is_active = True
            pars.MG.mu_0 = mu0
            pars.MG.sigma_0 = sig0
            results = camb.get_results(pars)
            out.append((pars, results, lab))
        return out
    return _multi


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    os.makedirs(FIGDIR, exist_ok=True)
    ensure_lcdm()
    print("LCDM baseline done.\n")

    models = [
        (setup_dmbaryon(), 'DM-Baryon Scattering', True),
        (setup_dmdr_ethos(), 'DM-Dark Radiation (ETHOS)', True),
        (setup_decayingdm(), 'Decaying Dark Matter', True),
        (setup_dmneutrino(), 'DM-Neutrino Scattering', True),
        (setup_fuzzydm(), 'Fuzzy Dark Matter', True),
        (setup_warmdm(), 'Warm Dark Matter', True),
        (setup_ethos_murgia(), 'ETHOS Transfer (Murgia)', True),
        (setup_ethos_physical(), 'ETHOS Transfer (Physical)', True),
        (setup_dmphoton(), 'DM-Photon Scattering', True),
        (setup_interactingde(), 'Interacting Dark Energy', True),
        (setup_horndeski(), 'Horndeski Gravity', True),
        (setup_musigma(), 'Mu-Sigma Modified Gravity', True),
        (setup_multi_idm(), 'Multi-Channel Interacting DM', True),
        (setup_fuzzydm_field(), 'Fuzzy DM Field', True),
    ]

    success, fail = 0, 0
    for setup_fn, name, is_multi in models:
        ok = run_model(setup_fn, name, multi=is_multi)
        if ok:
            success += 1
        else:
            fail += 1

    print(f"\nDone: {success} succeeded, {fail} failed.")
