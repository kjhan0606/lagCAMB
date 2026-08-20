Dark Matter models
==================================

This module provides several non-standard dark-matter models on top of the
standard CAMB CDM. Each subsection below shows the comparison of CMB TT and
matter :math:`P(k)` against :math:`\Lambda`\ CDM for representative parameter
choices (left/top: spectra; right/bottom: ratio to :math:`\Lambda`\ CDM).

.. autoclass:: camb.dark_matter.DarkMatterModel
   :members:

DM--Baryon Scattering
---------------------

.. autoclass:: camb.dark_matter.DMBaryonScattering
   :show-inheritance:
   :members:

.. image:: figures/dm_baryon_scattering.png
   :alt: DM-Baryon scattering vs LCDM
   :width: 100%

DM--Dark Radiation (ETHOS)
--------------------------

.. autoclass:: camb.dark_matter.DMDR_ETHOS
   :show-inheritance:
   :members:

.. image:: figures/dm_dark_radiation_ethos.png
   :alt: DMDR ETHOS vs LCDM
   :width: 100%

Decaying Dark Matter
--------------------

.. autoclass:: camb.dark_matter.DecayingDM
   :show-inheritance:
   :members:

.. image:: figures/decaying_dark_matter.png
   :alt: Decaying DM vs LCDM
   :width: 100%

DM--Neutrino Scattering
-----------------------

.. autoclass:: camb.dark_matter.DMNeutrinoScattering
   :show-inheritance:
   :members:

.. image:: figures/dm_neutrino_scattering.png
   :alt: DM-Neutrino scattering vs LCDM
   :width: 100%

Warm Dark Matter
----------------

.. autoclass:: camb.dark_matter.WarmDM
   :show-inheritance:
   :members:

.. image:: figures/warm_dark_matter.png
   :alt: WarmDM vs LCDM
   :width: 100%

Fuzzy Dark Matter
-----------------

``FuzzyDM`` evolves the quadratic Klein--Gordon background and two exact field
perturbations through the frozen and early oscillatory regimes.  At the common
physical surface ``m/H=match_ratio``, a Passaglia--Hu phase-envelope map
initializes density and heat-flux EFA variables; exact and EFA stress-energy
are then joined with a conservative C1 window.  Log-space shooting enforces
the requested present-day axion fraction of the total ``omch2`` budget.

Active models enforce ``1e-23 <= m_axion/eV <= 1e-21``,
``1e-3 <= f_axion <= 0.1``, ``50 <= match_ratio <= 75``, and a Lambda
dark-energy background.  The match-ratio scan changes the tested z=0 spectra
by at most ``5e-5`` over the validated k ranges.  At ``f_axion=0.05``, the
current AxiCLASS double-ratio residuals are about 1.15%, 1.11%, and 0.87% for
``m_axion=1e-23, 1e-22, 1e-21 eV`` over ``k <= 10, 30, 50 h/Mpc``.  These are
cross-code residuals, not an accuracy guarantee over the full abundance guard
range.  The older DarkEnergy-slot ``FuzzyDMField`` remains disabled.
The DarkMatter-slot result is conditionally validated only for the reported
``z=0``, ``f_axion=0.05`` linear-power benchmarks.  The abundance interval is
a guard rather than a precision domain, and CMB/lensing outputs have no
independent exact-solver comparison yet.

.. autoclass:: camb.dark_matter.FuzzyDM
   :show-inheritance:
   :members:

.. image:: figures/fuzzy_dark_matter.png
   :alt: FuzzyDM vs LCDM
   :width: 100%

DM--Photon Scattering
---------------------

.. autoclass:: camb.dark_matter.DMPhotonScattering
   :show-inheritance:
   :members:

.. image:: figures/dm_photon_scattering.png
   :alt: DM-Photon scattering vs LCDM
   :width: 100%

Multi-Channel Interacting DM
----------------------------

.. autoclass:: camb.dark_matter.MultiInteractingDM
   :show-inheritance:
   :members:

.. image:: figures/multi_channel_interacting_dm.png
   :alt: Multi-Channel IDM vs LCDM
   :width: 100%

ETHOS Transfer (Murgia fit)
---------------------------

Analytic surrogate for ETHOS small-scale suppression using the
Murgia et al. (2017) :math:`T(k) = [1 + (\alpha k)^\beta]^\gamma`
fitting form, applied directly to the matter power spectrum
(no CMB modification). Defaults: :math:`\beta=2.24,\ \gamma=-4.46`.

.. autoclass:: camb.dark_matter.ETHOSTransferMurgia
   :show-inheritance:
   :members:

.. image:: figures/ethos_transfer_murgia.png
   :alt: ETHOS Transfer (Murgia) vs LCDM
   :width: 100%

ETHOS Transfer (Physical surrogate)
-----------------------------------

Same Murgia :math:`T(k)` form, but with the suppression scale
:math:`\alpha` derived from physical ETHOS parameters
(:math:`a_{\rm dark,n}`, :math:`\xi_{\rm DR}`, :math:`\omega_{\rm dmdr}h^2`)
through a power-law fit whose exponents are user-tunable. Intended for
fast parameter scans; for the full self-consistent ETHOS Boltzmann
solution (DAO + CMB effects) use ``DMDR_ETHOS``.

.. autoclass:: camb.dark_matter.ETHOSTransferPhysical
   :show-inheritance:
   :members:

.. image:: figures/ethos_transfer_physical.png
   :alt: ETHOS Transfer (Physical) vs LCDM
   :width: 100%
