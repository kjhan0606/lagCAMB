import math

import numpy as np
import pytest

import camb
from camb.baseconfig import CAMBError
from camb.dark_matter import TransientGDM


def _power(model=None, maxkh=120.0):
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=68.4, ombh2=0.0225, omch2=0.118, mnu=0.06, tau=0.0575
    )
    pars.InitPower.set_params(As=math.exp(3.05) * 1.0e-10, ns=0.971)
    if model is not None:
        pars.DarkMatter = model
    pars.WantCls = False
    pars.set_matter_power(redshifts=[25.0, 12.0], kmax=maxkh * 1.05)
    results = camb.get_results(pars)
    return results.get_matter_power_spectrum(
        minkh=1.0e-3, maxkh=maxkh, npoints=240
    )


def _transfers(model=None, redshift=49.0, maxkh=60.0):
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=68.4, ombh2=0.0225, omch2=0.118, mnu=0.06, tau=0.0575
    )
    pars.InitPower.set_params(As=math.exp(3.05) * 1.0e-10, ns=0.971)
    if model is not None:
        pars.DarkMatter = model
    pars.WantCls = False
    pars.set_matter_power(redshifts=[redshift], kmax=maxkh)
    return camb.get_results(pars).get_matter_transfer_data().transfer_data


def test_transient_gdm_parameter_contract():
    with pytest.raises(CAMBError):
        TransientGDM().set_params(f_X=1.0)
    with pytest.raises(CAMBError):
        TransientGDM().set_params(w_p=-1.0)
    with pytest.raises(CAMBError):
        TransientGDM().set_params(a_i=1.0e-5, a_f=1.0e-6)


def test_transient_gdm_zero_fraction_is_native_cdm():
    k_native, z_native, p_native = _power()
    k_zero, z_zero, p_zero = _power(
        TransientGDM().set_params(f_X=0.0, w_p=-0.45)
    )
    assert np.array_equal(k_native, k_zero)
    assert np.array_equal(z_native, z_zero)
    assert np.array_equal(p_native, p_zero)


def test_transient_gdm_velocity_transfer_contract():
    assert camb.model.transfer_names[camb.model.Transfer_gdm_velocity - 1] == "v_newtonian_gdm"

    native = _transfers()
    zero = _transfers(TransientGDM().set_params(f_X=0.0, w_p=-0.45))
    ic_columns = [
        camb.model.Transfer_kh,
        camb.model.Transfer_cdm,
        camb.model.Transfer_b,
        camb.model.Transfer_tot,
        camb.model.Transfer_nonu,
        camb.model.Transfer_Newt_vel_cdm,
        camb.model.Transfer_Newt_vel_baryon,
        camb.model.Transfer_gdm_velocity,
    ]
    for column in ic_columns:
        assert np.array_equal(native[column - 1], zero[column - 1])
    # Three single-precision massive-neutrino entries can differ at roundoff
    # level when the otherwise inert DM object is allocated.
    assert np.allclose(native, zero, rtol=2.0e-7, atol=5.0e-13)

    set_ii = _transfers(
        TransientGDM().set_params(f_X=0.00155, w_p=-0.45, kappa0_hmpc=6.4)
    )
    v_total_dm = set_ii[camb.model.Transfer_Newt_vel_cdm - 1, :, 0]
    v_x = set_ii[camb.model.Transfer_gdm_velocity - 1, :, 0]
    delta_c = set_ii[camb.model.Transfer_cdm - 1, :, 0]
    delta_b = set_ii[camb.model.Transfer_b - 1, :, 0]
    delta_x = set_ii[camb.model.Transfer_dm_dr - 1, :, 0]
    delta_nonu = set_ii[camb.model.Transfer_nonu - 1, :, 0]
    f_x = 0.00155
    reconstructed_nonu = (
        0.0225 * delta_b
        + 0.118 * ((1.0 - f_x) * delta_c + f_x * delta_x)
    ) / (0.0225 + 0.118)
    assert np.all(np.isfinite(v_x))
    assert np.any(v_x != 0.0)
    assert np.any(v_total_dm != v_x)
    assert np.allclose(delta_nonu, reconstructed_nonu, rtol=3.0e-7, atol=0.0)


def test_transient_gdm_benchmark_scale_ordering_and_saturation():
    k, z, p_native = _power()
    k_i, z_i, p_i = _power(
        TransientGDM().set_params(f_X=0.00161, w_p=-0.99, kappa0_hmpc=3.6)
    )
    k_ii, z_ii, p_ii = _power(
        TransientGDM().set_params(f_X=0.00155, w_p=-0.45, kappa0_hmpc=6.4)
    )
    assert np.array_equal(k, k_i) and np.array_equal(k, k_ii)
    assert np.array_equal(z, z_i) and np.array_equal(z, z_ii)
    ratio_i = p_i / p_native
    ratio_ii = p_ii / p_native
    j10 = int(np.argmin(np.abs(k - 10.0)))
    j30 = int(np.argmin(np.abs(k - 30.0)))
    j100 = int(np.argmin(np.abs(k - 100.0)))
    assert np.all(ratio_i[:, j30] > ratio_i[:, j10])
    assert np.all(ratio_i[:, j10] > 1.0)
    assert np.all(ratio_ii[:, j30] > ratio_i[:, j30])
    assert np.all(ratio_ii[:, j100] > ratio_ii[:, j30])
    assert np.all(
        ratio_ii[:, j100] - ratio_ii[:, j30]
        < 2.0 * (ratio_ii[:, j30] - ratio_ii[:, j10])
    )
