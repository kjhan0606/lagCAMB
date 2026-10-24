"""Independent CAMB-iDEDM benchmark for the type-1 local-expansion fluid."""
import json
from pathlib import Path
import unittest

import numpy as np
import camb
from camb.dark_energy import InteractingDE
from camb.baseconfig import CAMBError


class IDEReferenceTest(unittest.TestCase):
    def test_nonzero_coupling_against_independent_reference(self):
        fixture = json.loads(Path(__file__).with_name("ide_reference_power.json").read_text())
        for case in fixture["cases"]:
            with self.subTest(case=case["name"]):
                p = camb.CAMBparams()
                p.set_cosmology(**fixture["settings"]["common"])
                p.WantCls = False
                p.InitPower.set_params(As=2.1e-9, ns=0.966, pivot_scalar=0.05)
                p.DarkEnergy = InteractingDE()
                p.DarkEnergy.set_params(w=case["w"], xi_ide=case["xi"], interaction_type=1, cs2_ide=1)
                p.set_accuracy(AccuracyBoost=2)
                p.set_matter_power(redshifts=[99, 1, 0], kmax=3)
                result = camb.get_results(p)
                kh, zs, power = result.get_matter_power_spectrum(minkh=0.001, maxkh=1, npoints=160)
                indices = fixture["indices"]
                np.testing.assert_allclose(kh[indices], case["k"], rtol=1e-12)
                np.testing.assert_array_equal(zs, [0, 1, 99])
                # 0.02% allows old/new CAMB baseline differences, but detects the
                # previous 0.13--0.66% source-term error without fitting A_s.
                np.testing.assert_allclose(power[:, indices], case["power"], rtol=2e-4)

    def test_nonzero_coupling_vacuum_fluid_is_rejected(self):
        p = camb.CAMBparams()
        p.set_cosmology(H0=67, ombh2=0.022, omch2=0.12)
        p.DarkEnergy = InteractingDE()
        p.DarkEnergy.set_params(w=-1, xi_ide=-0.01)
        with self.assertRaisesRegex(CAMBError, "Type-1 IDE fluid requires"):
            camb.get_background(p)


if __name__ == "__main__":
    unittest.main()
