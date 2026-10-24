"""Native CLI access to the existing InteractingDE implementation."""
import math
from pathlib import Path
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class InteractingDECLITest(unittest.TestCase):
    def calculate(self, folder, model, w=-1, xi=0):
        text = (ROOT / "inifiles/params.ini").read_text()
        values = {
            "output_root": "cli", "get_scalar_cls": "F", "get_transfer": "T", "do_lensing": "F",
            "dark_energy_model": model, "w": str(w), "wa": "0", "xi_ide": str(xi), "interaction_type": "1",
            "transfer_num_redshifts": "1", "transfer_redshift(1)": "0", "transfer_filename(1)": "transfer.dat",
            "transfer_matterpower(1)": "power.dat", "transfer_kmax": "1", "number_of_threads": "1",
        }
        lines = []
        for line in text.splitlines():
            key = line.split("=", 1)[0].strip() if "=" in line and not line.lstrip().startswith("#") else ""
            if key in values:
                lines.append(key + " = " + values.pop(key))
            else:
                lines.append(line)
        lines.extend(key + " = " + value for key, value in values.items())
        folder.mkdir()
        (folder / "params.ini").write_text("\n".join(lines) + "\n")
        run = subprocess.run([str(ROOT / "fortran/camb"), "params.ini"], cwd=folder, capture_output=True, text=True)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        lines = (folder / "cli_transfer.dat").read_text().splitlines()
        rows = [[float(value) for value in line.split()] for line in lines
                if line.strip() and not line.startswith("#")]
        self.assertTrue(all(math.isfinite(v) for row in rows for v in row))
        self.assertTrue(all(row[6] > 0 for row in rows))
        return rows

    def test_zero_coupling_matches_lcdm_and_nonzero_cli_is_finite(self):
        if not (ROOT / "fortran/camb").is_file():
            self.skipTest("Build the native CLI with make -j1 camb first")
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            standard = self.calculate(root / "lcdm", "fluid")
            zero = self.calculate(root / "zero", "InteractingDE")
            self.assertEqual(len(standard), len(zero))
            self.assertLess(max(abs(a[6]/b[6]-1) for a, b in zip(standard, zero)), 1e-6)
            self.calculate(root / "coupled", "InteractingDE", w=-0.9, xi=-0.01)


if __name__ == "__main__":
    unittest.main()
