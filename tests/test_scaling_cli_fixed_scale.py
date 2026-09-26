# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""rmc-autoscale fixed scaling: --scale/--offset must be finite, and the scale non-zero.

``--scale nan`` (or ``--offset nan``) exited 0 and wrote all nine files with
100 % NaN RMCProfile inputs, and ``--scale 0`` wrote a family from a discarded
data set. read_stog_inp already refuses a non-finite or zero yscale, and the
API a non-finite or zero ``a``; the flags bypassed both.
"""

from pathlib import Path
from tempfile import TemporaryDirectory
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from test_scaling_cli import CliSyntheticBase, run_cli  # noqa: E402
from rmc_toolkits.parsers import write_stog_xy  # noqa: E402


class FixedScaleFlagTests(CliSyntheticBase):
    def run_data(self, *flags):
        with TemporaryDirectory() as tmp:
            run = Path(tmp)
            write_stog_xy(run / "synth.dat", self.q, self.sq_meas, title="synthetic")
            code, out, err = run_cli([
                "--data", run / "synth.dat", "--qmin", "0.6", "--qmax", "29", "--rho0", "0.05",
                "--b-avg-sq", "0.02", *flags,
            ])
            written = sorted(path.name for path in (run / "autoscale").glob("*")) if (run / "autoscale").exists() else []
        return code, out, err, written

    def test_non_finite_or_zero_scale_and_offset_are_refused(self):
        for flags, fragment in (
            (["--scale", "nan"], "--scale must be a finite, non-zero number"),
            (["--scale", "inf"], "--scale must be a finite, non-zero number"),
            (["--scale", "0"], "--scale must be a finite, non-zero number"),
            (["--scale", "1", "--offset", "nan"], "--offset must be a finite number"),
            (["--scale", "1", "--offset=-inf"], "--offset must be a finite number"),
        ):
            with self.subTest(flags=flags):
                code, _, err, written = self.run_data(*flags)
                self.assertEqual(code, 2, err)
                self.assertIn(fragment, err)
                self.assertEqual(written, [])

    def test_finite_fixed_scaling_still_runs(self):
        code, _, err, written = self.run_data("--scale", "10", "--offset", "-9")
        self.assertEqual(code, 0, err)
        self.assertEqual(len(written), 9)


if __name__ == "__main__":
    unittest.main()
