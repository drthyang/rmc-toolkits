# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Auto StoG range checks on Q and on the low-r enforcement, in every entry point.

Before: ``--enforce-cutoff nan`` exited 0 reporting "hard-set below r = nan A"
while nothing was enforced; ``--enforce-cutoff 1000`` (rmax 50) replaced the
whole RMC G(r) by -<b>^2 with no warning; ``--peak-window 2.5 1.5`` was taken
silently; ``--qmin -3`` shifted the fit, and a NaN qmax passed ScalingConfig's
``qmax <= qmin`` test. ``ScalingConfig`` (JS ``makeConfig``) now requires a
finite qmin >= 0 and a finite qmax, and ``validate_enforcement``
(JS ``validateEnforcement``) a finite cutoff in [0, rmax) and a finite
first-peak window with rmin <= rmax -- the CLI, the API and the page's worker
all call it before computing.
"""

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock
import math
import os
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

from rmc_toolkits import scaling_cli  # noqa: E402
from rmc_toolkits.parsers import write_stog_xy  # noqa: E402
from rmc_toolkits.scaling import ScalingConfig, validate_enforcement  # noqa: E402
from test_scaling_cli import CliSyntheticBase, run_cli  # noqa: E402


def config(**overrides):
    base = dict(qmin=0.6, qmax=29.0, rho0=0.05, b_avg_sq=0.02, rmax=25.0, nr=1000)
    base.update(overrides)
    return ScalingConfig(**base)


class LibraryRangeTests(unittest.TestCase):
    def test_q_range_must_be_finite_and_non_negative(self):
        for overrides, fragment in (
            ({"qmin": math.nan}, "qmin must be finite and >= 0"),
            ({"qmin": -3.0}, "qmin must be finite and >= 0"),
            ({"qmax": math.nan}, "qmax must be finite"),
            ({"qmax": math.inf}, "qmax must be finite"),
        ):
            with self.subTest(**overrides):
                with self.assertRaisesRegex(ValueError, fragment):
                    config(**overrides)
        config(qmin=0.0)  # Q = 0 is allowed

    def test_enforcement_triple(self):
        validate_enforcement(2.48, 2.65, 3.1, rmax=25.0)
        validate_enforcement(0.0, 0.0, 0.0, rmax=25.0)
        for triple, fragment in (
            ((math.nan, 2.0, 2.0), "enforcement cutoff must be finite and >= 0"),
            ((-1.0, -1.0, -1.0), "enforcement cutoff must be finite and >= 0"),
            ((25.0, 25.0, 25.0), "must be below rmax"),
            ((1000.0, 1000.0, 1000.0), "must be below rmax"),
            ((2.0, 2.5, 1.5), "first-peak window must be finite with rmin <= rmax"),
            ((2.0, math.nan, 3.0), "first-peak window must be finite with rmin <= rmax"),
        ):
            with self.subTest(triple=triple):
                with self.assertRaisesRegex(ValueError, fragment):
                    validate_enforcement(*triple, rmax=25.0)


class CliRangeTests(CliSyntheticBase):
    def run_data(self, *flags):
        with TemporaryDirectory() as tmp:
            run = Path(tmp)
            write_stog_xy(run / "synth.dat", self.q, self.sq_meas, title="synthetic")
            with mock.patch.object(scaling_cli, "autoscale", wraps=scaling_cli.autoscale) as auto:
                code, out, err = run_cli([
                    "--data", run / "synth.dat", "--qmin", "0.6", "--qmax", "29", "--rho0", "0.05",
                    "--b-avg-sq", "0.02", "--rmax", "25", "--nr", "1000", *flags,
                ])
            written = sorted(p.name for p in (run / "autoscale").glob("*")) if (run / "autoscale").exists() else []
        return code, err, written, auto.called

    def test_bad_enforcement_and_q_flags_are_refused_before_computing(self):
        for flags, fragment in (
            (["--enforce-cutoff", "nan"], "enforcement cutoff must be finite and >= 0"),
            (["--enforce-cutoff", "1000"], "must be below rmax"),
            (["--enforce-cutoff", "2.0", "--peak-window", "2.5", "1.5"], "rmin <= rmax"),
            (["--qmin=-3"], "qmin must be finite and >= 0"),
        ):
            with self.subTest(flags=flags):
                code, err, written, computed = self.run_data(*flags)
                self.assertEqual(code, 2, err)
                self.assertIn(fragment, err)
                self.assertEqual(written, [])
                self.assertFalse(computed)

    def test_a_valid_explicit_cutoff_still_runs(self):
        code, err, written, _ = self.run_data("--enforce-cutoff", "2.0")
        self.assertEqual(code, 0, err)
        self.assertEqual(len(written), 9)


class ApiRangeTests(CliSyntheticBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        import app as backend_app

        backend_app.app.config.update(TESTING=True)
        cls.backend = backend_app
        cls.client = backend_app.app.test_client()
        cls.root = ROOT / "results" / "scaling_range_validation"
        shutil.rmtree(cls.root, ignore_errors=True)
        cls.root.mkdir(parents=True)
        write_stog_xy(cls.root / "synth.dat", cls.q, cls.sq_meas, title="synthetic")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.root, ignore_errors=True)

    def test_bad_ranges_are_400_before_computing(self):
        base = {"path": str((self.root / "synth.dat").relative_to(ROOT)), "qmin": 0.6, "qmax": 29,
                "rho0": 0.05, "bAvgSq": 0.02, "rmax": 25, "nr": 1000}
        for extra, fragment in (
            ({"enforceCutoff": 1000}, "must be below rmax"),
            ({"enforceCutoff": -1}, "enforcement cutoff must be finite and >= 0"),
            ({"enforceCutoff": 2.0, "peakWindow": [2.5, 1.5]}, "rmin <= rmax"),
            ({"qmin": -3}, "qmin must be finite and >= 0"),
        ):
            with self.subTest(extra=extra):
                with mock.patch.object(self.backend, "_compute_scaling") as compute:
                    response = self.client.post("/api/scaling/preview", json={**base, **extra})
                self.assertEqual(response.status_code, 400, response.get_data(as_text=True)[:300])
                self.assertIn(fragment, response.get_json()["error"])
                self.assertFalse(compute.called)


if __name__ == "__main__":
    unittest.main()
