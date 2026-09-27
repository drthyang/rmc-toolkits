# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""estimate_rho0 never accepts a physically impossible root (0.6.0 audit, stog-b).

On missing-low-Q data the density-limit amplitude a_density(rho0) can cross the
rho0-independent Faber-Ziman amplitude a second time at a density no solid has:
with the closest approach pinned (r0 = 2.6 A, as a stog.inp or MINIMUM_DISTANCES
header supplies it) the Mn3Sn 300 K run "converged" at 0.428 A^-3 (~50 g/cm^3,
6.8x the real 0.063) from the true density as seed, and the CLI adopted it. The
iterate is now confined to RHO0_PHYSICAL_RANGE and a concordant root must also
satisfy the density limit; otherwise converged=False with a reason, and the CLI
refuses the estimate.
"""

import contextlib
from dataclasses import replace
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from rmc_toolkits.parsers import read_stog_xy
from rmc_toolkits.scaling import (
    RHO0_PHYSICAL_RANGE,
    ScalingConfig,
    diagnostics_summary,
    estimate_rho0,
    scale_pipeline,
)
from rmc_toolkits.scaling_cli import main
from rmc_toolkits.scattering import faber_ziman
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

ROOT = Path(__file__).resolve().parents[1]
MN3SN = ROOT / "data" / "stog_tests"
MN3SN_RUNS = {
    "stog": "PG3_55537_rebin.sq",
    "stog_300K": "PG3_55526_SQ_rebin.sq",
    "stog_500K": "PG3_54139_SQ_rebin.dat",
    "stog_59438": "PG3_59438_SQ_rebin.dat",
}
HAVE_MN3SN = all((MN3SN / name / raw).exists() for name, raw in MN3SN_RUNS.items())
RHO0, B2 = 0.05, 0.02


def synthetic_sq():
    """Repo synthetic model (rho0 = 0.05) measured with a = 5, b = -4."""
    q = np.arange(20, 981) * 0.03
    r = np.arange(1, 12001) * 0.005
    g = 0.5 * (1.0 + np.tanh((r - 2.65) / 0.07)) + 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    sq_true = fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, RHO0), q))
    head = q <= q[0] + 1.0
    _, s_true_0 = np.polyfit(q[head], sq_true[head], 1)
    return q, (sq_true + 4.0) / 5.0, B2 * (1.0 - float(s_true_0))


def synthetic_config(**overrides):
    base = dict(
        qmin=0.6, qmax=30.0, rho0=0.02, b_avg_sq=B2, r0=2.5, r_fit_min=1.2,
        r_fit_max=2.25, rmax=25.0, nr=1000,
    )
    base.update(overrides)
    return ScalingConfig(**base)


class PhysicalRangeTests(unittest.TestCase):
    def test_default_range(self):
        self.assertEqual(RHO0_PHYSICAL_RANGE, (0.005, 0.25))

    def test_genuine_root_still_converges(self):
        q, sq, b_sq_avg = synthetic_sq()
        est = estimate_rho0(q, sq, synthetic_config(b_sq_avg=b_sq_avg))
        self.assertTrue(est["converged"])
        self.assertIsNone(est["reason"])
        self.assertLess(abs(est["rho0"] - RHO0) / RHO0, 0.05)

    def test_a_step_out_of_the_range_stops_with_a_reason(self):
        q, sq, b_sq_avg = synthetic_sq()
        est = estimate_rho0(q, sq, synthetic_config(b_sq_avg=b_sq_avg), rho_max=0.03)
        self.assertFalse(est["converged"])
        self.assertIn("physical density range", est["reason"])
        self.assertTrue(all(row[0] <= 0.03 for row in est["history"]))

    def test_the_range_is_validated(self):
        q, sq, b_sq_avg = synthetic_sq()
        with self.assertRaisesRegex(ValueError, "rho_min < rho_max"):
            estimate_rho0(q, sq, synthetic_config(b_sq_avg=b_sq_avg), rho_min=0.3)

    def test_concordant_root_that_fails_the_density_limit_is_spurious(self):
        # A pass whose amplitudes agree (a_fz = a) on a fit that leaves the
        # low-r window far from g = 0: concordance alone must not be accepted.
        q, sq, b_sq_avg = synthetic_sq()
        config = synthetic_config(b_sq_avg=b_sq_avg)
        # At the seed rho0 = 0.02 the density-limit scale is ~2; a = 1 leaves mean g ~0.4.
        bad = scale_pipeline(q, sq, config, 1.0, 0.0)
        self.assertFalse(diagnostics_summary(bad, config)["density_limit_satisfied"])
        bad = replace(bad, a_fz=bad.a)
        with mock.patch("rmc_toolkits.scaling.autoscale", return_value=bad):
            est = estimate_rho0(q, sq, config)
        self.assertFalse(est["converged"])
        self.assertIn("spurious root", est["reason"])


@unittest.skipUnless(HAVE_MN3SN, "Mn3Sn PG3 neutron runs not present")
class Mn3SnSpuriousRootTests(unittest.TestCase):
    """Real missing-low-Q data (Qmin 0.82): the density limit fails at every density."""

    @classmethod
    def setUpClass(cls):
        cls.fz = faber_ziman("Mn3Sn")
        cls.data = {
            name: read_stog_xy(MN3SN / name / raw)[:2] for name, raw in MN3SN_RUNS.items()
        }

    def config(self, seed, **overrides):
        return ScalingConfig(
            qmin=0.82, qmax=28.0, rho0=seed, b_avg_sq=self.fz.b_avg_sq_barn,
            b_sq_avg=self.fz.b_sq_avg_barn, **overrides,
        )

    def check(self, name, seed, **overrides):
        q, sq = self.data[name]
        try:
            est = estimate_rho0(q, sq, self.config(seed, **overrides))
        except ValueError:
            return  # the seed density itself cannot be fitted: nothing adopted
        self.assertFalse(est["converged"], est)
        self.assertTrue(est["reason"])
        low, high = RHO0_PHYSICAL_RANGE
        self.assertTrue(all(low <= row[0] <= high for row in est["history"]), est["history"])

    def test_pinned_closest_approach_never_adopts_an_impossible_density(self):
        for name in MN3SN_RUNS:
            for seed in (0.02, 0.05, 0.063049):
                with self.subTest(run=name, seed=seed):
                    self.check(name, seed, r0=2.6)

    def test_detected_window_never_adopts_an_impossible_density(self):
        for seed in (0.02, 0.05, 0.063049):
            with self.subTest(seed=seed):
                self.check("stog_300K", seed)

    def test_cli_refuses_the_spurious_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = main([
                    "--data", str(MN3SN / "stog_300K" / MN3SN_RUNS["stog_300K"]),
                    "--qmin", "0.82", "--qmax", "28", "--formula", "Mn3Sn",
                    "--rho0", "0.063049", "--r0", "2.6", "--estimate-rho0",
                    "--out-dir", str(Path(tmp) / "out"),
                ])
        self.assertEqual(code, 2, out.getvalue())
        self.assertIn("rho0 self-consistency did not converge", err.getvalue())
        self.assertNotIn("rho0 self-consistency:", out.getvalue())


if __name__ == "__main__":
    unittest.main()
