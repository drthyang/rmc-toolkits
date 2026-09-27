# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Despike runs exactly once in the auto path (0.6.0 audit, stog-b group).

``_autoscale_pass`` cropped + despiked the data, fitted (a, b) on it, then handed
the already-despiked arrays to ``scale_pipeline``, which cropped + despiked them a
SECOND time: the written files and fit diagnostics came from a smaller point set
than the one fitted (2161 vs 2318 points on the Mn3Sn 59438 run), and
``n_despiked`` counted only the second pass (157 of 539 removed). The JS engine
despikes once. Now both engines fit, write and report the same single pass.
"""

from dataclasses import replace
from pathlib import Path
import unittest

import numpy as np

from rmc_toolkits.parsers import read_stog_xy
from rmc_toolkits.scaling import ScalingConfig, autoscale, crop_sq, diagnostics_summary
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

ROOT = Path(__file__).resolve().parents[1]
RUN_59438 = ROOT / "data" / "stog_tests" / "stog_59438" / "PG3_59438_SQ_rebin.dat"
RHO0, B2 = 0.05, 0.02


def synthetic_g(r):
    """The shared synthetic model (tests/generate_autoscale_fixture.py)."""
    onset = 0.5 * (1.0 + np.tanh((r - 2.65) / 0.07))
    peak = 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    return onset + peak


def glitchy_sq():
    """Model S(Q) scaled by a = 10, b = -9, with noise and 12 tail glitches."""
    q = np.arange(20, 981) * 0.03
    r = np.arange(1, 12001) * 0.005
    sq_true = fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, synthetic_g(r), RHO0), q))
    rng = np.random.default_rng(3)
    sq = (sq_true + 9.0) / 10.0 + rng.normal(0.0, 2e-3, q.size)
    sq[rng.choice(np.where(q > 20)[0], 12, replace=False)] += 0.3
    return q, sq


class SingleDespikeTests(unittest.TestCase):
    def check_single_pass(self, q, sq, config):
        fitted_q, _, _ = crop_sq(q, sq, config)  # one despike pass
        raw_q, _, _ = crop_sq(q, sq, replace(config, despike=False))
        result = autoscale(q, sq, config)
        np.testing.assert_array_equal(result.q, fitted_q)
        self.assertEqual(result.provenance["n_q_points"], fitted_q.size)
        self.assertEqual(result.provenance["n_despiked"], raw_q.size - fitted_q.size)
        return result, raw_q.size - fitted_q.size

    def test_outputs_and_count_come_from_the_fitted_single_pass(self):
        q, sq = glitchy_sq()
        config = ScalingConfig(
            qmin=0.6, qmax=30.0, rho0=RHO0, b_avg_sq=B2, r0=2.5, rmax=25.0, nr=1000,
            despike=True,
        )
        result, removed = self.check_single_pass(q, sq, config)
        self.assertEqual(removed, 12)  # exactly the glitches
        self.assertLess(abs(result.a - 10.0) / 10.0, 0.02)

    def test_fz_mode_despikes_once_too(self):
        q, sq = glitchy_sq()
        config = ScalingConfig(
            qmin=0.6, qmax=30.0, rho0=RHO0, b_avg_sq=B2, b_sq_avg=0.0348, r0=2.5,
            rmax=25.0, nr=1000, despike=True, amplitude_criterion="fz",
        )
        self.check_single_pass(q, sq, config)

    @unittest.skipUnless(RUN_59438.exists(), "Mn3Sn 59438 sample data not present")
    def test_mn3sn_59438_reports_every_removed_point(self):
        q, sq = read_stog_xy(RUN_59438)[:2]
        config = ScalingConfig(
            qmin=1.0, qmax=28.0, rho0=0.063049, b_avg_sq=0.015407, r0=2.7,
            rmax=20.0, nr=2000, despike=True,
        )
        result, removed = self.check_single_pass(q, sq, config)
        # 2701 finite points in the window, 393 flagged by the single pass (the
        # second pass removed 172 more and n_despiked reported only those).
        self.assertEqual(removed, 393)
        self.assertEqual(result.q.size, 2308)
        summary = diagnostics_summary(result, config)
        self.assertEqual(summary["a"], result.a)


if __name__ == "__main__":
    unittest.main()
