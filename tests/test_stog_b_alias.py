# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The r range the discrete sine transform resolves: r < pi/dQ (1.0 audit, stog-b).

The docs said the aliasing period is 2 pi/dQ and called r_max = 50 A safe at
dQ = 0.1. On a uniform grid G(2 pi/dQ - r) = -G(r): beyond pi/dQ the output is a
negated mirror image (a shell at 20 A reappears inverted at 42.8 A). The engines
now report r_alias_limit = pi/dQ and flag r_max beyond it; the CLI warns.
"""

import contextlib
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import write_stog_xy
from rmc_toolkits.scaling import ScalingConfig, alias_limit, diagnostics_summary, scale_pipeline
from rmc_toolkits.scaling_cli import main
from rmc_toolkits.transforms import fq_to_gpdf, fq_to_sq, g_to_gpdf, gpdf_to_fq


def coarse_model(dq=0.1):
    q = np.arange(6, 295) * dq  # 0.6 .. 29.4 A^-1
    r = np.arange(1, 12001) * 0.005
    g = 0.5 * (1 + np.tanh((r - 2.65) / 0.07)) + 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    return q, fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, 0.05), q))


class AliasTests(unittest.TestCase):
    def test_transform_folds_at_pi_over_dq(self):
        q = np.arange(5, 301) * 0.1
        fq = np.exp(-0.5 * (0.1 * q) ** 2) * np.sin(20.0 * q)  # one shell at 20 A
        r = np.arange(1, 5001) * 0.01
        g = fq_to_gpdf(q, fq, r)
        ghost = int(np.argmin(g))
        self.assertAlmostEqual(r[ghost], 2 * np.pi / 0.1 - 20.0, places=1)  # 42.83 < 50
        self.assertGreater(-g[ghost], 0.99 * g.max())  # full-amplitude, inverted
        self.assertAlmostEqual(alias_limit(q), np.pi / 0.1, places=9)

    def test_summary_flags_rmax_beyond_the_limit(self):
        q, sq = coarse_model()
        for rmax, beyond in ((50.0, True), (25.0, False)):
            config = ScalingConfig(qmin=0.6, qmax=29.4, rho0=0.05, b_avg_sq=0.02, rmax=rmax, nr=int(rmax * 100))
            summary = diagnostics_summary(scale_pipeline(q, sq, config, 1.0, 0.0), config)
            with self.subTest(rmax=rmax):
                self.assertAlmostEqual(summary["r_alias_limit"], np.pi / 0.1, places=6)
                self.assertIs(summary["rmax_beyond_alias_limit"], beyond)

    def test_cli_warns(self):
        q, sq = coarse_model()
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "coarse.sq"
            write_stog_xy(data, q, sq)
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = main([
                    "--data", str(data), "--qmin", "0.6", "--qmax", "29.4", "--rho0", "0.05",
                    "--b-avg-sq", "0.02", "--manual", "--scale", "1", "--offset", "0",
                    "--out-dir", str(Path(tmp) / "out"),
                ])
        self.assertEqual(code, 0, err.getvalue())
        self.assertIn("exceeds the aliasing limit pi/dQ = 31.42 A", out.getvalue())


def log_binned_shells():
    """A log-binned grid (dQ/Q = 0.004, Q 0.5-30) and F(Q) with shells at 2.5 and 20 A."""
    q = 0.5 * 1.004 ** np.arange(1030)
    q = q[q <= 30.0]
    fq = np.exp(-0.5 * (0.08 * q) ** 2) * np.sin(2.5 * q) + 0.3 * np.exp(
        -0.5 * (0.05 * q) ** 2
    ) * np.sin(20.0 * q)
    return q, fq


class NonUniformGridTests(unittest.TestCase):
    """The limit is set by the coarsest step, not the median (1.0 review).

    pi/median(dQ) = 203 A on the log-binned grid, while its coarse high-Q
    steps (dQ = 0.119) already fail beyond pi/dQ_max = 26.4 A: against a
    0.001-spaced transform the error jumps from < 0.004 below 25 A to 0.07-0.09
    at 30-50 A, where the true |G| <= 0.016. r_max = 50 was not flagged.
    """

    def test_log_binned_grid_limit_is_pi_over_the_coarsest_step(self):
        q, fq = log_binned_shells()
        limit = alias_limit(q)
        self.assertAlmostEqual(limit, np.pi / np.diff(q).max(), places=9)
        self.assertLess(limit, 27.0)
        r = np.arange(1, 5001) * 0.01
        q_ref = np.arange(0.5, 30.0005, 0.001)
        fq_ref = np.exp(-0.5 * (0.08 * q_ref) ** 2) * np.sin(2.5 * q_ref) + 0.3 * np.exp(
            -0.5 * (0.05 * q_ref) ** 2
        ) * np.sin(20.0 * q_ref)
        error = np.abs(fq_to_gpdf(q, fq, r) - fq_to_gpdf(q_ref, fq_ref, r))
        self.assertLess(error[r < 0.95 * limit].max(), 0.005)
        self.assertGreater(error[(r > 30) & (r <= 50)].max(), 0.05)

    def test_summary_flags_rmax_50_on_the_log_grid(self):
        q, _ = log_binned_shells()
        _, sq = coarse_model()
        sq = np.interp(q, np.arange(6, 295) * 0.1, sq)
        config = ScalingConfig(qmin=0.5, qmax=30.0, rho0=0.05, b_avg_sq=0.02, rmax=50.0, nr=5000)
        summary = diagnostics_summary(scale_pipeline(q, sq, config, 1.0, 0.0), config)
        self.assertLess(summary["r_alias_limit"], 27.0)
        self.assertTrue(summary["rmax_beyond_alias_limit"])

    def test_despike_gaps_lower_the_limit(self):
        # Despiking deletes rows; the trapezoid then spans each gap with one
        # chord. On Mn3Sn 59438 (dQ = 0.01, despike on) the widest gap is 0.13,
        # and beyond pi/0.13 = 24 A the chords' error is 25-40 % rms of G(r).
        q = np.arange(50, 2801) * 0.01
        keep = np.ones(q.size, bool)
        keep[1000:1012] = False  # one 0.13-wide gap
        self.assertAlmostEqual(alias_limit(q), np.pi / 0.01, places=6)
        self.assertAlmostEqual(alias_limit(q[keep]), np.pi / 0.13, places=6)


if __name__ == "__main__":
    unittest.main()
