# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Conditioning of the Q->0 Faber-Ziman amplitude (0.6.0 audit, stog-b group).

a_fz = (s0 - 1)/(S_meas(0) - L) had no conditioning guard: when the extrapolated
head lands near the high-Q level, the denominator is a small difference of noisy
numbers (74 -> 98 -> 91 -> 141 -> 512 on the Mn3Sn 59438 run as Qmin moves by
0.02-0.05 A^-1), and nothing distinguished such a scale from a trustworthy one.
fz_limit_fit now returns the standard error of S_meas(0) - L and a ``reliable``
flag (relative error <= FZ_REL_SE_MAX); the summary, CLI, page and rho0 estimate
report it.
"""

import contextlib
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import read_stog_xy, write_stog_xy
from rmc_toolkits.scaling import (
    FZ_REL_SE_MAX,
    ScalingConfig,
    amplitude_from_fz_limit,
    autoscale,
    diagnostics_summary,
    fz_limit_fit,
    level_sweep,
)
from rmc_toolkits.scaling_cli import main
from rmc_toolkits.scattering import faber_ziman
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

ROOT = Path(__file__).resolve().parents[1]
MN3SN = ROOT / "data" / "stog_tests"
RHO0, B2 = 0.05, 0.02


def model():
    """The shared synthetic (a = 10, b = -9) and its true-S(0) <b^2>."""
    q = np.arange(20, 981) * 0.03
    r = np.arange(1, 12001) * 0.005
    g = 0.5 * (1 + np.tanh((r - 2.65) / 0.07)) + 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    sq_true = fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, RHO0), q))
    head = q <= q[0] + 1.0
    _, s_true_0 = np.polyfit(q[head], sq_true[head], 1)
    return q, (sq_true + 9.0) / 10.0, B2 * (1.0 - float(s_true_0))


def flat_head(q, sq):
    """The head replaced by noise just below the level: S_meas(0) - L ~ its error."""
    level = level_sweep(q, sq).level
    out = sq.copy()
    flat = q <= 1.6
    out[flat] = level - 0.08 + np.random.default_rng(11).normal(0.0, 0.05, int(flat.sum()))
    return out


def config(b_sq_avg, **overrides):
    return ScalingConfig(
        qmin=0.6, qmax=30.0, rho0=RHO0, b_avg_sq=B2, b_sq_avg=b_sq_avg, r0=2.5,
        r_fit_min=1.2, r_fit_max=2.25, rmax=25.0, nr=1000, **overrides,
    )


class FzConditioningTests(unittest.TestCase):
    def test_clean_head_is_reliable(self):
        q, sq, b_sq_avg = model()
        sweep = level_sweep(q, sq)
        fit = fz_limit_fit(q, sq, sweep.level, config(b_sq_avg),
                           level_uncertainty=sweep.level_uncertainty)
        self.assertTrue(fit["reliable"])
        self.assertLess(fit["a_fz_rel_se"], FZ_REL_SE_MAX)
        self.assertEqual(fit["a_fz"], amplitude_from_fz_limit(q, sq, sweep.level, config(b_sq_avg)))
        self.assertLess(abs(fit["a_fz"] - 10.0) / 10.0, 0.01)

    def test_head_within_noise_of_the_level_is_flagged(self):
        q, sq, b_sq_avg = model()
        bad = flat_head(q, sq)
        result = autoscale(q, bad, config(b_sq_avg, amplitude_criterion="fz"))
        summary = diagnostics_summary(result, config(b_sq_avg, amplitude_criterion="fz"))
        self.assertGreater(result.a, 30.0)  # a positive, confidently wrong scale (truth 10)
        self.assertFalse(summary["a_fz_reliable"])
        self.assertGreater(summary["a_fz_rel_se"], 0.5)
        # Density mode reports the same flag next to the concordance.
        density = autoscale(q, bad, config(b_sq_avg))
        self.assertFalse(diagnostics_summary(density, config(b_sq_avg))["a_fz_reliable"])

    def test_cli_warns(self):
        q, sq, b_sq_avg = model()
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "flat.sq"
            write_stog_xy(data, q, flat_head(q, sq))
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = main([
                    "--data", str(data), "--qmin", "0.6", "--qmax", "30", "--rho0", str(RHO0),
                    "--b-avg-sq", str(B2), "--b-sq-avg", str(b_sq_avg), "--r0", "2.5",
                    "--r-fit-min", "1.2", "--rmax", "25", "--nr", "1000", "--amplitude", "fz",
                    "--out-dir", str(Path(tmp) / "out"),
                ])
        self.assertEqual(code, 0, err.getvalue())
        self.assertIn("Faber-Ziman amplitude a_fz", out.getvalue())
        self.assertIn("ill-conditioned", out.getvalue())

    @unittest.skipUnless(
        all((MN3SN / name).exists() for name in ("stog_59438", "stog_300K")),
        "Mn3Sn PG3 runs not present",
    )
    def test_mn3sn_59438_is_ill_conditioned_300k_is_not(self):
        fz = faber_ziman("Mn3Sn")
        runs = {"stog_59438": ("PG3_59438_SQ_rebin.dat", False), "stog_300K": ("PG3_55526_SQ_rebin.sq", True)}
        for name, (raw, reliable) in runs.items():
            q, sq = read_stog_xy(MN3SN / name / raw)[:2]
            for qmin in (0.82, 1.0):
                cfg = ScalingConfig(
                    qmin=qmin, qmax=28.0, rho0=0.063049, b_avg_sq=fz.b_avg_sq_barn,
                    b_sq_avg=fz.b_sq_avg_barn, amplitude_criterion="fz",
                )
                summary = diagnostics_summary(autoscale(q, sq, cfg), cfg)
                with self.subTest(run=name, qmin=qmin):
                    self.assertIs(summary["a_fz_reliable"], reliable)


class ReliableIsNotSufficientTests(unittest.TestCase):
    """reliable=True only says the denominator is statistically resolved (0.6.0 review).

    On two of the three 'good' Mn3Sn runs the reliable-flagged a_fz still
    drifts ~45 % with Qmin (a systematic head bias), so the CLI says what else
    to check whenever it reports a reliable a_fz.
    """

    def run_cli(self, amplitude):
        q, sq, b_sq_avg = model()
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "clean.sq"
            write_stog_xy(data, q, sq)
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = main([
                    "--data", str(data), "--qmin", "0.6", "--qmax", "30", "--rho0", str(RHO0),
                    "--b-avg-sq", str(B2), "--b-sq-avg", str(b_sq_avg), "--r0", "2.5",
                    "--r-fit-min", "1.2", "--rmax", "25", "--nr", "1000",
                    "--amplitude", amplitude, "--out-dir", str(Path(tmp) / "out"),
                ])
        self.assertEqual(code, 0, err.getvalue())
        return out.getvalue()

    def test_cli_qualifies_a_reliable_a_fz(self):
        for amplitude in ("fz", "density"):
            with self.subTest(amplitude=amplitude):
                out = self.run_cli(amplitude)
                self.assertIn("Q->0 amplitude", out)
                self.assertIn("necessary, not sufficient", out)
                self.assertIn("--qmin", out)
                self.assertNotIn("ill-conditioned", out)

    @unittest.skipUnless(
        all((MN3SN / name).exists() for name in ("stog", "stog_500K")),
        "Mn3Sn PG3 runs not present",
    )
    def test_reliable_a_fz_still_drifts_with_qmin(self):
        fz = faber_ziman("Mn3Sn")
        runs = {"stog": "PG3_55537_rebin.sq", "stog_500K": "PG3_54139_SQ_rebin.dat"}
        for name, raw in runs.items():
            q, sq = read_stog_xy(MN3SN / name / raw)[:2]
            values = []
            for qmin in (0.82, 1.02, 1.05):
                cfg = ScalingConfig(
                    qmin=qmin, qmax=28.0, rho0=0.063049, b_avg_sq=fz.b_avg_sq_barn,
                    b_sq_avg=fz.b_sq_avg_barn,
                )
                crop = q >= qmin
                sweep = level_sweep(q[crop & (q <= 28.0)], sq[crop & (q <= 28.0)])
                fit = fz_limit_fit(
                    q[crop & (q <= 28.0)], sq[crop & (q <= 28.0)], sweep.level, cfg,
                    level_uncertainty=sweep.level_uncertainty,
                )
                with self.subTest(run=name, qmin=qmin):
                    self.assertTrue(fit["reliable"])
                values.append(fit["a_fz"])
            with self.subTest(run=name):  # 55537: 11.0 -> 6.2; 54139: 16.3 -> 23.7
                self.assertGreater(max(values) / min(values), 1.4)


if __name__ == "__main__":
    unittest.main()
