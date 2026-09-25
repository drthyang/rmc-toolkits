# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Automatic low-r enforcement must never remove first-shell signal from the RMC files.

Regression tests for the 1.0 audit (stog-a group). The automatic enforcement cutoff
used to be the detected first-shell onset itself -- a point ~35 % up the shell's
rising flank -- and first_peak_zero zeroed every r <= onset, deleting 6-9 % of the
first-shell pair density from <stem>_rmc.gr / _rmc.dr (the files RMCProfile fits).
"""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import read_stog_xy, write_stog_xy
from rmc_toolkits.scaling import ScalingConfig, scale_pipeline
from rmc_toolkits.scaling_cli import main
from rmc_toolkits.transforms import first_peak_zero, fq_to_sq, g_to_gpdf, gpdf_to_fq

ROOT = Path(__file__).resolve().parents[1]
STOG_TESTS = ROOT / "data" / "stog_tests"
FECOSN = STOG_TESTS / "199K" / "FeCoSn_199K_rebinned.sq"
MN3SN_59438 = STOG_TESTS / "stog_59438" / "PG3_59438_SQ_rebin.dat"

RHO0 = 0.05
R_TRUE = np.arange(1, 16001) * 0.005


def shell_g(r, sigma):
    """One Gaussian coordination shell at 2.8 A (+ a continuum from 3.4 A)."""
    return 1.6 * np.exp(-0.5 * ((r - 2.8) / sigma) ** 2) + 0.5 * (1.0 + np.tanh((r - 3.4) / 0.08))


def exact_sq(sigma, qmax):
    """S(Q) from Q = 0.01: nothing omitted, so the true scale is a = 1, b = 0."""
    q = np.arange(1, round(qmax / 0.01) + 1) * 0.01
    sq = fq_to_sq(q, gpdf_to_fq(R_TRUE, g_to_gpdf(R_TRUE, shell_g(R_TRUE, sigma), RHO0), q))
    return q, sq


def coordination(r, g, lo=2.0, hi=3.2):
    """First-shell coordination number 4 pi rho0 int r^2 g dr over [lo, hi]."""
    window = (r >= lo) & (r <= hi)
    return 4 * np.pi * RHO0 * np.trapezoid(r[window] ** 2 * g[window], r[window])


def run_cli(args):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main([str(arg) for arg in args])
    return code, out.getvalue(), err.getvalue()


def cli_outputs(tmp, stem):
    """(r, GK written for RMCProfile, GK before enforcement, provenance)."""
    rmc = read_stog_xy(Path(tmp) / f"{stem}_rmc.gr")
    ft = read_stog_xy(Path(tmp) / f"{stem}_ft.gr")  # g - 1 before enforcement
    provenance = json.loads((Path(tmp) / f"{stem}_provenance.json").read_text())
    return rmc[0], rmc[1], ft[1], provenance


class AutoEnforcementCoordinationTests(unittest.TestCase):
    def test_cli_data_mode_keeps_the_first_shell(self):
        # The default data-mode workflow: no --enforce flag, auto enforcement.
        q, sq = exact_sq(0.10, 26.0)
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "shell.sq"
            write_stog_xy(data, q, sq)
            code, out, err = run_cli([
                "--data", data, "--qmin", "0.01", "--qmax", "26", "--rho0", RHO0,
                "--b-avg-sq", "1.0", "--scale", "1", "--offset", "0",
                "--rmax", "20", "--nr", "2000", "--out-dir", Path(tmp) / "out",
            ])
            self.assertEqual(code, 0, err)
            r, gk_rmc, gm1_ft, provenance = cli_outputs(Path(tmp) / "out", "shell")
        enforcement = provenance["enforcement"]
        g_rmc, g_ft = gk_rmc + 1.0, gm1_ft + 1.0  # b_avg_sq = 1
        cn_rmc, cn_ft = coordination(r, g_rmc), coordination(r, g_ft)
        # Pre-1.0 the cutoff was the onset (~2.66 A) and CN dropped by 7.4 %.
        self.assertLess(abs(cn_rmc / cn_ft - 1.0), 0.005, (cn_rmc, cn_ft))
        self.assertLess(enforcement["cutoff"], 2.45)
        self.assertEqual(enforcement["source"], "auto (first-shell foot)")
        self.assertIn("foot of the first shell", out)

    def test_coordination_preserved_across_widths_and_resolution(self):
        from rmc_toolkits.scaling import auto_enforcement_cutoff, detect_first_peak_onset

        for sigma in (0.08, 0.10, 0.15):
            for qmax, lorch in ((26.0, False), (40.0, False), (26.0, True)):
                with self.subTest(sigma=sigma, qmax=qmax, lorch=lorch):
                    q, sq = exact_sq(sigma, qmax)
                    config = ScalingConfig(
                        qmin=0.01, qmax=qmax, rho0=RHO0, b_avg_sq=1.0, lorch=lorch,
                        rmax=20.0, nr=2000,
                    )
                    result = scale_pipeline(q, sq, config, 1.0, 0.0)
                    r, g = result.r, result.g_filtered
                    onset = detect_first_peak_onset(r, g, qmax, search_min=1.3)
                    cutoff = auto_enforcement_cutoff(r, g, config)
                    self.assertLess(cutoff, onset - 0.2)
                    enforced = first_peak_zero(r, g, cutoff=cutoff, peak_rmin=cutoff, peak_rmax=cutoff)
                    cn_pre, cn_enforced = coordination(r, g), coordination(r, enforced)
                    self.assertLess(abs(cn_enforced / cn_pre - 1.0), 0.005)
                    # And the shell itself: every point above the cutoff is untouched.
                    np.testing.assert_array_equal(enforced[r > cutoff], g[r > cutoff])
                    # The pre-1.0 cutoff (the onset) lost 6-9 % of the shell.
                    at_onset = first_peak_zero(r, g, cutoff=onset, peak_rmin=onset, peak_rmax=onset)
                    self.assertLess(coordination(r, at_onset) / cn_pre, 0.96)

    def test_foot_walks_to_the_sign_change_or_local_minimum(self):
        from rmc_toolkits.scaling import first_shell_foot

        r = np.arange(1, 501) * 0.01
        g = np.where(r < 2.0, -0.1 * np.sin(np.pi * (2.0 - r) / 0.2), 0.0)
        g = g + 2.0 * np.exp(-0.5 * ((r - 2.4) / 0.1) ** 2) * (r >= 2.0)
        # Starting on the flank at 2.25 A, the walk stops just across the sign
        # change at 2.0 A (g < 0 below it).
        self.assertLess(first_shell_foot(r, g, 2.25), 2.0)
        self.assertGreater(first_shell_foot(r, g, 2.25), 1.97)


@unittest.skipUnless(FECOSN.exists(), "FeCoSn 199K run not present")
class FeCoSnAutoEnforcementTests(unittest.TestCase):
    def test_first_peak_flank_survives(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, _, err = run_cli([
                "--data", FECOSN, "--qmin", "0.5", "--qmax", "26", "--rho0", "0.057329",
                "--b-avg-sq", "1.0", "--out-dir", tmp,
            ])
            self.assertEqual(code, 0, err)
            r, gk_rmc, gm1_ft, provenance = cli_outputs(tmp, "FeCoSn_199K_rebinned")
        cutoff = provenance["enforcement"]["cutoff"]
        # Pre-1.0: 2.53 A, on the flank of the 2.64 A shell (g(2.53) = 2.24 zeroed).
        self.assertLess(cutoff, 2.4)
        shell = (r > 2.4) & (r < 3.0)
        np.testing.assert_allclose(gk_rmc[shell], gm1_ft[shell], atol=1e-6)


@unittest.skipUnless(MN3SN_59438.exists(), "stog_59438 example run not present")
class Mn3SnAutoEnforcementTests(unittest.TestCase):
    def test_cutoff_matches_the_expert_below_the_inverted_first_shell(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, _, err = run_cli([
                "--data", MN3SN_59438, "--qmin", "1.0", "--qmax", "28",
                "--rho0", "0.063049", "--formula", "Mn3Sn", "--out-dir", tmp,
            ])
            self.assertEqual(code, 0, err)
            r, gk_rmc, gm1_ft, provenance = cli_outputs(tmp, "PG3_59438_SQ_rebin")
        cutoff = provenance["enforcement"]["cutoff"]
        # The expert's hand cutoff (rmccut) is 2.48 A with the first peak at
        # 2.65-3.1 A; pre-1.0 the auto cutoff was 3.49 A (the whole shell erased).
        self.assertLess(cutoff, 2.55)
        self.assertGreater(cutoff, 2.2)
        b2 = provenance["diagnostics"]["gk_low_r_theory"] * -1.0
        shell = (r > 2.55) & (r < 3.1)
        np.testing.assert_allclose(gk_rmc[shell], b2 * gm1_ft[shell], atol=1e-6)


if __name__ == "__main__":
    unittest.main()
