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

# numpy >= 2.0 renamed trapz; the package supports numpy >= 1.22.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz

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
    return 4 * np.pi * RHO0 * _trapezoid(r[window] ** 2 * g[window], r[window])


def run_cli(args):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main([str(arg) for arg in args])
    return code, out.getvalue(), err.getvalue()


def cli_outputs(tmp, stem):
    """(r, GK written for RMCProfile, GK before enforcement, provenance)."""
    rmc = read_stog_xy(Path(tmp) / f"{stem}_rmc.gr")
    ft = read_stog_xy(Path(tmp) / f"{stem}_ft.gr")  # classic scale_ft.gr: column 2 is g(r)
    provenance = json.loads((Path(tmp) / f"{stem}_provenance.json").read_text())
    return rmc[0], rmc[1], ft[1] - 1.0, provenance  # g - 1 before enforcement


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


class CliEnforcementHelpTests(unittest.TestCase):
    def test_help_describes_the_default_enforcement(self):
        from rmc_toolkits.scaling_cli import build_parser

        text = " ".join(build_parser().format_help().split())
        # Pre-fix: "(default: on in stog.inp mode, off in --data mode)", while
        # data mode actually enforced automatically.
        self.assertNotIn("off in --data mode", text)
        self.assertIn("foot of the detected first shell", text)
        self.assertIn("--no-enforce to disable", text)


class PinnedR0EnforcementTests(unittest.TestCase):
    """A given closest approach (--r0, MINIMUM_DISTANCES, stog.inp) is never overridden upward."""

    def run_synthetic(self, extra, header=""):
        q, sq = exact_sq(0.10, 26.0)  # the data's first shell starts at ~2.6 A
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "shell.dat"
            rows = "\n".join(f"{x:.6f} {y:.10f}" for x, y in zip(q, sq))
            data.write_text(header + rows + "\n")
            code, out, err = run_cli([
                "--data", data, "--qmin", "0.01", "--qmax", "26", "--rho0", RHO0,
                "--b-avg-sq", "1.0", "--scale", "1", "--offset", "0",
                "--rmax", "20", "--nr", "2000", "--out-dir", Path(tmp) / "out", *extra,
            ])
            self.assertEqual(code, 0, err)
            return cli_outputs(Path(tmp) / "out", "shell"), out

    def test_cli_r0_caps_the_automatic_cutoff(self):
        (r, gk_rmc, gm1_ft, provenance), out = self.run_synthetic(["--r0", "2.2"])
        # Pre-fix the cutoff came from the detected onset (~2.66 -> 2.40 A),
        # above the closest approach the user declared.
        self.assertLessEqual(provenance["enforcement"]["cutoff"], 2.2 - 0.25 + 1e-9)
        above = r > 2.2 - 0.25
        np.testing.assert_allclose(gk_rmc[above], gm1_ft[above], atol=1e-6)

    def test_minimum_distances_header_caps_the_automatic_cutoff(self):
        (r, gk_rmc, gm1_ft, provenance), _ = self.run_synthetic(
            [], header="MINIMUM_DISTANCES :: 2.3 2.2\n"
        )
        self.assertLessEqual(provenance["enforcement"]["cutoff"], 2.2 - 0.25 + 1e-9)

    def test_library_cap_and_conflict_flag(self):
        from rmc_toolkits.scaling import auto_enforcement_cutoff, diagnostics_summary

        q, sq = exact_sq(0.10, 26.0)
        base = dict(qmin=0.01, qmax=26.0, rho0=RHO0, b_avg_sq=1.0, rmax=20.0, nr=2000)
        free = scale_pipeline(q, sq, ScalingConfig(**base), 1.0, 0.0)
        cut_free = auto_enforcement_cutoff(free.r, free.g_filtered, ScalingConfig(**base))
        pinned_config = ScalingConfig(**base, r0=2.2)
        cut_pinned = auto_enforcement_cutoff(free.r, free.g_filtered, pinned_config)
        self.assertGreater(cut_free, 2.3)
        self.assertLessEqual(cut_pinned, 1.95 + 1e-12)
        # A pinned r0 ABOVE the detected shell is respected but flagged.
        high = ScalingConfig(**base, r0=3.0)
        result = scale_pipeline(q, sq, high, 1.0, 0.0)
        result.provenance["r0_detected"] = 2.66
        self.assertTrue(diagnostics_summary(result, high)["first_shell_below_r0"])
        low = ScalingConfig(**base, r0=2.6)
        result = scale_pipeline(q, sq, low, 1.0, 0.0)
        result.provenance["r0_detected"] = 2.66
        self.assertFalse(diagnostics_summary(result, low)["first_shell_below_r0"])


class PinnedR0WithoutDetectedShellTests(unittest.TestCase):
    """A given r0 anchors the automatic cutoff when no shell is detected (review follow-up).

    d6d646b let auto_enforcement_cutoff fall back to the given r0, but the report
    note still formatted the (None) detected onset: the CLI died with a TypeError
    before writing anything. A scale of 0.02 leaves every |g| feature below the
    detector's floor, so detection returns None.
    """

    def run_faint(self, extra, header=""):
        q, sq = exact_sq(0.10, 26.0)
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "faint.dat"
            rows = "\n".join(f"{x:.6f} {y:.10f}" for x, y in zip(q, sq))
            data.write_text(header + rows + "\n")
            code, out, err = run_cli([
                "--data", data, "--qmin", "0.01", "--qmax", "26", "--rho0", RHO0,
                "--b-avg-sq", "1.0", "--scale", "0.02", "--offset", "0",
                "--rmax", "20", "--nr", "2000", "--out-dir", Path(tmp) / "out", *extra,
            ])
            self.assertEqual(code, 0, err)
            written = sorted(p.name for p in (Path(tmp) / "out").glob("*"))
            provenance = json.loads((Path(tmp) / "out" / "faint_provenance.json").read_text())
        self.assertIn("faint_rmc.gr", written)
        self.assertIsNone(provenance["diagnostics"].get("r0_detected"))
        return provenance, out

    def test_cli_r0_anchors_the_cutoff_when_no_shell_is_detected(self):
        provenance, out = self.run_faint(["--r0", "2.4"])
        self.assertAlmostEqual(provenance["enforcement"]["cutoff"], 2.4 - 0.25, places=6)
        self.assertIn("given r0 2.4 A (no shell detected)", out)

    def test_minimum_distances_header_anchors_the_cutoff(self):
        provenance, out = self.run_faint([], header="MINIMUM_DISTANCES :: 2.4\n")
        self.assertAlmostEqual(provenance["enforcement"]["cutoff"], 2.4 - 0.25, places=6)
        self.assertIn("given r0 2.4 A (no shell detected)", out)

    def test_note_names_the_given_r0_when_it_caps_the_detected_onset(self):
        q, sq = exact_sq(0.10, 26.0)  # onset ~2.66 A, above the given 2.2
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "shell.sq"
            write_stog_xy(data, q, sq)
            code, out, err = run_cli([
                "--data", data, "--qmin", "0.01", "--qmax", "26", "--rho0", RHO0,
                "--b-avg-sq", "1.0", "--scale", "1", "--offset", "0", "--r0", "2.2",
                "--rmax", "20", "--nr", "2000", "--out-dir", Path(tmp) / "out",
            ])
        self.assertEqual(code, 0, err)
        self.assertIn("given r0 2.2 A", out)


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
        # Qmax 27: at the expert's 28 the unpinned density fit refuses since
        # 1.0 (tests/test_stog_a_detection.py), and writes nothing.
        with tempfile.TemporaryDirectory() as tmp:
            code, _, err = run_cli([
                "--data", MN3SN_59438, "--qmin", "1.0", "--qmax", "28",
                "--rho0", "0.063049", "--formula", "Mn3Sn", "--out-dir", tmp,
            ])
            self.assertEqual(code, 2)
            self.assertIn("could not locate the first coordination shell", err)
            self.assertEqual(list(Path(tmp).iterdir()), [])
        with tempfile.TemporaryDirectory() as tmp:
            code, _, err = run_cli([
                "--data", MN3SN_59438, "--qmin", "1.0", "--qmax", "27",
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
