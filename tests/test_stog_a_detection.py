# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""First-coordination-shell detection (Auto StoG r0): the FIRST shell, of either sign.

Regression tests for the 1.0 audit (stog-a group): the detector used to return the
flank of the *strongest* |g| feature, so a weak or inverted first shell (Ti-O in
titanates, Mn-Sn in Mn3Sn) was skipped and the low-r window / enforcement landed on
the real first shell.
"""

from pathlib import Path
import unittest

import numpy as np

from rmc_toolkits.parsers import read_stog_xy
from rmc_toolkits.scaling import (
    ScalingConfig,
    autoscale,
    detect_first_peak_onset,
    diagnostics_summary,
)
from rmc_toolkits.scattering import faber_ziman

ROOT = Path(__file__).resolve().parents[1]
STOG_59438 = ROOT / "data" / "stog_tests" / "stog_59438" / "PG3_59438_SQ_rebin.dat"

R = np.arange(1, 801) * 0.01  # 0.01 .. 8 A


def gauss(r, centre, sigma):
    return np.exp(-0.5 * ((r - centre) / sigma) ** 2)


def continuum(r, start=3.4):
    return 0.5 * (1.0 + np.tanh((r - start) / 0.08))


def inverted_first_shell_g(r):
    """SrTiO3-like neutron total: inverted Ti-O (depth 2.8) below a +10 A-O/O-O shell."""
    return -2.8 * gauss(r, 1.95, 0.07) + 10.0 * gauss(r, 2.76, 0.09) + continuum(r)


def weak_first_shell_g(r):
    """A weak positive first shell (0.6) 0.8 A below a strong one (5.6)."""
    return 0.6 * gauss(r, 2.1, 0.08) + 5.6 * gauss(r, 2.9, 0.09) + continuum(r, 3.6)


class FirstShellDetectorTests(unittest.TestCase):
    def test_inverted_first_shell_weaker_than_the_second_is_found(self):
        onset = detect_first_peak_onset(R, inverted_first_shell_g(R), 30.0, search_min=1.3)
        self.assertIsNotNone(onset)
        # Ti-O shell at 1.95 A (sigma 0.07): its 35%-height flank is ~1.85 A. The
        # old argmax-|g| detector returned the flank of the 2.76 A shell (~2.6 A).
        self.assertGreater(onset, 1.75)
        self.assertLess(onset, 1.95)

    def test_weak_first_shell_below_a_strong_second_is_found(self):
        onset = detect_first_peak_onset(R, weak_first_shell_g(R), 30.0, search_min=1.3)
        self.assertIsNotNone(onset)
        self.assertGreater(onset, 1.9)
        self.assertLess(onset, 2.1)

    def test_positive_single_shell_unchanged(self):
        g = 1.6 * gauss(R, 2.8, 0.15) + 0.5 * (1.0 + np.tanh((R - 2.65) / 0.07))
        onset = detect_first_peak_onset(R, g, 30.0, search_min=1.3)
        self.assertGreater(onset, 2.5)
        self.assertLess(onset, 2.7)

    def test_ripples_comparable_to_the_first_shell_are_not_a_shell(self):
        # A Mn3Sn-first-pass-like profile: a +-0.9 ripple field, then the inverted
        # first shell (-3.7) and a positive second shell. Neither a ripple crest nor
        # the (later, but not taller) second shell may be taken as the first shell.
        g = 0.9 * np.sin(2 * np.pi * (R - 1.3) / 0.24) * ((R > 1.3) & (R < 2.55))
        g = g - 3.7 * gauss(R, 2.84, 0.06) + 1.8 * gauss(R, 4.0, 0.08) + continuum(R, 4.4)
        onset = detect_first_peak_onset(R, g, 28.0, search_min=1.3)
        self.assertIsNotNone(onset)
        self.assertGreater(onset, 2.65)
        self.assertLess(onset, 2.84)

    def test_shell_at_the_search_start_is_not_replaced_by_a_later_one(self):
        # A B-O-like first shell at 1.37 A straddles search_min = 1.3: it cannot
        # be located, but the O-O shell at 2.4 A must not be reported as the
        # first shell either (the old detector returned its flank, ~2.3 A).
        g = 6.0 * gauss(R, 1.37, 0.045) + 3.0 * gauss(R, 2.4, 0.08) + continuum(R, 3.0)
        self.assertIsNone(detect_first_peak_onset(R, g, 30.0, search_min=1.3))

    def test_no_feature_returns_none(self):
        g = 0.2 * np.sin(R * 25.0)
        self.assertIsNone(detect_first_peak_onset(R, g, 25.0, search_min=1.3))


@unittest.skipUnless(STOG_59438.exists(), "stog_59438 example run not present")
class Mn3Sn59438DetectionTests(unittest.TestCase):
    """Real data: the 59438 run's inverted Mn-Sn first shell sits at 2.65-3.1 A."""

    def test_composition_only_run_detects_the_first_shell(self):
        data = read_stog_xy(STOG_59438)
        fz = faber_ziman("Mn3Sn")
        config = ScalingConfig(
            qmin=1.0, qmax=28.0, rho0=0.063049,
            b_avg_sq=fz.b_avg_sq_barn, b_sq_avg=fz.b_sq_avg_barn,
        )
        result = autoscale(data[0], data[1], config)
        summary = diagnostics_summary(result, config)
        # The old detector locked onto the second shell (r0 = 3.49 A) and refined
        # the C2 window to [1.2, 3.24] -- across the first shell.
        self.assertGreater(summary["r0_detected"], 2.4)
        self.assertLess(summary["r0_detected"], 2.9)
        self.assertLess(summary["r_fit_window"][1], 2.65)


if __name__ == "__main__":
    unittest.main()
