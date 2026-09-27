# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Auto StoG low-r window placement: confirm the first shell, never return a <= 0.

Regression tests for the 0.6.0 review of the stog-a group. On the real Mn3Sn 59438
run (Qmin 0.82 / 1.0, Qmax 24-30) the first window-placement loop

- accepted a sub-shell ripple lobe (~1.65 A, 2.9x its ripple field at 35 % of the
  range maximum) as the first shell,
- returned the refit below it although its scale was negative (a = -0.48 .. -0.72),
- "confirmed" it by a re-detected onset far ABOVE the one the window was built
  from, and reported that unrelated onset as r0_detected.

The placement logic is exercised on scripted passes (``SCENARIOS``, shared with the
JS port through tests/generate_autoscale_fixture.py) and end to end on the real run.
Set RMC_TOOLKITS_FULL_SWEEP=1 to sweep the whole Qmin {0.82, 1.0} x Qmax 24-30 grid
of the real run (~5 min); by default one configuration runs.
"""

import os
from pathlib import Path
from types import SimpleNamespace
import unittest

import numpy as np

from rmc_toolkits.parsers import read_stog_xy
from rmc_toolkits.scaling import (
    MAX_WINDOW_REFITS,
    ScalingConfig,
    _place_low_r_window,
    autoscale,
    detect_first_peak_onset,
    diagnostics_summary,
    first_shell_candidates,
)
from rmc_toolkits.scattering import faber_ziman

ROOT = Path(__file__).resolve().parents[1]
STOG_59438 = ROOT / "data" / "stog_tests" / "stog_59438" / "PG3_59438_SQ_rebin.dat"

R = np.arange(1, 801) * 0.01  # 0.01 .. 8 A


def gauss(r, centre, sigma):
    return np.exp(-0.5 * ((r - centre) / sigma) ** 2)


def continuum(r, start=3.4):
    return 0.5 * (1.0 + np.tanh((r - start) / 0.08))


def mn3sn_lobe_g(r):
    """The measured 59438 trial profile (Qmin 1.0, Qmax 25, trial window [1.2, 2.2]).

    A +-0.22 ripple field, a +0.643 lobe at 1.65 A (2.9x that field, 35 % of the
    range maximum: it passed the pre-review 2x / 35 % rule) and the inverted Mn-Sn
    first shell (-1.82 at 2.84 A).
    """
    ripple = 0.22 * np.sin(2 * np.pi * (r - 1.3) / 0.2) * ((r > 1.3) & (r < 1.5))
    return ripple + 0.643 * gauss(r, 1.65, 0.04) - 1.821 * gauss(r, 2.84, 0.06) + continuum(r)


#: Scripted-pass profiles g(r) on R (the placement scenarios below and the JS fixture).
PROFILES = {
    "flat": lambda r: 0.0 * r,
    "shell": lambda r: -3.0 * gauss(r, 2.84, 0.06) + continuum(r),
    "lobeOnly": lambda r: 1.2 * gauss(r, 1.62, 0.04) + 0.4 * continuum(r),
    "lobeShell": lambda r: 1.2 * gauss(r, 1.62, 0.04) - 3.0 * gauss(r, 2.84, 0.06) + continuum(r),
    "secondOnly": lambda r: 4.0 * gauss(r, 2.76, 0.08) + continuum(r),
    "twoShell": lambda r: 1.5 * gauss(r, 1.95, 0.06) + 4.0 * gauss(r, 2.76, 0.08) + continuum(r),
    "lobeSecond": lambda r: 1.2 * gauss(r, 1.62, 0.04) + 4.0 * gauss(r, 2.76, 0.08) + continuum(r),
    "short": lambda r: 2.0 * gauss(r, 1.58, 0.03) + continuum(r, 3.0),
    "veryShort": lambda r: 2.0 * gauss(r, 1.48, 0.02) + continuum(r, 3.0),
    "pair32": lambda r: 1.5 * gauss(r, 3.3, 0.06) + 4.0 * gauss(r, 3.9, 0.08) + continuum(r, 4.4),
    "pair29": lambda r: 1.5 * gauss(r, 3.0, 0.06) + 1.5 * gauss(r, 3.3, 0.06) + continuum(r, 4.4),
    "pair26": lambda r: 1.5 * gauss(r, 2.7, 0.06) + 1.5 * gauss(r, 3.0, 0.06) + continuum(r, 4.4),
    "pair23": lambda r: 1.5 * gauss(r, 2.4, 0.06) + 1.5 * gauss(r, 2.7, 0.06) + continuum(r, 4.4),
    "pair20": lambda r: 1.5 * gauss(r, 2.1, 0.06) + 1.5 * gauss(r, 2.4, 0.06) + continuum(r, 4.4),
}

#: Placement scenarios: trial passes keyed by window width, refits keyed by the r0
#: they are built from (matched within 0.05 A; ``default`` otherwise). ``expect``
#: holds what the placement must do.
SCENARIOS = {
    # The 59438 failure: the ripple lobe is proposed by the trial; its refit shows
    # only the real shell (far above): the lobe is dropped, the shell confirmed.
    "rippleDropped": {
        "qmax": 28.0,
        "trials": [[0.3, -0.25, "shell"], [1.0, 0.98, "lobeShell"]],
        "refits": [[1.57, 1.1, "shell"], [2.76, 1.2, "lobeShell"]],
        "default": [1.0, "flat"],
        "expect": {"a": 1.2, "r0": 2.76},
    },
    # A refit with a <= 0 is never returned (pre-review: a = -0.59 came back).
    "negativeRefit": {
        "qmax": 28.0,
        "trials": [[0.3, -0.1, "flat"], [1.0, 0.98, "lobeShell"]],
        "refits": [[1.57, -0.59, "shell"]],
        "default": [1.0, "flat"],
        "expect": {"error": "non-physical scale"},
    },
    # A refit across a lower shell uncovers it: the lower shell is confirmed.
    "lowerShellUncovered": {
        "qmax": 28.0,
        "trials": [[0.3, 5.0, "secondOnly"], [1.0, 5.0, "secondOnly"]],
        "refits": [[2.65, 8.0, "twoShell"], [1.87, 10.0, "twoShell"]],
        "default": [1.0, "flat"],
        "expect": {"a": 10.0, "r0": 1.87},
    },
    # The uncovered lower feature is not re-detected by its own refit: it is
    # dropped and the shell the loop came from is confirmed after all.
    "lowerDroppedThenBack": {
        "qmax": 28.0,
        "trials": [[0.3, 1.0, "secondOnly"], [1.0, 1.0, "secondOnly"]],
        "refits": [[2.65, 1.2, "lobeSecond"], [1.57, 1.1, "secondOnly"]],
        "default": [1.0, "flat"],
        "expect": {"a": 1.2, "r0": 2.65},
    },
    # Every refit uncovers another lower shell: the refit budget runs out.
    "budgetExhausted": {
        "qmax": 28.0,
        "trials": [[0.3, 1.0, "pair32"], [1.0, 1.0, "pair32"]],
        "refits": [[3.22, 1.0, "pair29"], [2.92, 1.0, "pair26"], [2.62, 1.0, "pair23"],
                   [2.32, 1.0, "pair20"]],
        "default": [1.0, "flat"],
        "expect": {"error": "no first-shell onset was confirmed within"},
    },
    # A confirmed shell too close to lo: the r_cutoff advice is given.
    "shortConfirmed": {
        "qmax": 28.0,
        "trials": [[0.3, 5.0, "short"], [1.0, -1.0, "flat"]],
        "refits": [[1.54, 6.0, "short"]],
        "default": [1.0, "flat"],
        "expect": {"error": "the first coordination shell starts at"},
    },
    # A candidate so close to lo that no window fits below it: conditional advice.
    "shortUnverifiable": {
        "qmax": 40.0,
        "trials": [[0.3, -1.5, "veryShort"], [1.0, -1.2, "veryShort"]],
        "refits": [],
        "default": [1.0, "flat"],
        "expect": {"error": "shell-like feature starts at"},
    },
    # A no-room candidate whose narrow verification fit is non-physical.
    "narrowNegative": {
        "qmax": 28.0,
        "trials": [[0.3, 5.0, "short"], [1.0, 5.0, "short"]],
        "refits": [[1.54, -0.3, "short"]],
        "default": [1.0, "flat"],
        "expect": {"error": "non-physical scale"},
    },
    # The only candidate vanishes on its own refit: fail, naming it.
    "notRedetected": {
        "qmax": 28.0,
        "trials": [[0.3, 1.0, "lobeOnly"], [1.0, 1.0, "lobeOnly"]],
        "refits": [[1.57, 1.0, "flat"]],
        "default": [1.0, "flat"],
        "expect": {"error": "could not locate the first coordination shell"},
    },
}

#: Error kinds, as phrases shared verbatim by the Python and JS messages.
ERROR_KINDS = (
    "non-physical scale",
    "no first-shell onset was confirmed within",
    "the first coordination shell starts at",
    "shell-like feature starts at",
    "could not locate the first coordination shell",
)


def scenario_config(scenario):
    return ScalingConfig(qmin=0.5, qmax=scenario["qmax"], rho0=0.05, b_avg_sq=1.0)


def scripted_pass(scenario, profiles, r=R):
    """A fake fit pass replaying ``scenario`` (same rules in the JS fixture test)."""
    calls = []

    def run(config):
        lo, hi = config.r_fit_window
        if hi - lo < 0.03:
            raise ValueError("fit windows contain fewer than 2 points")
        if config.r0 is None:
            a, name = next(
                (row[1], row[2]) for row in scenario["trials"] if abs(hi - lo - row[0]) < 1e-6
            )
        else:
            match = [row for row in scenario["refits"] if abs(row[0] - config.r0) <= 0.05]
            a, name = (match[0][1], match[0][2]) if match else tuple(scenario["default"])
        calls.append((config.r0, config.r_fit_max, a, name))
        return SimpleNamespace(
            a=float(a), r=r, g_filtered=profiles[name], provenance={"r_fit_window": [lo, hi]}
        )

    return run, calls


class FirstShellMarginTests(unittest.TestCase):
    """The acceptance rule needs a margin over real sub-shell ripple lobes."""

    def test_mn3sn_ripple_lobe_is_not_the_first_shell(self):
        g = mn3sn_lobe_g(R)
        # Pre-review rule (2x ripple at >= 35 % of the maximum, or 3x): the lobe.
        old = first_shell_candidates(R, g, 25.0, search_min=1.3, major=0.35, strong_prominence=3.0)
        self.assertLess(old[0], 1.65)
        onset = detect_first_peak_onset(R, g, 25.0, search_min=1.3)
        self.assertGreater(onset, 2.65)
        self.assertLess(onset, 2.84)

    def test_flank_reaching_the_search_start_is_not_a_shell(self):
        # A broad feature whose flank is still high at search_min: its onset is
        # not separable. Pre-review this returned search_min itself (1.30 A),
        # which then produced "lower r_cutoff to <= 0.75 A" on real Mn3Sn data.
        g = 3.0 * gauss(R, 1.58, 0.2) + continuum(R)
        self.assertIsNone(detect_first_peak_onset(R, g, 28.0, search_min=1.3))
        self.assertEqual(first_shell_candidates(R, g, 28.0, search_min=1.3), [])

    def test_candidates_list_every_shell_in_order(self):
        g = PROFILES["twoShell"](R)
        onsets = first_shell_candidates(R, g, 28.0, search_min=1.3)
        self.assertEqual(len(onsets), 2)
        self.assertLess(onsets[0], 1.95)
        self.assertGreater(onsets[1], 2.5)
        self.assertEqual(detect_first_peak_onset(R, g, 28.0, search_min=1.3), onsets[0])


class ScriptedPlacementTests(unittest.TestCase):
    """The trial / confirm loop on scripted passes."""

    def place(self, name):
        scenario = SCENARIOS[name]
        profiles = {key: build(R) for key, build in PROFILES.items()}
        run, calls = scripted_pass(scenario, profiles)
        return _place_low_r_window(run, scenario_config(scenario)), calls

    def test_expected_outcomes(self):
        for name, scenario in SCENARIOS.items():
            expect = scenario["expect"]
            with self.subTest(scenario=name):
                if "error" in expect:
                    with self.assertRaisesRegex(ValueError, expect["error"]):
                        self.place(name)
                    continue
                result, _ = self.place(name)
                self.assertGreater(result.a, 0)
                self.assertAlmostEqual(result.a, expect["a"])
                onset = result.provenance["r0_detected"]
                self.assertAlmostEqual(onset, expect["r0"], delta=0.02)
                # r0_detected is the onset the window was built from.
                self.assertAlmostEqual(result.provenance["r_fit_window"][1], onset - 0.25)
                self.assertTrue(result.provenance["window_refined"])

    def test_ripple_candidate_is_dropped_not_returned(self):
        result, calls = self.place("rippleDropped")
        # Refits at the lobe (1.57) and at the shell (2.76): the lobe fit (a = 1.1,
        # window [1.2, 1.31]) is not what comes back.
        self.assertEqual([round(call[0], 2) for call in calls if call[0] is not None], [1.57, 2.76])
        self.assertEqual(result.a, 1.2)

    def test_short_shell_advice_only_when_confirmed(self):
        with self.assertRaisesRegex(ValueError, r"starts at 1\.54 A.*Lower r_cutoff to <= 0\.95"):
            self.place("shortConfirmed")
        with self.assertRaisesRegex(ValueError, r"1\.46 A.*If it is the first coordination shell"):
            self.place("shortUnverifiable")

    def test_refit_budget_is_bounded(self):
        with self.assertRaises(ValueError):
            self.place("budgetExhausted")
        scenario = SCENARIOS["budgetExhausted"]
        profiles = {key: build(R) for key, build in PROFILES.items()}
        run, calls = scripted_pass(scenario, profiles)
        with self.assertRaises(ValueError):
            _place_low_r_window(run, scenario_config(scenario))
        refits = [call for call in calls if call[0] is not None]
        self.assertEqual(len(refits), MAX_WINDOW_REFITS)


MN3SN_SWEEP = [(qmin, qmax) for qmin in (0.82, 1.0) for qmax in range(24, 31)]


@unittest.skipUnless(STOG_59438.exists(), "stog_59438 example run not present")
class Mn3Sn59438PlacementTests(unittest.TestCase):
    """Real missing-low-Q data: every returned fit has a > 0 and a window below the shell."""

    @classmethod
    def setUpClass(cls):
        data = read_stog_xy(STOG_59438)
        cls.q, cls.sq = data[0], data[1]
        cls.fz = faber_ziman("Mn3Sn")

    def check(self, config):
        try:
            result = autoscale(self.q, self.sq, config)
        except ValueError as exc:  # raising is acceptable, a wrong fit is not
            self.assertNotIn("r_cutoff to", str(exc))
            return
        summary = diagnostics_summary(result, config)
        self.assertGreater(result.a, 0.0)
        hi = summary["r_fit_window"][1]
        self.assertLess(hi, 2.65)  # the inverted Mn-Sn shell spans 2.65-3.1 A
        self.assertAlmostEqual(hi, summary["r0_detected"] - 0.25, places=9)
        self.assertGreater(summary["r0_detected"], 2.4)
        self.assertLess(summary["r0_detected"], 2.9)
        # The honest verdict stays visible (a fit on a ripple had hidden it).
        self.assertFalse(summary["density_limit_satisfied"])

    def composition_config(self, qmin, qmax):
        return ScalingConfig(
            qmin=qmin, qmax=float(qmax), rho0=0.063049,
            b_avg_sq=self.fz.b_avg_sq_barn, b_sq_avg=self.fz.b_sq_avg_barn,
        )

    def test_qmax_29_returns_a_positive_scale_below_the_first_shell(self):
        # Pre-review: a = -0.592 on the window [1.2, 1.33], r0_detected 2.68.
        self.check(self.composition_config(1.0, 29))

    @unittest.skipUnless(os.environ.get("RMC_TOOLKITS_FULL_SWEEP"), "set RMC_TOOLKITS_FULL_SWEEP=1")
    def test_full_qmin_qmax_sweep(self):
        # Pre-review: a < 0 in 7 of 14 configurations, 3 raised with wrong advice.
        for qmin, qmax in MN3SN_SWEEP:
            with self.subTest(qmin=qmin, qmax=qmax):
                self.check(self.composition_config(qmin, qmax))
        # Without the composition (S(0) = 0 convention): pre-review "the first
        # coordination shell starts at 1.30 A ... Lower r_cutoff to <= 0.75 A".
        with self.subTest(qmin=1.0, qmax=25, composition=False):
            self.check(ScalingConfig(qmin=1.0, qmax=25.0, rho0=0.063049, b_avg_sq=0.015407))


if __name__ == "__main__":
    unittest.main()
