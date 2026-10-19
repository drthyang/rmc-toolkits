# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The displacement-modulation (superstructure) engine on synthetic configurations."""

from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cluster_configs import cluster_config_text, long_wave_config_text, wave_config_text  # noqa: E402
from rmc_toolkits.modulations import (  # noqa: E402
    LONG_WAVE_Q,
    LONG_WAVE_RATIO,
    MODE_KINDS,
    MODULATION_BECAUSE,
    MODULATION_RATIO,
    MODULATION_RULES,
    MODULATION_SIGMA,
    SIGMA_MAX,
    _folded,
    displacement_modulations,
    general_mode_sigma,
    is_long_wave,
    zone_boundary_points,
    zone_boundary_sigma,
)
from rmc_toolkits.pca_kde import load_site_displacements  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEMO_RMC6F = ROOT / "web_app" / "frontend" / "public" / "demo" / "GTS_250K.rmc6f"


def _load(mode: str, supercell=(6, 6, 6)):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / f"{mode}.rmc6f"
        path.write_text(cluster_config_text(mode, supercell=supercell), encoding="utf-8")
        return load_site_displacements(path)


def _element(result: dict, name: str) -> dict:
    return next(entry for entry in result["elements"] if entry["element"] == name)


class ZoneBoundaryTests(unittest.TestCase):
    def test_only_even_axes_can_be_halved(self):
        self.assertEqual(len(zone_boundary_points(np.array([4, 4, 4]))), 7)
        self.assertEqual(zone_boundary_points(np.array([3, 4, 5])), [(0, 1, 0)])


class ModulationTests(unittest.TestCase):
    def test_a_planted_cell_doubling_stands_out(self):
        result = displacement_modulations(_load("doubled"))
        self.assertEqual(result["supercell"], [6, 6, 6])
        self.assertEqual(result["wavevectors"], 215)
        nb = _element(result, "Nb")
        doubling = next(point for point in nb["zoneBoundary"] if point["q"] == [0.0, 0.0, 0.5])
        self.assertGreater(doubling["ratio"], 50)
        self.assertEqual(nb["strongest"][0]["q"], [0.0, 0.0, 0.5])
        # Compared: one mode per +q/-q pair, long waves out: 98 of the 215.
        self.assertEqual(nb["comparedModes"], 98)
        self.assertAlmostEqual(doubling["percentile"], 100 * 97 / 98, places=6)
        self.assertEqual(doubling["rank"], 1)
        # The other elements carry no order.
        for name in ("Ga", "Se"):
            self.assertTrue(all(point["ratio"] < 3.5 for point in _element(result, name)["zoneBoundary"]), name)

    def test_a_stronger_long_wave_does_not_hide_the_order(self):
        # Review 2026-10-12 #6: 55x the median at (0, 0, 1/2), but a long wave
        # and its -q partner above it put it at percentile 98.6 of all q.
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "long_wave.rmc6f"
            path.write_text(long_wave_config_text(), encoding="utf-8")
            ti = _element(displacement_modulations(load_site_displacements(path), "Ti"), "Ti")
        self.assertEqual(ti["strongest"][0]["q"], [1 / 6, 0.0, 0.0])
        order = next(point for point in ti["zoneBoundary"] if point["q"] == [0.0, 0.0, 0.5])
        self.assertGreater(order["ratio"], 50)
        self.assertEqual(order["rank"], 1)
        self.assertEqual(ti["comparedModes"], 98)

    def test_independent_displacements_have_no_order(self):
        result = displacement_modulations(_load("independent"))
        for entry in result["elements"]:
            self.assertTrue(all(point["ratio"] < 3.5 for point in entry["zoneBoundary"]), entry["element"])
            self.assertEqual(entry["missingCopies"], 0)

    def test_friedel_pairs_are_reported_once_as_their_smaller_index(self):
        result = displacement_modulations(_load("independent"), "Se", top=8)
        keys = set()
        for mode in _element(result, "Se")["strongest"]:
            index = tuple(round(value * 6) % 6 for value in mode["q"])
            partner = tuple((-value) % 6 for value in index)
            self.assertLessEqual(index, partner)
            self.assertNotIn(index, keys)
            keys.add(index)

    def test_an_odd_axis_leaves_points_untestable(self):
        result = displacement_modulations(_load("independent", supercell=(3, 4, 4)), "Nb")
        self.assertEqual([point["q"] for point in _element(result, "Nb")["zoneBoundary"]],
                         [[0.0, 0.0, 0.5], [0.0, 0.5, 0.0], [0.0, 0.5, 0.5]])
        self.assertEqual(result["untestable"], [[1, 0, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1]])

    def test_invalid_requests(self):
        sites = _load("independent", supercell=(4, 4, 4))
        with self.assertRaisesRegex(ValueError, "no site holds Xx"):
            displacement_modulations(sites, "Xx")
        with self.assertRaisesRegex(ValueError, "cell columns"):
            displacement_modulations(replace(sites, cell_indices=None))
        with self.assertRaisesRegex(ValueError, "too few wavevectors"):
            displacement_modulations(replace(sites, supercell=np.array([1, 1, 2])))


def _load_text(text: str):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "config.rmc6f"
        path.write_text(text, encoding="utf-8")
        return load_site_displacements(path)


class SignificanceTests(unittest.TestCase):
    """The calibrated significance against independent displacements (the tool's anchors, now the engine's)."""

    def test_zone_boundary_sigma(self):
        # Power 1 is the independent-displacement mean: nothing.
        self.assertLess(zone_boundary_sigma(1, 1, 7), 0)
        # One site, one point: P(chi2_1 >= 9) = 0.0027, 2.78 sigma; over seven points 1 - (1 - p)^7 = 0.019, 2.08 sigma.
        self.assertAlmostEqual(zone_boundary_sigma(9, 1, 1), 2.78, places=2)
        self.assertAlmostEqual(zone_boundary_sigma(9, 1, 7), 2.08, places=2)
        # More sites average the noise down: the same power is stronger evidence.
        self.assertGreater(zone_boundary_sigma(3, 16, 7), zone_boundary_sigma(3, 1, 7))
        self.assertGreater(zone_boundary_sigma(250, 1, 21), 8)
        for power in (None, 0.0, -1.0):
            self.assertIsNone(zone_boundary_sigma(power, 1, 7))

    def test_general_mode_sigma(self):
        # P >= chi2_2m / 2m. One site, one mode: P(chi2_2 >= 18) = e^-9 = 1.23e-4,
        # 3.67 sigma; corrected over 255 modes: 0.031, 1.87 sigma. Power 1 (the
        # independent mean) over 255 modes: 1 - (1 - e^-1)^255 rounds to 1,
        # sigma -inf, reported as None.
        self.assertLess(general_mode_sigma(1, 1, 7), 0)
        self.assertIsNone(general_mode_sigma(1, 1, 255))
        self.assertAlmostEqual(general_mode_sigma(9, 1, 1), 3.666, places=2)
        self.assertAlmostEqual(general_mode_sigma(9, 1, 255), 1.867, places=2)
        self.assertIsNone(general_mode_sigma(9, 1, 0))
        self.assertIsNone(general_mode_sigma(0.0, 1, 255))

    def test_infinite_sigma_is_reported_as_the_tool_showed_it(self):
        # An upper tail that underflows (sigma = +inf, shown '>8') is SIGMA_MAX,
        # which no finite double reaches; far below the null (sigma = -inf,
        # shown as no value) is None.
        self.assertEqual(zone_boundary_sigma(1e5, 4, 7), SIGMA_MAX)
        self.assertEqual(general_mode_sigma(1e5, 4, 99), SIGMA_MAX)
        self.assertLess(zone_boundary_sigma(650, 2, 1), SIGMA_MAX)  # p ~ 1e-282: finite, about 36 sigma
        self.assertIsNone(zone_boundary_sigma(0.01, 32, 21))
        self.assertEqual(SIGMA_MAX, 40)


class LongWavePredicateTests(unittest.TestCase):
    def test_one_predicate_on_every_grid_point(self):
        # The engine left out |q_i| <= LONG_WAVE_Q + 1e-12 from the compared
        # modes and the tool called a strongest mode a long wave at
        # |q_i| <= LONG_WAVE_Q: on a grid they never differ (k/N is 1/4 exactly
        # or at least 1/(4N) away), so one predicate replaces both. Each axis
        # is tested on its own, every N up to 256.
        for size in range(1, 257):
            for index in range(size):
                value = _folded(index, size)
                self.assertEqual(abs(value) <= LONG_WAVE_Q + 1e-12, abs(value) <= LONG_WAVE_Q, (index, size))
                self.assertEqual(is_long_wave([value, 0.0, 0.0]), abs(value) <= LONG_WAVE_Q, (index, size))
                # A zone-boundary component (0 or 1/2) by its index, as the
                # engine tests it, is the tool's 1e-9 test on the folded value.
                by_index = index == 0 or 2 * index == size
                by_value = abs(value) < 1e-9 or abs(abs(value) - 0.5) < 1e-9
                self.assertEqual(by_index, by_value, (index, size))
        self.assertTrue(is_long_wave([0.25, -0.25, 0.0]))
        self.assertFalse(is_long_wave([0.25, 0.0, 0.375]))


class ModulationVerdictTests(unittest.TestCase):
    """The engine's own verdict on configurations with known answers."""

    def check(self, result):
        """Rules and family echoed, every because in its closed set, ordered recomputable, strict JSON."""
        json.dumps(result, allow_nan=False)
        self.assertEqual(result["rules"], MODULATION_RULES)
        tested = [entry for entry in result["elements"] if entry["sites"]]
        self.assertEqual(result["family"], {
            "zoneBoundaryPoints": sum(len(entry["zoneBoundary"]) for entry in tested),
            "comparedModes": sum(entry["comparedModes"] for entry in tested),
        })
        for entry in result["elements"]:
            for name, finding in entry["verdict"].items():
                self.assertIn(finding["because"], MODULATION_BECAUSE[name], name)
            for point in entry.get("zoneBoundary", []) + entry.get("strongest", []):
                self.assertIn(point["kind"], MODE_KINDS)
                self.assertIn(point["because"], MODULATION_BECAUSE["point"])
                if point["bar"] is None:
                    self.assertEqual((point["sigma"], point["ordered"]), (None, None))
                    continue
                self.assertEqual(point["bar"], {"ratio": MODULATION_RATIO, "sigma": MODULATION_SIGMA})
                if point["because"] == "undefined_statistic":
                    self.assertIsNone(point["ordered"])
                else:
                    recomputed = point["ratio"] >= point["bar"]["ratio"] and point["sigma"] is not None \
                        and point["sigma"] >= point["bar"]["sigma"]
                    self.assertEqual(point["ordered"], recomputed)
                self.assertTrue(point["sigma"] is None or -SIGMA_MAX < point["sigma"] <= SIGMA_MAX)
            for point in entry.get("zoneBoundary", []):
                self.assertEqual(point["kind"], "zone_boundary")

    def element(self, result, name):
        return next(entry for entry in result["elements"] if entry["element"] == name)

    def test_a_planted_cell_doubling(self):
        result = displacement_modulations(_load("doubled"))
        self.check(result)
        nb = self.element(result, "Nb")
        self.assertEqual(nb["verdict"]["cellDoubling"], {
            "answer": True, "because": "ordered_point", "decidedBy": [[0.0, 0.0, 0.5]], "nearest": None,
        })
        self.assertEqual(nb["strongest"][0]["kind"], "zone_boundary")
        self.assertEqual(nb["strongest"][0]["because"], "zone_boundary")
        self.assertEqual(result["family"]["zoneBoundaryPoints"], 21)
        for name in ("Ga", "Se"):
            verdict = self.element(result, name)["verdict"]
            self.assertEqual(verdict["cellDoubling"]["answer"], False, name)
            self.assertEqual(verdict["cellDoubling"]["because"], "no_ordered_point", name)
            # The nearest miss is the point with the largest sigma.
            nearest = verdict["cellDoubling"]["nearest"]
            sigmas = [p["sigma"] for p in self.element(result, name)["zoneBoundary"] if p["sigma"] is not None]
            self.assertEqual(nearest["sigma"], max(sigmas))
            self.assertEqual(verdict["generalOrder"]["answer"], False, name)

    def test_order_under_a_stronger_long_wave(self):
        result = displacement_modulations(_load_text(long_wave_config_text()), "Ti")
        self.check(result)
        verdict = self.element(result, "Ti")["verdict"]
        self.assertEqual((verdict["cellDoubling"]["answer"], verdict["cellDoubling"]["decidedBy"]), (True, [[0.0, 0.0, 0.5]]))
        self.assertEqual(verdict["longWave"]["answer"], True)
        self.assertEqual(verdict["longWave"]["decidedBy"], [[1 / 6, 0.0, 0.0]])
        self.assertGreaterEqual(self.element(result, "Ti")["strongest"][0]["ratio"], LONG_WAVE_RATIO)
        # The family is the elements tested in this call: Ti alone, 7 points.
        self.assertEqual(result["family"]["zoneBoundaryPoints"], 7)
        self.assertEqual(displacement_modulations(_load_text(long_wave_config_text()))["family"]["zoneBoundaryPoints"], 14)

    def test_order_at_a_general_wavevector(self):
        # Ti displaced along z with q = (3/8, 0, 0): neither a zone-boundary
        # point nor a long wave -- correlation between cells, not a doubling.
        result = displacement_modulations(_load_text(wave_config_text()))
        self.check(result)
        ti = self.element(result, "Ti")
        self.assertEqual(ti["verdict"]["generalOrder"]["answer"], True)
        self.assertEqual(ti["verdict"]["generalOrder"]["decidedBy"], [[0.375, 0.0, 0.0]])
        self.assertEqual(ti["verdict"]["cellDoubling"]["answer"], False)
        self.assertEqual(ti["verdict"]["longWave"], {
            "answer": False, "because": "strongest_not_long_wave", "decidedBy": [], "nearest": None,
        })
        mode = ti["strongest"][0]
        self.assertEqual((mode["kind"], mode["ordered"], mode["because"]), ("general", True, "ordered"))
        self.assertGreater(mode["sigma"], MODULATION_SIGMA)
        o = self.element(result, "O")["verdict"]
        self.assertEqual((o["generalOrder"]["answer"], o["cellDoubling"]["answer"]), (False, False))

    def test_a_long_wave_and_general_order_together(self):
        result = displacement_modulations(_load_text(wave_config_text(long_amp=0.15, order_amp=0.12)))
        self.check(result)
        verdict = self.element(result, "Ti")["verdict"]
        self.assertEqual((verdict["longWave"]["answer"], verdict["longWave"]["decidedBy"]), (True, [[0.125, 0.0, 0.0]]))
        self.assertEqual((verdict["generalOrder"]["answer"], verdict["generalOrder"]["decidedBy"]), (True, [[0.375, 0.0, 0.0]]))
        self.assertEqual(verdict["cellDoubling"]["answer"], False)

    def test_independent_displacements_and_untestable_points(self):
        result = displacement_modulations(_load("independent", supercell=(3, 4, 4)))
        self.check(result)
        for entry in result["elements"]:
            self.assertEqual(entry["verdict"]["cellDoubling"]["answer"], False, entry["element"])
            self.assertNotEqual(entry["verdict"]["generalOrder"]["answer"], True, entry["element"])
        self.assertEqual(result["family"]["zoneBoundaryPoints"], 9)
        # No axis can be halved: the cell doubling cannot be tested.
        odd = displacement_modulations(_load("independent", supercell=(3, 3, 3)), "Nb")
        self.check(odd)
        self.assertEqual(odd["family"]["zoneBoundaryPoints"], 0)
        self.assertEqual(self.element(odd, "Nb")["verdict"]["cellDoubling"], {
            "answer": None, "because": "untestable", "decidedBy": [], "nearest": None,
        })

    def test_no_displacements(self):
        sites = _load("independent", supercell=(4, 4, 4))
        flat = replace(sites, displacements=np.zeros_like(sites.displacements))
        result = displacement_modulations(flat, "Nb")
        self.check(result)
        entry = self.element(result, "Nb")
        self.assertEqual(entry["sites"], 0)
        for finding in entry["verdict"].values():
            self.assertEqual((finding["answer"], finding["because"]), (None, "no_displacements"))
        self.assertEqual(result["family"], {"zoneBoundaryPoints": 0, "comparedModes": 0})


@unittest.skipUnless(DEMO_RMC6F.exists(), "the committed demo run is missing")
class DemoModulationTests(unittest.TestCase):
    def test_the_demo_has_no_cell_doubling(self):
        result = displacement_modulations(load_site_displacements(DEMO_RMC6F))
        for entry in result["elements"]:
            self.assertTrue(all(point["percentile"] < 99 for point in entry["zoneBoundary"]), entry["element"])

    def test_the_demo_verdict(self):
        # As the Copilot reports it: no cell doubling, and each element's
        # strongest mode a box-scale wave (no general mode among the five).
        result = displacement_modulations(load_site_displacements(DEMO_RMC6F))
        self.assertEqual(result["family"], {"zoneBoundaryPoints": 21, "comparedModes": 1323})
        for entry in result["elements"]:
            verdict = entry["verdict"]
            self.assertEqual((verdict["cellDoubling"]["answer"], verdict["cellDoubling"]["because"]),
                             (False, "no_ordered_point"), entry["element"])
            self.assertEqual(verdict["longWave"]["answer"], True, entry["element"])
            self.assertTrue(all(math.isfinite(p["ratio"]) for p in entry["zoneBoundary"]))


if __name__ == "__main__":
    unittest.main()
