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
from scipy.special import gammaincc

sys.path.insert(0, str(Path(__file__).resolve().parent))

import modulation_configs as configs  # noqa: E402
from cluster_configs import (  # noqa: E402
    POLARISED_CELL,
    cluster_config_text,
    long_wave_config_text,
    polarised_config_text,
    polarised_direction,
    wave_config_text,
)
from rmc_toolkits.clusters import evidence_clears  # noqa: E402
from rmc_toolkits.modulations import (  # noqa: E402
    AXIS_TIE_RTOL,
    DIFFUSE_CORRELATION,
    DIFFUSE_LAGS,
    DIFFUSE_PHASES,
    DIFFUSE_SIGMA,
    LONG_WAVE_Q,
    LONG_WAVE_RATIO,
    MODE_KINDS,
    MODULATION_BECAUSE,
    MODULATION_RATIO,
    MODULATION_RULES,
    MODULATION_SIGMA,
    SIGMA_MAX,
    _folded,
    diffuse_order,
    displacement_modulations,
    element_spectrum,
    general_mode_sigma,
    independent_cell_null,
    is_long_wave,
    mode_polarisation,
    null_tail,
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
        # The family is every element with displacements (F3): it holds every
        # element tested here, and each counts its points and compared modes.
        family = result["family"]
        self.assertEqual(set(family), {"zoneBoundaryPoints", "comparedModes", "diffuse", "elements"})
        self.assertEqual(family["elements"], sorted(family["elements"]))
        self.assertEqual(family["diffuse"], len(family["elements"]))
        self.assertTrue({entry["element"] for entry in tested} <= set(family["elements"]))
        for entry in tested:
            self.assertEqual(family["zoneBoundaryPoints"], len(entry["zoneBoundary"]) * len(family["elements"]))
            self.assertEqual(family["comparedModes"], entry["comparedModes"] * len(family["elements"]))
            # The diffuse finding: its answer is its one record's clears, recomputable.
            finding = entry["verdict"]["diffuseOrder"]
            diffuse = entry["diffuse"]
            self.assertEqual(len(diffuse["correlations"]), len(DIFFUSE_LAGS))
            if diffuse["df"] is None:
                self.assertEqual((finding["answer"], finding["because"], finding["evidence"]), (None, "too_few_modes", []))
            else:
                [record] = finding["evidence"]
                self.assertEqual(record["clears"], evidence_clears(record))
                self.assertEqual(record["bar"], {"z": DIFFUSE_SIGMA, "change": DIFFUSE_CORRELATION})
                if finding["answer"] is not None:
                    self.assertEqual(finding["answer"], record["clears"])
                if finding["answer"]:
                    self.assertIn(diffuse["phase"], DIFFUSE_PHASES)
                    self.assertEqual(finding["because"], f"correlated_{diffuse['phase']}")
            for point in entry["zoneBoundary"] + entry["strongest"]:
                pol = point["polarisation"]
                self.assertEqual(set(pol), {"direction", "uvw", "share"})
                self.assertAlmostEqual(float(np.linalg.norm(pol["direction"])), 1.0, places=12)
                self.assertAlmostEqual(max(abs(value) for value in pol["uvw"]), 1.0, places=12)
                self.assertTrue(1.0 / 3.0 - 1e-12 <= pol["share"] <= 1.0)
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
        # The family is every element of the configuration, whichever the call
        # names (F3): Ti alone counts O's 7 points too, as a call for both does.
        # (Until 2026-10-19 it was the elements of the call: 7.)
        self.assertEqual(result["family"], {"zoneBoundaryPoints": 14, "comparedModes": 196, "diffuse": 2, "elements": ["O", "Ti"]})
        self.assertEqual(displacement_modulations(_load_text(long_wave_config_text()))["family"], result["family"])

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
        self.assertEqual(result["family"], {"zoneBoundaryPoints": 0, "comparedModes": 0, "diffuse": 0, "elements": []})
        # An element that does not move is not in the family of the others either.
        still = replace(sites, displacements=np.where(
            (np.asarray(sites.elements)[sites.site_index] == "Nb")[:, None], 0.0, sites.displacements))
        others = displacement_modulations(still, "Ga")
        self.check(others)
        self.assertEqual(others["family"]["elements"], ["Ga", "Se"])
        self.assertEqual(displacement_modulations(still)["family"], others["family"])


def _angle_deg(a, b) -> float:
    """The angle between two axes (degrees; an axis has no sign)."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    return float(np.degrees(np.arccos(min(abs(a @ b) / (np.linalg.norm(a) * np.linalg.norm(b)), 1.0))))


#: The stated known answer of plan item 4.3: a planted linearly polarised mode
#: (0.1 A over 0.06 A of thermal noise, 2 sites x 128 cells) is found within
#: POLARISATION_DEG of its direction, with at least POLARISATION_SHARE of the
#: power along it. The scatter is the sampling error sigma / (delta sqrt(m n))
#: rad per component (0.0375 rad here: 5 degrees is its Rayleigh 93rd
#: percentile), which the sweep below checks.
POLARISATION_DEG = 5.0
POLARISATION_SHARE = 0.95


class PolarisationTests(unittest.TestCase):
    """Plan item 4.3: the displacement direction of a mode (mode_polarisation)."""

    def test_the_leading_axis_as_a_crystal_direction(self):
        # A rank-1 tensor along the direct-lattice [1 0 1] of a monoclinic
        # cell: uvw is [1 0 1], the direction its Cartesian vector a + c,
        # which is 14 degrees from the Cartesian (1, 0, 1).
        along = polarised_direction((1, 0, 1))
        self.assertGreater(_angle_deg(along, [1.0, 0.0, 1.0]), 13.0)
        for sign in (1.0, -1.0):
            pol = mode_polarisation(np.outer(sign * along, sign * along) * 2.5, POLARISED_CELL)
            # Either sign gives the sense that makes the closest crystal axis (c) acute.
            np.testing.assert_allclose(pol["direction"], along, atol=1e-12)
            np.testing.assert_allclose(pol["uvw"], [1.0, 0.0, 1.0], atol=1e-12)
            self.assertAlmostEqual(pol["share"], 1.0, places=12)

    def test_share_and_the_imaginary_part(self):
        cubic = np.eye(3) * 4.0
        pol = mode_polarisation(np.diag([0.2, 1.0, 0.05]), cubic)
        np.testing.assert_allclose(pol["uvw"], [0.0, 1.0, 0.0], atol=1e-12)
        self.assertAlmostEqual(pol["share"], 1.0 / 1.25)
        # A circular mode (a + ib with a, b orthogonal): its real part is
        # aa^T + bb^T, no axis is preferred in the plane (share 1/2).
        a, b = np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])
        amplitude = a + 1j * b
        self.assertAlmostEqual(mode_polarisation(np.outer(amplitude, amplitude.conj()).real, cubic)["share"], 0.5)
        self.assertIsNone(mode_polarisation(np.zeros((3, 3)), cubic))

    def test_the_sense_of_an_axis_equally_close_to_two(self):
        # |cos| to a and b tie (within AXIS_TIE_RTOL): the first, a, is made acute.
        cubic = np.eye(3) * 4.0
        for direction, want in (([1.0, 1.0, 0.0], [1.0, 1.0, 0.0]), ([-1.0, 1.0, 0.0], [1.0, -1.0, 0.0])):
            d = np.array(direction) / np.sqrt(2.0)
            np.testing.assert_allclose(mode_polarisation(np.outer(d, d), cubic)["uvw"], want, atol=1e-12)
        self.assertEqual(AXIS_TIE_RTOL, 1e-9)

    def test_a_planted_mode_in_a_monoclinic_cell(self):
        # Both Zr sites move +-0.1 A along [1 0 1] in alternate cells along a
        # (q = (1/2, 0, 0)), the second in antiphase with the first.
        result = displacement_modulations(_load_text(polarised_config_text()))
        zr = next(entry for entry in result["elements"] if entry["element"] == "Zr")
        point = next(point for point in zr["zoneBoundary"] if point["q"] == [0.5, 0.0, 0.0])
        self.assertTrue(point["ordered"])
        pol = point["polarisation"]
        want = polarised_direction((1, 0, 1))
        self.assertLess(_angle_deg(pol["direction"], want), POLARISATION_DEG)
        self.assertLess(_angle_deg(POLARISED_CELL.T @ np.array(pol["uvw"]), want), POLARISATION_DEG)
        self.assertGreater(pol["share"], POLARISATION_SHARE)
        self.assertGreater(float(np.dot(pol["direction"], want)), 0)  # the acute sense
        # The strongest mode is that point, with the same polarisation.
        self.assertEqual(zr["strongest"][0]["q"], [0.5, 0.0, 0.0])
        self.assertEqual(zr["strongest"][0]["polarisation"], pol)

    def test_the_scatter_is_the_sampling_error(self):
        # Over 60 seeds the polarisation is unbiased and scatters as the
        # sampling error: per component sigma / (delta sqrt(m n)) rad, so the
        # mean squared angle is twice its square (within 25 %: the mean of
        # 60 draws of a chi^2_2 has a relative sd of 13 %), none beyond 10
        # degrees, and the share always above POLARISATION_SHARE.
        sd = 0.06 / (0.1 * np.sqrt(2 * 128))
        want = polarised_direction((1, 0, 1))
        angles, shares = [], []
        for seed in range(1, 61):
            result = displacement_modulations(_load_text(polarised_config_text(seed=seed)), "Zr")
            pol = next(p for p in result["elements"][0]["zoneBoundary"] if p["q"] == [0.5, 0.0, 0.0])["polarisation"]
            angles.append(np.radians(_angle_deg(pol["direction"], want)))
            shares.append(pol["share"])
        ratio = float(np.mean(np.square(angles))) / (2 * sd ** 2)
        self.assertTrue(0.75 < ratio < 1.25, ratio)
        self.assertLess(np.degrees(max(angles)), 10.0)
        self.assertGreater(min(shares), POLARISATION_SHARE)

    def test_the_dev_style_orders(self):
        # The existing known answers: the doubled Nb along z, the long-wave
        # repro (order along z at (0, 0, 1/2), the long wave along y), and the
        # general-q order along z under a long wave along y.
        nb = _element(displacement_modulations(_load("doubled")), "Nb")
        doubling = next(point for point in nb["zoneBoundary"] if point["q"] == [0.0, 0.0, 0.5])
        self.assertLess(_angle_deg(doubling["polarisation"]["direction"], [0, 0, 1]), POLARISATION_DEG)
        self.assertGreater(doubling["polarisation"]["share"], POLARISATION_SHARE)
        ti = _element(displacement_modulations(_load_text(long_wave_config_text()), "Ti"), "Ti")
        self.assertLess(_angle_deg(ti["strongest"][0]["polarisation"]["direction"], [0, 1, 0]), POLARISATION_DEG)
        order = next(point for point in ti["zoneBoundary"] if point["q"] == [0.0, 0.0, 0.5])
        self.assertLess(_angle_deg(order["polarisation"]["direction"], [0, 0, 1]), POLARISATION_DEG)
        both = _element(displacement_modulations(_load_text(wave_config_text(long_amp=0.15, order_amp=0.12))), "Ti")
        by_kind = {mode["kind"]: mode for mode in both["strongest"][:2]}
        self.assertLess(_angle_deg(by_kind["long_wave"]["polarisation"]["direction"], [0, 1, 0]), POLARISATION_DEG)
        self.assertLess(_angle_deg(by_kind["general"]["polarisation"]["direction"], [0, 0, 1]), POLARISATION_DEG)


class FamilyTests(unittest.TestCase):
    """F3 (2026-10-19): the Sidak family is every element of the configuration, whichever the call names."""

    def test_a_point_just_below_the_bar_stays_below_it_when_its_element_is_named(self):
        # A weak doubling (Nb +-0.02 A along z in alternate c layers, under
        # 0.06 A of noise): its (0, 0, 1/2) point sits at 2.99 sigma over the
        # 21 zone-boundary points of Ga, Nb and Se (against the
        # independent-cell null, F1). Naming Nb used to shrink the family to
        # its own 7 points, which takes it to 3.31 sigma: "ordered". It is the
        # same point, the same data and the same verdict either way now.
        # (Under the null before F1 the same known answer was Nb +-0.03 A,
        # seed 26: 2.98 and 3.30 sigma.)
        sites = _load_text(cluster_config_text("doubled", supercell=(4, 4, 4), delta=0.06, seed=33))
        every = displacement_modulations(sites)
        named = displacement_modulations(sites, "Nb")
        self.assertEqual(every["family"], {"zoneBoundaryPoints": 21, "comparedModes": 66, "diffuse": 3, "elements": ["Ga", "Nb", "Se"]})
        self.assertEqual(named["family"], every["family"])
        point = next(p for p in _element(every, "Nb")["zoneBoundary"] if p["q"] == [0.0, 0.0, 0.5])
        again = next(p for p in _element(named, "Nb")["zoneBoundary"] if p["q"] == [0.0, 0.0, 0.5])
        self.assertEqual(again, point)
        self.assertTrue(2.95 < point["sigma"] < MODULATION_SIGMA)
        self.assertGreater(point["ratio"], MODULATION_RATIO)
        self.assertEqual((point["ordered"], point["because"]), (False, "below_sigma"))
        for result in (every, named):
            self.assertEqual(_element(result, "Nb")["verdict"]["cellDoubling"]["answer"], False)
        # The old family (Nb's 7 points) would have called it ordered.
        power = point["ratio"] * _element(named, "Nb")["medianPower"]
        self.assertGreaterEqual(zone_boundary_sigma(power, independent_cell_null(sites, "Nb"), 7), MODULATION_SIGMA)


@unittest.skipUnless(DEMO_RMC6F.exists(), "the committed demo run is missing")
class DemoModulationTests(unittest.TestCase):
    def test_the_demo_has_no_cell_doubling(self):
        result = displacement_modulations(load_site_displacements(DEMO_RMC6F))
        for entry in result["elements"]:
            self.assertTrue(all(point["percentile"] < 99 for point in entry["zoneBoundary"]), entry["element"])

    def test_the_demo_verdict(self):
        # As the Copilot reports it: no cell doubling, and each element's
        # strongest mode a box-scale wave (no general mode among the five).
        sites = load_site_displacements(DEMO_RMC6F)
        result = displacement_modulations(sites)
        self.assertEqual(result["family"], {"zoneBoundaryPoints": 21, "comparedModes": 1323, "diffuse": 3, "elements": ["Ga", "Se", "Ta"]})
        # F3: a call naming one element corrects over the same family.
        self.assertEqual(displacement_modulations(sites, "Ta")["family"], result["family"])
        for entry in result["elements"]:
            verdict = entry["verdict"]
            self.assertEqual((verdict["cellDoubling"]["answer"], verdict["cellDoubling"]["because"]),
                             (False, "no_ordered_point"), entry["element"])
            self.assertEqual(verdict["longWave"]["answer"], True, entry["element"])
            self.assertTrue(all(math.isfinite(p["ratio"]) for p in entry["zoneBoundary"]))
            # Plan item 4.4: the power beyond the long waves is spread toward
            # the zone centre -- neighbouring cells move together, in phase
            # (the tail of RMC's box-scale waves), never in antiphase.
            diffuse = entry["diffuse"]
            self.assertEqual(diffuse["phase"], "in_phase", entry["element"])
            self.assertGreater(diffuse["sigma"], 5, entry["element"])
            self.assertTrue(all(c > 0 for c in diffuse["correlations"][:3]), entry["element"])
        se = _element(result, "Se")
        self.assertGreater(se["diffuse"]["correlation"], 0.15)
        self.assertEqual(se["verdict"]["diffuseOrder"]["because"], "correlated_in_phase")


class NullTests(unittest.TestCase):
    """F1 (2026-10-19): the independent-cell null keeps every correlation inside a cell."""

    def test_the_weights_sum_to_the_mean_power(self):
        # Parseval with each site's mean subtracted: the power averages
        # N/(N - 1) over q != 0, and so do the null's weights.
        for make in configs.NULL_FAMILIES.values():
            sites = configs.site_displacements(make(5))
            n_q = int(np.prod(sites.supercell))
            for name in sorted(set(sites.elements)):
                power = element_spectrum(sites, name)[0].ravel()
                self.assertAlmostEqual(float(power[1:].mean()), n_q / (n_q - 1), places=12)
                self.assertAlmostEqual(float(independent_cell_null(sites, name).sum()), n_q / (n_q - 1), places=12)

    def test_a_shared_cell_variable_leaves_fewer_independent_components(self):
        # Independent isotropic sites: about three components per site. The
        # six F of a Jahn-Teller octahedron share its long axis: about three
        # in all, where the null before F1 assumed six (one per site).
        thermal = displacement_modulations(configs.site_displacements(configs.thermal(5)))
        self.assertAlmostEqual(_element(thermal, "Ti")["nullDof"], 3.0, delta=0.1)
        jt = displacement_modulations(configs.site_displacements(configs.jahn_teller(5)))
        self.assertEqual(_element(jt, "F")["sites"], 6)
        self.assertLess(_element(jt, "F")["nullDof"], 4.0)

    def test_equal_weights_are_the_exact_gamma(self):
        # An int m is m equal weights 1/m, the bound the engine used before F1.
        self.assertAlmostEqual(null_tail(9, 1, True), 0.0026997960632601866, places=15)
        self.assertAlmostEqual(null_tail(9, [0.5, 0.5], False), math.exp(-18) * 19, places=15)
        self.assertEqual(null_tail(9, 2, False), null_tail(9, [0.5, 0.5], False))
        self.assertIsNone(null_tail(0.0, 3, True))
        self.assertIsNone(null_tail(1.0, [], True))

    def test_the_saddlepoint_tail_is_never_anticonservative(self):
        # Against the exact tail (the Moschopoulos series of a sum of gammas)
        # for weight patterns the engine meets -- one dominant component over
        # several small ones, an anisotropic site, two equal ones over many --
        # the saddlepoint lies within -1 % and +8 % of it from p = 1e-2 to 1e-6.
        # (A three-moment gamma fell to 0.36 of it at 1e-6.)
        def exact_tail(weights, real):
            """P(sum lambda_k Y_k > x) as a function of x: a gamma mixture whose
            coefficients (they do not depend on x) are computed once."""
            lam = np.asarray(weights, dtype=float)
            alpha = np.full(lam.size, 0.5 if real else 1.0)
            beta = lam * (2.0 if real else 1.0)
            b1 = beta.min()
            rho = alpha.sum()
            ratio = 1.0 - b1 / beta
            terms = 2000  # the slowest pattern's coefficients fall below 1e-12 of the sum by about 900
            k = np.arange(1, terms)
            weighted = np.array([0.0] + [float(np.sum(alpha * ratio ** j)) for j in k])  # j gamma_j
            delta = np.zeros(terms)
            delta[0] = 1.0
            for j in k:
                delta[j] = float(weighted[1:j + 1] @ delta[j - 1::-1]) / j
            scale = math.exp(float(np.sum(alpha * np.log(b1 / beta))))
            orders = rho + np.arange(terms)
            return lambda x: scale * float(np.sum(delta * gammaincc(orders, x / b1)))

        patterns = [[0.429] + [0.0713] * 8, [0.6, 0.25, 0.15], [0.39, 0.39] + [0.01375] * 16]
        for weights in patterns:
            for real in (True, False):
                exact = exact_tail(weights, real)
                for target in (1e-2, 1e-4, 1e-6):
                    lo, hi = 0.0, 100.0
                    for _ in range(80):
                        mid = 0.5 * (lo + hi)
                        lo, hi = (mid, hi) if exact(mid) > target else (lo, mid)
                    ratio = null_tail(0.5 * (lo + hi), weights, real) / target
                    self.assertTrue(0.99 < ratio < 1.08, (weights[:2], real, target, ratio))

    def test_shared_variables_no_longer_read_as_order(self):
        # F1's repro: an F mode of Jahn-Teller octahedra and a W mode of
        # elongated octahedral clusters, five to six times the median, which
        # the null of independent sites put above 3 sigma (3.66 and 4.46): the
        # sites of each octahedron share one variable, so the element's power
        # has about three components, not six. Against the independent-cell
        # null both are below 2 sigma, and nothing is ordered.
        for make, seed, name, want_q in ((configs.jahn_teller, 11, "F", (0.5, 0.5, 1 / 3)),
                                         (configs.octahedral, 4, "W", (0.0, 0.375, 0.25))):
            sites = configs.site_displacements(make(seed))
            result = displacement_modulations(sites)
            entry = _element(result, name)
            mode = next(m for m in entry["strongest"] if np.allclose(m["q"], want_q))
            self.assertEqual(mode["kind"], "general")
            self.assertGreater(mode["ratio"], 4.5)
            old = general_mode_sigma(mode["ratio"] * entry["medianPower"], entry["sites"], result["family"]["comparedModes"])
            self.assertGreater(old, 3.5)
            self.assertLess(mode["sigma"], 2.0)
            self.assertEqual((mode["ordered"], mode["because"]), (False, "below_sigma"))
            for other in result["elements"]:
                for finding in ("cellDoubling", "generalOrder", "diffuseOrder"):
                    self.assertNotEqual(other["verdict"][finding]["answer"], True, (make.__name__, other["element"], finding))

    def test_independent_isotropic_sites_are_no_longer_held_to_the_bound(self):
        # One isotropic site has three components: the bound chi2_1 the null
        # took before F1 read a weak doubling (B +-0.018 A along z at
        # (1/2, 1/2, 0) over 0.06 A of thermal motion, the weak order the
        # validation's judgment pair had until protocol 2026-10-19.19 re-set
        # it) near 3 sigma; the null reads it above 5.
        result = displacement_modulations(configs.site_displacements(configs.zone_mode(3, delta=0.018, q=(0.5, 0.5, 0))))
        ti = _element(result, "Ti")
        point = next(p for p in ti["zoneBoundary"] if p["q"] == [0.5, 0.5, 0.0])
        old = zone_boundary_sigma(point["ratio"] * ti["medianPower"], 1, result["family"]["zoneBoundaryPoints"])
        self.assertLess(old, point["sigma"] - 1.5)
        self.assertTrue(point["ordered"])


class DiffuseOrderTests(unittest.TestCase):
    """Plan item 4.4: diffuse order between cells, spread over many q."""

    def test_short_range_order_in_antiphase_with_no_single_mode(self):
        # B displacements in neighbouring cells along a correlated by -0.27
        # (antiphase): the power lies on the whole sheet q_a = 1/2, no mode
        # stands out (no cell doubling, no general order), and the diffuse
        # test finds it, along [1 0 0], in antiphase.
        result = displacement_modulations(configs.site_displacements(configs.diffuse(1)))
        nb = _element(result, "Nb")
        verdict = nb["verdict"]
        self.assertEqual((verdict["cellDoubling"]["answer"], verdict["generalOrder"]["answer"]), (False, False))
        finding = verdict["diffuseOrder"]
        self.assertEqual((finding["answer"], finding["because"]), (True, "correlated_antiphase"))
        diffuse = nb["diffuse"]
        self.assertEqual((diffuse["lag"], diffuse["phase"], diffuse["df"], diffuse["modes"]), ([1, 0, 0], "antiphase", 13, 190))
        self.assertLess(diffuse["correlation"], -0.15)
        self.assertGreater(diffuse["sigma"], 5)
        [record] = finding["evidence"]
        self.assertEqual((record["observed"], record["z"], record["clears"]), (abs(diffuse["correlation"]), diffuse["sigma"], True))
        self.assertEqual(record["side"], "above")
        # The other elements move independently: no diffuse order there.
        for name in ("Ag", "O"):
            self.assertNotEqual(_element(result, name)["verdict"]["diffuseOrder"]["answer"], True, name)

    def test_in_phase_correlation(self):
        result = displacement_modulations(configs.site_displacements(configs.diffuse(1, sign=1.0, axis=2)))
        nb = _element(result, "Nb")
        self.assertEqual(nb["verdict"]["diffuseOrder"]["because"], "correlated_in_phase")
        self.assertEqual((nb["diffuse"]["lag"], nb["diffuse"]["phase"]), ([0, 0, 1], "in_phase"))
        self.assertGreater(nb["diffuse"]["correlation"], 0.15)

    def test_one_strong_mode_alone_is_not_diffuse_order(self):
        # A cell doubling (dev-zone-boundary's kind: B +-0.1 A at (0, 0, 1/2))
        # is the zone-boundary test's: its point is not among the diffuse
        # test's modes. An ordered general mode is the general test's: it is
        # left out of them. Ranks bound what any other one mode adds.
        doubled = displacement_modulations(configs.site_displacements(configs.zone_mode(7)))
        ti = _element(doubled, "Ti")
        self.assertTrue(ti["verdict"]["cellDoubling"]["answer"])
        self.assertEqual(ti["verdict"]["diffuseOrder"]["answer"], False)
        self.assertLess(ti["diffuse"]["sigma"], 2)
        wave = displacement_modulations(_load_text(wave_config_text()))
        ti = _element(wave, "Ti")
        self.assertTrue(ti["verdict"]["generalOrder"]["answer"])
        self.assertEqual(ti["verdict"]["diffuseOrder"]["answer"], False)
        self.assertEqual(ti["diffuse"]["modes"], _element(wave, "O")["diffuse"]["modes"] - 1)

    def test_the_statistic_reads_ranks(self):
        # Under independent cells the power at every compared mode is
        # identically distributed: the ranks are exchangeable, so the
        # statistic is the same for any increasing function of the power, and
        # one mode a million times stronger moves it no further than the top rank.
        rng = np.random.default_rng(4)
        q = np.array([(h / 8, k / 8, l / 8) for h in range(8) for k in range(-3, 5) for l in range(-3, 5)][:150], dtype=float)
        power = rng.exponential(1.0, q.shape[0])
        base = diffuse_order(power, q, 3)
        self.assertAlmostEqual(diffuse_order(np.log(power) + 10.0, q, 3)["chiSquare"], base["chiSquare"], places=9)
        spiked = power.copy()
        spiked[int(np.argmax(power))] = 1e6
        self.assertAlmostEqual(diffuse_order(spiked, q, 3)["chiSquare"], base["chiSquare"], places=9)

    def test_too_few_modes(self):
        # A 4 x 4 x 4 grid leaves 15 compared modes outside the long waves and
        # zone-boundary points: fewer than two per term of the regression.
        result = displacement_modulations(_load("independent", supercell=(4, 4, 4)))
        for entry in result["elements"]:
            self.assertEqual(entry["diffuse"]["modes"], 15)
            self.assertEqual(entry["verdict"]["diffuseOrder"], {"answer": None, "because": "too_few_modes", "evidence": []})


if __name__ == "__main__":
    unittest.main()
