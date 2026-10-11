# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The displacement-modulation (superstructure) engine on synthetic configurations."""

from __future__ import annotations

import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cluster_configs import cluster_config_text  # noqa: E402
from rmc_toolkits.modulations import displacement_modulations, zone_boundary_points  # noqa: E402
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
        self.assertAlmostEqual(doubling["percentile"], 100 * 214 / 215, places=6)
        # The other elements carry no order.
        for name in ("Ga", "Se"):
            self.assertTrue(all(point["ratio"] < 3.5 for point in _element(result, name)["zoneBoundary"]), name)

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


@unittest.skipUnless(DEMO_RMC6F.exists(), "the committed demo run is missing")
class DemoModulationTests(unittest.TestCase):
    def test_the_demo_has_no_cell_doubling(self):
        result = displacement_modulations(load_site_displacements(DEMO_RMC6F))
        for entry in result["elements"]:
            self.assertTrue(all(point["percentile"] < 99 for point in entry["zoneBoundary"]), entry["element"])


if __name__ == "__main__":
    unittest.main()
