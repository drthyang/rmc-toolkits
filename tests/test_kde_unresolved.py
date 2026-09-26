# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""A map whose kernel falls between all grid nodes is flagged, not contoured.

A needle kernel (sigma_minor far below the node spacing) can miss every grid
node: the map then holds only Gaussian tails, 1e-14 on the GaNb4Se8 AVERAGE Nb
layer, and a per-slice colour scale and eight contour levels would stretch
that round-off into a picture. The grid-summed linear density ("mass", about 1
for a resolved map, >= 0.13 for every drawn case of the parity fixture) tells
the two apart: below ``UNRESOLVED_MASS_LIMIT`` the engine attaches the
``unresolved`` warning and draws no contours, in both scales. The worker's twin
is ``workers/__tests__/unresolvedMap.test.js``.
"""

from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from rmc_toolkits.kde import (
    KDE_WARNINGS,
    KERNEL_MIN_SIGMA,
    UNRESOLVED_MASS_LIMIT,
    load_unit_cell_positions,
    oriented_kde_slice,
)

ROOT = Path(__file__).resolve().parents[1]
AVERAGE_5K = ROOT / "data" / "5K_try1" / "GaNb4Se8_5KAVERAGE.rmc6f"
C_SLICE = {"normal": np.array([0.0, 0.0, 1.0]), "u_axis": np.array([1.0, 0.0, 0.0]), "v_axis": np.array([0.0, 1.0, 0.0])}


def two_site_needle(spread: float) -> np.ndarray:
    """Two sites on y = 0.5, midway between the nodes of a 32-point grid.

    Each site is four atoms at (x0 +- 0.01, 0.5 +- spread): C is diagonal with
    a minor sigma of bw * spread, and every atom sits 1/62 - spread from the
    nearest node row.
    """
    return np.array(
        [(x0 + dx, 0.5 + dy, 0.5) for x0 in (0.2, 0.8) for dx in (-0.01, 0.01) for dy in (-spread, spread)]
    )


def grid_mass(result: dict) -> float:
    x0, x1, y0, y1 = result["extent"]
    step = (x1 - x0) / (result["grid"] - 1) * (y1 - y0) / (result["grid"] - 1)
    density = np.asarray(result["density"], dtype=float)
    if result["log"]:
        density = 10.0**density - 1e-12
    return float(density.sum() * step)


class UnresolvedMapTests(unittest.TestCase):
    def codes(self, result):
        return [warning["code"] for warning in result["warnings"]]

    def test_kernel_between_the_nodes_is_flagged_and_not_contoured(self):
        for log in (False, True):
            with self.subTest(log=log):
                result = oriented_kde_slice(two_site_needle(0.01), 0.5, 0.08, bw=0.076, grid=32, log=log, **C_SLICE)
                self.assertIsNone(result["message"])
                self.assertIsNotNone(result["kernel"])
                self.assertEqual(self.codes(result), ["subgrid", "unresolved"])
                self.assertEqual(result["warnings"][1]["message"], KDE_WARNINGS["unresolved"])
                self.assertEqual(result["contours"], [])
                # The density itself is still reported: tails, not zeros.
                self.assertGreater(result["vmax"], result["vmin"])
                self.assertLess(grid_mass(result), UNRESOLVED_MASS_LIMIT)

    def test_fully_underflowed_map_is_flagged_too(self):
        result = oriented_kde_slice(two_site_needle(0.004), 0.5, 0.08, bw=0.03, grid=32, **C_SLICE)
        self.assertEqual(result["vmax"], 0.0)
        self.assertEqual(self.codes(result), ["subgrid", "unresolved"])
        self.assertEqual(result["contours"], [])

    def test_a_kernel_below_the_evaluable_floor_is_a_flagged_zero_map(self):
        # bw = 1e-30: scipy's whitening of the absolute coordinates cannot
        # resolve the kernel (its value at an atom comes out 0 and the
        # _FixedCovarianceKDE self-check used to decline with the SciPy
        # "engine" message). bw = 1e-200: det H underflows as well. Both are
        # the zero map with the kernel summary, flagged -- never NaN, never a
        # decline that blames the SciPy release. The worker returns the same.
        for bw in (1e-30, 1e-200):
            with self.subTest(bw=bw):
                result = oriented_kde_slice(two_site_needle(0.01), 0.5, 0.08, bw=bw, grid=32, **C_SLICE)
                density = np.asarray(result["density"], dtype=float)
                self.assertTrue(np.all(density == 0.0))
                self.assertIsNone(result["message"])
                self.assertIsNotNone(result["kernel"])
                self.assertLess(result["kernel"]["sigmaMinor"], KERNEL_MIN_SIGMA)
                self.assertEqual(result["fitCount"], 8)
                self.assertEqual(self.codes(result), ["subgrid", "unresolved"])
                self.assertEqual(result["contours"], [])

    def test_the_floor_leaves_a_narrow_but_evaluable_kernel_alone(self):
        # sigma_minor = 1e-6 * 0.0115: far below the grid, above the floor.
        # An atom placed on a node shows the evaluated spike (scipy's sum).
        points = two_site_needle(0.01)
        points[0, :2] = (0.0, 0.0)
        result = oriented_kde_slice(points, 0.5, 0.08, bw=1e-6, grid=32, **C_SLICE)
        self.assertGreater(result["kernel"]["sigmaMinor"], KERNEL_MIN_SIGMA)
        self.assertGreater(result["vmax"], 1e9)
        self.assertTrue(np.all(np.isfinite(np.asarray(result["density"], dtype=float))))

    def test_the_same_layer_resolved_by_a_wider_kernel_is_not_flagged(self):
        # bw = 0.5 puts sigma_minor at 5e-3 against a node spacing of 0.032.
        result = oriented_kde_slice(two_site_needle(0.01), 0.5, 0.08, bw=0.5, grid=32, **C_SLICE)
        self.assertNotIn("unresolved", self.codes(result))
        self.assertGreater(grid_mass(result), 0.1)
        self.assertEqual(len(result["contours"]), 8)

    def test_declined_slab_carries_no_warning(self):
        result = oriented_kde_slice(two_site_needle(0.01), 0.5, 0.08, bw=0.0, grid=32, **C_SLICE)
        self.assertIsNotNone(result["message"])
        self.assertEqual(result["warnings"], [])

    @unittest.skipUnless(AVERAGE_5K.exists(), "GaNb4Se8 5K run not present in data/ (gitignored)")
    def test_average_nb_layer_needle_is_flagged(self):
        positions = load_unit_cell_positions(AVERAGE_5K, element="Nb").fractional_positions
        result = oriented_kde_slice(positions, 0.15, 0.08, bw=0.03, grid=120, log=True, **C_SLICE)
        self.assertLess(result["kernel"]["sigmaMinor"], 1e-5)
        self.assertEqual(self.codes(result), ["subgrid", "unresolved"])
        self.assertEqual(result["contours"], [])


if __name__ == "__main__":
    unittest.main()
