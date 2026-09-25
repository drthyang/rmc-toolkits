# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""A kernel narrower than half a grid step is flagged, not drawn silently.

``H = bw**2 * Cov(slab atoms)`` collapses for a slab of one or two sites (the
minor axis becomes bw times the thermal spread), and the map then samples a
needle far below the grid spacing: its peak and integral follow the grid, not
the atoms. Both runtimes attach the same ``subgrid`` warning
(``workers/__tests__/kernelDiagnostics.test.js`` is the worker's twin).
"""

from __future__ import annotations

import unittest

import numpy as np

from rmc_toolkits.kde import KDE_WARNINGS, oriented_kde_slice

C_SLICE = {
    "normal": np.array([0.0, 0.0, 1.0]),
    "u_axis": np.array([1.0, 0.0, 0.0]),
    "v_axis": np.array([0.0, 1.0, 0.0]),
}


class KernelDiagnosticsTests(unittest.TestCase):
    def test_two_site_layer_needle_is_flagged(self):
        rng = np.random.default_rng(5)
        sites = [np.array([0.25, 0.75, 0.5]), np.array([0.75, 0.25, 0.5])]
        points = np.vstack([site + np.array([0.006, 0.006, 0.0]) * rng.standard_normal((300, 3)) for site in sites])
        result = oriented_kde_slice(points, center=0.5, thickness=0.08, bw=0.03, grid=120, n_levels=0, **C_SLICE)
        self.assertLess(result["kernel"]["sigmaMinor"], 0.5 / 119)
        self.assertEqual(result["warnings"], [{"code": "subgrid", "message": KDE_WARNINGS["subgrid"]}])

    def test_cell_filling_slab_is_quiet_at_grid_120_but_not_at_grid_16(self):
        rng = np.random.default_rng(6)
        points = np.column_stack([rng.random(3000), rng.random(3000), np.full(3000, 0.5)])
        common = {"center": 0.5, "thickness": 0.08, "bw": 0.03, "n_levels": 0, **C_SLICE}
        self.assertEqual(oriented_kde_slice(points, grid=120, **common)["warnings"], [])
        coarse = oriented_kde_slice(points, grid=16, **common)
        self.assertEqual([warning["code"] for warning in coarse["warnings"]], ["subgrid"])

    def test_declined_slab_has_no_warnings(self):
        result = oriented_kde_slice(np.empty((0, 3)), center=0.5, thickness=0.08, grid=16, **C_SLICE)
        self.assertEqual(result["warnings"], [])
        self.assertIsNone(result["kernel"])


if __name__ == "__main__":
    unittest.main()
