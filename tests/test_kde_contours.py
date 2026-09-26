# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Log scale changes the contour levels, never whether there are contours.

The Python guard used to test ``density.max() > 0`` after the log10 transform,
so an oblique slice of a disordered model at bw >= 0.1 -- peak below 1 per unit
fractional area -- lost every contour on the SciPy path while the browser drew
all eight. Positivity is now judged on the linear density in both runtimes
(``workers/__tests__/logContours.test.js`` is the worker's twin).
"""

from __future__ import annotations

import unittest

import numpy as np

from rmc_toolkits.kde import oriented_kde_slice


def kronecker_points(count: int) -> np.ndarray:
    """A deterministic, uniform-looking cloud that JS reproduces bit for bit."""
    index = np.arange(1, count + 1, dtype=float)
    return np.column_stack([(index * 0.7548776662466927) % 1.0, (index * 0.5698402909980532) % 1.0, (index * 0.6180339887498949) % 1.0])


class LogContourTests(unittest.TestCase):
    def test_oblique_disordered_slice_keeps_its_contours_in_log_mode(self):
        points = kronecker_points(8000)
        common = {"center": 0.5, "thickness": 0.08, "normal": np.array([1.0, 1.0, 1.0]), "bw": 0.1, "grid": 48}
        log_map = oriented_kde_slice(points, log=True, **common)
        linear_map = oriented_kde_slice(points, log=False, **common)
        self.assertLess(log_map["vmax"], 0.0)  # the peak is below 1 per unit area
        self.assertEqual(len(linear_map["contours"]), 8)
        self.assertEqual(len(log_map["contours"]), 8)

    def test_declined_slab_has_no_contours_in_either_mode(self):
        for log in (False, True):
            result = oriented_kde_slice(
                kronecker_points(8000), center=0.5, thickness=0.08, normal=np.array([0.0, 0.0, 1.0]), bw=0.0, log=log
            )
            self.assertEqual(result["contours"], [])


if __name__ == "__main__":
    unittest.main()
