# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""kde_slice either draws gaussian_kde's exact kernel or declines with a reason.

The browser worker ports these rules one for one (``KDE_MESSAGES`` in
``localKdeWorker.js``); ``kdeParity.test.js`` checks the two against each other.
"""

from __future__ import annotations

import unittest

import numpy as np
from scipy.stats import gaussian_kde

from rmc_toolkits.kde import COVARIANCE_CONDITION_LIMIT, KDE_MESSAGES, kde_slice


def _slice(points, **kwargs):
    options = {"z_center": 0.5, "dz": 0.1, "xlim": (0.0, 1.0), "ylim": (0.0, 1.0), "grid": 16, "n_levels": 0}
    options.update(kwargs)
    return kde_slice(np.asarray(points, dtype=float), **options)


def _plane(xy):
    xy = np.asarray(xy, dtype=float)
    return np.column_stack([xy, np.full(len(xy), 0.5)])


class KdeDeclineTests(unittest.TestCase):
    def assertDeclined(self, result, reason):
        self.assertEqual(result["message"], KDE_MESSAGES[reason])
        self.assertEqual(result["fitCount"], 0)
        self.assertIsNone(result["kernel"])
        self.assertEqual(result["vmin"], 0.0)
        self.assertEqual(result["vmax"], 0.0)

    def test_empty_slab(self):
        result = _slice(_plane([[0.2, 0.3]]), z_center=0.9)
        self.assertEqual(result["slabCount"], 0)
        self.assertDeclined(result, "empty")

    def test_too_few_rows(self):
        self.assertDeclined(_slice(_plane([[0.2, 0.3], [0.5, 0.6], [0.4, 0.8]])), "too_few")

    def test_fewer_than_three_distinct_positions(self):
        self.assertDeclined(_slice(_plane([[0.3, 0.4]] * 10 + [[0.6, 0.7]] * 10)), "few_unique")

    def test_exactly_collinear_points(self):
        t = np.arange(16, 49) / 64.0  # dyadic: exactly on the line
        self.assertDeclined(_slice(_plane(np.column_stack([t, 0.25 + 0.5 * t]))), "collinear")

    def test_collinear_to_round_off_is_declined_not_drawn(self):
        # Coordinates written to 10 decimals (as a configuration file would
        # carry them) pass numpy's rank test (rank 2), but 1 - rho^2 is at
        # round-off level, where the sign of the Cholesky pivot is arbitrary.
        t = np.linspace(0.2, 0.8, 40)
        points = np.round(np.column_stack([t, 0.3 + 0.7 * t]), 10)
        self.assertEqual(np.linalg.matrix_rank(points - points.mean(axis=0)), 2)
        cov = np.cov(points, rowvar=False)
        self.assertLess(1 - cov[0, 1] ** 2 / (cov[0, 0] * cov[1, 1]), COVARIANCE_CONDITION_LIMIT)
        self.assertDeclined(_slice(_plane(points)), "singular")

    def test_non_positive_or_non_finite_bandwidth_is_declined(self):
        rng = np.random.default_rng(3)
        points = _plane(rng.random((50, 2)))
        for bw in (0.0, -0.03, float("nan"), float("inf")):
            with self.subTest(bw=bw):
                result = _slice(points, bw=bw)
                self.assertDeclined(result, "bandwidth")
                self.assertIsNone(result["bw"])

    def test_near_collinear_but_two_dimensional_slab_is_drawn_exactly(self):
        rng = np.random.default_rng(4)
        t = np.linspace(0.2, 0.8, 40)
        points = np.column_stack([t, 0.3 + 0.7 * t + 1e-4 * rng.standard_normal(40)])
        result = _slice(_plane(points), bw=0.03, grid=24)
        self.assertIsNone(result["message"])
        self.assertEqual(result["fitCount"], 40)
        # No images, no subsample: kde_slice's kernel is gaussian_kde's own. The
        # covariances differ only by round-off (np.cov vs scipy's weighted
        # cov), which this needle's conditioning (~1e6) amplifies to ~1e-10.
        reference = gaussian_kde(points.T, bw_method=0.03)
        np.testing.assert_allclose(result["kernel"]["covariance"], reference.covariance, rtol=1e-12)
        axis = np.linspace(0.0, 1.0, 24)
        mesh_x, mesh_y = np.meshgrid(axis, axis)
        expected = reference(np.vstack([mesh_x.ravel(), mesh_y.ravel()])).reshape(mesh_x.shape)
        np.testing.assert_allclose(result["density"], expected, rtol=1e-9, atol=1e-9 * expected.max())

    def test_kernel_summary_reports_principal_sigmas(self):
        rng = np.random.default_rng(5)
        points = np.column_stack([0.5 + 0.1 * rng.standard_normal(400), 0.5 + 0.02 * rng.standard_normal(400)])
        result = _slice(_plane(points), bw=0.05)
        eigenvalues = np.linalg.eigvalsh(np.asarray(result["kernel"]["covariance"]))
        self.assertAlmostEqual(result["kernel"]["sigmaMinor"], float(np.sqrt(eigenvalues[0])), places=12)
        self.assertAlmostEqual(result["kernel"]["sigmaMajor"], float(np.sqrt(eigenvalues[1])), places=12)


if __name__ == "__main__":
    unittest.main()
