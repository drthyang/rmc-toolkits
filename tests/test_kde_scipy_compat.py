# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""``_FixedCovarianceKDE`` evaluates the supplied kernel on every SciPy it can meet.

The subclass replaces the covariance that ``scipy.stats.gaussian_kde`` would
estimate. Which attribute SciPy's evaluator then reads depends on the release:
SciPy >= 1.10 whitens with ``cho_cov``, earlier compiled evaluators (1.8.1
checked) with ``inv_cov`` (the precision matrix, through its own Cholesky
factor), and the oldest pure-Python ``evaluate`` divides by ``_norm_factor``.
The pyproject pins no SciPy floor, so the kernel must be the supplied one
through each of those readers. The pre-1.10 evaluators are replicated here from
their SciPy sources (the suite runs on one SciPy); this file and the other KDE
tests also pass on a real SciPy 1.8.1.

When SciPy evaluates something else anyway, ``kde_slice`` declines with the
``engine`` message instead of raising, so ``/api/kde/slice`` answers 200.
"""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
from scipy.linalg import solve_triangular

from rmc_toolkits import kde as kde_module
from rmc_toolkits.kde import KDE_MESSAGES, _FixedCovarianceKDE, kde_slice

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"


def _slab_rows(seed: int = 3, n: int = 400) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n, 2)) @ np.array([[0.08, 0.0], [0.05, 0.02]]) + 0.5


def _source_covariance() -> np.ndarray:
    # Deliberately not the dataset's covariance: that is the point of the class.
    return np.array([[0.0061, 0.0019], [0.0019, 0.0013]])


def _scipy_pre_1_10_evaluate(kde, points: np.ndarray) -> np.ndarray:
    """scipy/stats/_kde.py evaluate() + _stats.pyx gaussian_kernel_estimate before SciPy 1.10.

    The data and the points are whitened by the Cholesky factor of the precision
    matrix ``inv_cov``; the normalisation is ``(2 pi)^(-d/2) prod(diag(W))``.
    """
    whitening = np.linalg.cholesky(kde.inv_cov)
    data = kde.dataset.T @ whitening
    xi = points.T @ whitening
    norm = (2 * np.pi) ** (-kde.d / 2) * np.prod(np.diag(whitening))
    values = np.empty(xi.shape[0])
    for index, point in enumerate(xi):
        values[index] = np.sum(np.exp(-0.5 * np.sum((data - point) ** 2, axis=1)) * kde.weights)
    return values * norm


def _scipy_pure_python_evaluate(kde, points: np.ndarray) -> np.ndarray:
    """scipy/stats/kde.py evaluate() from before the compiled gaussian_kernel_estimate."""
    values = np.empty(points.shape[1])
    for index in range(points.shape[1]):
        diff = kde.dataset - points[:, index, np.newaxis]
        energy = np.sum(diff * np.dot(kde.inv_cov, diff), axis=0) / 2.0
        values[index] = np.sum(np.exp(-energy) * kde.weights, axis=0)
    return values / kde._norm_factor


class FixedCovarianceKdeCompatTests(unittest.TestCase):
    def setUp(self):
        self.rows = _slab_rows()
        self.covariance = _source_covariance()
        self.bw = 0.5
        self.kde = _FixedCovarianceKDE(self.rows.T, self.covariance, self.bw)
        grid = np.linspace(0.35, 0.65, 9)
        mesh_x, mesh_y = np.meshgrid(grid, grid)
        self.points = np.vstack([mesh_x.ravel(), mesh_y.ravel()])

    def direct_sum(self, points: np.ndarray) -> np.ndarray:
        kernel = self.bw**2 * self.covariance
        inverse = np.linalg.inv(kernel)
        diff = self.rows[:, None, :] - points.T[None, :, :]
        energy = np.einsum("npi,ij,npj->np", diff, inverse, diff)
        return np.exp(-0.5 * energy).sum(axis=0) / (
            self.rows.shape[0] * 2 * np.pi * np.sqrt(np.linalg.det(kernel))
        )

    def assert_is_the_direct_sum(self, values):
        expected = self.direct_sum(self.points)
        np.testing.assert_allclose(values, expected, rtol=1e-12, atol=1e-13 * expected.max())

    def test_scipy_evaluates_the_supplied_kernel(self):
        self.assert_is_the_direct_sum(self.kde(self.points))

    def test_inv_cov_is_the_supplied_kernel_precision(self):
        # SciPy >= 1.10 exposes inv_cov as a property that re-estimates the
        # covariance from the dataset (and overwrites _data_covariance).
        expected = np.linalg.inv(self.bw**2 * self.covariance)
        np.testing.assert_allclose(self.kde.inv_cov, expected, rtol=1e-12)
        np.testing.assert_array_equal(self.kde._data_covariance, self.covariance)
        np.testing.assert_allclose(self.kde.covariance, self.bw**2 * self.covariance, rtol=1e-15)

    def test_pre_1_10_evaluator_sees_the_supplied_kernel(self):
        self.assert_is_the_direct_sum(_scipy_pre_1_10_evaluate(self.kde, self.points))

    def test_pure_python_evaluator_sees_the_supplied_kernel(self):
        self.assert_is_the_direct_sum(_scipy_pure_python_evaluate(self.kde, self.points))

    def test_the_construction_check_tolerates_a_needle_on_any_evaluator(self):
        # A kernel near COVARIANCE_CONDITION_LIMIT (cond(H) ~ 2e10). Correct
        # evaluators that whiten differently disagree at O(cond * eps) -- up to
        # 0.33 cond * eps measured for the pre-1.10 path against chol(H), i.e.
        # 1.7e-6 here -- which the construction check must not mistake for
        # ignored internals (a flat 1e-6 tolerance did). A kernel that is off
        # by 1e-3 must still be caught.
        rng = np.random.default_rng(38)
        t = rng.uniform(-0.3, 0.3, size=300)
        angle = rng.uniform(0, np.pi)
        along = np.array([np.cos(angle), np.sin(angle)])
        across = np.array([-np.sin(angle), np.cos(angle)])
        rows = 0.5 + t[:, None] * along + (1.2e-6 * rng.normal(size=t.size))[:, None] * across
        covariance = np.cov(rows, rowvar=False)
        rho2 = covariance[0, 1] ** 2 / (covariance[0, 0] * covariance[1, 1])
        self.assertGreater(1 - rho2, kde_module.COVARIANCE_CONDITION_LIMIT)
        self.assertLess(1 - rho2, 1e-9)

        reference = _FixedCovarianceKDE(rows.T, covariance, 0.03)
        cond_eps = reference._condition_number() * np.finfo(float).eps
        self.assertGreater(0.5 * cond_eps, 1e-6)
        original = _FixedCovarianceKDE.evaluate

        def skewed(factor):
            return lambda kde, points: original(kde, points) * factor

        with mock.patch.object(_FixedCovarianceKDE, "evaluate", skewed(1 + 0.5 * cond_eps)):
            _FixedCovarianceKDE(rows.T, covariance, 0.03)
        with mock.patch.object(
            _FixedCovarianceKDE,
            "evaluate",
            lambda kde, points: _scipy_pre_1_10_evaluate(kde, np.atleast_2d(points)),
        ):
            _FixedCovarianceKDE(rows.T, covariance, 0.03)
        with mock.patch.object(_FixedCovarianceKDE, "evaluate", skewed(1 + 1e-3)):
            with self.assertRaises(kde_module.ScipyKdeUnsupported):
                _FixedCovarianceKDE(rows.T, covariance, 0.03)


class KdeSliceDeclinesUnsupportedScipyTests(unittest.TestCase):
    """A SciPy whose evaluator ignores or cannot read the kernel declines the slab."""

    def slab(self):
        rng = np.random.default_rng(5)
        xy = rng.uniform(0.1, 0.9, size=(40, 2))
        return np.column_stack([xy, np.full(len(xy), 0.5)])

    def assert_engine_decline(self, result):
        self.assertEqual(result["message"], KDE_MESSAGES["engine"])
        self.assertIsNone(result["kernel"])
        self.assertEqual(result["fitCount"], 0)
        self.assertEqual(result["slabCount"], 40)
        self.assertEqual(result["contours"], [])
        self.assertEqual(result["vmax"], 0.0)

    def run_slice(self):
        return kde_slice(self.slab(), 0.5, 0.1, xlim=(0.0, 1.0), ylim=(0.0, 1.0), grid=16)

    def test_missing_attribute_declines(self):
        def missing(_kde, _points):
            raise AttributeError("'_FixedCovarianceKDE' object has no attribute 'inv_cov'")

        with mock.patch.object(_FixedCovarianceKDE, "evaluate", missing), self.assertLogs(
            "rmc_toolkits.kde", level="WARNING"
        ) as logs:
            self.assert_engine_decline(self.run_slice())
        self.assertIn("inv_cov", logs.output[0])

    def test_wrong_kernel_declines(self):
        original = _FixedCovarianceKDE.evaluate

        def ignores_the_supplied_kernel(kde, points):
            return 1.5 * original(kde, points)

        with mock.patch.object(_FixedCovarianceKDE, "evaluate", ignores_the_supplied_kernel), self.assertLogs(
            "rmc_toolkits.kde", level="WARNING"
        ):
            self.assert_engine_decline(self.run_slice())


@unittest.skipUnless((DEMO / "GTS_250K.rmc6f").exists(), "committed demo run missing")
class KdeSliceEndpointUnsupportedScipyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend = str(ROOT / "web_app" / "backend")
        if backend not in sys.path:
            sys.path.insert(0, backend)
        os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
        os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
        import app as backend_app

        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

    def test_endpoint_answers_with_a_message_not_a_500(self):
        def missing(_kde, _points):
            raise AttributeError("'_FixedCovarianceKDE' object has no attribute 'inv_cov'")

        with mock.patch.object(_FixedCovarianceKDE, "evaluate", missing), self.assertLogs(
            "rmc_toolkits.kde", level="WARNING"
        ):
            response = self.client.get(
                "/api/kde/slice?dir=web_app/frontend/public/demo&element=Ga&z=0.25&grid=32"
            )
        self.assertEqual(response.status_code, 200, response.get_json())
        payload = response.get_json()
        self.assertEqual(payload["message"], KDE_MESSAGES["engine"])
        self.assertIsNone(payload["kernel"])
        self.assertGreater(payload["slabCount"], 0)

    def test_endpoint_draws_on_this_scipy(self):
        response = self.client.get("/api/kde/slice?dir=web_app/frontend/public/demo&element=Ga&z=0.25&grid=32")
        self.assertEqual(response.status_code, 200, response.get_json())
        payload = response.get_json()
        self.assertIsNone(payload["message"])
        self.assertIsNotNone(payload["kernel"])


if __name__ == "__main__":
    unittest.main()
