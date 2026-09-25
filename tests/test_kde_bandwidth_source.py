# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The KDE kernel comes from the slab's source atoms, not from the periodic images.

``H = bw**2 * C`` with ``C`` the covariance of one row per source atom; the
periodic images (whose number depends on the margin
``min(0.5, max(0.1, 2 bw, dz))``) and the 6000-point subsample enter only the
density sum. So moving the thickness or bandwidth slider across a margin step,
or changing the subsample seed, must not change the kernel of a slab that holds
the same atoms.
"""

from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np

from rmc_toolkits.kde import MAX_KDE_FIT_POINTS, load_unit_cell_positions, oriented_kde_slice

ROOT = Path(__file__).resolve().parents[1]
GANB4SE8_5K = ROOT / "data" / "5K_try1" / "GaNb4Se8_5K.rmc6f"

C_SLICE = {
    "normal": np.array([0.0, 0.0, 1.0]),
    "u_axis": np.array([1.0, 0.0, 0.0]),
    "v_axis": np.array([0.0, 1.0, 0.0]),
}


def _two_site_layer(seed: int = 0, per_site: int = 400, sigma: float = 0.006) -> np.ndarray:
    """Two isotropic sites on the cell anti-diagonal at z = 0.25 (like Ga in GaNb4Se8)."""
    rng = np.random.default_rng(seed)
    sites = [np.array([0.25, 0.75, 0.25]), np.array([0.75, 0.25, 0.25])]
    points = np.vstack([site + sigma * rng.standard_normal((per_site, 3)) for site in sites])
    return points % 1.0


def _slice(points, **kwargs):
    options = {"center": 0.25, "thickness": 0.08, "bw": 0.03, "grid": 48, "n_levels": 0, **C_SLICE}
    options.update(kwargs)
    return oriented_kde_slice(points, **options)


class KdeBandwidthSourceTests(unittest.TestCase):
    def test_kernel_is_bw_squared_times_the_source_atom_covariance(self):
        points = _two_site_layer()
        result = _slice(points)
        self.assertEqual(result["slabCount"], points.shape[0])
        expected = 0.03**2 * np.cov(points[:, :2], rowvar=False)
        np.testing.assert_allclose(result["kernel"]["covariance"], expected, rtol=1e-12)

    def test_thickness_across_the_margin_step_leaves_the_kernel_alone(self):
        # dz = 0.08 -> margin 0.1; dz = 0.30 -> margin 0.3 admits the sites'
        # x/y images at 1.25 and -0.25. Same 800 atoms in both slabs.
        points = _two_site_layer()
        thin = _slice(points, thickness=0.08)
        thick = _slice(points, thickness=0.30)
        self.assertEqual(thin["slabCount"], thick["slabCount"])
        self.assertGreater(thick["fitCount"], thin["fitCount"])  # more image rows summed
        np.testing.assert_allclose(thick["kernel"]["covariance"], thin["kernel"]["covariance"], rtol=1e-12)
        # ... and so the site is drawn the same: equal density at the site centre.
        np.testing.assert_allclose(thick["vmax"], thin["vmax"], rtol=1e-6)

    def test_bandwidth_across_the_margin_step_scales_the_kernel_as_bw_squared(self):
        points = _two_site_layer()
        below = _slice(points, bw=0.12)  # margin 0.24
        above = _slice(points, bw=0.13)  # margin 0.26: the x/y images at 1.25 / -0.25 enter
        ratio = np.asarray(above["kernel"]["covariance"]) / np.asarray(below["kernel"]["covariance"])
        np.testing.assert_allclose(ratio, (0.13 / 0.12) ** 2, rtol=1e-12)

    def test_depth_wrapped_atoms_are_represented_by_their_nearest_image(self):
        rng = np.random.default_rng(2)
        points = np.column_stack([rng.random(60), rng.random(60), 0.97 + 0.02 * rng.random(60)])
        result = _slice(points, center=0.0, thickness=0.1)
        self.assertEqual(result["slabCount"], 60)
        expected = 0.03**2 * np.cov(points[:, :2], rowvar=False)
        np.testing.assert_allclose(result["kernel"]["covariance"], expected, rtol=1e-12)

    def test_kernel_does_not_depend_on_the_fit_subsample(self):
        rng = np.random.default_rng(3)
        count = MAX_KDE_FIT_POINTS + 2000
        points = np.column_stack([rng.random(count), rng.random(count), np.full(count, 0.5)])
        first = _slice(points, center=0.5, rng_seed=0)
        second = _slice(points, center=0.5, rng_seed=1)
        self.assertEqual(first["fitCount"], MAX_KDE_FIT_POINTS)
        self.assertEqual(first["kernel"], second["kernel"])
        expected = 0.03**2 * np.cov(points[:, :2], rowvar=False)
        np.testing.assert_allclose(first["kernel"]["covariance"], expected, rtol=1e-12)

    def test_fixed_covariance_kde_is_scipys_sum_and_guards_against_scipy_changes(self):
        from unittest import mock

        from rmc_toolkits.kde import _FixedCovarianceKDE

        rng = np.random.default_rng(4)
        data = rng.random((2, 300))
        covariance = np.array([[0.02, 0.005], [0.005, 0.01]])
        kde = _FixedCovarianceKDE(data, covariance, 0.05)
        np.testing.assert_allclose(kde.covariance, 0.05**2 * covariance, rtol=1e-15)
        # The density is the Gaussian mixture with that kernel, not with cov(data).
        point = np.array([[0.4], [0.6]])
        inverse = np.linalg.inv(kde.covariance)
        offsets = data - point
        quadratic = np.einsum("in,ij,jn->n", offsets, inverse, offsets)
        expected = np.exp(-0.5 * quadratic).sum() / (300 * 2 * np.pi * np.sqrt(np.linalg.det(kde.covariance)))
        self.assertAlmostEqual(float(kde(point)[0]) / expected, 1.0, places=10)
        # A scipy whose evaluate() ignored the supplied covariance must raise.
        with mock.patch.object(_FixedCovarianceKDE, "evaluate", lambda self, points: np.array([1.0])):
            with self.assertRaises(RuntimeError):
                _FixedCovarianceKDE(data, covariance, 0.05)

    @unittest.skipUnless(GANB4SE8_5K.exists(), "GaNb4Se8 5K run not present in data/ (gitignored)")
    def test_real_ga_layer_kernel_is_independent_of_the_slab_thickness(self):
        # The z = 0.75 Ga layer of GaNb4Se8; the z = 0.25 layer stays outside
        # every slab, so all four thicknesses select the same 2000 atoms.
        ga = load_unit_cell_positions(GANB4SE8_5K, element="Ga").fractional_positions
        results = [_slice(ga, center=0.75, thickness=dz, grid=24) for dz in (0.08, 0.2, 0.3, 0.5)]
        self.assertEqual({result["slabCount"] for result in results}, {2000})
        for result in results[1:]:
            np.testing.assert_allclose(
                result["kernel"]["covariance"], results[0]["kernel"]["covariance"], rtol=1e-12
            )


if __name__ == "__main__":
    unittest.main()
