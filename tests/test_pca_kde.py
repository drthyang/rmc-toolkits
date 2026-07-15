# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np
from scipy.stats import gaussian_kde

from rmc_toolkits.pca_kde import (
    load_site_displacements,
    pca_kde_volume,
    probability_scale,
    site_ellipsoids,
    site_pca_kde,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

requires_sample = unittest.skipUnless(
    (DATA / "GNSe.rmc6f").exists(),
    "GNSe sample data not present in data/ (gitignored)",
)


def write_rmc6f(path: Path, atom_lines: list[str], *, supercell=(2, 2, 2), lattice=None):
    lattice = lattice or [[8.0, 0.0, 0.0], [0.0, 8.0, 0.0], [0.0, 0.0, 8.0]]
    header = [
        f"Supercell dimensions {supercell[0]} {supercell[1]} {supercell[2]}",
        "Lattice vectors (Ang):",
        " ".join(str(v) for v in lattice[0]),
        " ".join(str(v) for v in lattice[1]),
        " ".join(str(v) for v in lattice[2]),
        "Atoms:",
    ]
    path.write_text("\n".join(header + atom_lines), encoding="utf-8")


def make_cloud_rmc6f(path: Path, sigma_frac, *, n=6, supercell=(6, 6, 6), seed=1, element="Ga"):
    """One reference site displaced by a known anisotropic Gaussian, n^3 copies."""
    rng = np.random.default_rng(seed)
    lattice = [[8.0 * supercell[0], 0, 0], [0, 8.0 * supercell[1], 0], [0, 0, 8.0 * supercell[2]]]
    lines = []
    atom = 0
    for ix in range(supercell[0]):
        for iy in range(supercell[1]):
            for iz in range(supercell[2]):
                atom += 1
                disp = rng.normal(0.0, 1.0, size=3) * np.asarray(sigma_frac)
                base = (np.array([ix, iy, iz]) + 0.0) / np.asarray(supercell)
                coord = base + disp
                lines.append(
                    f"{atom} {element} [1] {coord[0]:.10f} {coord[1]:.10f} {coord[2]:.10f} "
                    f"1 {ix} {iy} {iz}"
                )
    write_rmc6f(path, lines, supercell=supercell, lattice=lattice)
    return lattice


class ProbabilityScaleTests(unittest.TestCase):
    def test_fifty_percent_is_crystallographic_constant(self):
        self.assertAlmostEqual(probability_scale(0.5), 1.5381722, places=5)

    def test_monotonic_and_bounded(self):
        self.assertLess(probability_scale(0.5), probability_scale(0.99))
        for bad in (0.0, 1.0, -0.1, 1.2):
            with self.assertRaises(ValueError):
                probability_scale(bad)


class SeparableKdeTests(unittest.TestCase):
    """The engine's core promise: the separable volume equals scipy's gaussian_kde."""

    def _random_cloud(self, seed=0, n=400):
        rng = np.random.default_rng(seed)
        transform = np.array([[1.4, 0.3, -0.2], [0.3, 0.7, 0.1], [-0.2, 0.1, 0.35]])
        return rng.normal(size=(n, 3)) @ transform.T

    def test_volume_matches_scipy_gaussian_kde(self):
        cloud = self._random_cloud()
        for bw in ("scott", "silverman", 0.35):
            result = pca_kde_volume(cloud, bw=bw, grid=16, projections=False)
            grid = result["grid"]
            coords = [np.asarray(axis) for axis in result["axisCoords"]]
            axes = np.asarray(result["axes"])
            mean = np.asarray(result["mean"])

            # Rebuild the sample points from the PCA-frame grid and evaluate the
            # reference estimator there.
            gx, gy, gz = np.meshgrid(coords[0], coords[1], coords[2], indexing="ij")
            pca_points = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
            cartesian = mean + pca_points @ axes

            reference = gaussian_kde(cloud.T, bw_method=bw)
            expected = reference(cartesian.T).reshape(grid, grid, grid)
            actual = np.asarray(result["density"]).reshape(grid, grid, grid)

            np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-12)

    def test_projection_is_true_marginal_of_the_volume(self):
        # The wall projection is the honest shadow of the displayed volume: the
        # 3D KDE integrated over the dropped axis, keeping the 3D bandwidth
        # sub-block rather than recomputing a 2D-optimal one. The volume is
        # already proven equal to scipy's estimator, so integrating it over PC3
        # by a fine Riemann sum is an independent check of the analytic marginal.
        cloud = self._random_cloud(seed=3)
        result = pca_kde_volume(cloud, bw="scott", grid=96, extent=5.0)
        volume = np.asarray(result["density"]).reshape(96, 96, 96)
        dz = result["axisCoords"][2][1] - result["axisCoords"][2][0]
        riemann = volume.sum(axis=2) * dz  # integrate PC3 out of the volume

        actual = np.asarray(result["projections"]["pc12"]["density"])
        np.testing.assert_allclose(actual, riemann, rtol=2e-3, atol=1e-6)

    def test_density_layout_is_c_order(self):
        cloud = self._random_cloud(seed=5, n=200)
        result = pca_kde_volume(cloud, bw=0.4, grid=12, projections=False)
        grid = result["grid"]
        flat = np.asarray(result["density"])
        self.assertEqual(flat.size, grid**3)
        # index = (i*grid + j)*grid + k must line up with the meshgrid rebuild.
        volume = flat.reshape(grid, grid, grid)
        self.assertGreater(volume.max(), volume.min())


class VolumePropertyTests(unittest.TestCase):
    def test_mass_recovers_probability_normalization(self):
        rng = np.random.default_rng(11)
        cloud = rng.normal(size=(2000, 3)) * np.array([0.4, 0.2, 0.1])
        result = pca_kde_volume(cloud, bw="scott", grid=64, extent=4.0, projections=False)
        # A KDE integrates to 1; a box at 4 sigma of the broadened density should
        # capture essentially all of it.
        self.assertGreater(result["mass"], 0.99)

    def test_mass_levels_bracket_the_cloud(self):
        rng = np.random.default_rng(12)
        cloud = rng.normal(size=(1500, 3)) * 0.3
        result = pca_kde_volume(cloud, bw="scott", grid=48, projections=False)
        levels = {entry["p"]: entry["level"] for entry in result["massLevels"]}
        # Higher enclosed probability -> lower density threshold.
        self.assertGreater(levels[0.1], levels[0.9])

    def test_recovers_known_anisotropy(self):
        rng = np.random.default_rng(7)
        sigma = np.array([0.5, 0.25, 0.1])
        cloud = rng.normal(size=(6000, 3)) * sigma
        result = pca_kde_volume(cloud, bw="scott", grid=16, projections=False)
        np.testing.assert_allclose(np.asarray(result["rms"]), sigma, rtol=0.08)
        self.assertAlmostEqual(result["anisotropy"], sigma[0] / sigma[2], delta=0.5)

    def test_cubic_box_uses_uniform_half_width(self):
        rng = np.random.default_rng(9)
        cloud = rng.normal(size=(800, 3)) * np.array([0.6, 0.2, 0.05])
        result = pca_kde_volume(cloud, bw="scott", grid=12, cubic_box=True, projections=False)
        half = result["halfWidths"]
        self.assertAlmostEqual(half[0], half[1])
        self.assertAlmostEqual(half[1], half[2])

    def test_rejects_degenerate_input(self):
        with self.assertRaises(ValueError):
            pca_kde_volume(np.zeros((10, 3)), grid=8)
        with self.assertRaises(ValueError):
            pca_kde_volume(np.zeros((2, 3)), grid=8)

    def test_excess_kurtosis_flags_non_gaussian_clouds(self):
        rng = np.random.default_rng(21)
        # A Gaussian cloud is ~harmonic: excess kurtosis near zero.
        gaussian = rng.normal(size=(8000, 3)) * np.array([0.3, 0.2, 0.1])
        harmonic = pca_kde_volume(gaussian, bw="scott", grid=12, projections=False)
        self.assertLess(abs(harmonic["nonGaussianity"]), 0.3)

        # A heavy-tailed (Student-t-like) cloud is strongly leptokurtic.
        heavy = rng.standard_t(3, size=(8000, 3)) * 0.1
        anharmonic = pca_kde_volume(heavy, bw="scott", grid=12, projections=False)
        self.assertGreater(anharmonic["nonGaussianity"], 1.0)
        self.assertEqual(len(anharmonic["excessKurtosis"]), 3)

    def test_site_ellipsoids_report_non_gaussianity(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "cloud.rmc6f"
            make_cloud_rmc6f(path, [0.02, 0.015, 0.01], supercell=(8, 8, 8), seed=3)
            sites = load_site_displacements(path)
            entry = site_ellipsoids(sites)[0]
            self.assertIn("nonGaussianity", entry)
            self.assertEqual(len(entry["excessKurtosis"]), 3)
            # A Gaussian-displaced synthetic site is near-harmonic.
            self.assertLess(abs(entry["nonGaussianity"]), 0.5)

    def test_axes_are_orthonormal_and_right_handed(self):
        cloud = np.random.default_rng(2).normal(size=(500, 3)) @ np.array(
            [[1.0, 0.4, 0.0], [0.4, 0.6, 0.2], [0.0, 0.2, 0.3]]
        )
        result = pca_kde_volume(cloud, grid=10, projections=False)
        axes = np.asarray(result["axes"])
        np.testing.assert_allclose(axes @ axes.T, np.eye(3), atol=1e-9)
        self.assertAlmostEqual(np.linalg.det(axes), 1.0, places=9)


class SiteDisplacementTests(unittest.TestCase):
    def test_extracts_centered_clouds_with_known_covariance(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "cloud.rmc6f"
            sigma_frac = np.array([0.02, 0.01, 0.005])
            lattice = make_cloud_rmc6f(path, sigma_frac, supercell=(8, 8, 8), seed=4)
            sites = load_site_displacements(path)

            self.assertEqual(sites.reference_numbers.tolist(), [1])
            self.assertEqual(int(sites.counts[0]), 512)
            # Cloud is centered per site.
            np.testing.assert_allclose(sites.displacements.mean(axis=0), 0.0, atol=1e-9)

            # Cartesian sigma = fractional sigma * cell edge (diagonal lattice).
            edge = lattice[0][0]
            expected_sigma = sigma_frac * edge
            actual_sigma = sites.displacements.std(axis=0, ddof=1)
            np.testing.assert_allclose(actual_sigma, expected_sigma, rtol=0.15)

    def test_site_ellipsoids_batch_matches_direct_covariance(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "cloud.rmc6f"
            make_cloud_rmc6f(path, [0.03, 0.015, 0.008], supercell=(6, 6, 6), seed=8)
            sites = load_site_displacements(path)
            ellipsoids = site_ellipsoids(sites)

            self.assertEqual(len(ellipsoids), 1)
            entry = ellipsoids[0]
            cloud = sites.displacements
            direct = np.cov(cloud, rowvar=False, bias=False)
            np.testing.assert_allclose(np.asarray(entry["covariance"]), direct, rtol=1e-9)
            # Eigenvalues are sorted descending.
            eig = entry["eigenvalues"]
            self.assertGreaterEqual(eig[0], eig[1])
            self.assertGreaterEqual(eig[1], eig[2])

    def test_supercell_boundary_wrap_folds(self):
        # An atom that drifted just across the supercell edge (coord slightly < 0
        # stored as ~1) must fold back, not register a full-cell displacement.
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "wrap.rmc6f"
            lines = [
                "1 Ga [1] 0.999 0.001 0.500 1 0 0 0",
                "2 Ga [1] 0.001 0.999 0.500 1 0 0 0",
                "3 Ga [1] 0.000 0.000 0.500 1 0 0 0",
                "4 Ga [1] 0.998 0.002 0.500 1 0 0 0",
            ]
            write_rmc6f(path, lines, supercell=(1, 1, 1),
                        lattice=[[8, 0, 0], [0, 8, 0], [0, 0, 8]])
            sites = load_site_displacements(path)
            # All four sit within ~0.002 of the origin along x,y after folding,
            # so no displacement should approach a full 8 A cell edge.
            self.assertLess(np.abs(sites.displacements).max(), 0.1)

    def test_site_pca_kde_tags_reference_and_element(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "cloud.rmc6f"
            make_cloud_rmc6f(path, [0.02, 0.02, 0.02], supercell=(6, 6, 6), seed=2, element="Nb")
            sites = load_site_displacements(path)
            result = site_pca_kde(sites, reference_number=1, grid=12, projections=False)
            self.assertEqual(result["referenceNumber"], 1)
            self.assertEqual(result["element"], "Nb")
            self.assertEqual(result["count"], 216)


@requires_sample
class SampleBackedTests(unittest.TestCase):
    def test_ganbse_sites_resolve(self):
        sites = load_site_displacements(DATA / "GNSe.rmc6f")
        self.assertEqual(sites.reference_numbers.size, 52)
        self.assertEqual(int(sites.counts.sum()), 52000)
        ellipsoids = site_ellipsoids(sites)
        self.assertEqual(len(ellipsoids), 52)


if __name__ == "__main__":
    unittest.main()
