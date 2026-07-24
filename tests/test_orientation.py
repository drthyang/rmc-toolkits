# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from rmc_toolkits.orientation import (
    MAX_FREQUENCY,
    MIN_FREQUENCY,
    assign_cells,
    goldberg_tiling,
    orientation_histogram,
    recommended_frequency,
    site_orientation_histogram,
)
from rmc_toolkits.pca_kde import load_site_displacements


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


def make_cloud_rmc6f(path: Path, sigma_frac, *, supercell=(6, 6, 6), seed=1, element="Ga"):
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


class TilingTests(unittest.TestCase):
    """The Goldberg tessellation invariants the histogram relies on."""

    def test_cell_count_and_pentagon_count(self):
        for frequency in (1, 2, 3, 5, 8):
            tiling = goldberg_tiling(frequency)
            self.assertEqual(tiling.cell_count, 10 * frequency**2 + 2)
            self.assertEqual(int(tiling.is_pentagon.sum()), 12)
            self.assertTrue(np.isin(tiling.sizes, (5, 6)).all())

    def test_centers_are_unit_vectors(self):
        tiling = goldberg_tiling(6)
        np.testing.assert_allclose(np.linalg.norm(tiling.centers, axis=1), 1.0, atol=1e-12)

    def test_areas_are_positive_and_sum_to_full_sphere(self):
        for frequency in (1, 4, 9):
            tiling = goldberg_tiling(frequency)
            self.assertGreater(tiling.areas.min(), 0.0)
            self.assertAlmostEqual(tiling.areas.sum(), 4.0 * np.pi, places=9)

    def test_area_spread_is_modest(self):
        # The whole point of the Goldberg grid: cells are nearly equal-area,
        # unlike a (theta, phi) grid whose polar cells collapse to slivers.
        tiling = goldberg_tiling(8)
        self.assertLess(tiling.areas.max() / tiling.areas.min(), 2.0)

    def test_neighbors_are_symmetric(self):
        tiling = goldberg_tiling(4)
        for cell in range(tiling.cell_count):
            for neighbor in tiling.neighbors[cell]:
                if neighbor < 0:
                    continue
                self.assertIn(cell, tiling.neighbors[neighbor])

    def test_antipode_is_exact(self):
        tiling = goldberg_tiling(5)
        np.testing.assert_allclose(
            tiling.centers[tiling.antipode], -tiling.centers, atol=1e-12
        )
        # Folding is an involution.
        np.testing.assert_array_equal(tiling.antipode[tiling.antipode], np.arange(tiling.cell_count))

    def test_polygons_enclose_their_center(self):
        # Each cell's polygon vertices should all be closer to their own center
        # than to any other center is far -- specifically the center should have
        # positive dot with every vertex (cells are small caps).
        tiling = goldberg_tiling(6)
        for cell in (0, 13, 101, tiling.cell_count - 1):
            size = tiling.sizes[cell]
            dots = tiling.polygons[cell, :size] @ tiling.centers[cell]
            self.assertGreater(dots.min(), 0.9)

    def test_rejects_out_of_range_frequency(self):
        with self.assertRaises(ValueError):
            goldberg_tiling(MIN_FREQUENCY - 1)
        with self.assertRaises(ValueError):
            goldberg_tiling(MAX_FREQUENCY + 1)


class AssignmentTests(unittest.TestCase):
    def test_matches_brute_force_nearest_center(self):
        rng = np.random.default_rng(3)
        directions = rng.normal(size=(4000, 3))
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        for frequency in (2, 6, 11):
            tiling = goldberg_tiling(frequency)
            assigned = assign_cells(tiling, directions)
            brute = np.argmax(directions @ tiling.centers.T, axis=1)
            # Boundary ties are the only legitimate disagreements; random
            # directions almost surely avoid them entirely.
            agreement = float((assigned == brute).mean())
            self.assertGreater(agreement, 0.9999)

    def test_cell_centers_assign_to_themselves(self):
        tiling = goldberg_tiling(7)
        assigned = assign_cells(tiling, tiling.centers)
        np.testing.assert_array_equal(assigned, np.arange(tiling.cell_count))

    def test_accepts_unnormalized_directions(self):
        tiling = goldberg_tiling(3)
        direction = np.array([[0.0, 0.0, 5.0]])
        self.assertEqual(
            assign_cells(tiling, direction)[0],
            assign_cells(tiling, direction / 5.0)[0],
        )


class RecommendedFrequencyTests(unittest.TestCase):
    def test_scales_with_point_count(self):
        self.assertLessEqual(
            recommended_frequency(100), recommended_frequency(10000)
        )
        self.assertEqual(recommended_frequency(0), MIN_FREQUENCY)

    def test_respects_target_per_cell(self):
        frequency = recommended_frequency(12000, target_per_cell=12)
        cells = 10 * frequency**2 + 2
        per_cell = 12000 / cells
        self.assertGreater(per_cell, 6)  # never wildly over-binned

    def test_exact_values_pin_cross_engine_parity(self):
        # Hard-coded expectations shared verbatim with the JS suite, so the
        # browser and server always pick the same tiling for the same data.
        # 774 points put the sqrt at exactly 2.5: Python rounds half to even
        # (nu=2), and the JS port must match, not Math.round's half-up 3.
        self.assertEqual(recommended_frequency(774), 2)
        self.assertEqual(recommended_frequency(300), 2)
        self.assertEqual(recommended_frequency(12000), 10)


class HistogramTests(unittest.TestCase):
    def _isotropic(self, n=20000, seed=0):
        rng = np.random.default_rng(seed)
        return rng.normal(size=(n, 3)) * 0.1

    def _lobed(self, n=20000, seed=1, sigma=(0.5, 0.1, 0.1)):
        rng = np.random.default_rng(seed)
        return rng.normal(size=(n, 3)) * np.asarray(sigma)

    def test_polygon_areas_match_assignment_voronoi(self):
        # The normalization divides by _spherical_polygon_areas, and the honesty
        # of the whole map rests on those analytic areas matching the regions
        # assign_cells actually routes directions into. Cross-check them by
        # Monte Carlo: the fraction of uniform random directions landing in a
        # cell, times 4*pi, must reproduce the cell's analytic solid angle.
        # (A plain sum(density * areas) == 1 check would be circular -- density
        # is defined by dividing by these same areas.)
        rng = np.random.default_rng(17)
        directions = rng.normal(size=(200000, 3))
        tiling = goldberg_tiling(3)
        cells = assign_cells(tiling, directions)
        fractions = np.bincount(cells, minlength=tiling.cell_count) / directions.shape[0]
        monte_carlo_areas = fractions * 4.0 * np.pi
        # Per-cell Poisson noise is ~1/sqrt(200000/92) = 2%; 12% is >5 sigma.
        np.testing.assert_allclose(monte_carlo_areas, tiling.areas, rtol=0.12)

    def test_isotropic_cloud_is_flat(self):
        result = orientation_histogram(self._isotropic(), frequency=8)
        enhancement = np.asarray(result["enhancement"])
        # ~31 points per cell -> Poisson scatter ~18%; the mean must be very
        # close to 1 and no cell should stray far.
        self.assertAlmostEqual(float(np.mean(enhancement)), 1.0, delta=0.05)
        self.assertLess(float(np.abs(enhancement - 1.0).max()), 1.0)
        # Overall departure from isotropy consistent with pure Poisson noise.
        self.assertLess(result["significance"], 1.5)
        self.assertLess(abs(result["orientationAnisotropy"]), 0.05)

    def test_lobed_cloud_peaks_on_the_long_axis(self):
        result = orientation_histogram(self._lobed(), frequency=8)
        peak = np.asarray(result["peakDirection"])
        self.assertGreater(abs(peak[0]), 0.9)
        self.assertGreater(result["peakEnhancement"], 2.0)
        self.assertGreater(result["significance"], 3.0)
        self.assertGreater(result["orientationAnisotropy"], 0.5)

    def test_one_sided_cloud_is_antipodally_asymmetric(self):
        # Displacements biased into the +x hemisphere: the map must show it --
        # this odd-anharmonicity signature is the purpose of the view, and the
        # asymmetry readout must flag it well above the Poisson noise floor.
        rng = np.random.default_rng(5)
        cloud = rng.normal(size=(8000, 3)) * 0.05
        cloud[:, 0] = np.abs(cloud[:, 0]) + 0.05
        result = orientation_histogram(cloud, frequency=6)
        tiling = goldberg_tiling(6)
        enhancement = np.asarray(result["enhancement"])
        peak = int(result["peakCell"])
        self.assertGreater(enhancement[peak], 4.0 * enhancement[tiling.antipode[peak]])
        self.assertGreater(result["antipodalAsymmetry"], 0.5)
        self.assertGreater(
            result["antipodalAsymmetry"], 3.0 * result["antipodalAsymmetryNull"]
        )

    def test_symmetric_cloud_asymmetry_sits_at_the_noise_floor(self):
        # A centrosymmetric (Gaussian) cloud has no real +u/-u imbalance: the
        # asymmetry readout should land near its predicted Poisson null, far
        # below the one-sided regime.
        result = orientation_histogram(self._lobed(n=20000, seed=13), frequency=6)
        null = result["antipodalAsymmetryNull"]
        self.assertLess(result["antipodalAsymmetry"], 2.0 * null)
        self.assertGreater(result["antipodalAsymmetry"], 0.2 * null)

    def test_smoothing_conserves_mass(self):
        cloud = self._lobed(n=5000)
        raw = orientation_histogram(cloud, frequency=8, smoothing=0)
        smooth = orientation_histogram(cloud, frequency=8, smoothing=3)
        self.assertAlmostEqual(
            float(np.sum(smooth["mass"])), float(np.sum(raw["mass"])), places=9
        )
        # Smoothing lowers the peak but cannot move the lobe.
        self.assertLess(smooth["peakEnhancement"], raw["peakEnhancement"])
        self.assertGreater(abs(np.asarray(smooth["peakDirection"])[0]), 0.85)

    def test_cell_mean_amplitude_tracks_per_direction_travel(self):
        # Two one-sided populations along +/-x with 4x different amplitudes:
        # the relief quantity must read each direction's own mean |dr|, and an
        # empty direction (+/-z here is populated, but a sparse tiling cell far
        # from x may be empty) must report 0, not NaN.
        rng = np.random.default_rng(23)
        plus = np.abs(rng.normal(size=(4000, 1))) * 0.05 + 0.35
        minus = -(np.abs(rng.normal(size=(4000, 1))) * 0.012 + 0.09)
        lateral = rng.normal(size=(8000, 2)) * 0.004
        cloud = np.column_stack([np.vstack([plus, minus]), lateral])
        result = orientation_histogram(cloud, frequency=4)
        tiling = goldberg_tiling(4)
        mean_amplitude = np.asarray(result["cellMeanAmplitude"])
        plus_cell = assign_cells(tiling, np.array([[1.0, 0.0, 0.0]]))[0]
        minus_cell = assign_cells(tiling, np.array([[-1.0, 0.0, 0.0]]))[0]
        self.assertGreater(mean_amplitude[plus_cell], 3.0 * mean_amplitude[minus_cell])
        self.assertAlmostEqual(mean_amplitude[plus_cell], 0.39, delta=0.05)
        # Empty cells report exactly 0 and stay finite.
        counts = np.asarray(result["counts"])
        self.assertTrue(np.all(mean_amplitude[counts == 0] == 0.0))
        self.assertTrue(np.all(np.isfinite(mean_amplitude)))

    def test_cell_mean_amplitude_smoothing_stays_bounded(self):
        # Smoothing the ratio's numerator and denominator (not the ratio) keeps
        # every smoothed mean inside the raw data's amplitude range.
        rng = np.random.default_rng(29)
        cloud = rng.normal(size=(5000, 3)) * np.array([0.4, 0.1, 0.1])
        result = orientation_histogram(cloud, frequency=5, smoothing=2)
        mean_amplitude = np.asarray(result["cellMeanAmplitude"])
        amplitudes = np.linalg.norm(cloud, axis=1)
        self.assertTrue(np.all(mean_amplitude <= amplitudes.max() + 1e-12))
        self.assertTrue(np.all(mean_amplitude >= 0.0))
        # Diffusion populates the neighbourhood: far fewer zero cells than raw.
        raw = orientation_histogram(cloud, frequency=5, smoothing=0)
        self.assertLess(
            np.sum(np.asarray(result["cellMeanAmplitude"]) == 0.0),
            max(np.sum(np.asarray(raw["cellMeanAmplitude"]) == 0.0), 1),
        )

    def test_amplitude_weighting_changes_the_map(self):
        # Two populations: many small displacements along x, few large along z.
        rng = np.random.default_rng(9)
        small = rng.normal(size=(9000, 3)) * np.array([0.08, 0.01, 0.01])
        large = rng.normal(size=(1000, 3)) * np.array([0.02, 0.02, 0.8])
        cloud = np.vstack([small, large])
        by_count = orientation_histogram(cloud, frequency=6, weight="count")
        by_amp2 = orientation_histogram(cloud, frequency=6, weight="amplitude2")
        # Counting heads, x wins; weighting by |dr|^2, z wins.
        self.assertGreater(abs(np.asarray(by_count["peakDirection"])[0]), 0.9)
        self.assertGreater(abs(np.asarray(by_amp2["peakDirection"])[2]), 0.9)

    def test_amplitude_cutoffs_drop_points(self):
        rng = np.random.default_rng(4)
        cloud = rng.normal(size=(1000, 3)) * 0.1
        result = orientation_histogram(cloud, frequency=3, min_amplitude_quantile=0.25)
        self.assertEqual(result["totalPoints"], 1000)
        self.assertAlmostEqual(result["usedPoints"], 750, delta=2)
        self.assertEqual(result["usedPoints"] + result["rejectedPoints"], 1000)

    def test_pca_frame_aligns_lobe_with_first_axis(self):
        # A cloud stretched along an oblique direction must, in the PCA frame,
        # peak along +/- PC1 = the local x axis.
        rng = np.random.default_rng(6)
        base = rng.normal(size=(15000, 3)) * np.array([0.5, 0.08, 0.08])
        angle = 0.7
        rotation = np.array(
            [
                [np.cos(angle), -np.sin(angle), 0.0],
                [np.sin(angle), np.cos(angle), 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        cloud = base @ rotation.T
        result = orientation_histogram(cloud, frequency=8, frame="pca")
        peak = np.asarray(result["peakDirection"])
        self.assertGreater(abs(peak[0]), 0.9)

    def test_auto_frequency_tracks_point_count(self):
        small = orientation_histogram(self._isotropic(n=300, seed=2))
        large = orientation_histogram(self._isotropic(n=30000, seed=2))
        self.assertLessEqual(small["frequency"], large["frequency"])
        self.assertEqual(small["frequency"], small["recommendedFrequency"])

    def test_orientation_tensor_matches_direct_computation(self):
        cloud = self._lobed(n=3000, seed=8)
        result = orientation_histogram(cloud, frequency=4)
        directions = cloud / np.linalg.norm(cloud, axis=1, keepdims=True)
        direct = directions.T @ directions / directions.shape[0]
        np.testing.assert_allclose(np.asarray(result["orientationTensor"]), direct, atol=1e-9)
        # Eigenvalues of <u u^T> sum to 1.
        self.assertAlmostEqual(float(np.sum(result["orientationEigenvalues"])), 1.0, places=9)

    def test_geometry_toggle(self):
        cloud = self._isotropic(n=500, seed=1)
        with_geometry = orientation_histogram(cloud, frequency=3, geometry=True)
        without = orientation_histogram(cloud, frequency=3, geometry=False)
        self.assertIn("polygons", with_geometry)
        self.assertNotIn("polygons", without)
        # Pentagons carry 5 vertices, hexagons 6.
        sizes = with_geometry["sizes"]
        for polygon, size in zip(with_geometry["polygons"], sizes):
            self.assertEqual(len(polygon), size)

    def test_rejects_bad_input(self):
        with self.assertRaises(ValueError):
            orientation_histogram(np.zeros((5, 2)))
        with self.assertRaises(ValueError):
            orientation_histogram(np.zeros((5, 3)))  # all-zero: no directions
        with self.assertRaises(ValueError):
            orientation_histogram(np.ones((5, 3)), weight="nope")
        with self.assertRaises(ValueError):
            orientation_histogram(np.ones((5, 3)), frame="lab")
        with self.assertRaises(ValueError):
            orientation_histogram(np.ones((5, 3)), min_amplitude_quantile=1.0)


class SiteHistogramTests(unittest.TestCase):
    def test_tags_reference_and_element(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "cloud.rmc6f"
            make_cloud_rmc6f(path, [0.03, 0.01, 0.01], supercell=(8, 8, 8), seed=2, element="Nb")
            sites = load_site_displacements(path)
            result = site_orientation_histogram(
                sites, reference_number=1, frequency=4, geometry=False
            )
            self.assertEqual(result["referenceNumber"], 1)
            self.assertEqual(result["element"], "Nb")
            self.assertEqual(result["totalPoints"], 512)
            # The cloud was drawn 3x wider along x: the lobe must sit there.
            self.assertGreater(abs(np.asarray(result["peakDirection"])[0]), 0.8)


if __name__ == "__main__":
    unittest.main()
