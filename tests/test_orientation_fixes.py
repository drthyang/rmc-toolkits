# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regression tests for the 1.0 audit of the displacement-direction engine.

Each class pins one root cause found by the audit (the finding ids are in
the class docstrings). The JS twin of every engine-level assertion lives in
web_app/frontend/src/workers/__tests__/orientationFixes.test.js, with the
same inputs and the same expected values, so the two engines cannot drift.
"""

import unittest

import numpy as np

from rmc_toolkits.orientation import (
    MIN_FREQUENCY,
    assign_cells,
    goldberg_tiling,
    orientation_histogram,
    recommended_frequency,
)


def _cloud(n=200, seed=0):
    return np.random.default_rng(seed).normal(size=(n, 3))


class NonFiniteInputTests(unittest.TestCase):
    """orientation.numerics.5/.15/.29, orientation.parity.12/.18.

    A NaN/inf row used to reach np.cov and fail with LAPACK's 'Eigenvalues did
    not converge' (even in the cartesian frame), while the JS port silently
    returned NaN PCA axes or crashed with a TypeError. Both engines now reject
    such input up front with the same, clear message.
    """

    def test_nan_row_is_rejected_with_a_clear_message(self):
        for frame in ("cartesian", "pca"):
            vectors = _cloud()
            vectors[5, 0] = np.nan
            with self.assertRaisesRegex(ValueError, r"non-finite.*row 5"):
                orientation_histogram(vectors, frame=frame, frequency=3)

    def test_inf_row_is_rejected_with_a_clear_message(self):
        for frame in ("cartesian", "pca"):
            vectors = _cloud()
            vectors[7, 2] = -np.inf
            with self.assertRaisesRegex(ValueError, r"1 non-finite.*row 7"):
                orientation_histogram(vectors, frame=frame, frequency=3)

    def test_non_finite_options_are_rejected(self):
        vectors = _cloud()
        for options in (
            {"smoothing": float("nan")},
            {"smoothing": -1},
            {"target_per_cell": float("nan")},
            {"target_per_cell": 0},
            {"min_amplitude": float("nan")},
            {"frequency": float("nan")},
        ):
            with self.subTest(options=options):
                with self.assertRaises(ValueError):
                    orientation_histogram(vectors, **options)

    def test_recommended_frequency_rejects_invalid_bounds(self):
        with self.assertRaises(ValueError):
            recommended_frequency(1000, max_frequency=0)
        with self.assertRaises(ValueError):
            recommended_frequency(float("nan"))
        with self.assertRaises(ValueError):
            recommended_frequency(1000, target_per_cell=float("inf"))
        self.assertEqual(recommended_frequency(1000, max_frequency=MIN_FREQUENCY), MIN_FREQUENCY)


# Shared verbatim with orientationFixes.test.js: (N, recommended frequency).
# The boundaries sit exactly where 12 * (10 nu^2 + 2) == N.
RECOMMENDED_FREQUENCY_PINS = [
    (0, 1),
    (294, 1),
    (300, 1),
    (503, 1),
    (504, 2),
    (774, 2),
    (1000, 2),
    (1103, 2),
    (1104, 3),
    (12000, 9),
    (12023, 9),
    (12024, 10),
    (10_000_000, 24),
]


class RecommendedFrequencyFloorTests(unittest.TestCase):
    """orientation.physics.10/.25, orientation.numerics.30.

    The docstring promised the largest frequency whose cells still average
    target_per_cell points, but the code rounded to the nearest frequency and
    returned cells averaging as few as 7 points (N=294 -> nu=2, 42 cells).
    """

    def test_never_drops_below_the_target_occupancy(self):
        for n in range(1, 20000, 7):
            frequency = recommended_frequency(n)
            cells = 10 * frequency**2 + 2
            if frequency > MIN_FREQUENCY:
                self.assertGreaterEqual(n / cells, 12, msg=f"N={n} nu={frequency}")
            if frequency < 24:
                finer = 10 * (frequency + 1) ** 2 + 2
                self.assertLess(n / finer, 12, msg=f"N={n}: nu+1 would still hold 12/cell")

    def test_pinned_values_shared_with_the_js_engine(self):
        for n, expected in RECOMMENDED_FREQUENCY_PINS:
            self.assertEqual(recommended_frequency(n), expected, msg=f"N={n}")
        self.assertEqual(recommended_frequency(1000, target_per_cell=5), 4)
        self.assertEqual(recommended_frequency(5000, max_frequency=3), 3)


def _axis_cloud(copies=500, length=0.1):
    """Exactly centrosymmetric: `copies` atoms at each of +/-x, +/-y, +/-z."""
    axes = np.vstack([np.eye(3), -np.eye(3)]) * length
    return np.repeat(axes, copies, axis=0)


def _body_diagonal_cloud(copies=200, length=0.1):
    signs = np.array([[a, b, c] for a in (1, -1) for b in (1, -1) for c in (1, -1)], float)
    return np.repeat(signs * length / np.sqrt(3.0), copies, axis=0)


class CentrosymmetricTieBreakTests(unittest.TestCase):
    """orientation.numerics.28.

    The icosahedron puts 2-fold axes on x, y, z, so at odd nu the directions
    +/-x, +/-y, +/-z sit exactly on a Voronoi boundary (and <111> on triple
    points when 3 does not divide nu). The first-maximum tie-breaks resolved
    +u and -u to cells that are not antipodes, so an exactly centrosymmetric
    cloud reported antipodalAsymmetry = 1.000 and tripped the red flag.
    Assignment is now exactly inversion-equivariant.
    """

    def test_axis_cloud_is_antipodally_symmetric_at_every_frequency(self):
        for frequency in (1, 2, 3, 5, 6, 9, 10, 11):
            with self.subTest(frequency=frequency):
                result = orientation_histogram(_axis_cloud(), frequency=frequency, geometry=False)
                counts = np.asarray(result["counts"])
                np.testing.assert_array_equal(counts, counts[np.asarray(result["antipode"])])
                self.assertEqual(result["antipodalAsymmetry"], 0.0)

    def test_body_diagonal_cloud_is_antipodally_symmetric(self):
        for frequency in (4, 5, 7, 11):
            with self.subTest(frequency=frequency):
                result = orientation_histogram(_body_diagonal_cloud(), frequency=frequency, geometry=False)
                self.assertEqual(result["antipodalAsymmetry"], 0.0)

    def test_assignment_is_inversion_equivariant(self):
        rng = np.random.default_rng(31)
        tiling = goldberg_tiling(5)
        # Random directions plus every exact tie family: the axes, the body
        # diagonals, and the cell-polygon vertices (Voronoi triple points).
        directions = np.vstack([
            rng.normal(size=(2000, 3)),
            np.eye(3),
            _body_diagonal_cloud(copies=1),
            tiling.polygons[:, 0, :],
            tiling.centers,
        ])
        plus = assign_cells(tiling, directions)
        minus = assign_cells(tiling, -directions)
        np.testing.assert_array_equal(minus, tiling.antipode[plus])

    def test_axis_cells_pinned_for_cross_engine_parity(self):
        # Shared verbatim with orientationFixes.test.js.
        tiling = goldberg_tiling(3)
        axes = np.vstack([np.eye(3), -np.eye(3)])
        self.assertEqual(assign_cells(tiling, axes).tolist(), AXIS_CELLS_NU3)


# Cells of +x, +y, +z, -x, -y, -z at nu = 3 (the -u cells are the antipodes).
AXIS_CELLS_NU3 = [36, 15, 16, 86, 71, 62]


if __name__ == "__main__":
    unittest.main()
