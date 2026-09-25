# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regression tests for the 1.0 audit of the displacement-direction engine.

Each class pins one root cause found by the audit (the finding ids are in
the class docstrings). The JS twin of every engine-level assertion lives in
web_app/frontend/src/workers/__tests__/orientationFixes.test.js, with the
same inputs and the same expected values, so the two engines cannot drift.
"""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from rmc_toolkits.orientation import (
    MIN_FREQUENCY,
    assign_cells,
    goldberg_tiling,
    orientation_histogram,
    recommended_frequency,
    site_orientation_histogram,
)
from rmc_toolkits.pca_kde import load_site_displacements


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


# Shared verbatim with orientationFixes.test.js:
# sum over cells c and slots k of (k+1) * (neighbors[c][k]+1) * (c % 97 + 1).
NEIGHBOR_CHECKSUMS = {4: 13107810, 6: 68604861, 10: 524199205, 38: 108147532718}
# First polygon vertex of the two nu = 38 cells whose start used to differ.
POLYGON_STARTS_NU38 = {
    2052: [0.008313449390497447, 0.01625704168659895, 0.9998332836802503],
    4275: [0.9998332836802503, 0.00831344939049745, -0.016257041686598955],
}


def neighbor_checksum(neighbors):
    neighbors = np.asarray(neighbors)
    cells = np.arange(neighbors.shape[0])
    weights = (np.arange(6) + 1)[None, :] * ((cells % 97) + 1)[:, None]
    return int((weights * (neighbors + 1)).sum())


class CanonicalCyclicOrderTests(unittest.TestCase):
    """orientation.parity.17.

    Each cell's neighbours (and polygon vertices) were sorted by atan2 about
    the cell centre; a neighbour on the -e1 ray sits at +pi in one engine and
    -pi in the other (a 1e-17 round-off sign), so the exported ``neighbors``
    rows were cyclically rotated between Python and JS (185 of 1002 rows at
    nu=10), and two polygons at nu=38 started on a different vertex.
    """

    def test_neighbor_rows_start_at_their_smallest_index(self):
        for frequency in (3, 6, 10):
            tiling = goldberg_tiling(frequency)
            for row in tiling.neighbors:
                valid = row[row >= 0]
                self.assertEqual(valid[0], valid.min())

    def test_neighbor_checksums_pinned_for_cross_engine_parity(self):
        for frequency, expected in NEIGHBOR_CHECKSUMS.items():
            self.assertEqual(neighbor_checksum(goldberg_tiling(frequency).neighbors), expected)

    def test_polygon_starts_pinned_for_cross_engine_parity(self):
        tiling = goldberg_tiling(38)
        for cell, vertex in POLYGON_STARTS_NU38.items():
            np.testing.assert_allclose(tiling.polygons[cell, 0], vertex, atol=1e-12)

    def test_cycle_is_still_counter_clockwise(self):
        tiling = goldberg_tiling(6)
        for cell in range(tiling.cell_count):
            size = tiling.sizes[cell]
            polygon = tiling.polygons[cell, :size]
            for i in range(size):
                turn = np.cross(polygon[i], polygon[(i + 1) % size]) @ tiling.centers[cell]
                self.assertGreater(turn, 0.0)


class ElementAllTests(unittest.TestCase):
    """orientation.parity.19 (Python side; the worker fix is in the JS twin)."""

    def test_pooled_result_carries_no_element_key(self):
        rng = np.random.default_rng(41)
        lines = [
            "Supercell dimensions 4 4 4",
            "Lattice vectors (Ang):",
            "32 0 0",
            "0 32 0",
            "0 0 32",
            "Atoms:",
        ]
        atom = 0
        for element, reference, offset in (("Ga", 1, 0.0), ("Se", 2, 0.5)):
            for ix in range(4):
                for iy in range(4):
                    for iz in range(4):
                        atom += 1
                        coord = (np.array([ix, iy, iz]) + offset) / 4 + rng.normal(size=3) * 0.003
                        lines.append(
                            f"{atom} {element} [{reference}] "
                            + " ".join(f"{value:.10f}" for value in coord)
                            + f" {reference} {ix} {iy} {iz}"
                        )
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "two.rmc6f"
            path.write_text("\n".join(lines), encoding="utf-8")
            sites = load_site_displacements(path)
        for element in ("all", "", None):
            result = site_orientation_histogram(sites, element=element, frequency=2, geometry=False)
            self.assertEqual(result["totalPoints"], 128)
            self.assertNotIn("element", result)
        self.assertEqual(
            site_orientation_histogram(sites, element="Se", frequency=2, geometry=False)["element"], "Se"
        )


class TiedPeakTests(unittest.TestCase):
    """orientation.numerics.2/.14, orientation.parity.11/.16, orientation.physics.21.

    Symmetry-equivalent Goldberg cells have mathematically equal solid angles
    that differ in the last bits -- differently in NumPy and JS. With equal
    counts the plain argmax picked a peak by round-off, so the Flask and
    browser engines reported peak directions up to 180 degrees apart. Both
    now take the lowest-index cell within a relative 1e-9 of the maximum and
    report how many cells tie.
    """

    def _tied_cloud(self):
        # nu = 3 cells 2 and 10 are symmetry-equivalent; in both engines the
        # computed area of cell 2 is a few ulp larger, so a plain argmax picks
        # cell 10 when both hold the same count.
        tiling = goldberg_tiling(3)
        return np.vstack(
            [np.repeat(tiling.centers[[2, 10]], 20, axis=0), tiling.centers[[50, 60, 70]]]
        )

    def test_tied_maximum_resolves_to_the_lowest_index(self):
        result = orientation_histogram(self._tied_cloud(), frequency=3, geometry=False)
        tiling = goldberg_tiling(3)
        self.assertEqual(result["peakCell"], 2)
        self.assertEqual(result["peakTieCount"], 2)
        self.assertEqual(result["peakDirection"], tiling.centers[2].tolist())

    def test_untied_maximum_reports_a_single_cell(self):
        rng = np.random.default_rng(1)
        cloud = rng.normal(size=(20000, 3)) * np.array([0.5, 0.1, 0.1])
        result = orientation_histogram(cloud, frequency=8, geometry=False)
        self.assertEqual(result["peakTieCount"], 1)
        self.assertEqual(result["peakCell"], int(np.argmax(result["enhancement"])))


def golden_cloud(n_sphere=900, n_lobe=60, lobe=(0.3, -0.5, 0.81)):
    """Deterministic, RNG-free cloud built identically in the JS suite.

    A Fibonacci-spiral sphere (radius modulated by index) plus a tight lobe of
    `n_lobe` points around `lobe` -- a cross-engine golden input for the
    significance statistics.
    """
    i = np.arange(n_sphere, dtype=float)
    z = 1.0 - (2.0 * i + 1.0) / n_sphere
    r = np.sqrt(1.0 - z * z)
    phi = i * 2.399963229728653
    radius = 0.05 + 0.03 * ((i * 7) % 11) / 11.0
    sphere = np.column_stack([r * np.cos(phi), r * np.sin(phi), z]) * radius[:, None]
    j = np.arange(n_lobe, dtype=float)
    lobe = np.asarray(lobe) / np.linalg.norm(lobe)
    jitter = 0.02 * np.column_stack([np.cos(j * 1.3), np.sin(j * 1.7), np.cos(j * 0.9)])
    lobe_points = (lobe[None, :] + jitter) * 0.12
    return np.vstack([sphere, lobe_points])


def _isotropic_units(rng, n):
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


class PeakSignificanceTests(unittest.TestCase):
    """orientation.numerics.1/.26, orientation.physics.20.

    peakZScore is the local Gaussian z of the chosen cell: the largest of ~C
    correlated values (no look-elsewhere correction) and Gaussian-read at
    expected counts of ~0.2-1, so an exactly isotropic cloud printed a median
    'z = 3.8' at the UI defaults. peakSignificance is the exact Poisson upper
    tail of the peak cell's raw count, Sidak-corrected over all C cells and
    reported as a one-sided normal deviate.
    """

    def test_isotropic_clouds_read_as_not_significant_at_the_ui_defaults(self):
        rng = np.random.default_rng(2024)
        for n in (216, 1000):
            values = [
                orientation_histogram(_isotropic_units(rng, n), frequency=10, smoothing=2, geometry=False)
                for _ in range(150)
            ]
            significance = np.array([v["peakSignificance"] for v in values])
            local = np.array([v["peakZScore"] for v in values])
            # The old readout: most noise maps printed a >= 3 sigma peak.
            self.assertGreater(np.mean(local >= 3), 0.5)
            # The calibrated one: <= 2.3% expected above 2 sigma (Sidak is
            # conservative); allow binomial scatter over 150 draws.
            self.assertLessEqual(np.mean(significance > 2), 0.05, msg=f"N={n}")
            self.assertLessEqual(np.mean(significance > 3), 0.02, msg=f"N={n}")

    def test_value_is_the_sidak_corrected_poisson_tail(self):
        from scipy.stats import norm, poisson

        result = orientation_histogram(golden_cloud(), frequency=6, smoothing=1, geometry=False)
        peak = result["peakCell"]
        count = result["counts"][peak]
        expected = result["expected"][peak]
        self.assertEqual(result["peakCount"], count)
        self.assertAlmostEqual(result["peakExpected"], expected, places=12)
        local = poisson.sf(count - 1, expected)
        self.assertAlmostEqual(result["peakLocalPValue"] / local, 1.0, places=10)
        corrected = -np.expm1(result["cellCount"] * np.log1p(-local))
        self.assertAlmostEqual(result["peakPValue"] / corrected, 1.0, places=8)
        self.assertAlmostEqual(result["peakSignificance"], norm.isf(corrected), places=8)

    def test_a_real_lobe_is_significant(self):
        result = orientation_histogram(golden_cloud(), frequency=6, geometry=False)
        self.assertGreater(result["peakSignificance"], 5.0)
        self.assertLess(result["peakPValue"], 1e-6)

    def test_golden_values_shared_with_the_js_engine(self):
        assert_golden(self, GOLDEN_PEAK)


def assert_golden(case, table):
    for (frequency, smoothing, n_lobe), expected in table.items():
        result = orientation_histogram(
            golden_cloud(n_lobe=n_lobe), frequency=frequency, smoothing=smoothing, geometry=False
        )
        for key, value in expected.items():
            message = f"nu={frequency} s={smoothing} lobe={n_lobe} {key}"
            if isinstance(value, (bool, int)) or value is None:
                case.assertEqual(result[key], value, msg=message)
            else:
                # Relative tolerance -- p-values here reach 1e-142, so an
                # absolute bound would check nothing.
                case.assertLessEqual(abs(result[key] - value), 1e-9 * abs(value), msg=message)


# Shared verbatim with orientationFixes.test.js (GOLDEN_PEAK there):
# (frequency, smoothing, n_lobe) -> expected peak-test fields.
GOLDEN_PEAK = {
    (6, 1, 60): {"peakCell": 181, "peakCount": 39, "peakExpected": 2.5649892741371745,
                 "peakLocalPValue": 3.6267987030691646e-32, "peakPValue": 1.3129011305110377e-29,
                 "peakSignificance": 11.238918217145391},
    (10, 2, 60): {"peakCell": 480, "peakCount": 45, "peakExpected": 0.8960392805524379,
                  "peakLocalPValue": 2.4905767221833515e-59, "peakPValue": 2.4955578756277183e-56,
                  "peakSignificance": 15.770175177825124},
    (None, 0, 60): {"peakCell": 10, "peakCount": 72, "peakExpected": 23.631937108780438,
                    "peakLocalPValue": 1.0240227941836545e-15, "peakPValue": 4.300895735571259e-14,
                    "peakSignificance": 7.460770577845533},
    (10, 2, 0): {"peakCell": 436, "peakCount": 2, "peakExpected": 0.9029912627020413,
                 "peakLocalPValue": 0.22861236738599036, "peakPValue": 1.0,
                 "peakSignificance": -22.629329003077444},
    (2, 0, 6): {"peakCell": 10, "peakCount": 27, "peakExpected": 22.30264064641154,
                "peakLocalPValue": 0.18477705092278982, "peakPValue": 0.999812237582693,
                "peakSignificance": -3.5567124319438927},
}


def _hemisphere_cloud(rng, n):
    u = _isotropic_units(rng, n)
    u[:, 0] = np.abs(u[:, 0])
    return u


class MapSignificanceTests(unittest.TestCase):
    """orientation.numerics.3, orientation.physics.7.

    'map significance N.N sigma' was the RMS of the per-cell z: its null is
    1.00 +/- 1/sqrt(2C) (0.02 at the UI default), not a sigma level, so 1.1
    meant ~4 sigma and a cloud with every atom in one hemisphere printed
    '1.4 sigma'. mapSignificance is Pearson's X^2 = sum z^2 against chi^2 with
    C-1 degrees of freedom, as a one-sided normal deviate.
    """

    def test_value_is_the_pearson_chi_square_tail(self):
        from scipy.stats import chi2, norm

        result = orientation_histogram(golden_cloud(), frequency=6, geometry=False)
        z = np.asarray(result["zScore"])
        self.assertAlmostEqual(result["mapChiSquare"] / float(np.sum(z * z)), 1.0, places=12)
        self.assertEqual(result["mapDegreesOfFreedom"], result["cellCount"] - 1)
        p = chi2.sf(result["mapChiSquare"], result["cellCount"] - 1)
        self.assertAlmostEqual(result["mapPValue"] / p, 1.0, places=9)
        self.assertAlmostEqual(result["mapSignificance"], norm.isf(p), places=8)

    def test_isotropic_null_reads_as_noise_at_the_ui_defaults(self):
        rng = np.random.default_rng(77)
        for n in (216, 1000):
            values = np.array([
                orientation_histogram(_isotropic_units(rng, n), frequency=10, smoothing=2, geometry=False)["mapSignificance"]
                for _ in range(150)
            ])
            self.assertLess(abs(values.mean()), 0.3, msg=f"N={n}")
            self.assertLessEqual(np.mean(values > 2), 0.06, msg=f"N={n}")
            self.assertLessEqual(np.mean(values > 3), 0.02, msg=f"N={n}")

    def test_a_one_sided_cloud_is_overwhelmingly_significant(self):
        rng = np.random.default_rng(5)
        result = orientation_histogram(_hemisphere_cloud(rng, 1000), frequency=10, smoothing=2, geometry=False)
        # The old RMS readout printed this as '1.4 sigma'.
        self.assertLess(result["significance"], 1.5)
        self.assertGreater(result["mapSignificance"], 10.0)

    def test_golden_values_shared_with_the_js_engine(self):
        assert_golden(self, GOLDEN_MAP)


# Shared verbatim with GOLDEN_MAP in orientationFixes.test.js.
GOLDEN_MAP = {
    (6, 1, 60): {"mapChiSquare": 800.1657022469657, "mapDegreesOfFreedom": 361,
                 "mapPValue": 2.594835912106325e-35, "mapSignificance": 12.344903050139939},
    (10, 2, 60): {"mapChiSquare": 2596.980147082603, "mapDegreesOfFreedom": 1001,
                  "mapPValue": 5.119975741811708e-142, "mapSignificance": 25.344856220908564},
    (None, 0, 60): {"mapChiSquare": 109.05527451398147, "mapDegreesOfFreedom": 41,
                    "mapPValue": 4.324724387217442e-08, "mapSignificance": 5.353027939305108},
    (10, 2, 0): {"mapChiSquare": 206.08354289942974, "mapDegreesOfFreedom": 1001,
                 "mapPValue": 1.0, "mapSignificance": -28.039678814235533},
    (2, 0, 6): {"mapChiSquare": 2.8616952082936087, "mapDegreesOfFreedom": 41,
                "mapPValue": 1.0, "mapSignificance": -8.344487440440679},
}


def _skewed_cloud(rng, n, fraction=0.2, shift=0.4, sigma=0.05):
    """Partial one-sided off-centring about the (re-centred) site mean."""
    v = rng.normal(size=(n, 3)) * sigma
    moved = rng.random(n) < fraction
    v[moved, 2] += shift
    return v - v.mean(axis=0)


class AntipodalNullTests(unittest.TestCase):
    """orientation.numerics.4/.13, orientation.physics.6.

    antipodalAsymmetryNull was sqrt(C/(pi N)): a Gaussian-limit, isotropic,
    unconditional mean. It exceeded the statistic's own maximum of 1 at the
    UI default (1.215 for 216 copies), overstated the floor of anisotropic
    sites by 12-68%, and the UI flag compared A with 3x that *mean* instead of
    its spread, so 5-20 sigma asymmetries were never flagged. The null is now
    conditional on the observed antipodal-pair totals T: under inversion
    symmetry each pair splits Bin(T, 1/2), with exact E|2X - T| and variance,
    and the flag is A > null + 3 null SDs.
    """

    def test_null_matches_the_exact_binomial_split(self):
        from scipy.stats import binom

        result = orientation_histogram(golden_cloud(), frequency=4, geometry=False)
        counts = np.asarray(result["counts"])
        antipode = np.asarray(result["antipode"])
        mean = variance = 0.0
        for cell in range(counts.size):
            if cell < antipode[cell]:
                total = counts[cell] + counts[antipode[cell]]
                x = np.arange(total + 1)
                weights = binom.pmf(x, total, 0.5)
                split = np.abs(2 * x - total)
                mean += float((weights * split).sum())
                variance += float((weights * split**2).sum()) - float((weights * split).sum()) ** 2
        used = result["usedPoints"]
        self.assertAlmostEqual(result["antipodalAsymmetryNull"], mean / used, places=12)
        self.assertAlmostEqual(result["antipodalAsymmetryNullSd"], np.sqrt(variance) / used, places=12)
        self.assertAlmostEqual(
            result["antipodalAsymmetryZ"],
            (result["antipodalAsymmetry"] - mean / used) / (np.sqrt(variance) / used),
            places=9,
        )

    def test_null_never_exceeds_the_statistic_bound(self):
        rng = np.random.default_rng(3)
        result = orientation_histogram(_isotropic_units(rng, 216), frequency=10, smoothing=2, geometry=False)
        self.assertLessEqual(result["antipodalAsymmetryNull"], 1.0)
        # The old Gaussian-limit formula gave sqrt(1002 / (pi * 216)) = 1.215.
        self.assertLess(result["antipodalAsymmetryNull"], 0.9)

    def test_centrosymmetric_clouds_are_rarely_flagged(self):
        rng = np.random.default_rng(11)
        for frequency in (10, 5):
            z_values, flags = [], []
            for _ in range(150):
                v = rng.normal(size=(1000, 3)) * np.array([0.2, 0.03, 0.03])
                result = orientation_histogram(v, frequency=frequency, geometry=False)
                z_values.append(result["antipodalAsymmetryZ"])
                flags.append(result["antipodalAsymmetrySignificant"])
            z_values = np.asarray(z_values)
            self.assertLess(abs(z_values.mean()), 0.3, msg=f"nu={frequency}")
            self.assertLess(abs(z_values.std() - 1.0), 0.25, msg=f"nu={frequency}")
            self.assertLessEqual(np.mean(flags), 0.02, msg=f"nu={frequency}")

    def test_real_asymmetry_is_flagged_at_the_ui_default(self):
        rng = np.random.default_rng(12)
        one_sided = orientation_histogram(_hemisphere_cloud(rng, 1000), frequency=10, smoothing=2, geometry=False)
        self.assertTrue(one_sided["antipodalAsymmetrySignificant"])
        self.assertGreater(one_sided["antipodalAsymmetryZ"], 10.0)
        skewed = orientation_histogram(_skewed_cloud(rng, 1000), frequency=10, smoothing=2, geometry=False)
        self.assertTrue(skewed["antipodalAsymmetrySignificant"])
        self.assertGreater(skewed["antipodalAsymmetryZ"], 5.0)
        # The old floor sqrt(C/(pi N)) = 0.565 made both invisible (3x > 1).
        self.assertLess(skewed["antipodalAsymmetry"], 3.0 * np.sqrt(1002 / (np.pi * 1000)))

    def test_all_pairs_singly_occupied_has_no_defined_z(self):
        # One atom in each of 20 cells whose antipodes are empty: every pair
        # total is 1, so |2X - 1| = 1 is deterministic -- A equals its null
        # exactly, the spread is 0 and there is no evidence either way.
        tiling = goldberg_tiling(4)
        cells = np.flatnonzero(np.arange(tiling.cell_count) < tiling.antipode)[:20]
        result = orientation_histogram(tiling.centers[cells] * 0.1, frequency=4, geometry=False)
        self.assertEqual(result["antipodalAsymmetryNullSd"], 0.0)
        self.assertIsNone(result["antipodalAsymmetryZ"])
        self.assertFalse(result["antipodalAsymmetrySignificant"])
        self.assertEqual(result["antipodalAsymmetry"], 1.0)
        self.assertEqual(result["antipodalAsymmetryNull"], 1.0)

    def test_golden_values_shared_with_the_js_engine(self):
        assert_golden(self, GOLDEN_ASYMMETRY)


# Shared verbatim with GOLDEN_ASYMMETRY in orientationFixes.test.js.
GOLDEN_ASYMMETRY = {
    (6, 1, 60): {"antipodalAsymmetry": 0.2, "antipodalAsymmetryNull": 0.3401168600567151,
                 "antipodalAsymmetryNullSd": 0.019345746381926036,
                 "antipodalAsymmetryZ": -7.2427735425923245, "antipodalAsymmetrySignificant": False},
    (None, 0, 60): {"antipodalAsymmetry": 0.08333333333333333,
                    "antipodalAsymmetryNull": 0.11738132743663726,
                    "antipodalAsymmetryNullSd": 0.01946344130131603,
                    "antipodalAsymmetryZ": -1.7493306335813166,
                    "antipodalAsymmetrySignificant": False},
    (10, 2, 0): {"antipodalAsymmetry": 0.19111111111111112,
                 "antipodalAsymmetryNull": 0.5672222222222222,
                 "antipodalAsymmetryNullSd": 0.020957040126511128,
                 "antipodalAsymmetryZ": -17.94676675907692, "antipodalAsymmetrySignificant": False},
    (4, 0, 400): {"antipodalAsymmetry": 0.3415384615384615,
                  "antipodalAsymmetryNull": 0.1752318529176018,
                  "antipodalAsymmetryNullSd": 0.01673666421475673,
                  "antipodalAsymmetryZ": 9.936663990320547, "antipodalAsymmetrySignificant": True},
}


class AnisotropyNullTests(unittest.TestCase):
    """orientation.physics.9/.24, orientation.numerics.27.

    orientationAnisotropy = 3 lambda_1 - 1 was described as '0 for an isotropic
    distribution', but the largest eigenvalue of a finite sample tensor is
    biased upward (0.11 +/- 0.04 for 216 isotropic copies) and no noise
    reference was reported. The engines now report the isotropic expectation
    9 / sqrt(10 pi N_eff) and Bingham's test S = (15 N_eff / 2) sum (l_i - 1/3)^2
    ~ chi^2_5, with its p-value and normal deviate.
    """

    def test_isotropic_expectation_matches_monte_carlo(self):
        rng = np.random.default_rng(21)
        for n in (216, 1000):
            values = [
                orientation_histogram(_isotropic_units(rng, n), frequency=2, geometry=False)
                for _ in range(200)
            ]
            anisotropy = np.array([v["orientationAnisotropy"] for v in values])
            expected = values[0]["orientationAnisotropyNull"]
            self.assertAlmostEqual(expected, 9.0 / np.sqrt(10.0 * np.pi * n), places=12)
            self.assertLess(abs(anisotropy.mean() / expected - 1.0), 0.1, msg=f"N={n}")
            significance = np.array([v["orientationAnisotropySignificance"] for v in values])
            self.assertLessEqual(np.mean(significance > 2), 0.05, msg=f"N={n}")
            self.assertLessEqual(np.mean(significance > 3), 0.015, msg=f"N={n}")

    def test_bingham_statistic_and_p_value(self):
        from scipy.stats import chi2, norm

        rng = np.random.default_rng(22)
        vectors = rng.normal(size=(500, 3)) * np.array([0.12, 0.1, 0.1])
        for weight in ("count", "amplitude2"):
            result = orientation_histogram(vectors, frequency=3, weight=weight, geometry=False)
            amplitude = np.linalg.norm(vectors, axis=1)
            w = np.ones(500) if weight == "count" else amplitude**2
            n_eff = w.sum() ** 2 / (w * w).sum()
            self.assertAlmostEqual(result["orientationEffectivePoints"] / n_eff, 1.0, places=12)
            eigenvalues = np.asarray(result["orientationEigenvalues"])
            statistic = 7.5 * n_eff * float(np.sum((eigenvalues - 1.0 / 3.0) ** 2))
            self.assertAlmostEqual(result["orientationBinghamStatistic"] / statistic, 1.0, places=9)
            p = chi2.sf(statistic, 5)
            self.assertAlmostEqual(result["orientationBinghamPValue"] / p, 1.0, places=8)
            self.assertAlmostEqual(result["orientationAnisotropySignificance"], norm.isf(p), places=7)

    def test_a_modest_harmonic_anisotropy_is_detected(self):
        # sigma ratio 1.25 at 1000 copies: the old map RMS read '1.0 sigma'.
        rng = np.random.default_rng(23)
        vectors = rng.normal(size=(1000, 3)) * np.array([0.125, 0.1, 0.1])
        result = orientation_histogram(vectors, frequency=10, smoothing=2, geometry=False)
        self.assertGreater(result["orientationAnisotropySignificance"], 4.0)
        self.assertGreater(result["orientationAnisotropy"], 3 * result["orientationAnisotropyNull"])

    def test_golden_values_shared_with_the_js_engine(self):
        assert_golden(self, GOLDEN_ANISOTROPY)
        weighted = orientation_histogram(
            golden_cloud(), frequency=6, smoothing=1, weight="amplitude2", geometry=False
        )
        for key, value in GOLDEN_ANISOTROPY_AMPLITUDE2.items():
            self.assertLessEqual(abs(weighted[key] - value), 1e-9 * abs(value), msg=key)


# Shared verbatim with GOLDEN_ANISOTROPY in orientationFixes.test.js.
GOLDEN_ANISOTROPY = {
    (6, 1, 60): {"orientationAnisotropy": 0.12495614840544977, "orientationEffectivePoints": 960.0,
                 "orientationAnisotropyNull": 0.05182412242070032,
                 "orientationBinghamStatistic": 18.73684781045718,
                 "orientationBinghamPValue": 0.0021515346560243643,
                 "orientationAnisotropySignificance": 2.855045282725758},
    (10, 2, 0): {"orientationEffectivePoints": 900.0, "orientationAnisotropyNull": 0.05352372348458313,
                 "orientationBinghamPValue": 0.9999999999999987,
                 "orientationAnisotropySignificance": -7.9084691847240896},
    (4, 0, 400): {"orientationAnisotropy": 0.6150392816231949, "orientationEffectivePoints": 1300.0,
                  "orientationAnisotropyNull": 0.04453442987940475,
                  "orientationBinghamStatistic": 614.6941425286471,
                  "orientationBinghamPValue": 1.3514132534260834e-130,
                  "orientationAnisotropySignificance": 24.286800999672977},
    (2, 0, 6): {"orientationAnisotropy": 0.013271734268148538, "orientationEffectivePoints": 906.0,
                "orientationAnisotropyNull": 0.05334619820786263,
                "orientationBinghamStatistic": 0.19947859577278887,
                "orientationBinghamPValue": 0.9991194624727917,
                "orientationAnisotropySignificance": -3.127820820781734},
}
GOLDEN_ANISOTROPY_AMPLITUDE2 = {
    "orientationEffectivePoints": 725.3165446406116,
    "orientationAnisotropyNull": 0.05962162122187657,
    "orientationBinghamStatistic": 129.85364394492348,
    "orientationBinghamPValue": 2.5564292957448824e-26,
    "orientationAnisotropySignificance": 10.549388830489368,
}


if __name__ == "__main__":
    unittest.main()
