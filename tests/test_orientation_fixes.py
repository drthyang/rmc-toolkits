# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regression tests for the 0.6.0 audit of the displacement-direction engine.

Each class pins one root cause found by the audit (the finding ids are in
the class docstrings). The JS twin of each rule lives in
web_app/frontend/src/workers/__tests__/orientationFixes.test.js. The values
the two engines must share -- the recommendedFrequency table, the tie
families' cells, the neighbour checksums and polygon starts, the 8-well
peaks and every significance statistic on the golden cloud -- are no longer
copied into both suites: tests/generate_orientation_fixture.py writes them
once (orientation_fixture.json, checked here by test_orientation_fixture.py
and in the browser by src/__tests__/orientationParity.test.js).
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

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_orientation_fixture import (
        PEARSON_MOMENTS as PEARSON_MOMENT_CASES,
        golden_cloud,
        pinned_probabilities,
        special_directions,
    )
except ModuleNotFoundError:  # `python -m unittest tests.test_orientation_fixes`
    from tests.generate_orientation_fixture import (
        PEARSON_MOMENTS as PEARSON_MOMENT_CASES,
        golden_cloud,
        pinned_probabilities,
        special_directions,
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


class RecommendedFrequencyFloorTests(unittest.TestCase):
    """orientation.physics.10/.25, orientation.numerics.30.

    The docstring promised the largest frequency whose cells still average
    target_per_cell points, but the code rounded to the nearest frequency and
    returned cells averaging as few as 7 points (N=294 -> nu=2, 42 cells).
    The boundary table both engines share (where 12 (10 nu^2 + 2) == N) is in
    orientation_fixture.json.
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


def _lowest_index_rule(tiling, directions, tolerance=1e-12):
    """Brute-force reference for the tie rule, over all cells.

    Fold into the canonical hemisphere (first non-zero component positive),
    take the lowest cell index whose centre is within ``tolerance`` of the
    best dot product, and map folded directions back through the antipode.
    """
    u = directions / np.linalg.norm(directions, axis=1, keepdims=True)
    x, y, z = u[:, 0], u[:, 1], u[:, 2]
    flip = (x < 0) | ((x == 0) & ((y < 0) | ((y == 0) & (z < 0))))
    folded = np.where(flip[:, None], -u, u)
    dots = folded @ tiling.centers.T
    tied = dots >= dots.max(axis=1, keepdims=True) - tolerance
    cells = np.argmax(tied, axis=1)  # first True = lowest index
    return np.where(flip, tiling.antipode[cells], cells)


class ExactTieResolutionTests(unittest.TestCase):
    """review:4.

    The centrosymmetric fold made each engine inversion-equivariant, but a
    direction on a Voronoi edge or vertex -- <111> when 3 does not divide nu,
    some <110> -- still landed in whichever tied cell the greedy walk reached
    first, which depends on last-bit round-off that differs between NumPy and
    V8: <111> went to different cells in the two engines at 18 of 31
    frequencies, <110> at 9, and an 8-well <111> model had a different peak
    direction in Flask and in the browser. Ties are now resolved
    combinatorially: the lowest cell index among the cells within
    ASSIGN_TIE_TOL of the best dot product. (+x and +y, exact two-way ties
    at nu = 3, now go to the lower index of their pair.) The cells both
    engines give the 26 directions, and the 8-well model's peaks, are in
    orientation_fixture.json.
    """

    def test_special_directions_follow_the_lowest_index_rule(self):
        directions = special_directions()
        for frequency in list(range(1, 31)) + [38]:
            tiling = goldberg_tiling(frequency)
            np.testing.assert_array_equal(
                assign_cells(tiling, directions),
                _lowest_index_rule(tiling, directions),
                err_msg=f"nu={frequency}",
            )

    def test_voronoi_vertices_follow_the_lowest_index_rule(self):
        # Every polygon vertex is equidistant from three cell centres.
        for frequency in (3, 5, 10):
            tiling = goldberg_tiling(frequency)
            vertices = tiling.polygons.reshape(-1, 3)
            np.testing.assert_array_equal(
                assign_cells(tiling, vertices),
                _lowest_index_rule(tiling, vertices),
                err_msg=f"nu={frequency}",
            )

    def test_generic_directions_are_the_plain_nearest_cell(self):
        rng = np.random.default_rng(12)
        directions = rng.normal(size=(20000, 3))
        tiling = goldberg_tiling(10)
        unit = directions / np.linalg.norm(directions, axis=1, keepdims=True)
        np.testing.assert_array_equal(
            assign_cells(tiling, directions), np.argmax(unit @ tiling.centers.T, axis=1)
        )


class CanonicalCyclicOrderTests(unittest.TestCase):
    """orientation.parity.17.

    Each cell's neighbours (and polygon vertices) were sorted by atan2 about
    the cell centre; a neighbour on the -e1 ray sits at +pi in one engine and
    -pi in the other (a 1e-17 round-off sign), so the exported ``neighbors``
    rows were cyclically rotated between Python and JS (185 of 1002 rows at
    nu=10), and two polygons at nu=38 started on a different vertex. The
    neighbour checksums and those two polygon starts, as both engines give
    them, are in orientation_fixture.json.
    """

    def test_neighbor_rows_start_at_their_smallest_index(self):
        for frequency in (3, 6, 10):
            tiling = goldberg_tiling(frequency)
            for row in tiling.neighbors:
                valid = row[row >= 0]
                self.assertEqual(valid[0], valid.min())

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


# golden_cloud() (tests/generate_orientation_fixture.py): a deterministic,
# RNG-free Fibonacci-spiral sphere plus a tight lobe of n_lobe points -- the
# golden fixture's input for the significance statistics, whose values on it
# both engines are held to.


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
            significance = np.array([v["peakSignificance"] for v in values], dtype=float)
            local = np.array([v["peakZScore"] for v in values], dtype=float)
            # review:3 -- a NaN (or missing) field would make every rate check
            # below pass vacuously; require real numbers first.
            self.assertTrue(np.isfinite(significance).all(), msg=f"N={n}")
            self.assertTrue(np.isfinite(local).all(), msg=f"N={n}")
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


def _hemisphere_cloud(rng, n):
    u = _isotropic_units(rng, n)
    u[:, 0] = np.abs(u[:, 0])
    return u


class MapSignificanceTests(unittest.TestCase):
    """orientation.numerics.3, orientation.physics.7, review:1.

    'map significance N.N sigma' was the RMS of the per-cell z: its null is
    1.00 +/- 1/sqrt(2C) (0.02 at the UI default), not a sigma level, so 1.1
    meant ~4 sigma and a cloud with every atom in one hemisphere printed
    '1.4 sigma'. mapSignificance is now Pearson's X^2 = sum z^2 as a one-sided
    normal deviate.

    review:1 -- the first fix referred X^2 to chi^2_{C-1}. Below ~0.1 atom per
    cell X^2 mostly counts atom pairs that share a cell, a skewed, Poisson-like
    quantity, and chi^2_{C-1} (nearly normal at large C) put pure noise above
    3 sigma 15-40x more often than nominal at UI-reachable settings (27
    copies at nu=10, 64 at nu=24). The reference is now a gamma matched to
    the exact multinomial mean, variance and third moment of X^2, and a map
    expecting fewer than MAP_TEST_MIN_PAIRS coincident pairs reports no value.
    """

    def test_null_moments_match_exact_multinomial_enumeration(self):
        from rmc_toolkits.orientation import _pearson_null_moments

        rng = np.random.default_rng(3)
        for used, cells in ((1, 3), (2, 3), (3, 4), (5, 4), (6, 5), (4, 6), (7, 3), (8, 5), (9, 4)):
            p = rng.uniform(0.5, 1.5, cells)
            p /= p.sum()
            exact = _exact_pearson_moments(used, p)
            moments = _pearson_null_moments(used, p)
            np.testing.assert_allclose(moments, exact, rtol=1e-10, atol=1e-10, err_msg=f"N={used} C={cells}")

    def test_null_moments_match_enumeration_on_the_golden_cases(self):
        # The (N, C) cases orientation_fixture.json holds both engines to.
        from rmc_toolkits.orientation import _pearson_null_moments

        for used, cells in PEARSON_MOMENT_CASES:
            p = pinned_probabilities(cells)
            np.testing.assert_allclose(_pearson_null_moments(used, p), _exact_pearson_moments(used, p), rtol=1e-12)

    def test_value_is_the_three_moment_gamma_tail(self):
        from scipy.stats import gamma, norm, pearson3

        tiling = goldberg_tiling(6)
        p = tiling.areas / (4.0 * np.pi)
        for n_lobe in (60, 0):
            result = orientation_histogram(golden_cloud(n_lobe=n_lobe), frequency=6, geometry=False)
            z = np.asarray(result["zScore"])
            self.assertAlmostEqual(result["mapChiSquare"] / float(np.sum(z * z)), 1.0, places=12)
            used = result["usedPoints"]
            self.assertEqual(result["mapDegreesOfFreedom"], result["cellCount"] - 1)
            # Haldane's exact variance, from the areas directly.
            variance = 2.0 * (p.size - 1) + (np.sum(1.0 / p) - p.size**2 - 2 * p.size + 2) / used
            self.assertAlmostEqual(result["mapNullSd"] / np.sqrt(variance), 1.0, places=12)
            self.assertAlmostEqual(
                result["mapExpectedPairs"] / (0.5 * used * (used - 1) * np.sum(p * p)), 1.0, places=12
            )
            # scipy's pearson3 is parameterized by exactly (skewness, mean, SD). Its
            # far upper tail underflows to 0 on older SciPy (1.8), so take the tail
            # from the equivalent shifted gamma: shape 4/g^2, scale sd*g/2, loc
            # mean - 2*sd/g (g > 0 here), which stays accurate on every supported SciPy.
            skew, sd = result["mapNullSkewness"], result["mapNullSd"]
            self.assertGreater(skew, 0.0)
            reference = pearson3(skew, loc=p.size - 1, scale=sd)
            equivalent = gamma(4.0 / skew**2, loc=p.size - 1 - 2.0 * sd / skew, scale=sd * skew / 2.0)
            tail = equivalent.sf(result["mapChiSquare"])
            self.assertAlmostEqual(equivalent.cdf(result["mapChiSquare"]) / reference.cdf(result["mapChiSquare"]), 1.0, places=10)
            self.assertAlmostEqual(result["mapPValue"] / tail, 1.0, places=8)
            # One-sided deviate, from whichever tail is accurate.
            deviate = norm.isf(tail) if tail <= 0.5 else norm.ppf(reference.cdf(result["mapChiSquare"]))
            self.assertAlmostEqual(result["mapSignificance"] / deviate, 1.0, places=7)

    def test_reference_tends_to_chi_square_when_cells_are_well_filled(self):
        from scipy.stats import chi2

        rng = np.random.default_rng(8)
        result = orientation_histogram(_isotropic_units(rng, 20000), frequency=2, geometry=False)
        tail = chi2.sf(result["mapChiSquare"], result["mapDegreesOfFreedom"])
        self.assertLess(abs(result["mapPValue"] - tail), 0.01)
        self.assertLess(abs(result["mapNullSkewness"] / np.sqrt(8.0 / 41.0) - 1.0), 0.01)

    def test_isotropic_null_reads_as_noise_at_the_ui_defaults(self):
        rng = np.random.default_rng(77)
        for n in (216, 1000):
            values = np.array([
                orientation_histogram(_isotropic_units(rng, n), frequency=10, smoothing=2, geometry=False)["mapSignificance"]
                for _ in range(150)
            ])
            self.assertTrue(np.isfinite(values).all(), msg=f"N={n}")
            self.assertLess(abs(values.mean()), 0.3, msg=f"N={n}")
            self.assertLessEqual(np.mean(values > 2), 0.06, msg=f"N={n}")
            self.assertLessEqual(np.mean(values > 3), 0.02, msg=f"N={n}")

    def test_sparse_isotropic_maps_are_calibrated_in_the_engine(self):
        # The reviewer's UI-reachable cases: chi^2_{C-1} put 2.1% (N=27, nu=10)
        # and 2.9% (N=64, nu=24) of pure-noise maps above 3 sigma.
        rng = np.random.default_rng(2718)
        for n, frequency in ((27, 10), (64, 24)):
            values = np.array([
                orientation_histogram(_isotropic_units(rng, n), frequency=frequency, smoothing=2, geometry=False)["mapSignificance"]
                for _ in range(600)
            ])
            self.assertTrue(np.isfinite(values).all(), msg=f"N={n} nu={frequency}")
            self.assertLessEqual(np.sum(values > 3), 4, msg=f"N={n} nu={frequency}")
            self.assertLessEqual(np.mean(values > 2), 0.04, msg=f"N={n} nu={frequency}")

    def test_sparse_reference_is_calibrated_on_exact_multinomial_draws(self):
        # Under isotropy the raw counts are exactly Multinomial(N, Omega/4pi)
        # (assignment is the exact Voronoi partition), so the reference can be
        # checked on direct draws -- 12000 per case, nominal 2.28% / 0.135%
        # (chi^2_{C-1} gives 0.34-2.9% above 3 sigma on these cases).
        from rmc_toolkits.orientation import _map_test

        def threshold(level, n, p):
            # The deviate is monotone in X^2: bisect for the X^2 that reads
            # `level` sigma rather than evaluating the reference per draw.
            low = p.size - 1.0
            high = low + 80.0 * _map_test(low, n, p)["nullSd"]
            for _ in range(80):
                middle = 0.5 * (low + high)
                if _map_test(middle, n, p)["deviate"] > level:
                    high = middle
                else:
                    low = middle
            return high

        rng = np.random.default_rng(99)
        for n, frequency in ((27, 10), (64, 10), (216, 10), (64, 24), (20, 4)):
            p = goldberg_tiling(frequency).areas / (4.0 * np.pi)
            expected = n * p
            chi_square = np.concatenate([
                np.sum((rng.multinomial(n, p / p.sum(), size=2000) - expected) ** 2 / expected, axis=1)
                for _ in range(6)
            ])
            above_three = np.mean(chi_square > threshold(3.0, n, p))
            above_two = np.mean(chi_square > threshold(2.0, n, p))
            self.assertLessEqual(above_three, 0.0025, msg=f"N={n} nu={frequency}")
            self.assertLess(abs(above_two - 0.0228), 0.006, msg=f"N={n} nu={frequency}")

    def test_too_sparse_map_reports_no_significance(self):
        from rmc_toolkits.orientation import MAP_TEST_MIN_PAIRS

        rng = np.random.default_rng(4)
        for n in (1, 2, 8):
            result = orientation_histogram(_isotropic_units(rng, n), frequency=10, geometry=False)
            self.assertLess(result["mapExpectedPairs"], MAP_TEST_MIN_PAIRS)
            self.assertIsNone(result["mapPValue"], msg=f"N={n}")
            self.assertIsNone(result["mapSignificance"], msg=f"N={n}")
            self.assertTrue(np.isfinite(result["mapChiSquare"]))
        populated = orientation_histogram(_isotropic_units(rng, 27), frequency=10, geometry=False)
        self.assertGreaterEqual(populated["mapExpectedPairs"], MAP_TEST_MIN_PAIRS)
        self.assertTrue(np.isfinite(populated["mapSignificance"]))

    def test_a_one_sided_cloud_is_overwhelmingly_significant(self):
        rng = np.random.default_rng(5)
        result = orientation_histogram(_hemisphere_cloud(rng, 1000), frequency=10, smoothing=2, geometry=False)
        # The old RMS readout printed this as '1.4 sigma'.
        self.assertLess(result["significance"], 1.5)
        self.assertGreater(result["mapSignificance"], 10.0)
        # Also in the sparse regime (0.03 atom per cell): 27 atoms confined
        # to a 10-degree cap pile up in a few cells.
        cap = _isotropic_units(rng, 40000)
        cap = cap[cap[:, 2] > np.cos(np.radians(10.0))][:27]
        sparse = orientation_histogram(cap, frequency=10, geometry=False)
        self.assertGreater(sparse["mapSignificance"], 8.0)


def _exact_pearson_moments(used, p):
    """Mean, variance, third central moment of X^2 by full enumeration."""
    from math import factorial

    p = np.asarray(p, dtype=float)

    def compositions(total, parts):
        if parts == 1:
            yield (total,)
            return
        for first in range(total + 1):
            for rest in compositions(total - first, parts - 1):
                yield (first,) + rest

    values, weights = [], []
    for counts in compositions(used, p.size):
        counts = np.asarray(counts)
        coefficient = factorial(used) / np.prod([factorial(int(c)) for c in counts])
        weights.append(coefficient * np.prod(p**counts))
        values.append(np.sum((counts - used * p) ** 2 / (used * p)))
    values, weights = np.asarray(values), np.asarray(weights)
    mean = np.sum(weights * values)
    return (
        mean,
        np.sum(weights * (values - mean) ** 2),
        np.sum(weights * (values - mean) ** 3),
    )


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
            significance = np.array([v["orientationAnisotropySignificance"] for v in values], dtype=float)
            self.assertTrue(np.isfinite(significance).all(), msg=f"N={n}")
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


if __name__ == "__main__":
    unittest.main()
