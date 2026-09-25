# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

import math
import unittest
from itertools import combinations, product
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from rmc_toolkits.triplets import (
    bond_angle_distribution,
    bond_angle_summary,
    bond_angles_from_rmc6f,
)
from rmc_toolkits.triplets_cli import main as triplets_main


ROOT = Path(__file__).resolve().parents[1]
SAMPLE_RMC6F = ROOT / "data" / "5K_try1" / "GaNb4Se8_5K.rmc6f"
# RMCProfile's own TRIPLETS output for the same 5 K configuration.
SAMPLE_TRIPLETS = ROOT / "data" / "5K_try1" / "bonds_hist.pct"

requires_sample = unittest.skipUnless(
    SAMPLE_RMC6F.exists(), "GaNb4Se8 sample data not present in data/ (gitignored)"
)

CUBIC_10 = np.diag([10.0, 10.0, 10.0])


def brute_force_angles(fractional, elements, lattice, triplet, window12, window23, span=2):
    """Independent O(N^2 * images) reference: enumerate every periodic image."""
    end1, apex, end2 = triplet
    fractional = np.asarray(fractional, dtype=float) % 1.0
    lattice = np.asarray(lattice, dtype=float)
    images = [np.asarray(m) for m in product(range(-span, span + 1), repeat=3)]
    angles = []
    for b, element in enumerate(elements):
        if element != apex:
            continue
        center = fractional[b] @ lattice
        bonds1, bonds2 = [], []
        for j, other in enumerate(elements):
            for image in images:
                if j == b and not image.any():
                    continue
                vector = (fractional[j] + image) @ lattice - center
                length = float(np.linalg.norm(vector))
                if other == end1 and window12[0] <= length <= window12[1]:
                    bonds1.append((j, tuple(image), vector, length))
                if other == end2 and window23[0] <= length <= window23[1]:
                    bonds2.append((j, tuple(image), vector, length))
        if end1 == end2:
            # Same end element: one physical triplet {x, B, y} per unordered
            # pair of distinct bond images, kept when either assignment puts
            # one bond in each window (all pairs when the windows are equal).
            in12 = {(bond[0], bond[1]) for bond in bonds1}
            in23 = {(bond[0], bond[1]) for bond in bonds2}
            union = {(bond[0], bond[1]): bond for bond in bonds1 + bonds2}
            keys = sorted(union)
            pairs = [
                (union[x], union[y])
                for x, y in combinations(keys, 2)
                if (x in in12 and y in in23) or (y in in12 and x in in23)
            ]
        else:
            pairs = [(one, two) for one in bonds1 for two in bonds2]
        for one, two in pairs:
            cosine = np.dot(one[2], two[2]) / (one[3] * two[3])
            angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
    return np.sort(np.asarray(angles))


def place(cartesian, lattice=CUBIC_10):
    """Cartesian coordinates -> supercell fractions for test fixtures."""
    return np.asarray(cartesian, dtype=float) @ np.linalg.inv(lattice)


class SingleAngleTests(unittest.TestCase):
    def test_one_triplet_exact_angle(self):
        theta = math.radians(60.0)
        positions = place(
            [[5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.0 + math.cos(theta), 5.0 + math.sin(theta), 5.0]]
        )
        result = bond_angle_distribution(
            positions,
            ["Nb", "Se", "Se"],
            CUBIC_10,
            triplet=("Se", "Nb", "Se"),
            bond12=(0.5, 1.5),
            collect_angles=True,
        )
        self.assertEqual(result.angle_count, 1)
        self.assertAlmostEqual(result.angles[0], 60.0, places=9)
        self.assertEqual(result.apex_count, 1)
        self.assertEqual(result.bond12_count, 2)
        self.assertAlmostEqual(result.mean_length12, 1.0, places=12)

    def test_window_excludes_everything(self):
        positions = place([[5.0, 5.0, 5.0], [6.0, 5.0, 5.0]])
        result = bond_angle_distribution(
            positions,
            ["Nb", "Se"],
            CUBIC_10,
            triplet=("Se", "Nb", "Se"),
            bond12=(2.0, 3.0),
        )
        self.assertEqual(result.angle_count, 0)
        self.assertIsNone(result.mean_angle)
        self.assertIsNone(result.mean_length12)
        self.assertTrue(np.all(result.counts == 0))
        self.assertTrue(np.all(result.density == 0))
        self.assertTrue(np.all(result.sin_corrected == 0))


class OctahedronTests(unittest.TestCase):
    def test_unordered_pair_counting(self):
        center = np.array([5.0, 5.0, 5.0])
        offsets = np.array(
            [[2, 0, 0], [-2, 0, 0], [0, 2, 0], [0, -2, 0], [0, 0, 2], [0, 0, -2]],
            dtype=float,
        )
        positions = place(np.vstack([center, center + offsets]))
        result = bond_angle_distribution(
            positions,
            ["Nb"] + ["O"] * 6,
            CUBIC_10,
            triplet=("O", "Nb", "O"),
            bond12=(1.0, 3.0),
            collect_angles=True,
        )
        self.assertEqual(result.angle_count, 15)
        angles = np.sort(result.angles)
        np.testing.assert_allclose(angles[:12], 90.0, atol=1e-9)
        np.testing.assert_allclose(angles[12:], 180.0, atol=1e-9)
        self.assertAlmostEqual(float(result.density.sum()) * (180.0 / result.counts.size), 1.0)


class PeriodicBoundaryTests(unittest.TestCase):
    def test_bond_through_the_wall(self):
        # Central atom flanked across the supercell boundary: 0.2 A inside the
        # wall on one side, neighbours 0.3 A beyond it on both axes' images.
        positions = np.array(
            [[0.98, 0.5, 0.5], [0.02, 0.5, 0.5], [0.90, 0.5, 0.5]]
        )
        result = bond_angle_distribution(
            positions,
            ["Nb", "Se", "Se"],
            CUBIC_10,
            triplet=("Se", "Nb", "Se"),
            bond12=(0.1, 1.0),
            collect_angles=True,
        )
        # Neighbours at +0.4 A (through the wall) and -0.8 A: one straight angle.
        self.assertEqual(result.angle_count, 1)
        self.assertAlmostEqual(result.angles[0], 180.0, places=9)

    def test_wrap_invariance(self):
        rng = np.random.default_rng(7)
        positions = rng.uniform(size=(30, 3))
        elements = ["Nb" if index % 3 else "Se" for index in range(30)]
        reference = bond_angle_distribution(
            positions, elements, CUBIC_10, triplet=("Se", "Nb", "Se"), bond12=(1.0, 4.0)
        )
        shifted = positions + rng.integers(-3, 4, size=(30, 3))
        moved = bond_angle_distribution(
            shifted, elements, CUBIC_10, triplet=("Se", "Nb", "Se"), bond12=(1.0, 4.0)
        )
        np.testing.assert_array_equal(reference.counts, moved.counts)


class SmallBoxImageTests(unittest.TestCase):
    def test_multiple_images_of_one_neighbour(self):
        # 4 A box, neighbour at +1 A: its -x image sits at 3 A, so a window
        # catching both makes a straight A-B-A' angle through the images.
        lattice = np.diag([4.0, 4.0, 4.0])
        positions = np.array([[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]])
        result = bond_angle_distribution(
            positions,
            ["Nb", "Se"],
            lattice,
            triplet=("Se", "Nb", "Se"),
            bond12=(0.5, 3.5),
            collect_angles=True,
        )
        self.assertEqual(result.bond12_count, 2)
        self.assertEqual(result.angle_count, 1)
        self.assertAlmostEqual(result.angles[0], 180.0, places=9)

    def test_self_image_bonds_of_the_central_element(self):
        # A single atom bonded to its own periodic images: 6 image bonds in a
        # cubic box, forming 90/180 degree angles like an octahedron.
        lattice = np.diag([3.0, 3.0, 3.0])
        positions = np.array([[0.1, 0.2, 0.3]])
        result = bond_angle_distribution(
            positions,
            ["Se"],
            lattice,
            triplet=("Se", "Se", "Se"),
            bond12=(2.5, 3.5),
            collect_angles=True,
        )
        self.assertEqual(result.bond12_count, 6)
        self.assertEqual(result.angle_count, 15)
        angles = np.sort(result.angles)
        np.testing.assert_allclose(angles[:12], 90.0, atol=1e-9)
        np.testing.assert_allclose(angles[12:], 180.0, atol=1e-9)


class BondCountTests(unittest.TestCase):
    """'Bonds' is the physical (undirected) number of bonds.

    triplets.physics.8/20, numerics.15/32: the engines count bonds from each
    central atom, so when the end element is the central element every bond
    is found from both of its ends and the directed count is twice the
    number of bonds.
    """

    def _tetrahedron(self):
        # One Nb4 tetrahedron (edge 3 A) and a Se above each face, in a big box.
        corners = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
        nb = 5.0 + corners * (3.0 / (2 * math.sqrt(2)))
        se = 5.0 - corners * 1.8
        return place(np.vstack([nb, se])), ["Nb"] * 4 + ["Se"] * 4

    def test_homonuclear_bonds_are_counted_once(self):
        positions, elements = self._tetrahedron()
        result = bond_angle_distribution(
            positions, elements, CUBIC_10, triplet=("Nb", "Nb", "Nb"), bond12=(2.5, 3.5)
        )
        self.assertEqual(result.bond12_count, 12)  # B-centred: 4 atoms x 3 neighbours
        self.assertEqual(result.unique_bonds12, 6)  # the tetrahedron's six edges
        self.assertEqual(result.unique_bonds23, 6)
        summary = bond_angle_summary(
            positions, elements, CUBIC_10, triplet=("Nb", "Nb", "Nb"), bond12=(2.5, 3.5)
        )
        self.assertEqual(summary["lengths12"]["count"], 12)
        self.assertEqual(summary["lengths12"]["uniqueBonds"], 6)

    def test_heteronuclear_and_mixed_triplets(self):
        positions, elements = self._tetrahedron()
        # A != B: each Se-Nb bond is found once, from its Nb end.
        hetero = bond_angle_distribution(
            positions, elements, CUBIC_10, triplet=("Se", "Nb", "Se"), bond12=(1.5, 3.2)
        )
        self.assertEqual(hetero.unique_bonds12, hetero.bond12_count)
        # A = B, C != B: halve the Nb-Nb side only.
        mixed = bond_angle_summary(
            positions, elements, CUBIC_10, triplet=("Nb", "Nb", "Se"),
            bond12=(2.5, 3.5), bond23=(1.5, 3.2),
        )
        self.assertEqual(mixed["lengths12"]["uniqueBonds"], 6)
        self.assertEqual(mixed["lengths23"]["uniqueBonds"], mixed["lengths23"]["count"])

    def test_self_image_bonds_pair_up(self):
        # One atom, its six face images in a 3 A cube: three physical bonds
        # (to +a and -a are the same periodic bond), six B-centred vectors.
        result = bond_angle_distribution(
            np.array([[0.1, 0.2, 0.3]]), ["Se"], np.diag([3.0, 3.0, 3.0]),
            triplet=("Se", "Se", "Se"), bond12=(2.5, 3.5),
        )
        self.assertEqual(result.bond12_count, 6)
        self.assertEqual(result.unique_bonds12, 3)

    def test_directed_bonds_come_in_exact_pairs_at_a_window_bound(self):
        # Ideal lattice, bound exactly on a shell: rounding once put a bond
        # inside the window from one end and outside it from the other, so
        # the directed count could be odd (2675 here). Bond vectors are now
        # exact negatives of each other, so the halving is exact.
        cells, a = 5, 3.9
        grid = np.array(list(product(range(cells), repeat=3)), float) / cells
        positions = np.concatenate([grid, (grid + 0.5 / cells) % 1.0])
        lattice = np.diag([a * cells] * 3)
        for rmax in (a, a * math.sqrt(3) / 2, a * math.sqrt(2)):
            summary = bond_angle_summary(
                positions, ["Se"] * len(positions), lattice,
                triplet=("Se", "Se", "Se"), bond12=(0.1, rmax),
            )
            directed = summary["lengths12"]["count"]
            self.assertEqual(directed % 2, 0, f"rmax {rmax}")
            self.assertEqual(summary["lengths12"]["uniqueBonds"] * 2, directed)


class CoincidentAtomTests(unittest.TestCase):
    def test_zero_length_pair_is_never_a_bond(self):
        # Two distinct atoms at bitwise-identical positions with rmin = 0: the
        # zero-length pair has no direction and must be dropped, not fed into
        # a 0/0 NaN angle that silently escapes the histogram.
        positions = np.array(
            [[0.5, 0.5, 0.5], [0.5, 0.5, 0.5], [0.6, 0.5, 0.5], [0.4, 0.5, 0.5]]
        )
        with np.errstate(invalid="raise", divide="raise"):
            result = bond_angle_distribution(
                positions,
                ["Nb", "Se", "Se", "Se"],
                CUBIC_10,
                triplet=("Se", "Nb", "Se"),
                bond12=(0.0, 1.5),
                collect_angles=True,
            )
        self.assertEqual(result.bond12_count, 2)
        self.assertEqual(result.angle_count, 1)
        self.assertEqual(int(result.counts.sum()), result.angle_count)
        self.assertAlmostEqual(result.angles[0], 180.0, places=9)


class SameElementDistinctWindowTests(unittest.TestCase):
    def test_no_zero_angle_from_a_bond_with_itself(self):
        positions = place([[5.0, 5.0, 5.0], [6.0, 5.0, 5.0]])
        result = bond_angle_distribution(
            positions,
            ["Nb", "Se"],
            CUBIC_10,
            triplet=("Se", "Nb", "Se"),
            bond12=(0.5, 1.5),
            bond23=(0.8, 2.0),
        )
        self.assertEqual(result.bond12_count, 1)
        self.assertEqual(result.bond23_count, 1)
        self.assertEqual(result.angle_count, 0)

    def test_overlapping_windows_count_each_triplet_once(self):
        # Both neighbours fall in both windows. {Se, Nb, Se'} is one physical
        # triplet, so it counts once -- not once per (1->2, 2->3) assignment,
        # which used to double every overlap triplet (triplets.physics.21).
        theta = math.radians(120.0)
        positions = place(
            [[5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.0 + math.cos(theta), 5.0 + math.sin(theta), 5.0]]
        )
        result = bond_angle_distribution(
            positions,
            ["Nb", "Se", "Se"],
            CUBIC_10,
            triplet=("Se", "Nb", "Se"),
            bond12=(0.5, 1.5),
            bond23=(0.6, 1.6),
            collect_angles=True,
        )
        self.assertEqual(result.angle_count, 1)
        np.testing.assert_allclose(result.angles, [120.0], atol=1e-9)

    def test_count_is_continuous_as_the_windows_meet(self):
        # Nudging the B-C bound by 1e-4 A past every bond must not change a
        # thing: the distinct-window rule reduces to the shared-window one.
        center = np.array([5.0, 5.0, 5.0])
        offsets = np.array(
            [[2, 0, 0], [-2, 0, 0], [0, 2, 0], [0, -2, 0], [0, 0, 2], [0, 0, -2]],
            dtype=float,
        )
        positions = place(np.vstack([center, center + offsets]))
        kwargs = dict(triplet=("O", "Nb", "O"), bond12=(1.0, 3.0), bin_width=5.0)
        shared = bond_angle_distribution(positions, ["Nb"] + ["O"] * 6, CUBIC_10, **kwargs)
        nudged = bond_angle_distribution(
            positions, ["Nb"] + ["O"] * 6, CUBIC_10, bond23=(1.0, 3.0001), **kwargs
        )
        self.assertEqual(shared.angle_count, 15)
        self.assertEqual(nudged.angle_count, 15)
        np.testing.assert_array_equal(nudged.counts, shared.counts)

    def test_disjoint_windows_pair_across_shells(self):
        # Short and long Nb-Se bonds (1 A and 2 A): with disjoint windows only
        # short-long pairs qualify, each counted once.
        positions = place(
            [[5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.0, 6.0, 5.0], [3.0, 5.0, 5.0], [5.0, 5.0, 7.0]]
        )
        result = bond_angle_distribution(
            positions,
            ["Nb", "Se", "Se", "Se", "Se"],
            CUBIC_10,
            triplet=("Se", "Nb", "Se"),
            bond12=(0.5, 1.5),
            bond23=(1.5, 2.5),
            collect_angles=True,
        )
        self.assertEqual(result.bond12_count, 2)
        self.assertEqual(result.bond23_count, 2)
        self.assertEqual(result.angle_count, 4)
        np.testing.assert_allclose(np.sort(result.angles), [90.0, 90.0, 90.0, 180.0], atol=1e-9)


class TriclinicBruteForceTests(unittest.TestCase):
    """Exact agreement with an independent all-images reference in a skewed cell."""

    LATTICE = np.array([[6.0, 0.0, 0.0], [3.0, 5.0, 0.0], [1.0, 1.0, 7.0]])

    def _random_configuration(self, seed=11, count=40):
        rng = np.random.default_rng(seed)
        positions = rng.uniform(size=(count, 3))
        elements = ["Se" if value < 0.6 else "Nb" for value in rng.uniform(size=count)]
        return positions, elements

    def assert_matches_brute_force(self, triplet, window12, window23):
        positions, elements = self._random_configuration()
        result = bond_angle_distribution(
            positions,
            elements,
            self.LATTICE,
            triplet=triplet,
            bond12=window12,
            bond23=window23,
            collect_angles=True,
        )
        expected = brute_force_angles(
            positions,
            elements,
            self.LATTICE,
            triplet,
            window12,
            window23 or window12,
        )
        self.assertEqual(result.angle_count, expected.size)
        np.testing.assert_allclose(np.sort(result.angles), expected, atol=1e-8)

    def test_shared_end_windows(self):
        self.assert_matches_brute_force(("Se", "Nb", "Se"), (1.0, 3.4), None)

    def test_distinct_end_windows(self):
        self.assert_matches_brute_force(("Se", "Nb", "Nb"), (1.0, 3.0), (1.5, 3.4))

    def test_same_element_everywhere(self):
        self.assert_matches_brute_force(("Se", "Se", "Se"), (1.0, 3.2), None)

    def test_same_end_element_distinct_windows(self):
        self.assert_matches_brute_force(("Se", "Nb", "Se"), (1.0, 3.0), (1.5, 3.4))

    def test_same_element_everywhere_distinct_windows(self):
        self.assert_matches_brute_force(("Se", "Se", "Se"), (1.0, 3.0), (2.0, 3.4))


class IsotropicReferenceTests(unittest.TestCase):
    def test_sin_correction_is_flat_for_random_directions(self):
        rng = np.random.default_rng(3)
        count = 1200
        directions = rng.normal(size=(count, 3))
        directions /= np.linalg.norm(directions, axis=1)[:, None]
        lattice = np.diag([200.0, 200.0, 200.0])
        cartesian = np.vstack([[100.0, 100.0, 100.0], 100.0 + 2.0 * directions])
        positions = cartesian / 200.0  # diagonal lattice; avoids an Accelerate matmul quirk
        result = bond_angle_distribution(
            positions,
            ["Nb"] + ["Se"] * count,
            lattice,
            triplet=("Se", "Nb", "Se"),
            bond12=(1.0, 3.0),
            bin_width=15.0,
        )
        self.assertEqual(result.angle_count, count * (count - 1) // 2)
        np.testing.assert_allclose(result.sin_corrected, 1.0, atol=0.05)


def ideal_perovskite(cells=3, a=3.905, shift=0.0):
    """Undisplaced cubic SrTiO3 supercell (fractions (i + x) / n, as a CIF-built
    RMCProfile start configuration stores them): every angle is a symmetry
    angle, 60/90/120/180 deg up to float noise."""
    basis = [("Sr", (0, 0, 0)), ("Ti", (0.5, 0.5, 0.5)), ("O", (0.5, 0.5, 0)),
             ("O", (0.5, 0, 0.5)), ("O", (0, 0.5, 0.5))]
    positions, elements = [], []
    for i, j, k in product(range(cells), repeat=3):
        for element, (x, y, z) in basis:
            positions.append([(i + x) / cells + shift, (j + y) / cells + shift, (k + z) / cells + shift])
            elements.append(element)
    return np.asarray(positions), elements, np.diag([a * cells] * 3)


class IdealConfigurationTests(unittest.TestCase):
    """Symmetry angles sit exactly on bin edges; float noise must not split them.

    triplets.physics.6, numerics.12, parity.24, numerics.30: an undisplaced
    configuration puts each symmetry angle a few ulp either side of an edge,
    so a whole symmetry class used to split between two bins -- differently
    per engine (libm vs V8 acos) and under a rigid shift.
    """

    def test_each_symmetry_angle_lands_in_one_bin(self):
        positions, elements, lattice = ideal_perovskite()
        for width in (1.0, 0.5, 5.0):
            result = bond_angle_distribution(
                positions, elements, lattice, triplet=("O", "Sr", "O"),
                bond12=(2.0, 3.0), bin_width=width,
            )
            occupied = {round(float(result.bin_edges[k]), 6): int(c)
                        for k, c in enumerate(result.counts) if c}
            # Half-open bins: an exact edge angle belongs to the bin it starts;
            # 180 deg to the last bin.
            last = round(180.0 - float(result.bin_edges[1]), 6)
            self.assertEqual(
                occupied, {60.0: 648, 90.0: 324, 120.0: 648, last: 162}, f"width {width}"
            )

    def test_rigid_shift_changes_nothing(self):
        reference = None
        for shift in (0.0, 0.001, 0.0123, 0.5):
            positions, elements, lattice = ideal_perovskite(shift=shift)
            counts = bond_angle_distribution(
                positions, elements, lattice, triplet=("O", "O", "O"), bond12=(2.0, 3.0)
            ).counts
            if reference is None:
                reference = counts
            np.testing.assert_array_equal(counts, reference, f"shift {shift}")

    def test_summary_bins_like_the_distribution(self):
        positions, elements, lattice = ideal_perovskite()
        kwargs = dict(triplet=("O", "Sr", "O"), bond12=(2.0, 3.0))
        summary = bond_angle_summary(positions, elements, lattice, **kwargs)
        result = bond_angle_distribution(positions, elements, lattice, **kwargs)
        self.assertEqual(summary["counts"], result.counts.tolist())

    def test_bins_match_numpy_histogram_off_the_edges(self):
        from rmc_toolkits.triplets import _angle_bins

        rng = np.random.default_rng(1)
        for nbins in (1, 7, 23, 180, 3600):
            edges = np.linspace(0.0, 180.0, nbins + 1)
            angles = np.concatenate([
                rng.uniform(0.0, 180.0, 5000),
                np.nextafter(edges, -np.inf)[1:],  # one ulp below each edge
                edges[1:-1] + 1e-6,                # just above, beyond the snap
                [0.0, 180.0],
            ])
            angles = np.clip(angles, 0.0, 180.0)
            expected, _ = np.histogram(angles, bins=nbins, range=(0.0, 180.0))
            ulp_below = np.nextafter(edges, -np.inf)[1:]
            snapped = np.isin(angles, ulp_below[:-1])
            # Away from the edges: numpy's assignment, exactly.
            got = np.bincount(_angle_bins(angles[~snapped], nbins), minlength=nbins)
            reference, _ = np.histogram(angles[~snapped], bins=nbins, range=(0.0, 180.0))
            np.testing.assert_array_equal(got, reference, f"nbins {nbins}")
            # One ulp below an interior edge: snapped up into that edge's bin.
            np.testing.assert_array_equal(
                _angle_bins(ulp_below[:-1], nbins), np.arange(1, nbins)
            )
            self.assertEqual(int(expected.sum()), angles.size)

    def test_snap_is_limited_to_float_noise(self):
        # 1e-7 deg off an edge is a real (if tiny) displacement, not noise:
        # it bins by its value, below the edge.
        theta = math.radians(60.0 - 1e-7)
        positions = place(
            [[5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.0 + math.cos(theta), 5.0 + math.sin(theta), 5.0]]
        )
        result = bond_angle_distribution(
            positions, ["Nb", "Se", "Se"], CUBIC_10, triplet=("Se", "Nb", "Se"), bond12=(0.5, 1.5)
        )
        self.assertEqual(int(np.flatnonzero(result.counts)[0]), 59)


class SinCorrectionIdentityTests(unittest.TestCase):
    """What sin_corrected is (triplets.physics.19, physics.2).

    The bin-integral reference (cos a - cos b) / 2 equals sin(c) sin(w/2)
    exactly, so the curve is the bin-centre 1/sin(c) correction times the
    constant 1/sin(w/2) -- finite either way at the end bins -- and
    RMCProfile's TRIPLETS norm/sin(theta) is it times sin(w/2)/w_deg.
    """

    def test_bin_integral_is_centre_sine_times_a_constant(self):
        rng = np.random.default_rng(4)
        positions = rng.uniform(size=(60, 3))
        elements = ["Se" if index % 3 else "Nb" for index in range(60)]
        for width in (1.0, 5.0, 15.0):
            result = bond_angle_distribution(
                positions, elements, CUBIC_10, triplet=("Se", "Nb", "Se"),
                bond12=(0.5, 5.0), bin_width=width,
            )
            half = math.radians(width) / 2
            centre_form = (result.counts / result.angle_count) / (
                np.sin(np.radians(result.bin_centers)) * math.sin(half)
            )
            np.testing.assert_allclose(result.sin_corrected, centre_form, rtol=1e-12)
            # Both forms are finite in the 0 and 180 degree bins.
            self.assertTrue(np.all(np.isfinite(centre_form)))
            # Density / sin(centre) -- RMCProfile's norm/sin(theta) -- is the
            # same curve times sin(w/2) / w, i.e. ~pi/360 at small widths.
            rmcprofile = result.density / np.sin(np.radians(result.bin_centers))
            np.testing.assert_allclose(
                rmcprofile, result.sin_corrected * math.sin(half) / width, rtol=1e-12
            )


@unittest.skipUnless(
    SAMPLE_TRIPLETS.exists() and SAMPLE_RMC6F.exists(),
    "RMCProfile TRIPLETS output for the 5 K sample not present in data/ (gitignored)",
)
class RmcProfileTripletsTests(unittest.TestCase):
    """Cross-check against RMCProfile's TRIPLETS on the same configuration.

    bonds_hist.pct: rmax 3.5 A for every pair, 1000 bins of 0.18 deg,
    columns theta, norm/sin(theta), norm, un_norm; type 1 = Ga, 2 = Nb,
    3 = Se.
    """

    WIDTH = 0.18

    @classmethod
    def setUpClass(cls):
        cls.lines = SAMPLE_TRIPLETS.read_text(encoding="utf-8").splitlines()

    def section(self, tag):
        start = next(n for n, line in enumerate(self.lines) if line.strip() == tag)
        rows = self.lines[start + 4:start + 4 + 1000]
        return np.array([[float(value) for value in row.split()] for row in rows])

    def test_totals_density_and_scale_match(self):
        factor = math.sin(math.radians(self.WIDTH) / 2) / self.WIDTH
        for tag, triplet, total in [
            ("b323", ("Se", "Nb", "Se"), 239326),
            ("b222", ("Nb", "Nb", "Nb"), 47078),
            ("b313", ("Se", "Ga", "Se"), 24132),
            ("b232", ("Nb", "Se", "Nb"), 95731),
            ("b322", ("Se", "Nb", "Nb"), 284483),
        ]:
            reference = self.section(tag)
            result = bond_angles_from_rmc6f(
                SAMPLE_RMC6F, triplet=triplet, bond12=(0.0, 3.5), bin_width=self.WIDTH
            )
            self.assertEqual(result.angle_count, total, tag)
            self.assertEqual(int(reference[:, 3].sum()), total, tag)
            # RMCProfile bins in single precision: a few angles sit across a
            # neighbouring edge, never further.
            drift = np.abs(np.cumsum(result.counts) - np.cumsum(reference[:, 3])).max()
            self.assertLessEqual(drift, 5, tag)
            same = (result.counts == reference[:, 3]) & (result.counts > 0)
            np.testing.assert_allclose(result.density[same], reference[same, 2], rtol=1e-5)
            np.testing.assert_allclose(
                reference[same, 1], result.sin_corrected[same] * factor, rtol=1e-4
            )


class HistogramConventionTests(unittest.TestCase):
    def test_bin_count_rounds_half_up_like_the_js_port(self):
        # 180/8 = 22.5 exactly: banker's rounding would give 22 bins while
        # JavaScript's Math.round gives 23 — the engines must agree, so the
        # rule is floor(x + 0.5) in both (see _bin_count).
        positions = place([[5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.0, 6.0, 5.0]])
        for width, expected in [(8.0, 23), (40.0, 5), (1.0, 180)]:
            result = bond_angle_distribution(
                positions,
                ["Nb", "Se", "Se"],
                CUBIC_10,
                triplet=("Se", "Nb", "Se"),
                bond12=(0.5, 1.5),
                bin_width=width,
            )
            self.assertEqual(result.counts.size, expected, f"width {width}")

    def test_bin_width_rounds_to_exact_tiling(self):
        positions = place([[5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.0, 6.0, 5.0]])
        result = bond_angle_distribution(
            positions,
            ["Nb", "Se", "Se"],
            CUBIC_10,
            triplet=("Se", "Nb", "Se"),
            bond12=(0.5, 1.5),
            bin_width=2.5,
        )
        self.assertEqual(result.counts.size, 72)
        self.assertAlmostEqual(result.bin_edges[1] - result.bin_edges[0], 2.5)
        self.assertEqual(int(result.counts.sum()), result.angle_count)

    def test_isotropic_fractions_sum_to_one(self):
        positions = place([[5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.0, 6.0, 5.0]])
        result = bond_angle_distribution(
            positions,
            ["Nb", "Se", "Se"],
            CUBIC_10,
            triplet=("Se", "Nb", "Se"),
            bond12=(0.5, 1.5),
        )
        edges = np.radians(result.bin_edges)
        fractions = (np.cos(edges[:-1]) - np.cos(edges[1:])) / 2.0
        self.assertAlmostEqual(float(fractions.sum()), 1.0, places=12)
        # The one angle sits at 90 degrees; its sin-corrected value is the
        # count fraction over the isotropic fraction of that bin.
        bin_index = int(np.flatnonzero(result.counts)[0])
        self.assertAlmostEqual(
            result.sin_corrected[bin_index], 1.0 / fractions[bin_index], places=9
        )


class SummaryTests(unittest.TestCase):
    """The JSON payload agrees with the dataclass result and its own totals."""

    LATTICE = np.array([[6.0, 0.0, 0.0], [3.0, 5.0, 0.0], [1.0, 1.0, 7.0]])

    def _configuration(self):
        rng = np.random.default_rng(11)
        positions = rng.uniform(size=(40, 3))
        elements = ["Se" if value < 0.6 else "Nb" for value in rng.uniform(size=40)]
        return positions, elements

    def test_summary_matches_distribution(self):
        positions, elements = self._configuration()
        kwargs = dict(triplet=("Se", "Nb", "Nb"), bond12=(1.0, 3.0), bond23=(1.5, 3.4))
        summary = bond_angle_summary(positions, elements, self.LATTICE, **kwargs)
        result = bond_angle_distribution(positions, elements, self.LATTICE, **kwargs)
        np.testing.assert_array_equal(summary["counts"], result.counts)
        np.testing.assert_allclose(summary["density"], result.density)
        np.testing.assert_allclose(summary["sinCorrected"], result.sin_corrected)
        self.assertEqual(summary["angleCount"], result.angle_count)
        self.assertEqual(summary["apexCount"], result.apex_count)
        self.assertEqual(summary["lengths12"]["count"], result.bond12_count)
        self.assertEqual(summary["lengths23"]["count"], result.bond23_count)
        self.assertAlmostEqual(summary["lengths12"]["meanLength"], result.mean_length12)
        self.assertFalse(summary["sharedEnds"])

    def test_summary_internal_totals(self):
        positions, elements = self._configuration()
        summary = bond_angle_summary(
            positions, elements, self.LATTICE, triplet=("Se", "Nb", "Se"), bond12=(1.0, 3.2)
        )
        self.assertTrue(summary["sharedEnds"])
        self.assertIsNone(summary["lengths23"])
        # Length-histogram bins tile the window, so their counts sum to the
        # bond count; the coordination histogram partitions the apex atoms.
        self.assertEqual(sum(summary["lengths12"]["counts"]), summary["lengths12"]["count"])
        self.assertEqual(sum(summary["coordination"]), summary["apexCount"])
        expected_pairs = sum(
            count * n * (n - 1) // 2 for n, count in enumerate(summary["coordination"])
        )
        self.assertEqual(summary["angleCount"], expected_pairs)

    def test_summary_is_json_safe(self):
        import json

        positions, elements = self._configuration()
        summary = bond_angle_summary(
            positions, elements, self.LATTICE, triplet=("Se", "Nb", "Se"), bond12=(1.0, 3.2)
        )
        json.dumps(summary)  # raises on any lingering numpy type


class WorkBudgetTests(unittest.TestCase):
    """The angle list is streamed, and the app budget is exact and up front.

    triplets.physics.1/7/18, numerics.11/31, parity.28, backend.api.4 and
    backend.cache.15: the app's rmax cap bounds the neighbour search, but the
    angle count grows ~rmax^6 and both engines used to hold every angle.
    """

    def _cloud(self, count=1000, box=20.8, seed=5):
        rng = np.random.default_rng(seed)
        return rng.uniform(size=(count, 3)), ["Se"] * count, np.diag([box] * 3)

    def test_budget_rejects_before_forming_angles(self):
        positions, elements, lattice = self._cloud()
        kwargs = dict(triplet=("Se", "Se", "Se"), bond12=(0.5, 3.0))
        exact = bond_angle_summary(positions, elements, lattice, **kwargs)["angleCount"]
        self.assertGreater(exact, 1000)
        # At the budget: allowed. One angle over: a clear error naming both.
        at_budget = bond_angle_summary(positions, elements, lattice, max_angles=exact, **kwargs)
        self.assertEqual(at_budget["angleCount"], exact)
        for engine in (bond_angle_summary, bond_angle_distribution):
            with self.assertRaisesRegex(ValueError, rf"{exact:,} angles.*{exact - 1:,}"):
                engine(positions, elements, lattice, max_angles=exact - 1, **kwargs)

    def test_budget_counts_exactly_for_every_counting_rule(self):
        positions, elements, lattice = self._cloud(count=300, box=12.0)
        elements = ["Se" if index % 3 else "Nb" for index in range(300)]
        for triplet, bond12, bond23 in [
            (("Se", "Nb", "Se"), (1.0, 3.0), None),
            (("Se", "Nb", "Nb"), (1.0, 3.0), (1.5, 3.4)),
            (("Se", "Nb", "Se"), (1.0, 3.0), (1.5, 3.4)),
            (("Se", "Se", "Se"), (1.0, 2.0), (2.5, 3.4)),
        ]:
            kwargs = dict(triplet=triplet, bond12=bond12, bond23=bond23)
            exact = bond_angle_distribution(positions, elements, lattice, **kwargs).angle_count
            bond_angle_distribution(positions, elements, lattice, max_angles=exact, **kwargs)
            with self.assertRaises(ValueError):
                bond_angle_distribution(
                    positions, elements, lattice, max_angles=exact - 1, **kwargs
                )

    def test_app_budget_is_the_shared_constant(self):
        from rmc_toolkits import triplets

        self.assertEqual(triplets.APP_MAX_ANGLES, 50_000_000)

    def test_memory_is_bounded_by_the_chunk_not_the_angle_count(self):
        import tracemalloc

        # ~1.6e6 angles: materializing them (the old engine) peaks at
        # ~200 B/angle, i.e. over 300 MB; streamed it stays near one chunk.
        positions, elements, lattice = self._cloud()
        tracemalloc.start()
        try:
            summary = bond_angle_summary(
                positions, elements, lattice, triplet=("Se", "Se", "Se"), bond12=(0.5, 5.0)
            )
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        self.assertGreater(summary["angleCount"], 1_000_000)
        self.assertLess(peak, 96 * 2**20)

    def test_refusal_stores_neither_bonds_nor_angles(self):
        import tracemalloc

        # ~8e7 angles from ~4e5 bonds: the budgeted path counts them with a
        # search that stores nothing, so refusing costs one search block.
        positions, elements, lattice = self._cloud()
        tracemalloc.start()
        try:
            with self.assertRaisesRegex(ValueError, "over the limit of 50,000,000"):
                bond_angle_summary(
                    positions, elements, lattice, triplet=("Se", "Se", "Se"),
                    bond12=(0.5, 9.5), max_angles=50_000_000,
                )
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        self.assertLess(peak, 48 * 2**20)

    def test_search_blocking_is_invisible(self):
        from unittest import mock

        from rmc_toolkits import triplets

        positions, elements, lattice = self._cloud(count=200, box=10.0)
        elements = ["Se" if index % 3 else "Nb" for index in range(200)]
        for triplet, bond23 in [(("Se", "Nb", "Se"), None), (("Se", "Nb", "Se"), (1.5, 3.4)),
                                (("Se", "Nb", "Nb"), (1.5, 3.4))]:
            kwargs = dict(triplet=triplet, bond12=(1.0, 3.0), bond23=bond23)
            whole = bond_angle_summary(positions, elements, lattice, **kwargs)
            with mock.patch.object(triplets, "SEARCH_CHUNK", 5):
                blocked = bond_angle_summary(
                    positions, elements, lattice, max_angles=whole["angleCount"], **kwargs
                )
            self.assertEqual(blocked["counts"], whole["counts"])
            self.assertEqual(blocked["coordination"], whole["coordination"])
            self.assertEqual(blocked["lengths12"]["counts"], whole["lengths12"]["counts"])
            self.assertAlmostEqual(blocked["meanAngle"], whole["meanAngle"], places=9)

    def test_chunking_is_invisible(self):
        from unittest import mock

        from rmc_toolkits import triplets

        positions, elements, lattice = self._cloud(count=200, box=10.0)
        elements = ["Se" if index % 3 else "Nb" for index in range(200)]
        for triplet, bond23 in [(("Se", "Nb", "Se"), None), (("Se", "Nb", "Se"), (1.5, 3.4)),
                                (("Se", "Nb", "Nb"), (1.5, 3.4))]:
            kwargs = dict(triplet=triplet, bond12=(1.0, 3.0), bond23=bond23, collect_angles=True)
            whole = bond_angle_distribution(positions, elements, lattice, **kwargs)
            with mock.patch.object(triplets, "PAIR_CHUNK", 7):
                chunked = bond_angle_distribution(positions, elements, lattice, **kwargs)
            self.assertGreater(whole.angle_count, 100)
            self.assertEqual(chunked.angle_count, whole.angle_count)
            np.testing.assert_array_equal(chunked.counts, whole.counts)
            np.testing.assert_array_equal(np.sort(chunked.angles), np.sort(whole.angles))
            self.assertAlmostEqual(chunked.mean_angle, float(np.mean(whole.angles)), places=9)
            self.assertAlmostEqual(chunked.std_angle, float(np.std(whole.angles)), places=9)


class ValidationTests(unittest.TestCase):
    def setUp(self):
        self.positions = place([[5.0, 5.0, 5.0], [6.0, 5.0, 5.0]])
        self.elements = ["Nb", "Se"]

    def test_unknown_element_lists_available(self):
        with self.assertRaisesRegex(ValueError, "available: Nb, Se"):
            bond_angle_distribution(
                self.positions,
                self.elements,
                CUBIC_10,
                triplet=("O", "Nb", "Se"),
                bond12=(1.0, 2.0),
            )

    def test_bad_windows(self):
        for window in [(2.0, 1.0), (-1.0, 2.0), (1.0, 1.0), (1.0, float("inf"))]:
            with self.assertRaises(ValueError):
                bond_angle_distribution(
                    self.positions,
                    self.elements,
                    CUBIC_10,
                    triplet=("Se", "Nb", "Se"),
                    bond12=window,
                )

    def test_singular_lattice(self):
        with self.assertRaisesRegex(ValueError, "singular"):
            bond_angle_distribution(
                self.positions,
                self.elements,
                np.array([[1.0, 0, 0], [2.0, 0, 0], [0, 0, 1.0]]),
                triplet=("Se", "Nb", "Se"),
                bond12=(1.0, 2.0),
            )

    def test_shape_and_length_mismatches(self):
        with self.assertRaises(ValueError):
            bond_angle_distribution(
                self.positions,
                ["Nb"],
                CUBIC_10,
                triplet=("Se", "Nb", "Se"),
                bond12=(1.0, 2.0),
            )
        with self.assertRaises(ValueError):
            bond_angle_distribution(
                self.positions,
                self.elements,
                CUBIC_10,
                triplet=("Se", "Nb"),
                bond12=(1.0, 2.0),
            )
        with self.assertRaises(ValueError):
            bond_angle_distribution(
                self.positions,
                self.elements,
                CUBIC_10,
                triplet=("Se", "Nb", "Se"),
                bond12=(1.0, 2.0),
                bin_width=0.0,
            )

    def test_case_insensitive_symbols(self):
        result = bond_angle_distribution(
            self.positions,
            ["NB", "se"],
            CUBIC_10,
            triplet=("SE", "nb", "Se"),
            bond12=(0.5, 1.5),
        )
        self.assertEqual(result.triplet, ("Se", "Nb", "Se"))
        self.assertEqual(result.bond12_count, 1)


def write_rmc6f(path: Path, atom_lines, *, supercell=(1, 1, 1), lattice=CUBIC_10):
    header = [
        f"Supercell dimensions: {supercell[0]} {supercell[1]} {supercell[2]}",
        "Lattice vectors (Ang):",
        " ".join(f"{value:.6f}" for value in lattice[0]),
        " ".join(f"{value:.6f}" for value in lattice[1]),
        " ".join(f"{value:.6f}" for value in lattice[2]),
        "Atoms:",
    ]
    path.write_text("\n".join(header + atom_lines) + "\n", encoding="utf-8")


class LoaderTests(unittest.TestCase):
    def test_loader_matches_direct_call(self):
        theta = math.radians(75.0)
        cartesian = np.array(
            [[5.0, 5.0, 5.0], [6.2, 5.0, 5.0], [5.0 + 1.2 * math.cos(theta), 5.0 + 1.2 * math.sin(theta), 5.0]]
        )
        positions = place(cartesian)
        lines = [
            f"{index + 1} {element} [1] {p[0]:.12f} {p[1]:.12f} {p[2]:.12f} {index + 1} 0 0 0"
            for index, (element, p) in enumerate(zip(["Nb", "Se", "Se"], positions))
        ]
        with TemporaryDirectory() as scratch:
            config = Path(scratch) / "tiny.rmc6f"
            write_rmc6f(config, lines)
            from_file = bond_angles_from_rmc6f(
                config, triplet=("Se", "Nb", "Se"), bond12=(0.5, 1.5), collect_angles=True
            )
        direct = bond_angle_distribution(
            positions,
            ["Nb", "Se", "Se"],
            CUBIC_10,
            triplet=("Se", "Nb", "Se"),
            bond12=(0.5, 1.5),
            collect_angles=True,
        )
        np.testing.assert_array_equal(from_file.counts, direct.counts)
        np.testing.assert_allclose(from_file.angles, direct.angles, atol=1e-9)
        self.assertAlmostEqual(from_file.angles[0], 75.0, places=6)


class CliTests(unittest.TestCase):
    def test_end_to_end_csv_and_angles(self):
        positions = place([[5.0, 5.0, 5.0], [6.0, 5.0, 5.0], [5.0, 6.0, 5.0]])
        lines = [
            f"{index + 1} {element} [1] {p[0]:.12f} {p[1]:.12f} {p[2]:.12f} {index + 1} 0 0 0"
            for index, (element, p) in enumerate(zip(["Nb", "Se", "Se"], positions))
        ]
        with TemporaryDirectory() as scratch:
            scratch = Path(scratch)
            config = scratch / "tiny.rmc6f"
            write_rmc6f(config, lines)
            output = scratch / "out.csv"
            angles_path = scratch / "angles.txt"
            code = triplets_main(
                [
                    str(config),
                    "--triplet",
                    "Se",
                    "Nb",
                    "Se",
                    "--bond12",
                    "0.5",
                    "1.5",
                    "--output",
                    str(output),
                    "--dump-angles",
                    str(angles_path),
                ]
            )
            self.assertEqual(code, 0)
            content = output.read_text(encoding="utf-8")
            self.assertIn("angle_deg,counts,density_per_deg,sin_corrected", content)
            # The exact factor to RMCProfile's norm/sin(theta) at 1-degree bins.
            self.assertIn(
                "RMCProfile TRIPLETS norm/sin(theta) = sin_corrected * 0.008726535498", content
            )
            data_lines = [
                line for line in content.splitlines() if line and not line.startswith("#")
            ]
            self.assertEqual(len(data_lines), 181)  # header + 180 bins
            angles = [float(line) for line in angles_path.read_text().split()]
            self.assertEqual(len(angles), 1)
            self.assertAlmostEqual(angles[0], 90.0, places=5)

    def test_plot_and_overwrite_protection(self):
        # An empty window exercises write_plot's zero-histogram guard; the
        # second run must refuse to overwrite until --force is passed.
        positions = place([[5.0, 5.0, 5.0], [6.0, 5.0, 5.0]])
        lines = [
            f"{index + 1} {element} [1] {p[0]:.12f} {p[1]:.12f} {p[2]:.12f} {index + 1} 0 0 0"
            for index, (element, p) in enumerate(zip(["Nb", "Se"], positions))
        ]
        with TemporaryDirectory() as scratch:
            scratch = Path(scratch)
            config = scratch / "tiny.rmc6f"
            write_rmc6f(config, lines)
            output = scratch / "out.csv"
            plot = scratch / "out.png"
            argv = [
                str(config),
                "--triplet", "Se", "Nb", "Se",
                "--bond12", "3.0", "4.0",
                "--output", str(output),
                "--plot", str(plot),
            ]
            self.assertEqual(triplets_main(argv), 0)
            self.assertTrue(output.exists())
            self.assertTrue(plot.exists())
            self.assertEqual(triplets_main(argv), 1)
            self.assertEqual(triplets_main(argv + ["--force"]), 0)

    def test_csv_and_stdout_report_physical_bonds(self):
        import contextlib
        import io

        # Nb4 tetrahedron: 6 Nb-Nb bonds, found 12 times from the Nb centers.
        corners = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], float)
        positions = place(5.0 + corners * (3.0 / (2 * math.sqrt(2))))
        lines = [
            f"{index + 1} Nb [1] {p[0]:.12f} {p[1]:.12f} {p[2]:.12f} {index + 1} 0 0 0"
            for index, p in enumerate(positions)
        ]
        with TemporaryDirectory() as scratch:
            scratch = Path(scratch)
            config = scratch / "tetra.rmc6f"
            write_rmc6f(config, lines)
            output = scratch / "out.csv"
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                code = triplets_main(
                    [str(config), "--triplet", "Nb", "Nb", "Nb", "--bond12", "2.5", "3.5",
                     "--output", str(output)]
                )
            self.assertEqual(code, 0)
            header = output.read_text(encoding="utf-8")
        self.assertIn("# bonds in window12: 6 physical bonds", header)
        self.assertIn("12 bond vectors counted from the central atoms", header)
        self.assertIn("bonds 1-2:     6 ", stdout.getvalue())
        self.assertIn("3.00 per central atom", stdout.getvalue())

    def test_run_folder_picks_the_run_configuration(self):
        # triplets.parity.26: RMCProfile folders often hold the input
        # supercell next to the refined one; '.' sorts before '_', so the
        # first sorted .rmc6f was the start configuration, not the model the
        # app shows. The run's own outputs name the refined one.
        from rmc_toolkits.triplets_cli import resolve_config

        with TemporaryDirectory() as scratch:
            run = Path(scratch)
            for name in ("GaNb4Se8.rmc6f", "GaNb4Se8_5K.rmc6f", "GaNb4Se8_5K-00.log",
                         "GaNb4Se8_5K_PDFpartials.csv"):
                (run / name).write_text("", encoding="utf-8")
            self.assertEqual(resolve_config(run), run / "GaNb4Se8_5K.rmc6f")
            # No run outputs: the first sorted file, as before.
            bare = run / "bare"
            bare.mkdir()
            for name in ("b.rmc6f", "a.rmc6f"):
                (bare / name).write_text("", encoding="utf-8")
            self.assertEqual(resolve_config(bare), bare / "a.rmc6f")

    def test_documented_flags_exist(self):
        # triplets.physics.5 et al.: the algorithm reference named a flag
        # (--angles-out) the parser never had.
        import re

        from rmc_toolkits.triplets_cli import build_parser

        doc = (ROOT / "docs" / "algorithms" / "bond-geometry.md").read_text(encoding="utf-8")
        section = doc.split("### The `rmc-triplets` CLI", 1)[1].split("\n### ", 1)[0]
        flags = set(re.findall(r"(?<![\w-])--[a-z][a-z0-9-]*", section))
        self.assertIn("--dump-angles", flags)
        known = set(build_parser()._option_string_actions)
        self.assertEqual(flags - known, set())

    def test_missing_config_fails_cleanly(self):
        code = triplets_main(
            [
                "/nonexistent/path.rmc6f",
                "--triplet",
                "Se",
                "Nb",
                "Se",
                "--bond12",
                "1.0",
                "2.0",
            ]
        )
        self.assertEqual(code, 1)

    def test_truncated_header_fails_cleanly(self):
        with TemporaryDirectory() as scratch:
            config = Path(scratch) / "broken.rmc6f"
            config.write_text(
                "Supercell dimensions: 1 1 1\nLattice vectors (Ang):\n", encoding="utf-8"
            )
            code = triplets_main(
                [str(config), "--triplet", "Se", "Nb", "Se", "--bond12", "1.0", "2.0"]
            )
        self.assertEqual(code, 1)


@requires_sample
class SampleDataTests(unittest.TestCase):
    def test_se_nb_se_octahedra(self):
        result = bond_angles_from_rmc6f(
            SAMPLE_RMC6F, triplet=("Se", "Nb", "Se"), bond12=(2.2, 2.9)
        )
        self.assertEqual(result.apex_count, 16000)
        self.assertGreater(result.angle_count, 10000)
        self.assertEqual(int(result.counts.sum()), result.angle_count)
        self.assertAlmostEqual(
            float(result.density.sum()) * (180.0 / result.counts.size), 1.0, places=12
        )
        self.assertTrue(2.3 < result.mean_length12 < 2.8)


if __name__ == "__main__":
    unittest.main()
