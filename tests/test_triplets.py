# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

import math
import unittest
from itertools import combinations, product
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from rmc_toolkits.triplets import bond_angle_distribution, bond_angles_from_rmc6f
from rmc_toolkits.triplets_cli import main as triplets_main


ROOT = Path(__file__).resolve().parents[1]
SAMPLE_RMC6F = ROOT / "data" / "5K_try1" / "GaNb4Se8_5K.rmc6f"

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
        if end1 == end2 and tuple(window12) == tuple(window23):
            pairs = [(bonds1[i], bonds1[j]) for i, j in combinations(range(len(bonds1)), 2)]
        else:
            pairs = [
                (one, two)
                for one in bonds1
                for two in bonds2
                if not (one[0] == two[0] and one[1] == two[1])
            ]
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

    def test_overlapping_windows_use_ordered_counting(self):
        # Both neighbours fall in both windows: the A-B-C convention counts the
        # (1->2, 2->3) assignment both ways, so the same geometric angle
        # appears twice. This documents the ordered-counting convention.
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
        self.assertEqual(result.angle_count, 2)
        np.testing.assert_allclose(result.angles, [120.0, 120.0], atol=1e-9)


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


class HistogramConventionTests(unittest.TestCase):
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
