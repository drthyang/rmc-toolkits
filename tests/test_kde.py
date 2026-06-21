from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from rmc_toolkits.kde import MAX_KDE_FIT_POINTS, kde_slice, load_unit_cell_positions, oriented_kde_slice


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

# The GNSe example dataset is gitignored (see README), so sample-backed tests
# skip rather than fail when it is not present locally / in CI.
requires_sample = unittest.skipUnless(
    (DATA / "GNSe.rmc6f").exists(),
    "GNSe sample data not present in data/ (gitignored)",
)


class KdeTests(unittest.TestCase):
    @requires_sample
    def test_load_unit_cell_positions_filters_by_element(self):
        all_positions = load_unit_cell_positions(DATA / "GNSe.rmc6f")
        ga_positions = load_unit_cell_positions(DATA / "GNSe.rmc6f", element="Ga")

        self.assertEqual(all_positions.positions.shape, (52000, 3))
        self.assertEqual(ga_positions.positions.shape, (4000, 3))
        np.testing.assert_allclose(all_positions.cell_lengths, [10.4116, 10.4116, 10.4116])
        self.assertTrue(np.all(ga_positions.positions >= 0.0))
        self.assertTrue(np.all(ga_positions.positions <= 10.4116))

    def test_load_unit_cell_positions_preserves_nonorthogonal_basis(self):
        with TemporaryDirectory() as tmpdir:
            rmc6f = Path(tmpdir) / "skewed.rmc6f"
            rmc6f.write_text(
                "\n".join(
                    [
                        "Supercell dimensions 1 1 1",
                        "Lattice vectors",
                        "2.0 0.0 0.0",
                        "1.0 3.0 0.0",
                        "0.5 0.25 4.0",
                        "Atoms:",
                        "1 Ga Ga 0.25 0.50 0.75 1 0 0 0",
                    ]
                ),
                encoding="utf-8",
            )

            positions = load_unit_cell_positions(rmc6f)

        np.testing.assert_allclose(
            positions.unit_vectors,
            [[2.0, 0.0, 0.0], [1.0, 3.0, 0.0], [0.5, 0.25, 4.0]],
        )
        np.testing.assert_allclose(positions.fractional_positions, [[0.25, 0.5, 0.75]])
        np.testing.assert_allclose(positions.positions, [[1.375, 1.6875, 3.0]])

    def test_kde_slice_returns_density_grid_and_contours(self):
        positions = np.array(
            [
                [0.2, 0.2, 0.48],
                [0.3, 0.4, 0.49],
                [0.45, 0.2, 0.5],
                [0.55, 0.7, 0.51],
                [0.7, 0.55, 0.52],
                [0.8, 0.8, 0.5],
                [0.1, 0.9, 0.9],
            ]
        )

        result = kde_slice(
            positions,
            z_center=0.5,
            dz=0.08,
            xlim=(0.0, 1.0),
            ylim=(0.0, 1.0),
            bw=0.2,
            grid=32,
            n_levels=3,
        )

        self.assertEqual(result["grid"], 32)
        self.assertEqual(result["slabCount"], 6)
        self.assertEqual(result["fitCount"], 6)
        self.assertEqual(len(result["density"]), 32)
        self.assertEqual(len(result["density"][0]), 32)
        self.assertGreater(result["vmax"], result["vmin"])
        self.assertLessEqual(len(result["contours"]), 3)

    def test_kde_slice_clamps_grid_and_empty_slab(self):
        result = kde_slice(
            np.empty((0, 3)),
            z_center=0.5,
            dz=0.1,
            xlim=(0.0, 1.0),
            ylim=(0.0, 1.0),
            grid=2,
        )

        self.assertEqual(result["grid"], 16)
        self.assertEqual(result["slabCount"], 0)
        self.assertEqual(result["fitCount"], 0)
        self.assertEqual(result["vmin"], 0.0)
        self.assertEqual(result["vmax"], 0.0)

    def test_kde_slice_handles_nonempty_structure_with_empty_slab(self):
        result = kde_slice(
            np.array([[0.1, 0.2, 0.1], [0.3, 0.4, 0.2]]),
            z_center=0.9,
            dz=0.05,
            xlim=(0.0, 1.0),
            ylim=(0.0, 1.0),
            grid=16,
        )

        self.assertEqual(result["slabCount"], 0)
        self.assertEqual(result["fitCount"], 0)
        self.assertEqual(result["vmin"], 0.0)
        self.assertEqual(result["vmax"], 0.0)

    def test_kde_slice_handles_degenerate_slab_without_error(self):
        positions = np.array(
            [
                [0.1, 0.1, 0.5],
                [0.2, 0.2, 0.5],
                [0.3, 0.3, 0.5],
                [0.4, 0.4, 0.5],
                [0.5, 0.5, 0.5],
            ]
        )

        result = kde_slice(
            positions,
            z_center=0.5,
            dz=0.1,
            xlim=(0.0, 1.0),
            ylim=(0.0, 1.0),
            grid=16,
        )

        self.assertEqual(result["slabCount"], 5)
        self.assertEqual(result["fitCount"], 0)
        self.assertEqual(result["vmin"], 0.0)
        self.assertEqual(result["vmax"], 0.0)

    def test_kde_slice_reports_subsampled_fit_count(self):
        rng = np.random.default_rng(0)
        positions = np.column_stack(
            [
                rng.random(MAX_KDE_FIT_POINTS + 10),
                rng.random(MAX_KDE_FIT_POINTS + 10),
                np.full(MAX_KDE_FIT_POINTS + 10, 0.5),
            ]
        )

        result = kde_slice(
            positions,
            z_center=0.5,
            dz=0.1,
            xlim=(0.0, 1.0),
            ylim=(0.0, 1.0),
            bw=0.2,
            grid=16,
            n_levels=0,
        )

        self.assertEqual(result["slabCount"], MAX_KDE_FIT_POINTS + 10)
        self.assertEqual(result["fitCount"], MAX_KDE_FIT_POINTS)

    def test_oriented_kde_slice_supports_axis_and_custom_normals(self):
        positions = np.array(
            [
                [0.48, 0.2, 0.2],
                [0.49, 0.3, 0.4],
                [0.5, 0.45, 0.2],
                [0.51, 0.55, 0.7],
                [0.52, 0.7, 0.55],
                [0.5, 0.8, 0.8],
                [0.9, 0.1, 0.9],
            ]
        )

        axis_result = oriented_kde_slice(
            positions,
            center=0.5,
            thickness=0.08,
            normal=np.array([1.0, 0.0, 0.0]),
            u_axis=np.array([0.0, 1.0, 0.0]),
            v_axis=np.array([0.0, 0.0, 1.0]),
            bw=0.2,
            grid=24,
            n_levels=2,
        )
        custom_result = oriented_kde_slice(
            positions,
            center=0.5,
            thickness=0.2,
            normal=np.array([1.0, 1.0, 0.0]),
            bw=0.2,
            grid=24,
            n_levels=2,
        )

        self.assertEqual(axis_result["slabCount"], 6)
        self.assertEqual(axis_result["normal"], [1.0, 0.0, 0.0])
        self.assertEqual(axis_result["grid"], 24)
        self.assertTrue(axis_result["planePolygon"])
        self.assertGreaterEqual(custom_result["slabCount"], 1)
        self.assertAlmostEqual(np.linalg.norm(custom_result["normal"]), 1.0)
        self.assertTrue(custom_result["planeVertices"])


if __name__ == "__main__":
    unittest.main()
