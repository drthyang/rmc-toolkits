# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Atoms exactly on a slab face are in the slab, in every runtime.

An ideal (unrelaxed) configuration puts whole sites on coordinates such as
1/8, which slider positions (z_c on a 0.001 grid, dz on a 0.01 grid) hit
exactly. Both runtimes and the Slab-In-Cell highlight share one inclusive test,
``|d - z_c| <= dz/2 + SLAB_FACE_TOLERANCE`` on the normalised depth; the
worker's twin of this file is ``workers/__tests__/slabFaces.test.js``.
"""

from __future__ import annotations

import unittest

import numpy as np

from rmc_toolkits.kde import SLAB_FACE_TOLERANCE, oriented_kde_slice

C_SLICE = {
    "normal": np.array([0.0, 0.0, 1.0]),
    "u_axis": np.array([1.0, 0.0, 0.0]),
    "v_axis": np.array([0.0, 1.0, 0.0]),
}


def _ideal_layers() -> np.ndarray:
    """16 atoms on each of the layers z = k/8 (exact binary fractions)."""
    grid = [(0.1 + 0.2 * i, 0.15 + 0.2 * j) for i in range(4) for j in range(4)]
    return np.array([(x, y, k / 8) for k in range(8) for (x, y) in grid])


def _exact_count(center_milli: int, thickness_milli: int) -> int:
    """Atoms in the slab by exact integer arithmetic in thousandths (with wrap)."""
    count = 0
    for k in range(8):
        depth = 125 * k
        if any(2 * abs(depth + shift * 1000 - center_milli) <= thickness_milli for shift in (-1, 0, 1)):
            count += 16
    return count


class SlabFaceTests(unittest.TestCase):
    def test_face_atoms_of_an_ideal_configuration_are_included(self):
        points = _ideal_layers()
        for thickness_milli in (80, 100, 150, 250):
            for center_milli in range(0, 1001, 5):
                with self.subTest(center=center_milli / 1000, thickness=thickness_milli / 1000):
                    result = oriented_kde_slice(
                        points,
                        center=center_milli / 1000,
                        thickness=thickness_milli / 1000,
                        bw=0.05,
                        grid=16,
                        n_levels=0,
                        **C_SLICE,
                    )
                    self.assertEqual(result["slabCount"], _exact_count(center_milli, thickness_milli))

    def test_face_tolerance_admits_round_off_but_not_more(self):
        rng = np.random.default_rng(0)
        xy = rng.random((10, 2))
        inside = np.column_stack([xy, np.full(10, 0.54 + 0.5 * SLAB_FACE_TOLERANCE)])
        outside = np.column_stack([xy, np.full(10, 0.54 + 2.0 * SLAB_FACE_TOLERANCE)])
        for layer, expected in ((inside, 10), (outside, 0)):
            result = oriented_kde_slice(layer, center=0.5, thickness=0.08, bw=0.05, grid=16, n_levels=0, **C_SLICE)
            self.assertEqual(result["slabCount"], expected)

    def test_payload_z_and_dz_echo_the_slider_fractions(self):
        rng = np.random.default_rng(1)
        result = oriented_kde_slice(
            rng.random((200, 3)), center=0.5, thickness=0.08, normal=np.array([1.0, -1.0, 0.0]), grid=16, n_levels=0
        )
        self.assertEqual(result["z"], 0.5)
        self.assertEqual(result["dz"], 0.08)
        self.assertAlmostEqual(result["depth"], 0.0, places=12)  # d_min + 0.5 * sqrt(2), d_min = -1/sqrt(2)
        self.assertAlmostEqual(result["depthThickness"], 0.08 * np.sqrt(2.0), places=12)


if __name__ == "__main__":
    unittest.main()
