# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regression tests for the PCA-ellipsoid engine (rmc_toolkits.pca_kde).

Each class pins one defect found by the 1.0 audit; the class docstring names it.
Fixtures are synthetic ``.rmc6f`` files written the way RMCProfile writes them:
coordinates are supercell fractions wrapped into [0, 1), and the cell indices
name the box copy the atom belongs to.
"""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from rmc_toolkits.pca_kde import (
    load_site_displacements,
    pca_kde_volume,
    site_ellipsoids,
    site_pca_kde,
)


ROOT = Path(__file__).resolve().parents[1]
AVERAGE_RMC6F = ROOT / "data" / "5K_try1" / "GaNb4Se8_5KAVERAGE.rmc6f"


def write_rmc6f(path: Path, atom_lines, *, supercell, lattice):
    header = [
        f"Supercell dimensions {supercell[0]} {supercell[1]} {supercell[2]}",
        "Lattice vectors (Ang):",
        *(" ".join(repr(float(v)) for v in row) for row in lattice),
        "Atoms:",
    ]
    path.write_text("\n".join(header + list(atom_lines)) + "\n", encoding="utf-8")


def wrapped_site_lines(site, sigma_angstrom, *, supercell, cell_edge, seed,
                       element="Se", reference=1, start=1):
    """One site at unit-cell fraction ``site`` with isotropic Gaussian spread.

    Every box copy (ix, iy, iz) gets one atom; its stored coordinate is the
    supercell fraction wrapped into [0, 1) exactly as RMCProfile writes it.
    """
    rng = np.random.default_rng(seed)
    supercell = np.asarray(supercell, dtype=float)
    lines = []
    atom = start
    for ix in range(int(supercell[0])):
        for iy in range(int(supercell[1])):
            for iz in range(int(supercell[2])):
                cell = np.array([ix, iy, iz], dtype=float)
                delta = rng.normal(0.0, sigma_angstrom, size=3) / cell_edge  # unit-cell fraction
                coord = np.mod((cell + np.asarray(site) + delta) / supercell, 1.0)
                lines.append(
                    f"{atom} {element} [1] {coord[0]:.12f} {coord[1]:.12f} {coord[2]:.12f} "
                    f"{reference} {ix} {iy} {iz}"
                )
                atom += 1
    return lines


class SiteCentredFoldTests(unittest.TestCase):
    """pca.parity.1 / physics.9 / numerics.13 / parity.25 / numerics.31 / physics.38.

    The half-box fold ``o -= round(o)`` is taken about 0, so a site whose own
    position s/N_i sits near half a supercell period -- s ~ 1/2 with N_i = 1, or
    s ~ 1 with N_i = 2 -- is torn into two lumps a whole supercell edge apart.
    """

    def _uiso_and_site(self, supercell, site, *, cell_edge=10.0, sigma=0.08, seed=3):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "fold.rmc6f"
            lattice = np.diag(np.asarray(supercell, dtype=float) * cell_edge)
            lines = wrapped_site_lines(site, sigma, supercell=supercell,
                                       cell_edge=cell_edge, seed=seed)
            write_rmc6f(path, lines, supercell=supercell, lattice=lattice)
            sites = load_site_displacements(path)
            entry = site_ellipsoids(sites)[0]
            return entry, sites

    def assert_intact(self, supercell, site, sigma=0.08):
        entry, sites = self._uiso_and_site(supercell, site, sigma=sigma)
        # The true U is sigma^2 on every axis; a torn cloud is ~10^3-10^4 larger.
        np.testing.assert_allclose(np.asarray(entry["rms"]), sigma, rtol=0.3)
        self.assertLess(np.abs(sites.displacements).max(), 8 * sigma)
        # The marker sits at the site, not between the two halves (circular distance).
        delta = np.asarray(entry["siteFractional"]) - np.asarray(site)
        delta -= np.round(delta)
        self.assertLess(np.abs(delta).max(), 0.01)

    def test_single_cell_axis_site_at_one_half(self):
        self.assert_intact((6, 6, 1), (0.3, 0.3, 0.5))

    def test_two_cell_axis_site_near_one(self):
        self.assert_intact((6, 6, 2), (0.3, 0.3, 0.99))

    def test_single_cell_axis_with_a_large_spread(self):
        # A large cell refined in a 1xNxN box, site on the special position x = 1/2.
        self.assert_intact((1, 8, 8), (0.5, 0.5, 0.5), sigma=0.2)

    def test_regular_boxes_are_unchanged(self):
        # N_i >= 3 never reached the old fold threshold; the fix must not move them.
        self.assert_intact((6, 6, 6), (0.5, 0.99, 0.0))


class ZeroSpreadTests(unittest.TestCase):
    """pca.parity.3 / parity.24 / parity.30.

    In an *AVERAGE.rmc6f (or an ideal start) every copy of a site sits at the
    same offset; the covariance is ~1e-28 A^2 of round-off. Nothing caught it:
    anisotropy, kurtosis and axes were noise that differed between engines, and
    the KDE drew a femtometre 'cloud'.
    """

    def _write_frozen(self, path, supercell=(6, 6, 6), cell_edge=10.0):
        lattice = np.diag(np.asarray(supercell, dtype=float) * cell_edge)
        frozen = wrapped_site_lines((0.25, 0.5, 0.75), 0.0, supercell=supercell,
                                    cell_edge=cell_edge, seed=1, element="Ga", reference=1)
        moving = wrapped_site_lines((0.6, 0.1, 0.3), 0.07, supercell=supercell,
                                    cell_edge=cell_edge, seed=2, element="Se", reference=2,
                                    start=len(frozen) + 1)
        write_rmc6f(path, frozen + moving, supercell=supercell, lattice=lattice)

    def test_frozen_site_is_flagged_and_its_noise_suppressed(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "frozen.rmc6f"
            self._write_frozen(path)
            sites = load_site_displacements(path)
            frozen, moving = site_ellipsoids(sites)
        self.assertTrue(frozen["zeroSpread"])
        self.assertTrue(frozen["degenerate"])
        self.assertIsNone(frozen["anisotropy"])
        self.assertIsNone(frozen["nonGaussianity"])
        self.assertEqual(frozen["excessKurtosis"], [None, None, None])
        self.assertIsNone(frozen["axes"])
        self.assertLess(frozen["uIso"], 1e-12)
        # The moving site next to it is untouched.
        self.assertFalse(moving["zeroSpread"])
        self.assertFalse(moving["degenerate"])
        self.assertIsNotNone(moving["axes"])
        self.assertTrue(all(np.isfinite(moving["excessKurtosis"])))

    def test_kde_refuses_a_zero_spread_cloud(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "frozen.rmc6f"
            self._write_frozen(path)
            sites = load_site_displacements(path)
        with self.assertRaisesRegex(ValueError, "zero spread"):
            site_pca_kde(sites, reference_number=1, grid=12, projections=False)
        # A cloud of pure round-off (the femtometre case) is refused too.
        rng = np.random.default_rng(0)
        with self.assertRaisesRegex(ValueError, "zero spread"):
            pca_kde_volume(rng.normal(size=(200, 3)) * 1e-14, grid=8)

    def test_collapsed_axis_has_no_kurtosis(self):
        # A planar cloud: z identically 0. kappa along PC3 is 0/0, not a number.
        rng = np.random.default_rng(5)
        cloud = np.column_stack([rng.normal(size=2000) * 0.1, rng.normal(size=2000) * 0.07,
                                 np.zeros(2000)])
        result = pca_kde_volume(cloud, grid=16, projections=False)
        self.assertTrue(result["degenerate"])
        self.assertIsNone(result["excessKurtosis"][2])
        self.assertTrue(np.isfinite(result["excessKurtosis"][0]))
        self.assertTrue(np.isfinite(result["nonGaussianity"]))

    @unittest.skipUnless(AVERAGE_RMC6F.exists(), "GaNb4Se8 AVERAGE sample not present in data/")
    def test_real_average_configuration_is_all_zero_spread(self):
        sites = load_site_displacements(AVERAGE_RMC6F)
        entries = site_ellipsoids(sites)
        self.assertEqual(len(entries), 52)
        self.assertTrue(all(entry["zeroSpread"] and entry["degenerate"] for entry in entries))
        self.assertTrue(all(entry["nonGaussianity"] is None for entry in entries))


if __name__ == "__main__":
    unittest.main()
