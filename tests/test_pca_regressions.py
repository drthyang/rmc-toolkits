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
    displacement_cloud,
    load_site_displacements,
    pca_kde_volume,
    probability_scale,
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


class RotationInvariantKurtosisTests(unittest.TestCase):
    """pca.physics.8 / physics.40.

    nonGaussianity was the mean of three marginal kurtoses along the site's own
    PCA axes. For a (near-)isotropic site those axes are set by sampling noise
    and correlate with the fourth moments, so the headline number depended on
    an arbitrary frame. It is now Mardia's b2, normalised as (b2 - 15)/5 --
    affine invariant, and equal to the marginal excess kurtosis of any
    elliptical distribution -- and per-axis kappa carries a resolution flag.
    """

    @staticmethod
    def _scale_mixture(n, seed):
        # Spherical scale mixture of normals (elliptical): scale 1 (80%) or 2 (20%).
        # Marginal excess kurtosis along EVERY direction: 3*E[s^4]/E[s^2]^2 - 3.
        rng = np.random.default_rng(seed)
        scales = np.where(rng.random(n) < 0.8, 1.0, 2.0)
        return rng.normal(size=(n, 3)) * scales[:, None]

    def test_non_gaussianity_is_affine_invariant(self):
        cloud = self._scale_mixture(3000, 1) * 0.05
        transform = np.array([[1.0, 0.4, -0.3], [0.0, 0.8, 0.5], [0.2, 0.0, 0.6]])
        base = pca_kde_volume(cloud, grid=8, projections=False)["nonGaussianity"]
        for matrix in (transform, np.diag([1.0, 1.0 + 1e-7, 1.0 + 2e-7]), np.diag([1.0 + 2e-7, 1.0, 1.0 + 1e-7])):
            moved = pca_kde_volume(cloud @ matrix, grid=8, projections=False)["nonGaussianity"]
            self.assertAlmostEqual(moved, base, places=9)

    def test_non_gaussianity_is_the_elliptical_marginal_kurtosis(self):
        truth = 3.0 * (0.8 + 0.2 * 16.0) / (0.8 + 0.2 * 4.0) ** 2 - 3.0  # 1.6875
        cloud = self._scale_mixture(40000, 2) @ np.diag([0.12, 0.08, 0.05])
        result = pca_kde_volume(cloud, grid=8, projections=False, max_fit_points=40000)
        self.assertAlmostEqual(result["nonGaussianity"], truth, delta=0.12)

    def test_degenerate_axes_are_flagged_unresolved(self):
        rng = np.random.default_rng(4)
        isotropic = pca_kde_volume(rng.normal(size=(1000, 3)) * 0.1, grid=8, projections=False)
        self.assertEqual(isotropic["axisResolved"], [False, False, False])
        uniaxial = pca_kde_volume(rng.normal(size=(1000, 3)) * [0.2, 0.1, 0.1], grid=8, projections=False)
        self.assertEqual(uniaxial["axisResolved"], [True, False, False])
        triaxial = pca_kde_volume(rng.normal(size=(1000, 3)) * [0.2, 0.1, 0.05], grid=8, projections=False)
        self.assertEqual(triaxial["axisResolved"], [True, True, True])

    def test_site_table_carries_the_same_statistics(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "iso.rmc6f"
            supercell = (10, 10, 10)
            lattice = np.diag(np.asarray(supercell, dtype=float) * 8.0)
            lines = wrapped_site_lines((0.25, 0.25, 0.25), 0.08, supercell=supercell,
                                       cell_edge=8.0, seed=9)
            write_rmc6f(path, lines, supercell=supercell, lattice=lattice)
            sites = load_site_displacements(path)
        entry = site_ellipsoids(sites)[0]
        volume = pca_kde_volume(sites.displacements, grid=8, projections=False)
        self.assertEqual(entry["axisResolved"], [False, False, False])
        self.assertAlmostEqual(entry["nonGaussianity"], volume["nonGaussianity"], places=9)

    def test_symmetric_split_site_is_platykurtic(self):
        # A symmetric double well along x: excess kurtosis -2 d^4/(s^2+d^2)^2 < 0.
        rng = np.random.default_rng(6)
        n, d, s = 8000, 0.15, 0.08
        x = rng.choice([-d, d], size=n) + rng.normal(size=n) * s
        cloud = np.column_stack([x, rng.normal(size=n) * s, rng.normal(size=n) * s])
        result = pca_kde_volume(cloud, grid=8, projections=False)
        analytic = -2 * d**4 / (s**2 + d**2) ** 2
        self.assertTrue(result["axisResolved"][0])
        self.assertAlmostEqual(result["excessKurtosis"][0], analytic, delta=0.08)
        self.assertLess(result["nonGaussianity"], 0.0)


def mixed_site_file(path, *, supercell=(4, 4, 4), cell_edge=8.0, minority_every=4,
                    majority="Ga", minority="In", other_token="SE"):
    """Reference 1 is mixed (majority first in the file, minority LAST, the way
    RMCProfile groups atoms by type); reference 2 is a single species written
    with an upper-case token."""
    lattice = np.diag(np.asarray(supercell, dtype=float) * cell_edge)
    lines = wrapped_site_lines((0.25, 0.25, 0.25), 0.06, supercell=supercell,
                               cell_edge=cell_edge, seed=11, element=majority, reference=1)
    # Re-label every `minority_every`-th copy and move those lines to the end.
    major = [line for k, line in enumerate(lines) if k % minority_every]
    minor = [line.replace(f" {majority} ", f" {minority} ", 1)
             for k, line in enumerate(lines) if not k % minority_every]
    other = wrapped_site_lines((0.6, 0.1, 0.3), 0.07, supercell=supercell, cell_edge=cell_edge,
                               seed=12, element=other_token, reference=2)
    write_rmc6f(path, major + other + minor, supercell=supercell, lattice=lattice)
    return len(major), len(minor)


class MixedOccupancySiteTests(unittest.TestCase):
    """pca.parity.4 / numerics.17 / parity.20 / parity.26 / numerics.33 / physics.39.

    A reference number carried by several species (a solid solution, swap
    moves) was labelled by the LAST atom in Python and the FIRST in JS, the
    minority species vanished from the element list, and element pooling
    selected whole sites by that label.
    """

    def test_majority_label_with_composition(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "mixed.rmc6f"
            n_major, n_minor = mixed_site_file(path)
            sites = load_site_displacements(path)
        mixed, pure = site_ellipsoids(sites)
        self.assertEqual(mixed["element"], "Ga")          # majority, not the last atom
        self.assertTrue(mixed["mixed"])
        self.assertEqual(mixed["elementCounts"], {"Ga": n_major, "In": n_minor})
        self.assertEqual(pure["element"], "Se")           # 'SE' normalised
        self.assertFalse(pure["mixed"])
        self.assertEqual(pure["elementCounts"], {"Se": 64})
        self.assertEqual(sites.species, ["Ga", "In", "Se"])

    def test_element_pooling_selects_atoms_by_their_own_element(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "mixed.rmc6f"
            n_major, n_minor = mixed_site_file(path)
            sites = load_site_displacements(path)
        self.assertEqual(len(displacement_cloud(sites, element="In")), n_minor)
        self.assertEqual(len(displacement_cloud(sites, element="ga")), n_major)
        self.assertEqual(len(displacement_cloud(sites, element="Se")), 64)
        with self.assertRaisesRegex(ValueError, "Unknown element"):
            displacement_cloud(sites, element="Nb")

    def test_tie_goes_to_the_alphabetically_first_species(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "tie.rmc6f"
            mixed_site_file(path, minority_every=2, majority="Fe", minority="Co")
            sites = load_site_displacements(path)
        entry = site_ellipsoids(sites)[0]
        self.assertEqual(entry["elementCounts"], {"Co": 32, "Fe": 32})
        self.assertEqual(entry["element"], "Co")


class NonFiniteInputTests(unittest.TestCase):
    """pca.parity.6 / parity.23 / parity.27 / numerics.37.

    NaN/inf passed every ``<= 0`` guard: a NaN extent, bw_scale or bw returned
    an all-NaN volume (HTTP 200 with bare NaN tokens), and one NaN coordinate
    failed the whole batched eigensolve with 'Eigenvalues did not converge'.
    """

    def setUp(self):
        self.cloud = np.random.default_rng(3).normal(size=(300, 3)) * 0.1

    def test_non_finite_points_are_rejected(self):
        for bad in (np.nan, np.inf):
            cloud = self.cloud.copy()
            cloud[17, 1] = bad
            with self.assertRaisesRegex(ValueError, "non-finite"):
                pca_kde_volume(cloud, grid=8)

    def test_non_finite_parameters_are_rejected(self):
        for kwargs in ({"extent": np.nan}, {"extent": np.inf}, {"bw_scale": np.nan},
                       {"bw_scale": np.inf}, {"bw": np.nan}, {"bw": np.inf}, {"grid": np.nan}):
            with self.subTest(**{k: str(v) for k, v in kwargs.items()}):
                with self.assertRaises(ValueError):
                    pca_kde_volume(self.cloud, **{"grid": 8, **kwargs})

    def test_a_nan_coordinate_is_skipped_and_the_warning_names_the_line(self):
        # Integrated 1.0 rule (parsers group, both runtimes): a .rmc6f atom line
        # with a non-finite coordinate is skipped and counted, never silently --
        # SiteDisplacements.parse_warning names the first such line (the
        # browser worker's parseWarning carries the same text).
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "nan.rmc6f"
            supercell = (4, 4, 4)
            lines = wrapped_site_lines((0.25, 0.25, 0.25), 0.06, supercell=supercell,
                                       cell_edge=8.0, seed=1)
            clean = load_site_displacements(_written(Path(tmp) / "clean.rmc6f", lines, supercell))
            self.assertIsNone(clean.parse_warning)

            parts = lines[9].split()
            parts[3] = "NaN"
            lines[9] = " ".join(parts)
            sites = load_site_displacements(_written(path, lines, supercell))
            self.assertEqual(int(sites.counts.sum()), 64 - 1)
            self.assertTrue(np.all(np.isfinite(sites.displacements)))
            self.assertEqual(
                sites.parse_warning,
                f"1 atom lines skipped for non-finite coordinates (first: '{lines[9]}')",
            )
            self.assertIn("'10 Se [1] NaN ", sites.parse_warning)
            # The skipped atom does not move the site: its mean is that of the
            # other 63 copies, so the ADP differs from the clean file only by
            # one atom's worth.
            self.assertLess(
                abs(site_ellipsoids(sites)[0]["uIso"] - site_ellipsoids(clean)[0]["uIso"])
                / site_ellipsoids(clean)[0]["uIso"],
                0.1,
            )

    def test_no_parseable_atom_is_an_error_that_says_why(self):
        with TemporaryDirectory() as tmp:
            supercell = (2, 2, 2)
            lines = wrapped_site_lines((0.25, 0.25, 0.25), 0.06, supercell=supercell,
                                       cell_edge=8.0, seed=1)
            lines = [" ".join(line.split()[:3] + ["inf"] + line.split()[4:]) for line in lines]
            with self.assertRaisesRegex(
                ValueError, "no atoms could be parsed — 8 atom lines skipped for non-finite"
            ):
                load_site_displacements(_written(Path(tmp) / "all_nan.rmc6f", lines, supercell))


def _written(path, lines, supercell):
    write_rmc6f(path, lines, supercell=supercell,
                lattice=np.diag(np.asarray(supercell, dtype=float) * 8.0))
    return path


class CubicDisplayBoxTests(unittest.TestCase):
    """pca.numerics.14 / numerics.32 / physics.44.

    With cubic_box every axis was sampled over PC1's half-width on the same G
    nodes, so a thin PC3 kernel fell between nodes: captured mass and the
    enclosed-mass iso levels collapsed (0% for a planar cloud on an even grid,
    ~228000% on an odd one). The cube is now a display box only; the volume is
    always sampled on the per-axis box, where the separable estimator is exact.
    """

    def test_planar_cloud_keeps_unit_mass_on_even_and_odd_grids(self):
        rng = np.random.default_rng(8)
        cloud = np.column_stack([rng.normal(size=1000) * 0.1, rng.normal(size=1000) * 0.07,
                                 np.zeros(1000)])
        for grid in (40, 41):
            with self.subTest(grid=grid):
                result = pca_kde_volume(cloud, grid=grid, extent=4.0, cubic_box=True, projections=False)
                self.assertTrue(result["degenerate"])
                self.assertGreater(result["mass"], 0.99)
                self.assertLess(result["mass"], 1.01)
                self.assertEqual(len(result["massLevels"]), 101)

    def test_thin_cloud_volume_does_not_depend_on_the_display_box(self):
        rng = np.random.default_rng(9)
        cloud = rng.normal(size=(1000, 3)) * np.array([0.3, 0.1, 0.015])  # sigma1/sigma3 = 20
        cubic = pca_kde_volume(cloud, grid=40, extent=4.0, cubic_box=True, projections=False)
        plain = pca_kde_volume(cloud, grid=40, extent=4.0, cubic_box=False, projections=False)
        self.assertGreater(cubic["mass"], 0.99)
        np.testing.assert_allclose(cubic["density"], plain["density"], rtol=0, atol=0)
        self.assertEqual(cubic["massLevels"], plain["massLevels"])
        # The display box is the cube; the sampled box stays per axis.
        self.assertEqual(len(set(cubic["boxHalfWidths"])), 1)
        self.assertAlmostEqual(cubic["boxHalfWidths"][0], max(cubic["halfWidths"]))
        self.assertEqual(plain["boxHalfWidths"], plain["halfWidths"])


class KernelBroadeningGuidanceTests(unittest.TestCase):
    """pca.physics.12 / physics.41 -- the physics the page's guidance now states.

    At EQUAL levels the p% mass isosurface of a Gaussian cloud sits sqrt(1+f^2)
    outside the p% ellipsoid (kernel broadening), so only a surface inside the
    ellipsoid signals anharmonicity; at the old defaults (25% surface, 50%
    ellipsoid) a Gaussian site sat ~16% inside by construction.
    """

    def _surface_radius_along_pc1(self, result, level):
        grid = result["grid"]
        centre = grid // 2
        density = np.asarray(result["density"]).reshape(grid, grid, grid)[:, centre, centre]
        coords = np.asarray(result["axisCoords"][0])
        outer = np.arange(centre, grid - 1)
        crossing = outer[(density[outer] >= level) & (density[outer + 1] < level)][0]
        t = (density[crossing] - level) / (density[crossing] - density[crossing + 1])
        return coords[crossing] + t * (coords[crossing + 1] - coords[crossing])

    def test_equal_level_gaussian_surface_sits_sqrt_one_plus_f2_outside(self):
        rng = np.random.default_rng(12)
        cloud = rng.normal(size=(20000, 3)) * np.array([0.12, 0.10, 0.08])
        # A wide fixed kernel (f = 0.6, broadening 1.166) makes the offset unmistakable.
        result = pca_kde_volume(cloud, bw=0.6, grid=41, extent=4.0, projections=False,
                                max_fit_points=20000)
        broadening = np.sqrt(1.0 + result["factor"] ** 2)
        levels = {round(entry["p"], 2): entry["level"] for entry in result["massLevels"]}
        sigma1 = result["rms"][0]
        at_50 = self._surface_radius_along_pc1(result, levels[0.5]) / (probability_scale(0.5) * sigma1)
        self.assertAlmostEqual(at_50, broadening, delta=0.03)
        # The pre-1.0 default pairing: a 25% surface against the 50% ellipsoid.
        at_25 = self._surface_radius_along_pc1(result, levels[0.25]) / (probability_scale(0.5) * sigma1)
        self.assertLess(at_25, 0.9)


if __name__ == "__main__":
    unittest.main()
