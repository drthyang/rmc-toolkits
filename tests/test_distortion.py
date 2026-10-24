# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Known answers for the polyhedral-distortion engine (rmc_toolkits/distortion.py).

Baur's DI, Robinson's quadratic elongation <lambda> and bond-angle variance
sigma^2_theta on polyhedra whose values follow in closed form; the average
structure of configurations whose average is known; the thermal expectation
of DI; and the refusals. The browser port is held to the same engine by the
parity fixture (tests/test_distortion_fixture.py).
"""

from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    import distortion_configs as configs
except ModuleNotFoundError:  # `python -m unittest tests.test_distortion`
    from tests import distortion_configs as configs

from rmc_toolkits import polyhedral_distortion as exported, polyhedral_distortion_from_file
from rmc_toolkits.distortion import (
    APP_MAX_BONDS,
    OCTAHEDRON_MATCHINGS,
    PAIRS,
    TETRAHEDRAL_ANGLE_DEG,
    distortion_configuration,
    polyhedral_distortion,
    polyhedron_indices,
    site_species,
    summarize,
)
from rmc_toolkits.triplets import APP_MAX_ANGLES


def run(built, centre="Cu", ligand="F", window=(1.0, 3.5), **options):
    atoms, lattice, supercell = built
    return polyhedral_distortion(distortion_configuration(atoms, lattice, supercell), centre, ligand, *window, **options)


def regular_octahedron_values(result, block="configuration"):
    stats = result[block]
    return (
        stats["distortionIndex"]["mean"],
        stats["quadraticElongation"]["mean"],
        stats["angleVariance"]["mean"],
    )


class ConstantsTests(unittest.TestCase):
    def test_the_octahedron_matchings_are_the_15_perfect_matchings(self):
        self.assertEqual(len(OCTAHEDRON_MATCHINGS), 15)
        seen = set()
        for matching in OCTAHEDRON_MATCHINGS:
            vertices = sorted(v for k in matching for v in PAIRS[6][k])
            self.assertEqual(vertices, list(range(6)))
            seen.add(frozenset(matching))
        self.assertEqual(len(seen), 15)
        # The tie order: the first vertex pairs with 1, then 2, ...
        self.assertEqual([PAIRS[6][k] for k in OCTAHEDRON_MATCHINGS[0]], [(0, 1), (2, 3), (4, 5)])
        self.assertEqual([PAIRS[6][k] for k in OCTAHEDRON_MATCHINGS[-1]], [(0, 5), (1, 4), (2, 3)])

    def test_the_tetrahedral_angle_is_exact(self):
        self.assertAlmostEqual(TETRAHEDRAL_ANGLE_DEG, 109.47122063449069, places=12)
        self.assertEqual(round(TETRAHEDRAL_ANGLE_DEG, 2), 109.47)

    def test_summaries_interpolate_between_order_statistics(self):
        summary = summarize([4.0, 1.0, 3.0, 2.0, 5.0])
        self.assertEqual(summary["count"], 5)
        self.assertEqual(summary["mean"], 3.0)
        self.assertAlmostEqual(summary["sd"], math.sqrt(2.5), places=14)
        self.assertAlmostEqual(summary["p10"], 1.4, places=14)
        self.assertEqual(summary["median"], 3.0)
        self.assertAlmostEqual(summary["p90"], 4.6, places=14)
        self.assertIsNone(summarize([]))
        self.assertIsNone(summarize([1.0])["sd"])
        self.assertEqual(summarize([np.nan, 2.0])["count"], 1)

    def test_a_site_is_its_majority_species_ties_to_the_code_point_first(self):
        self.assertEqual(site_species({"Ti": 7, "Zr": 3}), "Ti")
        self.assertEqual(site_species({"Zr": 3, "Ti": 7}), "Ti")
        self.assertEqual(site_species({"Zr": 5, "Ti": 5}), "Ti")
        self.assertEqual(site_species({"b": 1, "B": 1}), "B")

    def test_the_package_exports_the_engine(self):
        self.assertIs(exported, polyhedral_distortion)
        self.assertLess(APP_MAX_BONDS, APP_MAX_ANGLES)


class RegularPolyhedraTests(unittest.TestCase):
    def test_a_regular_octahedron_is_0_1_0(self):
        result = run(configs.isolated_polyhedron(configs.octahedron()))
        self.assertEqual(result["polyhedron"], {"coordination": 6, "shape": "octahedron", "idealAngle": 90.0, "cisAngles": 12})
        self.assertEqual(result["configuration"]["polyhedra"], 8)
        di, elongation, variance = regular_octahedron_values(result)
        self.assertLess(abs(di), 1e-14)
        self.assertAlmostEqual(elongation, 1.0, places=13)
        self.assertLess(variance, 1e-20)
        self.assertAlmostEqual(result["configuration"]["transAngle"]["mean"], 180.0, places=10)
        # The average structure is the same octahedron.
        self.assertEqual(result["average"]["polyhedra"], 1)
        self.assertIsNone(result["average"]["reason"])
        di, elongation, variance = regular_octahedron_values(result, "average")
        self.assertLess(abs(di), 1e-14)
        self.assertAlmostEqual(elongation, 1.0, places=13)
        self.assertLess(variance, 1e-20)

    def test_a_regular_tetrahedron_is_0_1_0(self):
        result = run(configs.isolated_polyhedron(configs.TETRAHEDRON * 2.0))
        self.assertEqual(result["polyhedron"]["shape"], "tetrahedron")
        self.assertEqual(result["polyhedron"]["cisAngles"], 6)
        self.assertAlmostEqual(result["polyhedron"]["idealAngle"], 109.47122063449069, places=12)
        di, elongation, variance = regular_octahedron_values(result)
        self.assertLess(abs(di), 1e-14)
        self.assertAlmostEqual(elongation, 1.0, places=13)
        self.assertLess(variance, 1e-20)
        self.assertIsNone(result["configuration"]["transAngle"])

    def test_a_rotated_and_sheared_cell_changes_nothing(self):
        """The indices are geometric: a rigid rotation of the whole configuration
        (a triclinic-looking lattice) gives the same values."""
        atoms, lattice, supercell = configs.isolated_polyhedron(configs.trigonal_octahedron(48.0, 2.1))
        theta = 0.7
        rotation = np.array([[math.cos(theta), -math.sin(theta), 0], [math.sin(theta), math.cos(theta), 0], [0, 0, 1]])
        tilt = np.array([[1, 0, 0], [0, math.cos(0.4), -math.sin(0.4)], [0, math.sin(0.4), math.cos(0.4)]])
        turned = lattice @ (tilt @ rotation).T
        first = run((atoms, lattice, supercell))
        second = run((atoms, turned, supercell))
        for key in ("distortionIndex", "quadraticElongation", "angleVariance", "transAngle"):
            self.assertAlmostEqual(first["configuration"][key]["mean"], second["configuration"][key]["mean"], places=10)
            self.assertAlmostEqual(first["average"][key]["mean"], second["average"][key]["mean"], places=10)


class ClosedFormTests(unittest.TestCase):
    def test_a_jahn_teller_4_plus_2_octahedron(self):
        """4 x 1.93 + 2 x 2.25 A at right angles: DI = (4|1.93 - <l>| + 2|2.25 - <l>|) / (6 <l>),
        <l> = 12.22 / 6 = 2.03667 (DI = 0.0698); <lambda> = (4 1.93^2 + 2 2.25^2) / (6 l_0^2) with
        l_0^3 = 1.93^2 2.25 (V = (4/3) 1.93^2 2.25); sigma^2 = 0."""
        result = run(configs.isolated_polyhedron(configs.octahedron((1.93, 1.93, 1.93, 1.93, 2.25, 2.25))))
        mean = (4 * 1.93 + 2 * 2.25) / 6
        self.assertAlmostEqual(mean, 2.0366666666666666, places=14)
        di = (4 * abs(1.93 - mean) + 2 * abs(2.25 - mean)) / (6 * mean)
        self.assertAlmostEqual(di, 0.06983087834151, places=12)
        elongation = (4 * 1.93**2 + 2 * 2.25**2) / (6 * (1.93**2 * 2.25) ** (2.0 / 3.0))
        stats = result["configuration"]
        self.assertAlmostEqual(stats["distortionIndex"]["mean"], di, places=12)
        self.assertAlmostEqual(stats["quadraticElongation"]["mean"], elongation, places=12)
        self.assertAlmostEqual(stats["quadraticElongation"]["mean"], 1.0108448, places=6)
        self.assertLess(stats["angleVariance"]["mean"], 1e-20)
        self.assertAlmostEqual(stats["bondLength"]["mean"], mean, places=12)

    def test_a_trigonally_compressed_octahedron(self):
        """Six bonds of one length R at polar angle psi from the threefold axis:
        the 6 cis angles inside each triangle are theta_t = arccos((3 cos^2 psi - 1) / 2),
        the 6 between the triangles 180 - theta_t, so with delta = 90 - theta_t,
        sigma^2 = 12 delta^2 / 11; DI = 0; <lambda> = ((3 sqrt(3) / 2) sin^2 psi cos psi)^(-2/3)."""
        for psi in (50.0, 54.7356103172, 60.0):
            with self.subTest(psi=psi):
                result = run(configs.isolated_polyhedron(configs.trigonal_octahedron(psi, 2.0)))
                p = math.radians(psi)
                delta = 90.0 - math.degrees(math.acos((3 * math.cos(p) ** 2 - 1) / 2))
                stats = result["configuration"]
                self.assertLess(abs(stats["distortionIndex"]["mean"]), 1e-14)
                self.assertAlmostEqual(stats["angleVariance"]["mean"], 12 * delta**2 / 11, places=9)
                expected = ((3 * math.sqrt(3) / 2) * math.sin(p) ** 2 * math.cos(p)) ** (-2.0 / 3.0)
                self.assertAlmostEqual(stats["quadraticElongation"]["mean"], expected, places=12)
                # arccos near -1 turns the cosine's last bit into ~1e-6 deg.
                self.assertAlmostEqual(stats["transAngle"]["mean"], 180.0, places=5)
        # psi = 50 deg: delta = 6.88 deg, sigma^2 = 51.6 deg^2.
        self.assertAlmostEqual(12 * (90.0 - math.degrees(math.acos((3 * math.cos(math.radians(50)) ** 2 - 1) / 2))) ** 2 / 11, 51.6, places=1)

    def test_a_flattened_tetrahedron(self):
        """(±s, 0, c) and (0, ±s, -c) with s^2 + c^2 = 1, bonds R: two angles
        arccos(c^2 - s^2), four arccos(-c^2); V = (4/3) c s^2 R^3, so
        <lambda> = ((3 sqrt(3) / 2) c s^2)^(-2/3); DI = 0."""
        c = 0.5
        s2 = 1 - c * c
        result = run(configs.isolated_polyhedron(configs.flattened_tetrahedron(c, 2.0)))
        wide = math.degrees(math.acos(c * c - s2))
        narrow = math.degrees(math.acos(-c * c))
        self.assertAlmostEqual(wide, 120.0, places=10)
        variance = (2 * (wide - TETRAHEDRAL_ANGLE_DEG) ** 2 + 4 * (narrow - TETRAHEDRAL_ANGLE_DEG) ** 2) / 5
        stats = result["configuration"]
        self.assertAlmostEqual(stats["angleVariance"]["mean"], variance, places=9)
        self.assertAlmostEqual(stats["quadraticElongation"]["mean"], ((3 * math.sqrt(3) / 2) * c * s2) ** (-2.0 / 3.0), places=12)
        self.assertLess(abs(stats["distortionIndex"]["mean"]), 1e-14)

    def test_an_off_centred_atom_in_a_regular_octahedron(self):
        """The centre delta off the middle along [0 0 1]: <lambda> = 1 + (delta/R)^2
        (the volume does not move), and DI and sigma^2 from the six bonds
        (R -+ delta along the axis, sqrt(R^2 + delta^2) across it)."""
        radius, delta = 2.0, 0.2
        result = run(configs.isolated_polyhedron(configs.octahedron(), offset=(0.0, 0.0, delta)))
        lengths = [math.hypot(radius, delta)] * 4 + [radius - delta, radius + delta]
        mean = sum(lengths) / 6
        di = sum(abs(length - mean) for length in lengths) / 6 / mean
        hyp = math.hypot(radius, delta)
        equatorial = math.degrees(math.acos(delta**2 / hyp**2))
        upper = math.degrees(math.acos(-delta / hyp))
        lower = math.degrees(math.acos(delta / hyp))
        variance = 4 * ((equatorial - 90) ** 2 + (upper - 90) ** 2 + (lower - 90) ** 2) / 11
        stats = result["configuration"]
        self.assertAlmostEqual(stats["quadraticElongation"]["mean"], 1 + (delta / radius) ** 2, places=12)
        self.assertAlmostEqual(stats["distortionIndex"]["mean"], di, places=12)
        self.assertAlmostEqual(stats["angleVariance"]["mean"], variance, places=9)

    def test_the_per_polyhedron_function_on_arrays(self):
        vectors = np.array([configs.octahedron(), configs.octahedron((1.93,) * 4 + (2.25,) * 2)])
        lengths = np.linalg.norm(vectors, axis=2)
        indices = polyhedron_indices(vectors, lengths)
        self.assertEqual(indices["distortionIndex"].shape, (2,))
        self.assertLess(indices["distortionIndex"][0], 1e-14)
        self.assertAlmostEqual(indices["distortionIndex"][1], 0.06983087834151, places=12)
        # A triangle (n = 3) has DI only.
        three = polyhedron_indices(vectors[:, :3], lengths[:, :3])
        self.assertIsNone(three["quadraticElongation"])
        self.assertIsNone(three["angleVariance"])

    def test_a_flat_polyhedron_has_no_quadratic_elongation(self):
        """A square plane read as a tetrahedron: zero volume, so no <lambda>
        (left out of its summary), and the angle variance of 90/90/90/90/180/180."""
        square = [np.array(v) * 2.0 for v in ((1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0))]
        result = run(configs.isolated_polyhedron(square))
        stats = result["configuration"]
        self.assertEqual(result["polyhedron"]["shape"], "tetrahedron")
        self.assertIsNone(stats["quadraticElongation"])
        expected = (4 * (90 - TETRAHEDRAL_ANGLE_DEG) ** 2 + 2 * (180 - TETRAHEDRAL_ANGLE_DEG) ** 2) / 5
        self.assertAlmostEqual(stats["angleVariance"]["mean"], expected, places=9)
        self.assertIsNone(result["average"]["sites"][0]["quadraticElongation"])


class AverageStructureTests(unittest.TestCase):
    def test_orientationally_disordered_jahn_teller_octahedra_are_regular_on_average(self):
        """Every octahedron 4 + 2 with its long axis along x, y or z, each F site
        long in a third of the cells: the configuration's DI is the 4 + 2 value
        everywhere, the average octahedron is regular at <l> = 2.03667 A."""
        result = run(configs.jahn_teller_network(), window=(1.6, 2.6))
        stats = result["configuration"]
        self.assertEqual(stats["polyhedra"], 27)
        self.assertAlmostEqual(stats["distortionIndex"]["mean"], 0.06983087834151, places=10)
        self.assertLess(stats["distortionIndex"]["sd"], 1e-12)
        self.assertLess(stats["angleVariance"]["p90"], 1e-12)
        average = result["average"]
        self.assertEqual(average["polyhedra"], 1)
        self.assertLess(average["distortionIndex"]["mean"], 1e-12)
        self.assertAlmostEqual(average["quadraticElongation"]["mean"], 1.0, places=12)
        self.assertLess(average["angleVariance"]["mean"], 1e-16)
        np.testing.assert_allclose(average["sites"][0]["bondLengths"], [(2.25 + 2 * 1.93) / 3] * 6, atol=1e-12)
        # The first peak alone: 4 ligands in the configuration, 6 in the average.
        partial = run(configs.jahn_teller_network(), window=(1.6, 2.1))
        self.assertEqual(partial["polyhedron"]["coordination"], 4)
        self.assertEqual(partial["average"]["reason"], "coordination_differs")
        self.assertEqual(partial["average"]["polyhedra"], 0)
        self.assertIsNone(partial["average"]["distortionIndex"])
        self.assertEqual(partial["average"]["sites"][0]["coordination"], 6)

    def test_a_thermal_only_network_is_regular_on_average_and_its_di_follows_the_bond_spread(self):
        """Corner-sharing TiO6 octahedra, only the O displaced (isotropic Gaussian,
        0.05 A per axis, each site's displacements centred over its 1000 copies).

        The average octahedron is regular exactly. The six bonds of an
        octahedron are six different O atoms, so they are independent, each
        l_i = <l> + eps_i to first order in s/<l> with eps_i ~ N(0, s^2), s the
        bond-length sd: |l_i - mean_k l_k| has mean sqrt(2/pi) sqrt((n-1)/n) s,
        so E[DI] = sqrt(2/pi) sqrt(5/6) s / <l> = 0.7284 s / <l>. Tolerance:
        four standard errors of the mean DI over the polyhedra (sd / sqrt(M))
        plus 1 % for the second-order terms ((s/<l>)^2 ~ 7e-4).

        The angles follow the same displacements: a cis angle moves by the two
        tangential components (u_i . e_j + u_j . e_i) / <l>, so
        E[(theta - 90)^2] = 2 (s/<l>)^2 rad^2 and E[sigma^2_theta] =
        (12 / 11) 2 (s/<l>)^2 (180/pi)^2 (4.71 deg^2 here), to the same
        tolerance (four standard errors plus 2 %)."""
        result = run(configs.perovskite(supercell=(10, 10, 10), sigma_o=0.05, centred=True), "Ti", "O", (1.5, 2.5))
        stats = result["configuration"]
        self.assertEqual(stats["polyhedra"], 1000)
        spread, mean = stats["bondLength"]["sd"], stats["bondLength"]["mean"]
        expected = math.sqrt(2 / math.pi) * math.sqrt(5 / 6) * spread / mean
        tolerance = 4 * stats["distortionIndex"]["sd"] / math.sqrt(stats["polyhedra"]) + 0.01 * expected
        self.assertLess(abs(stats["distortionIndex"]["mean"] - expected), tolerance)
        self.assertLess(tolerance, 0.06 * expected)
        variance = (12 / 11) * 2 * (spread / mean) ** 2 * (180 / math.pi) ** 2
        self.assertAlmostEqual(variance, 4.7, delta=0.15)
        tolerance = 4 * stats["angleVariance"]["sd"] / math.sqrt(stats["polyhedra"]) + 0.02 * variance
        self.assertLess(abs(stats["angleVariance"]["mean"] - variance), tolerance)
        # The configuration's polyhedra are distorted; their average is regular.
        self.assertGreater(stats["distortionIndex"]["mean"], 0.01)
        self.assertGreater(stats["angleVariance"]["mean"], 1.0)
        average = result["average"]
        self.assertEqual(average["polyhedra"], 1)
        self.assertLess(average["distortionIndex"]["mean"], 1e-12)
        self.assertAlmostEqual(average["quadraticElongation"]["mean"], 1.0, places=12)
        self.assertLess(average["angleVariance"]["mean"], 1e-16)

    def test_uncentred_thermal_motion_leaves_the_average_nearly_regular(self):
        """Without the centring, each site mean carries its sampling error
        (0.05 / sqrt(64) A per axis over 64 copies): the average octahedron's DI
        is of the configuration's divided by sqrt(64), not the configuration's."""
        result = run(configs.perovskite(supercell=(4, 4, 4), sigma_o=0.05, seed=5), "Ti", "O", (1.5, 2.5))
        self.assertLess(result["average"]["distortionIndex"]["mean"], 0.3 * result["configuration"]["distortionIndex"]["mean"])

    def test_a_static_distortion_shows_in_the_average_too(self):
        """A sheared (triclinic) cell distorts every octahedron alike: the average
        polyhedron carries the configuration's distortion."""
        cell = [[3.9, 0.0, 0.0], [0.25, 3.9, 0.0], [0.15, 0.2, 3.95]]
        result = run(configs.perovskite(unit_cell=cell, sigma_o=0.02, seed=3), "Ti", "O", (1.5, 2.5))
        config_di = result["configuration"]["distortionIndex"]["mean"]
        average_di = result["average"]["distortionIndex"]["mean"]
        self.assertGreater(average_di, 0.5 * config_di)
        self.assertGreater(result["average"]["angleVariance"]["mean"], 1.0)

    def test_a_mixed_site_counts_as_its_majority_species(self):
        b_site = lambda rng, cell: "Zr" if (cell[0] + 2 * cell[1] + 3 * cell[2]) % 4 == 0 else "Ti"  # noqa: E731
        built = configs.perovskite(b_site=b_site, sigma_o=0.03, seed=2)
        titanium = run(built, "Ti", "O", (1.5, 2.5))
        self.assertEqual(titanium["average"]["polyhedra"], 1)
        self.assertTrue(titanium["average"]["sites"][0]["mixed"])
        self.assertEqual(titanium["average"]["sites"][0]["element"], "Ti")
        self.assertEqual(titanium["average"]["sites"][0]["copies"], 64)
        zirconium = run(built, "Zr", "O", (1.5, 2.5))
        self.assertEqual(zirconium["configuration"]["polyhedra"], 16)
        self.assertEqual(zirconium["average"]["reason"], "no_centre_sites")
        self.assertEqual(zirconium["average"]["sites"], [])

    def test_a_coords_only_file_has_no_average_structure(self):
        atoms, lattice, supercell = configs.perovskite(sigma_o=0.03)
        text = configs.rmc6f_text(atoms, lattice, supercell, coords_only=True)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "run.rmc6f"
            path.write_text(text, encoding="utf-8")
            result = polyhedral_distortion_from_file(path, "ti", "o", 1.5, 2.5)
        self.assertEqual(result["average"]["reason"], "no_basis")
        self.assertEqual(result["average"]["polyhedra"], 0)
        self.assertEqual(result["configuration"]["polyhedra"], 64)
        self.assertIsNone(result["parseWarning"])

    def test_a_box_smaller_than_the_window(self):
        """One perovskite cell (3.9 A) holds each TiO6 only through periodic
        images: the six O are three sites, each twice."""
        result = run(configs.perovskite(supercell=(1, 1, 1)), "Ti", "O", (1.5, 2.5))
        self.assertEqual(result["configuration"]["polyhedra"], 1)
        self.assertLess(result["configuration"]["distortionIndex"]["mean"], 1e-14)
        self.assertEqual(result["average"]["sites"][0]["coordination"], 6)


class CoordinationTests(unittest.TestCase):
    def test_the_most_common_coordination_is_measured_and_the_rest_counted(self):
        """The O at (0, 1/2, 1/2) of every fourth cell along a moved 0.4 A along
        -x: 2.35 A from its own cell's Ti, beyond the window (5 O), 1.55 A from
        the previous cell's (still 6). The 48 octahedra are measured and the
        16 five-coordinated Ti counted."""
        def displace(rng, cell, site, element):
            if site == 4 and cell[0] % 4 == 0:
                return np.array([-0.4, 0.0, 0.0])
            return np.zeros(3)

        built = configs.supercell_atoms(np.eye(3) * 3.9, configs.PEROVSKITE_SITES, (4, 4, 4), displace=displace)
        result = run(built, "Ti", "O", (1.4, 2.3))
        self.assertEqual(result["centres"], 64)
        self.assertEqual(result["coordination"], [0, 0, 0, 0, 0, 16, 48])
        self.assertEqual(result["polyhedron"]["coordination"], 6)
        self.assertEqual(result["configuration"]["polyhedra"], 48)
        self.assertEqual(result["skipped"], 16)
        self.assertEqual(result["coordinationSummary"]["modal"], 6)
        self.assertEqual(result["coordinationSummary"]["modalShare"], 0.75)

    def test_a_tie_goes_to_the_lowest_count(self):
        sites = [("Cu", (0.25, 0.5, 0.5)), ("Cu", (0.75, 0.5, 0.5))]
        sites += [("F", (0.25 + 0.2 * v[0], 0.5 + 0.2 * v[1], 0.5 + 0.2 * v[2])) for v in configs.AXES[:4]]
        sites += [("F", (0.75 + 0.2 * v[0], 0.5 + 0.2 * v[1], 0.5 + 0.2 * v[2])) for v in configs.AXES[:2]]
        built = configs.supercell_atoms(np.eye(3) * 10.0, sites, (1, 1, 1))
        result = run(built, window=(1.0, 2.5))
        self.assertEqual(result["coordination"], [0, 0, 1, 0, 1])
        self.assertEqual(result["polyhedron"]["coordination"], 2)


class RefusalTests(unittest.TestCase):
    built = configs.isolated_polyhedron(configs.octahedron())

    def test_a_window_must_be_ordered(self):
        for window in ((3.0, 2.0), (2.0, 2.0), (-0.1, 2.0), (float("nan"), 2.0), (1.0, float("inf"))):
            with self.subTest(window=window):
                with self.assertRaisesRegex(ValueError, "window"):
                    run(self.built, window=window)

    def test_an_unknown_element_names_the_elements(self):
        with self.assertRaisesRegex(ValueError, r"unknown element 'Xx' for the centre: the configuration holds Cu, F"):
            run(self.built, centre="xx")
        with self.assertRaisesRegex(ValueError, "element 'Yy' for the ligand"):
            run(self.built, ligand="YY")

    def test_no_polyhedra(self):
        with self.assertRaisesRegex(ValueError, r"no polyhedra: no Cu atom has a F neighbour within \[0.5, 1\] A"):
            run(self.built, window=(0.5, 1.0))
        one = configs.isolated_polyhedron([np.array([2.0, 0.0, 0.0])])
        with self.assertRaisesRegex(ValueError, "no polyhedra: the most common Cu coordination .* is one F"):
            run(one)

    def test_the_bond_budget_is_checked_before_any_bond_is_stored(self):
        self.assertEqual(run(self.built, max_bonds=48)["configuration"]["polyhedra"], 8)
        with self.assertRaisesRegex(ValueError, "holds 48 bonds, over the limit of 47 for one request; narrow the window"):
            run(self.built, max_bonds=47)

    def test_elements_and_window_are_echoed_as_read(self):
        result = run(self.built, centre=" cu ", ligand="f", window=(1, 3.5))
        self.assertEqual((result["centre"], result["ligand"], result["window"]), ("Cu", "F", [1.0, 3.5]))


if __name__ == "__main__":
    unittest.main()
