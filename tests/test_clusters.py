# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The cluster-distortion engine on synthetic configurations with planted distortions."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cluster_configs import chains_config_text, cluster_config_text, origin_dimer_config_text, with_vacancies  # noqa: E402
from rmc_toolkits.clusters import (  # noqa: E402
    cluster_distortion,
    cluster_distortion_from_file,
    lattice_families,
    mix32,
    shuffle_permutation,
)
from rmc_toolkits.pca_kde import load_site_displacements  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEMO_RMC6F = ROOT / "web_app" / "frontend" / "public" / "demo" / "GTS_250K.rmc6f"


class _Configs:
    """Each synthetic mode written once to a temporary .rmc6f and loaded."""

    _dir = tempfile.TemporaryDirectory()
    _cache: dict[str, object] = {}

    @classmethod
    def sites(cls, mode: str):
        if mode not in cls._cache:
            path = Path(cls._dir.name) / f"{mode}.rmc6f"
            path.write_text(cluster_config_text(mode, supercell=(6, 6, 6)), encoding="utf-8")
            cls._cache[mode] = load_site_displacements(path)
        return cls._cache[mode]


class ShuffleHashTests(unittest.TestCase):
    def test_mix32_is_the_lowbias32_hash(self):
        # Reference values of lowbias32 (x ^= x>>16; x *= 0x7feb352d; x ^= x>>15; x *= 0x846ca68b; x ^= x>>16).
        self.assertEqual(int(mix32(np.uint32(0))), 0)
        self.assertEqual(int(mix32(np.uint32(1))), 0x688990C0)
        self.assertEqual(int(mix32(np.uint32(12345))), 0x912EFCF7)
        values = mix32(np.arange(1000, dtype=np.uint32))
        self.assertEqual(len(set(values.tolist())), 1000)

    def test_permutations_are_deterministic_and_distinct(self):
        first = shuffle_permutation(500, 0, 0, 0)
        self.assertEqual(sorted(first.tolist()), list(range(500)))
        np.testing.assert_array_equal(first, shuffle_permutation(500, 0, 0, 0))
        for other in (shuffle_permutation(500, 1, 0, 0), shuffle_permutation(500, 0, 1, 0), shuffle_permutation(500, 0, 0, 1)):
            self.assertFalse(np.array_equal(first, other))


class LatticeFamilyTests(unittest.TestCase):
    def test_cubic_families(self):
        families = lattice_families(np.eye(3) * 8.0)
        self.assertEqual({name: len(v) for name, v in families.items()}, {"100": 3, "110": 6, "111": 4})
        np.testing.assert_allclose(families["111"][0], np.ones(3) / np.sqrt(3.0))
        for vectors in families.values():
            np.testing.assert_allclose(np.linalg.norm(vectors, axis=1), 1.0)


class ClusterDistortionTests(unittest.TestCase):
    def test_templates_from_the_average_structure(self):
        result = cluster_distortion(_Configs.sites("independent"), "nb", (2.5, 3.5))
        self.assertEqual(result["element"], "Nb")
        self.assertFalse(result["extended"])
        self.assertEqual(result["kind"], {"members": 4, "edges": 6, "clusters": 216, "templates": 1})
        template = result["templates"][0]
        self.assertEqual(template["referenceNumbers"], [1, 2, 3, 4])
        self.assertEqual(template["edges"], [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
        self.assertEqual(template["incomplete"], 0)
        np.testing.assert_allclose(template["idealEdges"], 3.0, atol=0.03)
        self.assertEqual(result["outliers"]["extraContacts"], 0)

    def test_independent_atoms_match_the_null(self):
        result = cluster_distortion(_Configs.sites("independent"), "Nb", (2.5, 3.5))
        z = result["z"]
        for key in ("ladderWidth", "shapeRms", "breathingSd", "eRms", "t2Rms", "t2Trigonality"):
            self.assertLess(abs(z[key]), 3.0, key)
        self.assertLess(result["ideal"]["spreadOverSe"], 5.0)

    def test_a_local_trigonal_distortion_is_found_and_classified(self):
        # One random vertex per cluster pushed out along its 3-fold axis: the
        # average tetrahedron stays regular, the clusters do not.
        result = cluster_distortion(_Configs.sites("trigonal"), "Nb", (2.5, 3.5))
        z = result["z"]
        self.assertLess(result["ideal"]["spreadOverSe"], 5.0)
        self.assertLess(z["t2Trigonality"], -3.0)
        self.assertLess(result["observed"]["t2Trigonality"], result["nullMean"]["t2Trigonality"])
        self.assertGreater(z["longAxis"]["111"], 3.0)
        # Exactly one vertex out per cluster: less size variation than a
        # binomial number of them (the null).
        self.assertLess(z["breathingSd"], -3.0)

    def test_a_static_distortion_shows_in_the_average_cluster_only(self):
        result = cluster_distortion(_Configs.sites("static"), "Nb", (2.5, 3.5))
        self.assertGreater(result["ideal"]["spreadOverSe"], 10.0)
        self.assertGreater(result["ideal"]["spread"], 0.1)
        self.assertLess(abs(result["z"]["ladderWidth"]), 3.0)

    def test_dimers(self):
        result = cluster_distortion(_Configs.sites("independent"), "Ga", (2.0, 3.0))
        self.assertEqual(result["kind"]["members"], 2)
        self.assertEqual(result["kind"]["edges"], 1)
        self.assertNotIn("t2Rms", result["observed"])
        self.assertNotIn("shortAxis", result["observed"])
        # The bond lies along x: every dimer's axis is near <100>.
        self.assertEqual(result["observed"]["longAxis"]["100"], 1.0)

    def test_a_site_at_fraction_zero_pairs_atoms_of_one_cell(self):
        # Review 2026-10-12 #1: a site mean a hair below 0 wrapped to 0.9999 and
        # the template paired each atom with its partner of the next cell: the
        # rigid dimers' correlation vanished (z -32 -> -0.6) and edges broke.
        results = {}
        with tempfile.TemporaryDirectory() as directory:
            for bias in (1e-4, -1e-4):
                path = Path(directory) / f"origin{bias}.rmc6f"
                path.write_text(origin_dimer_config_text(bias), encoding="utf-8")
                sites = load_site_displacements(path)
                self.assertAlmostEqual(sites.site_cell_fractional[0][0], bias / 5.0, places=9)
                results[bias] = cluster_distortion(sites, "V", (2.0, 2.5))
        for bias, result in results.items():
            self.assertEqual(result["templates"][0]["offsets"], [[0, 0, 0], [0, 0, 0]], bias)
            self.assertEqual(result["outliers"]["brokenEdgeFraction"], 0.0, bias)
            self.assertLess(result["z"]["breathingSd"], -10.0, bias)
        self.assertAlmostEqual(results[-1e-4]["z"]["breathingSd"], results[1e-4]["z"]["breathingSd"], delta=0.01)

    def test_templates_of_one_kind_are_measured_with_their_own_edges(self):
        # Review 2026-10-12 #2: chain B (bonded 0-1, 0-2) was measured with
        # chain A's edges (0-1, 1-2), so one of its "edges" was its 4 A end to
        # end: ladderWidth 1.03 A for chains whose bonds differ by 0.01 A.
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "chains.rmc6f"
            path.write_text(chains_config_text(), encoding="utf-8")
            result = cluster_distortion(load_site_displacements(path), "Mo", (1.8, 2.2))
        self.assertEqual([t["edges"] for t in result["templates"]], [[[0, 1], [1, 2]], [[0, 1], [0, 2]]])
        self.assertEqual(result["kind"], {"members": 3, "edges": 2, "clusters": 128, "templates": 2})
        self.assertLess(result["observed"]["ladderWidth"], 0.2)
        self.assertAlmostEqual(result["observed"]["meanEdge"], 2.0, delta=0.01)
        self.assertEqual(result["outliers"]["brokenEdgeFraction"], 0.0)

    def test_missing_atoms_are_not_contacts_between_clusters(self):
        # Review 2026-10-12 #16: the bonds of a cluster that lost an atom were
        # counted as contacts between clusters (102 of them, 0.25 per atom).
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "vacancies.rmc6f"
            path.write_text(with_vacancies(cluster_config_text("independent", supercell=(6, 6, 6)), "Nb", 0.05), encoding="utf-8")
            result = cluster_distortion(load_site_displacements(path), "Nb", (2.5, 3.5))
        self.assertGreater(result["templates"][0]["incomplete"], 0)
        self.assertEqual(result["outliers"]["extraContacts"], 0)

    def test_an_extended_network_and_an_empty_window(self):
        extended = cluster_distortion(_Configs.sites("independent"), "Se", (7.5, 8.0))
        self.assertTrue(extended["extended"])
        self.assertEqual(extended["templates"], [])
        self.assertNotIn("observed", extended)
        empty = cluster_distortion(_Configs.sites("independent"), "Nb", (0.5, 1.0))
        self.assertFalse(empty["extended"])
        self.assertEqual(empty["templates"], [])
        self.assertEqual(empty["isolatedSites"], 4)

    def test_invalid_requests(self):
        sites = _Configs.sites("independent")
        with self.assertRaisesRegex(ValueError, "no site holds Xx"):
            cluster_distortion(sites, "Xx", (2.5, 3.5))
        with self.assertRaisesRegex(ValueError, "window must be"):
            cluster_distortion(sites, "Nb", (3.5, 2.5))
        with self.assertRaisesRegex(ValueError, "shuffles"):
            cluster_distortion(sites, "Nb", (2.5, 3.5), shuffles=1)
        from dataclasses import replace
        with self.assertRaisesRegex(ValueError, "cell columns"):
            cluster_distortion(replace(sites, cell_indices=None), "Nb", (2.5, 3.5))


@unittest.skipUnless(DEMO_RMC6F.exists(), "the committed demo run is missing")
class DemoRunTests(unittest.TestCase):
    def test_the_demo_ta4_clusters(self):
        result = cluster_distortion_from_file(str(DEMO_RMC6F), "Ta", 2.65, 3.33)
        self.assertEqual(result["kind"], {"members": 4, "edges": 6, "clusters": 4000, "templates": 4})
        # The RMC clusters are more rigid than independent atoms would make them.
        self.assertLess(result["z"]["ladderWidth"], -10)
        self.assertLess(result["observed"]["ladderWidth"], 0.9 * result["nullMean"]["ladderWidth"])


if __name__ == "__main__":
    unittest.main()
