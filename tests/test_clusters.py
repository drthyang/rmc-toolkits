# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The cluster-distortion engine on synthetic configurations with planted distortions."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cluster_configs import (  # noqa: E402
    chains_config_text,
    cluster_config_text,
    octahedra_config_text,
    origin_dimer_config_text,
    rigid_dimer_config_text,
    tie_config_text,
    with_vacancies,
)
from rmc_toolkits.clusters import (  # noqa: E402
    AXIS_ENHANCEMENT,
    CHANGE_KINDS,
    CLUSTER_BECAUSE,
    CLUSTER_EFFECT,
    CLUSTER_RULES,
    CLUSTER_Z,
    EVIDENCE_AGAINST,
    FAMILIES,
    SPLIT_SITE_KURTOSIS,
    SPLIT_SITE_SE,
    cluster_distortion,
    cluster_distortion_from_file,
    evidence_clears,
    is_split_site,
    lattice_families,
    mix32,
    shuffle_permutation,
    split_site_evidence,
)
from rmc_toolkits.pca_kde import load_site_displacements, site_ellipsoids  # noqa: E402

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
        # The shape null (each cluster's breathing kept) matches too: the
        # radial noise its per-cluster fit takes out is put back.
        shape = result["shapeNull"]
        for key in ("ladderWidth", "shapeRms", "eRms", "t2Rms", "t2Trigonality"):
            self.assertLess(abs(shape["z"][key]), 3.0, key)
            self.assertLess(abs(result["observed"][key] / shape["mean"][key] - 1.0), 0.03, key)
        self.assertNotIn("breathingSd", shape["mean"])
        self.assertNotIn("meanEdge", shape["z"])

    def test_breathing_is_size_not_shape(self):
        # Review 2026-10-14 L2: every tetrahedron scaled by 1 + eps, sd(eps)
        # 0.04, its shape regular. Against independent atoms the size
        # fluctuates more, and the shape reads "less distorted" with a
        # spurious axis preference (the null spreads the breathing over
        # independent atoms); against the shape null the shape is regular.
        result = cluster_distortion(_Configs.sites("breathing"), "Nb", (2.5, 3.5))
        observed, null, z, shape = result["observed"], result["nullMean"], result["z"], result["shapeNull"]
        self.assertGreater(z["breathingSd"], 3.0)
        self.assertGreater(observed["breathingSd"] / null["breathingSd"] - 1.0, 0.05)
        self.assertLess(z["ladderWidth"], -3.0)
        self.assertLess(observed["ladderWidth"] / null["ladderWidth"] - 1.0, -0.05)
        for key in ("ladderWidth", "shapeRms", "t2Rms", "t2Trigonality"):
            self.assertLess(abs(shape["z"][key]), 3.0, key)
            self.assertLess(abs(observed[key] / shape["mean"][key] - 1.0), 0.05, key)
        for axis in ("longAxis", "shortAxis"):
            for family in FAMILIES:
                self.assertLess(abs(shape["z"][axis][family] or 0.0), 3.0, (axis, family))

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
        # At the clusters' own size the vertex sites' two positions hold the
        # edge spread; the pattern still shows: trigonal, long axis <111>.
        shape = result["shapeNull"]
        self.assertLess(shape["z"]["t2Trigonality"], -3.0)
        self.assertGreater(shape["z"]["longAxis"]["111"], 3.0)
        self.assertLess(abs(shape["z"]["ladderWidth"]), 3.0)

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

    def test_a_tie_lists_the_analysed_kind_first(self):
        # Review 2026-10-12 #17: 64 triangles and 64 dimers; the dimers are
        # analysed (the smaller kind wins a tie), so they must lead the list.
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tie.rmc6f"
            path.write_text(tie_config_text(), encoding="utf-8")
            result = cluster_distortion(load_site_displacements(path), "Mo", (1.8, 2.2))
        self.assertEqual(result["kind"]["members"], 2)
        self.assertEqual(result["kinds"], [
            {"members": 2, "edges": 1, "clusters": 64},
            {"members": 3, "edges": 3, "clusters": 64},
        ])

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


def _load_text(text: str, name: str = "config"):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / f"{name}.rmc6f"
        path.write_text(text, encoding="utf-8")
        return load_site_displacements(path)


def _answers(verdict: dict) -> dict:
    """Every finding as (answer, because), the window's as window.<name>."""
    out = {}
    for name, finding in verdict.items():
        if name == "rules":
            continue
        if name == "window":
            out.update({f"window.{key}": (value["answer"], value["because"]) for key, value in finding.items()})
        else:
            out[name] = (finding["answer"], finding["because"])
    return out


def _records(verdict: dict):
    """Every Evidence record of a verdict, with the finding it belongs to."""
    for name, finding in verdict.items():
        if name == "rules":
            continue
        findings = finding.items() if name == "window" else [(name, finding)]
        for key, item in findings:
            for record in item["evidence"]:
                yield (key if name != "window" else f"window.{key}"), record


def _clears_as_documented(record: dict) -> bool:
    """The Evidence rule as the docs state it, written out apart from evidence_clears."""
    side, strict = record.get("side"), record.get("strict", False)
    if not record["bar"]:
        return False
    for field, bar in record["bar"].items():
        value = record[field]
        if value is None:
            return False
        if side == "above" and not (value > bar if strict else value >= bar):
            return False
        if side == "below" and not (value < -bar if strict else value <= -bar):
            return False
        if side is None and not (abs(value) > bar if strict else abs(value) >= bar):
            return False
    return True


class ClusterVerdictTests(unittest.TestCase):
    """The engine's own verdict, on configurations with known answers."""

    def verdict(self, mode, element="Nb", window=(2.5, 3.5)):
        return cluster_distortion(_Configs.sites(mode), element, window)["verdict"]

    def check_records(self, verdict):
        for name, record in _records(verdict):
            with self.subTest(finding=name, statistic=record["statistic"]):
                self.assertIn(record["against"], EVIDENCE_AGAINST)
                self.assertIn(record["changeKind"], CHANGE_KINDS + (None,))
                self.assertEqual(record["clears"], _clears_as_documented(record))
                self.assertEqual(record["clears"], evidence_clears(record))
                # change is the record's own observed and null, by its kind.
                if record["changeKind"] == "ratio" and record["change"] is not None:
                    self.assertEqual(record["change"], record["observed"] / record["null"] - 1.0)
                if record["changeKind"] == "difference" and record["change"] is not None:
                    self.assertEqual(record["change"], record["observed"] - record["null"])

    def check_shape(self, verdict):
        """Rules echoed, every because in its closed set, the clearing evidence first, strict JSON."""
        self.assertEqual(verdict["rules"], CLUSTER_RULES)
        for name, (answer, because) in _answers(verdict).items():
            self.assertIn(answer, (True, False, None), name)
            self.assertIn(because, CLUSTER_BECAUSE[name], name)
        for family in verdict["favouredAxes"]["families"]:
            self.assertIn(family["because"], CLUSTER_BECAUSE["favouredAxes.families"])
        for name, finding in verdict.items():
            if name in ("rules", "window"):
                continue
            clears = [record["clears"] for record in finding["evidence"]]
            self.assertEqual(clears, sorted(clears, reverse=True), name)
        json.dumps(verdict, allow_nan=False)
        self.check_records(verdict)

    def test_independent_atoms_find_nothing(self):
        verdict = self.verdict("independent")
        self.check_shape(verdict)
        answers = _answers(verdict)
        for name in ("distorted", "lessDistorted", "breathing", "rigid", "trigonalLean", "favouredAxes"):
            self.assertEqual(answers[name], (False, "within_null"), name)
        self.assertEqual(answers["staticDistortion"], (False, "within_sampling_error"))
        self.assertEqual(answers["splitMemberSites"], (False, "none_split"))
        self.assertEqual(answers["window.extendedNetwork"], (False, "finite_clusters"))
        for name in ("mixedSizes", "bondsOutsideWindow", "contactsBetweenClusters"):
            self.assertEqual(answers[f"window.{name}"], (False, "within_bar"), name)
        # A negative keeps every statistic it consulted: the edge spread, the
        # trigonality and the six axis families.
        self.assertEqual(
            [record["statistic"] for record in verdict["distorted"]["evidence"]],
            ["ladderWidth", "t2Trigonality", "longAxis.100", "longAxis.110", "longAxis.111",
             "shortAxis.100", "shortAxis.110", "shortAxis.111"],
        )

    def test_breathing_is_size_not_shape(self):
        # Against independent atoms the edge spread is far below (z < -3, the
        # breathing spread over independent atoms), against the shape null it
        # is not: the verdict reads the shape null (review 2026-10-14 L2).
        result = cluster_distortion(_Configs.sites("breathing"), "Nb", (2.5, 3.5))
        verdict = result["verdict"]
        self.check_shape(verdict)
        self.assertLess(result["z"]["ladderWidth"], -3.0)
        width = verdict["lessDistorted"]["evidence"][0]
        self.assertEqual((width["statistic"], width["against"]), ("ladderWidth", "shapeNull"))
        self.assertLess(abs(width["z"]), 3.0)
        answers = _answers(verdict)
        self.assertEqual(answers["breathing"], (True, "above_null"))
        self.assertEqual(answers["distorted"], (False, "within_null"))
        self.assertEqual(answers["lessDistorted"], (False, "within_null"))
        self.assertEqual(answers["rigid"], (False, "within_null"))
        self.assertEqual(answers["favouredAxes"], (False, "within_null"))
        breathing = verdict["breathing"]["evidence"][0]
        self.assertEqual((breathing["against"], breathing["side"], breathing["clears"]), ("independentNull", "above", True))

    def test_a_trigonal_distortion_is_a_pattern(self):
        verdict = self.verdict("trigonal")
        self.check_shape(verdict)
        answers = _answers(verdict)
        self.assertEqual(answers["distorted"], (True, "axis_family"))
        self.assertEqual(answers["trigonalLean"], (True, "leans_trigonal"))
        self.assertEqual(answers["staticDistortion"], (False, "within_sampling_error"))
        favoured = [(f["axis"], f["family"]) for f in verdict["favouredAxes"]["families"] if f["answer"]]
        self.assertIn(("long", "111"), favoured)
        # The clearing evidence leads.
        self.assertTrue(verdict["distorted"]["evidence"][0]["clears"])

    def test_a_pinched_edge_is_distorted_beyond_independent_atoms(self):
        # Vertices 0 and 1 move along their edge by +t and -t: the spread of
        # the edges is above the shape null, and edge 0-1 differs from its
        # opposite (the T2 pattern of one edge, not three at a vertex).
        verdict = self.verdict("pinched")
        self.check_shape(verdict)
        answers = _answers(verdict)
        self.assertEqual(answers["distorted"], (True, "edge_spread"))
        self.assertEqual(answers["lessDistorted"], (False, "above_null"))
        self.assertEqual(answers["trigonalLean"], (True, "leans_one_edge"))
        width = verdict["distorted"]["evidence"][0]
        self.assertEqual(width["statistic"], "ladderWidth")
        self.assertGreaterEqual(width["z"], CLUSTER_Z)
        self.assertGreaterEqual(width["change"], CLUSTER_EFFECT)

    def test_rigid_clusters_are_less_distorted(self):
        verdict = self.verdict("rigid")
        self.check_shape(verdict)
        answers = _answers(verdict)
        self.assertEqual(answers["lessDistorted"], (True, "below_null"))
        self.assertEqual(answers["distorted"], (False, "less_distorted"))
        self.assertEqual(answers["rigid"], (True, "less_distorted"))
        self.assertEqual(answers["breathing"], (False, "below_null"))
        self.assertEqual([r["statistic"] for r in verdict["rigid"]["evidence"]], ["ladderWidth", "breathingSd"])

    def test_a_static_distortion_is_the_average_cluster(self):
        verdict = self.verdict("static")
        self.check_shape(verdict)
        answers = _answers(verdict)
        self.assertEqual(answers["staticDistortion"], (True, "above_sampling_error"))
        self.assertEqual(answers["distorted"], (False, "within_null"))
        static = verdict["staticDistortion"]["evidence"][0]
        self.assertEqual((static["against"], static["bar"]), ("samplingError", {"z": 5}))
        self.assertGreater(static["z"], 10.0)

    def test_dimers(self):
        verdict = self.verdict("independent", "Ga", (2.0, 3.0))
        self.check_shape(verdict)
        answers = _answers(verdict)
        self.assertEqual(answers["lessDistorted"], (None, "single_edge"))
        self.assertEqual(answers["trigonalLean"], (None, "not_tetrahedra"))
        self.assertEqual(answers["breathing"], (False, "within_null"))
        # Every dimer lies along x: the axis families have no spread, no test.
        self.assertEqual(answers["favouredAxes"], (None, "undefined_statistic"))
        self.assertEqual(answers["distorted"], (None, "undefined_statistic"))

    def test_rigid_dimers(self):
        verdict = cluster_distortion(_load_text(rigid_dimer_config_text()), "V", (2.0, 2.5))["verdict"]
        self.check_shape(verdict)
        answers = _answers(verdict)
        self.assertEqual(answers["rigid"], (True, "stiffer_bonds"))
        self.assertEqual(answers["breathing"], (False, "below_null"))
        breathing = verdict["rigid"]["evidence"][0]
        self.assertEqual((breathing["statistic"], breathing["side"], breathing["clears"]), ("breathingSd", "below", True))
        self.assertLess(breathing["change"], -0.5)

    def test_split_octahedra(self):
        # Jahn-Teller F6 octahedra, long axis random: every F site has two
        # positions, the null keeps them (so the edge spread reads LESS
        # distorted), and the pattern shows the distortion (axis families).
        sites = _load_text(octahedra_config_text())
        result = cluster_distortion(sites, "F", (2.5, 3.25))
        verdict = result["verdict"]
        self.check_shape(verdict)
        self.assertEqual(result["kind"], {"members": 6, "edges": 12, "clusters": 216, "templates": 1})
        answers = _answers(verdict)
        self.assertEqual(answers["splitMemberSites"], (True, "half_or_more_split"))
        self.assertEqual(verdict["splitMemberSites"]["sites"], [2, 3, 4, 5, 6, 7])
        self.assertEqual(verdict["splitMemberSites"]["memberSites"], [2, 3, 4, 5, 6, 7])
        self.assertEqual(answers["distorted"], (True, "axis_family"))
        self.assertEqual(answers["lessDistorted"], (True, "below_null"))
        self.assertEqual(answers["rigid"], (False, "split_sites"))
        # The per-site evidence is the one rule list_sites applies.
        rows = {row["referenceNumber"]: row for row in site_ellipsoids(sites)}
        for record in verdict["splitMemberSites"]["evidence"]:
            self.assertEqual(record["clears"], is_split_site(rows[record["site"]]))
            self.assertEqual(record, split_site_evidence(rows[record["site"]]))

    def test_the_window(self):
        tie = cluster_distortion(_load_text(tie_config_text()), "Mo", (1.8, 2.2))["verdict"]
        self.check_shape(tie)
        self.assertEqual(_answers(tie)["window.mixedSizes"], (True, "above_bar"))
        share = tie["window"]["mixedSizes"]["evidence"][0]
        self.assertEqual((share["statistic"], share["observed"], share["strict"]), ("otherKindsShare", 0.5, True))
        # A window narrower than the bonds' thermal spread: bonds fall outside.
        narrow = self.verdict("independent", "Nb", (2.95, 3.05))
        self.check_shape(narrow)
        self.assertEqual(_answers(narrow)["window.bondsOutsideWindow"], (True, "above_bar"))

    def test_an_extended_network_and_no_clusters(self):
        extended = self.verdict("independent", "Se", (7.5, 8.0))
        empty = self.verdict("independent", "Nb", (0.5, 1.0))
        for verdict, because in ((extended, "extended_network"), (empty, "no_clusters")):
            self.check_shape(verdict)
            for name, (answer, code) in _answers(verdict).items():
                if name == "window.extendedNetwork":
                    self.assertEqual((answer, code), (because == "extended_network", because))
                else:
                    self.assertEqual((answer, code), (None, because), name)
            self.assertEqual(verdict["splitMemberSites"]["sites"], [])
            self.assertEqual(verdict["favouredAxes"]["families"], [])


class SplitRuleTests(unittest.TestCase):
    """The split rule (one definition: list_sites and the split member sites)."""

    @staticmethod
    def row(kappa, count=4000, resolved=True):
        return {"referenceNumber": 1, "count": count, "axisResolved": [resolved, True, True],
                "excessKurtosis": [kappa, 0.0, 0.0]}

    def test_the_bound_and_its_sampling_errors(self):
        self.assertEqual((SPLIT_SITE_KURTOSIS, SPLIT_SITE_SE), (-0.5, 3))
        self.assertTrue(is_split_site(self.row(SPLIT_SITE_KURTOSIS - 0.01)))
        self.assertFalse(is_split_site(self.row(SPLIT_SITE_KURTOSIS + 0.01)))
        # 216 copies: three sampling errors are sqrt(24/216) * 3 = 1.0.
        self.assertFalse(is_split_site(self.row(-0.99, count=216)))
        self.assertTrue(is_split_site(self.row(-1.01, count=216)))
        record = split_site_evidence(self.row(-1.01, count=216))
        self.assertAlmostEqual(record["z"], -1.01 / (24 / 216) ** 0.5, places=12)
        self.assertEqual((record["against"], record["side"], record["bar"]), ("siteKurtosis", "below", {"z": 3, "change": 0.5}))

    def test_an_unresolved_or_undefined_axis_is_never_split(self):
        self.assertFalse(is_split_site(self.row(-1.5, resolved=False)))
        self.assertIsNone(split_site_evidence(self.row(-1.5, resolved=False))["z"])
        self.assertFalse(is_split_site(self.row(None)))
        self.assertFalse(is_split_site(self.row(-1.5, count=0)))

    def test_the_rule_list_sites_had(self):
        # Away from the bounds, z = kappa / sqrt(24/n) <= -3 is kappa <= -3 sqrt(24/n),
        # the comparison the Copilot's isSplitSite made before the rule moved here.
        for kappa in (-2.0, -1.2, -0.8, -0.51, -0.49, -0.2, 0.3):
            for count in (64, 216, 512, 4000):
                old = kappa <= SPLIT_SITE_KURTOSIS and kappa <= -SPLIT_SITE_SE * (24 / count) ** 0.5
                self.assertEqual(is_split_site(self.row(kappa, count=count)), old, (kappa, count))


@unittest.skipUnless(DEMO_RMC6F.exists(), "the committed demo run is missing")
class DemoRunTests(unittest.TestCase):
    def test_the_demo_ta4_clusters(self):
        result = cluster_distortion_from_file(str(DEMO_RMC6F), "Ta", 2.65, 3.33)
        self.assertEqual(result["kind"], {"members": 4, "edges": 6, "clusters": 4000, "templates": 4})
        # The RMC clusters are more rigid than independent atoms would make them.
        self.assertLess(result["z"]["ladderWidth"], -10)
        self.assertLess(result["observed"]["ladderWidth"], 0.9 * result["nullMean"]["ladderWidth"])

    def test_the_demo_verdict(self):
        # As the Copilot's analyze_clusters reports it (its model-visible
        # snapshot): 20 % less distorted at the clusters' own size, yet more
        # clusters elongated and compressed along <111>; breathing; the
        # average cluster statically distorted (5.5 SE); 8 % of the bonds
        # outside the window; no split member site.
        verdict = cluster_distortion_from_file(str(DEMO_RMC6F), "Ta", 2.65, 3.33)["verdict"]
        answers = _answers(verdict)
        self.assertEqual(answers["distorted"], (True, "axis_family"))
        self.assertEqual(answers["lessDistorted"], (True, "below_null"))
        self.assertEqual(answers["breathing"], (True, "above_null"))
        self.assertEqual(answers["rigid"], (False, "pattern"))
        self.assertEqual(answers["staticDistortion"], (True, "above_sampling_error"))
        self.assertEqual(answers["splitMemberSites"], (False, "none_split"))
        self.assertEqual(answers["trigonalLean"], (False, "within_null"))
        self.assertEqual(answers["window.bondsOutsideWindow"], (True, "above_bar"))
        self.assertEqual(
            [(f["axis"], f["family"]) for f in verdict["favouredAxes"]["families"] if f["answer"]],
            [("long", "111"), ("short", "111")],
        )
        self.assertEqual(len(verdict["splitMemberSites"]["memberSites"]), 16)
        self.assertEqual(AXIS_ENHANCEMENT - 1.0, verdict["favouredAxes"]["evidence"][0]["bar"]["change"])


if __name__ == "__main__":
    unittest.main()
