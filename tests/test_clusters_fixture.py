# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The committed cluster golden must still be what clusters.py computes.

``web_app/frontend/src/__tests__/fixtures/clusters_fixture.json`` holds the
``cluster_distortion`` payloads, verdicts included, that ``workers/clusters.js``
is held to (``src/workers/__tests__/clusters.test.js``). This test re-runs the
generator in memory on the committed configurations and fails when the engine
has drifted from the file, i.e. when
``PYTHONPATH=. python tests/generate_clusters_fixture.py`` needs to be re-run
(and the JS port re-checked against the new numbers). The comparison is
tests/golden_compare.py, the same rule the JS test applies.
"""

from __future__ import annotations

import json
import unittest

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_clusters_fixture import CASES, OUT, build_fixture
    from golden_compare import assert_same
except ModuleNotFoundError:  # `python -m unittest tests.test_clusters_fixture`
    from tests.generate_clusters_fixture import CASES, OUT, build_fixture
    from tests.golden_compare import assert_same

STALE = "clusters_fixture.json is stale: run `PYTHONPATH=. python tests/generate_clusters_fixture.py`"


class ClustersFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.committed = json.loads(OUT.read_text(encoding="utf-8"))
        regenerated = build_fixture({case["name"]: case["text"] for case in cls.committed["cases"]})
        cls.regenerated = json.loads(json.dumps(regenerated, allow_nan=False))

    def test_cases(self):
        got_cases = self.regenerated["cases"]
        want_cases = self.committed["cases"]
        self.assertEqual([case["name"] for case in got_cases], [case["name"] for case in want_cases], STALE)
        for got, want in zip(got_cases, want_cases):
            with self.subTest(case=want["name"]):
                self.assertEqual(got["queries"], want["queries"], f"{STALE} (queries)")
                assert_same(self, got["expected"], want["expected"], "expected", message=STALE)

    def test_hashes_and_enums(self):
        assert_same(self, self.regenerated["hashes"], self.committed["hashes"], "hashes", message=STALE)
        self.assertEqual(self.regenerated["enums"], self.committed["enums"], STALE)

    def test_the_fixture_covers_every_kind_of_verdict(self):
        # Each answer the renderer says something for appears at least once,
        # so the JS port is held to it.
        answers = set()
        for case in self.committed["cases"]:
            for result in case["expected"]:
                verdict = result["verdict"]
                for name, finding in verdict.items():
                    if name == "rules":
                        continue
                    items = finding.items() if name == "window" else [(name, finding)]
                    answers.update((key, item["because"]) for key, item in items)
        for wanted in [
            ("distorted", "edge_spread"), ("distorted", "axis_family"), ("distorted", "less_distorted"),
            ("distorted", "within_null"), ("distorted", "extended_network"), ("distorted", "no_clusters"),
            ("lessDistorted", "below_null"), ("lessDistorted", "single_edge"),
            ("breathing", "above_null"), ("breathing", "below_null"),
            ("rigid", "less_distorted"), ("rigid", "stiffer_bonds"), ("rigid", "split_sites"),
            ("staticDistortion", "above_sampling_error"), ("splitMemberSites", "half_or_more_split"),
            ("trigonalLean", "leans_one_edge"), ("trigonalLean", "not_tetrahedra"),
            ("favouredAxes", "enhanced"), ("extendedNetwork", "extended_network"),
            ("mixedSizes", "above_bar"),
        ]:
            self.assertIn(wanted, answers)
        self.assertEqual([name for name, _, _ in CASES], [case["name"] for case in self.committed["cases"]])


if __name__ == "__main__":
    unittest.main()
