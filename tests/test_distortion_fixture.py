# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The committed polyhedral-distortion golden must still be what distortion.py computes.

``web_app/frontend/src/__tests__/fixtures/distortion_fixture.json`` holds the
``polyhedral_distortion`` payloads (and refusals) that ``workers/distortion.js``
is held to (``src/workers/__tests__/distortion.test.js``). This test re-runs the
generator in memory on the committed configurations and fails when the engine
has drifted from the file, i.e. when
``PYTHONPATH=. python tests/generate_distortion_fixture.py`` needs to be re-run
(and the JS port re-checked against the new numbers). The comparison is
tests/golden_compare.py, the same rule the JS test applies.
"""

from __future__ import annotations

import json
import unittest

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_distortion_fixture import CASES, OUT, build_fixture, constants
    from golden_compare import assert_same
except ModuleNotFoundError:  # `python -m unittest tests.test_distortion_fixture`
    from tests.generate_distortion_fixture import CASES, OUT, build_fixture, constants
    from tests.golden_compare import assert_same

STALE = "distortion_fixture.json is stale: run `PYTHONPATH=. python tests/generate_distortion_fixture.py`"


class DistortionFixtureTests(unittest.TestCase):
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

    def test_constants(self):
        assert_same(self, json.loads(json.dumps(constants())), self.committed["constants"], "constants", message=STALE)

    def test_the_fixture_covers_every_path(self):
        results = [result for case in self.committed["cases"] for result in case["expected"]]
        shapes = {result["polyhedron"]["shape"] for result in results if "error" not in result}
        self.assertEqual(shapes, {"octahedron", "tetrahedron", None})
        reasons = {result["average"]["reason"] for result in results if "error" not in result}
        self.assertEqual(reasons, {None, "no_basis", "no_centre_sites", "coordination_differs"})
        errors = " | ".join(result["error"] for result in results if "error" in result)
        for phrase in ("no polyhedra: no", "is one Cu", "over the limit of 100", "unknown element", "window must be"):
            self.assertIn(phrase, errors)
        # A flat-ish polyhedron, a mixed site and a box smaller than the window.
        self.assertTrue(any(
            "error" not in result and (result["configuration"]["quadraticElongation"] or {}).get("mean", 0) > 2
            for result in results
        ))
        self.assertTrue(any(site["mixed"] for result in results if "error" not in result for site in result["average"]["sites"]))
        self.assertEqual([name for name, _, _ in CASES], [case["name"] for case in self.committed["cases"]])


if __name__ == "__main__":
    unittest.main()
