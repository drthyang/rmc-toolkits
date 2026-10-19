# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The committed displacement-direction golden must still be what orientation.py computes.

``web_app/frontend/src/__tests__/fixtures/orientation_fixture.json`` holds the
Python engine's payloads that ``workers/orientation.js`` is held to
(``src/__tests__/orientationParity.test.js``). This test re-runs the generator
in memory on the committed inputs and fails when the engine has drifted from
the file, i.e. when ``PYTHONPATH=. python tests/generate_orientation_fixture.py``
needs to be re-run (and the JS port re-checked against the new numbers). The
comparison is tests/golden_compare.py, the same rule the JS test applies.
"""

from __future__ import annotations

import json
import unittest

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_orientation_fixture import OUT, build_fixture
    from golden_compare import assert_same
except ModuleNotFoundError:  # `python -m unittest tests.test_orientation_fixture`
    from tests.generate_orientation_fixture import OUT, build_fixture
    from tests.golden_compare import assert_same

STALE = "orientation_fixture.json is stale: run `PYTHONPATH=. python tests/generate_orientation_fixture.py`"
RELATIVE = frozenset({"peakLocalPValue", "peakPValue", "mapPValue", "orientationBinghamPValue"})


class OrientationFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.committed = json.loads(OUT.read_text(encoding="utf-8"))
        regenerated = build_fixture(cls.committed["inputs"])
        cls.regenerated = json.loads(json.dumps(regenerated, allow_nan=False))

    def test_sections(self):
        self.assertEqual(sorted(self.regenerated), sorted(self.committed), STALE)

    def test_histogram_payloads(self):
        got_cases = self.regenerated["histograms"]
        want_cases = self.committed["histograms"]
        self.assertEqual([case["name"] for case in got_cases], [case["name"] for case in want_cases], STALE)
        for got, want in zip(got_cases, want_cases):
            with self.subTest(case=want["name"]):
                for field in ("input", "options", "cells", "keys"):
                    self.assertEqual(got[field], want[field], f"{STALE} ({field})")
                assert_same(self, got["spectra"], want["spectra"], "spectra", message=STALE)
                assert_same(
                    self, got["expected"], want["expected"], "expected",
                    relative=RELATIVE, bases=want["spectra"], message=STALE,
                )

    def test_tilings_and_tables(self):
        for section in ("tilings", "assignments", "tilingChecks", "recommendedFrequency", "pearsonMoments"):
            with self.subTest(section=section):
                assert_same(self, self.regenerated[section], self.committed[section], section, message=STALE)

    def test_the_fixture_covers_every_kind_of_case(self):
        cases = [case["expected"] for case in self.committed["histograms"]]
        self.assertTrue(any(case["mapSignificance"] is None for case in cases))
        self.assertTrue(any(case["antipodalAsymmetryZ"] is None for case in cases))
        self.assertTrue(any(case["antipodalAsymmetrySignificant"] for case in cases))
        self.assertTrue(any("referenceNumber" in case for case in cases))
        self.assertTrue(any("polygons" in case for case in cases))
        self.assertEqual({case["weight"] for case in cases}, {"count", "amplitude", "amplitude2"})
        self.assertEqual({case["frame"] for case in cases}, {"cartesian", "pca"})
        self.assertGreaterEqual(len({case["frequency"] for case in cases}), 3)


if __name__ == "__main__":
    unittest.main()
