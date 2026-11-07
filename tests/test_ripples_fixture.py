# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The committed termination-ripple golden must still be what ripples.py computes.

``web_app/frontend/src/__tests__/fixtures/ripples_fixture.json`` holds run
folders, the calls made on them and what ``ripples_for_run`` returns (what
Flask's ``/api/ripples`` sends); the browser's port is held to it by
``src/__tests__/ripplesParity.test.js``. This test re-runs the engine on the
committed files and fails when it has drifted, i.e. when
``PYTHONPATH=. python tests/generate_ripples_fixture.py`` needs to be re-run
(and the JS parity re-checked). The comparison is tests/golden_compare.py at
the ripple golden's tolerance (RIPPLES_RTOL, the port's: an iterative fit).
"""

from __future__ import annotations

import json
import unittest

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_ripples_fixture import OUT, build_fixture, committed_cases, near_thresholds
    from golden_compare import assert_same
except ModuleNotFoundError:  # `python -m unittest tests.test_ripples_fixture`
    from tests.generate_ripples_fixture import OUT, build_fixture, committed_cases, near_thresholds
    from tests.golden_compare import assert_same

RIPPLES_RTOL = 1e-7
STALE = "ripples_fixture.json is stale: run `PYTHONPATH=. python tests/generate_ripples_fixture.py`"


class RipplesFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.committed = json.loads(OUT.read_text(encoding="utf-8"))
        regenerated = build_fixture(committed_cases(cls.committed))
        cls.regenerated = json.loads(json.dumps(regenerated, allow_nan=False, default=float))

    def test_the_results(self):
        for got, want in zip(self.regenerated["demo"]["calls"], self.committed["demo"]["calls"]):
            assert_same(self, got["expected"], want["expected"], "demo", message=STALE, rtol=RIPPLES_RTOL)
        got_cases = self.regenerated["cases"]
        want_cases = self.committed["cases"]
        self.assertEqual([case["name"] for case in got_cases], [case["name"] for case in want_cases], STALE)
        for got, want in zip(got_cases, want_cases):
            with self.subTest(case=want["name"]):
                for got_call, want_call in zip(got["calls"], want["calls"]):
                    assert_same(self, got_call["expected"], want_call["expected"], "expected", message=STALE,
                                rtol=RIPPLES_RTOL)

    def test_no_reading_sits_on_a_threshold(self):
        self.assertEqual(near_thresholds(self.committed), [])

    def test_the_demo_is_a_matched_real_space_fit(self):
        first = self.committed["demo"]["calls"][0]["expected"]
        xpdf = next(d for d in first["datasets"] if d["kind"] == "xpdf")
        self.assertEqual((xpdf["status"], xpdf["fit"]["risk"], xpdf["transform"]["source"]), ("assessed", "matched", "run-fq"))
        self.assertTrue(first["found"])


if __name__ == "__main__":
    unittest.main()
