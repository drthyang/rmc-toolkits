# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The committed G(r)-ribbon golden must still be what ft_uncertainty.py computes.

``web_app/frontend/src/__tests__/fixtures/ft_uncertainty_fixture.json`` holds
run folders and the ``uncertainty`` payload ``ft_uncertainty_for_run`` makes of
each (what Flask's ``/api/plot/data`` adds to an xPDF plot); the browser's port
is held to it by ``src/__tests__/ftUncertaintyParity.test.js``. This test
re-runs the engine on the committed files and fails when it has drifted, i.e.
when ``PYTHONPATH=. python tests/generate_ft_uncertainty_fixture.py`` needs to
be re-run (and the JS parity re-checked). The comparison is
tests/golden_compare.py, the rule the other goldens use.
"""

from __future__ import annotations

import json
import unittest

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_ft_uncertainty_fixture import OUT, build_fixture, committed_cases
    from golden_compare import assert_same
except ModuleNotFoundError:  # `python -m unittest tests.test_ft_uncertainty_fixture`
    from tests.generate_ft_uncertainty_fixture import OUT, build_fixture, committed_cases
    from tests.golden_compare import assert_same

STALE = (
    "ft_uncertainty_fixture.json is stale: run "
    "`PYTHONPATH=. python tests/generate_ft_uncertainty_fixture.py`"
)


class FtUncertaintyFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.committed = json.loads(OUT.read_text(encoding="utf-8"))
        regenerated = build_fixture(committed_cases(cls.committed))
        cls.regenerated = json.loads(json.dumps(regenerated, allow_nan=False))

    def test_the_payloads(self):
        assert_same(self, self.regenerated["demo"], self.committed["demo"], "demo", message=STALE)
        got_cases = self.regenerated["cases"]
        want_cases = self.committed["cases"]
        self.assertEqual([case["name"] for case in got_cases], [case["name"] for case in want_cases], STALE)
        for got, want in zip(got_cases, want_cases):
            with self.subTest(case=want["name"]):
                assert_same(self, got["expected"], want["expected"], "expected", message=STALE)

    def test_the_fixture_covers_every_outcome(self):
        statuses = {case["name"]: case["expected"]["status"] for case in self.committed["cases"]}
        self.assertEqual(sorted(set(statuses.values())), ["failed", "propagated", "unavailable"])
        self.assertEqual(self.committed["demo"]["expected"]["status"], "unavailable")
        self.assertTrue(self.committed["demo"]["expected"]["reproduction"]["gr"]["pass"])
        propagated = [case for case in self.committed["cases"] if case["expected"]["status"] == "propagated"]
        self.assertTrue(any(case["expected"]["convolutionCheck"] for case in propagated))
        self.assertTrue(any(case["expected"]["convolutionCheck"] is None for case in propagated))
        reasons = [case["expected"]["reason"] for case in self.committed["cases"] if case["expected"]["status"] == "failed"]
        self.assertEqual(len(reasons), 3)


if __name__ == "__main__":
    unittest.main()
