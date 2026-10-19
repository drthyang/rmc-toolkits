# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The committed site-statistics golden must still be what pca_kde.py computes.

``web_app/frontend/src/__tests__/fixtures/pca_sites_fixture.json`` holds the
``site_ellipsoids`` rows that ``workers/pcaKde.js`` is held to
(``src/__tests__/pcaSitesParity.test.js``). This test re-runs the generator in
memory on the committed configurations and fails when the engine has drifted
from the file, i.e. when ``PYTHONPATH=. python tests/generate_pca_sites_fixture.py``
needs to be re-run (and the JS port re-checked against the new numbers). The
comparison is tests/golden_compare.py, the same rule the JS test applies.
"""

from __future__ import annotations

import json
import unittest

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_pca_sites_fixture import OUT, build_fixture
    from golden_compare import assert_same
except ModuleNotFoundError:  # `python -m unittest tests.test_pca_sites_fixture`
    from tests.generate_pca_sites_fixture import OUT, build_fixture
    from tests.golden_compare import assert_same

STALE = "pca_sites_fixture.json is stale: run `PYTHONPATH=. python tests/generate_pca_sites_fixture.py`"


class PcaSitesFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.committed = json.loads(OUT.read_text(encoding="utf-8"))
        regenerated = build_fixture(cls.committed["synthetic"])
        cls.regenerated = json.loads(json.dumps(regenerated, allow_nan=False))

    def test_site_rows(self):
        got_cases = self.regenerated["cases"]
        want_cases = self.committed["cases"]
        self.assertEqual([case["name"] for case in got_cases], [case["name"] for case in want_cases], STALE)
        for got, want in zip(got_cases, want_cases):
            with self.subTest(case=want["name"]):
                for field in ("source", "probability", "referenceNumbers", "species", "parseWarning"):
                    self.assertEqual(got[field], want[field], f"{STALE} ({field})")
                self.assertEqual(len(got["expected"]), len(want["expected"]), STALE)
                for got_row, want_row in zip(got["expected"], want["expected"]):
                    assert_same(
                        self, got_row, want_row, f"site {want_row['referenceNumber']}",
                        bases={"axes": want_row["eigenvalues"]}, message=STALE,
                    )

    def test_the_fixture_covers_every_kind_of_site(self):
        rows = [row for case in self.committed["cases"] for row in case["expected"]]
        self.assertTrue(any(row["zeroSpread"] and row["axes"] is None for row in rows))
        self.assertTrue(any(row["mixed"] for row in rows))
        self.assertTrue(any((row["nonGaussianity"] or 0.0) < 0.0 for row in rows))
        self.assertTrue(any(all(row["axisResolved"]) for row in rows))
        self.assertTrue(any(not row["zeroSpread"] and not any(row["axisResolved"]) for row in rows))
        self.assertLess(max(row["count"] for row in rows), 20000)


if __name__ == "__main__":
    unittest.main()
