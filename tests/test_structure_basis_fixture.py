# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The committed structure-basis golden must still be what structure_basis.py computes.

``web_app/frontend/src/__tests__/fixtures/structure_basis_fixture.json`` holds
the ``basis`` rows that ``/api/structure`` returns and that
``structureFromRmc6f()`` (``browserData.js``) is held to
(``src/__tests__/structureBasisParity.test.js``); the frontend's Flask-mode tests
also use the demo run's rows as the server's payload. This test re-runs the
generator in memory on the committed configurations and fails when the engine
has drifted from the file, i.e. when
``PYTHONPATH=. python tests/generate_structure_basis_fixture.py`` needs to be
re-run (and the JS parity re-checked against the new numbers). The comparison
is tests/golden_compare.py, the rule the other goldens use.
"""

from __future__ import annotations

import json
import unittest

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_structure_basis_fixture import OUT, build_fixture
    from golden_compare import assert_same
except ModuleNotFoundError:  # `python -m unittest tests.test_structure_basis_fixture`
    from tests.generate_structure_basis_fixture import OUT, build_fixture
    from tests.golden_compare import assert_same

STALE = (
    "structure_basis_fixture.json is stale: run "
    "`PYTHONPATH=. python tests/generate_structure_basis_fixture.py`"
)


class StructureBasisFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.committed = json.loads(OUT.read_text(encoding="utf-8"))
        regenerated = build_fixture(cls.committed["synthetic"], cls.committed["coordsOnly"])
        cls.regenerated = json.loads(json.dumps(regenerated, allow_nan=False))

    def test_the_cap_is_the_engines(self):
        # The JS test holds symmetryModel.js MAX_SYMMETRY_SITES to this value.
        self.assertEqual(self.regenerated["maxBasisSites"], self.committed["maxBasisSites"], STALE)

    def test_basis_rows(self):
        got_cases = self.regenerated["cases"]
        want_cases = self.committed["cases"]
        self.assertEqual([case["name"] for case in got_cases], [case["name"] for case in want_cases], STALE)
        for got, want in zip(got_cases, want_cases):
            with self.subTest(case=want["name"]):
                self.assertEqual(got["source"], want["source"], STALE)
                assert_same(self, got["basis"], want["basis"], "basis", message=STALE)

    def test_the_fixture_covers_every_kind_of_site(self):
        demo, synthetic, coords_only = (case["basis"] for case in self.committed["cases"])
        self.assertEqual(len(demo), 52)
        self.assertEqual(coords_only, [])
        by_ref = {site["referenceNumber"]: site for site in synthetic}
        # Reference numbers out of file order come back sorted.
        self.assertEqual(list(by_ref), sorted(by_ref))
        # A site whose copies straddle the cell edge, its centre clear of the
        # edge itself (so neither runtime can land on the other side of it).
        edge = by_ref[3]
        self.assertTrue(0.0 < edge["frac"][0] < 0.01 and 0.99 < edge["frac"][2] < 1.0)
        self.assertTrue(min(edge["frac"][0], 1.0 - edge["frac"][2]) > 1e-6)
        # A mixed site named by its first atom, the minority species.
        mixed = by_ref[1]
        self.assertEqual(mixed["el"], "Ta")
        self.assertEqual(mixed["elementCounts"], {"Ta": 12, "Nb": 20})
        # A split site whose arithmetic mean leaves the circular-mean centre.
        self.assertGreater(abs(by_ref[7]["mean"][0] - by_ref[7]["frac"][0]), 1e-5)
        # A correlated site (off-diagonal covariance), a still site and a one-copy site.
        self.assertGreater(abs(by_ref[4]["covFrac"][0][1]), 1e-5)
        self.assertEqual((by_ref[12]["dispA"], by_ref[12]["count"]), (0.0, 32))
        self.assertEqual((by_ref[9]["dispA"], by_ref[9]["count"]), (0.0, 1))
        # The three coords-only atoms are not in any site.
        self.assertEqual(sum(site["count"] for site in synthetic) + 3, int(self.committed["synthetic"].split("Number of atoms:")[1].split()[0]))


if __name__ == "__main__":
    unittest.main()
