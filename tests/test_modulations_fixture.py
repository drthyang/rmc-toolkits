# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The committed modulation golden must still be what modulations.py computes.

``web_app/frontend/src/__tests__/fixtures/modulations_fixture.json`` holds the
``displacement_modulations`` payloads, verdicts included, that
``workers/modulations.js`` is held to (``src/workers/__tests__/modulations.test.js``).
This test re-runs the generator in memory on the committed configurations and
fails when the engine has drifted from the file, i.e. when
``PYTHONPATH=. python tests/generate_modulations_fixture.py`` needs to be
re-run. The comparison is tests/golden_compare.py.
"""

from __future__ import annotations

import json
import unittest

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_modulations_fixture import OUT, build_fixture
    from golden_compare import assert_same
except ModuleNotFoundError:  # `python -m unittest tests.test_modulations_fixture`
    from tests.generate_modulations_fixture import OUT, build_fixture
    from tests.golden_compare import assert_same

STALE = "modulations_fixture.json is stale: run `PYTHONPATH=. python tests/generate_modulations_fixture.py`"


class ModulationsFixtureTests(unittest.TestCase):
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
                self.assertEqual(got["elements"], want["elements"], f"{STALE} (elements)")
                assert_same(self, got["expected"], want["expected"], "expected", message=STALE)

    def test_enums(self):
        self.assertEqual(self.regenerated["enums"], self.committed["enums"], STALE)

    def test_the_fixture_covers_every_kind_of_verdict(self):
        seen = set()
        for case in self.committed["cases"]:
            for result in case["expected"]:
                for entry in result["elements"]:
                    seen.update((name, finding["answer"], finding["because"]) for name, finding in entry["verdict"].items())
                    seen.update(("point", point["kind"], point["because"])
                                for point in entry.get("zoneBoundary", []) + entry.get("strongest", []))
        for wanted in [
            ("cellDoubling", True, "ordered_point"), ("cellDoubling", False, "no_ordered_point"),
            ("generalOrder", True, "ordered_mode"), ("generalOrder", False, "no_ordered_mode"),
            ("longWave", True, "long_wave_dominates"), ("longWave", False, "strongest_not_long_wave"),
            ("point", "zone_boundary", "ordered"), ("point", "general", "ordered"),
            ("point", "general", "below_sigma"), ("point", "zone_boundary", "below_ratio"),
            ("point", "long_wave", "long_wave"), ("point", "zone_boundary", "zone_boundary"),
            ("diffuseOrder", True, "correlated_antiphase"), ("diffuseOrder", True, "correlated_in_phase"),
            ("diffuseOrder", False, "below_sigma"), ("diffuseOrder", None, "too_few_modes"),
        ]:
            self.assertIn(wanted, seen)


if __name__ == "__main__":
    unittest.main()
