# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""A legacy coords-only ``.rmc6f`` reads the same in Flask/CLI as in the browser.

``web_app/frontend/src/__tests__/fixtures/rmc6f_coords_only_fixture.json``
holds a coords-only configuration (``id element [label] x y z``, no site or
cell columns) and what the Python FILE loaders make of it; the browser side is
pinned to the same golden by ``coordsOnlyParity.test.js``. This test recomputes
the Python side from the fixture's own text and fails when it drifts, i.e. when
``tests/generate_coords_only_fixture.py`` needs to be re-run (and the JS port
re-checked). The KDE loader used to yield zero positions for such a file, the
Flask ``/api/kde/slice`` therefore declined an empty slab where the browser
folded every atom.
"""

from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.kde import load_unit_cell_positions
from rmc_toolkits.parsers import Rmc6fParseReport, iter_rmc6f_atoms
from rmc_toolkits.triplets import bond_angle_summary_from_file

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_coords_only_fixture import OUT
except ModuleNotFoundError:  # `python -m unittest tests.test_coords_only_fixture`
    from tests.generate_coords_only_fixture import OUT


class CoordsOnlyFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fixture = json.loads(OUT.read_text(encoding="utf-8"))
        cls._tmp = tempfile.TemporaryDirectory()
        cls.path = Path(cls._tmp.name) / "coords_only.rmc6f"
        cls.path.write_text(cls.fixture["rmc6f"], encoding="utf-8")

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def test_the_fixture_really_is_coords_only(self):
        report = Rmc6fParseReport()
        full_layout = list(iter_rmc6f_atoms(self.path, report=report))
        self.assertEqual(full_layout, [])
        self.assertEqual(report.coords_only_atoms, self.fixture["kdePositions"]["count"])
        self.assertEqual(report.declared_atoms, report.coords_only_atoms)
        self.assertIsNone(report.warning())

    def test_kde_positions_match_the_golden(self):
        positions = load_unit_cell_positions(self.path).fractional_positions
        expected = np.asarray(self.fixture["kdePositions"]["fractional"], dtype=float)
        self.assertGreater(len(expected), 0)
        self.assertEqual(positions.shape, expected.shape)
        np.testing.assert_allclose(positions, expected, rtol=0, atol=1e-15)

    def test_bond_angle_summaries_match_the_golden(self):
        for case in self.fixture["triplets"]:
            spec, expected = case["request"], case["summary"]
            with self.subTest(triplet=(spec["end1"], spec["apex"], spec["end2"])):
                actual = bond_angle_summary_from_file(
                    self.path,
                    spec["end1"], spec["apex"], spec["end2"],
                    spec["r12Min"], spec["r12Max"],
                    spec.get("r23Min", spec["r12Min"]), spec.get("r23Max", spec["r12Max"]),
                    spec["binWidth"],
                )
                self.assertGreater(actual["angleCount"], 0)
                self.assertEqual(actual["angleCount"], expected["angleCount"])
                self.assertEqual(actual["counts"], expected["counts"])
                self.assertEqual(actual["coordination"], expected["coordination"])
                self.assertEqual(actual["lengths12"]["counts"], expected["lengths12"]["counts"])
                self.assertIsNone(actual["parseWarning"])
                self.assertAlmostEqual(actual["meanAngle"], expected["meanAngle"], places=9)


if __name__ == "__main__":
    unittest.main()
