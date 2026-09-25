# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The committed browser-parity golden must still be what kde.py computes.

``web_app/frontend/src/__tests__/fixtures/kde_parity_fixture.json`` pins the
browser KDE worker to the Python reference (``kdeParity.test.js``). This test
recomputes every case whose input is available and fails when the Python engine
has drifted from the golden, i.e. when ``tests/generate_kde_fixture.py`` needs
to be re-run (and the JS port re-checked against the new numbers).
"""

from __future__ import annotations

import json
import unittest

import numpy as np

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_kde_fixture import OUT, case_positions, compute_case
except ModuleNotFoundError:  # `python -m unittest tests.test_kde_parity_fixture`
    from tests.generate_kde_fixture import OUT, case_positions, compute_case


class KdeParityFixtureTests(unittest.TestCase):
    def test_fixture_matches_the_python_engine(self):
        fixture = json.loads(OUT.read_text(encoding="utf-8"))
        checked = 0
        for case in fixture["cases"]:
            positions = case_positions(case)
            if positions is None:
                self.assertTrue(case["requiresData"], f"{case['name']}: committed input is missing")
                continue
            with self.subTest(case=case["name"]):
                expected = case["expected"]
                actual = compute_case(case, positions)
                self.assertEqual(actual["slabCount"], expected["slabCount"])
                self.assertEqual(actual["fitCount"], expected["fitCount"])
                self.assertEqual(actual["message"], expected["message"])
                self.assertEqual(actual["warnings"], expected["warnings"])
                if expected["kernel"] is None:
                    self.assertIsNone(actual["kernel"])
                else:
                    np.testing.assert_allclose(
                        actual["kernel"]["covariance"], expected["kernel"]["covariance"], rtol=1e-12
                    )
                    # The peak is SciPy's compiled sum, which whitens through
                    # chol(H) (>= 1.10) or chol(H^-1) (older) with release-
                    # dependent round-off, amplified by the kernel's condition
                    # number: near-collinear (cond 6.6e6) moves by 6e-11 on
                    # SciPy 1.17 and 1.7e-9 on 1.8 against the 1.13 golden;
                    # the well-conditioned cases by <= 1e-11.
                    kernel = expected["kernel"]
                    condition = (kernel["sigmaMajor"] / kernel["sigmaMinor"]) ** 2
                    rtol = 1e-10 + 20 * condition * np.finfo(float).eps
                    np.testing.assert_allclose(actual["vmax"], expected["vmax"], rtol=rtol)
                np.testing.assert_allclose(
                    actual["densityOverPeak"], expected["densityOverPeak"], rtol=0, atol=1e-8
                )
                checked += 1
        self.assertGreater(checked, 0)


if __name__ == "__main__":
    unittest.main()
