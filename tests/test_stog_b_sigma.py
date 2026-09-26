# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The sigma-column guard shared by the CLI, the API and (ported) the page (stog-b).

``usable_sigma`` drops the whole uncertainty column when any row with finite Q and
S has a zero, negative or non-finite sigma; the JS port ``usableSigma`` now guards
the browser path the same way (it used to feed a 1e12 weight or NaN to the fit).
"""

import unittest

import numpy as np

from rmc_toolkits.scaling_cli import usable_sigma


class UsableSigmaTests(unittest.TestCase):
    def setUp(self):
        self.q = np.linspace(0.5, 30.0, 100)
        self.sq = np.ones_like(self.q)
        self.sigma = np.full_like(self.q, 1e-3)

    def test_clean_column_is_kept(self):
        self.assertIs(usable_sigma(self.q, self.sq, self.sigma), self.sigma)
        self.assertIsNone(usable_sigma(self.q, self.sq, None))

    def test_any_bad_value_on_a_usable_row_drops_the_column(self):
        for bad in (0.0, -1e-3, np.nan, np.inf):
            sigma = self.sigma.copy()
            sigma[90] = bad
            with self.subTest(bad=bad):
                self.assertIsNone(usable_sigma(self.q, self.sq, sigma))

    def test_rows_without_data_do_not_count(self):
        sq = self.sq.copy()
        sq[3] = np.nan
        sigma = self.sigma.copy()
        sigma[3] = np.nan
        self.assertIs(usable_sigma(self.q, sq, sigma), sigma)


if __name__ == "__main__":
    unittest.main()
