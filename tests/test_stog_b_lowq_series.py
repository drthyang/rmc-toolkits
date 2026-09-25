# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Unwindowed omitted-low-Q basis at small v = Q0 r (stog-b, found verifying r = 0).

The closed forms (2v sin v - (v^2 - 2) cos v - 2)/r^3 and (sin v - v cos v)/r^2
cancel O(1) terms down to O(v^3)/O(v^4): for v below ~1e-3 (Q0 = 0.01 on the
0.01 A grid, or fine r grids) the coefficient was wrong by 100 % or more. Below
|v| = 0.5 the Taylor series is used; checked against quadrature.
"""

import unittest

import numpy as np
from scipy.integrate import quad

from rmc_toolkits.transforms import low_q_correction_basis


class UnwindowedSeriesTests(unittest.TestCase):
    def test_moments_match_quadrature_at_every_v(self):
        for q0 in (0.01, 0.5, 1.0, 2.0):
            q = np.linspace(q0, 28.0, 300)
            r = np.array([1e-6, 1e-4, 1e-3, 0.002, 0.01, 0.1, 0.49 / q0, 0.51 / q0, 1.0, 5.0])
            coef, const = low_q_correction_basis(q, r)
            for index, radius in enumerate(r):
                ref_coef = (2 / np.pi) * quad(
                    lambda x: x * (x / q0) * np.sin(x * radius), 0, q0, epsabs=0, epsrel=1e-13,
                )[0]
                ref_const = (2 / np.pi) * quad(
                    lambda x: x * np.sin(x * radius), 0, q0, epsabs=0, epsrel=1e-13,
                )[0]
                with self.subTest(q0=q0, r=radius):
                    self.assertLess(abs(coef[index] - ref_coef) / abs(ref_coef), 1e-12)
                    self.assertLess(abs(const[index] - ref_const) / abs(ref_const), 1e-12)

    def test_r_zero_stays_zero(self):
        coef, const = low_q_correction_basis(np.linspace(0.5, 28.0, 100), np.array([0.0]))
        self.assertEqual(coef[0], 0.0)
        self.assertEqual(const[0], 0.0)


if __name__ == "__main__":
    unittest.main()
