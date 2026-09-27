# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Lorch omitted-low-Q basis without catastrophic cancellation (0.6.0 audit, stog-b).

The Lorch branch evaluated (v sin v + cos v - 1)/(r - a)^2, v = Q0 (r - a), with
cos v - 1 by subtraction: just outside the 1e-9 patch around the removable
singularity r = a = pi/Qmax the coefficient lost every digit (80x off at
r - a = 2e-9, sign-flipped at r = 0.11 for Qmax = 28.56 on the default grid).
The basis is now written with sinc(v) and cos v - 1 = -2 sin^2(v/2), which needs
no patch; this test sweeps the whole band |r - a| in 1e-9 .. 1e-5 (and r = a,
r = 0) against adaptive quadrature of the defining integrals.
"""

import unittest

import numpy as np
from scipy.integrate import quad

from rmc_toolkits.transforms import low_q_correction_basis


def reference(q0, qmax, r, s0_target=0.0):
    """coef/const of the Lorch-windowed linear extrapolation, by quadrature."""
    a = np.pi / qmax

    def window(q):
        return np.sinc(a * q / np.pi)  # sin(aQ)/(aQ)

    opts = dict(epsabs=0.0, epsrel=1e-13, limit=200)
    coef = (2 / np.pi) * quad(lambda q: q * (q / q0) * window(q) * np.sin(q * r), 0, q0, **opts)[0]
    const = (2 / np.pi) * quad(lambda q: q * window(q) * np.sin(q * r), 0, q0, **opts)[0]
    if s0_target:
        const = (1 - s0_target) * const + s0_target * coef
    return coef, const


class LorchBasisTests(unittest.TestCase):
    def check(self, q0, qmax, r_values, s0_target=0.0):
        q = np.linspace(q0, qmax, 400)
        coef, const = low_q_correction_basis(q, r_values, lorch=True, s0_target=s0_target)
        for index, r in enumerate(r_values):
            ref_coef, ref_const = reference(q0, qmax, r, s0_target)
            scale = max(abs(ref_coef), 1e-300)
            with self.subTest(q0=q0, qmax=qmax, r=r):
                self.assertLess(abs(coef[index] - ref_coef) / scale, 1e-11)
                self.assertLess(abs(const[index] - ref_const), 1e-11 * max(abs(ref_const), 1e-3))

    def test_band_around_the_removable_singularity(self):
        for q0, qmax in ((1.0, 28.0), (0.5, 28.56), (1.0, 52.36), (0.5, 26.18)):
            a = np.pi / qmax
            deltas = np.concatenate([[0.0], np.logspace(-9, -5, 9)])
            r_values = np.concatenate([a - deltas[1:], a + deltas])
            self.check(q0, qmax, r_values)

    def test_default_grid_points_near_pi_over_qmax(self):
        # Realistic Qmax values put a 0.01 A grid point 1e-7..3e-7 A from pi/Qmax.
        r_grid = np.arange(1, 31) * 0.01
        for qmax in (22.44, 26.18, 28.56, 39.27, 44.88, 52.36):
            self.check(0.5, qmax, r_grid)

    def test_r_zero_and_s0_target(self):
        self.check(1.0, 28.0, np.array([0.0, 0.05, np.pi / 28.0, 0.5, 2.0]), s0_target=-12.06)


if __name__ == "__main__":
    unittest.main()
