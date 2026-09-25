# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""fourier_filter on an r grid that starts at r = 0 (1.0 audit, stog-b group).

The filter divided G_PDF by 4 pi rho0 r to get g and multiplied back for the
section integrand: at r = 0 that is 0/0 = NaN, and the NaN reached every output Q
through the back transform — 100 % NaN sq_filtered / sq_ft / g_filtered for the
common np.linspace(0, rmax, n) grid. The integrand is now G_PDF + 4 pi rho0 r (the
same number, no division) and g_filtered(0) is the continuous extension
1 + G_PDF'(0)/(4 pi rho0) (gpdf_slope_at_zero).
"""

import unittest

import numpy as np

from rmc_toolkits.transforms import (
    fourier_filter,
    fq_to_gpdf,
    fq_to_sq,
    g_to_gpdf,
    gpdf_slope_at_zero,
    gpdf_to_fq,
    sq_to_fq,
)

RHO0 = 0.05
OPTIONS = [
    dict(lorch=lorch, low_q_correction=lqc, s0_target=s0)
    for lorch in (False, True) for lqc in (False, True) for s0 in (0.0, -12.06)
]


def model():
    q = np.arange(20, 981) * 0.03
    r = np.arange(1, 12001) * 0.005
    g = 0.5 * (1 + np.tanh((r - 2.65) / 0.07)) + 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    return q, fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, RHO0), q))


class RZeroTests(unittest.TestCase):
    def test_grid_from_zero_gives_finite_outputs(self):
        q, sq = model()
        r0 = np.linspace(0.0, 25.0, 2501)
        for options in OPTIONS:
            with self.subTest(**options):
                sq_filtered, sq_ft, g_filtered = fourier_filter(
                    q, sq, r0, rho0=RHO0, cutoff=1.0, **options
                )
                for array in (sq_filtered, sq_ft, g_filtered):
                    self.assertTrue(np.all(np.isfinite(array)))
                # g(0) continues g(r): even in r, so (4 g(dr) - g(2 dr))/3 ~ g(0).
                extrapolated = (4.0 * g_filtered[1] - g_filtered[2]) / 3.0
                self.assertLess(abs(g_filtered[0] - extrapolated), 1e-4)
                # r = 0 adds only the [0, dr] panel to the section integral.
                ref = fourier_filter(q, sq, r0[1:], rho0=RHO0, cutoff=1.0, **options)
                self.assertLess(np.max(np.abs(sq_filtered - ref[0])), 1e-6)

    def test_slope_at_zero_matches_richardson(self):
        q, sq = model()
        fq = sq_to_fq(q, sq)
        h = 0.002
        for options in OPTIONS:
            if options["low_q_correction"] and not options["lorch"]:
                continue  # its closed form cancels at r ~ 1e-3; the moments are exact
            with self.subTest(**options):
                g = fq_to_gpdf(q, fq, np.array([h, 2 * h, 4 * h]), **options)
                s1 = (8 * g[0] - g[1]) / (6 * h)
                s2 = (8 * g[1] - g[2]) / (12 * h)
                richardson = (16 * s1 - s2) / 15
                slope = gpdf_slope_at_zero(q, fq, **options)
                self.assertLess(abs(slope - richardson) / abs(richardson), 1e-9)

    def test_negative_r_is_rejected(self):
        q, sq = model()
        with self.assertRaisesRegex(ValueError, "non-negative"):
            fourier_filter(q, sq, np.linspace(-1.0, 5.0, 61), rho0=RHO0, cutoff=1.0)


if __name__ == "__main__":
    unittest.main()
