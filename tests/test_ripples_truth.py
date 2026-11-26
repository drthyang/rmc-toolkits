# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Known-truth acceptance test for check_ripples' first principles (rmc_toolkits/ripples.py), synthetic data only.

What check_ripples ships besides RMCProfile's fit space (tests/test_ripples.py)
is the r-resolution of each dataset's transform and the rule its checks state:
ripples repeat every 2π/Q_max, a peak at r₀ has them near r₀ ± (n + ½)π/Q_max,
and the kernel's FWHM is 3.79/Q_max for a rectangle from ~0. Both are tested
here against numpy alone, independent of the engine's own helpers:

* ``resolution()`` on the rectangle Q = ΔQ … Q_max: ripplePeriod = 2π/Q_max and
  fwhm = 3.79/Q_max within 1 %, for Q_max 15, 20 and 25 Å⁻¹;
* a single Gaussian peak A[g(r − r₀) − g(r + r₀)] in D(r), its exact sine
  transform cut at Q_max and transformed back by the rectangle rule
  (2/π) Σ_k ΔQ F(Q_k) sin(Q_k r), Q_k = ΔQ … Q_max, ΔQ = 0.01; the true
  termination artefact is that less the untruncated peak. Its extrema lie
  within 0.05 Å of r₀ ± (n + ½)π/Q_max for n = 1 … 4, for σ 0.03, 0.07 and
  0.12 Å and Q_max 15, 20 and 30 Å⁻¹.

The verdict cases (groups A, B and C: whether a given feature is read as a
ripple) and their KNOWN_FAILURES live on the branch wip/joint-fit with the
strong-peak fit they test; that reading is not shipped.
"""

from __future__ import annotations

import math
import unittest

import numpy as np

from rmc_toolkits import ripples as R

DQ = 0.01
R0 = 5.0
R_FINE = R0 + np.round(np.arange(-1500, 1501) * 0.001, 10)  # r₀ ± 1.5 Å at 0.001 Å
POSITION_TOLERANCE = 0.05  # Å


def rectangle_q(qmax):
    return np.arange(1, int(round(qmax / DQ)) + 1) * DQ


def artefact(sigma, qmax):
    """The true termination artefact of one Gaussian peak at R0 (unit area): truncated less untruncated."""
    q = rectangle_q(qmax)
    f = np.exp(-q * q * sigma * sigma / 2) * np.sin(q * R0)
    truncated = (2 / math.pi) * DQ * (np.sin(np.outer(R_FINE, q)) @ f)
    c = 1.0 / (sigma * math.sqrt(2 * math.pi))
    untruncated = c * (np.exp(-(R_FINE - R0) ** 2 / (2 * sigma * sigma)) - np.exp(-(R_FINE + R0) ** 2 / (2 * sigma * sigma)))
    return truncated - untruncated


def extrema(y):
    i = np.arange(1, y.size - 1)
    return i[((y[i] > y[i - 1]) & (y[i] >= y[i + 1])) | ((y[i] < y[i - 1]) & (y[i] <= y[i + 1]))]


class Resolution(unittest.TestCase):
    def test_ripple_period_and_fwhm_of_a_rectangle_from_zero(self):
        for qmax in (15.0, 20.0, 25.0):
            with self.subTest(qmax=qmax):
                q = rectangle_q(qmax)
                res = R.resolution(q, R.transform_weights(q, DQ, "rectangle", "none"))
                self.assertAlmostEqual(res["ripplePeriod"], 2 * math.pi / qmax, places=12)
                self.assertLess(abs(res["fwhm"] / (3.79 / qmax) - 1), 0.01, res["fwhm"])


class RipplePositions(unittest.TestCase):
    def test_a_peaks_ripples_lie_at_r0_plus_minus_n_and_a_half_pi_over_qmax(self):
        # The rule follows from the extrema of the truncation kernel
        # sin(Q_max·Δr)/Δr: they are where tan x = x (x = Q_max·Δr), x ≈ (n + ½)π
        # for n ≥ 1. Their sign is set by the peak's width (σ·Q_max), so only
        # the positions are tested, not the sign.
        for sigma in (0.03, 0.07, 0.12):
            for qmax in (15.0, 20.0, 30.0):
                with self.subTest(sigma=sigma, qmax=qmax):
                    found = R_FINE[extrema(artefact(sigma, qmax))]
                    misses = []
                    for side in (-1, 1):
                        for n in range(1, 5):
                            expected = R0 + side * (n + 0.5) * math.pi / qmax
                            nearest = float(found[np.argmin(np.abs(found - expected))])
                            if abs(nearest - expected) > POSITION_TOLERANCE:
                                misses.append((side * n, round(expected, 3), round(nearest, 3)))
                    self.assertEqual(misses, [], "ripple extrema off r₀ ± (n + ½)π/Q_max")


if __name__ == "__main__":
    unittest.main()
