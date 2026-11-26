# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Known-truth acceptance test for check_ripples' first principles (rmc_toolkits/ripples.py), synthetic data only.

What check_ripples ships besides RMCProfile's fit space (tests/test_ripples.py)
is the r-resolution of each dataset's transform and the rule its checks state:
ripples repeat every 2π/Q_max, a peak at r₀ has them between about r₀ ± (n + ½)π/Q_max
and r₀ ± (n + 1)π/Q_max, and the kernel's FWHM is 3.79/Q_max for a rectangle from ~0. Both are tested
here against numpy alone, independent of the engine's own helpers:

* ``resolution()`` on the rectangle Q = ΔQ … Q_max: ripplePeriod = 2π/Q_max and
  fwhm = 3.79/Q_max within 1 %, for Q_max 15, 20 and 25 Å⁻¹;
* a single Gaussian peak A[g(r − r₀) − g(r + r₀)] in D(r), its exact sine
  transform cut at Q_max and transformed back by the rectangle rule
  (2/π) Σ_k ΔQ F(Q_k) sin(Q_k r), Q_k = ΔQ … Q_max, ΔQ = 0.01; the true
  termination artefact is that less the untruncated peak. For n = 1 … 4 on
  both sides, σ 0.03, 0.07 and 0.12 Å and Q_max 15, 20 and 30 Å⁻¹, its extrema
  (x = Q_max·Δr) lie between the root of tan x = x near (n + ½)π and (n + 1)π, to
  INTERVAL_TOLERANCE, and at the root of tan(x − φ) = x, tan φ = Q_max·σ²/Δr, to
  PHASE_TOLERANCE; the rule Δr = nπ/Q_max, a negative control, must fail that at
  every Q_max.

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
INTERVAL_TOLERANCE = 0.01  # Å
PHASE_TOLERANCE = 0.015  # Å
NS = range(1, 5)


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


def kernel_root(n, qmax):
    """Δr of the sin x/x kernel's n-th extremum: tan x = x, x = Q_max·Δr, just inside (n + ½)π."""
    x = (n + 0.5) * math.pi
    for _ in range(200):
        x = n * math.pi + math.atan(x)
    return x / qmax


def phase_rule(n, sigma, qmax):
    """Δr of a Gaussian peak's n-th ripple extremum: tan(x − φ) = x, x = Q_max·Δr, tan φ = Q_max·σ²/Δr.

    A Gaussian of width σ cut at Q_max leaves Δ ∝ −∫_{Q_max}^∞ e^{−Q²σ²/2} cos(QΔr) dQ
    ≈ e^{−Q_max²σ²/2}[sin(Q_maxΔr)/Δr − Q_maxσ² cos(Q_maxΔr)/Δr² + …] (by parts:
    the next asymptotic term, from the envelope's slope at Q_max) ∝ sin(x − φ)/Δr,
    whose extrema solve tan(x − φ) = x: the kernel's own (tan x = x) shifted
    outward by φ ∈ [0, π/2), toward (n + 1)π for a broad peak.
    """
    x = (n + 0.5) * math.pi
    for _ in range(200):
        x = n * math.pi + math.atan(x) + math.atan(qmax * qmax * sigma * sigma / x)
    return x / qmax


def nearest_offsets(sigma, qmax):
    """For each side and n: (side, n, the nearest extremum's |Δr|)."""
    found = R_FINE[extrema(artefact(sigma, qmax))]
    out = []
    for side in (-1, 1):
        for n in NS:
            guess = phase_rule(n, sigma, qmax)
            nearest = float(found[np.argmin(np.abs(found - (R0 + side * guess)))])
            out.append((side, n, abs(nearest - R0)))
    return out


SIGMAS = (0.03, 0.07, 0.12)
QMAXES = (15.0, 20.0, 30.0)


class RipplePositions(unittest.TestCase):
    """Where a peak's ripples fall (the r-resolution line of check_ripples)."""

    def test_extrema_lie_between_the_kernel_root_and_n_plus_one_pi_over_qmax(self):
        for sigma in SIGMAS:
            for qmax in QMAXES:
                with self.subTest(sigma=sigma, qmax=qmax):
                    outside = [(side * n, round(x, 3)) for side, n, x in nearest_offsets(sigma, qmax)
                               if not (kernel_root(n, qmax) - INTERVAL_TOLERANCE
                                       <= x <= (n + 1) * math.pi / qmax + INTERVAL_TOLERANCE)]
                    self.assertEqual(outside, [], "extrema outside [tan x = x root, (n + 1)π]/Q_max")

    def test_extrema_follow_the_phase_rule(self):
        for sigma in SIGMAS:
            for qmax in QMAXES:
                with self.subTest(sigma=sigma, qmax=qmax):
                    misses = [(side * n, round(x, 3), round(phase_rule(n, sigma, qmax), 3))
                              for side, n, x in nearest_offsets(sigma, qmax)
                              if abs(x - phase_rule(n, sigma, qmax)) > PHASE_TOLERANCE]
                    self.assertEqual(misses, [], "extrema off the phase rule")

    def test_the_rule_n_pi_over_qmax_fails_at_every_qmax(self):
        # Negative control: the same check with the wrong rule Δr = nπ/Q_max
        # must miss by more than PHASE_TOLERANCE somewhere at every Q_max.
        for qmax in QMAXES:
            with self.subTest(qmax=qmax):
                worst = max(abs(x - n * math.pi / qmax) for sigma in SIGMAS
                            for _side, n, x in nearest_offsets(sigma, qmax))
                self.assertGreater(worst, PHASE_TOLERANCE, f"nπ/Q_max not rejected at Q_max {qmax}")


if __name__ == "__main__":
    unittest.main()
