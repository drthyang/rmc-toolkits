# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The Auto StoG robust fits are Huber's M-estimator (0.6.0 integration).

``_solve_affine`` (the (a, b) fit) and ``fz_limit_fit`` (the Q->0 head fit)
re-weight by Huber IRLS with ``c = 1.345`` and a MAD scale. A weighted
least-squares pass minimises ``sum w r^2``, so its rows must be scaled by
``sqrt(w)``; before 0.6.0 both engines scaled them by ``w``, an effective weight
``w^2`` -- a redescending estimator, not the documented Huber one (95 %
Gaussian efficiency). The reference below is independent of the engine's
row-scaling ``lstsq``: it solves the weighted normal equations
``X^T W X beta = X^T W y`` with ``W = diag(w)``. Run for the engine's own pass
count from the same (unweighted) start, it must reproduce the engine to
round-off on a contaminated regression, where the old ``w^2`` weighting is far
off; iterated to convergence it satisfies Huber's estimating equation.
The JS engine (workers/autoScale.js) matches Python through the parity fixture.
"""

import unittest
from unittest import mock

import numpy as np

from rmc_toolkits import scaling
from rmc_toolkits.scaling import ScalingConfig, _fit_windows, _solve_affine, fz_limit_fit
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

HUBER_C = 1.345


def mad_scale(residuals):
    return 1.4826 * np.median(np.abs(residuals - np.median(residuals)))


def huber_weights(residuals):
    scale = mad_scale(residuals)
    if scale <= 1e-14:
        return np.ones_like(residuals)
    return np.minimum(1.0, HUBER_C * scale / np.maximum(np.abs(residuals), 1e-14 * scale))


def reference_irls(x, y, passes, blocks=None, row_power=1.0):
    """Huber IRLS by the weighted normal equations, from the unweighted solve.

    ``blocks`` lists (start, stop, min_rows) row ranges that get their own MAD
    scale (rows outside every block, or of a block with fewer than min_rows
    rows, keep weight 1). ``row_power = 2`` reproduces the pre-0.6.0 engines
    (effective weight w^2). Returns (beta, weights of the last solve).
    """
    blocks = blocks or [(0, y.size, 0)]
    weights = np.ones_like(y)
    beta = np.linalg.solve(x.T @ x, x.T @ y)
    for _ in range(passes):
        residuals = x @ beta - y
        weights = np.ones_like(y)
        for start, stop, min_rows in blocks:
            if stop - start >= min_rows:
                weights[start:stop] = huber_weights(residuals[start:stop])
        effective = weights**row_power
        beta = np.linalg.solve(x.T @ (effective[:, None] * x), x.T @ (effective * y))
    return beta, weights


def contaminated_head(seed=5):
    """A Q head with Gaussian noise and one-sided Bragg-like spikes on 12 % of rows."""
    rng = np.random.default_rng(seed)
    q = np.arange(0.8, 1.8, 0.01)
    spikes = (rng.random(q.size) < 0.12) * rng.uniform(0.05, 0.25, q.size)
    return q, 0.4 + 0.3 * q + rng.normal(0.0, 0.01, q.size) + spikes


class FzHeadFitIsHuberTests(unittest.TestCase):
    CONFIG = ScalingConfig(qmin=0.8, qmax=30.0, rho0=0.05, b_avg_sq=1.0, b_sq_avg=2.0)

    def design(self, q):
        return np.column_stack([np.ones_like(q), q - q.mean()])

    def test_matches_the_reference_irls_pass_for_pass(self):
        q, sq = contaminated_head()
        fit = fz_limit_fit(np.r_[q, 30.0], np.r_[sq, 1.0], 1.0, self.CONFIG)
        # fz_limit_fit: one unweighted solve and three re-weighted ones.
        beta, _ = reference_irls(self.design(q), sq, passes=3)
        reference_s0 = beta[0] - beta[1] * q.mean()
        self.assertAlmostEqual(fit["s_meas_0"], reference_s0, delta=1e-12)
        # The old row scaling is a different estimator on this head.
        old, _ = reference_irls(self.design(q), sq, passes=3, row_power=2.0)
        self.assertGreater(abs((old[0] - old[1] * q.mean()) - reference_s0), 1e-3)

    def test_converged_reference_solves_hubers_estimating_equation(self):
        q, sq = contaminated_head()
        x = self.design(q)
        beta, _ = reference_irls(x, sq, passes=200)
        residuals = x @ beta - sq
        # Scale fixed at the converged MAD: sum psi(r_i) x_i = 0.
        psi = np.clip(residuals, -HUBER_C * mad_scale(residuals), HUBER_C * mad_scale(residuals))
        np.testing.assert_allclose(x.T @ psi, 0.0, atol=1e-12)
        # The engine's four passes land within a small fraction of the
        # intercept's standard error of that fixed point.
        fit = fz_limit_fit(np.r_[q, 30.0], np.r_[sq, 1.0], 1.0, self.CONFIG)
        self.assertLess(abs(fit["s_meas_0"] - (beta[0] - beta[1] * q.mean())), 0.1 * fit["s_meas_0_se"])


class AffineFitIsHuberTests(unittest.TestCase):
    """The (a, b) solve re-weights each block (C1 tail, C2 low-r) by its own MAD."""

    RHO0 = 0.05

    def data(self):
        q = np.arange(20, 981) * 0.03
        r = np.arange(1, 12001) * 0.005
        g = 0.5 * (1 + np.tanh((r - 2.65) / 0.07)) + 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
        sq = fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, self.RHO0), q))
        rng = np.random.default_rng(7)
        # Bragg-like one-sided spikes on 10 % of the high-Q tail, plus noise.
        tail = q > 22.0
        spikes = tail * (rng.random(q.size) < 0.10) * rng.uniform(0.02, 0.1, q.size)
        return q, (sq + 9.0) / 10.0 + rng.normal(0.0, 0.002, q.size) + spikes

    def test_matches_the_reference_irls_on_the_captured_design(self):
        q, sq = self.data()
        config = ScalingConfig(
            qmin=0.6, qmax=29.4, rho0=self.RHO0, b_avg_sq=0.02, r0=2.5,
            r_fit_min=1.2, r_fit_max=2.25, rmax=25.0, nr=1000,
        )
        r = config.r_grid
        tail, window = _fit_windows(q, r, config)
        calls = []
        real_lstsq = np.linalg.lstsq

        def capture(design, rhs, rcond=None):
            calls.append((np.array(design), np.array(rhs)))
            return real_lstsq(design, rhs, rcond=rcond)

        with mock.patch.object(scaling.np.linalg, "lstsq", capture):
            a, b = _solve_affine(q, sq, np.zeros_like(q), r, tail, window, config)
        self.assertEqual(len(calls), 4)  # the unweighted solve + 3 IRLS passes
        design, rhs = calls[0]
        n1 = int(tail.sum())
        blocks = [(0, n1, 0), (n1, rhs.size, 4)]
        beta, weights = reference_irls(design, rhs, passes=3, blocks=blocks)
        self.assertLess(weights[:n1].min(), 0.5)  # the spikes are down-weighted
        np.testing.assert_allclose([a, b], beta, rtol=1e-9)
        old, _ = reference_irls(design, rhs, passes=3, blocks=blocks, row_power=2.0)
        self.assertGreater(abs(old[0] - beta[0]) / abs(beta[0]), 1e-4)


if __name__ == "__main__":
    unittest.main()
