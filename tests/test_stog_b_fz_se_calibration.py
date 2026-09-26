# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The Q->0 intercept's standard error matches its scatter (1.0 review, stog-b).

fz_limit_fit reported the naive weighted-least-squares error of the Huber head
fit, sigma^2 = sum(w^2 r^2)/dof over N = D^T W^2 D. With the Huber-clipped
residuals that is biased low for an M-estimator: 0.78-0.80 of the empirical
scatter of S_meas(0) on a Gaussian head, so the FZ_REL_SE_MAX = 0.2 gate
(|denominator| >= 5 sigma) was effectively ~4 sigma. It now uses Huber's
sandwich (Huber 1981, Eq. 7.10; statsmodels RLM 'H1'):
K^2 [sum psi^2/(n-p)] / mean(psi')^2 (D^T D)^-1.
"""

import unittest

import numpy as np

from rmc_toolkits.scaling import ScalingConfig, fz_limit_fit

Q = np.arange(0.8, 30.0, 0.02)
CONFIG = ScalingConfig(qmin=0.8, qmax=30.0, rho0=0.05, b_avg_sq=1.0, b_sq_avg=2.0)


def calibration(noise, realizations=300, seed=1):
    """mean reported se(S_meas(0)) / empirical sd(S_meas(0)) over realizations."""
    rng = np.random.default_rng(seed)
    values, errors = [], []
    for _ in range(realizations):
        sq = 0.7 + 0.05 * Q + noise(rng, Q.size)
        fit = fz_limit_fit(Q, sq, 1.0, CONFIG)
        values.append(fit["s_meas_0"])
        errors.append(fit["s_meas_0_se"])
    return float(np.mean(errors) / np.std(values))


class InterceptErrorCalibrationTests(unittest.TestCase):
    # 300 realizations: the empirical sd itself is known to ~4 %.
    def test_gaussian_head(self):
        for sigma, seed in ((0.01, 1), (0.05, 4)):
            with self.subTest(sigma=sigma):
                ratio = calibration(lambda rng, n: rng.normal(0.0, sigma, n), seed=seed)
                self.assertGreater(ratio, 0.9)
                self.assertLess(ratio, 1.1)

    def test_heavy_tailed_head(self):
        ratio = calibration(lambda rng, n: 0.01 * rng.standard_t(3, n), seed=2)
        self.assertGreater(ratio, 0.9)
        self.assertLess(ratio, 1.1)

    def test_head_with_one_sided_spikes(self):
        def spiky(rng, n):
            spikes = (rng.random(n) < 0.05) * rng.uniform(0.05, 0.3, n)
            return rng.normal(0.0, 0.01, n) + spikes

        ratio = calibration(spiky, seed=3)
        self.assertGreater(ratio, 0.9)
        self.assertLess(ratio, 1.15)


if __name__ == "__main__":
    unittest.main()
