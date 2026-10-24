# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""RMCProfile's transform of F(Q) and the propagation of F(Q)'s σ through it.

Gates, strongest first (rmc_toolkits/ft_uncertainty.py):

1. the reproduction tests on real RMCProfile output: the committed demo run's
   G(r) column from its F(Q) column (CI), and the two local runs' F(Q) columns
   from their input files through the renormalisation and the box convolution
   (data-gated: ``data/`` is gitignored);
2. a brute-force Jacobian: every operator column from the value path, so the σ
   and covariance paths must reproduce J·Σ·Jᵀ;
3. Monte Carlo: seeded noise realisations, transformed, give the propagated
   variances and correlations (no formula shared with 2), also through the
   box convolution;
4. closed forms: white noise against the continuum integral, the Nyquist-grid
   orthogonality, the continuum correlation law;
5. external goldens adapted to RMCProfile's convention: pystog's (MATERIA's
   copy) and RMCProfile's StoG on FeCoSn (data-gated).
"""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

from rmc_toolkits.ft_uncertainty import (
    CONVOLUTION_CHECK_LIMIT,
    RMC_FT_TOLERANCE,
    box_convolution,
    convolution_check_points,
    convolution_sigma_change,
    correlation_summary,
    file_sigma_column,
    assess_ft_uncertainty,
    ft_uncertainty_for_run,
    noise_autocorrelation,
    reproduce_input,
    reproduce_transform,
    rmc_ft,
    rmc_ft_coefficients,
    rmc_ft_covariance,
    rmc_ft_sigma,
    run_datasets,
    uniform_step,
    xray_dataset,
)
from rmc_toolkits.parsers import read_rmc_csv, read_stog_xy
from rmc_toolkits.transforms import fq_to_gpdf, omitted_low_q_correction

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_ft_uncertainty_fixture import STEM, synthetic_run, toy_f
except ModuleNotFoundError:  # `python -m unittest tests.test_ft_uncertainty`
    from tests.generate_ft_uncertainty_fixture import STEM, synthetic_run, toy_f

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"
DATA = ROOT / "data"
PYSTOG = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "pystog_golden.json"


def _columns(path: Path):
    series = read_rmc_csv(path)
    return series.data


def _dense_box_kernel(q: np.ndarray, radius: float) -> np.ndarray:
    """K as a dense matrix from the exact Q differences (no Toeplitz shortcut)."""
    step = uniform_step(q)
    diff = q[:, None] - q[None, :]
    safe = np.where(diff == 0, 1.0, diff)
    first = np.where(diff == 0, radius, np.sin(diff * radius) / safe)
    second = np.sin((q[:, None] + q[None, :]) * radius) / (q[:, None] + q[None, :])
    return (first - second) * (step * q[None, :]) / (math.pi * q[:, None])


class TransformTests(unittest.TestCase):
    def test_the_value_is_the_rectangle_rule_without_two_over_pi(self):
        q = 0.8 + 0.05 * np.arange(300)
        f = toy_f(q)
        r = np.array([0.3, 1.7, 2.51, 9.9])
        got = rmc_ft(q, f, r)
        for j, at in enumerate(r):
            want = sum(0.05 * qk * fk * math.sin(qk * at) for qk, fk in zip(q, f))
            self.assertAlmostEqual(got[j], want, delta=1e-13)

    def test_a_non_uniform_or_unsorted_grid_is_refused(self):
        with self.assertRaisesRegex(ValueError, "not uniform"):
            rmc_ft_coefficients(np.array([0.5, 0.6, 0.71, 0.8]))
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            rmc_ft_coefficients(np.array([0.5, 0.6, 0.6, 0.8]))
        with self.assertRaisesRegex(ValueError, "at least two"):
            rmc_ft_coefficients(np.array([0.5]))


class PropagationTests(unittest.TestCase):
    q = 0.6 + 0.31 * np.arange(70)
    sigma = 0.004 + 0.0007 * q + 0.003 * np.exp(-q)
    r_uniform = 0.05 + 0.137 * np.arange(24)
    r_ragged = np.array([0.02, 0.4, 0.41, 1.3, 2.75, 2.8, 6.1, 9.0])

    def jacobian(self, r: np.ndarray) -> np.ndarray:
        """J[i, k] = ∂G(r_i)/∂F_k from the value path alone."""
        return np.stack([rmc_ft(self.q, np.eye(self.q.size)[k], r) for k in range(self.q.size)], axis=1)

    def test_sigma_is_the_jacobian_law(self):
        jac = self.jacobian(self.r_uniform)
        want = np.sqrt((jac ** 2) @ (self.sigma ** 2))
        got = rmc_ft_sigma(self.q, self.sigma, self.r_uniform)
        np.testing.assert_allclose(got, want, rtol=1e-12, atol=0)

    def test_covariance_is_j_sigma_jt_on_uniform_and_ragged_grids(self):
        for r in (self.r_uniform, self.r_ragged):
            jac = self.jacobian(r)
            want = (jac * self.sigma ** 2) @ jac.T
            got = rmc_ft_covariance(self.q, self.sigma, r)
            scale = float(np.max(np.diag(want)))
            self.assertLess(float(np.max(np.abs(got - want))), 1e-12 * scale)

    def test_the_covariance_diagonal_is_sigma_squared_bit_for_bit(self):
        cov = rmc_ft_covariance(self.q, self.sigma, self.r_uniform)
        sigma = rmc_ft_sigma(self.q, self.sigma, self.r_uniform)
        self.assertTrue(np.array_equal(np.sqrt(np.diag(cov)), sigma))

    def test_the_covariance_is_positive_semidefinite(self):
        cov = rmc_ft_covariance(self.q, self.sigma, self.r_uniform)
        jitter = 1e-12 * float(np.max(np.diag(cov)))
        np.linalg.cholesky(cov + jitter * np.eye(cov.shape[0]))

    def test_white_noise_gives_the_analytic_sigma(self):
        # σ_k = σ0: σ(r)² = σ0² ΔQ² Σ Q_k² sin²(Q_k r) exactly, and the rectangle
        # sum is the midpoint rule for ΔQ σ0² ∫ Q² sin²(Qr) dQ over
        # [Q_0 − ΔQ/2, Q_last + ΔQ/2] (antiderivative below), to O((ΔQ r)²).
        dq, sigma0 = 0.01, 0.003
        q = 0.8 + dq * np.arange(2619)
        r = np.array([0.5, 2.34, 5.0, 10.0])
        got = rmc_ft_sigma(q, np.full(q.size, sigma0), r)
        explicit = np.sqrt([sigma0 ** 2 * dq ** 2 * float(np.sum(q ** 2 * np.sin(q * at) ** 2)) for at in r])
        np.testing.assert_allclose(got, explicit, rtol=1e-12)

        def antiderivative(x, at):
            return (x ** 3 / 6 - x ** 2 / (4 * at) * math.sin(2 * x * at) - x / (4 * at ** 2) * math.cos(2 * x * at)
                    + math.sin(2 * x * at) / (8 * at ** 3))

        lo, hi = q[0] - dq / 2, q[-1] + dq / 2
        continuum = np.sqrt([dq * sigma0 ** 2 * (antiderivative(hi, at) - antiderivative(lo, at)) for at in r])
        np.testing.assert_allclose(got, continuum, rtol=1e-3)

    def test_monte_carlo_reproduces_the_propagated_covariance(self):
        q = 0.5 + 0.12 * np.arange(160)
        f0 = toy_f(q)
        sigma = 0.01 * (1 + 0.08 * q)
        r = np.array([0.15, 1.0, 1.08, 2.5, 7.7])
        cov = rmc_ft_covariance(q, sigma, r)
        rng = np.random.default_rng(20261024)
        draws = 20000
        samples = np.stack([rmc_ft(q, f0 + sigma * rng.standard_normal(q.size), r) for _ in range(draws)])
        empirical = np.cov(samples, rowvar=False, bias=True)
        sd, emp_sd = np.sqrt(np.diag(cov)), np.sqrt(np.diag(empirical))
        # SE of a sample variance σ²·√(2/M) ≈ 1 %: 5 % is 5 SE.
        np.testing.assert_allclose(np.diag(empirical), np.diag(cov), rtol=0.05)
        rho = cov / np.outer(sd, sd)
        emp_rho = empirical / np.outer(emp_sd, emp_sd)
        # SE of a sample correlation ≤ 1/√M ≈ 0.007: 0.035 is 5 SE.
        self.assertLess(float(np.max(np.abs(emp_rho - rho))), 0.035)
        # 0.08 Å apart, half of π/Q_max ≈ 0.16 Å: correlated (the continuum law
        # for σ ∝ 1 + 0.08 Q, u = Q_max·Δ ≈ 1.6, gives ≈ 0.3), never independent.
        self.assertGreater(rho[1, 2], 0.25)

    def test_monte_carlo_through_the_box_convolution(self):
        # The run transforms s·K·F_input: the exact σ is through B = A·K; the
        # propagation leaves K out, which convolution_sigma_change measures.
        q = 0.8 + 0.05 * np.arange(300)
        f_input = toy_f(q)
        sigma = 0.002 * (1 + 0.15 * q)
        radius, scale = 25.0, 1.15
        r = np.array([2.0, 6.0, 9.9])
        rng = np.random.default_rng(7)
        draws = 20000
        samples = np.stack([
            rmc_ft(q, scale * box_convolution(q, f_input + sigma * rng.standard_normal(q.size), radius), r)
            for _ in range(draws)
        ])
        empirical_sd = samples.std(axis=0)
        approx = rmc_ft_sigma(q, scale * sigma, r)
        exact = approx * (1 + convolution_sigma_change(q, scale * sigma, r, radius))
        np.testing.assert_allclose(empirical_sd, exact, rtol=0.025)
        np.testing.assert_allclose(empirical_sd, approx, rtol=0.03)

    def test_white_integrand_noise_on_the_nyquist_grid_is_uncorrelated(self):
        # Q_k = kΔQ from 0, σ_k ∝ 1/Q_k (white noise in Q·F, the integrand):
        # Σ_k sin(πkj/K) sin(πkl/K) = (K/2) δ_jl on r_j = jπ/Q_max.
        big_k, dq, sigma0 = 400, 0.06, 0.02
        q = dq * np.arange(big_k + 1)
        sigma = np.empty_like(q)
        sigma[0] = 1.0  # c_0 = ΔQ·0 = 0: never enters
        sigma[1:] = sigma0 / q[1:]
        r = (np.arange(60) + 1) * math.pi / q[-1]
        cov = rmc_ft_covariance(q, sigma, r)
        expected = (dq * sigma0) ** 2 * big_k / 2
        np.testing.assert_allclose(np.diag(cov) / expected, 1.0, atol=1e-12)
        off = cov - np.diag(np.diag(cov))
        self.assertLess(float(np.max(np.abs(off))), 1e-12 * expected)

    def test_white_f_noise_follows_the_continuum_correlation_law(self):
        # σ_k const: integrand noise ∝ Q, so ρ(u) = 3[(u² − 2) sin u + 2u cos u]/u³
        # (u = Q_max Δ), ρ(π) = −6/π² — MATERIA's "white S(Q)" law.
        qmax = 25.0
        q = qmax * np.arange(5001) / 5000
        sigma = np.full(q.size, 0.01)
        for u in (0.3, 1.0, math.pi, 5.0):
            p0, pu = noise_autocorrelation(q, sigma, np.array([0.0, u / qmax]))
            law = 3 * ((u * u - 2) * math.sin(u) + 2 * u * math.cos(u)) / u ** 3
            self.assertAlmostEqual(pu / p0, law, places=3)
        p0, pnyq = noise_autocorrelation(q, sigma, np.array([0.0, math.pi / qmax]))
        self.assertAlmostEqual(pnyq / p0, -6 / math.pi ** 2, places=3)

    def test_the_correlation_summary_counts_independent_points(self):
        q = 0.01 * np.arange(2001)  # 0 … 20 Å⁻¹
        r = 10 + 0.02 * np.arange(1001)  # 10 … 30 Å
        nyquist = (30 - 10) * 20 / math.pi
        white_integrand = np.empty_like(q)
        white_integrand[0] = 1.0
        white_integrand[1:] = 0.01 / q[1:]
        summary = correlation_summary(q, white_integrand, r)
        self.assertAlmostEqual(summary["nyquistPoints"], nyquist, places=9)
        self.assertGreater(summary["effectivePoints"] / nyquist, 0.95)
        self.assertLess(summary["effectivePoints"] / nyquist, 1.05)
        # σ const (the integrand's noise ∝ Q): ≈ 5/9 of the Nyquist count.
        flat = correlation_summary(q, np.full(q.size, 0.01), r)
        self.assertGreater(flat["effectivePoints"] / nyquist, 5 / 9 * 0.95)
        self.assertLess(flat["effectivePoints"] / nyquist, 5 / 9 * 1.05)
        self.assertEqual(summary["points"], r.size)
        self.assertAlmostEqual(summary["gridStep"], 0.02, places=12)
        self.assertAlmostEqual(summary["nyquistStep"], math.pi / 20, places=12)


class BoxConvolutionTests(unittest.TestCase):
    q = 0.8 + 0.05 * np.arange(300)

    def test_the_toeplitz_form_is_the_dense_kernel(self):
        f = toy_f(self.q)
        got = box_convolution(self.q, f, 25.0)
        want = _dense_box_kernel(self.q, 25.0) @ f
        self.assertLess(float(np.max(np.abs(got - want))), 1e-13 * float(np.max(np.abs(want))))

    def test_the_sigma_change_is_the_exact_chain(self):
        sigma = 0.002 * (1 + 0.15 * self.q)
        r = np.array([1.5, 6.0, 9.95, 24.0])
        step = uniform_step(self.q)
        a = np.sin(np.outer(r, self.q)) * (step * self.q)
        b = a @ _dense_box_kernel(self.q, 25.0)
        want = np.sqrt((b ** 2) @ (sigma ** 2)) / np.sqrt((a ** 2) @ (sigma ** 2)) - 1
        got = convolution_sigma_change(self.q, sigma, r, 25.0)
        np.testing.assert_allclose(got, want, rtol=0, atol=1e-12)
        # Inside R it is small; at the edge it is not.
        self.assertLess(abs(got[0]), CONVOLUTION_CHECK_LIMIT)
        self.assertGreater(abs(got[-1]), CONVOLUTION_CHECK_LIMIT)

    def test_the_checked_points_span_the_last_nyquist_step(self):
        r = 1.5 + 0.04 * np.arange(213)
        points = convolution_check_points(r, 15.75)
        self.assertEqual(points.size, 9)
        self.assertEqual(points[0], r[0])
        self.assertEqual(points[-1], r[-1])
        self.assertAlmostEqual(points[-1] - points[1], math.pi / 15.75, places=12)

    def test_the_input_is_reproduced_through_the_chain(self):
        f_input = toy_f(self.q)
        f_run = 1.15 * box_convolution(self.q, f_input, 25.0) - 2.7e-4
        result = reproduce_input(self.q, f_run, f_input, radius=25.0)
        self.assertTrue(result["pass"])
        self.assertAlmostEqual(result["scale"], 1.15, places=12)
        self.assertAlmostEqual(result["offset"], -2.7e-4, places=14)
        # Without the convolution it does not reproduce.
        self.assertFalse(reproduce_input(self.q, f_run, f_input)["pass"])


class DemoReproductionTests(unittest.TestCase):
    """The committed demo run (GaTa4Se8 250 K): RMCProfile's own files."""

    @classmethod
    def setUpClass(cls):
        cls.q, _, cls.f = _columns(DEMO / "GTS_250K_XFQ1.csv")
        cls.r, _, cls.g = _columns(DEMO / "GTS_250K_FT_XFQ1.csv")

    def test_the_transform_reproduces_the_runs_g_of_r(self):
        result = reproduce_transform(self.q, self.f, self.r, self.g)
        self.assertTrue(result["pass"])
        self.assertLess(result["residual"], 1e-13)  # measured 1.5e-14 of max |G|
        self.assertAlmostEqual(result["offset"], 2.7257839810e-4, delta=1e-12)

    def test_the_7_decimal_fq1_file_reproduces_it_within_the_tolerance(self):
        q, _, f = _columns(DEMO / "GTS_250K_FQ1.csv")
        result = reproduce_transform(q, f, self.r, self.g)
        self.assertTrue(result["pass"])
        self.assertLess(result["residual"], 1e-7)  # measured 4.5e-8: the file's 7 decimals

    def test_every_other_convention_fails(self):
        dq = uniform_step(self.q)
        sines = np.sin(np.outer(self.r, self.q))
        weights = dq * self.q
        scale = float(np.max(np.abs(self.g)))

        def misfit(model):
            return float(np.max(np.abs(self.g - model))) / scale

        offset = reproduce_transform(self.q, self.f, self.r, self.g)["offset"]
        rect = sines @ (weights * (self.f + offset))
        self.assertLess(misfit(rect), RMC_FT_TOLERANCE)
        self.assertGreater(misfit(sines @ (weights * self.f)), 100 * RMC_FT_TOLERANCE)  # the offset left out
        trapezoid = weights.copy()
        trapezoid[[0, -1]] *= 0.5
        self.assertGreater(misfit(sines @ (trapezoid * (self.f + offset))), 100 * RMC_FT_TOLERANCE)
        self.assertGreater(misfit((2 / math.pi) * rect), 0.1)  # D(r) = (2/π)·this
        lorch = np.sinc(self.q / self.q[-1])  # sin(πQ/Qmax)/(πQ/Qmax)
        self.assertGreater(misfit(sines @ (weights * lorch * (self.f + offset))), 0.01)
        rho0 = 0.046690
        keen = rect / (2 * math.pi ** 2 * rho0 * self.r)  # Keen's G(r) = D(r)/(4πρ0 r), no 2/π
        self.assertGreater(misfit(keen), 0.1)

    def test_the_demo_has_no_ribbon_and_says_which_file_is_missing(self):
        payload = ft_uncertainty_for_run(DEMO / "GTS_250K_FT_XFQ1.csv")
        self.assertEqual(payload["status"], "unavailable")
        self.assertEqual(
            payload["reason"],
            "F(Q) carries no uncertainties: scale_ft_rmc.fq (named in GTS_250K.dat) is not in the run folder",
        )
        self.assertIsNone(payload["sigma"])
        self.assertTrue(payload["reproduction"]["gr"]["pass"])
        self.assertIsNone(payload["reproduction"]["fq"])
        self.assertEqual(payload["transform"], {
            "qMin": 0.8, "qMax": 26.98, "dQ": payload["transform"]["dQ"], "nQ": 2619,
            "rMin": 2.34, "rMax": 20.0, "nR": 1767, "convolve": True, "radius": 51.82,
        })
        self.assertAlmostEqual(payload["transform"]["dQ"], 0.01, places=12)
        self.assertEqual(payload["sources"], {
            "transform": "GTS_250K_FT_XFQ1.csv", "fq": "GTS_250K_XFQ1.csv",
            "settings": "GTS_250K.dat", "partials": "GTS_250K_PDFpartials.csv",
        })
        json.dumps(payload, allow_nan=False)

    def test_the_run_control_file_names_the_x_ray_dataset(self):
        settings = run_datasets((DEMO / "GTS_250K.dat").read_text())
        self.assertEqual(settings["rSpacing"], 0.01)
        self.assertEqual(settings["numberDensity"], 0.04669)
        dataset = xray_dataset(settings, 1)
        self.assertEqual(dataset["file"], "scale_ft_rmc.fq")
        self.assertTrue(dataset["convolve"] and dataset["fittedScale"] and dataset["fittedOffset"])
        self.assertEqual(dataset["reciprocalParameters"], [1, 2619, 1])
        self.assertEqual(dataset["realParameters"], [234, 2000, 1])
        # The G(r) file's r grid is REAL_SPACE_PARAMETERS × R_SPACING.
        first, last, _ = dataset["realParameters"]
        np.testing.assert_allclose(self.r, settings["rSpacing"] * np.arange(first, last + 1), atol=1e-12)
        self.assertIsNone(xray_dataset(settings, 2))


@unittest.skipUnless((DATA / "250K_try1" / "scale_ft_rmc.fq").is_file() and (DATA / "5K_try1" / "scale_ft_rmc.fq").is_file(),
                     "the local GTS 250 K / GNS 5 K runs (data/, gitignored)")
class LocalRunChainTests(unittest.TestCase):
    """RMCProfile 6.7.9.5's chain on two real runs: input → F(Q)_Expt → G(r)."""

    RUNS = (("250K_try1", "GTS_250K", 1.1475350510), ("5K_try1", "GaNb4Se8_5K", 1.2831308333))

    def test_the_input_reproduces_the_runs_f_of_q_and_g_of_r(self):
        for folder, stem, scale in self.RUNS:
            with self.subTest(run=folder):
                payload = ft_uncertainty_for_run(DATA / folder / f"{stem}_FT_XFQ1.csv")
                self.assertEqual(payload["reason"], "F(Q) carries no uncertainties: scale_ft_rmc.fq has no error column")
                gr, fq = payload["reproduction"]["gr"], payload["reproduction"]["fq"]
                self.assertLess(gr["residual"], 1e-13)
                self.assertLess(fq["residual"], 1e-13)
                self.assertTrue(fq["convolve"])
                self.assertAlmostEqual(fq["scale"], scale, delta=1e-9)
                # The F(Q)_Expt column is s·K·F_input − c and the transform adds c back.
                self.assertLess(abs(gr["offset"] + fq["offset"]), 1e-14)

    def test_leaving_the_box_convolution_out_changes_sigma_by_under_0_2_percent(self):
        for folder, stem, _ in self.RUNS:
            with self.subTest(run=folder):
                q, _, _ = _columns(DATA / folder / f"{stem}_XFQ1.csv")
                r = _columns(DATA / folder / f"{stem}_FT_XFQ1.csv")[0]
                radius = float(_columns(DATA / folder / f"{stem}_PDFpartials.csv")[0][-1])
                for label, sigma in (("flat", np.full(q.size, 0.01)), ("rising", 0.002 * (1 + 0.15 * q))):
                    change = convolution_sigma_change(q, sigma, r[::20], radius)
                    self.assertLess(float(np.max(np.abs(change))), 0.002, label)


class ExternalGoldenTests(unittest.TestCase):
    """MATERIA's pystog golden (S(Q) → G(r) with dG), adapted to RMCProfile's convention."""

    @classmethod
    def setUpClass(cls):
        cls.golden = json.loads(PYSTOG.read_text(encoding="utf-8"))
        cls.q = 0.4 + 0.02 * np.arange(1181)
        cls.r = 0.1 * np.arange(1, 201)
        s = np.ones_like(cls.q)
        for d, n, u in ((2.5, 0.6, 0.006), (3.6, 0.45, 0.009), (4.4, 0.35, 0.011)):
            s = s + (n * np.sin(cls.q * d)) / (cls.q * d) * np.exp(-0.5 * u * cls.q * cls.q * d)
        cls.f = s - 1.0  # RMCProfile's F(Q) = S(Q) − 1
        cls.dsq = 0.002 * (1 + 0.15 * cls.q)

    def endpoint(self, values: np.ndarray) -> np.ndarray:
        q, r = self.q, self.r
        return values[0] * np.sin(q[0] * r), values[-1] * np.sin(q[-1] * r)

    def test_the_value_is_pi_over_2_pystog_plus_the_endpoint_halves(self):
        # pystog: (2/π)·trapezoid; RMCProfile: rectangle, no 2/π. So RMC =
        # (π/2)·pystog + (ΔQ/2)·(the two endpoint terms of Q·F·sin(Qr)).
        first, last = self.endpoint(self.q * self.f)
        adapted = 0.5 * math.pi * np.array(self.golden["plain"]["g"]) + 0.01 * (first + last)
        got = rmc_ft(self.q, self.f, self.r)
        self.assertLess(float(np.max(np.abs(got - adapted))), 1e-12 * float(np.max(np.abs(got))))

    def test_sigma_is_pystogs_once_its_endpoint_weights_are_reconciled(self):
        # pystog integrates the variance with its own trapezoid (endpoint weight
        # ½ΔQ²); RMCProfile's rectangle gives every point ΔQ². So σ_RMC² =
        # (π/2)²·dG² + ΔQ²·(e_first + e_last)/2, e = (Q dS sin Qr)².
        first, last = self.endpoint(self.q * self.dsq)
        dg = np.array(self.golden["plain"]["dg"])
        adapted = np.sqrt((0.5 * math.pi * dg) ** 2 + 0.02 ** 2 * (first ** 2 + last ** 2) / 2)
        got = rmc_ft_sigma(self.q, self.dsq, self.r)
        np.testing.assert_allclose(got, adapted, rtol=1e-12)
        raw = float(np.max(np.abs(got / (0.5 * math.pi * dg) - 1)))
        self.assertGreater(raw, 1e-4)  # the convention difference is real …
        self.assertLess(raw, 2e-3)  # … and an O(1/N_Q) endpoint effect


STOG = DATA / "stog_tests" / "199K"


@unittest.skipUnless(all((STOG / name).is_file() for name in ("scale.fq", "scale.gr", "stog_input.dat")),
                     "the FeCoSn 199 K StoG run (data/stog_tests, gitignored)")
class StogGoldenTests(unittest.TestCase):
    """RMCProfile's StoG (FeCoSn 199 K) against transforms.py, and RMCProfile's FT against both."""

    def test_stog_is_the_trapezoid_with_the_low_q_term_and_rmcprofile_is_not(self):
        deck = [line.strip() for line in (STOG / "stog_input.dat").read_text().splitlines()]
        qmin, qmax = (float(value) for value in deck[2].split())
        rho0 = float(deck[10])
        q_all, s_all = read_stog_xy(STOG / "scale.fq")[:2]
        window = (q_all >= qmin) & (q_all <= qmax)
        q, s = q_all[window], s_all[window]
        reference = read_stog_xy(STOG / "scale.gr")
        r, g_ref = reference[0], reference[1]
        keep = r >= 0.05
        fq = q * (s - 1.0)
        stog = fq_to_gpdf(q, fq, r, low_q_correction=True)
        g = 1.0 + stog / (4 * math.pi * rho0 * r)
        self.assertLess(float(np.max(np.abs(g[keep] - g_ref[keep]))), 1e-10)
        # RMCProfile's operator on the same S(Q): (π/2)·(StoG − its low-Q term)
        # plus the rectangle's endpoint halves.
        first = q[0] * (s[0] - 1.0) * np.sin(q[0] * r)
        last = q[-1] * (s[-1] - 1.0) * np.sin(q[-1] * r)
        dq = uniform_step(q)
        adapted = 0.5 * math.pi * (stog - omitted_low_q_correction(q, fq, r)) + 0.5 * dq * (first + last)
        got = rmc_ft(q, s - 1.0, r)
        self.assertLess(float(np.max(np.abs(got - adapted))), 1e-10 * float(np.max(np.abs(got))))


class RunFolderTests(unittest.TestCase):
    """ft_uncertainty_for_run on synthetic run folders (the fixture's construction)."""

    def run_folder(self, files: dict, resolve=None) -> dict:
        with tempfile.TemporaryDirectory() as folder:
            for name, text in files.items():
                (Path(folder) / name).write_text(text, encoding="utf-8")
            return ft_uncertainty_for_run(Path(folder) / f"{STEM}_FT_XFQ1.csv", resolve=resolve)

    def test_a_clean_error_column_gives_the_ribbon(self):
        payload = self.run_folder(synthetic_run())
        self.assertEqual(payload["status"], "propagated")
        self.assertIsNone(payload["reason"])
        q = 0.8 + 0.05 * np.arange(300)
        r = 1.5 + 0.04 * np.arange(213)
        want = rmc_ft_sigma(q, 1.15 * 0.002 * (1 + 0.15 * q), r)
        np.testing.assert_allclose(payload["sigma"], want, rtol=1e-9)
        self.assertLessEqual(payload["convolutionCheck"]["maxChange"], CONVOLUTION_CHECK_LIMIT)
        self.assertGreater(payload["correlation"]["neighbour"], 0.5)

    def test_the_resolver_guards_the_input_path(self):
        # A path the caller's guard refuses is never read: as if it were absent.
        def refuse(path):
            raise PermissionError(f"{path.name} is outside the data root")

        payload = self.run_folder(synthetic_run(), resolve=refuse)
        self.assertEqual(payload["status"], "unavailable")
        self.assertEqual(
            payload["reason"],
            "F(Q) carries no uncertainties: scale_ft_rmc.fq (named in rmc_run.dat) is not in the run folder",
        )
        self.assertNotIn("input", payload["sources"])

    def test_a_file_that_is_not_a_transform_is_refused(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "run_FQ1.csv"
            path.write_text("Q, a, b\n1, 2, 3\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "not an RMCProfile transform file"):
                ft_uncertainty_for_run(path)

    def test_the_messages_name_what_is_missing(self):
        full = synthetic_run()
        cases = {
            "scale_ft_rmc.fq": "F(Q) carries no uncertainties: scale_ft_rmc.fq (named in rmc_run.dat) is not in the run folder",
            "rmc_run.dat": "no run-control file rmc_run.dat beside rmc_run_FT_XFQ1.csv",
            "rmc_run_XFQ1.csv": "no F(Q) file (rmc_run_XFQ1.csv or rmc_run_FQ1.csv) to reproduce the transform from",
            "rmc_run_PDFpartials.csv": ("rmc_run_PDFpartials.csv is not in the run folder: the radius of the box "
                                        "convolution RMCProfile applied is unknown"),
        }
        for missing, reason in cases.items():
            with self.subTest(missing=missing):
                payload = self.run_folder({name: text for name, text in full.items() if name != missing})
                self.assertEqual((payload["status"], payload["reason"]), ("unavailable", reason))
                self.assertIsNone(payload["sigma"])


class SmallPieceTests(unittest.TestCase):
    def test_a_sigma_column_is_found_by_its_header(self):
        self.assertIsNone(file_sigma_column(["r(A)", "X_ray-calc", "X_ray_exp_renorm"]))
        self.assertEqual(file_sigma_column(["r(A)", "X_ray-calc", "X_ray_exp_renorm", "sigma"]), 3)
        self.assertEqual(file_sigma_column(["r", "calc", "exp", "dG err"]), 3)
        self.assertIsNone(file_sigma_column(["sigma", "calc", "exp"]))  # never the x column

    def test_run_datasets_reads_flags_with_and_without_a_value(self):
        settings = run_datasets(
            "R_SPACING :: 0.02 Angstrom\nXRAY_RECIPROCAL_SPACE_DATA ::\n  > FILENAME :: a.fq\n  > CONVOLVE\n"
            "NEUTRON_REAL_SPACE_DATA ::\n  > FILENAME :: b.gr\nXRAY_RECIPROCAL_SPACE_DATA ::\n  > FILENAME :: c.fq\n"
            "END ::\nXRAY_RECIPROCAL_SPACE_DATA ::\n"
        )
        self.assertEqual(settings["rSpacing"], 0.02)
        self.assertEqual([entry["file"] for entry in settings["datasets"]], ["a.fq", "b.gr", "c.fq"])
        self.assertEqual(xray_dataset(settings, 1)["file"], "a.fq")
        self.assertTrue(xray_dataset(settings, 1)["convolve"])
        self.assertEqual(xray_dataset(settings, 2)["file"], "c.fq")
        self.assertFalse(xray_dataset(settings, 2)["convolve"])

    def test_the_payload_without_an_f_of_q_file(self):
        r = np.linspace(1, 2, 11)
        payload = assess_ft_uncertainty(r=r, g=np.sin(r), names={"transform": "x_FT_XFQ1.csv"})
        self.assertEqual(payload["status"], "unavailable")
        self.assertEqual(payload["reason"], "no run-control .dat beside x_FT_XFQ1.csv")


if __name__ == "__main__":
    unittest.main()
