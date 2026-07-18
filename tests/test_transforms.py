# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

from pathlib import Path
import unittest

import numpy as np

from rmc_toolkits.parsers import StogInput, read_stog_xy
from rmc_toolkits.transforms import (
    density_line,
    enforce_low_r,
    first_peak_zero,
    fourier_filter,
    fq_to_gpdf,
    fq_to_sq,
    g_to_gk,
    g_to_gpdf,
    gk_to_dr,
    gk_to_g,
    gpdf_to_fq,
    gpdf_to_g,
    lorch_window,
    low_q_correction_basis,
    omitted_low_q_correction,
    sq_to_fk,
    sq_to_fq,
)

ROOT = Path(__file__).resolve().parents[1]
STOG_RUN = ROOT / "data" / "stog_tests" / "stog_59438"

# stog_59438: a complete Fortran stog run of Mn3Sn (POWGEN PG3), local-only and
# gitignored with the rest of data/. The run's parameters are recovered from its
# own outputs, so the sample-backed tests need only the rebinned data file, not
# the interactive stog.inp (reading a real stog.inp is exercised by the x-ray
# run in tests/test_scaling.py). scale.fq vs the rebinned input pins a = 10,
# b = -9 and Qmin = 1.0 (its first finite Q; the other Mn3Sn runs start at 0.82),
# and the flat -<b>^2 stretch of scale_ft_rmc.gr pins b_avg_sq = 0.015407 and the
# 2.48 A enforcement cutoff. Kept in sync with tests/test_scaling.py.
STOG_59438 = StogInput(
    n_files=1, data_file="PG3_59438_SQ_rebin.dat", qmin=1.0, qmax=28.0,
    yoffset=-9.0, yscale=0.1, qoffset=0.0, out_sq="scale.fq", out_gr="scale.gr",
    rmax=50.0, nr=5000, lorch=False, rho0=0.063049, yoffset2=0.0,
    try_again=False, use_filter=True, r_cutoff=1.0, out_ft_sq="scale_ft.sq",
    out_ft_gr="scale_ft.gr", b_avg_sq=0.015407, out_rmc_fq="scale_ft_rmc.fq",
    out_rmc_gr="scale_ft_rmc.gr", out_rmc_dr="scale_ft_rmc.dr",
    peak_cutoff=2.48, peak_rmin=2.65, peak_rmax=3.1,
)

# Gitignored, so sample-backed tests skip rather than fail when the run is
# absent; the rebinned data file alone is enough (no stog.inp required).
requires_stog_run = unittest.skipUnless(
    (STOG_RUN / STOG_59438.data_file).exists(), "stog_59438 example run not present"
)

RHO0 = 0.05
B2 = 0.02


def synthetic_g(r: np.ndarray) -> np.ndarray:
    """Physically consistent g(r): zero below ~2.5 A, first shell at 2.8 A, -> 1."""
    onset = 0.5 * (1.0 + np.tanh((r - 2.65) / 0.07))
    peak = 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    return onset + peak


def synthetic_sq(q: np.ndarray, rho0: float = RHO0) -> np.ndarray:
    """Forward-model S(Q) from the synthetic g(r) on a dense r grid."""
    r = np.arange(1, 12001) * 0.005  # 0.005 .. 60 A
    gpdf = g_to_gpdf(r, synthetic_g(r), rho0)
    return fq_to_sq(q, gpdf_to_fq(r, gpdf, q))


class ConversionTests(unittest.TestCase):
    def setUp(self):
        self.q = np.linspace(0.5, 30.0, 500)
        self.r = np.linspace(0.1, 20.0, 400)
        self.sq = 1.0 + 0.3 * np.sin(self.q)
        self.g = 1.0 + 0.2 * np.cos(self.r)

    def test_sq_fq_round_trip(self):
        fq = sq_to_fq(self.q, self.sq)
        np.testing.assert_allclose(fq_to_sq(self.q, fq), self.sq, atol=1e-12)

    def test_sq_fk_round_trip(self):
        fk = sq_to_fk(self.sq, B2)
        np.testing.assert_allclose(fk / B2 + 1.0, self.sq, atol=1e-12)
        # Keen Eq. 19: FK = <b>^2 [S - 1]
        np.testing.assert_allclose(fk, B2 * (self.sq - 1.0), atol=1e-15)

    def test_g_gpdf_round_trip(self):
        gpdf = g_to_gpdf(self.r, self.g, RHO0)
        np.testing.assert_allclose(gpdf_to_g(self.r, gpdf, RHO0), self.g, atol=1e-12)

    def test_g_gk_round_trip_and_dr(self):
        gk = g_to_gk(self.g, B2)
        np.testing.assert_allclose(gk_to_g(gk, B2), self.g, atol=1e-12)
        # Keen Eq. 29: D(r) = 4 pi rho0 r GK(r)
        np.testing.assert_allclose(
            gk_to_dr(self.r, gk, RHO0),
            4.0 * np.pi * RHO0 * self.r * gk,
            atol=1e-15,
        )

    def test_density_line_matches_low_r_theory(self):
        # D(r < r0) = -4 pi rho0 <b>^2 r  (Keen Eqs. 29 + 15)
        line = density_line(self.r, RHO0, B2)
        gk_flat = np.full_like(self.r, -B2)
        np.testing.assert_allclose(line, gk_to_dr(self.r, gk_flat, RHO0), atol=1e-15)

    def test_lorch_window_limits(self):
        q = np.array([1e-12, 15.0, 30.0])
        window = lorch_window(q, 30.0)
        self.assertAlmostEqual(window[0], 1.0, places=9)
        self.assertAlmostEqual(window[1], np.sin(np.pi / 2) / (np.pi / 2), places=12)
        self.assertAlmostEqual(window[2], 0.0, places=12)

    def test_enforce_low_r(self):
        gk = np.linspace(-1.0, 1.0, self.r.size)
        enforced = enforce_low_r(self.r, gk, cutoff=5.0, b_avg_sq=B2)
        below = self.r <= 5.0
        np.testing.assert_allclose(enforced[below], -B2, atol=1e-15)
        np.testing.assert_allclose(enforced[~below], gk[~below], atol=1e-15)

    def test_first_peak_zero_fortran_semantics(self):
        # stog_new3.f90: zero g where r <= cutoff AND outside [rmin, rmax].
        r = np.array([0.5, 1.5, 2.7, 3.5, 4.5, 6.0])
        g = np.ones_like(r)
        out = first_peak_zero(r, g, cutoff=5.0, peak_rmin=2.0, peak_rmax=3.0)
        np.testing.assert_allclose(out, [0.0, 0.0, 1.0, 0.0, 0.0, 1.0])
        # Peak window entirely above cutoff (the 59438 layout: 2.48/2.65/3.1)
        # degenerates to zero-everything-below-cutoff == enforce_low_r.
        out = first_peak_zero(r, g, cutoff=2.48, peak_rmin=2.65, peak_rmax=3.1)
        np.testing.assert_allclose(out, [0.0, 0.0, 1.0, 1.0, 1.0, 1.0])

    def test_fourier_filter_rejects_nonpositive_q(self):
        q = np.arange(0, 300) * 0.1  # starts at Q = 0
        with self.assertRaises(ValueError):
            fourier_filter(q, np.ones_like(q), self.r, rho0=RHO0, cutoff=1.0)


class FourierTransformTests(unittest.TestCase):
    def test_round_trip_recovers_g(self):
        # r -> Q -> r with wide Q coverage recovers the synthetic model.
        q = np.arange(1, 3501) * 0.01  # 0.01 .. 35
        sq = synthetic_sq(q)
        r = np.linspace(0.5, 8.0, 751)
        g_back = gpdf_to_g(r, fq_to_gpdf(q, sq_to_fq(q, sq), r), RHO0)
        np.testing.assert_allclose(g_back, synthetic_g(r), atol=2e-2)

    def test_high_q_asymptote_and_low_r_zero(self):
        q = np.arange(1, 3501) * 0.01
        sq = synthetic_sq(q)
        # Keen Eq. 21: S(Q) -> 1 at high Q.
        self.assertLess(abs(np.mean(sq[q > 25.0]) - 1.0), 5e-3)
        # g(r) = 0 below the closest approach.
        r = np.linspace(1.0, 2.2, 121)
        g_low = gpdf_to_g(r, fq_to_gpdf(q, sq_to_fq(q, sq), r), RHO0)
        self.assertLess(np.sqrt(np.mean(g_low**2)), 3e-2)

    def test_low_q_correction_improves_truncated_transform(self):
        q_full = np.arange(1, 3501) * 0.01
        sq_full = synthetic_sq(q_full)
        keep = q_full >= 0.9  # truncate low-Q data, as measured datasets are
        q, sq = q_full[keep], sq_full[keep]
        r = np.linspace(1.0, 2.2, 121)
        fq = sq_to_fq(q, sq)
        g_plain = gpdf_to_g(r, fq_to_gpdf(q, fq, r), RHO0)
        g_corr = gpdf_to_g(
            r, fq_to_gpdf(q, fq, r, low_q_correction=True), RHO0
        )
        rms_plain = np.sqrt(np.mean(g_plain**2))
        rms_corr = np.sqrt(np.mean(g_corr**2))
        self.assertLess(rms_corr, rms_plain)

    def test_low_q_correction_basis_consistency(self):
        q = np.arange(90, 3000) * 0.01
        fq = np.sin(q)
        r = np.linspace(0.5, 5.0, 50)
        coef, const = low_q_correction_basis(q, r)
        s0 = fq[0] / q[0] + 1.0
        np.testing.assert_allclose(
            omitted_low_q_correction(q, fq, r), coef * s0 - const, atol=1e-14
        )

    def test_low_q_correction_zero_when_data_start_at_zero(self):
        # Data starting at Q = 0 omit nothing (pystog parity: F1 = F2 = 0).
        q = np.arange(0, 3000) * 0.01
        r = np.linspace(0.5, 5.0, 20)
        coef, const = low_q_correction_basis(q, r)
        np.testing.assert_allclose(coef, 0.0, atol=0.0)
        np.testing.assert_allclose(const, 0.0, atol=0.0)
        np.testing.assert_allclose(
            omitted_low_q_correction(q, np.sin(q), r), 0.0, atol=0.0
        )

    def test_low_q_correction_lorch_finite_at_singular_r(self):
        # The Lorch branch has a removable singularity at r = pi/Qmax.
        q = np.arange(90, 3000) * 0.01
        r_sing = np.pi / q[-1]
        r = np.array([r_sing - 0.01, r_sing, r_sing + 0.01])
        coef, const = low_q_correction_basis(q, r, lorch=True)
        self.assertTrue(np.all(np.isfinite(coef)))
        self.assertTrue(np.all(np.isfinite(const)))
        # The limit value sits between its neighbors (smooth function).
        self.assertLess(min(coef[0], coef[2]) - 1e-6, coef[1])
        self.assertGreater(max(coef[0], coef[2]) + 1e-6, coef[1])

    def test_fourier_filter_identity_and_cleanup(self):
        q = np.arange(50, 3001) * 0.01
        sq = synthetic_sq(q)
        # Contaminate with a smooth background that produces sub-cutoff artifacts.
        sq_dirty = sq + 0.4 * np.exp(-0.5 * (q / 4.0) ** 2)
        r = np.arange(1, 3001) * 0.01
        sq_filtered, sq_ft, g_filtered = fourier_filter(
            q, sq_dirty, r, rho0=RHO0, cutoff=1.5
        )
        # Defining relation of the classic stog ft.dat correction file.
        np.testing.assert_allclose(sq_filtered, sq_dirty - (sq_ft - 1.0), atol=1e-10)
        # The filter must reduce unphysical sub-cutoff structure.
        g_dirty = gpdf_to_g(r, fq_to_gpdf(q, sq_to_fq(q, sq_dirty), r), RHO0)
        below = r <= 1.5
        self.assertLess(
            np.sqrt(np.mean(g_filtered[below] ** 2)),
            np.sqrt(np.mean(g_dirty[below] ** 2)),
        )


@requires_stog_run
class StogRunParityTests(unittest.TestCase):
    """V1: reproduce the classic Fortran stog intermediates from the example run."""

    @classmethod
    def setUpClass(cls):
        cls.inp = STOG_59438
        data = read_stog_xy(STOG_RUN / cls.inp.data_file)
        q, sq = data[0], data[1]
        keep = np.isfinite(sq) & (q >= cls.inp.qmin - 1e-9) & (q <= cls.inp.qmax + 1e-9)
        cls.q, cls.sq_raw = q[keep], sq[keep]
        cls.sq_scaled = cls.inp.a * cls.sq_raw + cls.inp.b
        cls.r = np.arange(1, cls.inp.nr + 1) * (cls.inp.rmax / cls.inp.nr)

    def test_scaling_reproduces_scale_fq_exactly(self):
        ref = read_stog_xy(STOG_RUN / "scale.fq")
        ours = np.interp(ref[0], self.q, self.sq_scaled)
        np.testing.assert_allclose(ours, ref[1], atol=1e-9)

    def test_forward_transform_matches_scale_gr(self):
        # scale.gr stores g(r) - 1 (pinned by the pystog cross-run).
        ref = read_stog_xy(STOG_RUN / "scale.gr")
        gpdf = fq_to_gpdf(self.q, sq_to_fq(self.q, self.sq_scaled), self.r)
        ours = gpdf_to_g(self.r, gpdf, self.inp.rho0) - 1.0
        ref_on_ours = np.interp(self.r, ref[0], ref[1])
        rms = np.sqrt(np.mean((ours - ref_on_ours) ** 2))
        self.assertLess(rms / np.abs(ref_on_ours).max(), 5e-3)

    def test_fourier_filter_matches_ft_dat_and_scale_ft_sq(self):
        sq_filtered, sq_ft, _ = fourier_filter(
            self.q,
            self.sq_scaled,
            self.r,
            rho0=self.inp.rho0,
            cutoff=self.inp.r_cutoff,
        )
        ft_ref = read_stog_xy(STOG_RUN / "ft.dat")
        ours_ft = np.interp(ft_ref[0], self.q, sq_ft)
        keep = (ft_ref[0] >= self.inp.qmin) & (ft_ref[0] <= self.inp.qmax)
        rms_ft = np.sqrt(np.mean((ours_ft[keep] - ft_ref[1][keep]) ** 2))
        self.assertLess(rms_ft, 2e-3)

        sq_ref = read_stog_xy(STOG_RUN / "scale_ft.sq")
        ours_sq = np.interp(sq_ref[0], self.q, sq_filtered)
        rms_sq = np.sqrt(np.mean((ours_sq - sq_ref[1]) ** 2))
        self.assertLess(rms_sq, 2e-3)

    def test_enforced_rmc_outputs_match_theory_and_file(self):
        _, _, g_filtered = fourier_filter(
            self.q,
            self.sq_scaled,
            self.r,
            rho0=self.inp.rho0,
            cutoff=self.inp.r_cutoff,
        )
        gk = g_to_gk(g_filtered, self.inp.b_avg_sq)
        gk_enf = enforce_low_r(
            self.r, gk, cutoff=self.inp.peak_cutoff, b_avg_sq=self.inp.b_avg_sq
        )
        dr_enf = gk_to_dr(self.r, gk_enf, self.inp.rho0)

        gr_ref = read_stog_xy(STOG_RUN / "scale_ft_rmc.gr")
        dr_ref = read_stog_xy(STOG_RUN / "scale_ft_rmc.dr")
        below = self.r <= self.inp.peak_cutoff
        # Below the cutoff the Fortran files hold the exact theoretical values.
        np.testing.assert_allclose(
            gk_enf[below], np.interp(self.r[below], gr_ref[0], gr_ref[1]), atol=1e-9
        )
        np.testing.assert_allclose(
            dr_enf[below], np.interp(self.r[below], dr_ref[0], dr_ref[1]), atol=1e-9
        )
        # And they follow Keen Eqs. 15 / 29+15 by construction.
        np.testing.assert_allclose(gk_enf[below], -self.inp.b_avg_sq, atol=1e-12)
        np.testing.assert_allclose(
            dr_enf[below],
            density_line(self.r[below], self.inp.rho0, self.inp.b_avg_sq),
            atol=1e-12,
        )


if __name__ == "__main__":
    unittest.main()
