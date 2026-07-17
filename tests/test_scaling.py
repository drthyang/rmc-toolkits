# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import (
    read_dat_header,
    read_stog_inp,
    read_stog_xy,
    write_stog_xy,
)
from rmc_toolkits.scaling import (
    ScalingConfig,
    autoscale,
    diagnostics_summary,
    scale_pipeline,
)
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

ROOT = Path(__file__).resolve().parents[1]
STOG_RUN = ROOT / "data" / "stog_tests" / "stog_59438"
XRAY_RUN = ROOT / "data" / "stog_tests" / "100K"

requires_stog_run = unittest.skipUnless(
    (STOG_RUN / "stog.inp").exists(), "stog_59438 example run not present"
)
requires_xray_run = unittest.skipUnless(
    (XRAY_RUN / "stog_input.dat").exists(), "FeCoSn 100K x-ray run not present"
)

RHO0 = 0.05
B2 = 0.02


def synthetic_g(r: np.ndarray) -> np.ndarray:
    """Same synthetic model as tests/test_transforms.py (kept in sync)."""
    onset = 0.5 * (1.0 + np.tanh((r - 2.65) / 0.07))
    peak = 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    return onset + peak


def synthetic_sq(q: np.ndarray, rho0: float = RHO0) -> np.ndarray:
    r = np.arange(1, 12001) * 0.005
    gpdf = g_to_gpdf(r, synthetic_g(r), rho0)
    return fq_to_sq(q, gpdf_to_fq(r, gpdf, q))


def synthetic_config(**overrides) -> ScalingConfig:
    values = dict(
        qmin=0.6,
        qmax=30.0,
        rho0=RHO0,
        b_avg_sq=B2,
        r_cutoff=1.0,
        r0=2.5,
        r_fit_min=1.2,
        r_fit_max=2.25,
        rmax=30.0,
        nr=3000,
    )
    values.update(overrides)
    return ScalingConfig(**values)


class AutoscaleSyntheticTests(unittest.TestCase):
    """V4: known-scale recovery and physics limits on synthetic data."""

    @classmethod
    def setUpClass(cls):
        cls.q = np.arange(60, 3001) * 0.01  # 0.60 .. 30.00
        cls.sq_true = synthetic_sq(cls.q)

    def test_recovers_known_scale_and_offset(self):
        a_true, b_true = 10.0, -9.0
        sq_meas = (self.sq_true - b_true) / a_true
        result = autoscale(self.q, sq_meas, synthetic_config())
        self.assertTrue(result.converged)
        # The omitted-low-Q correction (on by default) is what makes this
        # tight: without it the Qmin=0.6 truncation biases a by ~8%.
        self.assertLess(abs(result.a - a_true) / a_true, 0.02)
        self.assertLess(abs(result.b - b_true) / abs(b_true), 0.02)

    def test_recovers_scale_with_offset_disabled(self):
        a_true = 2.5
        sq_meas = self.sq_true / a_true
        # With b frozen at 0 the model is exact, so recovery is tight.
        result = autoscale(self.q, sq_meas, synthetic_config(fit_offset=False))
        self.assertEqual(result.b, 0.0)
        self.assertLess(abs(result.a - a_true) / a_true, 0.02)

    def test_result_satisfies_physics_limits(self):
        sq_meas = (self.sq_true + 4.0) / 5.0
        result = autoscale(self.q, sq_meas, synthetic_config())
        # C1 (Keen Eq. 21): filtered S(Q) tail oscillates about 1.
        self.assertLess(abs(result.c1_tail_mean - 1.0), 0.02)
        # C2 (Keen Eq. 15 in g-space): g ~ 0 below the first peak.
        self.assertLess(result.low_r_rms, 0.05)
        # GK sits on the -<b>^2 level in the window (loose: ripples remain).
        window = (result.r >= 1.2) & (result.r <= 2.25)
        self.assertLess(abs(result.gk[window].mean() + B2), 0.15 * B2)

    def test_noise_robustness(self):
        rng = np.random.default_rng(1234)
        sq_meas = (self.sq_true - (-2.0)) / 4.0 + rng.normal(0.0, 0.01, self.q.size)
        result = autoscale(self.q, sq_meas, synthetic_config())
        self.assertTrue(result.converged)
        self.assertLess(abs(result.a - 4.0) / 4.0, 0.1)

    def test_sigma_weighting_accepted(self):
        rng = np.random.default_rng(5)
        sq_meas = (self.sq_true + 1.0) / 2.0 + rng.normal(0.0, 0.005, self.q.size)
        sigma = np.full(self.q.size, 0.005)
        result = autoscale(self.q, sq_meas, synthetic_config(), sigma=sigma)
        self.assertTrue(result.converged)
        self.assertLess(abs(result.a - 2.0) / 2.0, 0.05)

    def test_despike_restores_recovery_under_tail_glitches(self):
        # Detector-glitch spikes ring through the transform into the C2 window
        # — a channel Huber IRLS cannot reject (measured: ~80% scale error).
        # Opt-in despiking removes them before any transform and restores
        # clean recovery. Despike stays OFF by default because it also flags
        # real Bragg maxima on crystalline data (12% of the 59438 points).
        rng = np.random.default_rng(7)
        sq_meas = (self.sq_true + 9.0) / 10.0
        idx = rng.choice(np.where(self.q > 25.5)[0], 8, replace=False)
        sq_meas = sq_meas.copy()
        sq_meas[idx] += 2.0

        spiked = autoscale(self.q, sq_meas, synthetic_config())
        despiked = autoscale(self.q, sq_meas, synthetic_config(despike=True))
        self.assertLess(abs(despiked.a - 10.0) / 10.0, 0.02)
        self.assertLess(
            abs(despiked.a - 10.0), 0.2 * abs(spiked.a - 10.0)
        )
        self.assertGreaterEqual(despiked.provenance["n_despiked"], 8)

    def test_slope_nuisance_config_accepted(self):
        result = autoscale(
            self.q, (self.sq_true + 1.0) / 2.0,
            synthetic_config(c1_slope_nuisance=True),
        )
        self.assertTrue(result.converged)
        self.assertLess(abs(result.a - 2.0) / 2.0, 0.02)

    def test_history_and_provenance(self):
        result = autoscale(self.q, (self.sq_true + 1.0) / 2.0, synthetic_config())
        self.assertGreaterEqual(len(result.history), 1)
        self.assertEqual(result.provenance["mode"], "auto")
        self.assertEqual(result.provenance["model"], "S_corr = a*S_meas + b")
        summary = diagnostics_summary(result, synthetic_config())
        self.assertIn("low_r_rms_pre_enforcement", summary)

    def test_enforcement_stage(self):
        config = synthetic_config(enforce_cutoff=2.4)
        result = scale_pipeline(self.q, self.sq_true, config, 1.0, 0.0)
        below = result.r <= 2.4
        np.testing.assert_allclose(result.gk_enforced[below], -B2, atol=1e-12)
        # Above the cutoff the enforced array equals the plain one.
        np.testing.assert_allclose(
            result.gk_enforced[~below], result.gk[~below], atol=1e-15
        )

    def test_rejects_bad_config(self):
        with self.assertRaises(ValueError):
            synthetic_config(rho0=-1.0)
        with self.assertRaises(ValueError):
            synthetic_config(qmin=30.0, qmax=1.0)
        with self.assertRaises(ValueError):
            synthetic_config(r_fit_min=2.0, r_fit_max=1.0).r_fit_window
        with self.assertRaises(ValueError):
            synthetic_config(nr=0)
        with self.assertRaises(ValueError):
            synthetic_config(rmax=-5.0)

    def test_q_zero_rows_are_dropped(self):
        # A grid containing Q = 0 must not poison the pipeline with NaN.
        q = np.concatenate([[0.0], self.q])
        sq = np.concatenate([[1.0], self.sq_true])
        result = scale_pipeline(q, sq, synthetic_config(), 1.0, 0.0)
        self.assertGreater(result.q[0], 0.0)
        self.assertTrue(np.all(np.isfinite(result.gk)))
        self.assertTrue(np.all(np.isfinite(result.sq_filtered)))


class StogFileIoTests(unittest.TestCase):
    def test_read_stog_xy_skips_headers_and_keeps_nan(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sample.dat"
            path.write_text(
                "        3401\n0.956998\n"
                "  0.50 NaN\n  0.51 1.25\n  0.52 1.30\nsome title text\n"
            )
            data = read_stog_xy(path)
            self.assertEqual(data.shape, (2, 3))
            self.assertTrue(np.isnan(data[1][0]))
            self.assertAlmostEqual(data[1][2], 1.30)

    def test_write_read_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "out.sq"
            x = np.linspace(1.0, 2.0, 11)
            y = np.sin(x)
            write_stog_xy(path, x, y, title="rmc S(Q) test")
            data = read_stog_xy(path)
            np.testing.assert_allclose(data[0], x, atol=1e-15)
            np.testing.assert_allclose(data[1], y, atol=1e-15)

    def test_read_dat_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "meta.dat"
            path.write_text(
                "TITLE :: GaTa4Se8 250K\n"
                "NUMBER_DENSITY :: 0.046690 Angstrom^(-3)\n"
                "MINIMUM_DISTANCES :: 2.4 2.2 2.6 Angstrom\n\n1.0 2.0\n"
            )
            header = read_dat_header(path)
            self.assertEqual(header["title"], "GaTa4Se8 250K")
            self.assertAlmostEqual(header["number_density"], 0.046690)
            self.assertAlmostEqual(header["min_distance"], 2.2)

    def test_read_stog_inp_variants(self):
        base = (
            "1\ndata.dat\n1.0 28.0\n{scale_line}\n0\nscale.fq\nscale.gr\n50\n5000\nN\n"
            "0.063049\n0\n{try_again}\nY\n1.0\nscale_ft.sq\nscale_ft.gr\n0.015407\n"
            "rmc.fq\nrmc.gr\nrmc.dr\n2.48 2.65 3.1\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            good = Path(tmp) / "stog.inp"
            good.write_text(base.format(scale_line="-9 0.1", try_again="N"))
            inp = read_stog_inp(good)
            self.assertAlmostEqual(inp.a, 10.0)
            self.assertAlmostEqual(inp.b, -9.0)
            self.assertFalse(inp.lorch)

            bad = Path(tmp) / "loop.inp"
            bad.write_text(base.format(scale_line="-9 0.1", try_again="Y"))
            with self.assertRaises(NotImplementedError):
                read_stog_inp(bad)

            zero_scale = Path(tmp) / "zero.inp"
            zero_scale.write_text(base.format(scale_line="-9 0", try_again="N"))
            with self.assertRaises(ValueError):
                read_stog_inp(zero_scale)

    def test_read_stog_xy_modal_columns(self):
        # A two-token numeric header must not become the column template.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hdr.dat"
            path.write_text(
                "3401 0.956998\n  0.51 1.25 0.01\n  0.52 1.30 0.01\n  0.53 1.28 0.01\n"
            )
            data = read_stog_xy(path)
            self.assertEqual(data.shape, (3, 3))
            self.assertAlmostEqual(data[0][0], 0.51)


@requires_stog_run
class AutoscaleExampleRunTests(unittest.TestCase):
    """V3: the auto-fit must do at least as well as the expert's hand tuning."""

    @classmethod
    def setUpClass(cls):
        cls.inp = read_stog_inp(STOG_RUN / "stog.inp")
        data = read_stog_xy(STOG_RUN / cls.inp.data_file)
        cls.q, cls.sq = data[0], data[1]
        cls.config = ScalingConfig(
            qmin=cls.inp.qmin,
            qmax=cls.inp.qmax,
            rho0=cls.inp.rho0,
            b_avg_sq=cls.inp.b_avg_sq,
            r_cutoff=cls.inp.r_cutoff,
            r0=cls.inp.peak_cutoff,
            r_fit_min=1.2,
            r_fit_max=cls.inp.peak_cutoff - 0.25,
            rmax=cls.inp.rmax,
            nr=cls.inp.nr,
            enforce_cutoff=cls.inp.peak_cutoff,
        )

    def test_autoscale_beats_hand_tuning(self):
        manual = scale_pipeline(self.q, self.sq, self.config, self.inp.a, self.inp.b)
        auto = autoscale(self.q, self.sq, self.config)
        self.assertTrue(auto.converged)
        self.assertLessEqual(auto.low_r_rms, manual.low_r_rms * 1.001)
        self.assertGreater(auto.a, 0.0)
        self.assertLess(abs(auto.c1_tail_mean - 1.0), 0.05)
        # This dataset misses O(1) structure below Qmin = 1.0 (its filtered
        # outputs violate the Krogh-Moe sum rule 26x even at the expert's hand
        # scaling), so the absolute scale is NOT recoverable from
        # self-consistency. The diagnostic is ONE-SIDED: False proves
        # non-recoverability (this case); True would not certify the scale.
        summary = diagnostics_summary(auto, self.config)
        self.assertFalse(summary["density_limit_satisfied"])

    def test_enforced_outputs_ready_for_rmcprofile(self):
        auto = autoscale(self.q, self.sq, self.config)
        below = auto.r <= self.inp.peak_cutoff
        np.testing.assert_allclose(
            auto.gk_enforced[below], -self.inp.b_avg_sq, atol=1e-12
        )


@requires_xray_run
class Xray100KTests(unittest.TestCase):
    """FeCoSn 100 K x-ray run: normalized S(Q) (<b>^2 = 1), Qmin = 0.5.

    Unlike the neutron 59438 case, this dataset CAN satisfy the density limit
    (low enough Qmin), so the auto-scaler should land near the expert's hand
    values and the one-sided diagnostic should report True.
    """

    @classmethod
    def setUpClass(cls):
        cls.inp = read_stog_inp(XRAY_RUN / "stog_input.dat")
        data = read_stog_xy(XRAY_RUN / cls.inp.data_file)
        cls.q, cls.sq = data[0], data[1]
        cls.config = ScalingConfig(
            qmin=cls.inp.qmin,
            qmax=cls.inp.qmax,
            rho0=cls.inp.rho0,
            b_avg_sq=cls.inp.b_avg_sq,
            r_cutoff=cls.inp.r_cutoff,
            r0=2.4,
            r_fit_min=1.2,
            r_fit_max=2.15,
            rmax=cls.inp.rmax,
            nr=cls.inp.nr,
            enforce_cutoff=cls.inp.peak_cutoff,
        )

    def test_inp_decodes_xray_conventions(self):
        self.assertAlmostEqual(self.inp.a, 1.0 / 0.9)
        self.assertAlmostEqual(self.inp.b, -0.111)
        self.assertAlmostEqual(self.inp.b_avg_sq, 1.0)  # normalized x-ray S(Q)
        self.assertAlmostEqual(self.inp.rho0, 0.057329)
        self.assertAlmostEqual(self.inp.peak_cutoff, 1.0)

    def test_manual_parity_with_fortran(self):
        manual = scale_pipeline(self.q, self.sq, self.config, self.inp.a, self.inp.b)
        ref = read_stog_xy(XRAY_RUN / "scale.fq")
        ours = np.interp(ref[0], manual.q, manual.sq_scaled)
        np.testing.assert_allclose(ours, ref[1], atol=1e-9)

        ft = read_stog_xy(XRAY_RUN / "ft.dat")
        keep = (ft[0] >= self.inp.qmin) & (ft[0] <= self.inp.qmax)
        ours_ft = np.interp(ft[0][keep], manual.q, manual.sq_ft)
        self.assertLess(
            float(np.sqrt(np.mean((ours_ft - ft[1][keep]) ** 2))), 5e-5
        )

        gr = read_stog_xy(XRAY_RUN / "scale_ft_rmc.gr")
        ours_gr = np.interp(gr[0], manual.r, manual.gk_enforced)
        self.assertLess(float(np.sqrt(np.mean((ours_gr - gr[1]) ** 2))), 2e-3)
        # x-ray normalized flat level: -<b>^2 = -1 below the 1.0 A cutoff.
        below = manual.r <= self.inp.peak_cutoff
        np.testing.assert_allclose(manual.gk_enforced[below], -1.0, atol=1e-12)

    def test_autoscale_succeeds_on_satisfiable_data(self):
        manual = scale_pipeline(self.q, self.sq, self.config, self.inp.a, self.inp.b)
        auto = autoscale(self.q, self.sq, self.config)
        self.assertTrue(auto.converged)
        self.assertLessEqual(auto.low_r_rms, manual.low_r_rms * 1.001)
        # Lands in the expert's neighborhood (hand a = 1.111): the density
        # limit is satisfiable here, unlike the neutron 59438 case.
        self.assertLess(abs(auto.a - self.inp.a) / self.inp.a, 0.15)
        self.assertLess(abs(auto.c1_tail_mean - 1.0), 0.02)
        summary = diagnostics_summary(auto, self.config)
        self.assertTrue(summary["density_limit_satisfied"])


if __name__ == "__main__":
    unittest.main()
