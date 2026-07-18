# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import (
    StogInput,
    read_dat_header,
    read_stog_inp,
    read_stog_xy,
    write_stog_xy,
)
from rmc_toolkits.scaling import (
    ScalingConfig,
    autoscale,
    detect_first_peak_onset,
    diagnostics_summary,
    level_sweep,
    scale_pipeline,
)
from rmc_toolkits.scattering import (
    faber_ziman,
    mass_density_from_number_density,
    number_density_from_mass_density,
)
from rmc_toolkits.transforms import first_peak_zero, fq_to_sq, g_to_gpdf, gpdf_to_fq

ROOT = Path(__file__).resolve().parents[1]
STOG_RUN = ROOT / "data" / "stog_tests" / "stog_59438"


def _first_existing_run(*names: str) -> Path:
    candidates = [ROOT / "data" / "stog_tests" / name for name in names]
    for candidate in candidates:
        if (candidate / "stog_input.dat").exists():
            return candidate
    return candidates[0]


# The FeCoSn x-ray series shares one stog parameterization across temperatures
# (same 0.5-26 window, hand scaling -0.111/0.9, rho0, and 1.0 0 0 enforcement),
# so any available temperature exercises the same assertions.
XRAY_RUN = _first_existing_run("100K", "199K")

# stog_59438: a complete Fortran stog run of Mn3Sn (POWGEN PG3), local-only. Its
# parameters are recovered from the run's own outputs (scale.fq vs the rebinned
# input pins a = 10, b = -9 and Qmin = 1.0; the flat -<b>^2 stretch of
# scale_ft_rmc.gr pins b_avg_sq = 0.015407 and the 2.48 A enforcement cutoff),
# so the sample-backed tests need only the rebinned data file, not the
# interactive stog.inp. Kept in sync with tests/test_transforms.py.
STOG_59438 = StogInput(
    n_files=1, data_file="PG3_59438_SQ_rebin.dat", qmin=1.0, qmax=28.0,
    yoffset=-9.0, yscale=0.1, qoffset=0.0, out_sq="scale.fq", out_gr="scale.gr",
    rmax=50.0, nr=5000, lorch=False, rho0=0.063049, yoffset2=0.0,
    try_again=False, use_filter=True, r_cutoff=1.0, out_ft_sq="scale_ft.sq",
    out_ft_gr="scale_ft.gr", b_avg_sq=0.015407, out_rmc_fq="scale_ft_rmc.fq",
    out_rmc_gr="scale_ft_rmc.gr", out_rmc_dr="scale_ft_rmc.dr",
    peak_cutoff=2.48, peak_rmin=2.65, peak_rmax=3.1,
)

requires_stog_run = unittest.skipUnless(
    (STOG_RUN / STOG_59438.data_file).exists(), "stog_59438 example run not present"
)
requires_xray_run = unittest.skipUnless(
    (XRAY_RUN / "stog_input.dat").exists(), "FeCoSn x-ray run not present"
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
        # Scale-only fitting requires joint mode: in the default sweep mode
        # the offset is tied to the measured level (b = 1 - a*L), not frozen.
        result = autoscale(
            self.q, sq_meas, synthetic_config(c1_mode="joint", fit_offset=False)
        )
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


class LevelSweepTests(unittest.TestCase):
    """The criterion-driven two-sided high-Q window search."""

    @classmethod
    def setUpClass(cls):
        cls.q = np.arange(60, 3001) * 0.01
        cls.sq_true = synthetic_sq(cls.q)

    def test_recovers_known_level(self):
        # Corruption (S - b)/a has asymptote (1 - b)/a.
        a_t, b_t = 4.0, -2.0
        sq_meas = (self.sq_true - b_t) / a_t
        result = level_sweep(self.q, sq_meas)
        self.assertTrue(result.asymptote_found)
        expected = (1.0 - b_t) / a_t  # 0.75
        self.assertLess(abs(result.level - expected), 0.01)
        self.assertGreater(result.n_admissible, 10)

    def test_excludes_dead_tail(self):
        # Detector rolloff: the tail dies above Q = 24 — any admissible
        # window must end below the artifact.
        sq_meas = self.sq_true.copy()
        bad = self.q > 24.0
        sq_meas[bad] *= np.exp(-(self.q[bad] - 24.0) / 3.0)
        result = level_sweep(self.q, sq_meas)
        self.assertTrue(result.asymptote_found)
        self.assertLessEqual(result.q_hi, 24.5)
        self.assertLess(abs(result.level - 1.0), 0.02)

    def test_reports_no_asymptote_on_sloped_data(self):
        # A uniform slope everywhere: no statistically flat window exists.
        rng = np.random.default_rng(3)
        sq_meas = 1.0 + 0.05 * self.q + rng.normal(0, 1e-4, self.q.size)
        result = level_sweep(self.q, sq_meas)
        self.assertFalse(result.asymptote_found)
        self.assertEqual(result.n_admissible, 0)

    def test_sweep_mode_recovery_and_provenance(self):
        a_t, b_t = 10.0, -9.0
        sq_meas = (self.sq_true - b_t) / a_t
        result = autoscale(self.q, sq_meas, synthetic_config())
        self.assertEqual(result.provenance["c1_mode_effective"], "sweep")
        self.assertIsNotNone(result.sweep)
        self.assertLess(abs(result.a - a_t) / a_t, 0.02)
        # b is tied to the level: b = 1 - a * L.
        self.assertAlmostEqual(result.b, 1.0 - result.a * result.sweep.level, places=10)

    def test_fz_amplitude_mode_round_trip(self):
        # amplitude_criterion="fz": subtract the measured level, pin Q->0 on
        # the Faber-Ziman limit, shift the level back to 1. Choosing b_sq_avg
        # so the FZ target equals the synthetic's own extrapolated S(0) makes
        # this a true round trip (affine corruption cancels exactly).
        head = self.q <= self.q[0] + 1.0
        _, s_true_0 = np.polyfit(self.q[head], self.sq_true[head], 1)
        config = synthetic_config(
            b_sq_avg=B2 * (1.0 - float(s_true_0)), amplitude_criterion="fz"
        )
        a_t, b_t = 3.0, -1.5
        result = autoscale(self.q, (self.sq_true - b_t) / a_t, config)
        self.assertEqual(result.iterations, 0)  # closed form, no loop
        self.assertEqual(result.provenance["config"]["amplitude_criterion"], "fz")
        self.assertLess(abs(result.a - a_t) / a_t, 0.05)
        self.assertAlmostEqual(result.b, 1.0 - result.a * result.sweep.level, places=10)
        summary = diagnostics_summary(result, config)
        self.assertIn("a_fz", summary)
        # a IS a_fz here: no self-concordance is reported.
        self.assertNotIn("amplitude_concordance", summary)

    def test_fz_mode_config_validation(self):
        with self.assertRaises(ValueError):
            synthetic_config(amplitude_criterion="fz")  # b_sq_avg missing
        with self.assertRaises(ValueError):
            synthetic_config(
                amplitude_criterion="fz", b_sq_avg=B2, c1_mode="joint"
            )
        with self.assertRaises(ValueError):
            synthetic_config(amplitude_criterion="nope")

    def test_fz_amplitude_and_concordance_mechanics(self):
        # With b_sq_avg supplied, the independent Q->0 amplitude and the
        # concordance diagnostic are reported. The synthetic model is not
        # compressibility-consistent, so we assert mechanics, not agreement.
        a_t, b_t = 2.0, -1.0
        sq_meas = (self.sq_true - b_t) / a_t
        config = synthetic_config(b_sq_avg=B2)  # monatomic-like: s0_target = 0
        result = autoscale(self.q, sq_meas, config)
        self.assertIsNotNone(result.a_fz)
        summary = diagnostics_summary(result, config)
        self.assertIn("amplitude_concordance", summary)
        self.assertIn("amplitudes_concordant", summary)
        self.assertGreater(summary["a_fz"], 0.0)


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
        cls.inp = STOG_59438
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
class XrayFeCoSnTests(unittest.TestCase):
    """FeCoSn x-ray run (any available temperature): normalized S(Q), Qmin = 0.5.

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


#: Complete Fortran stog runs of Mn3Sn (POWGEN PG3), local-only. Parameters
#: recovered from the outputs themselves (max residual ~1e-12): hand (a, b)
#: from scale.fq vs the rebinned input, the enforcement cutoff from the flat
#: -<b>^2 stretch of scale_ft_rmc.gr.
MN3SN_RUNS = {
    "stog": ("PG3_55537_rebin.sq", 2.5, -1.5, 2.60),
    "stog_300K": ("PG3_55526_SQ_rebin.sq", 2.050021, -1.05, 2.68),
    "stog_500K": ("PG3_54139_SQ_rebin.dat", 10.0, -9.0, 2.40),
}
MN3SN_ROOT = ROOT / "data" / "stog_tests"

requires_mn3sn = unittest.skipUnless(
    all((MN3SN_ROOT / name / raw).exists() for name, (raw, _, _, _) in MN3SN_RUNS.items()),
    "Mn3Sn PG3 neutron runs not present",
)


@requires_mn3sn
class Mn3SnNeutronTests(unittest.TestCase):
    """Negative-b neutron material: <b^2>/<b>^2 = 13.06, S(0) = -12.06.

    The composition carries the absolute-scale information here: the density
    limit alone is degenerate (one-sided flag False on every run) and the
    historical hand scalings are mutually inconsistent (2.5x / 2.05x / 10x
    for the same material). Classic-parity comparisons keep b_sq_avg unset so
    the omitted-low-Q correction stays at the Fortran's S(0) = 0 convention.
    """

    QMIN, QMAX, RHO0, B2 = 0.82, 28.0, 0.063049, 0.015407

    @classmethod
    def setUpClass(cls):
        cls.fz = faber_ziman("Mn3Sn")
        cls.data = {
            name: read_stog_xy(MN3SN_ROOT / name / raw)
            for name, (raw, _, _, _) in MN3SN_RUNS.items()
        }

    def parity_config(self):
        return ScalingConfig(
            qmin=self.QMIN, qmax=self.QMAX, rho0=self.RHO0, b_avg_sq=self.B2
        )

    def composition_config(self, **overrides):
        return ScalingConfig(
            qmin=self.QMIN, qmax=self.QMAX, rho0=self.RHO0,
            b_avg_sq=self.fz.b_avg_sq_barn, b_sq_avg=self.fz.b_sq_avg_barn,
            **overrides,
        )

    def test_composition_reproduces_the_stog_coefficients(self):
        self.assertAlmostEqual(self.fz.b_avg_sq_barn, 0.015407, places=6)
        self.assertAlmostEqual(self.fz.b_sq_avg_barn, 0.201223, places=5)
        # ADDIE-convention mass <-> number density round trip.
        mass_density = mass_density_from_number_density("Mn3Sn", self.RHO0)
        self.assertAlmostEqual(mass_density, 7.4209, places=3)
        self.assertAlmostEqual(
            number_density_from_mass_density("Mn3Sn", mass_density), self.RHO0, places=9
        )

    def test_manual_parity_with_fortran(self):
        for name, (_, a, b, cutoff) in MN3SN_RUNS.items():
            with self.subTest(run=name):
                q, sq = self.data[name][0], self.data[name][1]
                manual = scale_pipeline(q, sq, self.parity_config(), a, b)

                ref = read_stog_xy(MN3SN_ROOT / name / "scale_ft.sq")
                ours = np.interp(ref[0], manual.q, manual.sq_filtered)
                rms = float(np.sqrt(np.mean((ours - ref[1]) ** 2)))
                self.assertLess(rms, 2e-3)

                gr = read_stog_xy(MN3SN_ROOT / name / "scale_ft_rmc.gr")
                enforced = first_peak_zero(
                    manual.r, manual.g_filtered,
                    cutoff=cutoff, peak_rmin=0.0, peak_rmax=0.0,
                )
                gk = self.B2 * (enforced - 1.0)
                ours_gr = np.interp(gr[0], manual.r, gk)
                rms_gr = float(np.sqrt(np.mean((ours_gr - gr[1]) ** 2)))
                self.assertLess(rms_gr, 5e-3)
                below = gr[0] <= cutoff - 0.02
                np.testing.assert_allclose(gr[1][below], -self.B2, atol=1e-6)

    def test_autoscale_composition_only_detects_first_shell(self):
        q, sq = self.data["stog_300K"][0], self.data["stog_300K"][1]
        result = autoscale(q, sq, self.composition_config())
        summary = diagnostics_summary(result, self.composition_config())
        self.assertTrue(result.converged)
        # The inverted Mn-Sn first shell sits at ~2.7-2.8 A; detection must
        # land on its left flank and refine the fit window there.
        self.assertIsNotNone(summary.get("r0_detected"))
        self.assertGreater(summary["r0_detected"], 2.4)
        self.assertLess(summary["r0_detected"], 3.0)
        self.assertTrue(summary.get("window_refined"))
        # Honest one-sided verdict: the low-Q hole is O(S(0)) = O(-12) deep,
        # so the absolute scale is NOT recoverable from self-consistency.
        self.assertFalse(summary["density_limit_satisfied"])

    def test_fz_amplitude_uses_the_composition(self):
        # The FZ criterion injects S(0) = 1 - <b^2>/<b>^2 = -12.06 and lands
        # at O(10) amplitudes — the only route consistent across runs (the
        # density-limit amplitudes sit at O(1) and the hand values disagree
        # with each other by 5x).
        for name in MN3SN_RUNS:
            with self.subTest(run=name):
                q, sq = self.data[name][0], self.data[name][1]
                result = autoscale(
                    q, sq, self.composition_config(amplitude_criterion="fz")
                )
                self.assertGreater(result.a, 5.0)
                self.assertLess(result.a, 25.0)


class DetectionAndDensityTests(unittest.TestCase):
    """Synthetic checks for the r0 detector and the s0-aware correction."""

    @classmethod
    def setUpClass(cls):
        cls.q = np.arange(60, 3001) * 0.01
        cls.sq_true = synthetic_sq(cls.q)

    def test_detects_first_shell_and_refines_window(self):
        config = synthetic_config(r0=None, r_fit_max=None)
        result = autoscale(self.q, (self.sq_true + 4.0) / 5.0, config)
        summary = diagnostics_summary(result, config)
        onset = summary.get("r0_detected")
        self.assertIsNotNone(onset)
        # Synthetic first shell: onset 2.65, peak 2.8.
        self.assertGreater(onset, 2.3)
        self.assertLess(onset, 2.85)
        self.assertLess(abs(result.a - 5.0) / 5.0, 0.03)

    def test_detector_ignores_subcutoff_ripples(self):
        rng = np.random.default_rng(11)
        r = np.arange(1, 1001) * 0.01
        g = 0.3 * np.sin(r * 28.0) * rng.uniform(0.5, 1.0, r.size)
        peak = 3.0 * np.exp(-0.5 * ((r - 2.8) / 0.12) ** 2)
        onset = detect_first_peak_onset(r, g + peak, 28.0, search_min=1.3)
        self.assertIsNotNone(onset)
        self.assertGreater(onset, 2.3)
        self.assertLess(onset, 2.8)

    def test_s0_target_changes_only_the_constant_term(self):
        from rmc_toolkits.transforms import low_q_correction_basis

        q = np.arange(82, 2801) * 0.01
        r = np.arange(1, 501) * 0.01
        coef0, const0 = low_q_correction_basis(q, r)
        coef1, const1 = low_q_correction_basis(q, r, s0_target=-12.06)
        np.testing.assert_allclose(coef1, coef0, atol=1e-15)
        np.testing.assert_allclose(
            const1, (1.0 - (-12.06)) * const0 + (-12.06) * coef0, atol=1e-12
        )

    def test_effective_s0_target_resolution(self):
        base = synthetic_config()
        self.assertEqual(base.effective_s0_target, 0.0)
        with_b2 = synthetic_config(b_sq_avg=2.0 * B2)
        self.assertAlmostEqual(with_b2.effective_s0_target, -1.0)
        pinned = synthetic_config(b_sq_avg=2.0 * B2, s0_target=0.25)
        self.assertAlmostEqual(pinned.effective_s0_target, 0.25)


if __name__ == "__main__":
    unittest.main()
