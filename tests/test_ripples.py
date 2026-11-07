# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Fourier termination ripples (rmc_toolkits/ripples.py): known answers.

1. How RMCProfile makes the model G(r) of a real-space fit, from its own files:
   the committed Demo (an X-ray F(Q) fitted in real space: the model column
   is the same transform of the model F(Q)) and RMCProfile's output on a
   synthetic configuration (tests/fixtures/rmcprofile_model_gr.json: a neutron
   G(r) without CONVOLVE is the histogram; with CONVOLVE the sinc in
   D-space, the baseline below R_CUTOFF; an X-ray G(r) the rectangle-rule
   transform on its QMIN–QMAX grid).
2. Synthetic G(r)s of Gaussian peaks, transformed to F(Q) exactly, cut at
   Q_max and transformed back with RMCProfile's rectangle rule: the predicted
   ripples are the actual ones; a planted shoulder above the resolution is
   not a ripple; a ripple is ripple-consistent; two Q_max separate a ripple
   from a feature; Lorch suppresses as predicted; the low-r level is the
   planted noise and ripple.
3. The Demo end to end.
"""

from __future__ import annotations

import json
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

from rmc_toolkits import ripples as R
from rmc_toolkits.ft_uncertainty import rmc_ft
from rmc_toolkits.parsers import read_rmc_csv

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"
EVIDENCE = Path(__file__).resolve().parent / "fixtures" / "rmcprofile_model_gr.json"

DQ = 0.01
R_GRID = np.arange(1, 1201) * 0.01
PEAKS = (
    {"amplitude": 8.0, "r": 2.5, "sigma": 0.06},
    {"amplitude": 6.0, "r": 3.6, "sigma": 0.08},
    {"amplitude": 10.0, "r": 4.4, "sigma": 0.1},
)


def synthetic(qmax, *, window="none", extra=(), noise=0.0, slope=0.0, seed=1):
    """A D-space function of Gaussian peaks through RMCProfile's transform: Σ ΔQ w_k S(Q_k) sin(Q_k r) · 2/π.

    The Q grid is ΔQ, 2ΔQ, … Q_max (RMCProfile's rectangle rule), so the only
    cut is at Q_max. Returns (y, the actual termination artefact Δ).
    """
    q = R.q_grid(DQ, DQ, int(round(qmax / DQ)))
    weights = R.transform_weights(q, DQ, "rectangle", window)
    peaks = list(PEAKS) + list(extra)
    spectrum = sum(R.gaussian_sine(q, peak) for peak in peaks)
    clean = R.apply_transform(q, weights, spectrum, R_GRID)
    truth = clean - sum(R.gaussian_d(R_GRID, peak) for peak in peaks)
    y = clean - slope * R_GRID
    if noise:
        y = y + np.random.default_rng(seed).normal(0.0, noise, R_GRID.size)
    return y, truth


def analyse(qmax, **kwargs):
    y, truth = synthetic(qmax, **kwargs)
    analysis = R.analyse_series(R_GRID, y, space="D", q_min=DQ, q_max=qmax, dq=DQ, rule="rectangle",
                                window=kwargs.get("window", "none"), r_first=2.3)
    return analysis, truth


class RmcProfileModelGr(unittest.TestCase):
    """What RMCProfile compares a real-space fit with, from its own output files."""

    @classmethod
    def setUpClass(cls):
        cls.evidence = json.loads(EVIDENCE.read_text())
        cls.r = np.array(cls.evidence["r"])
        partials = cls.evidence["partials"]
        # A neutron G(r) with NEUTRON_COEFFICIENTS 1.0: Σ w_ij (g_ij − 1) = g − 1.
        cls.histogram = np.interp(cls.r, np.array(partials["r"]), np.array(partials["g"]) - 1.0)
        cls.rp = np.array(partials["r"])
        cls.gp = np.array(partials["g"]) - 1.0

    def sinc_d_space(self, qmax):
        """(1/r)(1/π) Σ_j Δr r_j G(r_j) [sin(Q(r − r_j))/(r − r_j) − sin(Q(r + r_j))/(r + r_j)]."""
        a = self.r[:, None] - self.rp[None, :]
        b = self.r[:, None] + self.rp[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            first = np.where(a == 0, qmax, np.sin(qmax * a) / a)
        kernel = (first - np.sin(qmax * b) / b) / math.pi
        step = self.rp[1] - self.rp[0]
        return (kernel @ (self.rp * self.gp)) * step / self.r

    def test_without_convolve_the_model_is_the_histogram(self):
        calc = np.array(self.evidence["neutron"]["2"]["calc"])
        self.assertLess(np.max(np.abs(calc - self.histogram)), 1e-12)

    def test_with_convolve_the_model_is_the_d_space_sinc_above_r_cutoff(self):
        sinc = self.sinc_d_space(20.0)
        for key, cutoff in (("1", 2.0), ("3", 0.5)):
            calc = np.array(self.evidence["neutron"][key]["calc"])
            above = self.r >= cutoff - 1e-9
            self.assertLess(np.max(np.abs(calc[above] - sinc[above])), 1e-12, key)
            # Below R_CUTOFF (strictly): the baseline, −Σ w_ij = −1 here.
            self.assertTrue(np.all(calc[~above] == -1.0), key)

    def test_the_sinc_in_g_space_is_not_it(self):
        # Convolving G(r) itself (no r factor) misses by percent: the D-space
        # form above is RMCProfile's.
        a = self.r[:, None] - self.rp[None, :]
        b = self.r[:, None] + self.rp[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            first = np.where(a == 0, 20.0, np.sin(20.0 * a) / a)
        g_space = ((first - np.sin(20.0 * b) / b) / math.pi) @ self.gp * (self.rp[1] - self.rp[0])
        calc = np.array(self.evidence["neutron"]["1"]["calc"])
        self.assertGreater(np.max(np.abs(calc - g_space)[self.r > 2.0]), 1e-2)

    def test_xray_real_space_is_the_rectangle_transform_on_its_grid(self):
        xray = self.evidence["xray"]
        q = np.arange(1, int(round(xray["qMax"] / xray["qStep"])) + 1) * xray["qStep"]
        step = self.rp[1] - self.rp[0]
        qf = np.sin(np.outer(q, self.rp)) @ (step * self.rp * self.gp)
        model = (2.0 / math.pi) * (np.sin(np.outer(self.r, q)) @ (xray["qStep"] * qf)) / self.r
        calc = np.array(xray["calc"])
        # The file has 7 decimals; its largest |value| is 4.8.
        self.assertLess(np.max(np.abs(calc - model)) / np.max(np.abs(calc)), 2e-5)

    def test_the_demo_model_g_is_the_same_transform_of_the_model_f(self):
        ft = read_rmc_csv(DEMO / "GTS_250K_FT_XFQ1.csv")
        fq = read_rmc_csv(DEMO / "GTS_250K_XFQ1.csv")
        model = rmc_ft(fq.data[0], fq.data[1], ft.data[0])
        # 1.6e-5 of the largest |G| (the experiment reproduces to 1.5e-14,
        # test_ft_uncertainty); the termination ripples are percents.
        self.assertLess(np.max(np.abs(model - ft.data[1])) / np.max(np.abs(ft.data[1])), 3e-5)

    def test_fit_space_risks(self):
        self.assertEqual(R.fit_space_risk("xpdf", {}, None)["risk"], "matched")
        self.assertEqual(R.fit_space_risk("reciprocal", {}, None)["risk"], "none")
        self.assertEqual(R.fit_space_risk("npdf", {"convolve": False}, None)["risk"], "high")
        convolved = {"convolve": True, "convolveQ": 28.0}
        self.assertEqual(R.fit_space_risk("npdf", convolved, None)["risk"], "unknown")
        known = {"source": "stog", "qMax": 28.0, "window": "none"}
        self.assertEqual(R.fit_space_risk("npdf", convolved, known)["risk"], "matched")
        self.assertEqual(R.fit_space_risk("npdf", convolved, {**known, "window": "lorch"})["risk"], "high")
        self.assertEqual(R.fit_space_risk("npdf", convolved, {**known, "qMax": 32.0})["risk"], "high")
        self.assertEqual(R.fit_space_risk("npdf", {"convolve": True}, known)["because"], "sinc-no-q")
        self.assertEqual(R.fit_space_risk("xray-real", {"qMax": 20.0}, None)["risk"], "unknown")
        self.assertEqual(R.fit_space_risk("xray-real", {"qMax": 20.0}, {**known, "qMax": 20.0})["risk"], "matched")


class Resolution(unittest.TestCase):
    def test_rectangle_and_lorch_kernels(self):
        q = R.q_grid(DQ, DQ, 2500)  # to 25 Å⁻¹
        rect = R.resolution(q, R.transform_weights(q, DQ, "rectangle", "none"))
        self.assertAlmostEqual(rect["nyquistStep"], math.pi / 25.0)
        self.assertAlmostEqual(rect["ripplePeriod"], 2 * math.pi / 25.0)
        # sin(x)/x = 1/2 at x = 1.8955: FWHM 3.791/Q_max.
        self.assertAlmostEqual(rect["fwhm"] * 25.0, 3.7909, delta=0.003)
        self.assertAlmostEqual(rect["mainLobe"], math.pi / 25.0, delta=2e-4)
        lorch = R.resolution(q, R.transform_weights(q, DQ, "rectangle", "lorch"))
        self.assertGreater(lorch["fwhm"], 1.4 * rect["fwhm"])
        self.assertGreater(lorch["mainLobe"], 1.4 * rect["mainLobe"])

    def test_ripple_pattern_is_the_transform_less_the_uncut_one(self):
        # Δ = T[d] − T_[Q_min,∞)[d], against the brute-force sum to 300 Å⁻¹.
        q = R.q_grid(0.5, DQ, 1951)
        peak = PEAKS[0]
        for rule in ("rectangle", "trapezoid"):
            weights = R.transform_weights(q, DQ, rule, "none")
            wide = R.q_grid(0.5, DQ, int(round((300.0 - 0.5) / DQ)) + 1)
            wide_weights = np.full(wide.size, DQ)
            if rule == "trapezoid":
                wide_weights[0] *= 0.5
            brute = (R.apply_transform(q, weights, R.gaussian_sine(q, peak), R_GRID)
                     - R.apply_transform(wide, wide_weights, R.gaussian_sine(wide, peak), R_GRID))
            short = R.ripple_pattern(q, weights, DQ, rule, peak, R_GRID)
            self.assertLess(np.max(np.abs(short - brute)), 1e-12, rule)


class KnownAnswers(unittest.TestCase):
    def test_predicted_ripples_are_the_actual_ones(self):
        analysis, truth = analyse(20.0)
        self.assertEqual([round(peak["r"], 2) for peak in analysis["peaks"]], [2.5, 3.6, 4.4])
        for peak, true_peak in zip(analysis["peaks"], PEAKS):
            self.assertAlmostEqual(peak["sigma"], true_peak["sigma"], delta=0.002)
            self.assertAlmostEqual(peak["amplitude"], true_peak["amplitude"], delta=0.02 * true_peak["amplitude"])
            single = R.ripple_pattern(R.q_grid(DQ, DQ, 2000), R.transform_weights(R.q_grid(DQ, DQ, 2000), DQ, "rectangle", "none"),
                                      DQ, "rectangle", true_peak, R_GRID)
            self.assertEqual(len(peak["ripples"]), 4)
            for ripple in peak["ripples"]:
                index = int(round(ripple["r"] / 0.01)) - 1
                extrema = [i for i in R.local_extrema(single) if abs(i - index) <= 1]
                self.assertTrue(extrema, ripple)
                self.assertAlmostEqual(ripple["amplitude"], single[extrema[0]], delta=0.02 * abs(single[extrema[0]]))
        # The summed prediction is the actual artefact beyond the peaks' cores.
        predicted = analysis["_arrays"]["predicted"]
        cores = np.any([np.abs(R_GRID - p["r"]) < 2.5 * p["sigma"] for p in PEAKS], axis=0)
        self.assertLess(np.max(np.abs(predicted - truth)[~cores]), 0.02 * max(p["amplitude"] for p in PEAKS))

    def test_a_ripple_is_ripple_consistent(self):
        analysis, _ = analyse(20.0)
        first = analysis["peaks"][0]["ripples"]
        for ripple in first:
            reading = R.classify_feature(analysis, ripple["r"])
            self.assertEqual(reading["class"], "ripple-consistent", ripple)

    def test_a_planted_shoulder_above_the_resolution_is_not_a_ripple(self):
        # 15 % of the first peak, 0.42 Å out (the FWHM is 0.19 Å): a peak of its own.
        shoulder = {"amplitude": 1.2, "r": 2.92, "sigma": 0.06}
        analysis, _ = analyse(20.0, extra=[shoulder])
        reading = R.classify_feature(analysis, 2.92)
        self.assertIn(reading["class"], ("resolved-peak", "larger-than-ripple"))
        # 8 % of the 4.4 Å peak, 0.46 Å out, where its ripples are 3 %.
        shoulder = {"amplitude": 0.8, "r": 4.86, "sigma": 0.06}
        analysis, _ = analyse(20.0, extra=[shoulder])
        reading = R.classify_feature(analysis, 4.86)
        self.assertIn(reading["class"], ("larger-than-ripple", "no-ripple-here"))
        # Three times the same-sign ripple it sits on: larger than a ripple.
        bump = {"amplitude": 1.5, "r": 2.11, "sigma": 0.04}
        analysis, _ = analyse(20.0, extra=[bump])
        self.assertEqual(R.classify_feature(analysis, 2.11)["class"], "larger-than-ripple")

    def test_a_ripple_filled_in_is_smaller_than_the_ripple(self):
        # A dip of the other sign on the first peak's +5 % ripple at 2.11 Å:
        # the predicted ripple is not there in full.
        dip = {"amplitude": -0.8, "r": 2.11, "sigma": 0.04}
        analysis, _ = analyse(20.0, extra=[dip])
        self.assertEqual(R.classify_feature(analysis, 2.11)["class"], "smaller-than-ripple")

    def test_two_qmax_separate_a_ripple_from_a_feature(self):
        small = {"amplitude": 0.6, "r": 5.3, "sigma": 0.07}
        low, _ = analyse(20.0, extra=[small])
        high, _ = analyse(30.0, extra=[small])
        ripple_at = [ripple["r"] for ripple in low["peaks"][0]["ripples"] if ripple["r"] > 2.5][0]
        readings = R.classify_features([low, high], [ripple_at, 5.3])
        self.assertEqual(readings[0]["datasets"][0]["class"], "ripple-consistent")
        self.assertEqual(readings[0]["across"]["verdict"], "moves-with-qmax")
        self.assertTrue(readings[0]["across"]["discriminating"])
        self.assertEqual(readings[1]["across"]["verdict"], "in-every-dataset")
        self.assertEqual([d["class"] for d in readings[1]["datasets"]], ["no-ripple-here", "no-ripple-here"])

    def test_lorch_suppresses_the_ripples_as_predicted(self):
        rect, _ = analyse(20.0)
        lorch, truth = analyse(20.0, window="lorch")
        for peak_rect, peak_lorch, true_peak in zip(rect["peaks"], lorch["peaks"], PEAKS):
            self.assertAlmostEqual(peak_lorch["sigma"], true_peak["sigma"], delta=0.002)
            largest_rect = max(abs(r["amplitude"]) for r in peak_rect["ripples"])
            # Beyond Lorch's own broadening (2σ and its kernel's main lobe): the sidelobes.
            beyond = [r for r in peak_lorch["ripples"] if abs(r["r"] - peak_lorch["r"]) > lorch["resolution"]["fwhm"]]
            self.assertTrue(beyond)
            self.assertLess(max(abs(r["amplitude"]) for r in beyond), 0.15 * largest_rect)
            q = R.q_grid(DQ, DQ, 2000)
            single = R.ripple_pattern(q, R.transform_weights(q, DQ, "rectangle", "lorch"), DQ, "rectangle", true_peak, R_GRID)
            for ripple in peak_lorch["ripples"]:
                actual = single[int(round(ripple["r"] / 0.01)) - 1]
                self.assertAlmostEqual(ripple["amplitude"], actual, delta=0.03 * abs(actual) + 1e-3)

    def test_the_low_r_level_is_the_planted_noise_and_ripple(self):
        noise = 0.02
        analysis, truth = analyse(20.0, noise=noise, slope=0.4)
        low = analysis["lowR"]
        self.assertEqual(low["state"], "measured")
        lo, hi = low["range"]
        window = (R_GRID >= lo - 1e-9) & (R_GRID <= hi + 1e-9)
        planted = truth[window] + np.random.default_rng(1).normal(0.0, noise, R_GRID.size)[window]
        rr = R_GRID[window]
        about = planted - rr * (np.dot(rr, planted) / np.dot(rr, rr))
        self.assertAlmostEqual(low["rms"], float(np.sqrt(np.mean(about * about))), delta=0.05 * low["rms"])
        self.assertAlmostEqual(low["predictedRms"], float(np.sqrt(np.mean(truth[window] ** 2))), delta=0.05 * low["predictedRms"])
        self.assertAlmostEqual(low["residualRms"], noise, delta=0.25 * noise)
        self.assertAlmostEqual(low["baseline"]["value"], -0.4, delta=0.02)

    def test_an_enforced_or_filtered_low_r_is_detected(self):
        y, truth = synthetic(20.0, noise=0.0, slope=0.4)
        enforced = y.copy()
        enforced[R_GRID <= 1.0] = -0.4 * R_GRID[R_GRID <= 1.0]
        analysis = R.analyse_series(R_GRID, enforced, space="D", q_min=DQ, q_max=20.0, dq=DQ, rule="rectangle",
                                    window="none", r_first=2.3)
        self.assertAlmostEqual(analysis["lowR"]["enforcedBelow"], 1.0, places=6)
        self.assertEqual(analysis["lowR"]["range"][0], 1.01)
        filtered = y.copy()
        below = R_GRID <= 2.11
        filtered[below] = -0.4 * R_GRID[below] + np.random.default_rng(3).normal(0.0, 0.005, int(below.sum()))
        analysis = R.analyse_series(R_GRID, filtered, space="D", q_min=DQ, q_max=20.0, dq=DQ, rule="rectangle",
                                    window="none", r_first=2.3)
        self.assertEqual(analysis["lowR"]["state"], "filtered")

    def test_a_g_space_function_reads_alike(self):
        # The same peaks as a G(r) (y = d / r): the peaks and ripples in G's units.
        y, _ = synthetic(20.0)
        g = y / R_GRID
        analysis = R.analyse_series(R_GRID, g, space="G", q_min=DQ, q_max=20.0, dq=DQ, rule="rectangle",
                                    window="none", r_first=2.3)
        first = analysis["peaks"][0]
        self.assertAlmostEqual(first["sigma"], 0.06, delta=0.002)
        self.assertAlmostEqual(first["height"], 8.0 / first["r"], delta=0.02 * 8.0 / 2.5)
        self.assertEqual(R.classify_feature(analysis, first["ripples"][1]["r"])["class"], "ripple-consistent")

    def test_a_split_below_the_resolution(self):
        low, _ = analyse(20.0)
        reading = R.split_reading([low], 3.55, 3.68)
        self.assertTrue(reading["datasets"][0]["belowResolution"])
        self.assertFalse(R.split_reading([low], 3.4, 3.7)["datasets"][0]["belowResolution"])


class RunFiles(unittest.TestCase):
    def test_settings_and_data_processing(self):
        text = (
            "MINIMUM_DISTANCES :: 2.2 2.0 Angstrom\n"
            "NEUTRON_REAL_SPACE_DATA :: 1\n  > FILENAME :: a.gr\n  > FIT_TYPE :: G(r)\n  > CONVOLVE :: 28.5\n"
            "  > R_CUTOFF :: 1.8\n  > NEUTRON_COEFFICIENTS :: 0.07 -0.08 0.02\n"
            "XRAY_REAL_SPACE_DATA ::\n  > FILENAME :: b.gr\n  > QMIN :: 0.0\n  > QMAX :: 20\n  > QSTEP :: 0.01\n"
            "XRAY_RECIPROCAL_SPACE_DATA ::\n  > FILENAME :: c.fq\n  > CONVOLVE ::\n  > REAL_SPACE_FIT :: 1 5000 1\nEND ::\n"
        )
        settings = R.ripple_settings(text)
        self.assertEqual(settings["minimumDistances"], [2.2, 2.0])
        neutron, xray, reciprocal = settings["datasets"]
        self.assertEqual((neutron["convolveQ"], neutron["rCutoff"], neutron["coefficients"]), (28.5, 1.8, [0.07, -0.08, 0.02]))
        self.assertEqual((xray["qMin"], xray["qMax"], xray["qStep"]), (0.0, 20.0, 0.01))
        self.assertTrue(reciprocal["convolve"] and reciprocal["convolveQ"] is None and reciprocal["realSpaceFit"])
        self.assertTrue(R.is_total(neutron))
        self.assertFalse(R.is_total({"coefficients": [0.0, -0.08, 0.0]}))
        self.assertAlmostEqual(R.expected_baseline(neutron, "G"), -0.01)
        provenance = {"tool": "rmc-autoscale", "outputs": {"rmc_gr": "/x/autoscale/scale_ft_rmc.gr"},
                      "enforcement": {"cutoff": 1.5}, "provenance": {"n_q_points": 2551,
                      "config": {"qmin": 0.5, "qmax": 26.0, "lorch": True, "r_cutoff": 1.0}}}
        found = R.provenance_transform(provenance, "scale_ft_rmc.gr")
        self.assertEqual((found["qMax"], found["window"], found["filterCutoff"], found["enforceCutoff"]), (26.0, "lorch", 1.0, 1.5))
        self.assertAlmostEqual(found["dQ"], 0.01)
        self.assertIsNone(R.provenance_transform(provenance, "other.gr"))
        self.assertTrue(R.auto_stog_title("        5000\nrmc-autoscale 0.4.0: a=1 b=0\n 0.01 -1\n"))
        self.assertFalse(R.auto_stog_title(" 5000\nrmc G(r)\n"))

    def test_a_neutron_g_with_its_stog_input(self):
        # A run fitting a G(r) made by stog (Lorch on, Q_max 25): read through
        # the stog input beside it; the model is not truncated (no CONVOLVE).
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            q = R.q_grid(DQ, DQ, 2500)
            weights = R.transform_weights(q, DQ, "trapezoid", "lorch")
            d = R.apply_transform(q, weights, sum(R.gaussian_sine(q, p) for p in PEAKS), R_GRID)
            g = d / R_GRID - 0.05
            (folder / "run.dat").write_text(
                "MINIMUM_DISTANCES :: 2.2 Angstrom\nNEUTRON_REAL_SPACE_DATA :: 1\n  > FILENAME :: scale_ft_rmc.gr\n"
                "  > FIT_TYPE :: G(r)\n  > NEUTRON_COEFFICIENTS :: 0.05\nEND ::\n")
            (folder / "scale_ft_rmc.gr").write_text(
                f"{R_GRID.size}\nstog G(r)\n" + "".join(f"{x:.6f} {v:.12f}\n" for x, v in zip(R_GRID, g)))
            inp = ["1", "data.sq", "0.01 25.0", "0 1", "0", "scale.fq", "scale.gr", "12", "1200", "Y", "0.09", "0",
                   "N", "Y", "1.0", "scale_ft.sq", "scale_ft.gr", "1.0", "scale_ft_rmc.fq", "scale_ft_rmc.gr",
                   "scale_ft_rmc.dr", "1.0 0 0"]
            (folder / "stog.inp").write_text("\n".join(inp) + "\n")
            keep = R_GRID >= 1.5
            (folder / "run_PDF1.csv").write_text("r (Ang), PDF (RMC), PDF (Expt)\n" + "".join(
                f"{float(x)!r} , {float(v)!r} , {float(v)!r}\n" for x, v in zip(R_GRID[keep], g[keep])))
            out = R.ripples_for_run([folder / "run_PDF1.csv"], features=[2.26])
        dataset = out["datasets"][0]
        self.assertEqual(dataset["status"], "assessed")
        self.assertEqual(dataset["fit"]["risk"], "high")
        self.assertEqual(dataset["transform"]["source"], "stog")
        self.assertAlmostEqual(dataset["transform"]["qMax"], 25.0)
        self.assertEqual(dataset["transform"]["window"], "lorch")
        self.assertEqual(dataset["range"][0], 0.01)  # the input file reaches below the fitted range
        self.assertEqual(dataset["processing"]["filterCutoff"], 1.0)  # the input's Fourier filter, Y at 1.0 Å
        self.assertAlmostEqual(dataset["lowR"]["baseline"]["value"], -0.05, delta=0.005)
        self.assertAlmostEqual(dataset["lowR"]["baseline"]["expected"], -0.05)
        self.assertAlmostEqual(dataset["peaks"][0]["sigma"], 0.06, delta=0.003)
        self.assertEqual(len(out["features"]), 1)


class DemoRun(unittest.TestCase):
    def test_the_demo_end_to_end(self):
        paths = sorted(p for p in DEMO.glob("GTS_250K_*.csv") if R.dataset_file(p.name))
        out = R.ripples_for_run(paths, r_range=[2.0, 6.0], split=[2.45, 2.6])
        by_kind = {dataset["kind"]: dataset for dataset in out["datasets"]}
        xpdf = by_kind["xpdf"]
        self.assertEqual(xpdf["status"], "assessed")
        self.assertEqual(xpdf["fit"]["risk"], "matched")
        self.assertEqual((xpdf["transform"]["qMin"], xpdf["transform"]["qMax"], xpdf["transform"]["nQ"]), (0.8, 26.98, 2619))
        self.assertLess(xpdf["transform"]["reproduction"]["residual"], 1e-12)
        self.assertAlmostEqual(xpdf["resolution"]["fwhm"], 3.79 / 26.98, delta=0.003)
        self.assertEqual(xpdf["rFirst"], {"value": 2.2, "source": "partials"})
        self.assertAlmostEqual(xpdf["peaks"][0]["r"], 2.51, delta=0.02)
        self.assertEqual(xpdf["lowR"]["state"], "measured")
        # The model is truncated alike: its own low-r ripples are those the
        # strong peaks predict, within the Gaussian approximation.
        self.assertLess(abs(xpdf["lowR"]["modelRms"] / xpdf["lowR"]["predictedRms"] - 1), 0.5)
        self.assertEqual(by_kind["xray-reciprocal"]["fit"]["risk"], "none")
        self.assertTrue(out["found"])
        self.assertTrue(1 <= len(out["features"]) <= R.MAX_FEATURES)
        self.assertTrue(out["split"]["datasets"][0]["belowResolution"] is False)


if __name__ == "__main__":
    unittest.main()
