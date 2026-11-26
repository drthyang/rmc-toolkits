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
2. The r-resolution of the run's own kernel (rectangle and Lorch).
3. The run's files: the .dat settings, an Auto StoG provenance, a stog input
   beside a neutron G(r); the Demo end to end.

The tests of the strong-peak fit, the predicted ripples, the low-r level and
the feature readings are on the branch wip/joint-fit with that code.
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


def transformed_peaks(q, weights):
    """The peaks' D-space function through the transform: (2/π) Σ_k w_k S(Q_k) sin(Q_k r), S their exact sine transform."""
    spectrum = sum(p["amplitude"] * math.sqrt(2.0 * math.pi) * p["sigma"] * np.exp(-0.5 * q * q * p["sigma"] ** 2)
                   * np.sin(q * p["r"]) for p in PEAKS)
    return np.sin(np.outer(R_GRID, q)) @ ((2.0 / math.pi) * weights * spectrum)


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
            d = transformed_peaks(q, weights)
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
            out = R.ripples_for_run([folder / "run_PDF1.csv"])
        dataset = out["datasets"][0]
        self.assertEqual(dataset["status"], "assessed")
        self.assertEqual(dataset["fit"]["risk"], "high")
        self.assertEqual(dataset["transform"]["source"], "stog")
        self.assertAlmostEqual(dataset["transform"]["qMax"], 25.0)
        self.assertEqual(dataset["transform"]["window"], "lorch")
        self.assertIn("input", dataset["transform"])  # the input file maps onto the fitted one
        self.assertEqual(dataset["processing"]["filterCutoff"], 1.0)  # the input's Fourier filter, Y at 1.0 Å


class DemoRun(unittest.TestCase):
    def test_the_demo_end_to_end(self):
        paths = sorted(p for p in DEMO.glob("GTS_250K_*.csv") if R.dataset_file(p.name))
        out = R.ripples_for_run(paths)
        by_kind = {dataset["kind"]: dataset for dataset in out["datasets"]}
        xpdf = by_kind["xpdf"]
        self.assertEqual(xpdf["status"], "assessed")
        self.assertEqual(xpdf["fit"]["risk"], "matched")
        self.assertEqual((xpdf["transform"]["qMin"], xpdf["transform"]["qMax"], xpdf["transform"]["nQ"]), (0.8, 26.98, 2619))
        self.assertLess(xpdf["transform"]["reproduction"]["residual"], 1e-12)
        self.assertAlmostEqual(xpdf["resolution"]["fwhm"], 3.79 / 26.98, delta=0.003)
        self.assertEqual(by_kind["xray-reciprocal"]["fit"]["risk"], "none")


if __name__ == "__main__":
    unittest.main()
