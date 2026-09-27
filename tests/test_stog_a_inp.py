# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Classic stog.inp first-peak line (line 22) semantics in the Auto StoG front ends.

Regression tests for the 0.6.0 audit (stog-a group): the r0 fallback took
max(peak_cutoff, peak_rmin), which for a line whose first-peak window starts below the
cleanup cutoff (the very case the window exists for) put the C2 window over the first
peak -- '2.3 1.6 2.2' on a synthetic with its first shell at 1.7 A gave a = -4.9.
"""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import StogInput, write_stog_xy
from rmc_toolkits.scaling_cli import main
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

A_TRUE, B_TRUE = 10.0, -9.0
RHO0 = 0.05

INP = (
    "1\nsample.sq\n0.6 29.4\n-9 0.1\n0\nscale.fq\nscale.gr\n25\n2500\nN\n"
    "0.05\n0\nN\nY\n1.0\nscale_ft.sq\nscale_ft.gr\n0.02\n"
    "scale_ft_rmc.fq\nscale_ft_rmc.gr\nscale_ft_rmc.dr\n{peak_line}\n"
)


def short_onset_sq():
    """The repo synthetic moved down: first-shell onset 1.70 A (peak 1.85 A)."""
    q = np.arange(20, 981) * 0.03
    r = np.arange(1, 12001) * 0.005
    g = 0.5 * (1.0 + np.tanh((r - 1.70) / 0.07)) + 1.6 * np.exp(-0.5 * ((r - 1.85) / 0.15) ** 2)
    sq_true = fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, RHO0), q))
    return q, (sq_true - B_TRUE) / A_TRUE


def write_run(directory, peak_line):
    q, sq = short_onset_sq()
    write_stog_xy(Path(directory) / "sample.sq", q, sq, title="synthetic onset 1.70")
    (Path(directory) / "stog.inp").write_text(INP.format(peak_line=peak_line))
    return Path(directory) / "stog.inp"


def inp_with(peak_cutoff, peak_rmin, peak_rmax):
    return StogInput(
        n_files=1, data_file="x", qmin=0.6, qmax=30, yoffset=0, yscale=1, qoffset=0,
        out_sq="a", out_gr="b", rmax=25, nr=2500, lorch=False, rho0=0.05, yoffset2=0,
        try_again=False, use_filter=True, r_cutoff=1.0, out_ft_sq="c", out_ft_gr="d",
        b_avg_sq=0.02, out_rmc_fq="e", out_rmc_gr="f", out_rmc_dr="g",
        peak_cutoff=peak_cutoff, peak_rmin=peak_rmin, peak_rmax=peak_rmax,
    )


class StogInpClosestApproachTests(unittest.TestCase):
    def test_cli_uses_the_first_peak_window_start(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp = write_run(tmp, "2.3 1.6 2.2")
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = main([str(inp), "--out-dir", str(Path(tmp) / "out")])
            self.assertEqual(code, 0, err.getvalue())
            provenance = json.loads((Path(tmp) / "out" / "stog_provenance.json").read_text())
        # Pre-fix: r0 = max(2.3, 1.6) = 2.3, window [1.2, 2.05] across the
        # first peak, a = -4.92.
        self.assertAlmostEqual(provenance["provenance"]["config"]["r0"], 1.6)
        a = provenance["diagnostics"]["a"]
        self.assertLess(abs(a / A_TRUE - 1.0), 0.03, a)

    def test_rule(self):
        from rmc_toolkits.scaling_cli import stog_inp_closest_approach

        # A first-peak window starting inside the cleanup radius: its start.
        self.assertEqual(stog_inp_closest_approach(inp_with(2.3, 1.6, 2.2), 1.0), 1.6)
        # Window outside [0, cutoff] (the Mn3Sn 59438 line): the cutoff.
        self.assertEqual(stog_inp_closest_approach(inp_with(2.48, 2.65, 3.1), 1.0), 2.48)
        # An empty/inverted window keeps nothing: the cutoff.
        self.assertEqual(stog_inp_closest_approach(inp_with(2.7, 2.3, 2.2), 1.0), 2.7)
        # '1.0 0 0' (FeCoSn): too low for a fit window -> detect from the data.
        self.assertIsNone(stog_inp_closest_approach(inp_with(1.0, 0.0, 0.0), 1.0))
        # A sliver is no window: the proxy must leave >= MIN_AUTO_WINDOW (0.1 A)
        # above r_cutoff + 0.2, else r0 is detected. '1.46 0 0' at r_cutoff 1.0
        # left [1.2, 1.21] (FeCoSn a 15 % low); '1.0 0 0' at r_cutoff 0.5 left
        # [0.7, 0.75] (a 43 % low), both reported converged.
        self.assertIsNone(stog_inp_closest_approach(inp_with(1.46, 0.0, 0.0), 1.0))
        self.assertIsNone(stog_inp_closest_approach(inp_with(1.5, 0.0, 0.0), 1.0))
        self.assertIsNone(stog_inp_closest_approach(inp_with(1.0, 0.0, 0.0), 0.5))
        self.assertIsNone(stog_inp_closest_approach(inp_with(1.0, 0.0, 0.0), 0.54))
        self.assertEqual(stog_inp_closest_approach(inp_with(1.55, 0.0, 0.0), 1.0), 1.55)
        self.assertEqual(stog_inp_closest_approach(inp_with(1.0, 0.0, 0.0), 0.4), 1.0)


XRAY_RUN = Path(__file__).resolve().parents[1] / "data" / "stog_tests" / "199K"


@unittest.skipUnless((XRAY_RUN / "stog_input.dat").exists(), "FeCoSn 199K run not present")
class StogInpSliverWindowRealDataTests(unittest.TestCase):
    """FeCoSn 199 K: a line-22 cutoff just above r_cutoff + 0.45 no longer pins a sliver."""

    def run_cli(self, tmp, peak_line, extra=()):
        lines = (XRAY_RUN / "stog_input.dat").read_text().splitlines()
        lines[1] = str(XRAY_RUN / lines[1].strip())
        lines[21] = peak_line
        inp = Path(tmp) / "stog.inp"
        inp.write_text("\n".join(lines) + "\n")
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            code = main([str(inp), "--out-dir", str(Path(tmp) / "out"), *extra])
        self.assertEqual(code, 0, err.getvalue())
        return json.loads((Path(tmp) / "out" / "stog_provenance.json").read_text())

    def test_sliver_line_falls_back_to_detection(self):
        with tempfile.TemporaryDirectory() as tmp:
            reference = self.run_cli(tmp, "1.0 0 0")["diagnostics"]["a"]
        for peak_line, extra in (("1.46 0 0", ()), ("1.0 0 0", ("--r-cutoff", "0.54"))):
            with self.subTest(peak_line=peak_line, extra=extra):
                with tempfile.TemporaryDirectory() as tmp:
                    payload = self.run_cli(tmp, peak_line, extra)
                # r0 came from the data (placement), not from line 22.
                self.assertTrue(payload["provenance"].get("window_refined"))
                self.assertNotEqual(payload["provenance"]["r0_detected"], float(peak_line.split()[0]))
                lo, hi = payload["diagnostics"]["r_fit_window"]
                self.assertGreaterEqual(hi - lo, 0.1)
                if not extra:
                    # Pre-fix: [1.2, 1.21], a = 1.0013 (15 % low).
                    self.assertLess(abs(payload["diagnostics"]["a"] / reference - 1), 0.01)


if __name__ == "__main__":
    unittest.main()
