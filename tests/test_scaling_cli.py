# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from rmc_toolkits.parsers import read_stog_xy, write_stog_xy
from rmc_toolkits.scaling_cli import main
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

ROOT = Path(__file__).resolve().parents[1]
STOG_RUN = ROOT / "data" / "stog_tests" / "stog_59438"
# stog_59438's parameters are recovered from its outputs (see tests/test_scaling.py),
# so the CLI parity runs in --data mode from the rebinned file -- no stog.inp needed.
STOG_DATA = "PG3_59438_SQ_rebin.dat"


def _first_existing_run(*names: str) -> Path:
    """Same FeCoSn-series candidate logic as tests/test_scaling.py (kept in sync)."""
    candidates = [ROOT / "data" / "stog_tests" / name for name in names]
    for candidate in candidates:
        if (candidate / "stog_input.dat").exists():
            return candidate
    return candidates[0]


XRAY_RUN = _first_existing_run("100K", "199K")

requires_stog_run = unittest.skipUnless(
    (STOG_RUN / STOG_DATA).exists(), "stog_59438 example run not present"
)
requires_xray_run = unittest.skipUnless(
    (XRAY_RUN / "stog_input.dat").exists(), "FeCoSn x-ray run not present"
)

RHO0 = 0.05
B2 = 0.02
A_TRUE, B_TRUE = 10.0, -9.0

# Classic 23-line stog.inp exercising the validated single-dataset layout;
# the -9 0.1 scale line encodes the "hand" correction a=10, b=-9.
INP_TEMPLATE = (
    "1\n{data}\n0.60 30.0\n-9 0.1\n0\nscale.fq\nscale.gr\n{rmax}\n{nr}\nN\n"
    "0.05\n0\nN\nY\n1.0\nscale_ft.sq\nscale_ft.gr\n0.02\n"
    "scale_ft_rmc.fq\nscale_ft_rmc.gr\nscale_ft_rmc.dr\n2.48 2.65 3.1\n"
)


def synthetic_g(r: np.ndarray) -> np.ndarray:
    """Same synthetic model as tests/test_scaling.py (kept in sync)."""
    onset = 0.5 * (1.0 + np.tanh((r - 2.65) / 0.07))
    peak = 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    return onset + peak


def run_cli(args):
    """Invoke the CLI in-process, capturing stdout/stderr."""
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main([str(arg) for arg in args])
    return code, out.getvalue(), err.getvalue()


class CliSyntheticBase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.q = np.arange(20, 981) * 0.03  # 0.60 .. 29.40
        r = np.arange(1, 12001) * 0.005
        gpdf = g_to_gpdf(r, synthetic_g(r), RHO0)
        cls.sq_true = fq_to_sq(cls.q, gpdf_to_fq(r, gpdf, cls.q))
        cls.sq_meas = (cls.sq_true - B_TRUE) / A_TRUE

    def make_run_dir(self, tmp: str) -> Path:
        run = Path(tmp)
        write_stog_xy(run / "synth.dat", self.q, self.sq_meas, title="synthetic")
        (run / "stog.inp").write_text(
            INP_TEMPLATE.format(data="synth.dat", rmax=25, nr=1000)
        )
        return run


class StogInpModeTests(CliSyntheticBase):
    def test_auto_end_to_end(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.make_run_dir(tmp)
            code, out, err = run_cli([run / "stog.inp"])
            self.assertEqual(code, 0, msg=err)
            self.assertIn("hand values: a = 10", out)

            out_dir = run / "autoscale"
            for name in (
                "scale.fq", "scale.gr", "scale_ft.sq", "scale_ft.gr",
                "scale_ft_rmc.fq", "scale_ft_rmc.gr", "scale_ft_rmc.dr",
                "ft.dat", "stog_provenance.json",
            ):
                self.assertTrue((out_dir / name).exists(), name)

            # Classic identity: sq_filtered = sq_scaled - (ft - 1).
            ft = read_stog_xy(out_dir / "ft.dat")
            filtered = read_stog_xy(out_dir / "scale_ft.sq")
            scaled_q = read_stog_xy(out_dir / "scale.fq")
            np.testing.assert_allclose(
                filtered[1], scaled_q[1] - (ft[1] - 1.0), atol=1e-10
            )

            payload = json.loads((out_dir / "stog_provenance.json").read_text())
            a, b = payload["diagnostics"]["a"], payload["diagnostics"]["b"]
            self.assertTrue(payload["diagnostics"]["converged"])
            self.assertEqual(payload["provenance"]["mode"], "auto")
            self.assertEqual(payload["stog_inp_reference"]["a"], A_TRUE)
            self.assertLess(abs(a - A_TRUE) / A_TRUE, 0.05)
            self.assertLess(abs(b - B_TRUE) / abs(B_TRUE), 0.05)

            # Written scaled S(Q) is exactly a*S_meas + b on the cropped grid.
            scaled = read_stog_xy(out_dir / "scale.fq")
            expected = a * np.interp(scaled[0], self.q, self.sq_meas) + b
            np.testing.assert_allclose(scaled[1], expected, atol=1e-10)

            # Classic enforcement (on by default in stog.inp mode): GK is the
            # flat -<b>^2 below the cutoff, and D(r) is on the theory line.
            gk = read_stog_xy(out_dir / "scale_ft_rmc.gr")
            below = gk[0] <= 2.48
            np.testing.assert_allclose(gk[1][below], -B2, atol=1e-12)
            dr = read_stog_xy(out_dir / "scale_ft_rmc.dr")
            np.testing.assert_allclose(
                dr[1][below],
                -4.0 * np.pi * RHO0 * B2 * dr[0][below],
                atol=1e-12,
            )

            # The filtered .gr carries the D(r) companion column.
            ftgr = read_stog_xy(out_dir / "scale_ft.gr")
            self.assertEqual(ftgr.shape[0], 3)
            np.testing.assert_allclose(
                ftgr[2], 4.0 * np.pi * RHO0 * ftgr[0] * ftgr[1], rtol=1e-10
            )

    def test_manual_mode_reproduces_hand_scaling(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.make_run_dir(tmp)
            code, out, err = run_cli([run / "stog.inp", "--manual"])
            self.assertEqual(code, 0, msg=err)
            self.assertIn("manual (fixed a, b)", out)
            scaled = read_stog_xy(run / "autoscale" / "scale.fq")
            expected = A_TRUE * np.interp(scaled[0], self.q, self.sq_meas) + B_TRUE
            np.testing.assert_allclose(scaled[1], expected, atol=1e-12)
            payload = json.loads(
                (run / "autoscale" / "stog_provenance.json").read_text()
            )
            self.assertEqual(payload["provenance"]["mode"], "manual")

    def test_no_clobber_guard_and_force(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.make_run_dir(tmp)
            code, _, _ = run_cli([run / "stog.inp", "--manual"])
            self.assertEqual(code, 0)
            code, _, err = run_cli([run / "stog.inp", "--manual"])
            self.assertEqual(code, 2)
            self.assertIn("refusing to overwrite", err)
            code, _, err = run_cli([run / "stog.inp", "--manual", "--force"])
            self.assertEqual(code, 0, msg=err)

    def test_no_enforce_keeps_honest_low_r(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.make_run_dir(tmp)
            code, _, err = run_cli([run / "stog.inp", "--manual", "--no-enforce"])
            self.assertEqual(code, 0, msg=err)
            gk = read_stog_xy(run / "autoscale" / "scale_ft_rmc.gr")
            below = gk[0] <= 2.0
            # Without enforcement the low-r region keeps its real ripples.
            self.assertGreater(float(np.abs(gk[1][below] + B2).max()), 1e-6)


class DataModeTests(CliSyntheticBase):
    def test_header_metadata_and_formula(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "sample.dat"
            with data.open("w", encoding="utf-8") as handle:
                handle.write(
                    "TITLE :: synthetic nickel\n"
                    "NUMBER_DENSITY :: 0.05 Angstrom^(-3)\n"
                    "MINIMUM_DISTANCES :: 2.65 3.2 Angstrom\n"
                )
                for qi, si in zip(self.q, self.sq_meas):
                    handle.write(f"{qi:.6f} {si:.10f}\n")
            code, out, err = run_cli(
                [
                    "--data", data, "--qmin", "0.6", "--qmax", "30",
                    "--formula", "Ni", "--rmax", "25", "--nr", "1000",
                ]
            )
            self.assertEqual(code, 0, msg=err)
            out_dir = Path(tmp) / "autoscale"
            self.assertTrue((out_dir / "sample.sq").exists())
            self.assertTrue((out_dir / "sample_rmc.gr").exists())
            payload = json.loads((out_dir / "sample_provenance.json").read_text())
            config = payload["provenance"]["config"]
            self.assertAlmostEqual(config["rho0"], 0.05)  # NUMBER_DENSITY header
            self.assertAlmostEqual(config["r0"], 2.65)  # MINIMUM_DISTANCES header
            self.assertAlmostEqual(config["b_avg_sq"], 1.0609, delta=0.01)  # Ni
            # <b^2> from the formula enables the Q->0 amplitude diagnostic.
            self.assertIn("amplitude_concordance", payload["diagnostics"])

    def test_sigma_column_and_fixed_scaling(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "sample.dat"
            sigma = np.full(self.q.size, 0.005)
            write_stog_xy(data, self.q, self.sq_meas, title="s", extra=sigma)
            code, _, err = run_cli(
                [
                    "--data", data, "--qmin", "0.6", "--qmax", "30",
                    "--b-avg-sq", B2, "--rho0", RHO0, "--r0", "2.65",
                    "--rmax", "25", "--nr", "1000",
                    "--scale", A_TRUE, "--offset", B_TRUE,
                ]
            )
            self.assertEqual(code, 0, msg=err)
            scaled = read_stog_xy(Path(tmp) / "autoscale" / "sample.sq")
            expected = A_TRUE * np.interp(scaled[0], self.q, self.sq_meas) + B_TRUE
            np.testing.assert_allclose(scaled[1], expected, atol=1e-12)

    def test_fz_amplitude_mode(self):
        head = self.q <= self.q[0] + 1.0
        _, s_true_0 = np.polyfit(self.q[head], self.sq_true[head], 1)
        b_sq_avg = B2 * (1.0 - float(s_true_0))
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "sample.dat"
            write_stog_xy(data, self.q, self.sq_meas, title="s")
            base = [
                "--data", data, "--qmin", "0.6", "--qmax", "30",
                "--b-avg-sq", B2, "--rho0", RHO0, "--r0", "2.65",
                "--rmax", "25", "--nr", "1000",
            ]
            code, out, err = run_cli(base + ["--b-sq-avg", b_sq_avg, "--amplitude", "fz"])
            self.assertEqual(code, 0, msg=err)
            self.assertIn("FZ-limit amplitude", out)
            payload = json.loads((Path(tmp) / "autoscale" / "sample_provenance.json").read_text())
            self.assertEqual(payload["provenance"]["config"]["amplitude_criterion"], "fz")
            self.assertLess(abs(payload["diagnostics"]["a"] - A_TRUE) / A_TRUE, 0.05)

            code, _, err = run_cli(base + ["--amplitude", "fz", "--force"])
            self.assertEqual(code, 2)
            self.assertIn("requires <b^2>", err)

            code, _, err = run_cli(
                base + ["--b-sq-avg", b_sq_avg, "--amplitude", "fz",
                        "--manual", "--scale", "10", "--force"]
            )
            self.assertEqual(code, 2)
            self.assertIn("cannot be combined", err)

    def test_estimate_rho0_flag(self):
        head = self.q <= self.q[0] + 1.0
        _, s_true_0 = np.polyfit(self.q[head], self.sq_true[head], 1)
        b_sq_avg = B2 * (1.0 - float(s_true_0))
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "sample.dat"
            write_stog_xy(data, self.q, self.sq_meas, title="s")
            base = [
                "--data", data, "--qmin", "0.6", "--qmax", "30",
                "--b-avg-sq", B2, "--r0", "2.65",
                "--rmax", "25", "--nr", "1000",
            ]
            # Seeded at the wrong density: the estimate must land on RHO0 and
            # replace the seed for the run (config + provenance record it).
            code, out, err = run_cli(
                base + ["--rho0", "0.02", "--b-sq-avg", b_sq_avg, "--estimate-rho0"]
            )
            self.assertEqual(code, 0, msg=err)
            self.assertIn("rho0 self-consistency", out)
            payload = json.loads(
                (Path(tmp) / "autoscale" / "sample_provenance.json").read_text()
            )
            estimate = payload["rho0_estimate"]
            self.assertTrue(estimate["converged"])
            self.assertLess(abs(estimate["rho0"] - RHO0) / RHO0, 0.05)
            self.assertAlmostEqual(
                payload["provenance"]["config"]["rho0"], estimate["rho0"]
            )

            code, _, err = run_cli(base + ["--rho0", "0.02", "--estimate-rho0", "--force"])
            self.assertEqual(code, 2)
            self.assertIn("requires <b^2>", err)

            code, _, err = run_cli(
                base + ["--rho0", "0.02", "--b-sq-avg", b_sq_avg, "--estimate-rho0",
                        "--manual", "--scale", "10", "--force"]
            )
            self.assertEqual(code, 2)
            self.assertIn("cannot be combined", err)

    def test_error_cases(self):
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "sample.dat"
            write_stog_xy(data, self.q, self.sq_meas)
            cases = [
                ([], "exactly one input"),
                (["--data", data, "--qmax", "30", "--b-avg-sq", "0.02",
                  "--rho0", "0.05"], "requires --qmin"),
                (["--data", data, "--qmin", "0.6", "--qmax", "30",
                  "--rho0", "0.05"], "--b-avg-sq or --formula"),
                (["--data", data, "--qmin", "0.6", "--qmax", "30",
                  "--b-avg-sq", "0.02"], "number density unknown"),
                (["--data", data, "--qmin", "0.6", "--qmax", "30",
                  "--b-avg-sq", "0.02", "--rho0", "0.05", "--manual"],
                 "requires --scale"),
                ([Path(tmp) / "missing.inp"], "not found"),
            ]
            for args, needle in cases:
                code, _, err = run_cli(args)
                self.assertEqual(code, 2, msg=f"{args} -> {err}")
                self.assertIn(needle, err)

    def test_module_entry_point(self):
        env = dict(os.environ, MPLCONFIGDIR=os.environ.get(
            "MPLCONFIGDIR", tempfile.gettempdir() + "/rmc_toolkits_matplotlib"
        ))
        result = subprocess.run(
            [sys.executable, "-m", "rmc_toolkits.scaling_cli", "--version"],
            capture_output=True, text=True, cwd=ROOT, env=env,
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn("rmc-autoscale", result.stdout)


@requires_stog_run
class ExampleRunCliTests(unittest.TestCase):
    """CLI-level parity against the complete classic Fortran stog run."""

    def test_manual_run_matches_fortran_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            # Recovered stog_59438 parameters (see tests/test_scaling.py), driven
            # through --data mode with the run's hand scaling (a = 10, b = -9).
            code, _, err = run_cli(
                [
                    "--data", STOG_RUN / STOG_DATA,
                    "--qmin", "1.0", "--qmax", "28.0",
                    "--rho0", "0.063049", "--b-avg-sq", "0.015407",
                    "--r-cutoff", "1.0", "--rmax", "50.0", "--nr", "5000",
                    "--manual", "--scale", "10.0", "--offset", "-9.0",
                    "--enforce", "--enforce-cutoff", "2.48",
                    "--peak-window", "2.65", "3.1",
                    "--out-stem", "scale", "--out-dir", tmp,
                ]
            )
            self.assertEqual(code, 0, msg=err)

            # --out-stem names the scaled S(Q) scale.sq; it holds the same scaled
            # function as the Fortran scale.fq.
            ours = read_stog_xy(Path(tmp) / "scale.sq")
            ref = read_stog_xy(STOG_RUN / "scale.fq")
            interp = np.interp(ref[0], ours[0], ours[1])
            np.testing.assert_allclose(interp, ref[1], atol=1e-8)

            ours_ft = read_stog_xy(Path(tmp) / "scale_ft.sq")
            ref_ft = read_stog_xy(STOG_RUN / "scale_ft.sq")
            interp = np.interp(ref_ft[0], ours_ft[0], ours_ft[1])
            rms = float(np.sqrt(np.mean((interp - ref_ft[1]) ** 2)))
            self.assertLess(rms, 1e-3)

            ours_gr = read_stog_xy(Path(tmp) / "scale_rmc.gr")
            ref_gr = read_stog_xy(STOG_RUN / "scale_ft_rmc.gr")
            interp = np.interp(ref_gr[0], ours_gr[0], ours_gr[1])
            rms = float(np.sqrt(np.mean((interp - ref_gr[1]) ** 2)))
            self.assertLess(rms, 2e-3)


@requires_xray_run
class XrayRunCliTests(unittest.TestCase):
    """CLI parity + auto-scale against a complete FeCoSn x-ray Fortran run."""

    def test_manual_parity_and_auto_smoke(self):
        with tempfile.TemporaryDirectory() as tmp:
            code, _, err = run_cli(
                [XRAY_RUN / "stog_input.dat", "--manual", "--out-dir", tmp]
            )
            self.assertEqual(code, 0, msg=err)

            ours = read_stog_xy(Path(tmp) / "scale.fq")
            ref = read_stog_xy(XRAY_RUN / "scale.fq")
            np.testing.assert_allclose(
                np.interp(ref[0], ours[0], ours[1]), ref[1], atol=1e-8
            )

            ours_ft = read_stog_xy(Path(tmp) / "scale_ft.sq")
            ref_ft = read_stog_xy(XRAY_RUN / "scale_ft.sq")
            rms = float(np.sqrt(np.mean(
                (np.interp(ref_ft[0], ours_ft[0], ours_ft[1]) - ref_ft[1]) ** 2
            )))
            self.assertLess(rms, 1e-3)

            # Classic ft.dat parity on the overlapping (>= Qmin) grid.
            ours_c = read_stog_xy(Path(tmp) / "ft.dat")
            ref_c = read_stog_xy(XRAY_RUN / "ft.dat")
            keep = ref_c[0] >= ours_c[0][0] - 1e-9
            rms = float(np.sqrt(np.mean(
                (np.interp(ref_c[0][keep], ours_c[0], ours_c[1]) - ref_c[1][keep]) ** 2
            )))
            self.assertLess(rms, 1e-3)

            ours_gr = read_stog_xy(Path(tmp) / "scale_ft_rmc.gr")
            ref_gr = read_stog_xy(XRAY_RUN / "scale_ft_rmc.gr")
            rms = float(np.sqrt(np.mean(
                (np.interp(ref_gr[0], ours_gr[0], ours_gr[1]) - ref_gr[1]) ** 2
            )))
            self.assertLess(rms, 2e-3)
            # Normalized x-ray enforcement: flat -1 below the 1.0 A cutoff.
            below = ours_gr[0] <= 1.0
            np.testing.assert_allclose(ours_gr[1][below], -1.0, atol=1e-12)

        with tempfile.TemporaryDirectory() as tmp:
            code, out, err = run_cli([XRAY_RUN / "stog_input.dat", "--out-dir", tmp])
            self.assertEqual(code, 0, msg=err)
            payload = json.loads(
                Path(tmp, "stog_input_provenance.json").read_text()
            )
            self.assertTrue(payload["diagnostics"]["converged"])
            # Hand values are a = 1/0.9, b = -0.111; the auto-fit should land
            # in that neighborhood on this density-limit-satisfiable dataset.
            self.assertLess(abs(payload["diagnostics"]["a"] - 1.0 / 0.9) * 0.9, 0.15)


MN3SN_500K_DATA = ROOT / "data" / "stog_tests" / "stog_500K" / "PG3_54139_SQ_rebin.dat"


@unittest.skipUnless(MN3SN_500K_DATA.exists(), "Mn3Sn PG3 500K run not present")
class EstimateRho0GuardCliTests(unittest.TestCase):
    """--estimate-rho0 must refuse data whose criteria are discordant."""

    def test_unconverged_estimate_is_refused(self):
        # Mn3Sn misses O(S(0)) = O(-12) structure below Qmin: the density-limit
        # amplitude is broken, no rho0 reconciles it with the FZ amplitude, and
        # the CLI must error out instead of fitting with a garbage density.
        with tempfile.TemporaryDirectory() as tmp:
            code, _, err = run_cli([
                "--data", MN3SN_500K_DATA, "--qmin", "0.82", "--qmax", "28",
                "--formula", "Mn3Sn", "--rho0", "0.063049",
                "--estimate-rho0", "--out-dir", tmp,
            ])
            self.assertEqual(code, 2)
            self.assertIn("did not converge", err)
            self.assertIn("--amplitude fz", err)


if __name__ == "__main__":
    unittest.main()
