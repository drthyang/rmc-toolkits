# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The classic-named outputs hold the classic stog functions (1.0 audit, stog-b group).

rmc-autoscale writes under the stog.inp-declared names (scale.gr, scale_ft.gr) as
a drop-in replacement for a classic session, but wrote g(r) - 1 in column 2 and
4 pi rho0 r [g - 1] (= G_PDF) in scale_ft.gr's column 3. Every Fortran run in
data/stog_tests has g(r) in column 2 (oscillating about 1) and exactly
r [g(r) - 1] in column 3. The writer (CLI, API) and the page now use the Fortran
conventions.
"""

import contextlib
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import read_stog_inp, read_stog_xy, write_stog_xy
from rmc_toolkits.scaling_cli import main
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

ROOT = Path(__file__).resolve().parents[1]
FECOSN = ROOT / "data" / "stog_tests" / "199K"
FORTRAN_RUNS = [ROOT / "data" / "stog_tests" / name for name in (
    "199K", "stog", "stog_300K", "stog_500K", "stog_59438",
)]


def run_cli(args):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = main([str(arg) for arg in args])
    return code, out.getvalue(), err.getvalue()


class ClassicConventionTests(unittest.TestCase):
    @unittest.skipUnless(all((run / "scale_ft.gr").exists() for run in FORTRAN_RUNS),
                         "Fortran stog reference runs not present")
    def test_fortran_reference_conventions(self):
        for run in FORTRAN_RUNS:
            with self.subTest(run=run.name):
                gr = read_stog_xy(run / "scale.gr")
                ft = read_stog_xy(run / "scale_ft.gr")
                self.assertAlmostEqual(gr[1][gr[0] >= 20].mean(), 1.0, delta=0.03)
                self.assertAlmostEqual(ft[1][ft[0] >= 20].mean(), 1.0, delta=0.03)
                np.testing.assert_array_equal(ft[2], ft[0] * (ft[1] - 1.0))

    def test_written_files_follow_the_classic_conventions(self):
        q = np.arange(20, 981) * 0.03
        r = np.arange(1, 12001) * 0.005
        g = 0.5 * (1 + np.tanh((r - 2.65) / 0.07)) + 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
        sq = (fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, 0.05), q)) + 9.0) / 10.0
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "model.sq"
            write_stog_xy(data, q, sq)
            code, _, err = run_cli([
                "--data", data, "--qmin", "0.6", "--qmax", "29.4", "--rho0", "0.05",
                "--b-avg-sq", "0.02", "--manual", "--scale", "10", "--offset", "-9",
                "--rmax", "25", "--nr", "1000", "--out-dir", Path(tmp) / "out",
            ])
            self.assertEqual(code, 0, err)
            gr = read_stog_xy(Path(tmp) / "out" / "model.gr")
            ft = read_stog_xy(Path(tmp) / "out" / "model_ft.gr")
            gk = read_stog_xy(Path(tmp) / "out" / "model_rmc.gr")
        self.assertAlmostEqual(gr[1][gr[0] >= 10].mean(), 1.0, delta=0.02)
        self.assertAlmostEqual(ft[1][ft[0] >= 10].mean(), 1.0, delta=0.02)
        np.testing.assert_allclose(ft[2], ft[0] * (ft[1] - 1.0), rtol=1e-12, atol=1e-15)
        # The RMC file stays Keen's G_K = <b>^2 (g - 1) (enforced below the cutoff).
        above = gk[0] > 3.0
        np.testing.assert_allclose(gk[1][above], 0.02 * (ft[1][above] - 1.0), atol=1e-12)

    @unittest.skipUnless((FECOSN / "stog_input.dat").exists(), "FeCoSn 199K run not present")
    def test_manual_run_matches_the_fortran_files(self):
        inp = read_stog_inp(FECOSN / "stog_input.dat")
        with tempfile.TemporaryDirectory() as tmp:
            code, _, err = run_cli([FECOSN / "stog_input.dat", "--manual", "--out-dir", tmp])
            self.assertEqual(code, 0, err)
            ours_gr = read_stog_xy(Path(tmp) / inp.out_gr)
            ours_ft = read_stog_xy(Path(tmp) / inp.out_ft_gr)
        ref_gr = read_stog_xy(FECOSN / "scale.gr")
        ref_ft = read_stog_xy(FECOSN / "scale_ft.gr")
        rms = lambda a, b: float(np.sqrt(np.mean((a - b) ** 2)))  # noqa: E731
        self.assertLess(rms(np.interp(ref_gr[0], ours_gr[0], ours_gr[1]), ref_gr[1]), 0.01)
        for column in (1, 2):
            ours = np.interp(ref_ft[0], ours_ft[0], ours_ft[column])
            self.assertLess(rms(ours, ref_ft[column]), 0.01)


if __name__ == "__main__":
    unittest.main()
