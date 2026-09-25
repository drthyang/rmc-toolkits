# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""STOG readers accept the same inputs as the browser ports (1.0 audit, stog-b group)."""

from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import read_dat_header, read_stog_inp, read_stog_xy

INP = [
    "1", "data.dat", "1.0 28.0", "-9 0.1", "0", "scale.fq", "scale.gr", "50", "5000", "N",
    "0.063049", "0", "N", "Y", "1.0", "scale_ft.sq", "scale_ft.gr", "0.015407",
    "rmc.fq", "rmc.gr", "rmc.dr", "2.48 2.65 3.1",
]
XY = ["        3", "Q S(Q)", " 0.50 1.10 0.01", " 0.51 1.25 0.01", " 0.52 1.30 0.01"]
DAT = [
    "TITLE :: FeCoSn 199K", "NUMBER_DENSITY :: 0.057329 Angstrom^(-3)",
    "MINIMUM_DISTANCES :: 2.4 2.2", " 0.5 1.0",
]


class LineEndingTests(unittest.TestCase):
    """LF, CRLF and CR-only files parse identically (the JS readers now split alike)."""

    def write(self, directory, name, lines, eol):
        path = Path(directory) / name
        path.write_bytes((eol.join(lines) + eol).encode("utf-8"))
        return path

    def test_all_line_endings(self):
        for eol in ("\n", "\r\n", "\r"):
            with self.subTest(eol=repr(eol)), tempfile.TemporaryDirectory() as tmp:
                inp = read_stog_inp(self.write(tmp, "stog.inp", INP, eol))
                self.assertEqual(inp.data_file, "data.dat")
                self.assertAlmostEqual(inp.peak_rmax, 3.1)
                xy = read_stog_xy(self.write(tmp, "data.dat", XY, eol))
                self.assertEqual(xy.shape, (3, 3))
                self.assertEqual(list(xy[1]), [1.1, 1.25, 1.3])
                header = read_dat_header(self.write(tmp, "h.dat", DAT, eol))
                self.assertEqual(header["title"], "FeCoSn 199K")
                self.assertEqual(header["min_distance"], 2.2)


class EncodingTests(unittest.TestCase):
    """Non-UTF-8 bytes in text lines and a UTF-8 BOM are tolerated like the browser.

    The browser decodes uploads with File.text() (UTF-8, invalid bytes replaced,
    BOM stripped); read_stog_xy / read_stog_inp used strict UTF-8, so one latin-1
    'A-ring' in a title line made the CLI and the API refuse a file the page
    scaled. Only numeric tokens are consumed, so replacement characters in text
    lines are harmless.
    """

    def test_latin1_title_line(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "latin1.sq"
            path.write_bytes(
                b"        3\n# Q (\xc5^-1)   S(Q)\n 0.50 1.10\n 0.51 1.25\n 0.52 1.30\n"
            )
            xy = read_stog_xy(path)
            self.assertEqual(list(xy[1]), [1.1, 1.25, 1.3])
            inp = Path(tmp) / "stog.inp"
            lines = list(INP)
            inp.write_bytes("\n".join(lines).encode("utf-8").replace(b"scale.gr", b"scale\xe9.gr"))
            self.assertEqual(read_stog_inp(inp).data_file, "data.dat")
            header = Path(tmp) / "h.dat"
            header.write_bytes(b"TITLE :: 5 \xc5 run\nNUMBER_DENSITY :: 0.05\n")
            self.assertEqual(read_dat_header(header)["number_density"], 0.05)

    def test_utf8_bom(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "stog.inp"
            inp.write_bytes(b"\xef\xbb\xbf" + "\n".join(INP).encode("utf-8"))
            self.assertEqual(read_stog_inp(inp).n_files, 1)
            header = Path(tmp) / "h.dat"
            header.write_bytes(b"\xef\xbb\xbfTITLE :: run\nNUMBER_DENSITY :: 0.05\n")
            self.assertEqual(read_dat_header(header)["title"], "run")

    def test_cli_scales_a_file_with_a_latin1_title(self):
        import contextlib
        import io

        import numpy as np

        from rmc_toolkits.parsers import write_stog_xy
        from rmc_toolkits.scaling_cli import main
        from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

        q = np.arange(20, 981) * 0.03
        r = np.arange(1, 12001) * 0.005
        g = 0.5 * (1 + np.tanh((r - 2.65) / 0.07)) + 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
        sq = (fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, 0.05), q)) + 9.0) / 10.0
        with tempfile.TemporaryDirectory() as tmp:
            clean = Path(tmp) / "clean.sq"
            write_stog_xy(clean, q, sq, title="Q (A^-1)  S(Q)")
            data = Path(tmp) / "latin1.sq"
            data.write_bytes(clean.read_bytes().replace(b"Q (A^-1)", b"Q (\xc5^-1)"))
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = main([
                    "--data", str(data), "--qmin", "0.6", "--qmax", "29.4", "--rho0", "0.05",
                    "--b-avg-sq", "0.02", "--r0", "2.5", "--rmax", "25", "--nr", "1000",
                    "--out-dir", str(Path(tmp) / "out"),
                ])
        self.assertEqual(code, 0, err.getvalue())
        self.assertIn("961 S(Q) points used", out.getvalue())


class NumericTokenTests(unittest.TestCase):
    """read_stog_xy accepts exactly the JS port's numeric tokens (readStogXy NUMERIC_TOKEN)."""

    def test_fortran_d_exponents_and_rejected_forms(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "d.dat"
            path.write_text(
                "3\ntitle\n 0.5D+00 1.1D+00\n 0.51d0 1.25E0\n 5.2E-01 1.3\n"
                " 1_0 2\n 0x10 3\n"
            )
            xy = read_stog_xy(path)
        # D exponents parse (as in the browser); '1_0' and '0x10' are not rows.
        self.assertEqual(list(xy[0]), [0.5, 0.51, 0.52])
        self.assertEqual(list(xy[1]), [1.1, 1.25, 1.3])

    def test_nan_and_inf_spellings_are_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "n.dat"
            path.write_text(" 0.5 NaN\n 0.6 -Infinity\n 0.7 +inf\n")
            xy = read_stog_xy(path)
        self.assertTrue(np.isnan(xy[1][0]))
        self.assertEqual(list(xy[1][1:]), [-np.inf, np.inf])


if __name__ == "__main__":
    unittest.main()
