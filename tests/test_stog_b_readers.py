# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""STOG readers accept the same inputs as the browser ports (1.0 audit, stog-b group)."""

from pathlib import Path
import tempfile
import unittest

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


if __name__ == "__main__":
    unittest.main()
