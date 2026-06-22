from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import (
    frac_lines_from_rmc6f,
    read_atom_indices,
    read_cell_vectors,
    read_chi,
    read_exafs_csv,
    read_rmc_csv,
    read_structure,
    related_r_value_logs,
    rwp,
    write_frac_from_rmc6f,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

# The GNSe example dataset is gitignored (see README), so sample-backed tests
# skip rather than fail when it is not present locally / in CI.
requires_sample = unittest.skipUnless(
    (DATA / "GNSe.rmc6f").exists(),
    "GNSe sample data not present in data/ (gitignored)",
)


class ParserTests(unittest.TestCase):
    @requires_sample
    def test_read_rmc_csv_loads_labels_and_numeric_columns(self):
        series = read_rmc_csv(DATA / "GNSe_FQ1.csv")

        self.assertEqual(series.labels, ["Q", "F(Q)_RMC", "F(Q)_Expt"])
        self.assertEqual(series.data.shape, (3, 2649))
        np.testing.assert_allclose(series.data[:, 0], [0.5, -1.0350516, -0.6322244])

    def test_read_exafs_csv_skips_q_output_title_row(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "Nb-EXAFS-1_Q_OUTPUT.csv"
            path.write_text(
                " EXAFS #1,   chi(k)*k^2\n"
                "      k    ,  calculated  ,  experiment\n"
                " 3.300 ,    -0.32999 ,    -0.23780\n"
                " 3.350 ,    -0.62595 ,    -0.33942\n",
                encoding="utf-8",
            )

            series = read_exafs_csv(path)

        self.assertEqual(series.labels, ["k", "calculated", "experiment"])
        self.assertEqual(series.data.shape, (3, 2))
        np.testing.assert_allclose(series.data[:, 0], [3.3, -0.32999, -0.2378])

    def test_read_exafs_csv_accepts_r_output_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "Nb-EXAFS-1_R_OUTPUT.csv"
            path.write_text(
                "     r   ,   Re_Calc  ,  Im_Calc  ,  Mod_Calc  ,   Re_Ex   ,   Im_Ex  ,   Mod_Ex\n"
                "   0.25000 ,    0.04898 ,   -0.29264 ,    0.29671 ,    0.04412 ,   -0.27550 ,    0.27901\n",
                encoding="utf-8",
            )

            series = read_exafs_csv(path)

        self.assertEqual(series.labels, ["r", "Re_Calc", "Im_Calc", "Mod_Calc", "Re_Ex", "Im_Ex", "Mod_Ex"])
        self.assertEqual(series.data.shape, (7, 1))
        np.testing.assert_allclose(series.data[:, 0], [0.25, 0.04898, -0.29264, 0.29671, 0.04412, -0.2755, 0.27901])

    @requires_sample
    def test_read_chi_extracts_q_and_r_columns(self):
        chi_q, chi_r = read_chi([DATA / "GNSe-02.log"])

        self.assertEqual(len(chi_q), 571)
        self.assertEqual(len(chi_r), 571)
        self.assertAlmostEqual(float(chi_q[0]), 0.00541)
        self.assertAlmostEqual(float(chi_r[-1]), 0.00405)

    def test_related_r_value_logs_use_numeric_suffix_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            for name in ("run-10.log", "run-02.log", "run-01.log", "other-01.log"):
                (directory / name).write_text("header\nheader\n1 0.1 0.2\n", encoding="utf-8")

            logs = related_r_value_logs(directory / "run-02.log")

            self.assertEqual([path.name for path in logs], ["run-01.log", "run-02.log", "run-10.log"])

    def test_rwp_uses_observed_series_as_denominator(self):
        value = rwp(
            np.array([0.0, 1.0, 2.0]),
            np.array([2.0, 4.0, 4.0]),
            np.array([1.0, 5.0, 6.0]),
        )

        self.assertAlmostEqual(value, np.sqrt(6.0 / 36.0))

    @requires_sample
    def test_read_rmc6f_metadata_and_atom_indices(self):
        atom_indices = read_atom_indices(DATA / "GNSe.rmc6f")
        lattice_vectors, supercell = read_cell_vectors(DATA / "GNSe.rmc6f")

        self.assertEqual(atom_indices["Ga"], [4, 8, 12, 16])
        self.assertEqual(len(atom_indices["Nb"]), 16)
        self.assertEqual(len(atom_indices["Se"]), 32)
        np.testing.assert_allclose(supercell, [10.0, 10.0, 10.0])
        np.testing.assert_allclose(np.diag(lattice_vectors), [104.116, 104.116, 104.116])

    @requires_sample
    def test_frac_conversion_matches_sample_shape_and_first_atom(self):
        lines = frac_lines_from_rmc6f(DATA / "GNSe.rmc6f")

        self.assertEqual(len(lines), 52005)
        self.assertEqual(lines[5], "  4    0.07507    0.07511    0.07275  0  0  0\n")

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "Frac_coord_GNSe.txt"
            written = write_frac_from_rmc6f(DATA / "GNSe.rmc6f", output, overwrite=False)
            self.assertEqual(written, output)
            self.assertEqual(output.read_text(encoding="utf-8").splitlines()[:8], [line.rstrip("\n") for line in lines[:8]])

    @requires_sample
    def test_read_structure_loads_full_folded_unit_cell(self):
        structure = read_structure(DATA)

        self.assertEqual(structure.positions.shape, (52000, 3))
        self.assertEqual(sorted(structure.atom_indices), ["Ga", "Nb", "Se"])
        self.assertEqual(len(set(structure.atom_types)), 52)
        self.assertTrue(np.all(structure.positions >= 0.0))
        self.assertTrue(np.all(structure.positions <= 10.4116))

    @requires_sample
    def test_read_structure_rejects_unknown_element_and_mode(self):
        with self.assertRaisesRegex(ValueError, "Unknown element"):
            read_structure(DATA, element="Zn")

        with self.assertRaisesRegex(ValueError, "mode must be"):
            read_structure(DATA, mode="bad-mode")


if __name__ == "__main__":
    unittest.main()
