from pathlib import Path
import os
import tempfile
import unittest

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

from rmc_toolkits.plots import close_plot, detect_plot_kind, make_plot, plot_to_png


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"

# The GNSe example dataset is gitignored (see README), so sample-backed tests
# skip rather than fail when it is not present locally / in CI.
requires_sample = unittest.skipUnless(
    (DATA / "GNSe.rmc6f").exists(),
    "GNSe sample data not present in data/ (gitignored)",
)


class PlotTests(unittest.TestCase):
    def test_detect_plot_kind_for_supported_outputs(self):
        cases = {
            "GNSe_FT_XFQ1.csv": "xpdf",
            "GNSe_FT_XFQ2.csv": "xpdf",
            "GNSe_FQ1.csv": "xray_sq",
            "GNSe_bragg.csv": "bragg",
            "GNSe_PDFpartials.csv": "pdf_partials",
            "Nb-EXAFS-1_Q_OUTPUT.csv": "exafs_q",
            "Nb-EXAFS-1_R_OUTPUT.csv": "exafs_r",
            "Nb-EXAFS-1_OUTPUT.csv": None,
            "GNSe-02.log": "r_value",
            "GNSe-123.log": "r_value",
            "GNSe.log": None,
            "run-info.log": None,
            "scale_ft.gr": "stog",
            "notes.txt": None,
        }

        for filename, expected in cases.items():
            with self.subTest(filename=filename):
                self.assertEqual(detect_plot_kind(filename), expected)

    def test_make_plot_supports_exafs_q_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "Nb-EXAFS-1_Q_OUTPUT.csv"
            path.write_text(
                " EXAFS #1,   chi(k)*k^2\n"
                "      k    ,  calculated  ,  experiment\n"
                " 3.300 ,    -0.32999 ,    -0.23780\n"
                " 3.350 ,    -0.62595 ,    -0.33942\n",
                encoding="utf-8",
            )

            result = make_plot(path)
            try:
                self.assertEqual(result.kind, "exafs_q")
                self.assertEqual(result.title, "EXAFS Q-space")
                png = plot_to_png(result, dpi=72)
                self.assertTrue(png.startswith(b"\x89PNG\r\n\x1a\n"))
            finally:
                close_plot(result)

    def test_make_plot_supports_exafs_r_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "Nb-EXAFS-1_R_OUTPUT.csv"
            path.write_text(
                "     r   ,   Re_Calc  ,  Im_Calc  ,  Mod_Calc  ,   Re_Ex   ,   Im_Ex  ,   Mod_Ex\n"
                "   0.25000 ,    0.04898 ,   -0.29264 ,    0.29671 ,    0.04412 ,   -0.27550 ,    0.27901\n"
                "   0.26000 ,    0.06389 ,   -0.28345 ,    0.29057 ,    0.05774 ,   -0.26676 ,    0.27294\n",
                encoding="utf-8",
            )

            result = make_plot(path)
            try:
                self.assertEqual(result.kind, "exafs_r")
                self.assertEqual(result.title, "EXAFS R-space")
            finally:
                close_plot(result)

    @requires_sample
    def test_make_plot_returns_metadata_and_png_bytes(self):
        result = make_plot(DATA / "GNSe_FQ1.csv")
        try:
            self.assertEqual(result.kind, "xray_sq")
            self.assertEqual(result.title, "S(Q) (x-ray)")
            self.assertIn("rwp", result.metrics)
            self.assertGreater(result.metrics["rwp"], 0.0)

            png = plot_to_png(result, dpi=72)
            self.assertTrue(png.startswith(b"\x89PNG\r\n\x1a\n"))
            self.assertGreater(len(png), 1000)
        finally:
            close_plot(result)

    @requires_sample
    def test_log_plot_reports_final_chi(self):
        result = make_plot(DATA / "GNSe-02.log")
        try:
            self.assertEqual(result.kind, "r_value")
            self.assertAlmostEqual(result.metrics["final_chi_r"], 0.00405)
        finally:
            close_plot(result)

    def test_log_plot_combines_related_logs_in_numeric_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            (directory / "run-10.log").write_text("header\nheader\n1 0.1 10.0\n", encoding="utf-8")
            (directory / "run-01.log").write_text("header\nheader\n1 0.1 1.0\n", encoding="utf-8")
            (directory / "run-02.log").write_text("header\nheader\n1 0.1 2.0\n", encoding="utf-8")

            result = make_plot(directory / "run-01.log")
            try:
                self.assertEqual(result.kind, "r_value")
                self.assertAlmostEqual(result.metrics["final_chi_r"], 10.0)
            finally:
                close_plot(result)


if __name__ == "__main__":
    unittest.main()
