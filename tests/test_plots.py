from pathlib import Path
import os
import tempfile
import unittest

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

from rmc_toolkits.plots import close_plot, detect_plot_kind, make_plot, plot_to_png


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


class PlotTests(unittest.TestCase):
    def test_detect_plot_kind_for_supported_outputs(self):
        cases = {
            "GNSe_FT_XFQ1.csv": "xpdf",
            "GNSe_FQ1.csv": "xray_sq",
            "GNSe_bragg.csv": "bragg",
            "GNSe_PDFpartials.csv": "pdf_partials",
            "GNSe-02.log": "r_value",
            "scale_ft.gr": "stog",
            "notes.txt": None,
        }

        for filename, expected in cases.items():
            with self.subTest(filename=filename):
                self.assertEqual(detect_plot_kind(filename), expected)

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

    def test_log_plot_reports_final_chi(self):
        result = make_plot(DATA / "GNSe-02.log")
        try:
            self.assertEqual(result.kind, "r_value")
            self.assertAlmostEqual(result.metrics["final_chi_r"], 0.00405)
        finally:
            close_plot(result)


if __name__ == "__main__":
    unittest.main()
