# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

from pathlib import Path
import os
import tempfile
import unittest

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

from rmc_toolkits.plots import bragg_is_tof, close_plot, detect_plot_kind, make_plot, plot_to_png


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
            "GNSe_bragg_1.csv": "bragg",
            "GNSe_braggish.csv": None,
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


def _write_fit_csv(directory: Path, name: str, header: str, rows: list[tuple[float, ...]]) -> Path:
    path = directory / name
    path.write_text(
        header + "\n" + "".join(", ".join(f"{value:.7f}" for value in row) + "\n" for row in rows),
        encoding="utf-8",
    )
    return path


class RwpColumnRoleTests(unittest.TestCase):
    """The dashboard R-factor is normalized by the EXPERIMENT, whatever the column order.

    RMCProfile writes fit CSVs as (x, calculated, experimental). With the calculated
    curve a uniform 0.7 x the experiment, the conventional R = ||calc - expt|| / ||expt||
    is exactly 0.3; normalizing by the calculated column instead gives 0.3 / 0.7.
    """

    EXPT = (1.0, -2.0, 3.0, 0.5, -1.5)

    def _rwp(self, name: str, header: str, calc_first: bool = True) -> float:
        rows = []
        for index, expt in enumerate(self.EXPT):
            calc = 0.7 * expt
            rows.append((0.1 * index, calc, expt) if calc_first else (0.1 * index, expt, calc))
        with tempfile.TemporaryDirectory() as tmpdir:
            result = make_plot(_write_fit_csv(Path(tmpdir), name, header, rows))
            try:
                return result.metrics["rwp"]
            finally:
                close_plot(result)

    def test_rmcprofile_order_divides_by_the_experiment(self):
        for name, header in (
            ("run_FQ1.csv", "Q, F(Q)_RMC, F(Q)_Expt"),
            ("run_FT_XFQ1.csv", "r(A), X_ray-calc, X_ray_exp_renorm"),
            ("run_PDF1.csv", "r, G(r)_RMC, G(r)_Expt"),
            ("run_SQ1.csv", "Q, S(Q)_RMC, S(Q)_Expt"),
            ("run_bragg.csv", "Flight time (us), Calculated, Experiment"),
        ):
            with self.subTest(name=name):
                self.assertAlmostEqual(self._rwp(name, header), 0.3, places=6)

    def test_unlabelled_columns_follow_the_rmcprofile_positional_order(self):
        # No role names in the header: column 2 is the calculation, column 3 the data.
        self.assertAlmostEqual(self._rwp("run_FQ1.csv", "Q, a, b"), 0.3, places=6)

    def test_header_roles_override_the_positional_order(self):
        # A file written (x, experimental, calculated) is still normalized by the data.
        self.assertAlmostEqual(
            self._rwp("run_FQ1.csv", "Q, F(Q)_Expt, F(Q)_RMC", calc_first=False), 0.3, places=6
        )


DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"


class ChiHistoryLabelTests(unittest.TestCase):
    """The log series is ONE column's chi^2, named by its header — not a total "R-value"."""

    def test_demo_run_logs_are_labelled_by_their_last_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            for name in ("GTS_250K-00.log", "GTS_250K-01.log", "GTS_250K-02.log"):
                (directory / name).write_bytes((DEMO / name).read_bytes())
            header = (DEMO / "GTS_250K-00.log").read_text(encoding="utf-8").splitlines()[0].split()
            last = (DEMO / "GTS_250K-02.log").read_text(encoding="utf-8").splitlines()[-1].split()[-1]

            result = make_plot(directory / "GTS_250K-01.log")
            try:
                self.assertEqual(result.kind, "r_value")
                self.assertEqual(result.title, f"χ² history: {header[-1]}")
                self.assertEqual(result.title, "χ² history: X_ray_(R)1")
                self.assertEqual(result.metrics["final_chi_r"], float(last))
                self.assertEqual(result.figure.axes[0].get_legend_handles_labels()[1], ["X_ray_(R)1"])
            finally:
                close_plot(result)

    def test_headerless_log_says_it_is_the_last_column(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "run-00.log"
            path.write_text("header\nheader\n1 0.1 2.0\n", encoding="utf-8")
            result = make_plot(path)
            try:
                self.assertEqual(result.title, "χ² history: last log column")
            finally:
                close_plot(result)


class BraggAxisTests(unittest.TestCase):
    def test_time_of_flight_headers(self):
        for header in ("Flight time (us)", "TOF,ms", "Time"):
            self.assertTrue(bragg_is_tof(header))

    def test_reciprocal_space_is_not_tof(self):
        for header in ("Q or theta", "2-theta, deg", "", None):
            self.assertFalse(bragg_is_tof(header))


if __name__ == "__main__":
    unittest.main()
