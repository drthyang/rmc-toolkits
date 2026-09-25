# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Parser, plot and API tests on the COMMITTED demo run (web_app/frontend/public/demo).

The GNSe reference dataset the sample-backed tests use is gitignored, so those
tests skip in CI. These are their counterparts on real RMCProfile files that are
always present (GaTa4Se8 at 250 K: GTS_250K.rmc6f, three -NN.log restarts, the
F(Q), xPDF and partials CSVs). Every expected value is derived here, from the
files' own text with the csv module / plain string splitting — never copied from
the parsers under test.
"""

from __future__ import annotations

import csv
import math
import os
from pathlib import Path
import re
import shutil
import sys
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import (
    frac_lines_from_rmc6f,
    iter_rmc6f_atoms,
    read_atom_indices,
    read_cell_vectors,
    read_chi,
    read_moves_metadata,
    read_rmc_csv,
    read_structure,
    related_r_value_logs,
    write_frac_from_rmc6f,
)
from rmc_toolkits.plots import close_plot, make_plot, plot_to_png

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"
RMC6F = DEMO / "GTS_250K.rmc6f"
LOGS = [DEMO / f"GTS_250K-0{index}.log" for index in range(3)]


# --- independent readers of the demo files -----------------------------------

def csv_rows(path: Path) -> tuple[list[str], list[list[float]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.reader(handle) if any(cell.strip() for cell in row)]
    header = [cell.strip() for cell in rows[0]]
    data = [[float(cell) for cell in row if cell.strip()] for row in rows[1:]]
    return header, data


def header_value(key: str) -> list[str]:
    for line in RMC6F.read_text(encoding="utf-8").splitlines():
        if line.startswith(key):
            return line.split(":", 1)[1].split()
    raise KeyError(key)


def atom_lines() -> list[list[str]]:
    lines = RMC6F.read_text(encoding="utf-8").splitlines()
    start = lines.index("Atoms:") + 1
    return [line.split() for line in lines[start:] if line.strip()]


def log_last_column(path: Path) -> list[float]:
    return [float(line.split()[-1]) for line in path.read_text(encoding="utf-8").splitlines()[2:] if line.strip()]


def conventional_r(calc: list[float], expt: list[float]) -> float:
    return math.sqrt(sum((c - e) ** 2 for c, e in zip(calc, expt)) / sum(e * e for e in expt))


ATOM_TYPES = header_value("Atom types present")                     # Ga Ta Se
TYPE_COUNTS = [int(v) for v in header_value("Number of each atom type")]
DECLARED = int(header_value("Number of atoms")[0])
SUPERCELL = [float(v) for v in header_value("Supercell dimensions")]


class DemoRunParserTests(unittest.TestCase):
    def test_fit_csv_labels_shape_and_first_row(self):
        header, rows = csv_rows(DEMO / "GTS_250K_FQ1.csv")
        series = read_rmc_csv(DEMO / "GTS_250K_FQ1.csv")
        self.assertEqual(series.labels, header)
        self.assertEqual(series.labels, ["Q", "F(Q)_RMC", "F(Q)_Expt"])
        self.assertEqual(series.data.shape, (3, len(rows)))
        np.testing.assert_allclose(series.data[:, 0], rows[0])
        np.testing.assert_allclose(series.data[:, -1], rows[-1])

    def test_partials_with_trailing_commas(self):
        header, rows = csv_rows(DEMO / "GTS_250K_PDFpartials.csv")
        series = read_rmc_csv(DEMO / "GTS_250K_PDFpartials.csv")
        self.assertEqual(series.labels, header)
        self.assertEqual(series.data.shape, (len(header), len(rows)))

    def test_log_restarts_are_concatenated_in_order(self):
        logs = related_r_value_logs(LOGS[1])
        self.assertEqual([path.name for path in logs], [path.name for path in LOGS])
        chi_q, chi_r = read_chi(logs)
        expected = [value for path in LOGS for value in log_last_column(path)]
        np.testing.assert_array_equal(chi_r, expected)
        self.assertEqual(len(chi_q), len(expected))

    def test_rmc6f_metadata_atom_indices_and_counts(self):
        lattice, supercell = read_cell_vectors(RMC6F)
        np.testing.assert_array_equal(supercell, SUPERCELL)
        cell = [float(v) for v in header_value("Cell (Ang/deg)")]
        np.testing.assert_allclose(np.linalg.norm(lattice, axis=1), cell[:3])

        sites: dict[str, set[int]] = {}
        for parts in atom_lines():
            sites.setdefault(parts[1], set()).add(int(parts[6]))
        self.assertEqual(read_atom_indices(RMC6F), {el: sorted(refs) for el, refs in sites.items()})

        atoms = list(iter_rmc6f_atoms(RMC6F))
        self.assertEqual(len(atoms), DECLARED)
        counts = {el: sum(1 for atom in atoms if atom["element"] == el) for el in ATOM_TYPES}
        self.assertEqual([counts[el] for el in ATOM_TYPES], TYPE_COUNTS)

    def test_moves_metadata_from_the_header(self):
        moves = read_moves_metadata(RMC6F)
        self.assertEqual(moves["generated"], float(header_value("Number of moves generated")[0]))
        self.assertEqual(moves["tried"], float(header_value("Number of moves tried")[0]))
        self.assertEqual(moves["accepted"], float(header_value("Number of moves accepted")[0]))

    def test_frac_conversion_shape_and_first_atom(self):
        lines = frac_lines_from_rmc6f(RMC6F)
        self.assertEqual(len(lines), 5 + DECLARED)
        first = atom_lines()[0]
        box = [float(v) for v in first[3:6]]
        cells = [int(v) for v in first[7:10]]
        reduced = [box[i] - cells[i] / SUPERCELL[i] for i in range(3)]
        expected = (
            f"{int(first[6]):3d}    {reduced[0]:.5f}    {reduced[1]:.5f}    {reduced[2]:.5f}  "
            f"{cells[0]:d}  {cells[1]:d}  {cells[2]:d}\n"
        )
        self.assertEqual(lines[5], expected)

    def test_read_structure_on_a_copy_with_its_frac_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            shutil.copy(RMC6F, directory / RMC6F.name)
            write_frac_from_rmc6f(directory / RMC6F.name)
            structure = read_structure(directory)
            gallium = read_structure(directory, element="Ga")
        edge = float(header_value("Cell (Ang/deg)")[0])
        self.assertEqual(structure.positions.shape, (DECLARED, 3))
        self.assertEqual(len(set(structure.atom_types)), sum(len(v) for v in read_atom_indices(RMC6F).values()))
        self.assertTrue(np.all(structure.positions >= 0.0))
        self.assertTrue(np.all(structure.positions <= edge / SUPERCELL[0] + 1e-9))
        self.assertEqual(len(gallium.positions), TYPE_COUNTS[ATOM_TYPES.index("Ga")])


class DemoRunPlotTests(unittest.TestCase):
    def test_fq_plot_metadata_rwp_and_png(self):
        _, rows = csv_rows(DEMO / "GTS_250K_FQ1.csv")
        expected = conventional_r([row[1] for row in rows], [row[2] for row in rows])
        result = make_plot(DEMO / "GTS_250K_FQ1.csv")
        try:
            self.assertEqual((result.kind, result.title), ("xray_sq", "F(Q)"))
            self.assertAlmostEqual(result.metrics["rwp"], expected, places=12)
            png = plot_to_png(result, dpi=72)
            self.assertTrue(png.startswith(b"\x89PNG\r\n\x1a\n"))
        finally:
            close_plot(result)

    def test_xpdf_rwp_uses_the_experiment(self):
        header, rows = csv_rows(DEMO / "GTS_250K_FT_XFQ1.csv")
        self.assertEqual(header, ["r(A)", "X_ray-calc", "X_ray_exp_renorm"])
        expected = conventional_r([row[1] for row in rows], [row[2] for row in rows])
        result = make_plot(DEMO / "GTS_250K_FT_XFQ1.csv")
        try:
            self.assertEqual(result.kind, "xpdf")
            self.assertAlmostEqual(result.metrics["rwp"], expected, places=12)
        finally:
            close_plot(result)

    def test_log_plot_final_chi_is_the_last_row_of_the_last_restart(self):
        result = make_plot(LOGS[0])
        try:
            self.assertEqual(result.kind, "r_value")
            self.assertEqual(result.metrics["final_chi_r"], log_last_column(LOGS[-1])[-1])
        finally:
            close_plot(result)


if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402

DEMO_REL = DEMO.relative_to(ROOT).as_posix()


class DemoRunApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

    def test_files_lists_the_demo_outputs_with_their_kinds(self):
        response = self.client.get("/api/files", query_string={"dir": DEMO_REL})
        self.assertEqual(response.status_code, 200)
        kinds = {item["name"]: item["plotKind"] for item in response.get_json()["files"]}
        self.assertEqual(kinds["GTS_250K_FQ1.csv"], "xray_sq")
        self.assertEqual(kinds["GTS_250K_FT_XFQ1.csv"], "xpdf")
        self.assertEqual(kinds["GTS_250K_PDFpartials.csv"], "pdf_partials")
        self.assertEqual(kinds["GTS_250K-00.log"], "r_value")
        self.assertIsNone(kinds["GTS_250K_XFQ1.csv"])
        self.assertIsNone(kinds["GTS_250K_FQ1partials.csv"])

    def test_plot_metadata_and_data_for_the_fq_fit(self):
        _, rows = csv_rows(DEMO / "GTS_250K_FQ1.csv")
        expected = conventional_r([row[1] for row in rows], [row[2] for row in rows])
        path = f"{DEMO_REL}/GTS_250K_FQ1.csv"
        metadata = self.client.get("/api/plot/metadata", query_string={"path": path}).get_json()
        data = self.client.get("/api/plot/data", query_string={"path": path}).get_json()
        self.assertEqual((metadata["kind"], metadata["title"]), ("xray_sq", "F(Q)"))
        self.assertAlmostEqual(metadata["metrics"]["rwp"], expected, places=12)
        self.assertEqual(data["xLabel"], "Q (Å^{-1})")
        self.assertEqual(len(data["series"]), 2)
        self.assertEqual(len(data["series"][0]["x"]), len(rows))

    def test_log_data_concatenates_the_three_restarts(self):
        data = self.client.get("/api/plot/data", query_string={"path": f"{DEMO_REL}/GTS_250K-01.log"}).get_json()
        expected = [value for path in LOGS for value in log_last_column(path)]
        self.assertEqual(len(data["series"][0]["y"]), len(expected))
        self.assertAlmostEqual(data["series"][0]["y"][-1], math.log(expected[-1]), places=12)
        self.assertEqual(data["metrics"]["final_chi_r"], expected[-1])
        self.assertEqual(data["title"], "χ² history: X_ray_(R)1")

    def test_structure_endpoint_reports_the_header_composition(self):
        response = self.client.get("/api/structure", query_string={"dir": DEMO_REL, "maxPoints": 500})
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["totalAtoms"], DECLARED)
        self.assertEqual(payload["elements"], sorted(ATOM_TYPES))
        self.assertEqual({el: payload["elementCounts"][el] for el in ATOM_TYPES}, dict(zip(ATOM_TYPES, TYPE_COUNTS)))
        self.assertLessEqual(payload["sampledAtoms"], 500)
        self.assertIsNone(payload["parseWarning"])
        self.assertEqual(payload["moves"]["accepted"], float(header_value("Number of moves accepted")[0]))

    def test_convert_frac_writes_inside_the_data_root(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmpdir:
            output = Path(tmpdir) / "Frac_coord_demo.txt"
            response = self.client.post(
                "/api/convert/frac",
                json={"path": f"{DEMO_REL}/GTS_250K.rmc6f", "outputPath": str(output)},
            )
            self.assertEqual(response.status_code, 200)
            self.assertEqual(len(output.read_text(encoding="utf-8").splitlines()), 5 + DECLARED)
            self.assertTrue(re.match(r"^\s*\d+\s", output.read_text(encoding="utf-8").splitlines()[5]))


if __name__ == "__main__":
    unittest.main()
