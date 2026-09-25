# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The .rmc6f atom-line grammar, pinned on real-file variants.

Each variant is built from a real RMCProfile configuration — the GaNb4Se8 5 K run
in data/5K_try1 when it is present locally, else the committed demo run
web_app/frontend/public/demo/GTS_250K.rmc6f — cut to its first ATOMS atoms with
the header's ``Number of atoms`` adjusted. Expected values come from the base
file itself (the header count and the first line's own tokens), never from the
parser under test. The browser parser is pinned on the same variants with the
same expectations in web_app/frontend/src/__tests__/rmc6fGrammar.test.js.
"""

from __future__ import annotations

import os
from pathlib import Path
import re
import sys
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import (
    Rmc6fParseReport,
    iter_rmc6f_atoms,
    parse_fortran_number,
    parse_rmc6f_atoms,
    read_atom_indices,
    read_moves_metadata,
)

ROOT = Path(__file__).resolve().parents[1]
REAL_5K = ROOT / "data" / "5K_try1" / "GaNb4Se8_5K.rmc6f"
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo" / "GTS_250K.rmc6f"
ATOMS = 300


def _base() -> tuple[list[str], list[str]]:
    """(header lines incl. the Atoms marker, first ATOMS atom lines) of a real file."""
    source = REAL_5K if REAL_5K.exists() else DEMO
    header: list[str] = []
    atoms: list[str] = []
    in_atoms = False
    with source.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not in_atoms:
                header.append(re.sub(r"(Number of atoms:\s*)\d+", rf"\g<1>{ATOMS}", line))
                in_atoms = line.strip() == "Atoms:"
                continue
            atoms.append(line)
            if len(atoms) == ATOMS:
                break
    assert header[-1].strip() == "Atoms:" and len(atoms) == ATOMS
    return header, atoms


HEADER, ATOM_LINES = _base()
FIRST = ATOM_LINES[0].split()          # id el [1] x y z ref cx cy cz
FIRST_COORDS = [float(token) for token in FIRST[3:6]]
FIRST_REF = int(FIRST[6])
SUPERCELL = next(
    [int(float(v)) for v in line.split()[-3:]] for line in HEADER if line.startswith("Supercell")
)


def _with_atoms(atom_lines: list[str], header: list[str] | None = None, newline: str = "\n") -> str:
    return newline.join((header or HEADER) + atom_lines) + newline


def _edit_tokens(edit) -> list[str]:
    return ["   ".join(edit(line.split())) for line in ATOM_LINES]


def _fortran_d(value: str) -> str:
    mantissa, exponent = f"{float(value):.14E}".split("E")
    return f"{mantissa}D{exponent}"


VARIANTS: dict[str, str] = {
    "crlf": _with_atoms(ATOM_LINES, newline="\r\n"),
    "cr_only": _with_atoms(ATOM_LINES, newline="\r"),
    "tabs": _with_atoms(["\t".join(line.split()) for line in ATOM_LINES]),
    "bom": "\uFEFF" + _with_atoms(ATOM_LINES),
    "no_label": _with_atoms(_edit_tokens(lambda t: t[:2] + t[3:])),
    "split_label": _with_atoms(_edit_tokens(lambda t: t[:2] + ["[", t[2][1:]] + t[3:])),
    "e_notation": _with_atoms(_edit_tokens(lambda t: t[:3] + [f"{float(v):.15E}" for v in t[3:6]] + t[6:])),
    "fortran_d": _with_atoms(_edit_tokens(lambda t: t[:3] + [_fortran_d(v) for v in t[3:6]] + t[6:])),
    "trailing_blank_lines": _with_atoms(ATOM_LINES + ["", "   ", ""]),
    "marker_space": _with_atoms(ATOM_LINES, HEADER[:-1] + ["Atoms :"]),
    "marker_lower": _with_atoms(ATOM_LINES, HEADER[:-1] + ["atoms:"]),
    "marker_suffix": _with_atoms(ATOM_LINES, HEADER[:-1] + ["Atoms (fractional coordinates):"]),
    "upper_element": _with_atoms(_edit_tokens(lambda t: t[:1] + [t[1].upper()] + t[2:])),
    "extra_numeric_column": _with_atoms([line + "   2.500000" for line in ATOM_LINES]),
    "trailing_moment_token": _with_atoms([line + "   M:  2.500000" for line in ATOM_LINES]),
    "label_without_reference": _with_atoms(_edit_tokens(lambda t: t[:6] + t[7:])),
    "coords_only": _with_atoms(_edit_tokens(lambda t: t[:6])),
}

CLEAN = (
    "crlf", "cr_only", "tabs", "bom", "no_label", "split_label", "e_notation",
    "fortran_d", "trailing_blank_lines", "marker_space", "marker_lower",
    "marker_suffix", "upper_element",
)
ALL_INVALID = ("extra_numeric_column", "trailing_moment_token", "label_without_reference")


def _write(directory: Path, name: str, text: str) -> Path:
    path = directory / f"{name}.rmc6f"
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(text)
    return path


class Rmc6fGrammarVariantTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.directory = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _parse(self, name: str, text: str, **kwargs):
        report = Rmc6fParseReport()
        atoms = list(iter_rmc6f_atoms(_write(self.directory, name, text), report=report, **kwargs))
        return atoms, report

    def test_layout_variants_parse_every_atom_identically(self):
        for name in CLEAN:
            with self.subTest(variant=name):
                atoms, report = self._parse(name, VARIANTS[name])
                self.assertTrue(report.has_atoms_section)
                self.assertEqual(report.declared_atoms, ATOMS)
                self.assertEqual(report.parsed_atoms, ATOMS)
                self.assertEqual((report.invalid_lines, report.non_finite_lines), (0, 0))
                self.assertIsNone(report.warning())
                self.assertEqual(len(atoms), ATOMS)
                np.testing.assert_allclose(atoms[0]["coords"], FIRST_COORDS, rtol=0, atol=1e-12)
                self.assertEqual(atoms[0]["reference_number"], FIRST_REF)
                self.assertEqual(atoms[0]["element"], FIRST[1].capitalize())

    def test_unknown_layouts_are_reported_not_shifted_or_silently_dropped(self):
        # An extra trailing field used to shift every column in the browser (y,z -> x,y;
        # a cell index became the reference number) and drop every atom in Python while
        # read_atom_indices still reported cell indices as "sites".
        for name in ALL_INVALID:
            with self.subTest(variant=name):
                path = _write(self.directory, name, VARIANTS[name])
                atoms, report = self._parse(name, VARIANTS[name])
                self.assertEqual(atoms, [])
                self.assertEqual(report.invalid_lines, ATOMS)
                self.assertEqual(report.atom_lines, ATOMS)
                self.assertIn(f"{ATOMS} of {ATOMS} atom lines unparsed", report.warning())
                self.assertIn(f"parsed 0 of {ATOMS} atoms declared", report.warning())
                self.assertEqual(read_atom_indices(path), {})
                with self.assertRaisesRegex(ValueError, "no atoms could be parsed.*unparsed"):
                    parse_rmc6f_atoms(path)

    def test_coords_only_lines_are_recognized_and_opt_in(self):
        atoms, report = self._parse("coords_only", VARIANTS["coords_only"])
        self.assertEqual(atoms, [])  # default: full layout only (reference/cell always set)
        self.assertEqual(report.coords_only_atoms, ATOMS)
        self.assertIsNone(report.warning())
        atoms, _ = self._parse("coords_only", VARIANTS["coords_only"], include_coords_only=True)
        self.assertEqual(len(atoms), ATOMS)
        self.assertIsNone(atoms[0]["reference_number"])
        self.assertIsNone(atoms[0]["cell_indices"])
        np.testing.assert_allclose(atoms[0]["coords"], FIRST_COORDS, rtol=0, atol=1e-12)

    def test_truncated_file_reports_the_shortfall(self):
        # A Live Data read that lands mid-write: the header promises ATOMS atoms.
        text = _with_atoms(ATOM_LINES)
        cut = text[: int(len(text) * 0.6)]
        atoms, report = self._parse("truncated", cut)
        self.assertLess(len(atoms), ATOMS)
        self.assertIn(f"parsed {len(atoms)} of {ATOMS} atoms declared in the header", report.warning())

    def test_non_finite_coordinate_lines_are_skipped_and_counted(self):
        lines = list(ATOM_LINES)
        tokens = lines[5].split()
        tokens[3] = "NaN"
        lines[5] = "   ".join(tokens)
        tokens = lines[9].split()
        tokens[4] = "*********"
        lines[9] = "   ".join(tokens)
        atoms, report = self._parse("nan", _with_atoms(lines))
        self.assertEqual(len(atoms), ATOMS - 2)
        self.assertEqual(report.non_finite_lines, 2)
        self.assertTrue(all(np.all(np.isfinite(atom["coords"])) for atom in atoms))
        self.assertIn("2 atom lines skipped for non-finite coordinates", report.warning())

    def test_reference_and_cell_indices_are_validated(self):
        lines = list(ATOM_LINES)
        bad_cell = lines[1].split()
        bad_cell[7] = str(SUPERCELL[0])      # cell index N_x is outside [0, N_x)
        bad_ref = lines[2].split()
        bad_ref[6] = "0"                     # reference numbers start at 1
        bad_int = lines[3].split()
        bad_int[8] = "1.5"                   # cell indices are integers
        lines[1:4] = ["   ".join(bad_cell), "   ".join(bad_ref), "   ".join(bad_int)]
        atoms, report = self._parse("bad", _with_atoms(lines))
        self.assertEqual(len(atoms), ATOMS - 3)
        self.assertEqual(report.invalid_lines, 3)

    def test_moves_metadata_uses_the_same_marker_rule(self):
        path = _write(self.directory, "marker_lower", VARIANTS["marker_lower"])
        self.assertEqual(read_moves_metadata(path), read_moves_metadata(
            _write(self.directory, "crlf", VARIANTS["crlf"])
        ))


class FortranNumberTests(unittest.TestCase):
    def test_accepts_plain_e_and_d_exponents(self):
        self.assertEqual(parse_fortran_number("0.117D-03"), 0.117e-3)
        self.assertEqual(parse_fortran_number("0.117E-03"), 0.117e-3)
        self.assertEqual(parse_fortran_number("-.5"), -0.5)
        self.assertEqual(parse_fortran_number("12"), 12.0)

    def test_non_finite_and_invalid_tokens(self):
        for token in ("NaN", "nan", "Inf", "-Infinity", "********"):
            with self.subTest(token=token):
                self.assertTrue(np.isnan(parse_fortran_number(token)))
        for token in ("abc", "1_000", "0x10", "1.2.3", "", "M:"):
            with self.subTest(token=token):
                self.assertIsNone(parse_fortran_number(token))


# --- Flask /api/structure -----------------------------------------------------

if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402


class StructureEndpointReportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(dir=ROOT)
        self.directory = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _structure(self, text: str):
        _write(self.directory, "run", text)
        return self.client.get("/api/structure", query_string={"dir": str(self.directory)})

    def test_clean_file_reports_no_warning(self):
        response = self._structure(VARIANTS["fortran_d"])
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["totalAtoms"], ATOMS)
        self.assertIsNone(payload["parseWarning"])
        self.assertEqual(payload["parseReport"]["parsedAtoms"], ATOMS)

    def test_coords_only_file_counts_its_atoms_like_the_browser(self):
        response = self._structure(VARIANTS["coords_only"])
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["totalAtoms"], ATOMS)
        self.assertEqual(payload["atomIndices"], {})
        self.assertEqual(payload["parseReport"]["coordsOnlyAtoms"], ATOMS)

    def test_zero_parsed_atoms_is_an_error_that_says_why(self):
        response = self._structure(VARIANTS["extra_numeric_column"])
        self.assertEqual(response.status_code, 500)
        self.assertIn(f"{ATOMS} of {ATOMS} atom lines unparsed", response.get_json()["error"])

    def test_non_finite_line_is_skipped_and_reported(self):
        lines = list(ATOM_LINES)
        tokens = lines[0].split()
        tokens[5] = "NaN"
        lines[0] = "   ".join(tokens)
        response = self._structure(_with_atoms(lines))
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["totalAtoms"], ATOMS - 1)
        self.assertEqual(payload["parseReport"]["nonFiniteLines"], 1)
        self.assertIn("non-finite", payload["parseWarning"])


if __name__ == "__main__":
    unittest.main()
