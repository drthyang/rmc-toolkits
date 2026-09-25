# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Which configuration files a run folder is read from.

* ``read_structure()`` must pair ``Frac_coord_<stem>.txt`` with ``<stem>.rmc6f``:
  it used to take ``sorted(glob('Frac*.txt'))[0]`` and ``sorted(glob('*.rmc6f'))[0]``
  independently, so the supercell used for the fold and the element → reference
  map could come from a different configuration (data/250K_try1/supercell: half
  the atoms per element, wrong x fold).
* A 0-byte or marker-less ``.rmc6f`` (left by a killed run) must not hide valid
  configurations in the same folder (Flask ``_find_rmc6f``; the browser chooser
  is pinned in browserData.test.js).
"""

from __future__ import annotations

import os
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import read_structure, rmc6f_problem, write_frac_from_rmc6f

ROOT = Path(__file__).resolve().parents[1]


def write_rmc6f(path: Path, supercell, sites, edge: float = 5.0) -> Path:
    """A minimal valid configuration: ``sites`` = [(element, (u, v, w))], one copy per cell."""
    lines = []
    atom = 0
    for cx in range(supercell[0]):
        for cy in range(supercell[1]):
            for cz in range(supercell[2]):
                for ref, (element, frac) in enumerate(sites, start=1):
                    atom += 1
                    box = [(frac[i] + cell) / supercell[i] for i, cell in enumerate((cx, cy, cz))]
                    lines.append(
                        f"{atom:6d}  {element} [1]  {box[0]:.6f}  {box[1]:.6f}  {box[2]:.6f}"
                        f"  {ref}  {cx}  {cy}  {cz}"
                    )
    header = [
        "(Version 6f format configuration file)",
        f"Number of atoms:   {atom}",
        f"Supercell dimensions:   {supercell[0]}  {supercell[1]}  {supercell[2]}",
        "Lattice vectors (Ang):",
        f"  {edge * supercell[0]:.6f}  0.000000  0.000000",
        f"  0.000000  {edge * supercell[1]:.6f}  0.000000",
        f"  0.000000  0.000000  {edge * supercell[2]:.6f}",
        "Atoms:",
    ]
    path.write_text("\n".join(header + lines) + "\n", encoding="utf-8")
    return path


ZETA_SITES = [("Ga", (0.1, 0.2, 0.3)), ("Se", (0.6, 0.6, 0.6))]
ALPHA_SITES = [("Ga", (0.4, 0.4, 0.4))]


class ReadStructurePairingTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.directory = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _zeta(self) -> Path:
        path = write_rmc6f(self.directory / "zeta.rmc6f", (2, 1, 1), ZETA_SITES)
        write_frac_from_rmc6f(path)
        return path

    def test_frac_file_is_paired_with_its_own_configuration(self):
        # alpha.rmc6f sorts first but is a different configuration (3x1x1, Ga only).
        write_rmc6f(self.directory / "alpha.rmc6f", (3, 1, 1), ALPHA_SITES)
        self._zeta()

        selenium = read_structure(self.directory, element="Se", mode="fractional")

        np.testing.assert_array_equal(selenium.supercell, [2, 1, 1])
        self.assertEqual(len(selenium.positions), 2)          # one Se per cell, 2 cells
        np.testing.assert_allclose(selenium.positions, [[0.6, 0.6, 0.6]] * 2, atol=1e-4)
        gallium = read_structure(self.directory, element="Ga", mode="fractional")
        np.testing.assert_allclose(gallium.positions, [[0.1, 0.2, 0.3]] * 2, atol=1e-4)

    def test_empty_stem_match_is_skipped_for_the_next_pair(self):
        (self.directory / "aaa.rmc6f").write_text("", encoding="utf-8")
        (self.directory / "Frac_coord_aaa.txt").write_text("h\nh\nh\nh\nh\n 1 0.1 0.1 0.1 0 0 0\n", encoding="utf-8")
        self._zeta()

        structure = read_structure(self.directory, element="Se", mode="fractional")

        np.testing.assert_array_equal(structure.supercell, [2, 1, 1])
        self.assertEqual(len(structure.positions), 2)

    def test_ambiguous_folder_raises_instead_of_mixing_runs(self):
        write_rmc6f(self.directory / "alpha.rmc6f", (3, 1, 1), ALPHA_SITES)
        write_rmc6f(self.directory / "beta.rmc6f", (2, 1, 1), ZETA_SITES)
        (self.directory / "Frac_coord.txt").write_text("h\nh\nh\nh\nh\n 1 0.1 0.1 0.1 0 0 0\n", encoding="utf-8")

        with self.assertRaisesRegex(ValueError, "Cannot tell which configuration"):
            read_structure(self.directory)

    def test_explicit_paths_and_the_consistency_check(self):
        zeta = self._zeta()
        alpha = write_rmc6f(self.directory / "alpha.rmc6f", (3, 1, 1), ALPHA_SITES)
        frac = self.directory / "Frac_coord_zeta.txt"

        structure = read_structure(self.directory, element="Ga", frac_path=frac, rmc6f_path=zeta)
        self.assertEqual(len(structure.positions), 2)
        # zeta's Frac file carries reference number 2, which alpha does not have.
        with self.assertRaisesRegex(ValueError, "does not belong to alpha.rmc6f"):
            read_structure(self.directory, frac_path=frac, rmc6f_path=alpha)

    def test_cell_indices_beyond_the_supercell_are_rejected(self):
        zeta = self._zeta()
        small = write_rmc6f(self.directory / "small.rmc6f", (1, 1, 1), ZETA_SITES)
        with self.assertRaisesRegex(ValueError, "exceed the supercell"):
            read_structure(self.directory, frac_path=self.directory / "Frac_coord_zeta.txt", rmc6f_path=small)
        self.assertIsNotNone(zeta)


class Rmc6fProblemTests(unittest.TestCase):
    def test_reports_empty_and_marker_less_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            empty = directory / "empty.rmc6f"
            empty.write_text("", encoding="utf-8")
            header_only = directory / "header.rmc6f"
            header_only.write_text("(Version 6f format configuration file)\nNumber of atoms: 4\n", encoding="utf-8")
            good = write_rmc6f(directory / "good.rmc6f", (1, 1, 1), ZETA_SITES)

            self.assertEqual(rmc6f_problem(empty), "empty (0 bytes)")
            self.assertIn("no Atoms section", rmc6f_problem(header_only))
            self.assertIsNone(rmc6f_problem(good))


if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402


class FlaskStructureChoiceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

    def test_empty_stem_matched_rmc6f_does_not_hide_valid_configurations(self):
        # data/250K_try1/supercell: Frac_coord_new_x.txt stem-matches new_x.rmc6f (0 bytes).
        with tempfile.TemporaryDirectory(dir=ROOT) as tmpdir:
            directory = Path(tmpdir)
            (directory / "new_x.rmc6f").write_text("", encoding="utf-8")
            (directory / "Frac_coord_new_x.txt").write_text("h\n", encoding="utf-8")
            (directory / "Frac_coord_new_y.txt").write_text("h\n", encoding="utf-8")
            write_rmc6f(directory / "new_y.rmc6f", (2, 1, 1), ZETA_SITES)

            response = self.client.get("/api/structure", query_string={"dir": str(directory)})

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertTrue(payload["source"].endswith("new_y.rmc6f"))
        self.assertEqual(payload["totalAtoms"], 4)

    def test_no_usable_configuration_names_each_candidate(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmpdir:
            directory = Path(tmpdir)
            (directory / "run.rmc6f").write_text("", encoding="utf-8")

            response = self.client.get("/api/structure", query_string={"dir": str(directory)})

        self.assertEqual(response.status_code, 404)
        self.assertIn("run.rmc6f (empty (0 bytes))", response.get_json()["error"])


if __name__ == "__main__":
    unittest.main()
