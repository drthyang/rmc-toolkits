# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Flask edge cases left after the 1.0 gate.

Each is a clear 4xx with a JSON error naming the problem -- never a 200 with
an empty result, a half-written file, an HTML page or a 500. One class per
case; the class docstrings say what used to happen.
"""

from pathlib import Path
import json
import os
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402
from test_backend_validation import strict_json, write_synthetic_rmc6f  # noqa: E402

RUN = "results/backend_edges"


class _EdgeCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        cls.run_dir = ROOT / RUN
        shutil.rmtree(cls.run_dir, ignore_errors=True)
        cls.run_dir.mkdir(parents=True)
        write_synthetic_rmc6f(cls.run_dir / "synthetic.rmc6f", supercell=(3, 3, 3))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.run_dir, ignore_errors=True)

    def status_error(self, response):
        body = response.get_data(as_text=True)
        self.assertTrue(response.is_json, f"not JSON: {response.content_type} {body[:120]}")
        return response.status_code, strict_json(body).get("error", "")

    def folder_with(self, name: str, text: str) -> str:
        folder = self.run_dir / name
        folder.mkdir(exist_ok=True)
        (folder / "run.rmc6f").write_text(text, encoding="utf-8")
        return f"{RUN}/{name}"


EMPTY_ATOMS = (
    "Supercell dimensions: 1 1 1\nLattice vectors (Ang):\n8 0 0\n0 8 0\n0 0 8\nAtoms:\n"
)


class KdeSliceElementTests(_EdgeCase):
    """/api/kde/slice: an unknown element (or no parseable atom) was a 200 all-zero map
    reading "No atoms in this slab."; element case is now ignored as on /api/pca/*."""

    def slice(self, directory=RUN, **params):
        return self.client.get(
            "/api/kde/slice", query_string={"dir": directory, "grid": 32, "levels": 0, **params}
        )

    def test_unknown_element_is_a_400_naming_it(self):
        status, error = self.status_error(self.slice(element="Xx"))
        self.assertEqual(status, 400, error)
        self.assertEqual(error, "Unknown element 'Xx'; available: Nb, Se")

    def test_element_case_is_ignored(self):
        upper = self.slice(element="Se", z=0.5)
        lower = self.slice(element="se", z=0.5)
        self.assertEqual(upper.status_code, 200)
        self.assertEqual(lower.status_code, 200, lower.get_data(as_text=True)[:200])
        self.assertEqual(lower.get_json()["element"], "Se")
        self.assertEqual(lower.get_json()["slabCount"], upper.get_json()["slabCount"])
        self.assertGreater(upper.get_json()["slabCount"], 0)

    def test_a_file_without_parseable_atoms_is_a_400(self):
        folder = self.folder_with("no_atoms", EMPTY_ATOMS)
        for element in (None, "Se"):
            with self.subTest(element=element):
                params = {} if element is None else {"element": element}
                status, error = self.status_error(self.slice(folder, **params))
                self.assertEqual(status, 400, error)
                self.assertIn("no atoms could be parsed", error)


if __name__ == "__main__":
    unittest.main()
