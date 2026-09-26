# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""``.rmc6f`` header validation: lattice and supercell, one rule in both runtimes.

``read_cell_vectors`` (and the browser's ``readRmc6fCellVectors``) used to
accept ``Supercell dimensions: 0 0 0`` (every atom folded onto the origin), a
NaN, collinear or overflowing lattice (200 responses with nonsense or all-zero
results), and rejected Fortran ``D`` exponents with a raw float error although
the atom lines accept them. The cases live in a fixture shared with
``src/__tests__/rmc6fHeader.test.js`` so the two runtimes give the same values
and the same error text.
"""

from pathlib import Path
import json
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

from rmc_toolkits.parsers import read_cell_vectors, write_frac_from_rmc6f  # noqa: E402

FIXTURE = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "rmc6f_header_cases.json"
CASES = json.loads(FIXTURE.read_text(encoding="utf-8"))["cases"]


class ReadCellVectorsHeaderTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def write(self, case) -> Path:
        path = self.tmp / f"{case['name']}.rmc6f"
        path.write_text(case["text"], encoding="utf-8")
        return path

    def test_shared_cases(self):
        self.assertGreaterEqual(len(CASES), 15)
        for case in CASES:
            with self.subTest(case["name"]):
                path = self.write(case)
                if "error" in case:
                    expected = case["error"].replace("{name}", str(path))
                    with self.assertRaises(ValueError) as caught:
                        read_cell_vectors(path)
                    self.assertEqual(str(caught.exception), expected)
                else:
                    lattice, supercell = read_cell_vectors(path)
                    np.testing.assert_allclose(lattice, case["latticeVectors"], rtol=0, atol=1e-12)
                    np.testing.assert_array_equal(supercell, case["supercell"])
                    self.assertEqual(lattice.shape, (3, 3))
                    self.assertEqual(supercell.dtype, float)

    def test_frac_conversion_refuses_a_zero_supercell(self):
        case = next(case for case in CASES if case["name"] == "zero_supercell")
        path = self.write(case)
        out = self.tmp / "Frac_coord_zero.txt"
        with self.assertRaisesRegex(ValueError, "supercell dimensions must be three positive integers"):
            write_frac_from_rmc6f(path, out)
        self.assertFalse(out.exists())


class HeaderValidationRouteTests(unittest.TestCase):
    """Every analysis route answers a bad header with a 400 naming the problem."""

    @classmethod
    def setUpClass(cls):
        import app as backend_app

        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        cls.root = ROOT / "results" / "rmc6f_header_routes"
        shutil.rmtree(cls.root, ignore_errors=True)
        cls.root.mkdir(parents=True)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.root, ignore_errors=True)

    def folder(self, name: str) -> Path:
        case = next(case for case in CASES if case["name"] == name)
        folder = self.root / name
        folder.mkdir(exist_ok=True)
        (folder / "run.rmc6f").write_text(case["text"], encoding="utf-8")
        return folder

    def routes(self, folder: Path):
        rel = str(folder.relative_to(ROOT))
        return {
            "structure": f"/api/structure?dir={rel}",
            "kde": f"/api/kde/slice?dir={rel}&element=Se",
            "pca sites": f"/api/pca/sites?dir={rel}",
            "pca kde": f"/api/pca/kde?dir={rel}&referenceNumber=1",
            "orientation": f"/api/pca/orientation?dir={rel}&referenceNumber=1",
            "triplets": f"/api/triplets?dir={rel}&end1=Se&apex=Nb&end2=Se&r12Min=1&r12Max=5",
        }

    def test_bad_headers_are_400_on_every_route(self):
        for name, fragment in (
            ("zero_supercell", "supercell dimensions must be three positive integers"),
            ("nan_lattice", "lattice vector 1 must be three finite numbers"),
            ("collinear_lattice", "lattice vectors are singular"),
            ("huge_lattice", "non-finite cell volume"),
        ):
            folder = self.folder(name)
            for route, url in self.routes(folder).items():
                with self.subTest(case=name, route=route):
                    response = self.client.get(url)
                    self.assertEqual(response.status_code, 400, response.get_data(as_text=True))
                    self.assertIn(fragment, response.get_json()["error"])

    def test_fortran_d_exponent_header_is_read(self):
        folder = self.folder("fortran_d_exponents")
        response = self.client.get(f"/api/structure?dir={folder.relative_to(ROOT)}")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        payload = response.get_json()
        np.testing.assert_allclose(payload["latticeVectors"], np.eye(3) * 10.0)
        self.assertEqual(payload["supercell"], [2.0, 2.0, 2.0])

    def test_frac_conversion_of_a_bad_header_is_400(self):
        folder = self.folder("zero_supercell")
        response = self.client.post(
            "/api/convert/frac",
            json={"path": str((folder / "run.rmc6f").relative_to(ROOT))},
        )
        self.assertEqual(response.status_code, 400, response.get_data(as_text=True))
        self.assertFalse(any(folder.glob("Frac_coord*")))


if __name__ == "__main__":
    unittest.main()
