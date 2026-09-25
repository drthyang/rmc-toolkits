# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Documented contracts of the Flask API (docs/REFERENCE.md, notation.md §3b).

* ``/api/kde/slice``: ``z``/``dz`` are fractions of the unit cube's projection
  range along the slice normal and the whole KDE path stays fractional -- no
  conversion to Angstrom, so scaling the lattice changes nothing.
* ``/api/structure`` returns the ``.rmc6f`` move counters in Flask mode too.
* Every ``/api/*`` route is listed in the REFERENCE.md endpoint table.
"""

from pathlib import Path
import os
import re
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))

os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "rmc_toolkits_cache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import numpy as np  # noqa: E402

import app as backend_app  # noqa: E402


def rmc6f_text(cell: float, *, moves: bool = False) -> str:
    rng = np.random.default_rng(11)
    n = 4
    lines = []
    if moves:
        lines += [
            "Number of moves generated:           2689753",
            "Number of moves tried:               2685185",
            "Number of moves accepted:            586675",
        ]
    lines += [
        f"Supercell dimensions {n} {n} {n}",
        "Lattice vectors (Ang):",
        f"{cell * n} 0.0 0.0",
        f"0.0 {cell * n} 0.0",
        f"0.0 0.0 {cell * n}",
        "Atoms:",
    ]
    atom = 0
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                for reference, element, basis in ((1, "Nb", 0.0), (2, "Se", 0.5)):
                    atom += 1
                    coord = ((np.array([ix, iy, iz]) + basis) / n + rng.normal(size=3) * 0.004) % 1.0
                    lines.append(
                        f"{atom} {element} [{reference}] {coord[0]:.10f} {coord[1]:.10f} "
                        f"{coord[2]:.10f} {reference} {ix} {iy} {iz}"
                    )
    return "\n".join(lines) + "\n"


class ApiContractTests(unittest.TestCase):
    RUNS = ("results/contract_cell8", "results/contract_cell16")

    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        for run, cell in zip(cls.RUNS, (8.0, 16.0)):
            directory = ROOT / run
            directory.mkdir(parents=True, exist_ok=True)
            (directory / "run.rmc6f").write_text(rmc6f_text(cell, moves=True), encoding="utf-8")

    @classmethod
    def tearDownClass(cls):
        for run in cls.RUNS:
            shutil.rmtree(ROOT / run, ignore_errors=True)

    def slice(self, run, **params):
        response = self.client.get(
            "/api/kde/slice",
            query_string={"dir": run, "element": "Se", "z": 0.5, "dz": 0.1, "grid": 24, **params},
        )
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True)[:300])
        return response.get_json()

    def test_kde_slice_is_fractional_throughout(self):
        for params in ({"orientation": "c"}, {"orientation": "custom", "nx": 1, "ny": 1, "nz": 1}):
            with self.subTest(**params):
                small, large = (self.slice(run, **params) for run in self.RUNS)
                self.assertEqual(large["cellLengths"], [16.0, 16.0, 16.0])
                self.assertEqual(small["slabCount"], large["slabCount"])
                self.assertGreater(small["slabCount"], 0)
                np.testing.assert_allclose(small["density"], large["density"], rtol=1e-12)
                self.assertEqual(small["depthThickness"], large["depthThickness"])

    def test_dz_is_a_fraction_of_the_projection_range(self):
        cube = self.slice(self.RUNS[0], orientation="custom", nx=1, ny=1, nz=1)
        np.testing.assert_allclose(cube["depthRange"], [0.0, np.sqrt(3.0)], atol=1e-12)
        self.assertAlmostEqual(cube["depthThickness"], 0.1 * np.sqrt(3.0), places=12)
        preset = self.slice(self.RUNS[0], orientation="c")
        self.assertEqual(preset["depthRange"], [0.0, 1.0])
        self.assertAlmostEqual(preset["depthThickness"], 0.1, places=12)

    def test_structure_reports_move_counters(self):
        response = self.client.get("/api/structure", query_string={"dir": self.RUNS[0]})
        self.assertEqual(response.status_code, 200)
        moves = response.get_json()["moves"]
        self.assertEqual(moves["generated"], 2689753)
        self.assertEqual(moves["tried"], 2685185)
        self.assertEqual(moves["accepted"], 586675)

    def test_reference_lists_every_api_route(self):
        routes = sorted(
            {rule.rule for rule in backend_app.app.url_map.iter_rules() if rule.rule.startswith("/api/")}
        )
        self.assertEqual(len(routes), 15)
        reference = (ROOT / "docs" / "REFERENCE.md").read_text(encoding="utf-8")
        table = set(re.findall(r"^\| `(?:GET|POST) (/api/[^`?]+)`", reference, flags=re.MULTILINE))
        self.assertEqual(sorted(table), routes)


if __name__ == "__main__":
    unittest.main()
