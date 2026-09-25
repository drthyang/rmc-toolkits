# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""/api/triplets behaviour owned by the triplets engine (budget, windows).

Kept apart from test_backend_api.py so the triplets group's route tests do
not collide with other groups' edits there.
"""

import shutil
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

try:  # Flask is a backend-only dependency; skip cleanly without it.
    import web_app.backend.app as backend_app
except ImportError:  # pragma: no cover - exercised only without flask
    backend_app = None


def write_cloud_rmc6f(path: Path, count: int = 200, box: float = 10.0, seed: int = 3) -> None:
    rng = np.random.default_rng(seed)
    positions = rng.uniform(size=(count, 3))
    lines = [
        "Supercell dimensions: 1 1 1",
        "Lattice vectors (Ang):",
        f"{box:.6f} 0.000000 0.000000",
        f"0.000000 {box:.6f} 0.000000",
        f"0.000000 0.000000 {box:.6f}",
        "Atoms:",
    ]
    for index, p in enumerate(positions):
        element = "Nb" if index % 4 == 0 else "Se"
        lines.append(
            f"{index + 1} {element} [1] {p[0]:.12f} {p[1]:.12f} {p[2]:.12f} {index + 1} 0 0 0"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


@unittest.skipIf(backend_app is None, "Flask backend not importable")
class TripletsBudgetApiTests(unittest.TestCase):
    RUN = "results/triplets_budget_api_test"

    @classmethod
    def setUpClass(cls):
        cls.run_dir = ROOT / cls.RUN
        cls.run_dir.mkdir(parents=True, exist_ok=True)
        write_cloud_rmc6f(cls.run_dir / "cloud.rmc6f")
        cls.client = backend_app.app.test_client()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.run_dir, ignore_errors=True)

    def query(self, **extra):
        return {
            "dir": self.RUN,
            "end1": "Se",
            "apex": "Nb",
            "end2": "Se",
            "r12Min": 1.0,
            "r12Max": 4.0,
            **extra,
        }

    def test_over_budget_request_is_a_clear_400(self):
        allowed = self.client.get("/api/triplets", query_string=self.query())
        self.assertEqual(allowed.status_code, 200)
        exact = allowed.get_json()["angleCount"]
        self.assertGreater(exact, 10)
        with mock.patch.object(backend_app, "TRIPLETS_MAX_ANGLES", exact - 1):
            refused = self.client.get("/api/triplets", query_string=self.query(r12Min=1.0001))
        self.assertEqual(refused.status_code, 400)
        self.assertIn("angles", refused.get_json()["error"])
        self.assertIn("narrow", refused.get_json()["error"])

    def test_blank_bounds_are_missing_not_zero(self):
        # The same request shape the worker rejects: a cleared rmin box.
        for blank in ("", "   "):
            response = self.client.get("/api/triplets", query_string=self.query(r12Min=blank))
            self.assertEqual(response.status_code, 400, repr(blank))
            self.assertIn("required together; missing r12Min", response.get_json()["error"])
            half = self.client.get(
                "/api/triplets", query_string=self.query(r23Min=blank, r23Max=4.0)
            )
            self.assertEqual(half.status_code, 400, repr(blank))
            self.assertIn("missing r23Min", half.get_json()["error"])

    def test_route_uses_the_engine_budget_constant(self):
        from rmc_toolkits.triplets import APP_MAX_ANGLES

        self.assertEqual(backend_app.TRIPLETS_MAX_ANGLES, APP_MAX_ANGLES)


@unittest.skipIf(backend_app is None, "Flask backend not importable")
class RunConfigurationParityTests(unittest.TestCase):
    """rmc-triplets <run folder> analyses the configuration the app analyses."""

    LAYOUTS = [
        ["GaNb4Se8.rmc6f", "GaNb4Se8_5K.rmc6f", "GaNb4Se8_5K-00.log", "GaNb4Se8_5K_PDFpartials.csv"],
        ["start.rmc6f", "run.rmc6f", "run_FQ1.csv"],
        ["zeta.rmc6f", "alpha.rmc6f", "Frac_coord_zeta.txt"],
        ["one.rmc6f", "two.rmc6f", "two_bragg.csv", "one-01.log"],
        ["x.rmc6f", "y.rmc6f"],
    ]

    def test_cli_and_backend_choose_the_same_file(self):
        import tempfile

        from rmc_toolkits.triplets_cli import resolve_config

        for layout in self.LAYOUTS:
            with tempfile.TemporaryDirectory() as scratch:
                run = Path(scratch)
                for name in layout:
                    (run / name).write_text("", encoding="utf-8")
                self.assertEqual(
                    resolve_config(run), backend_app._find_rmc6f(run), layout
                )


if __name__ == "__main__":
    unittest.main()
