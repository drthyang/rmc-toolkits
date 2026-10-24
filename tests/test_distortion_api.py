# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""/api/distortion: the polyhedral-distortion route (rmc_toolkits/distortion.py)."""

import json
import shutil
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

import distortion_configs as configs  # noqa: E402
from rmc_toolkits.distortion import APP_MAX_BONDS, polyhedral_distortion_from_file  # noqa: E402

try:  # Flask is a backend-only dependency; skip cleanly without it.
    import web_app.backend.app as backend_app
except ImportError:  # pragma: no cover - exercised only without flask
    backend_app = None


@unittest.skipIf(backend_app is None, "Flask backend not importable")
class DistortionApiTests(unittest.TestCase):
    RUN = "results/distortion_api_test"

    @classmethod
    def setUpClass(cls):
        cls.run_dir = ROOT / cls.RUN
        cls.run_dir.mkdir(parents=True, exist_ok=True)
        cls.path = cls.run_dir / "network.rmc6f"
        cls.path.write_text(configs.rmc6f_text(*configs.jahn_teller_network(sigma=0.02)), encoding="utf-8")
        cls.client = backend_app.app.test_client()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.run_dir, ignore_errors=True)

    def setUp(self):
        backend_app._DISTORTION_CACHE.clear()

    def get(self, **params):
        return self.client.get("/api/distortion", query_string={"dir": self.RUN, **params})

    def test_the_payload_is_the_engine_s(self):
        response = self.get(centre="cu", ligand="F", rMin="1.6", rMax="2.6")
        self.assertEqual(response.status_code, 200, response.get_json())
        payload = response.get_json()
        self.assertTrue(payload.pop("source").endswith("network.rmc6f"))
        engine = polyhedral_distortion_from_file(self.path, "Cu", "F", 1.6, 2.6, APP_MAX_BONDS)
        self.assertEqual(payload, json.loads(json.dumps(engine, allow_nan=False)))
        self.assertEqual(payload["polyhedron"]["shape"], "octahedron")
        self.assertEqual(payload["configuration"]["polyhedra"], 27)
        self.assertGreater(payload["configuration"]["distortionIndex"]["mean"], 0.06)
        self.assertLess(payload["average"]["distortionIndex"]["mean"], 0.01)
        self.assertIsNone(payload["parseWarning"])

    def test_bad_requests(self):
        cases = [
            ({"centre": "Cu", "ligand": "F", "rMin": "1.6"}, "rMax required"),
            ({"ligand": "F", "rMin": "1.6", "rMax": "2.6"}, "centre required"),
            ({"centre": "Cu", "rMin": "1.6", "rMax": "2.6"}, "ligand required"),
            ({"centre": "Cu", "ligand": "F", "rMin": "1.6", "rMax": "9"}, "rMax must be <= 8"),
            ({"centre": "Cu", "ligand": "F", "rMin": "-1", "rMax": "2.6"}, "rMin must be >= 0"),
            ({"centre": "Cu", "ligand": "F", "rMin": "2.6", "rMax": "1.6"}, "window must be"),
            ({"centre": "Xx", "ligand": "F", "rMin": "1.6", "rMax": "2.6"}, "unknown element 'Xx'"),
            ({"centre": "Cu", "ligand": "F", "rMin": "0.5", "rMax": "1.0"}, "no polyhedra"),
        ]
        for params, message in cases:
            with self.subTest(params=params):
                response = self.get(**params)
                self.assertEqual(response.status_code, 400)
                self.assertIn(message, response.get_json()["error"])

    def test_the_bond_budget_is_the_engine_s_and_in_the_cache_key(self):
        self.assertEqual(backend_app.DISTORTION_MAX_BONDS, APP_MAX_BONDS)
        original = backend_app.DISTORTION_MAX_BONDS
        try:
            backend_app.DISTORTION_MAX_BONDS = 100
            response = self.get(centre="Cu", ligand="F", rMin="1.6", rMax="2.6")
            self.assertEqual(response.status_code, 400)
            self.assertIn("holds 162 bonds, over the limit of 100", response.get_json()["error"])
        finally:
            backend_app.DISTORTION_MAX_BONDS = original
        self.assertEqual(self.get(centre="Cu", ligand="F", rMin="1.6", rMax="2.6").status_code, 200)

    def test_a_rewritten_file_is_read_again(self):
        first = self.get(centre="Cu", ligand="F", rMin="1.6", rMax="2.6").get_json()
        text = self.path.read_text(encoding="utf-8")
        try:
            # The file signature (mtime, ctime, size, inode) keys the cache.
            self.path.write_text(configs.rmc6f_text(*configs.jahn_teller_network(sigma=0.04, seed=8)), encoding="utf-8")
            second = self.get(centre="Cu", ligand="F", rMin="1.6", rMax="2.6").get_json()
            self.assertNotEqual(first["configuration"]["distortionIndex"], second["configuration"]["distortionIndex"])
        finally:
            self.path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    unittest.main()
