# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""/api/clusters and /api/modulations: the cluster-distortion and modulation routes."""

import json
import shutil
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from cluster_configs import cluster_config_text  # noqa: E402
from rmc_toolkits.clusters import cluster_distortion  # noqa: E402
from rmc_toolkits.modulations import displacement_modulations  # noqa: E402
from rmc_toolkits.pca_kde import load_site_displacements  # noqa: E402

try:  # Flask is a backend-only dependency; skip cleanly without it.
    import web_app.backend.app as backend_app
except ImportError:  # pragma: no cover - exercised only without flask
    backend_app = None


@unittest.skipIf(backend_app is None, "Flask backend not importable")
class ClustersApiTests(unittest.TestCase):
    RUN = "results/clusters_api_test"

    @classmethod
    def setUpClass(cls):
        cls.run_dir = ROOT / cls.RUN
        cls.run_dir.mkdir(parents=True, exist_ok=True)
        (cls.run_dir / "trigonal.rmc6f").write_text(cluster_config_text("trigonal"), encoding="utf-8")
        cls.client = backend_app.app.test_client()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.run_dir, ignore_errors=True)

    def get(self, **params):
        return self.client.get("/api/clusters", query_string={"dir": self.RUN, **params})

    def test_the_payload(self):
        response = self.get(element="nb", rMin="2.5", rMax="3.5")
        self.assertEqual(response.status_code, 200, response.get_json())
        payload = response.get_json()
        self.assertEqual(payload["element"], "Nb")
        self.assertEqual(payload["kind"], {"members": 4, "edges": 6, "clusters": 64, "templates": 1})
        self.assertIn("t2Trigonality", payload["z"])
        self.assertTrue(payload["source"].endswith("trigonal.rmc6f"))
        self.assertIsNone(payload["parseWarning"])

    def test_bad_requests(self):
        cases = [
            ({"element": "Nb", "rMin": "2.5"}, "rMax required"),
            ({"rMin": "2.5", "rMax": "3.5"}, "element required"),
            ({"element": "Nb", "rMin": "2.5", "rMax": "20"}, "rMax must be <= 8"),
            ({"element": "Nb", "rMin": "3.5", "rMax": "2.5"}, "window must be"),
            ({"element": "Xx", "rMin": "2.5", "rMax": "3.5"}, "no site holds Xx"),
        ]
        for params, message in cases:
            with self.subTest(params=params):
                response = self.get(**params)
                self.assertEqual(response.status_code, 400)
                self.assertIn(message, response.get_json()["error"])

    def test_an_extended_network_is_reported_not_an_error(self):
        payload = self.get(element="Se", rMin="7.5", rMax="8").get_json()
        self.assertTrue(payload["extended"])
        self.assertEqual(payload["templates"], [])
        self.assertEqual(payload["verdict"]["window"]["extendedNetwork"]["because"], "extended_network")
        self.assertIsNone(payload["verdict"]["distorted"]["answer"])

    def test_the_route_returns_the_engine_verdict(self):
        sites = load_site_displacements(self.run_dir / "trigonal.rmc6f")
        for element, window in (("Nb", (2.5, 3.5)), ("Ga", (2.0, 3.0)), ("Se", (7.5, 8.0))):
            with self.subTest(element=element):
                payload = self.get(element=element, rMin=str(window[0]), rMax=str(window[1])).get_json()
                engine = json.loads(json.dumps(cluster_distortion(sites, element, window)["verdict"], allow_nan=False))
                self.assertEqual(payload["verdict"], engine)


if __name__ == "__main__":
    unittest.main()


@unittest.skipIf(backend_app is None, "Flask backend not importable")
class ModulationsApiTests(unittest.TestCase):
    RUN = "results/modulations_api_test"

    @classmethod
    def setUpClass(cls):
        cls.run_dir = ROOT / cls.RUN
        cls.run_dir.mkdir(parents=True, exist_ok=True)
        (cls.run_dir / "doubled.rmc6f").write_text(cluster_config_text("doubled"), encoding="utf-8")
        cls.client = backend_app.app.test_client()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.run_dir, ignore_errors=True)

    def get(self, **params):
        return self.client.get("/api/modulations", query_string={"dir": self.RUN, **params})

    def test_every_element_and_one(self):
        payload = self.get().get_json()
        self.assertEqual([entry["element"] for entry in payload["elements"]], ["Ga", "Nb", "Se"])
        self.assertEqual(len(self.get(element="all").get_json()["elements"]), 3)
        nb = self.get(element="nb").get_json()
        self.assertEqual(nb["elements"][0]["strongest"][0]["q"], [0.0, 0.0, 0.5])
        self.assertTrue(nb["source"].endswith("doubled.rmc6f"))

    def test_an_unknown_element_is_a_400(self):
        response = self.get(element="Xx")
        self.assertEqual(response.status_code, 400)
        self.assertIn("no site holds Xx", response.get_json()["error"])

    def test_the_route_returns_the_engine_verdict(self):
        sites = load_site_displacements(self.run_dir / "doubled.rmc6f")
        for element in (None, "Nb"):
            with self.subTest(element=element):
                payload = self.get(**({"element": element} if element else {})).get_json()
                engine = json.loads(json.dumps(displacement_modulations(sites, element), allow_nan=False))
                for key in ("family", "rules"):
                    self.assertEqual(payload[key], engine[key])
                self.assertEqual([entry["verdict"] for entry in payload["elements"]],
                                 [entry["verdict"] for entry in engine["elements"]])
                self.assertEqual([entry["zoneBoundary"] for entry in payload["elements"]],
                                 [entry["zoneBoundary"] for entry in engine["elements"]])
        nb = self.get(element="Nb").get_json()["elements"][0]
        self.assertEqual(nb["verdict"]["cellDoubling"]["decidedBy"], [[0.0, 0.0, 0.5]])
