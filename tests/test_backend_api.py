from pathlib import Path
import math
import os
import sys
import unittest
import tempfile


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))

os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "rmc_toolkits_cache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import app as backend_app


DATA = ROOT / "data"


class BackendApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

    def tearDown(self):
        generated = ROOT / "results" / "Frac_coord_api_test.txt"
        if generated.exists():
            generated.unlink()

    def test_health_reports_data_root(self):
        response = self.client.get("/api/health")

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["status"], "ok")
        self.assertEqual(payload["dataRoot"], str(ROOT))

    def test_files_lists_supported_outputs_with_plot_kind(self):
        response = self.client.get("/api/files", query_string={"dir": "data"})

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        files = {item["name"]: item for item in payload["files"]}
        self.assertEqual(payload["currentPath"], str(DATA))
        self.assertEqual(files["GNSe_FQ1.csv"]["plotKind"], "xray_sq")
        self.assertEqual(files["GNSe-02.log"]["plotKind"], "r_value")
        self.assertIn("GNSe.rmc6f", files)

    def test_find_rmc6f_prefers_output_stem_match(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            earlier = directory / "alpha.rmc6f"
            expected = directory / "beta.rmc6f"
            earlier.write_text("", encoding="utf-8")
            expected.write_text("", encoding="utf-8")
            (directory / "beta-01.log").write_text("", encoding="utf-8")

            self.assertEqual(backend_app._find_rmc6f(directory), expected)

    def test_rejects_paths_outside_data_root(self):
        response = self.client.get("/api/files", query_string={"dir": "/private/tmp"})

        self.assertEqual(response.status_code, 403)
        self.assertIn("outside configured data root", response.get_json()["error"])

    def test_folder_dialog_registers_selected_external_folder(self):
        original_choose_folder = backend_app._choose_folder
        with tempfile.TemporaryDirectory() as tmpdir:
            selected = Path(tmpdir).resolve()
            backend_app._choose_folder = lambda _initial_dir: selected
            try:
                response = self.client.post("/api/dialog/folder", json={"dir": "data"})
                self.assertEqual(response.status_code, 200)
                self.assertEqual(response.get_json()["path"], str(selected))

                files_response = self.client.get("/api/files", query_string={"dir": str(selected)})
                self.assertEqual(files_response.status_code, 200)
                self.assertEqual(files_response.get_json()["currentPath"], str(selected))
            finally:
                backend_app._choose_folder = original_choose_folder
                backend_app.SELECTED_DATA_ROOTS.discard(selected)

    def test_plot_metadata_and_data_endpoints(self):
        metadata_response = self.client.get(
            "/api/plot/metadata",
            query_string={"path": "data/GNSe_FQ1.csv"},
        )
        data_response = self.client.get(
            "/api/plot/data",
            query_string={"path": "data/GNSe_FQ1.csv"},
        )

        self.assertEqual(metadata_response.status_code, 200)
        metadata = metadata_response.get_json()
        self.assertEqual(metadata["kind"], "xray_sq")
        self.assertEqual(metadata["title"], "S(Q) (x-ray)")
        self.assertGreater(metadata["metrics"]["rwp"], 0.0)

        self.assertEqual(data_response.status_code, 200)
        payload = data_response.get_json()
        self.assertEqual(payload["xLabel"], "Q (Å^{-1})")
        self.assertEqual(len(payload["series"]), 2)
        self.assertEqual(len(payload["series"][0]["x"]), 2649)

    def test_r_value_data_endpoint_combines_related_logs_in_numeric_order(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmpdir:
            directory = Path(tmpdir)
            (directory / "run-10.log").write_text("header\nheader\n1 0.1 10.0\n", encoding="utf-8")
            (directory / "run-01.log").write_text("header\nheader\n1 0.1 1.0\n", encoding="utf-8")
            (directory / "run-02.log").write_text("header\nheader\n1 0.1 2.0\n", encoding="utf-8")

            response = self.client.get(
                "/api/plot/data",
                query_string={"path": str(directory / "run-01.log")},
            )

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["series"][0]["x"], [0, 1, 2])
        self.assertEqual(payload["series"][0]["y"], [0.0, math.log(2.0), math.log(10.0)])
        self.assertAlmostEqual(payload["metrics"]["final_chi_r"], 10.0)

    def test_convert_frac_writes_requested_output_inside_root(self):
        response = self.client.post(
            "/api/convert/frac",
            json={
                "path": "data/GNSe.rmc6f",
                "outputPath": "results/Frac_coord_api_test.txt",
            },
        )

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["name"], "Frac_coord_api_test.txt")
        self.assertTrue((ROOT / "results" / "Frac_coord_api_test.txt").exists())

    def test_structure_endpoint_samples_atoms_and_reports_metadata(self):
        response = self.client.get(
            "/api/structure",
            query_string={"dir": "data", "maxPoints": 500},
        )

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["totalAtoms"], 52000)
        self.assertLessEqual(payload["sampledAtoms"], 500)
        self.assertEqual(payload["elements"], ["Ga", "Nb", "Se"])
        self.assertEqual(payload["elementCounts"]["Ga"], 4000)
        self.assertEqual(len(payload["points"]), payload["sampledAtoms"])
        self.assertIn("latticeVectors", payload)

    def test_kde_slice_endpoint_returns_density_payload(self):
        response = self.client.get(
            "/api/kde/slice",
            query_string={
                "dir": "data",
                "element": "Ga",
                "z": 0.775,
                "dz": 0.08,
                "grid": 16,
                "levels": 2,
            },
        )

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["element"], "Ga")
        self.assertEqual(payload["grid"], 16)
        self.assertEqual(payload["cellLengths"], [10.4116, 10.4116, 10.4116])
        self.assertGreater(payload["slabCount"], 0)
        self.assertEqual(len(payload["density"]), 16)

    def test_kde_slice_endpoint_handles_empty_default_slab(self):
        response = self.client.get(
            "/api/kde/slice",
            query_string={
                "dir": "data",
                "orientation": "c",
                "z": 0.5,
                "dz": 0.08,
                "grid": 16,
                "levels": 1,
            },
        )

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["orientation"], "c")
        self.assertEqual(payload["slabCount"], 0)
        self.assertEqual(payload["fitCount"], 0)
        self.assertEqual(len(payload["density"]), 16)

    def test_kde_slice_endpoint_accepts_custom_orientation(self):
        response = self.client.get(
            "/api/kde/slice",
            query_string={
                "dir": "data",
                "orientation": "custom",
                "nx": 1,
                "ny": 1,
                "nz": 0,
                "z": 0.5,
                "dz": 0.12,
                "grid": 16,
                "levels": 1,
            },
        )

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["orientation"], "custom")
        self.assertEqual(payload["grid"], 16)
        self.assertEqual(len(payload["normal"]), 3)
        self.assertTrue(payload["planeVertices"])
        self.assertTrue(payload["slabVertices"])


if __name__ == "__main__":
    unittest.main()
