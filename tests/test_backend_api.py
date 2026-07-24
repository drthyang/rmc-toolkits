# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

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

# The GNSe example dataset is gitignored (see README), so sample-backed tests
# skip rather than fail when it is not present locally / in CI.
requires_sample = unittest.skipUnless(
    (DATA / "GNSe.rmc6f").exists(),
    "GNSe sample data not present in data/ (gitignored)",
)


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

    @requires_sample
    def test_files_lists_supported_outputs_with_plot_kind(self):
        response = self.client.get("/api/files", query_string={"dir": "data"})

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        files = {item["name"]: item for item in payload["files"]}
        self.assertEqual(payload["currentPath"], str(DATA))
        self.assertEqual(files["GNSe_FQ1.csv"]["plotKind"], "xray_sq")
        self.assertGreater(files["GNSe_FQ1.csv"]["modified"], 0)
        self.assertGreater(files["GNSe_FQ1.csv"]["size"], 0)
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

    def test_files_lists_exafs_outputs_with_plot_kind(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmpdir:
            directory = Path(tmpdir)
            (directory / "Nb-EXAFS-1_Q_OUTPUT.csv").write_text(
                " EXAFS #1,   chi(k)*k^2\n"
                "      k    ,  calculated  ,  experiment\n"
                " 3.300 ,    -0.32999 ,    -0.23780\n",
                encoding="utf-8",
            )
            (directory / "Nb-EXAFS-1_R_OUTPUT.csv").write_text(
                "     r   ,   Re_Calc  ,  Im_Calc  ,  Mod_Calc  ,   Re_Ex   ,   Im_Ex  ,   Mod_Ex\n"
                "   0.25000 ,    0.04898 ,   -0.29264 ,    0.29671 ,    0.04412 ,   -0.27550 ,    0.27901\n",
                encoding="utf-8",
            )

            response = self.client.get("/api/files", query_string={"dir": str(directory)})

        self.assertEqual(response.status_code, 200)
        files = {item["name"]: item for item in response.get_json()["files"]}
        self.assertEqual(files["Nb-EXAFS-1_Q_OUTPUT.csv"]["plotKind"], "exafs_q")
        self.assertEqual(files["Nb-EXAFS-1_R_OUTPUT.csv"]["plotKind"], "exafs_r")

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

    @requires_sample
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

    def test_exafs_q_data_endpoint_uses_exafs_labels(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmpdir:
            path = Path(tmpdir) / "Nb-EXAFS-1_Q_OUTPUT.csv"
            path.write_text(
                " EXAFS #1,   chi(k)*k^2\n"
                "      k    ,  calculated  ,  experiment\n"
                " 3.300 ,    -0.32999 ,    -0.23780\n"
                " 3.350 ,    -0.62595 ,    -0.33942\n",
                encoding="utf-8",
            )

            response = self.client.get("/api/plot/data", query_string={"path": str(path)})

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["kind"], "exafs_q")
        self.assertEqual(payload["title"], "EXAFS Q-space")
        self.assertEqual(payload["xLabel"], "k (Å^{-1})")
        self.assertEqual(payload["yLabel"], "χ(k) k²")
        self.assertEqual([series["label"] for series in payload["series"]], ["calculated", "experiment"])
        self.assertEqual(payload["series"][0]["x"], [3.3, 3.35])

    def test_exafs_r_data_endpoint_uses_exafs_labels(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmpdir:
            path = Path(tmpdir) / "Nb-EXAFS-1_R_OUTPUT.csv"
            path.write_text(
                "     r   ,   Re_Calc  ,  Im_Calc  ,  Mod_Calc  ,   Re_Ex   ,   Im_Ex  ,   Mod_Ex\n"
                "   0.25000 ,    0.04898 ,   -0.29264 ,    0.29671 ,    0.04412 ,   -0.27550 ,    0.27901\n",
                encoding="utf-8",
            )

            response = self.client.get("/api/plot/data", query_string={"path": str(path)})

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["kind"], "exafs_r")
        self.assertEqual(payload["title"], "EXAFS R-space")
        self.assertEqual(payload["xLabel"], "r (Å)")
        self.assertEqual(payload["yLabel"], "FT[χ(k) k²]")
        self.assertEqual([series["label"] for series in payload["series"]], ["Re_Calc", "Im_Calc", "Mod_Calc", "Re_Ex", "Im_Ex", "Mod_Ex"])

    @requires_sample
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

    @requires_sample
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

    @requires_sample
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

    @requires_sample
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

    @requires_sample
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


class OrientationApiTests(unittest.TestCase):
    """Displacement-direction histogram endpoint, on a synthetic one-site run.

    The run dir lives under results/ so API paths stay inside the data root.
    """

    @classmethod
    def setUpClass(cls):
        import numpy as np

        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

        cls.run_dir = ROOT / "results" / "orientation_api_test"
        cls.run_dir.mkdir(parents=True, exist_ok=True)
        rng = np.random.default_rng(2)
        supercell = (6, 6, 6)
        sigma = np.array([0.03, 0.01, 0.01])  # lobe along x
        lines = [
            "Supercell dimensions 6 6 6",
            "Lattice vectors (Ang):",
            "48.0 0.0 0.0",
            "0.0 48.0 0.0",
            "0.0 0.0 48.0",
            "Atoms:",
        ]
        atom = 0
        for ix in range(supercell[0]):
            for iy in range(supercell[1]):
                for iz in range(supercell[2]):
                    atom += 1
                    coord = np.array([ix, iy, iz]) / supercell + rng.normal(size=3) * sigma
                    lines.append(
                        f"{atom} Nb [1] {coord[0]:.10f} {coord[1]:.10f} {coord[2]:.10f} "
                        f"1 {ix} {iy} {iz}"
                    )
        (cls.run_dir / "synthetic.rmc6f").write_text("\n".join(lines), encoding="utf-8")

    @classmethod
    def tearDownClass(cls):
        import shutil

        shutil.rmtree(cls.run_dir, ignore_errors=True)

    def test_orientation_endpoint_returns_histogram(self):
        response = self.client.get(
            "/api/pca/orientation",
            query_string={
                "dir": "results/orientation_api_test",
                "referenceNumber": 1,
                "frequency": 4,
                "geometry": "false",
            },
        )

        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["referenceNumber"], 1)
        self.assertEqual(payload["element"], "Nb")
        self.assertEqual(payload["cellCount"], 162)
        self.assertEqual(payload["totalPoints"], 216)
        self.assertEqual(len(payload["enhancement"]), 162)
        self.assertNotIn("polygons", payload)
        # The cloud was drawn 3x wider along x: the lobe must sit there.
        self.assertGreater(abs(payload["peakDirection"][0]), 0.8)

    def test_orientation_endpoint_rejects_bad_weight(self):
        response = self.client.get(
            "/api/pca/orientation",
            query_string={
                "dir": "results/orientation_api_test",
                "referenceNumber": 1,
                "weight": "nope",
            },
        )
        self.assertEqual(response.status_code, 400)


class ScalingApiTests(unittest.TestCase):
    """Auto StoG endpoints, driven by a synthetic classic-stog run.

    Synthetic model kept in sync with tests/test_scaling_cli.py. The run dir
    lives under results/ so API paths stay inside the configured data root.
    """

    RHO0 = 0.05
    B2 = 0.02

    @classmethod
    def setUpClass(cls):
        import numpy as np
        from rmc_toolkits.parsers import write_stog_xy
        from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

        cls.run_dir = ROOT / "results" / "scaling_api_test"
        cls.run_dir.mkdir(parents=True, exist_ok=True)
        q = np.arange(20, 981) * 0.03
        r = np.arange(1, 12001) * 0.005
        onset = 0.5 * (1.0 + np.tanh((r - 2.65) / 0.07))
        peak = 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
        gpdf = g_to_gpdf(r, onset + peak, cls.RHO0)
        sq_true = fq_to_sq(q, gpdf_to_fq(r, gpdf, q))
        cls.sq_meas = (sq_true + 9.0) / 10.0  # true correction a=10, b=-9
        cls.q = q
        write_stog_xy(cls.run_dir / "synth.dat", q, cls.sq_meas, title="synthetic")
        (cls.run_dir / "stog.inp").write_text(
            "1\nsynth.dat\n0.60 30.0\n-9 0.1\n0\nscale.fq\nscale.gr\n25\n1000\nN\n"
            "0.05\n0\nN\nY\n1.0\nscale_ft.sq\nscale_ft.gr\n0.02\n"
            "scale_ft_rmc.fq\nscale_ft_rmc.gr\nscale_ft_rmc.dr\n2.48 2.65 3.1\n"
        )

    @classmethod
    def tearDownClass(cls):
        import shutil

        shutil.rmtree(cls.run_dir, ignore_errors=True)

    def test_inspect_reports_inp_and_header(self):
        response = self.client.post(
            "/api/scaling/preview",
            json={"path": "results/scaling_api_test/stog.inp", "inspect": True},
        )
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["kind"], "inp")
        self.assertAlmostEqual(payload["inp"]["a"], 10.0)
        self.assertAlmostEqual(payload["inp"]["b"], -9.0)
        self.assertTrue(payload["dataFile"].endswith("synth.dat"))

    def test_preview_auto_fits_and_returns_series(self):
        response = self.client.post(
            "/api/scaling/preview",
            json={"path": "results/scaling_api_test/stog.inp"},
        )
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertTrue(payload["result"]["converged"])
        self.assertLess(abs(payload["result"]["a"] - 10.0) / 10.0, 0.05)
        series = payload["series"]
        self.assertEqual(len(series["q"]), len(series["sqScaled"]))
        self.assertEqual(len(series["r"]), len(series["gk"]))
        self.assertAlmostEqual(payload["guides"]["gkLowR"], -self.B2)
        # stog.inp mode defaults to classic enforcement: flat -<b>^2 below 2.48.
        self.assertIsNotNone(payload["series"]["gkEnforced"])
        below = [
            value
            for radius, value in zip(series["r"], series["gkEnforced"])
            if radius <= 2.48
        ]
        for value in below[:: max(1, len(below) // 7)]:
            self.assertAlmostEqual(value, -self.B2, places=10)

    def test_preview_manual_data_mode(self):
        response = self.client.post(
            "/api/scaling/preview",
            json={
                "path": "results/scaling_api_test/synth.dat",
                "qmin": 0.6,
                "qmax": 30,
                "rho0": self.RHO0,
                "bAvgSq": self.B2,
                "r0": 2.65,
                "mode": "manual",
                "a": 10.0,
                "b": -9.0,
            },
        )
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["kind"], "data")
        self.assertEqual(payload["result"]["a"], 10.0)
        # Data mode auto-enforces at the data-derived first-shell onset
        # (synthetic model onset ~2.6 A); explicit enforce=False opts out.
        self.assertIsNotNone(payload["series"]["gkEnforced"])
        detected = payload["enforcement"]["cutoff"]
        self.assertGreater(detected, 2.3)
        self.assertLess(detected, 2.85)

        bare = self.client.post(
            "/api/scaling/preview",
            json={
                "path": "results/scaling_api_test/synth.dat",
                "qmin": 0.6, "qmax": 30, "rho0": self.RHO0,
                "bAvgSq": self.B2, "r0": 2.65,
                "mode": "manual", "a": 10.0, "b": -9.0,
                "enforce": False,
            },
        )
        self.assertIsNone(bare.get_json()["series"]["gkEnforced"])

    def test_preview_rejects_missing_qmin(self):
        response = self.client.post(
            "/api/scaling/preview",
            json={"path": "results/scaling_api_test/synth.dat", "qmax": 30},
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn("qmin", response.get_json()["error"])

    def test_run_writes_family_then_conflicts_then_force(self):
        body = {"path": "results/scaling_api_test/stog.inp", "mode": "manual"}
        response = self.client.post("/api/scaling/run", json=body)
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        out_dir = Path(payload["outDir"])
        self.assertEqual(out_dir, self.run_dir / "autoscale")
        for name in (
            "scale.fq", "scale_ft.sq", "scale_ft_rmc.gr", "ft.dat",
            "stog_provenance.json",
        ):
            self.assertTrue((out_dir / name).exists(), name)

        conflict = self.client.post("/api/scaling/run", json=body)
        self.assertEqual(conflict.status_code, 409)

        forced = self.client.post("/api/scaling/run", json={**body, "force": True})
        self.assertEqual(forced.status_code, 200)


if __name__ == "__main__":
    unittest.main()
