# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Flask /api/pca/sites and /api/pca/kde regressions (PCA-ellipsoid group)."""

from pathlib import Path
import json
import os
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
if str(ROOT / "tests") not in sys.path:
    sys.path.insert(0, str(ROOT / "tests"))

os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "rmc_toolkits_cache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402

from test_pca_regressions import mixed_site_file, wrapped_site_lines, write_rmc6f  # noqa: E402
import numpy as np  # noqa: E402


def strict_json(response):
    """Parse a body the way a browser's JSON.parse does: NaN/Infinity tokens are errors."""
    def reject(token):
        raise ValueError(f"non-standard JSON constant {token}")
    return json.loads(response.get_data(as_text=True), parse_constant=reject)


class PcaApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.run_dir = Path(self._tmp.name).resolve()
        backend_app.SELECTED_DATA_ROOTS.add(self.run_dir)

    def tearDown(self):
        backend_app.SELECTED_DATA_ROOTS.discard(self.run_dir)
        self._tmp.cleanup()

    def test_sites_lists_every_species_of_a_mixed_site(self):
        # pca.parity.4 / numerics.17 / parity.20 / parity.26 / numerics.33 / physics.39
        n_major, n_minor = mixed_site_file(self.run_dir / "mixed.rmc6f")
        response = self.client.get("/api/pca/sites", query_string={"dir": str(self.run_dir)})
        self.assertEqual(response.status_code, 200)
        payload = strict_json(response)
        self.assertEqual(payload["elements"], ["Ga", "In", "Se"])
        mixed = payload["sites"][0]
        self.assertEqual(mixed["element"], "Ga")
        self.assertTrue(mixed["mixed"])
        self.assertEqual(mixed["elementCounts"], {"Ga": n_major, "In": n_minor})

        kde = self.client.get("/api/pca/kde", query_string={
            "dir": str(self.run_dir), "element": "In", "grid": 8, "projections": "false"})
        self.assertEqual(kde.status_code, 200)
        self.assertEqual(strict_json(kde)["count"], n_minor)

    def _write_cloud(self, *, nan_line=None):
        supercell = (4, 4, 4)
        lines = wrapped_site_lines((0.25, 0.25, 0.25), 0.06, supercell=supercell, cell_edge=8.0, seed=1)
        if nan_line is not None:
            parts = lines[nan_line].split()
            parts[3] = "NaN"
            lines[nan_line] = " ".join(parts)
        write_rmc6f(self.run_dir / "cloud.rmc6f", lines, supercell=supercell,
                    lattice=np.diag(np.asarray(supercell, dtype=float) * 8.0))

    def test_non_finite_kde_parameters_are_a_400_not_nan_json(self):
        # pca.numerics.37
        self._write_cloud()
        for key in ("extent", "bwScale", "bw"):
            for value in ("nan", "inf"):
                with self.subTest(key=key, value=value):
                    response = self.client.get("/api/pca/kde", query_string={
                        "dir": str(self.run_dir), "referenceNumber": 1, "grid": 8, key: value})
                    self.assertEqual(response.status_code, 400)
                    self.assertIn("finite", strict_json(response)["error"])

    def test_a_nan_coordinate_is_a_clear_400(self):
        # pca.parity.6 / parity.23 / parity.27
        self._write_cloud(nan_line=9)
        response = self.client.get("/api/pca/sites", query_string={"dir": str(self.run_dir)})
        self.assertEqual(response.status_code, 400)
        self.assertRegex(strict_json(response)["error"], "atom 10 .*non-finite")

    def test_invalid_probability_is_a_400(self):
        self._write_cloud()
        response = self.client.get("/api/pca/sites", query_string={"dir": str(self.run_dir), "probability": "nan"})
        self.assertEqual(response.status_code, 400)


if __name__ == "__main__":
    unittest.main()
