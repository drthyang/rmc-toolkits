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

from test_pca_regressions import mixed_site_file  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
