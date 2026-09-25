# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Flask /api/plot/data payloads: labels name the function the file holds.

The committed demo run's *_FQ1.csv holds F(Q) (header ``F(Q)_RMC, F(Q)_Expt``;
-> 0 at high Q) and *_PDFpartials.csv the partial g_ij(r) (-> 1 at large r);
they used to be labelled S(Q) / G(r). The browser port is pinned on the same
files in web_app/frontend/src/__tests__/plotLabels.test.js.
"""

from __future__ import annotations

import os
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402


class PlotLabelApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

    def _data(self, path: Path) -> dict:
        response = self.client.get("/api/plot/data", query_string={"path": str(path)})
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return response.get_json()

    def test_demo_fq_is_f_of_q(self):
        payload = self._data(DEMO / "GTS_250K_FQ1.csv")
        self.assertEqual((payload["kind"], payload["title"], payload["yLabel"]), ("xray_sq", "F(Q)", "F(Q)"))
        self.assertEqual([series["label"] for series in payload["series"]], ["F(Q)_RMC", "F(Q)_Expt"])

    def test_demo_partials_are_g_of_r(self):
        payload = self._data(DEMO / "GTS_250K_PDFpartials.csv")
        self.assertEqual((payload["title"], payload["yLabel"], payload["xLabel"]), ("Partial g(r)", "g(r)", "r (Å)"))

    def test_second_reciprocal_dataset_is_charted(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmpdir:
            path = Path(tmpdir) / "run_SQ2.csv"
            path.write_text("Q, S(Q)_RMC, S(Q)_Expt\n1.0, 0.7, 1.0\n2.0, 1.4, 2.0\n", encoding="utf-8")
            listing = self.client.get("/api/files", query_string={"dir": tmpdir}).get_json()
            payload = self._data(path)
        self.assertEqual({item["name"]: item["plotKind"] for item in listing["files"]}["run_SQ2.csv"], "neutron_sq")
        self.assertEqual((payload["title"], payload["yLabel"]), ("S(Q) #2", "S(Q)"))
        self.assertAlmostEqual(payload["metrics"]["rwp"], 0.3)

    def test_stog_fq_is_f_of_q(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as tmpdir:
            path = Path(tmpdir) / "scale_ft_rmc.fq"
            path.write_text("2\ntitle\n0.5 -0.9\n1.0 -0.5\n", encoding="utf-8")
            payload = self._data(path)
        self.assertEqual(payload["yLabel"], "F(Q)")


if __name__ == "__main__":
    unittest.main()
