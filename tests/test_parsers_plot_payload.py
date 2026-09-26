# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Flask /api/plot/data payloads: labels name the function the file holds.

The committed demo run's *_FQ1.csv holds F(Q) (header ``F(Q)_RMC, F(Q)_Expt``;
-> 0 at high Q) and *_PDFpartials.csv the partial g_ij(r) (-> 1 at large r);
they used to be labelled S(Q) / G(r). The browser port is pinned on the same
files in web_app/frontend/src/__tests__/plotLabels.test.js.
"""

from __future__ import annotations

import json
import math
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import generate_plot_parity_fixture as parity  # noqa: E402


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


class PlotParityFixtureTests(unittest.TestCase):
    """Non-x-ray RMCProfile layouts: the committed Flask golden the browser is held to."""

    @classmethod
    def setUpClass(cls):
        cls.committed = json.loads(parity.OUT.read_text(encoding="utf-8"))["cases"]

    def test_committed_fixture_matches_the_current_flask_payloads(self):
        regenerated = json.loads(json.dumps(parity.flask_payloads(), ensure_ascii=False))
        # Floats compare to 1e-6 relative: a value rounded to 7 significant digits can land on
        # either side of a rounding boundary across numpy versions (summation order).
        stale = "plot_parity_fixture.json is stale: run `python tests/generate_plot_parity_fixture.py`"
        self._assert_same(regenerated, self.committed, stale, "$")

    def _assert_same(self, got, want, message, path):
        if isinstance(want, float) or isinstance(got, float):
            self.assertIsInstance(got, (int, float), f"{message} ({path})")
            self.assertTrue(math.isclose(got, want, rel_tol=1e-6, abs_tol=1e-12), f"{message} ({path}: {got} != {want})")
        elif isinstance(want, dict):
            self.assertIsInstance(got, dict, f"{message} ({path})")
            self.assertEqual(sorted(got), sorted(want), f"{message} ({path} keys)")
            for key in want:
                self._assert_same(got[key], want[key], message, f"{path}.{key}")
        elif isinstance(want, list):
            self.assertIsInstance(got, list, f"{message} ({path})")
            self.assertEqual(len(got), len(want), f"{message} ({path} length)")
            for index, (g, w) in enumerate(zip(got, want)):
                self._assert_same(g, w, message, f"{path}[{index}]")
        else:
            self.assertEqual(got, want, f"{message} ({path})")

    def test_rwp_column_roles_match_the_construction(self):
        # Independent of both implementations: the case records which column is the
        # calculation and which the experiment; compute R = ||calc-expt|| / ||expt||
        # over the rows finite in both, straight from the text.
        for case in self.committed:
            if case["truth"] is None:
                continue
            with self.subTest(name=case["name"]):
                rows = []
                for line in case["text"].splitlines():
                    cells = [cell.strip() for cell in line.split(",") if cell.strip()]
                    if not cells:
                        continue
                    try:
                        rows.append([float(cell) for cell in cells])
                    except ValueError:
                        continue
                calc = [row[case["truth"]["calculated"]] for row in rows]
                expt = [row[case["truth"]["experimental"]] for row in rows]
                pairs = [(c, e) for c, e in zip(calc, expt) if math.isfinite(c) and math.isfinite(e)]
                expected = math.sqrt(sum((c - e) ** 2 for c, e in pairs) / sum(e * e for _, e in pairs))
                self.assertAlmostEqual(case["expected"]["metrics"]["rwp"], expected, places=12)


if __name__ == "__main__":
    unittest.main()
