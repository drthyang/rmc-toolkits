# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Auto StoG scaling API (/api/scaling/*): 1.0 audit regressions, stog-b group."""

from pathlib import Path
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))

os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "rmc_toolkits_cache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402

from rmc_toolkits.parsers import write_stog_xy  # noqa: E402
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq  # noqa: E402

RHO0, B2 = 0.05, 0.02
RUN = ROOT / "results" / "stog_b_api_test"
REL = "results/stog_b_api_test"


def model_sq():
    """Repo synthetic model (true S(Q)) on the fixture Q grid."""
    q = np.arange(20, 981) * 0.03
    r = np.arange(1, 12001) * 0.005
    g = 0.5 * (1.0 + np.tanh((r - 2.65) / 0.07))
    g = g + 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    return q, fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, RHO0), q))


class ScalingApiStogBTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        RUN.mkdir(parents=True, exist_ok=True)
        q, sq_true = model_sq()
        write_stog_xy(RUN / "desc.dat", q[::-1], ((sq_true + 9.0) / 10.0)[::-1])
        write_stog_xy(RUN / "inverted.dat", q, 2.0 - sq_true)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(RUN, ignore_errors=True)

    def setUp(self):
        backend_app._SCALING_CACHE.clear()

    def body(self, name, **extra):
        return {
            "path": f"{REL}/{name}", "qmin": 0.6, "qmax": 30, "rho0": RHO0,
            "bAvgSq": B2, "r0": 2.5, "rmax": 25, "nr": 1000, **extra,
        }

    def test_descending_file_scales_like_the_ascending_one(self):
        response = self.client.post("/api/scaling/preview", json=self.body("desc.dat"))
        self.assertEqual(response.status_code, 200, response.get_json())
        result = response.get_json()["result"]
        self.assertAlmostEqual(result["a"], 9.969849853119817, places=9)  # the parity fixture's auto a
        self.assertTrue(result["converged"])

    def test_run_refuses_a_non_positive_scale(self):
        preview = self.client.post("/api/scaling/preview", json=self.body("inverted.dat"))
        self.assertEqual(preview.status_code, 200)
        body = preview.get_json()
        self.assertFalse(body["result"]["converged"])
        self.assertIn("non-physical scale", body["diagnostics"]["fit_failure"])
        out = RUN / "out_inverted"
        response = self.client.post(
            "/api/scaling/run", json=self.body("inverted.dat", outDir=f"{REL}/out_inverted"),
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn("auto-fit failed", response.get_json()["error"])
        self.assertFalse(out.exists() and any(out.iterdir()))

    def test_formula_b_sq_avg_is_not_paired_with_another_b_avg_sq(self):
        # x-ray style <b>^2 = 1 with a (neutron) formula: CLI parity, no mixed pair.
        body = self.body("desc.dat", bAvgSq=1.0, formula="FeCoSn")
        response = self.client.post("/api/scaling/preview", json=body)
        self.assertEqual(response.status_code, 200, response.get_json())
        config = response.get_json()["provenance"]["config"]
        self.assertEqual(config["b_avg_sq"], 1.0)
        self.assertIsNone(config["b_sq_avg"])
        # The dropped <b^2> is said, as the CLI prints it (integration item 4).
        (warning,) = response.get_json()["warnings"]
        self.assertIn("NOT the formula's <b^2>", warning)
        out = f"{REL}/out_mixed"
        run = self.client.post("/api/scaling/run", json={**body, "outDir": out})
        self.assertEqual(run.status_code, 200, run.get_json())
        self.assertEqual(run.get_json()["warnings"], [warning])
        clean = self.client.post("/api/scaling/preview", json=self.body("desc.dat"))
        self.assertEqual(clean.get_json()["warnings"], [])
        impossible = self.body("desc.dat", bAvgSq=1.0, bSqAvg=0.5)
        response = self.client.post("/api/scaling/preview", json=impossible)
        self.assertEqual(response.status_code, 400)
        self.assertIn("Cauchy-Schwarz", response.get_json()["error"])


if __name__ == "__main__":
    unittest.main()
