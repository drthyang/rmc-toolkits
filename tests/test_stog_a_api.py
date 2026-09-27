# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Auto StoG scaling API (/api/scaling/*): stog.inp r0, enforcement flag, cache hygiene.

Regression tests for the 0.6.0 audit (stog-a group).
"""

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
RUN = ROOT / "results" / "stog_a_api_test"
REL = "results/stog_a_api_test"

INP = (
    "1\n{data}\n0.6 29.4\n-9 0.1\n0\nscale.fq\nscale.gr\n25\n1000\nN\n"
    "0.05\n0\nN\nY\n1.0\nscale_ft.sq\nscale_ft.gr\n0.02\n"
    "scale_ft_rmc.fq\nscale_ft_rmc.gr\nscale_ft_rmc.dr\n{peak_line}\n"
)


def synthetic_sq(onset):
    """Repo synthetic model with its first-shell onset at `onset` (peak +0.15 A); a=10, b=-9."""
    q = np.arange(20, 981) * 0.03
    r = np.arange(1, 12001) * 0.005
    g = 0.5 * (1.0 + np.tanh((r - onset) / 0.07))
    g = g + 1.6 * np.exp(-0.5 * ((r - (onset + 0.15)) / 0.15) ** 2)
    sq_true = fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, RHO0), q))
    return q, (sq_true + 9.0) / 10.0


class ScalingApiStogATests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        RUN.mkdir(parents=True, exist_ok=True)
        q, sq = synthetic_sq(1.70)
        write_stog_xy(RUN / "short.sq", q, sq, title="synthetic onset 1.70")
        (RUN / "short.inp").write_text(INP.format(data="short.sq", peak_line="2.3 1.6 2.2"))
        q, sq = synthetic_sq(2.65)
        write_stog_xy(RUN / "synth.dat", q, sq, title="synthetic onset 2.65")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(RUN, ignore_errors=True)

    def setUp(self):
        backend_app._SCALING_CACHE.clear()

    def data_body(self, **extra):
        return {
            "path": f"{REL}/synth.dat", "qmin": 0.6, "qmax": 30, "rho0": RHO0,
            "bAvgSq": B2, "rmax": 25, "nr": 1000, **extra,
        }

    def test_inp_r0_uses_the_first_peak_window_start(self):
        response = self.client.post("/api/scaling/preview", json={"path": f"{REL}/short.inp"})
        self.assertEqual(response.status_code, 200, response.get_json())
        payload = response.get_json()
        # Pre-fix: r0 = max(2.3, 1.6) = 2.3 and a fit across the first peak.
        self.assertAlmostEqual(payload["provenance"]["config"]["r0"], 1.6)
        self.assertLess(abs(payload["result"]["a"] / 10.0 - 1.0), 0.03)

    def test_data_mode_honours_an_explicit_enforce_cutoff(self):
        # Pre-fix: without enforce:true the cutoff was discarded and the RMC G(r)
        # flattened out to the auto-detected onset instead.
        response = self.client.post("/api/scaling/preview", json=self.data_body(enforceCutoff=2.0))
        payload = response.get_json()
        self.assertEqual(response.status_code, 200, payload)
        self.assertEqual(payload["enforcement"], {"cutoff": 2.0, "peakRmin": 2.0, "peakRmax": 2.0})
        r = np.asarray(payload["series"]["r"])
        gk = np.asarray(payload["series"]["gk"])
        enforced = np.asarray(payload["series"]["gkEnforced"])
        above = r > 2.0
        np.testing.assert_allclose(enforced[above], gk[above], atol=1e-15)
        np.testing.assert_allclose(enforced[~above], -B2, atol=1e-15)

    def test_run_writes_the_explicit_cutoff(self):
        body = self.data_body(enforceCutoff=2.0, outDir=f"{REL}/out_cutoff", force=True)
        response = self.client.post("/api/scaling/run", json=body)
        self.assertEqual(response.status_code, 200, response.get_json())
        import json

        provenance = json.loads((RUN / "out_cutoff" / "synth_provenance.json").read_text())
        self.assertEqual(provenance["enforcement"]["cutoff"], 2.0)

    def test_string_and_numeric_false_disable_enforcement(self):
        # Pre-fix: _payload_bool read "false" as off for the explicit cutoff, but
        # the auto branch tested `is not False`, so "false"/0 re-enabled it.
        for flag in ("false", "0", 0, False, "no"):
            with self.subTest(enforce=flag):
                response = self.client.post("/api/scaling/preview", json=self.data_body(enforce=flag))
                payload = response.get_json()
                self.assertEqual(response.status_code, 200, payload)
                self.assertIsNone(payload["enforcement"])
                self.assertIsNone(payload["series"]["gkEnforced"])
        inp = self.client.post(
            "/api/scaling/preview", json={"path": f"{REL}/short.inp", "enforce": "false"}
        ).get_json()
        self.assertIsNone(inp["enforcement"])

    def test_default_and_true_still_enforce(self):
        for extra in ({}, {"enforce": True}, {"enforce": "true"}, {"enforce": ""}):
            with self.subTest(**extra):
                payload = self.client.post(
                    "/api/scaling/preview", json=self.data_body(**extra)
                ).get_json()
                self.assertIsNotNone(payload["enforcement"])

    def test_contradictory_flags_are_rejected(self):
        response = self.client.post(
            "/api/scaling/preview", json=self.data_body(enforce=False, enforceCutoff=2.0)
        )
        self.assertEqual(response.status_code, 400)
        response = self.client.post(
            "/api/scaling/preview", json=self.data_body(peakWindow=[2.1, 2.5])
        )
        self.assertEqual(response.status_code, 400)

    def test_cached_result_is_never_mutated(self):
        # Manual data mode with a pinned r0: scale_pipeline records no
        # r0_detected, and the auto-enforcement branch used to write it into the
        # shared lru-cached ScalingResult, leaking into later enforce=false calls.
        body = self.data_body(r0=2.65, mode="manual", a=10.0, b=-9.0)
        computed = []
        compute = backend_app._compute_scaling

        def counting_compute(*args):
            computed.append(args)
            return compute(*args)

        backend_app._compute_scaling = counting_compute
        self.addCleanup(setattr, backend_app, "_compute_scaling", compute)

        def diagnostics(**extra):
            payload = self.client.post("/api/scaling/preview", json={**body, **extra}).get_json()
            return payload["diagnostics"], payload["guides"]["r0Detected"], payload["provenance"]

        cold, cold_guide, cold_provenance = diagnostics(enforce=False)
        self.assertNotIn("r0_detected", cold)
        self.assertIsNone(cold_guide)
        auto, auto_guide, _ = diagnostics()
        self.assertIn("r0_detected", auto)
        self.assertIsNotNone(auto_guide)
        again, again_guide, again_provenance = diagnostics(enforce=False)
        self.assertEqual(again, cold)
        self.assertIsNone(again_guide)
        self.assertNotIn("r0_detected", again_provenance)
        self.assertEqual(len(computed), 1)  # the same cached object served all three


if __name__ == "__main__":
    unittest.main()
