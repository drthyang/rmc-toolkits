# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Auto StoG scaling API (/api/scaling/*): stog.inp r0, enforcement flag, cache hygiene.

Regression tests for the 1.0 audit (stog-a group).
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
        backend_app._cached_scaling.cache_clear()

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


if __name__ == "__main__":
    unittest.main()
