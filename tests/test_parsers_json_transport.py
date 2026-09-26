# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Every Flask JSON body must be strict JSON: no bare NaN / Infinity tokens.

RMCProfile outputs reach the readers with NaN in masked regions (and a blown-up
run can log NaN chi^2). Python's float() keeps those values; the default Flask
provider then wrote them as the bare tokens ``NaN`` / ``Infinity``, which
``JSON.parse`` rejects, so axios handed the chart a raw string and the Flask
dashboard lost the plot. The app-wide provider maps non-finite floats to null.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))

os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402


def _reject_constant(token: str):
    raise ValueError(f"non-standard JSON token {token}")


def strict_json(response) -> object:
    """Parse a response body the way a browser's JSON.parse does (no NaN/Infinity)."""
    return json.loads(response.get_data(as_text=True), parse_constant=_reject_constant)


LOG_HEADER = (
    "Time       moves_acc moves_gen     F(Q)_1  X_ray_(R)1\n"
    "h/m/s/.th    WEIGHT PARAMETERS   0.100E+01  0.100E+01\n"
)


class StrictJsonTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(dir=ROOT)
        self.directory = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def test_masked_csv_series_are_strict_json_with_null_gaps(self):
        path = self.directory / "run_FQ1.csv"
        rows = ["Q, F(Q)_RMC, F(Q)_Expt"]
        for index in range(6):
            calc = "NaN" if index in (2, 3) else f"{0.1 * index:.3f}"
            rows.append(f"{0.5 + index:.2f}, {calc}, {0.1 * index + 0.01:.3f}")
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")

        response = self.client.get("/api/plot/data", query_string={"path": str(path)})

        self.assertEqual(response.status_code, 200)
        payload = strict_json(response)
        calculated = payload["series"][0]["y"]
        self.assertIsNone(calculated[2])
        self.assertIsNone(calculated[3])
        self.assertAlmostEqual(calculated[1], 0.1)
        self.assertTrue(math.isfinite(payload["metrics"]["rwp"]))

    def test_fully_masked_column_gives_null_rwp_and_strict_series(self):
        path = self.directory / "run_FQ1.csv"
        path.write_text(
            "Q, F(Q)_RMC, F(Q)_Expt\n1.0, 0.1, NaN\n2.0, 0.2, NaN\n", encoding="utf-8"
        )

        for endpoint in ("/api/plot/data", "/api/plot/metadata"):
            with self.subTest(endpoint=endpoint):
                response = self.client.get(endpoint, query_string={"path": str(path)})
                self.assertEqual(response.status_code, 200)
                payload = strict_json(response)
                self.assertIsNone(payload["metrics"]["rwp"])

    def test_log_with_non_finite_chi_is_strict_json(self):
        path = self.directory / "run-00.log"
        path.write_text(
            LOG_HEADER
            + "1.0  10  20  0.300E-02  0.200E-03\n"
            + "2.0  20  40  0.300E-02  NaN\n",
            encoding="utf-8",
        )

        for endpoint in ("/api/plot/data", "/api/plot/metadata"):
            with self.subTest(endpoint=endpoint):
                response = self.client.get(endpoint, query_string={"path": str(path)})
                self.assertEqual(response.status_code, 200)
                payload = strict_json(response)
                self.assertIsNone(payload["metrics"]["final_chi_r"])

    def test_provider_maps_numpy_and_nested_non_finite_values_to_null(self):
        with backend_app.app.app_context():
            body = backend_app.app.json.dumps(
                {
                    "a": [1.0, float("nan"), float("inf"), -float("inf")],
                    "b": {"c": np.float64("nan"), "d": np.array([1.5, np.nan])},
                    "e": np.int64(3),
                }
            )
        parsed = json.loads(body, parse_constant=_reject_constant)
        self.assertEqual(parsed, {"a": [1.0, None, None, None], "b": {"c": None, "d": [1.5, None]}, "e": 3})


if __name__ == "__main__":
    unittest.main()
