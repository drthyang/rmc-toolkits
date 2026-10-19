# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The analysis request contract on the Flask adapter.

``web_app/frontend/src/__tests__/fixtures/analysis_requests.json`` says what each
kind of ``requestAnalysis`` (useSiteCloud.js) accepts and refuses. The JS test
(analysisRequests.test.js) replays it on the worker and in-process adapters;
this one replays every case through the Flask test client, sending the params
as axios serialises them: null dropped, booleans 'true'/'false', numbers as JS
writes them (a whole float as an integer). An accepted case needs 200 and its
``fields``; a refusal needs 400 with an ``error`` naming the parameter. Where
the adapters differ (``{flask, worker}``, each kept with its reason) the
``flask`` outcome applies. The cancel probes and an unknown kind are refused at
requestAnalysis before any request is sent, so they have no Flask replay.
"""

from pathlib import Path
import json
import os
import shutil
import sys
import tempfile
import unittest
import warnings


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))

os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "rmc_toolkits_cache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402


FIXTURE_PATH = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "analysis_requests.json"
FIXTURE = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))
RUN = "results/analysis_requests_contract"


def axios_value(value) -> str:
    """One query value as axios writes it (JS String(value))."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def axios_query(params: dict) -> dict:
    return {key: axios_value(value) for key, value in params.items() if value is not None}


def field(payload, path: str):
    """A payload value by dotted path (``.length`` of a list); a missing value is None."""
    value = payload
    for part in path.split("."):
        if isinstance(value, list):
            value = len(value) if part == "length" else value[int(part)]
        elif isinstance(value, dict):
            value = value.get(part)
        else:
            return None
    return value


def flask_outcome(case: dict):
    """('ok', fields) or ('error', name): the case's expectation on the Flask adapter."""
    expect = case["expect"]
    split = isinstance(expect, dict) and "flask" in expect
    expected = expect["flask"] if split else expect
    if expected == "ok":
        return "ok", ({} if split else case.get("fields", {}))
    if "fields" in expected:
        return "ok", expected["fields"]
    return "error", expected["error"]


def strict_json(body: str):
    """Parse like a browser's JSON.parse: NaN/Infinity tokens are an error."""

    def reject(token):
        raise ValueError(f"non-standard JSON token {token}")

    return json.loads(body, parse_constant=reject)


class AnalysisRequestContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        directory = ROOT / RUN
        directory.mkdir(parents=True, exist_ok=True)
        (directory / "run.rmc6f").write_text(FIXTURE["rmc6f"], encoding="utf-8")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(ROOT / RUN, ignore_errors=True)

    def test_every_route_exists(self):
        served = {rule.rule: rule.methods for rule in backend_app.app.url_map.iter_rules()}
        for kind, route in FIXTURE["routes"].items():
            with self.subTest(kind=kind):
                self.assertIn(route, served)
                self.assertIn("GET", served[route])

    def test_the_apps_copies_of_the_caps_equal_the_contracts(self):
        limits = {(limit["kind"], limit["param"]): limit for limit in FIXTURE["limits"]}
        self.assertEqual(tuple(limits[("kde", "grid")]["clamp"]), backend_app.PCA_GRID_CLAMP)
        self.assertEqual(limits[("orientation", "smoothing")]["le"], backend_app.MAX_ORIENTATION_SMOOTHING)
        self.assertEqual(limits[("clusters", "rMax")]["le"], backend_app.MAX_CLUSTER_WINDOW)

    def test_every_case_on_the_flask_adapter(self):
        replayed = 0
        for case in FIXTURE["cases"]:
            if case.get("probe") or case["kind"] not in FIXTURE["routes"]:
                # Refused at requestAnalysis: never reaches a route.
                self.assertTrue(case.get("probe") or case["expect"] == {"error": "kind"}, case["name"])
                continue
            with self.subTest(case=case["name"]):
                outcome, expected = flask_outcome(case)
                query = {**axios_query(case["params"]), "dir": RUN}
                with warnings.catch_warnings():
                    # Accelerate's BLAS warns spuriously inside matmul on some Macs.
                    warnings.simplefilter("ignore", RuntimeWarning)
                    response = self.client.get(FIXTURE["routes"][case["kind"]], query_string=query)
                body = response.get_data(as_text=True)
                if outcome == "ok":
                    self.assertEqual(response.status_code, 200, body[:300])
                    payload = strict_json(body)
                    for path, value in expected.items():
                        self.assertEqual(field(payload, path), value, path)
                else:
                    self.assertEqual(response.status_code, 400, body[:300])
                    self.assertIn(expected, strict_json(body)["error"])
                replayed += 1
        self.assertGreater(replayed, 100)


if __name__ == "__main__":
    unittest.main()
