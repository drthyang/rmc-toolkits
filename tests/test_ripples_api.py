# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Flask's /api/ripples is ripples_for_run: the golden the browser's port is held to.

A server folder and the same files picked in the browser give the Copilot's
check_ripples the same result (ripples_fixture.json); the arguments are
validated (a bad one is a 400 that names it); an input file is read only
inside the data roots.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_ripples_fixture import OUT, committed_cases
    from golden_compare import assert_same
except ModuleNotFoundError:  # `python -m unittest tests.test_ripples_api`
    from tests.generate_ripples_fixture import OUT, committed_cases
    from tests.golden_compare import assert_same

DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"
RIPPLES_RTOL = 1e-7


def query(paths, args: dict) -> list:
    out = [("path", str(path)) for path in paths]
    if args.get("features"):
        out.append(("features", ",".join(repr(value) for value in args["features"])))
    if args.get("rRange"):
        out += [("rMin", repr(args["rRange"][0])), ("rMax", repr(args["rRange"][1]))]
    if args.get("split"):
        out.append(("split", ",".join(repr(value) for value in args["split"])))
    for key in ("pair", "qmax", "qmin", "window"):
        if args.get(key) is not None:
            out.append((key, str(args[key])))
    return out


class RipplesApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        cls.fixture = json.loads(OUT.read_text(encoding="utf-8"))

    def setUp(self):
        backend_app._RIPPLES_CACHE.clear()
        self.folder = Path(tempfile.mkdtemp(dir=ROOT))
        self.addCleanup(shutil.rmtree, self.folder, ignore_errors=True)

    def get(self, paths, args):
        return self.client.get("/api/ripples", query_string=query(paths, args))

    def test_the_demo_is_the_golden(self):
        for call in self.fixture["demo"]["calls"]:
            response = self.get([DEMO / name for name in self.fixture["demo"]["fits"]], call["args"])
            self.assertEqual(response.status_code, 200, response.get_json())
            assert_same(self, response.get_json(), call["expected"], "demo", rtol=RIPPLES_RTOL)

    def test_every_synthetic_case_is_the_golden(self):
        for case, committed in zip(committed_cases(self.fixture), self.fixture["cases"]):
            folder = self.folder / f"case{len(list(self.folder.iterdir()))}"
            for name, text in case["files"].items():
                path = folder / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text, encoding="utf-8")
            for call in committed["calls"]:
                with self.subTest(case=case["name"], args=call["args"]):
                    response = self.get([folder / name for name in case["fits"]], call["args"])
                    self.assertEqual(response.status_code, 200, response.get_json())
                    got = response.get_json()
                    # A failure message names the file, never the server's folder.
                    for dataset in got["datasets"]:
                        self.assertNotIn(str(folder), dataset.get("reason") or "")
                    assert_same(self, got, call["expected"], "expected", rtol=RIPPLES_RTOL)

    def test_bad_arguments_are_400s_that_say_why(self):
        fit = DEMO / "GTS_250K_FT_XFQ1.csv"
        for args, needle in (
            ({}, None),
            ({"features": "x"}, "features"),
            ({"rMin": "3"}, "rMin and rMax"),
            ({"rMin": "3", "rMax": "2"}, "rMax"),
            ({"split": "2.5"}, "split"),
            ({"qmax": "-1"}, "qmax"),
            ({"window": "hann"}, "window"),
            ({"features": ",".join(["2"] * 9)}, "at most"),
        ):
            with self.subTest(args=args):
                response = self.client.get("/api/ripples", query_string=[("path", str(fit)), *args.items()])
                if needle is None:
                    self.assertEqual(response.status_code, 200)
                else:
                    self.assertEqual(response.status_code, 400)
                    self.assertIn(needle, response.get_json()["error"])
        self.assertEqual(self.client.get("/api/ripples").status_code, 400)
        not_fit = self.client.get("/api/ripples", query_string=[("path", str(DEMO / "GTS_250K.dat"))])
        self.assertEqual(not_fit.status_code, 400)
        outside = self.client.get("/api/ripples", query_string=[("path", "/etc/hosts")])
        self.assertIn(outside.status_code, (403, 404))

    def test_an_input_file_outside_the_data_roots_reads_as_absent(self):
        case = next(c for c in committed_cases(self.fixture) if c["name"].startswith("a neutron G(r) from a stog"))
        for name, text in case["files"].items():
            (self.folder / name).write_text(text, encoding="utf-8")
        dat = self.folder / "rmc_run.dat"
        dat.write_text(dat.read_text().replace("FILENAME :: scale_ft_rmc.gr", "FILENAME :: /etc/hosts"))
        response = self.get([self.folder / "rmc_run_PDF1.csv"], {"features": [2.26]})
        self.assertEqual(response.status_code, 200)
        dataset = response.get_json()["datasets"][0]
        self.assertNotIn("input", dataset["transform"] or {})


if __name__ == "__main__":
    unittest.main()
