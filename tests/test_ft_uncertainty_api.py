# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Flask's /api/plot/data gives an xPDF plot its error ribbon (``uncertainty``).

The payload is ``ft_uncertainty_for_run``'s, so it is the golden the browser's
port is held to (ft_uncertainty_fixture.json): a server folder and the same
files picked in the browser draw the same ribbon. A companion that cannot be
read is a failure payload on a plot that still loads; the input file is read
only inside the data roots; a file's own σ column is its experimental
series' error bar.
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
    from generate_ft_uncertainty_fixture import OUT, STEM, committed_cases
    from golden_compare import assert_same
except ModuleNotFoundError:  # `python -m unittest tests.test_ft_uncertainty_api`
    from tests.generate_ft_uncertainty_fixture import OUT, STEM, committed_cases
    from tests.golden_compare import assert_same

DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"


class FtUncertaintyApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        cls.fixture = json.loads(OUT.read_text(encoding="utf-8"))
        cls.cases = {case["name"]: case for case in committed_cases(cls.fixture)}
        cls.expected = {case["name"]: case["expected"] for case in cls.fixture["cases"]}

    def setUp(self):
        backend_app._FT_UNCERTAINTY_CACHE.clear()
        # Inside the data root (the repo root), as a user's run folder is.
        self.folder = Path(tempfile.mkdtemp(dir=ROOT))
        self.addCleanup(shutil.rmtree, self.folder, ignore_errors=True)

    def write(self, files: dict) -> Path:
        for name, text in files.items():
            (self.folder / name).write_text(text, encoding="utf-8")
        return self.folder / f"{STEM}_FT_XFQ1.csv"

    def plot(self, path: Path) -> dict:
        response = self.client.get("/api/plot/data", query_string={"path": str(path)})
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        return response.get_json()

    def test_the_demo_plot_says_why_it_has_no_ribbon(self):
        payload = self.plot(DEMO / "GTS_250K_FT_XFQ1.csv")
        assert_same(self, payload["uncertainty"], self.fixture["demo"]["expected"], "demo")
        self.assertEqual([series["label"] for series in payload["series"]], ["X_ray-calc", "X_ray_exp_renorm"])
        self.assertTrue(all("sigma" not in series for series in payload["series"]))

    def test_every_fixture_case_is_the_payload_flask_sends(self):
        for name, case in self.cases.items():
            with self.subTest(case=name):
                backend_app._FT_UNCERTAINTY_CACHE.clear()
                for stale in self.folder.iterdir():
                    stale.unlink()
                payload = self.plot(self.write(case["files"]))
                assert_same(self, payload["uncertainty"], self.expected[name], name)

    def test_other_plot_kinds_carry_no_uncertainty(self):
        for name in ("GTS_250K_FQ1.csv", "GTS_250K_PDFpartials.csv"):
            self.assertNotIn("uncertainty", self.plot(DEMO / name))

    def test_a_changed_companion_is_read_again(self):
        files = self.cases["propagated, CONVOLVE"]["files"]
        path = self.write(files)
        self.assertEqual(self.plot(path)["uncertainty"]["status"], "propagated")
        without_sigma = self.cases["no error column"]["files"]["scale_ft_rmc.fq"]
        (self.folder / "scale_ft_rmc.fq").write_text(without_sigma + "\n", encoding="utf-8")
        again = self.plot(path)["uncertainty"]
        self.assertEqual((again["status"], again["reason"]), (
            "unavailable", "F(Q) carries no uncertainties: scale_ft_rmc.fq has no error column"))

    def test_an_input_outside_the_data_root_reads_as_absent(self):
        files = dict(self.cases["propagated, CONVOLVE"]["files"])
        outside = Path(tempfile.mkdtemp())  # not under the repo root
        self.addCleanup(shutil.rmtree, outside, ignore_errors=True)
        (outside / "elsewhere.fq").write_text(files.pop("scale_ft_rmc.fq"), encoding="utf-8")
        files[f"{STEM}.dat"] = files[f"{STEM}.dat"].replace(
            "FILENAME :: scale_ft_rmc.fq", f"FILENAME :: {outside / 'elsewhere.fq'}")
        payload = self.plot(self.write(files))["uncertainty"]
        self.assertEqual(payload["status"], "unavailable")
        self.assertTrue(payload["reason"].startswith("F(Q) carries no uncertainties: "), payload["reason"])
        self.assertTrue(payload["reason"].endswith("is not in the run folder"), payload["reason"])

    def test_an_unreadable_companion_is_a_failure_on_a_plot_that_loads(self):
        files = dict(self.cases["propagated, CONVOLVE"]["files"])
        files[f"{STEM}_XFQ1.csv"] = "Q , F(Q)_RMC, F(Q)_Expt\n1.0 , oops , 2.0\n"
        payload = self.plot(self.write(files))
        self.assertEqual(payload["uncertainty"]["status"], "failed")
        self.assertTrue(payload["uncertainty"]["reason"].startswith("could not compute the error ribbon: "))
        self.assertIn("'oops' is not a number", payload["uncertainty"]["reason"])
        self.assertEqual(len(payload["series"]), 2)

    def test_a_sigma_column_is_the_experimental_series_error_bar(self):
        text = "r(A) , X_ray-calc , X_ray_exp_renorm , sigma\n1.0 , 0.1 , 0.2 , 0.01\n1.1 , 0.3 , 0.25 , 0.02\n"
        path = self.folder / "run_FT_XFQ1.csv"
        path.write_text(text, encoding="utf-8")
        payload = self.plot(path)
        self.assertEqual([series["label"] for series in payload["series"]], ["X_ray-calc", "X_ray_exp_renorm"])
        self.assertEqual(payload["series"][1]["sigma"], [0.01, 0.02])
        self.assertNotIn("sigma", payload["series"][0])
        # A column with a zero σ where G is measured is not used.
        path.write_text(text.replace("0.02\n", "0.0\n"), encoding="utf-8")
        payload = self.plot(path)
        self.assertEqual([series["label"] for series in payload["series"]], ["X_ray-calc", "X_ray_exp_renorm", "sigma"])
        self.assertTrue(all("sigma" not in series for series in payload["series"]))


if __name__ == "__main__":
    unittest.main()
