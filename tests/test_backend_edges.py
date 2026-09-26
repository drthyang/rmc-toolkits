# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Flask edge cases left after the 1.0 gate.

Each is a clear 4xx with a JSON error naming the problem -- never a 200 with
an empty result, a half-written file, an HTML page or a 500. One class per
case; the class docstrings say what used to happen.
"""

from pathlib import Path
import json
import os
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402
from test_backend_validation import strict_json, write_synthetic_rmc6f  # noqa: E402

RUN = "results/backend_edges"


class _EdgeCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        cls.run_dir = ROOT / RUN
        shutil.rmtree(cls.run_dir, ignore_errors=True)
        cls.run_dir.mkdir(parents=True)
        write_synthetic_rmc6f(cls.run_dir / "synthetic.rmc6f", supercell=(3, 3, 3))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.run_dir, ignore_errors=True)

    def status_error(self, response):
        body = response.get_data(as_text=True)
        self.assertTrue(response.is_json, f"not JSON: {response.content_type} {body[:120]}")
        return response.status_code, strict_json(body).get("error", "")

    def folder_with(self, name: str, text: str) -> str:
        folder = self.run_dir / name
        folder.mkdir(exist_ok=True)
        (folder / "run.rmc6f").write_text(text, encoding="utf-8")
        return f"{RUN}/{name}"


EMPTY_ATOMS = (
    "Supercell dimensions: 1 1 1\nLattice vectors (Ang):\n8 0 0\n0 8 0\n0 0 8\nAtoms:\n"
)


class KdeSliceElementTests(_EdgeCase):
    """/api/kde/slice: an unknown element (or no parseable atom) was a 200 all-zero map
    reading "No atoms in this slab."; element case is now ignored as on /api/pca/*."""

    def slice(self, directory=RUN, **params):
        return self.client.get(
            "/api/kde/slice", query_string={"dir": directory, "grid": 32, "levels": 0, **params}
        )

    def test_unknown_element_is_a_400_naming_it(self):
        status, error = self.status_error(self.slice(element="Xx"))
        self.assertEqual(status, 400, error)
        self.assertEqual(error, "Unknown element 'Xx'; available: Nb, Se")

    def test_element_case_is_ignored(self):
        upper = self.slice(element="Se", z=0.5)
        lower = self.slice(element="se", z=0.5)
        self.assertEqual(upper.status_code, 200)
        self.assertEqual(lower.status_code, 200, lower.get_data(as_text=True)[:200])
        self.assertEqual(lower.get_json()["element"], "Se")
        self.assertEqual(lower.get_json()["slabCount"], upper.get_json()["slabCount"])
        self.assertGreater(upper.get_json()["slabCount"], 0)

    def test_a_file_without_parseable_atoms_is_a_400(self):
        folder = self.folder_with("no_atoms", EMPTY_ATOMS)
        for element in (None, "Se"):
            with self.subTest(element=element):
                params = {} if element is None else {"element": element}
                status, error = self.status_error(self.slice(folder, **params))
                self.assertEqual(status, 400, error)
                self.assertIn("no atoms could be parsed", error)


class PcaKdeEmptyVolumeTests(_EdgeCase):
    """/api/pca/kde: a bandwidth or extent that leaves every kernel between the grid
    nodes gave a 200 all-zero volume (mass 0) with no warning."""

    def test_a_volume_that_captures_nothing_is_a_400(self):
        for params in ({"bw": "1e-6"}, {"bwScale": "1e-6"}, {"extent": "1e6"}):
            with self.subTest(**params):
                response = self.client.get(
                    "/api/pca/kde",
                    query_string={"dir": RUN, "referenceNumber": 1, "grid": 16, "projections": "false", **params},
                )
                status, error = self.status_error(response)
                self.assertEqual(status, 400, error)
                self.assertIn("captures less than 1e-6 of the density", error)

    def test_default_volume_is_unaffected(self):
        response = self.client.get(
            "/api/pca/kde", query_string={"dir": RUN, "referenceNumber": 1, "grid": 16, "projections": "false"}
        )
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True)[:200])
        self.assertGreater(response.get_json()["mass"], 0.9)


class PcaKdeVolumeLibraryTests(unittest.TestCase):
    """The engine refuses the empty volume itself (pcaKde.js: the same limit and message)."""

    def test_empty_volume_raises_and_a_normal_one_does_not(self):
        import numpy as np
        from rmc_toolkits.pca_kde import EMPTY_VOLUME_MESSAGE, pca_kde_volume

        rng = np.random.default_rng(4)
        points = rng.normal(size=(400, 3)) * np.array([0.1, 0.07, 0.05])
        self.assertGreater(pca_kde_volume(points, grid=16)["mass"], 0.9)
        for kwargs in ({"bw": 1e-6}, {"bw_scale": 1e-6}, {"extent": 1e6}):
            with self.subTest(**kwargs):
                with self.assertRaises(ValueError) as caught:
                    pca_kde_volume(points, grid=16, **kwargs)
                self.assertEqual(str(caught.exception), EMPTY_VOLUME_MESSAGE)


class FracConversionTests(_EdgeCase):
    """/api/convert/frac: zero parseable atoms wrote a header-only Frac file (200), the
    output could be the source itself, overwrite "false" counted as true, and an output
    path that is a directory was a 500."""

    def post(self, **payload):
        return self.client.post("/api/convert/frac", json=payload)

    def test_zero_parseable_atoms_writes_nothing(self):
        folder = self.folder_with("frac_no_atoms", EMPTY_ATOMS)
        status, error = self.status_error(self.post(path=f"{folder}/run.rmc6f"))
        self.assertEqual(status, 400, error)
        self.assertIn("no atoms could be parsed", error)
        self.assertEqual(list((ROOT / folder).glob("Frac_coord*")), [])

    def test_a_partial_file_reports_what_was_skipped(self):
        text = (self.run_dir / "synthetic.rmc6f").read_text(encoding="utf-8")
        lines = text.splitlines()
        marker = lines.index("Atoms:")
        lines[marker + 3] = lines[marker + 3].rsplit(" ", 1)[0]  # one torn atom line
        folder = self.folder_with("frac_partial", "\n".join(lines) + "\n")
        response = self.post(path=f"{folder}/run.rmc6f")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True)[:200])
        payload = response.get_json()
        self.assertIn("1 of 54 atom lines unparsed", payload["parseWarning"])
        written = (ROOT / folder / "Frac_coord_run.txt").read_text(encoding="utf-8").splitlines()
        self.assertEqual(len(written), 5 + 53)

    def test_clean_file_has_a_null_parse_warning(self):
        folder = self.folder_with("frac_clean", (self.run_dir / "synthetic.rmc6f").read_text())
        response = self.post(path=f"{folder}/run.rmc6f")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True)[:200])
        self.assertIsNone(response.get_json()["parseWarning"])

    def test_the_source_is_never_the_output(self):
        folder = self.folder_with("frac_self", (self.run_dir / "synthetic.rmc6f").read_text())
        before = (ROOT / folder / "run.rmc6f").read_bytes()
        status, error = self.status_error(
            self.post(path=f"{folder}/run.rmc6f", outputPath=f"{folder}/run.rmc6f", overwrite=True)
        )
        self.assertEqual(status, 400, error)
        self.assertIn("would overwrite its source", error)
        self.assertEqual((ROOT / folder / "run.rmc6f").read_bytes(), before)

    def test_overwrite_is_a_strict_boolean(self):
        folder = self.folder_with("frac_flag", (self.run_dir / "synthetic.rmc6f").read_text())
        first = self.post(path=f"{folder}/run.rmc6f")
        self.assertEqual(first.status_code, 200)
        status, _ = self.status_error(self.post(path=f"{folder}/run.rmc6f", overwrite="false"))
        self.assertEqual(status, 409)
        status, error = self.status_error(self.post(path=f"{folder}/run.rmc6f", overwrite="maybe"))
        self.assertEqual(status, 400, error)
        self.assertIn("overwrite must be a boolean", error)
        again = self.post(path=f"{folder}/run.rmc6f", overwrite="true")
        self.assertEqual(again.status_code, 200)

    def test_an_output_path_that_is_a_directory_is_a_400(self):
        folder = self.folder_with("frac_dir", (self.run_dir / "synthetic.rmc6f").read_text())
        (ROOT / folder / "out").mkdir(exist_ok=True)
        status, error = self.status_error(
            self.post(path=f"{folder}/run.rmc6f", outputPath=f"{folder}/out", overwrite=True)
        )
        self.assertEqual(status, 400, error)
        self.assertIn("is a directory", error)


class ScalingEdgeTests(_EdgeCase):
    """/api/scaling/*: a stog.inp naming '.' as its data file and a deeply nested JSON
    body were 500s; inspect "false" entered inspect mode; booleans were coerced."""

    def test_a_stog_inp_whose_data_file_is_a_folder_is_a_404(self):
        inp = self.run_dir / "dot.inp"
        inp.write_text(
            "1\n.\n0.60 30.0\n-9 0.1\n0\nscale.fq\nscale.gr\n25\n1000\nN\n0.05\n0\nN\nY\n1.0\n"
            "scale_ft.sq\nscale_ft.gr\n0.02\nscale_ft_rmc.fq\nscale_ft_rmc.gr\nscale_ft_rmc.dr\n2.48 2.65 3.1\n",
            encoding="utf-8",
        )
        for route in ("/api/scaling/preview", "/api/scaling/run"):
            with self.subTest(route=route):
                status, error = self.status_error(self.client.post(route, json={"path": f"{RUN}/dot.inp"}))
                self.assertEqual(status, 404, error)
                self.assertIn("not a file", error)

    def test_a_deeply_nested_body_is_a_400(self):
        body = "[" * 5000 + "]" * 5000
        for route in ("/api/scaling/preview", "/api/convert/frac"):
            with self.subTest(route=route):
                response = self.client.post(route, data=body, content_type="application/json")
                status, error = self.status_error(response)
                self.assertEqual(status, 400, error)

    def test_inspect_false_is_not_inspect_mode(self):
        response = self.client.post(
            "/api/scaling/preview", json={"path": f"{RUN}/synthetic.rmc6f", "inspect": "false"}
        )
        status, error = self.status_error(response)
        # Not the inspect reply ({"kind": "data", ...} with 200): a real preview,
        # which needs qmin/qmax for a data file.
        self.assertEqual(status, 400, error)
        self.assertIn("qmin", error)
        status, error = self.status_error(self.client.post(
            "/api/scaling/preview", json={"path": f"{RUN}/synthetic.rmc6f", "inspect": "perhaps"}))
        self.assertEqual(status, 400, error)
        self.assertIn("inspect must be a boolean", error)


if __name__ == "__main__":
    unittest.main()
