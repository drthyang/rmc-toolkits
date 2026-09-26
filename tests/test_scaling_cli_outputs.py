# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""rmc-autoscale / /api/scaling/run output integrity.

Before any computation or write, every output target must be a distinct file
(not a directory, not another target, not an input) whose parent can exist;
the family is then written to temporary names and renamed into place only after
every write succeeded. Before this, a stog.inp naming 'sub/rmc.gr' failed after
five files were written, a directory at out/ft.dat failed after seven, and a
stog.inp declaring the FK(Q) name as 'ft.dat' exited 0 with the RMCProfile
input silently replaced by the Fourier-filter correction.
"""

from pathlib import Path
from unittest import mock
import hashlib
import os
import shutil
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

from rmc_toolkits import scaling_cli  # noqa: E402
from rmc_toolkits.parsers import read_stog_xy, write_stog_xy  # noqa: E402
from test_scaling_cli import INP_TEMPLATE, CliSyntheticBase, run_cli  # noqa: E402

CLASSIC_NAMES = (
    "scale.fq", "scale.gr", "scale_ft.sq", "scale_ft.gr",
    "scale_ft_rmc.fq", "scale_ft_rmc.gr", "scale_ft_rmc.dr",
)


def inp_text(data: str = "synth.dat", **renames: str) -> str:
    """The synthetic stog.inp with some declared output names replaced."""
    text = INP_TEMPLATE.format(data=data, rmax=25, nr=1000)
    lines = text.split("\n")
    for old, new in renames.items():
        index = lines.index(old)
        lines[index] = new
    return "\n".join(lines)


def digests(folder: Path) -> dict[str, str]:
    return {
        name: hashlib.sha256((folder / name).read_bytes()).hexdigest() for name in files_under(folder)
    }


def files_under(folder: Path) -> list[str]:
    if not folder.exists():
        return []
    return sorted(str(path.relative_to(folder)) for path in folder.rglob("*") if path.is_file())


class OutputPreflightCliTests(CliSyntheticBase):
    def make_run(self, tmp: str, **renames: str) -> Path:
        run = Path(tmp)
        write_stog_xy(run / "synth.dat", self.q, self.sq_meas, title="synthetic")
        (run / "stog.inp").write_text(inp_text(**renames))
        return run

    def assert_refused_before_computing(self, args, fragment):
        with mock.patch.object(scaling_cli, "autoscale", wraps=scaling_cli.autoscale) as auto, \
                mock.patch.object(scaling_cli, "scale_pipeline", wraps=scaling_cli.scale_pipeline) as manual:
            code, _, err = run_cli(args)
        self.assertEqual(code, 2, msg=err)
        self.assertIn(fragment, err)
        self.assertFalse(auto.called, "autoscale ran before the output check")
        self.assertFalse(manual.called, "scale_pipeline ran before the output check")
        return err

    def test_two_declared_outputs_naming_one_file_are_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            # FK(Q) declared as the fixed-name Fourier-filter correction.
            run = self.make_run(tmp, **{"scale_ft_rmc.fq": "ft.dat"})
            err = self.assert_refused_before_computing([run / "stog.inp"], "same file")
            self.assertIn("ft.dat", err)
            self.assertEqual(files_under(run / "autoscale"), [])

            # Filtered S(Q) declared as the scaled S(Q), even with --force.
            run2 = Path(tmp) / "second"
            run2.mkdir()
            self.make_run(str(run2), **{"scale_ft.sq": "scale.fq"})
            self.assert_refused_before_computing([run2 / "stog.inp", "--force"], "same file")
            self.assertEqual(files_under(run2 / "autoscale"), [])

    def test_names_differing_only_in_case_are_refused(self):
        # One file on a case-insensitive filesystem (the macOS / Windows default).
        with tempfile.TemporaryDirectory() as tmp:
            run = self.make_run(tmp, **{"scale_ft.sq": "SCALE.fq"})
            self.assert_refused_before_computing([run / "stog.inp"], "same file")
            self.assertEqual(files_under(run / "autoscale"), [])

    def test_a_directory_at_an_output_path_is_refused_before_any_write(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.make_run(tmp)
            (run / "autoscale" / "ft.dat").mkdir(parents=True)
            err = self.assert_refused_before_computing([run / "stog.inp", "--force"], "is a directory")
            self.assertIn("ft.dat", err)
            self.assertEqual(files_under(run / "autoscale"), [])

    def test_a_parent_that_is_a_file_is_refused_before_any_write(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.make_run(tmp, **{"scale_ft_rmc.gr": "blocker/rmc.gr"})
            (run / "autoscale").mkdir()
            (run / "autoscale" / "blocker").write_text("not a folder\n")
            err = self.assert_refused_before_computing([run / "stog.inp"], "not a directory")
            self.assertIn("blocker", err)
            self.assertEqual(files_under(run / "autoscale"), ["blocker"])

    def test_declared_subfolders_are_created(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.make_run(tmp, **{"scale_ft_rmc.gr": "sub/rmc.gr"})
            code, _, err = run_cli([run / "stog.inp", "--manual"])
            self.assertEqual(code, 0, msg=err)
            self.assertTrue((run / "autoscale" / "sub" / "rmc.gr").is_file())
            self.assertEqual(len(files_under(run / "autoscale")), 9)

    def test_a_failed_write_leaves_no_family_and_no_temporary_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.make_run(tmp)
            real_write = scaling_cli.write_stog_xy
            calls = []

            def failing_write(path, *args, **kwargs):
                calls.append(path)
                if len(calls) == 5:
                    raise OSError(28, "No space left on device")
                return real_write(path, *args, **kwargs)

            with mock.patch.object(scaling_cli, "write_stog_xy", side_effect=failing_write):
                code, _, err = run_cli([run / "stog.inp", "--manual"])
            self.assertEqual(code, 2, msg=err)
            self.assertIn("No space left", err)
            self.assertEqual(files_under(run / "autoscale"), [])

    def test_a_failed_forced_rewrite_keeps_the_previous_family_intact(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.make_run(tmp)
            code, _, err = run_cli([run / "stog.inp", "--manual"])
            self.assertEqual(code, 0, msg=err)
            before = digests(run / "autoscale")
            self.assertEqual(len(before), 9)

            def failing_dump(*args, **kwargs):
                raise OSError(5, "Input/output error")

            with mock.patch.object(scaling_cli.json, "dump", side_effect=failing_dump):
                # A different scale, so a partial rewrite would change the bytes.
                code, _, err = run_cli([run / "stog.inp", "--force", "--scale", "9.5"])
            self.assertEqual(code, 2, msg=err)
            self.assertEqual(digests(run / "autoscale"), before)

    def test_successful_run_leaves_only_the_family(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = self.make_run(tmp)
            code, _, err = run_cli([run / "stog.inp", "--manual"])
            self.assertEqual(code, 0, msg=err)
            expected = sorted([*CLASSIC_NAMES, "ft.dat", "stog_provenance.json"])
            self.assertEqual(files_under(run / "autoscale"), expected)
            scaled = read_stog_xy(run / "autoscale" / "scale.fq")
            self.assertTrue(np.isfinite(scaled).all())


class OutputPreflightApiTests(CliSyntheticBase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        import app as backend_app

        backend_app.app.config.update(TESTING=True)
        cls.backend = backend_app
        cls.client = backend_app.app.test_client()
        cls.root = ROOT / "results" / "scaling_output_preflight"

    def setUp(self):
        shutil.rmtree(self.root, ignore_errors=True)
        self.root.mkdir(parents=True)
        self.backend._SCALING_CACHE.clear()

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def make_run(self, **renames: str) -> Path:
        write_stog_xy(self.root / "synth.dat", self.q, self.sq_meas, title="synthetic")
        (self.root / "stog.inp").write_text(inp_text(**renames))
        return self.root

    def post(self, **payload):
        return self.client.post("/api/scaling/run", json=payload)

    def test_colliding_declared_names_are_a_400_before_computing(self):
        run = self.make_run(**{"scale_ft_rmc.fq": "ft.dat"})
        with mock.patch.object(self.backend, "_compute_scaling") as compute:
            response = self.post(path=str((run / "stog.inp").relative_to(ROOT)), mode="manual")
        self.assertEqual(response.status_code, 400, response.get_data(as_text=True))
        self.assertIn("same file", response.get_json()["error"])
        self.assertFalse(compute.called)
        self.assertEqual(files_under(run / "autoscale"), [])

    def test_out_dir_naming_a_file_is_a_400_not_a_500(self):
        run = self.make_run()
        (run / "not_a_folder").write_text("x\n")
        with mock.patch.object(self.backend, "_compute_scaling") as compute:
            response = self.post(
                path=str((run / "stog.inp").relative_to(ROOT)),
                mode="manual",
                outDir=str((run / "not_a_folder").relative_to(ROOT)),
            )
        self.assertEqual(response.status_code, 400, response.get_data(as_text=True))
        self.assertIn("not a directory", response.get_json()["error"])
        self.assertFalse(compute.called)

    def test_a_failed_write_leaves_no_family(self):
        run = self.make_run()
        real_write = scaling_cli.write_stog_xy
        calls = []

        def failing_write(path, *args, **kwargs):
            calls.append(path)
            if len(calls) == 3:
                raise OSError(28, "No space left on device")
            return real_write(path, *args, **kwargs)

        with mock.patch.object(scaling_cli, "write_stog_xy", side_effect=failing_write):
            response = self.post(path=str((run / "stog.inp").relative_to(ROOT)), mode="manual")
        self.assertGreaterEqual(response.status_code, 400)
        self.assertEqual(files_under(run / "autoscale"), [])

    def test_success_writes_the_whole_family(self):
        run = self.make_run()
        response = self.post(path=str((run / "stog.inp").relative_to(ROOT)), mode="manual")
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        self.assertEqual(len(files_under(run / "autoscale")), 9)


if __name__ == "__main__":
    unittest.main()
