# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""rmc-triplets checks every destination before computing, and has --version.

An unsupported --plot format printed a traceback after the CSV had been
written; --plot or --dump-angles equal to --output silently replaced the
histogram (exit 0, both paths reported); a destination whose parent is a file
failed after the CSV was written; and --version was an argparse error.
"""

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import mock
import contextlib
import io
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from rmc_toolkits import __version__  # noqa: E402
from rmc_toolkits import triplets_cli  # noqa: E402
from rmc_toolkits.triplets_cli import main  # noqa: E402

TINY = """Supercell dimensions: 1 1 1
Lattice vectors (Ang):
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
Atoms:
1 Nb [1] 0.500000 0.500000 0.500000 1 0 0 0
2 Se [2] 0.600000 0.500000 0.500000 2 0 0 0
3 Se [2] 0.500000 0.600000 0.500000 2 0 0 0
"""


def run(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            code = main([str(arg) for arg in argv])
        except SystemExit as exit_:
            code = exit_.code
    return code, out.getvalue(), err.getvalue()


class TripletsCliDestinationTests(unittest.TestCase):
    def setUp(self):
        self._tmp = TemporaryDirectory()
        self.dir = Path(self._tmp.name)
        self.config = self.dir / "tiny.rmc6f"
        self.config.write_text(TINY, encoding="utf-8")
        self.base = [self.config, "--triplet", "Se", "Nb", "Se", "--bond12", "0.5", "1.5"]

    def tearDown(self):
        self._tmp.cleanup()

    def files(self):
        return sorted(path.name for path in self.dir.rglob("*") if path.is_file() and path != self.config)

    def refused_before_computing(self, argv, fragment):
        with mock.patch.object(triplets_cli, "bond_angles_from_rmc6f") as compute:
            code, _, err = run(self.base + argv)
        self.assertEqual(code, 1, err)
        self.assertIn(fragment, err)
        self.assertNotIn("Traceback", err)
        self.assertEqual(err.count("\n"), 1, err)  # one line
        self.assertFalse(compute.called, "computed before checking the destinations")
        self.assertEqual(self.files(), [])
        return err

    def test_unsupported_plot_format_is_refused_up_front(self):
        err = self.refused_before_computing(
            ["--output", self.dir / "o.csv", "--plot", self.dir / "p.xyz"], "unsupported plot format '.xyz'"
        )
        self.assertIn("png", err)

    def test_destinations_must_be_distinct(self):
        self.refused_before_computing(
            ["--output", self.dir / "o.csv", "--dump-angles", self.dir / "o.csv"], "must be different files"
        )
        self.refused_before_computing(
            ["--output", self.dir / "o.png", "--plot", self.dir / "O.png"], "must be different files"
        )

    def test_the_configuration_is_never_a_destination(self):
        before = self.config.read_bytes()
        self.refused_before_computing(["--output", self.config, "--force"], "is the input configuration")
        self.assertEqual(self.config.read_bytes(), before)

    def test_a_parent_that_is_a_file_is_refused_up_front(self):
        (self.dir / "blocker").write_text("x\n")
        with mock.patch.object(triplets_cli, "bond_angles_from_rmc6f") as compute:
            code, _, err = run(self.base + ["--output", self.dir / "o.csv", "--plot", self.dir / "blocker" / "p.png"])
        self.assertEqual(code, 1, err)
        self.assertIn("not a directory", err)
        self.assertFalse(compute.called)
        self.assertEqual(self.files(), ["blocker"])

    def test_a_directory_destination_is_refused(self):
        (self.dir / "out.csv").mkdir()
        with mock.patch.object(triplets_cli, "bond_angles_from_rmc6f") as compute:
            code, _, err = run(self.base + ["--output", self.dir / "out.csv", "--force"])
        self.assertEqual(code, 1, err)
        self.assertIn("is a directory", err)
        self.assertFalse(compute.called)

    def test_a_failed_plot_leaves_no_partial_outputs(self):
        def boom(*args, **kwargs):
            raise OSError(28, "No space left on device")

        with mock.patch.object(triplets_cli, "write_plot", side_effect=boom):
            code, _, err = run(self.base + ["--output", self.dir / "o.csv", "--plot", self.dir / "p.png",
                                            "--dump-angles", self.dir / "a.txt"])
        self.assertEqual(code, 1, err)
        self.assertIn("No space left", err)
        self.assertEqual(self.files(), [])

    def test_all_outputs_written_on_success(self):
        code, out, err = run(self.base + ["--output", self.dir / "o.csv", "--plot", self.dir / "p.png",
                                          "--dump-angles", self.dir / "a.txt"])
        self.assertEqual(code, 0, err)
        self.assertEqual(self.files(), ["a.txt", "o.csv", "p.png"])
        self.assertTrue((self.dir / "p.png").read_bytes().startswith(b"\x89PNG"))
        self.assertIn("plot:", out)

    def test_version(self):
        code, out, err = run(["--version"])
        self.assertEqual(code, 0)
        self.assertEqual((out + err).strip(), f"rmc-triplets {__version__}")


if __name__ == "__main__":
    unittest.main()
