# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Freshness of the backend's parsed-file caches.

The KDE-slice, PCA/orientation, triplets and scaling caches must notice any
rewrite of their source file -- including one that lands in the same
whole-second mtime, which is what sshfs/SFTP mounts, ``scp -p`` and rsync from
a coarse filesystem produce while an RMCProfile run is being watched -- and must
never keep a parse of a file that changed while it was being read.

Runs on synthetic fixtures under results/ (data root = repo root, as in
tests/test_backend_api.py), so it needs no sample data.
"""

from pathlib import Path
import json
import os
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))

os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "rmc_toolkits_cache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import numpy as np  # noqa: E402

import app as backend_app  # noqa: E402


# A whole-second mtime, as a coarse-resolution mount reports it.
WHOLE_SECOND_NS = 1_790_213_288 * 1_000_000_000


def rmc6f_lines(supercell=(6, 6, 6), seed=5) -> list[str]:
    """Two-site cubic run (Nb at the origin, Se at the body centre), 8 A cell."""
    rng = np.random.default_rng(seed)
    n1, n2, n3 = supercell
    lines = [
        f"Supercell dimensions {n1} {n2} {n3}",
        "Lattice vectors (Ang):",
        f"{8.0 * n1} 0.0 0.0",
        f"0.0 {8.0 * n2} 0.0",
        f"0.0 0.0 {8.0 * n3}",
        "Atoms:",
    ]
    atom = 0
    for ix in range(n1):
        for iy in range(n2):
            for iz in range(n3):
                for reference, element, basis in ((1, "Nb", 0.0), (2, "Se", 0.5)):
                    atom += 1
                    coord = (np.array([ix, iy, iz]) + basis) / np.asarray(supercell)
                    coord = (coord + rng.normal(size=3) * 0.003) % 1.0
                    lines.append(
                        f"{atom} {element} [{reference}] {coord[0]:.10f} {coord[1]:.10f} "
                        f"{coord[2]:.10f} {reference} {ix} {iy} {iz}"
                    )
    return lines


FULL_LINES = rmc6f_lines()
FULL_TEXT = "\n".join(FULL_LINES) + "\n"
FULL_ATOMS = len(FULL_LINES) - 6
# A half-written file: the header plus the first 55% of the atom lines.
PARTIAL_TEXT = "\n".join(FULL_LINES[: 6 + int(0.55 * FULL_ATOMS)]) + "\n"


def write_at(path: Path, text: str, mtime_ns: int = WHOLE_SECOND_NS) -> None:
    path.write_text(text, encoding="utf-8")
    os.utime(path, ns=(mtime_ns, mtime_ns))


class _CacheCase(unittest.TestCase):
    RUN = "results/backend_cache_test"
    PRISTINE = "results/backend_cache_test_pristine"

    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

    def setUp(self):
        self.run_dir = ROOT / self.RUN
        self.pristine_dir = ROOT / self.PRISTINE
        for directory in (self.run_dir, self.pristine_dir):
            shutil.rmtree(directory, ignore_errors=True)
            directory.mkdir(parents=True)
        self.path = self.run_dir / "run.rmc6f"
        # Same content in a folder the caches have never seen: the ground truth.
        write_at(self.pristine_dir / "run.rmc6f", FULL_TEXT)

    def tearDown(self):
        shutil.rmtree(self.run_dir, ignore_errors=True)
        shutil.rmtree(self.pristine_dir, ignore_errors=True)

    def get_json(self, url, directory, **params):
        response = self.client.get(url, query_string={"dir": directory, **params})
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True)[:300])
        return response.get_json()


class SameSecondRewriteTests(_CacheCase):
    KDE = {"element": "Nb", "z": 0.0, "dz": 0.1, "grid": 24, "levels": 0}
    TRIPLETS = {"end1": "Se", "apex": "Nb", "end2": "Se", "r12Min": 6.0, "r12Max": 7.5}

    def snapshot(self, directory):
        return {
            "sites": self.get_json("/api/pca/sites", directory)["totalAtoms"],
            "kde": self.get_json("/api/kde/slice", directory, **self.KDE)["slabCount"],
            "triplets": self.get_json("/api/triplets", directory, **self.TRIPLETS)["angleCount"],
            "orientation": self.get_json(
                "/api/pca/orientation", directory, referenceNumber=1, frequency=3, geometry="false"
            )["totalPoints"],
        }

    def test_completed_file_with_the_same_whole_second_mtime_is_reparsed(self):
        write_at(self.path, PARTIAL_TEXT)
        partial = self.snapshot(self.RUN)
        self.assertLess(partial["sites"], FULL_ATOMS)

        # The writer finishes within the same (whole) second.
        write_at(self.path, FULL_TEXT)
        self.assertEqual(os.stat(self.path).st_mtime_ns, WHOLE_SECOND_NS)

        truth = self.snapshot(self.PRISTINE)
        self.assertEqual(truth["sites"], FULL_ATOMS)
        self.assertEqual(self.snapshot(self.RUN), truth)


class RaceDuringReadTests(_CacheCase):
    def setUp(self):
        super().setUp()
        self.real_loader = backend_app.load_site_displacements
        self.calls = 0

    def tearDown(self):
        backend_app.load_site_displacements = self.real_loader
        super().tearDown()

    def test_a_file_that_changes_during_the_read_is_reread_not_cached(self):
        write_at(self.path, PARTIAL_TEXT)

        def loader(path):
            self.calls += 1
            sites = self.real_loader(path)
            if self.calls == 1:
                # The writer completes the file while the first parse runs.
                write_at(self.path, FULL_TEXT)
            return sites

        backend_app.load_site_displacements = loader
        first = self.get_json("/api/pca/sites", self.RUN)
        self.assertEqual(first["totalAtoms"], FULL_ATOMS)
        self.assertEqual(self.calls, 2)

        # The stable re-read was cached; the torn one was not.
        again = self.get_json("/api/pca/sites", self.RUN)
        self.assertEqual(again["totalAtoms"], FULL_ATOMS)
        self.assertEqual(self.calls, 2)

    def test_a_file_that_keeps_changing_is_a_409_not_a_torn_result(self):
        write_at(self.path, PARTIAL_TEXT)
        texts = [FULL_TEXT, PARTIAL_TEXT]

        def loader(path):
            self.calls += 1
            sites = self.real_loader(path)
            write_at(self.path, texts[self.calls % 2])
            return sites

        backend_app.load_site_displacements = loader
        response = self.client.get("/api/pca/sites", query_string={"dir": self.RUN})
        self.assertEqual(response.status_code, 409)
        self.assertIn("changed while it was being read", response.get_json()["error"])


class ScalingCacheTests(_CacheCase):
    def test_rewritten_data_file_with_the_same_mtime_is_reparsed(self):
        from rmc_toolkits.parsers import write_stog_xy

        q = np.arange(20, 981) * 0.03
        sq = 1.0 + 0.3 * np.exp(-q) * np.sin(2.6 * q)

        def write_data(path, rows, gain):
            write_stog_xy(path, q[:rows], gain * sq[:rows], title="synthetic")
            os.utime(path, ns=(WHOLE_SECOND_NS, WHOLE_SECOND_NS))

        def preview(directory):
            response = self.client.post(
                "/api/scaling/preview",
                json={
                    "path": f"{directory}/synth.dat",
                    "qmin": 0.6, "qmax": 26.0, "rho0": 0.05, "bAvgSq": 0.02, "r0": 2.65,
                    "mode": "manual", "a": 1.0, "b": 0.0, "enforce": False, "rmax": 20.0, "nr": 1000,
                },
            )
            self.assertEqual(response.status_code, 200, response.get_data(as_text=True)[:300])
            return response.get_json()["series"]["sqScaled"]

        # Half-written file, then the finished (different) data, same whole second.
        write_data(self.run_dir / "synth.dat", 900, 1.0)
        first = preview(self.RUN)
        write_data(self.run_dir / "synth.dat", q.size, 1.1)
        write_data(self.pristine_dir / "synth.dat", q.size, 1.1)

        truth = preview(self.PRISTINE)
        # Plain comparisons: a list diff of ~900 floats takes difflib seconds.
        self.assertTrue(first != truth, "fixture: the rewrite must change the data")
        self.assertTrue(preview(self.RUN) == truth, "stale cached scaling result served")


class FileSignatureTests(unittest.TestCase):
    def test_signature_sees_same_second_rewrites(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "run.rmc6f"
            write_at(path, PARTIAL_TEXT)
            before = backend_app._file_signature(path)
            write_at(path, FULL_TEXT)
            self.assertNotEqual(backend_app._file_signature(path), before)
            json.dumps(before)  # plain ints: usable in cache keys and payloads


if __name__ == "__main__":
    unittest.main()
