# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The committed PCA-KDE golden must still be what pca_kde.py computes.

``web_app/frontend/src/__tests__/fixtures/pca_kde_fixture.json`` holds the fit
subsample draw (mulberry32 stream, drawn indices) and two volumes above
``MAX_PCA_FIT_POINTS`` that ``workers/pcaKde.js`` is held to
(``src/__tests__/pcaKdeParity.test.js``). This test regenerates it in memory and
fails when the engine has drifted from the file, i.e. when
``PYTHONPATH=. python tests/generate_pca_kde_fixture.py`` needs to be re-run (and
the JS port re-checked against the new numbers). The draw is compared exactly;
the volumes with tests/golden_compare.py, the rule the JS test applies.
"""

from __future__ import annotations

import json
import unittest

import numpy as np

try:  # `python -m unittest discover -s tests` puts tests/ itself on sys.path
    from generate_pca_kde_fixture import OUT, build_fixture
    from golden_compare import assert_same
except ModuleNotFoundError:  # `python -m unittest tests.test_pca_kde_fixture`
    from tests.generate_pca_kde_fixture import OUT, build_fixture
    from tests.golden_compare import assert_same

from rmc_toolkits.pca_kde import MAX_PCA_FIT_POINTS, fit_subsample_indices, mulberry32_uint32

STALE = "pca_kde_fixture.json is stale: run `PYTHONPATH=. python tests/generate_pca_kde_fixture.py`"


class PcaKdeFixtureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.committed = json.loads(OUT.read_text(encoding="utf-8"))
        cls.regenerated = json.loads(json.dumps(build_fixture(), allow_nan=False))

    def test_the_draw_is_unchanged(self):
        for field in ("maxFitPoints", "stream", "draws", "synthetic"):
            self.assertEqual(self.regenerated[field], self.committed[field], f"{STALE} ({field})")

    def test_volumes(self):
        got_cases = self.regenerated["volumes"]
        want_cases = self.committed["volumes"]
        self.assertEqual([case["name"] for case in got_cases], [case["name"] for case in want_cases], STALE)
        for got, want in zip(got_cases, want_cases):
            with self.subTest(case=want["name"]):
                for field in ("source", "options", "fitChecksum"):
                    self.assertEqual(got[field], want[field], f"{STALE} ({field})")
                assert_same(self, got["cloud"], want["cloud"], "cloud", message=STALE)
                got_walls = got["expected"].pop("projections")
                want_walls = want["expected"].pop("projections")
                # The walls' `axes` are index pairs, not an eigenbasis: compared apart.
                assert_same(self, got["expected"], want["expected"], "volume",
                            bases={"axes": want["expected"]["eigenvalues"]}, message=STALE)
                assert_same(self, got_walls, want_walls, "projections", message=STALE)

    def test_every_volume_is_above_the_cap(self):
        for case in self.committed["volumes"]:
            self.assertGreater(case["cloud"]["count"], MAX_PCA_FIT_POINTS)
            self.assertEqual(case["expected"]["count"], case["cloud"]["count"])
            self.assertEqual(case["expected"]["fitCount"], MAX_PCA_FIT_POINTS)


class FitSubsampleDrawTests(unittest.TestCase):
    """Properties of the draw itself, beside the golden."""

    def test_mulberry32_reference_values(self):
        # mulberry32(0) as the JS reference implementation prints it (x 2^32).
        self.assertEqual(mulberry32_uint32(0, 3).tolist(), [1144304738, 1416247, 958946056])
        # seed >>> 0: -1 is 2^32 - 1, and 2^32 + 5 is 5.
        self.assertEqual(mulberry32_uint32(-1, 4).tolist(), mulberry32_uint32(2**32 - 1, 4).tolist())
        self.assertEqual(mulberry32_uint32(2**32 + 5, 4).tolist(), mulberry32_uint32(5, 4).tolist())
        self.assertEqual(mulberry32_uint32(7, 0).size, 0)

    def test_seed_must_be_an_integer(self):
        for seed in (0.5, "0", None, True):
            with self.subTest(seed=seed):
                with self.assertRaises(ValueError):
                    mulberry32_uint32(seed, 3)

    def test_indices_are_distinct_ascending_and_in_range(self):
        for count, limit, seed in ((21952, 20000, 0), (100, 3, 11), (5, 4, 2**31)):
            with self.subTest(count=count, limit=limit, seed=seed):
                indices = fit_subsample_indices(count, limit, seed)
                self.assertEqual(indices.size, limit)
                self.assertTrue(np.all(np.diff(indices) > 0))
                self.assertGreaterEqual(int(indices[0]), 0)
                self.assertLess(int(indices[-1]), count)

    def test_at_or_below_the_cap_nothing_is_drawn(self):
        for count in (1, 19999, 20000):
            self.assertEqual(fit_subsample_indices(count, 20000).tolist(), list(range(count)))

    def test_the_sparse_shuffle_is_a_full_fisher_yates(self):
        # The dict keeps only displaced entries; a full index array gives the same draw.
        count, limit, seed = 300, 120, 9
        uniforms = mulberry32_uint32(seed, limit).astype(float) / 4294967296.0
        indices = list(range(count))
        for index in range(limit):
            swap = index + int(np.floor(uniforms[index] * float(count - index)))
            indices[index], indices[swap] = indices[swap], indices[index]
        self.assertEqual(fit_subsample_indices(count, limit, seed).tolist(), sorted(indices[:limit]))


if __name__ == "__main__":
    unittest.main()
