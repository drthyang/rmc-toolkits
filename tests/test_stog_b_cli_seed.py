# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""``--estimate-rho0`` without any density source seeds the iteration (1.0 audit, stog-b).

The documented composition-only route (``rmc-autoscale --data ... --formula ...
--estimate-rho0``) exited with "number density unknown": _build_config demanded a
rho0 before the estimator (which only uses it as a seed) could run. The CLI now
seeds 0.05 A^-3 like the Auto StoG page and says so.
"""

import contextlib
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import write_stog_xy
from rmc_toolkits.scaling import RHO0_SEED
from rmc_toolkits.scaling_cli import main
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

# The model's true density is far from the 0.05 seed, so recovering it shows
# the self-consistency actually moved from the seed (1.0 review).
RHO0, B2 = 0.08, 0.02


def model():
    q = np.arange(20, 981) * 0.03
    r = np.arange(1, 12001) * 0.005
    g = 0.5 * (1 + np.tanh((r - 2.65) / 0.07)) + 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    sq_true = fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, RHO0), q))
    head = q <= q[0] + 1.0
    _, s_true_0 = np.polyfit(q[head], sq_true[head], 1)
    return q, (sq_true + 9.0) / 10.0, B2 * (1.0 - float(s_true_0))


class EstimateSeedTests(unittest.TestCase):
    def run_cli(self, *extra):
        q, sq, b_sq_avg = model()
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "model.sq"  # no NUMBER_DENSITY header
            write_stog_xy(data, q, sq)
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = main([
                    "--data", str(data), "--qmin", "0.6", "--qmax", "29.4",
                    "--b-avg-sq", str(B2), "--b-sq-avg", str(b_sq_avg), "--r0", "2.5",
                    "--r-fit-min", "1.2", "--rmax", "25", "--nr", "1000",
                    "--out-dir", str(Path(tmp) / "out"), *extra,
                ])
        return code, out.getvalue(), err.getvalue()

    def test_composition_only_estimate_runs_from_the_seed(self):
        self.assertEqual(RHO0_SEED, 0.05)
        code, out, err = self.run_cli("--estimate-rho0")
        self.assertEqual(code, 0, err)
        self.assertIn("seeding the self-consistency at 0.05", err)
        line = [row for row in out.splitlines() if row.startswith("rho0 self-consistency")][0]
        estimate = float(line.split(":")[1].split()[0])
        self.assertLess(abs(estimate - RHO0) / RHO0, 0.05)
        self.assertGreater(abs(estimate - RHO0_SEED) / RHO0_SEED, 0.4)

    def test_without_the_estimate_a_density_is_still_required(self):
        code, _, err = self.run_cli()
        self.assertEqual(code, 2)
        self.assertIn("number density unknown", err)


if __name__ == "__main__":
    unittest.main()
