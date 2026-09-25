# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""rho0 self-consistency: report WHY it stopped (1.0 review follow-up, stog-a group).

estimate_rho0 returns ``stopped`` when autoscale cannot fit a trial density, but the
CLI (like the page worker) still blamed the two amplitude criteria for disagreeing
"at every density" -- a different failure with a different remedy.
"""

import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from rmc_toolkits.parsers import write_stog_xy
from rmc_toolkits.scaling_cli import main
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

RHO0 = 0.05


def shell_sq():
    r = np.arange(1, 16001) * 0.005
    g = 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.1) ** 2) + 0.5 * (1.0 + np.tanh((r - 3.4) / 0.08))
    q = np.arange(50, 2601) * 0.01
    return q, fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, RHO0), q))


def estimate(stopped):
    return {
        "rho0": 0.063, "converged": False, "iterations": 2, "concordance": 9.87,
        "a_density": 1.0, "a_fz": 9.87, "extrapolated": False,
        "history": [[0.063, 1.0, 9.87, 9.87]], "stopped": stopped,
    }


class Rho0StoppedMessageTests(unittest.TestCase):
    def run_cli(self, stopped):
        q, sq = shell_sq()
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "shell.sq"
            write_stog_xy(data, q, sq)
            out, err = io.StringIO(), io.StringIO()
            with mock.patch(
                "rmc_toolkits.scaling_cli.estimate_rho0", return_value=estimate(stopped)
            ), contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = main([
                    "--data", str(data), "--qmin", "0.5", "--qmax", "26", "--rho0", str(RHO0),
                    "--b-avg-sq", "1.0", "--b-sq-avg", "2.0", "--estimate-rho0",
                    "--out-dir", str(Path(tmp) / "out"),
                ])
        return code, err.getvalue()

    def test_cli_reports_the_trial_density_that_could_not_be_fitted(self):
        reason = "autoscale failed at rho0 = 0.66: autoscale: could not locate the first shell"
        code, err = self.run_cli(reason)
        self.assertEqual(code, 2)
        self.assertIn(reason, err)
        self.assertNotIn("disagree at every density", err)

    def test_cli_keeps_the_discordance_message_otherwise(self):
        code, err = self.run_cli(None)
        self.assertEqual(code, 2)
        self.assertIn("disagree at every density", err)


if __name__ == "__main__":
    unittest.main()
