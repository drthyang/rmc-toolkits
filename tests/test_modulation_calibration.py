# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The modulation engine's false-positive rate under independent cells (F1, plan item 4.4).

Every family of ``tests/modulation_configs.py`` without order between cells --
thermal motion only, and configurations whose sites share a cell-level
variable (Jahn-Teller ligands, octahedral cluster vertices, breathing
triangles, rigid dimers) -- is run over many fresh seeds. Each of the three
tests (the zone-boundary points, the general modes, the diffuse test) is one
Sidak family per configuration, so its family-wise sigma clears a level s with
probability alpha(s) = 1 - Phi(s) (or less: the zone-boundary and general
tests read the strongest points, the diffuse test is a rank test). The counts
must lie within sampling error of that: at most the binomial 99.9 % quantile
of n alpha(s), for s = 3 (the bar), 2 and 1. Before F1 the null of independent
sites found order in 9-31 % of the Jahn-Teller and octahedral seeds.

By default a fast subset runs (``FAST_SEEDS`` per family: no configuration may
reach a verdict of order). ``RMC_TOOLKITS_FULL_CALIBRATION=1`` runs
``FULL_SEEDS`` per family and checks the counts at every level (about a
minute).
"""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.special import ndtr
from scipy.stats import binom

sys.path.insert(0, str(Path(__file__).resolve().parent))

import modulation_configs as configs  # noqa: E402
from rmc_toolkits.modulations import displacement_modulations  # noqa: E402

FAST_SEEDS = range(910001, 910016)
FULL_SEEDS = range(920001, 921001)
FULL = os.environ.get("RMC_TOOLKITS_FULL_CALIBRATION") == "1"
LEVELS = (3.0, 2.0, 1.0)
FINDINGS = ("cellDoubling", "generalOrder", "diffuseOrder")


def _maxima(result: dict) -> dict:
    """Each test's largest family-wise sigma over the configuration (-inf when none is defined)."""
    tested = [entry for entry in result["elements"] if entry["sites"]]

    def largest(values):
        values = [value for value in values if value is not None]
        return max(values) if values else -np.inf

    return {
        "cellDoubling": largest(point["sigma"] for entry in tested for point in entry["zoneBoundary"]),
        "generalOrder": largest(mode["sigma"] for entry in tested for mode in entry["strongest"] if mode["kind"] == "general"),
        "diffuseOrder": largest(entry["diffuse"]["sigma"] for entry in tested),
    }


def run_family(make, seeds) -> tuple[dict, dict]:
    """(sigma maxima per test, verdicts of order per test) over the seeds."""
    maxima = {name: [] for name in FINDINGS}
    orders = {name: 0 for name in FINDINGS}
    for seed in seeds:
        result = displacement_modulations(configs.site_displacements(make(seed)))
        for name, value in _maxima(result).items():
            maxima[name].append(value)
        for name in FINDINGS:
            orders[name] += any(entry["verdict"][name]["answer"] is True for entry in result["elements"] if entry["sites"])
    return {name: np.asarray(values) for name, values in maxima.items()}, orders


class FastCalibrationTests(unittest.TestCase):
    def test_no_configuration_without_order_between_cells_reads_as_order(self):
        for family, make in configs.NULL_FAMILIES.items():
            _, orders = run_family(make, FAST_SEEDS)
            self.assertEqual(orders, {name: 0 for name in FINDINGS}, family)


@unittest.skipUnless(FULL, "set RMC_TOOLKITS_FULL_CALIBRATION=1 for the full calibration (1000 seeds per family)")
class FullCalibrationTests(unittest.TestCase):
    def test_every_level_is_within_sampling_error(self):
        n = len(FULL_SEEDS)
        for family, make in configs.NULL_FAMILIES.items():
            maxima, orders = run_family(make, FULL_SEEDS)
            for name in FINDINGS:
                for level in LEVELS:
                    alpha = float(ndtr(-level))
                    count = int(np.sum(maxima[name] >= level))
                    bound = int(binom.ppf(0.999, n, alpha))
                    self.assertLessEqual(count, bound, f"{family} {name}: {count} of {n} at {level} sigma (alpha {alpha:.4g})")
                # A verdict of order also needs the effect size: never more often than its sigma alone.
                self.assertLessEqual(orders[name], int(np.sum(maxima[name] >= 3.0)), f"{family} {name}")


if __name__ == "__main__":
    unittest.main()
