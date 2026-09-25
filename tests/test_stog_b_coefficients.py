# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""<b>^2 and <b^2> come from one consistent source (1.0 audit, stog-b group).

The CLI (and the API and page) took <b^2> from --formula (Sears neutron values)
whenever --b-sq-avg was absent, even when <b>^2 came from a stog.inp or
--b-avg-sq in another unit system (1.0 for normalized x-ray data). The mixed
ratio fabricated an S(0) target (+0.55 for FeCoSn: impossible, <b^2> >= <b>^2)
that fed the low-Q correction, the FZ amplitude, the concordance verdict and
estimate_rho0 (a 60 %-low density "converged"). Now the formula's <b^2> is used
only with a <b>^2 that agrees with the formula's, S(0) > 0 is rejected by
ScalingConfig, and the CLI states the values in effect.
"""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import write_stog_xy
from rmc_toolkits.scaling import ScalingConfig
from rmc_toolkits.scaling_cli import main, resolve_coefficients
from rmc_toolkits.scattering import faber_ziman
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

RHO0 = 0.05


def model_sq():
    q = np.arange(20, 981) * 0.03
    r = np.arange(1, 12001) * 0.005
    g = 0.5 * (1.0 + np.tanh((r - 2.65) / 0.07)) + 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    return q, (fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, RHO0), q)) + 4.0) / 5.0


class CauchySchwarzTests(unittest.TestCase):
    def base(self, **overrides):
        values = dict(qmin=0.5, qmax=26.0, rho0=RHO0, b_avg_sq=1.0)
        values.update(overrides)
        return ScalingConfig(**values)

    def test_positive_s0_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Cauchy-Schwarz"):
            self.base(b_sq_avg=0.447511)  # FeCoSn neutron <b^2> against x-ray <b>^2 = 1

    def test_single_element_and_polyatomic_pass(self):
        self.assertEqual(self.base(b_sq_avg=1.0).effective_s0_target, 0.0)
        self.assertAlmostEqual(self.base(b_sq_avg=1.10426).effective_s0_target, -0.10426)

    def test_non_positive_or_nan_b_sq_avg_is_rejected(self):
        for value in (0.0, -1.0, float("nan")):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "b_sq_avg"):
                self.base(b_sq_avg=value)


class ResolveCoefficientsTests(unittest.TestCase):
    def test_formula_alone_gives_its_pair(self):
        fz = faber_ziman("Mn3Sn")
        out = resolve_coefficients(b_avg_sq=None, b_avg_sq_source=None, b_sq_avg=None, formula="Mn3Sn")
        self.assertEqual(out["b_avg_sq"], fz.b_avg_sq_barn)
        self.assertAlmostEqual(out["b_sq_avg"], fz.b_sq_avg_barn, places=15)
        self.assertEqual(out["warnings"], [])

    def test_agreeing_b_avg_sq_keeps_the_formula_ratio(self):
        fz = faber_ziman("Mn3Sn")
        out = resolve_coefficients(
            b_avg_sq=0.015407, b_avg_sq_source="stog.inp", b_sq_avg=None, formula="Mn3Sn",
        )
        self.assertAlmostEqual(
            out["b_sq_avg"] / out["b_avg_sq"], fz.b_sq_avg_barn / fz.b_avg_sq_barn, places=12,
        )

    def test_another_source_never_takes_the_formula_b_sq_avg(self):
        out = resolve_coefficients(
            b_avg_sq=1.0, b_avg_sq_source="stog.inp", b_sq_avg=None, formula="FeCoSn",
        )
        self.assertEqual(out["b_avg_sq"], 1.0)
        self.assertIsNone(out["b_sq_avg"])
        self.assertIn("NOT the formula's <b^2>", out["warnings"][0])
        explicit = resolve_coefficients(
            b_avg_sq=1.0, b_avg_sq_source="--b-avg-sq", b_sq_avg=1.10426, formula="FeCoSn",
        )
        self.assertEqual(explicit["b_sq_avg"], 1.10426)


class CliCoefficientTests(unittest.TestCase):
    def run_cli(self, *extra):
        q, sq = model_sq()
        with tempfile.TemporaryDirectory() as tmp:
            data = Path(tmp) / "model.sq"
            write_stog_xy(data, q, sq)
            out, err = io.StringIO(), io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = main([
                    "--data", str(data), "--qmin", "0.6", "--qmax", "29.4",
                    "--rmax", "25", "--nr", "1000", "--r0", "2.5",
                    "--out-dir", str(Path(tmp) / "out"), *extra,
                ])
            provenance = None
            path = Path(tmp) / "out" / "model_provenance.json"
            if path.exists():
                provenance = json.loads(path.read_text())
        return code, out.getvalue(), err.getvalue(), provenance

    def test_xray_b_avg_sq_with_a_formula_does_not_mix_sources(self):
        code, out, err, provenance = self.run_cli(
            "--rho0", str(RHO0), "--b-avg-sq", "1.0", "--formula", "FeCoSn",
        )
        self.assertEqual(code, 0, err)
        self.assertIn("NOT the formula's <b^2> = 0.447511 barn", err)
        self.assertIn("<b^2> = not set", out)
        self.assertIsNone(provenance["provenance"]["config"]["b_sq_avg"])

    def test_estimate_rho0_needs_a_consistent_b_sq_avg(self):
        code, _, err, _ = self.run_cli(
            "--rho0", str(RHO0), "--b-avg-sq", "1.0", "--formula", "FeCoSn", "--estimate-rho0",
        )
        self.assertEqual(code, 2)
        self.assertIn("--estimate-rho0 requires <b^2>", err)

    def test_explicit_impossible_pair_is_refused(self):
        code, _, err, _ = self.run_cli(
            "--rho0", str(RHO0), "--b-avg-sq", "1.0", "--b-sq-avg", "0.5",
        )
        self.assertEqual(code, 2)
        self.assertIn("Cauchy-Schwarz", err)

    def test_report_states_the_coefficients_in_effect(self):
        code, out, err, _ = self.run_cli("--rho0", str(RHO0), "--formula", "SrTiO3")
        self.assertEqual(code, 0, err)
        fz = faber_ziman("SrTiO3")
        self.assertIn(f"<b>^2 = {fz.b_avg_sq_barn:.6g} barn", out)
        self.assertIn(f"S(0) target = {1 - fz.b_sq_avg_barn / fz.b_avg_sq_barn:.4g}", out)


if __name__ == "__main__":
    unittest.main()
