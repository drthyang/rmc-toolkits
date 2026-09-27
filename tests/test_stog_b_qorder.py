# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Q-grid order (0.6.0 audit, stog-b group).

The trapezoid sine transform assumes a strictly increasing grid; with descending Q
every panel width is negative, G_PDF(r) comes out negated, Q0/Qmax swap in the
low-Q correction and the Lorch window, and autoscale "converged" to a negative
scale (a = -2.29 on Mn3Sn 59438) with no warning. Descending files are realistic
(Q = 2 pi/d keeps a d-ordered row order). Now crop_sq sorts to ascending Q and
rejects duplicate or overlapping (two concatenated banks) Q with a clear error,
the public transforms raise on non-increasing grids, and a fit with a <= 0 is
reported as failed, never converged.
"""

import contextlib
import io
from pathlib import Path
import tempfile
import unittest

import numpy as np

from rmc_toolkits.parsers import write_stog_xy
from rmc_toolkits.scaling import ScalingConfig, autoscale, crop_sq
from rmc_toolkits.scaling_cli import main
from rmc_toolkits.transforms import (
    fourier_filter,
    fq_to_gpdf,
    fq_to_sq,
    g_to_gpdf,
    gpdf_to_fq,
    low_q_correction_basis,
    sine_transform,
    sq_to_fq,
)

RHO0, B2 = 0.05, 0.02


def synthetic_g(r):
    """The shared synthetic model (tests/generate_autoscale_fixture.py)."""
    onset = 0.5 * (1.0 + np.tanh((r - 2.65) / 0.07))
    peak = 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    return onset + peak


def model_sq():
    q = np.arange(20, 981) * 0.03
    r = np.arange(1, 12001) * 0.005
    sq_true = fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, synthetic_g(r), RHO0), q))
    return q, sq_true


def config(**overrides):
    base = dict(
        qmin=0.6, qmax=30.0, rho0=RHO0, b_avg_sq=B2, r0=2.5, rmax=25.0, nr=1000,
    )
    base.update(overrides)
    return ScalingConfig(**base)


class CropOrderTests(unittest.TestCase):
    def test_descending_input_is_sorted_with_its_sigma(self):
        q, sq = model_sq()
        sigma = 1e-3 * (1.0 + q)
        asc = crop_sq(q, sq, config(), sigma)
        desc = crop_sq(q[::-1], sq[::-1], config(), sigma[::-1])
        for left, right in zip(asc, desc):
            np.testing.assert_array_equal(left, right)

    def test_non_overlapping_segments_in_any_order_are_sorted(self):
        q, sq = model_sq()
        half = q.size // 2
        order = np.r_[np.arange(half, q.size), np.arange(half)]  # high bank first
        cropped_q, cropped_sq, _ = crop_sq(q[order], sq[order], config())
        np.testing.assert_array_equal(cropped_q, q)
        np.testing.assert_array_equal(cropped_sq, sq)

    def test_overlapping_banks_are_rejected(self):
        q, sq = model_sq()
        bank1 = q <= 16.0
        bank2 = q >= 14.0  # overlaps bank 1 on [14, 16]
        q2 = np.r_[q[bank1], q[bank2] + 0.005]
        sq2 = np.r_[sq[bank1], sq[bank2]]
        with self.assertRaisesRegex(ValueError, "overlap"):
            crop_sq(q2, sq2, config())

    def test_duplicate_q_is_rejected(self):
        q, sq = model_sq()
        q2 = np.r_[q[:100], q[99:]]
        sq2 = np.r_[sq[:100], sq[99:]]
        with self.assertRaisesRegex(ValueError, "duplicate"):
            crop_sq(q2, sq2, config())


class DescendingAutoscaleTests(unittest.TestCase):
    def test_descending_file_gives_the_ascending_result(self):
        q, sq_true = model_sq()
        sq_meas = (sq_true + 9.0) / 10.0
        asc = autoscale(q, sq_meas, config())
        desc = autoscale(q[::-1], sq_meas[::-1], config())
        self.assertGreater(asc.a, 0)
        self.assertEqual(desc.a, asc.a)
        self.assertEqual(desc.b, asc.b)
        self.assertTrue(desc.converged)
        np.testing.assert_array_equal(desc.q, asc.q)

    def test_a_non_positive_scale_is_a_failed_fit(self):
        q, sq_true = model_sq()
        inverted = 2.0 - sq_true  # a sign-inverted measurement: the fit needs a < 0
        result = autoscale(q, inverted, config())
        self.assertLessEqual(result.a, 0)
        self.assertFalse(result.converged)
        self.assertIn("a <= 0", result.provenance["fit_failure"])

    def test_cli_descending_file_matches_and_negative_scale_is_refused(self):
        q, sq_true = model_sq()
        sq_meas = (sq_true + 9.0) / 10.0
        runs = {}
        with tempfile.TemporaryDirectory() as tmp:
            for name, (qq, ss) in {
                "asc": (q, sq_meas), "desc": (q[::-1], sq_meas[::-1]),
                "inverted": (q, 2.0 - sq_true),
            }.items():
                data = Path(tmp) / f"{name}.sq"
                write_stog_xy(data, qq, ss)
                out, err = io.StringIO(), io.StringIO()
                with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                    code = main([
                        "--data", str(data), "--qmin", "0.6", "--qmax", "30",
                        "--rho0", str(RHO0), "--b-avg-sq", str(B2), "--r0", "2.5",
                        "--rmax", "25", "--nr", "1000",
                        "--out-dir", str(Path(tmp) / f"out_{name}"),
                    ])
                runs[name] = (code, out.getvalue(), err.getvalue())
        self.assertEqual(runs["asc"][0], 0)
        self.assertEqual(runs["desc"][0], 0)
        line = [row for row in runs["asc"][1].splitlines() if "result" in row][0]
        self.assertIn(line, runs["desc"][1])
        code, _, err = runs["inverted"]
        self.assertEqual(code, 2)
        self.assertIn("non-physical scale", err)


class TransformGridTests(unittest.TestCase):
    def setUp(self):
        self.q, sq = model_sq()
        self.fq = sq_to_fq(self.q, sq)
        self.sq = sq
        self.r = np.arange(1, 501) * 0.02

    def test_descending_q_raises(self):
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            fq_to_gpdf(self.q[::-1], self.fq[::-1], self.r)
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            low_q_correction_basis(self.q[::-1], self.r)
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            fourier_filter(self.q[::-1], self.sq[::-1], self.r, rho0=RHO0, cutoff=1.0)

    def test_descending_r_raises(self):
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            gpdf_to_fq(self.r[::-1], np.ones_like(self.r), self.q)
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            fourier_filter(self.q, self.sq, self.r[::-1], rho0=RHO0, cutoff=1.0)

    def test_repeated_or_nan_grid_raises(self):
        q = self.q.copy()
        q[10] = q[9]
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            sine_transform(q, self.fq, self.r)
        q[10] = np.nan
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            sine_transform(q, self.fq, self.r)

    def test_output_grid_may_be_scalar(self):
        full = fq_to_gpdf(self.q, self.fq, self.r)
        single = fq_to_gpdf(self.q, self.fq, self.r[123])
        self.assertEqual(np.shape(single), ())
        self.assertAlmostEqual(float(single), full[123], places=12)

    def test_short_sections_still_integrate_to_zero(self):
        self.assertEqual(sine_transform(np.array([]), np.array([]), self.q).max(), 0.0)
        self.assertEqual(sine_transform(np.array([0.5]), np.array([1.0]), self.q).max(), 0.0)


if __name__ == "__main__":
    unittest.main()
