# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Auto StoG low-r window placement: never fit the density limit across the first shell.

Regression tests for the 1.0 audit (stog-a group). Without r0 / r_fit_max the first
pass used to fit g = 0 on the blind window [r_cutoff + 0.2, r_cutoff + 1.2] =
[1.2, 2.2] A, which contains the first shell of most oxides (Ti-O 1.95, Re-O 1.875,
Si-O 1.61 A); the scale came out 40-60 % low or negative, and short bonds were
refused refinement and returned as if fine.

The crystal models are physically consistent neutron totals: a periodic supercell with
Gaussian thermal displacements, partial g_ij from pair counts, Faber-Ziman weights from
the repo's Sears table, and S(Q) by the repo's own transform; measured = (S + 9) / 10,
so the true scale is a = 10.
"""

import unittest

import numpy as np

# numpy >= 2.0 renamed trapz; the package supports numpy >= 1.22.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz
from scipy.spatial import cKDTree

from rmc_toolkits.scaling import ScalingConfig, autoscale, diagnostics_summary
from rmc_toolkits.scattering import faber_ziman
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

A_TRUE, B_TRUE = 10.0, -9.0
Q = np.arange(50, 3001) * 0.01  # 0.5 .. 30 A^-1

PEROVSKITE = {"Sr": [(0, 0, 0)], "Ti": [(0.5, 0.5, 0.5)],
              "O": [(0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]}
REO3 = {"Re": [(0, 0, 0)], "O": [(0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5)]}


def crystal_sq(formula, sites, a0, n=6, u_iso=0.006, seed=7, dr=0.005, q=Q):
    """Measured S(Q) of an n^3 cubic supercell with Gaussian thermal displacements."""
    rng = np.random.default_rng(seed)
    grid = np.stack(np.meshgrid(*[np.arange(n)] * 3, indexing="ij"), -1).reshape(-1, 3)
    box = n * a0
    positions = {}
    for element, fractional in sites.items():
        cell = np.concatenate([(grid + np.array(site)) * a0 for site in fractional])
        positions[element] = np.mod(cell + rng.normal(0, np.sqrt(u_iso), cell.shape), box)
    fz = faber_ziman(formula)
    rho0 = sum(len(p) for p in positions.values()) / box**3
    rmax = box / 2 - 0.1
    edges = np.arange(0, rmax + dr, dr)
    r = 0.5 * (edges[1:] + edges[:-1])
    trees = {element: cKDTree(p, boxsize=box) for element, p in positions.items()}
    g = np.zeros_like(r)
    for i in positions:
        for j in positions:
            cumulative = trees[i].count_neighbors(trees[j], edges).astype(float)
            if i == j:
                cumulative -= len(positions[i])
            g_ij = np.diff(cumulative) / (
                len(positions[i]) * fz.fractions[j] * rho0 * 4 * np.pi * r**2 * dr
            )
            weight = fz.fractions[i] * fz.fractions[j] * fz.b_coh_fm[i] * fz.b_coh_fm[j]
            g += weight / fz.b_avg_sq_fm2 * g_ij
    taper = np.where(r < rmax / 2, 1.0, 0.5 * (1 + np.cos(np.pi * (r - rmax / 2) / (rmax / 2))))
    g = 1.0 + (g - 1.0) * taper
    sq_true = fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, rho0), q))
    config = dict(
        qmin=float(q[0]), qmax=30.0, rho0=rho0,
        b_avg_sq=fz.b_avg_sq_barn, b_sq_avg=fz.b_sq_avg_barn,
    )
    return (sq_true - B_TRUE) / A_TRUE, config


def shell_sq(formula, rho0, shells, r_continuum, q=Q):
    """Analytic Gaussian-shell FZ total (+ smooth continuum) and its consistent S(0)."""
    fz = faber_ziman(formula)
    r = np.arange(1, 12001) * 0.005
    g = 0.5 * (1.0 + np.tanh((r - r_continuum) / 0.25))
    for i, j, distance, cn, sigma in shells:
        amplitude = (1.0 if i == j else 2.0) * fz.fractions[i] * fz.b_coh_fm[i] * fz.b_coh_fm[j] * cn
        amplitude /= fz.b_avg_sq_fm2
        g += amplitude * np.exp(-0.5 * ((r - distance) / sigma) ** 2) / (
            4 * np.pi * r**2 * rho0 * sigma * np.sqrt(2 * np.pi)
        )
    s0 = 1.0 + 4 * np.pi * rho0 * _trapezoid(r**2 * (g - 1.0), r)
    sq_true = fq_to_sq(q, gpdf_to_fq(r, g_to_gpdf(r, g, rho0), q))
    config = dict(
        qmin=float(q[0]), qmax=30.0, rho0=rho0, b_avg_sq=fz.b_avg_sq_barn,
        b_sq_avg=fz.b_sq_avg_barn, s0_target=s0,
    )
    return (sq_true - B_TRUE) / A_TRUE, config


SIO2_GLASS = ("SiO2", 0.0663, [("Si", "O", 1.61, 4, 0.05), ("O", "O", 2.63, 6, 0.09),
                               ("Si", "Si", 3.08, 4, 0.1)], 3.6)
B2O3_LIKE = ("B2O3", 0.08, [("B", "O", 1.37, 3, 0.045), ("O", "O", 2.38, 4, 0.08),
                            ("B", "B", 2.45, 3, 0.09)], 3.4)


class OxideFirstShellWindowTests(unittest.TestCase):
    """Composition + Q window only: the window must end below the M-O first shell."""

    def check(self, formula, sites, a0, first_shell):
        sq, values = crystal_sq(formula, sites, a0)
        config = ScalingConfig(**values)
        result = autoscale(Q, sq, config)
        summary = diagnostics_summary(result, config)
        self.assertLess(abs(result.a / A_TRUE - 1.0), 0.02, f"a = {result.a}")
        self.assertLess(summary["r0_detected"], first_shell)
        self.assertGreater(summary["r0_detected"], first_shell - 0.3)
        self.assertLess(summary["r_fit_window"][1], first_shell - 0.25)
        self.assertTrue(summary["window_refined"])
        return result

    def test_srtio3_inverted_ti_o_shell(self):
        # Pre-1.0: window [1.2, 2.34] across Ti-O, a = 5.2 (-48 %).
        self.check("SrTiO3", PEROVSKITE, 3.905, 1.95)

    def test_reo3_short_m_o_shell(self):
        # Pre-1.0: window [1.2, 2.2] across Re-O, a < 0 (the data inverted).
        self.check("ReO3", REO3, 3.75, 1.875)

    def test_fz_mode_reports_the_first_shell(self):
        sq, values = crystal_sq("SrTiO3", PEROVSKITE, 3.905)
        config = ScalingConfig(**values, amplitude_criterion="fz")
        result = autoscale(Q, sq, config)
        self.assertLess(result.provenance["r0_detected"], 1.95)
        self.assertLess(result.provenance["r_fit_window"][1], 1.7)


class ShortBondTests(unittest.TestCase):
    """First shells too close to r_cutoff: fail loudly, never fit across them."""

    def test_short_bond_below_1_45_A_fails_loudly(self):
        formula, rho0, shells, r_continuum = B2O3_LIKE
        sq, values = shell_sq(formula, rho0, shells, r_continuum)
        # Pre-1.0: returned a = -0.65 fitted on [1.2, 2.2] with no error.
        with self.assertRaisesRegex(ValueError, "r_cutoff"):
            autoscale(Q, sq, ScalingConfig(**values))

    def test_si_o_shell_needs_a_lower_cutoff(self):
        formula, rho0, shells, r_continuum = SIO2_GLASS
        sq, values = shell_sq(formula, rho0, shells, r_continuum)
        with self.assertRaisesRegex(ValueError, r"starts at 1\.\d\d A.*r_cutoff to <= 0\.9"):
            autoscale(Q, sq, ScalingConfig(**values))
        config = ScalingConfig(**values, r_cutoff=0.7)
        result = autoscale(Q, sq, config)
        summary = diagnostics_summary(result, config)
        # A sharp (sigma 0.05 A) Si-O shell leaves strong termination ripples in
        # the narrow [0.9, 1.27] window; the point here is that the window sits
        # below the shell and the scale is sane (pre-1.0: a = 4.7 of 5 with the
        # window silently on the shell's foot, or a < 0 on [1.2, 2.2]).
        self.assertLess(abs(result.a / A_TRUE - 1.0), 0.08, f"a = {result.a}")
        self.assertLess(summary["r_fit_window"][1], 1.35)

    def test_unfittable_data_reports_the_real_error(self):
        # Both trial fits fail for a reason unrelated to the first shell: that
        # reason is reported, not "could not locate the first shell".
        formula, rho0, shells, r_continuum = SIO2_GLASS
        sq, values = shell_sq(formula, rho0, shells, r_continuum)
        with self.assertRaisesRegex(ValueError, "fewer than 16 usable"):
            autoscale(Q[:10], sq[:10], ScalingConfig(**values))

    def test_pinned_window_is_respected(self):
        formula, rho0, shells, r_continuum = SIO2_GLASS
        sq, values = shell_sq(formula, rho0, shells, r_continuum)
        config = ScalingConfig(**values, r_cutoff=0.7, r0=1.5)
        result = autoscale(Q, sq, config)
        np.testing.assert_allclose(result.provenance["r_fit_window"], [0.9, 1.25])
        self.assertNotIn("window_refined", result.provenance)

    def test_pinned_sliver_window_is_refused(self):
        # A pinned density-limit window narrower than MIN_AUTO_WINDOW (0.1 A)
        # is set by one truncation ripple: on FeCoSn 199 K a 0.01-0.05 A window
        # gave scales 11-43 % low, reported converged with no flag.
        formula, rho0, shells, r_continuum = SIO2_GLASS
        sq, values = shell_sq(formula, rho0, shells, r_continuum)
        for pins in ({"r0": 1.2}, {"r_fit_max": 0.95}, {"r_fit_min": 1.2, "r_fit_max": 1.25}):
            with self.subTest(**pins):
                config = ScalingConfig(**values, r_cutoff=0.7, **pins)
                with self.assertRaisesRegex(ValueError, "narrower than 0.1"):
                    autoscale(Q, sq, config)
        # The Q->0 amplitude does not depend on the window: fz still runs.
        fz = ScalingConfig(**values, r_cutoff=0.7, r0=1.2, amplitude_criterion="fz")
        self.assertTrue(np.isfinite(autoscale(Q, sq, fz).a))


class RCutoffValidationTests(unittest.TestCase):
    def test_r_cutoff_must_be_finite_and_non_negative(self):
        base = dict(qmin=0.5, qmax=30.0, rho0=0.05, b_avg_sq=0.02)
        for bad in (-1.0, -1e-9, float("nan"), float("inf")):
            with self.subTest(r_cutoff=bad):
                with self.assertRaisesRegex(ValueError, "r_cutoff"):
                    ScalingConfig(**base, r_cutoff=bad)
        ScalingConfig(**base, r_cutoff=0.0)  # no Fourier filter: allowed
        for name in ("r0", "r_fit_min", "r_fit_max"):
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, name):
                    ScalingConfig(**base, **{name: float("nan")})
        with self.assertRaisesRegex(ValueError, "r_fit_min"):
            ScalingConfig(**base, r_fit_min=-0.5)


if __name__ == "__main__":
    unittest.main()
