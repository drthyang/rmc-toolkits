# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Automatic scale/offset determination for total-scattering S(Q).

Replaces the classic stog "try again" manual loop with a minimizer inside a
self-consistent Fourier-filter loop. The correction model is the affine map

    ``S_corr(Q) = a * S_meas(Q) + b``            (multiply convention; the
                                                  classic stog input's
                                                  yoffset/yscale map to
                                                  ``a = 1/yscale, b = yoffset``)

fitted against two physics targets (Keen 2001):

- C1 — high-Q asymptote: ``S_corr -> 1`` over a tail window (Eq. 21).
- C2 — low-r density limit: the transformed ``g_corr(r) -> 0`` below the
  closest interatomic approach (Eqs. 15/29 in g-space).

Both residual blocks are affine in ``(a, b)`` because the sine transform is
linear in the data, so each inner fit is a closed-form linear least-squares
solve. The outer loop alternates fitting with the Fourier filter (whose
subtraction term is held fixed during each fit) until ``(a, b)`` converge.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .transforms import (
    enforce_low_r,
    fourier_filter,
    fq_to_gpdf,
    gk_to_dr,
    g_to_gk,
    low_q_correction_basis,
    sq_to_fk,
    sq_to_fq,
)


@dataclass(frozen=True)
class ScalingConfig:
    """Parameters for the auto-scaling pipeline.

    ``b_avg_sq`` is ``<b>^2 = (sum_i c_i b_i)^2`` in barns — the classic stog
    input's "Faber-Ziman coefficient" line. ``b_sq_avg`` is the *different*
    number ``<b^2> = sum_i c_i b_i^2`` (Keen Eq. 14 Q->0 limit), used only for
    diagnostics when provided.

    ``low_q_correction`` defaults on: measured data always omit ``[0, Qmin]``,
    and without the analytic correction that omission biases the fitted scale
    (8% on the synthetic benchmark vs 0.3% with it). Disable only for strict
    classic-stog parity. ``lorch=True`` additionally damps termination ripples
    and speeds loop convergence at the cost of real-space resolution.

    Scope note: the fit recovers the *true* scale only for data that can
    satisfy the physics targets (well-reduced, low enough Qmin). Data missing
    structure below Qmin violate the absolute sum rules; the fit still
    converges and minimizes unphysical low-r content, but a smooth low-Q
    deficiency is generically ABSORBED into a biased scale with all residual
    diagnostics clean — there is no internal signature that certifies the
    absolute scale. :func:`diagnostics_summary`'s ``density_limit_satisfied``
    is one-sided (False proves non-recoverability; True does not prove
    correctness); validate the absolute scale externally whenever Qmin is not
    small. ``c2_bins > 0`` switches C2 from pointwise residuals to binned mean
    levels; pointwise (default) is the stable choice on ripple-heavy data.
    """

    qmin: float
    qmax: float
    rho0: float
    b_avg_sq: float
    b_sq_avg: float | None = None
    r_cutoff: float = 1.0
    r0: float | None = None
    fit_offset: bool = True
    q_tail_frac: float = 0.15
    r_fit_min: float | None = None
    r_fit_max: float | None = None
    rmax: float = 50.0
    nr: int = 5000
    lorch: bool = False
    low_q_correction: bool = True
    c2_weight: float = 1.0
    c2_bins: int = 0
    max_iter: int = 50
    tol: float = 1.0e-6
    enforce_cutoff: float | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.rho0) or self.rho0 <= 0:
            raise ValueError(f"rho0 must be finite and positive, got {self.rho0}")
        if not np.isfinite(self.b_avg_sq) or self.b_avg_sq <= 0:
            raise ValueError(f"b_avg_sq must be finite and positive, got {self.b_avg_sq}")
        if self.qmax <= self.qmin:
            raise ValueError("qmax must exceed qmin")
        if int(self.nr) != self.nr or self.nr <= 0:
            raise ValueError(f"nr must be a positive integer, got {self.nr}")
        if not np.isfinite(self.rmax) or self.rmax <= 0:
            raise ValueError(f"rmax must be finite and positive, got {self.rmax}")

    @property
    def r_grid(self) -> np.ndarray:
        step = self.rmax / self.nr
        return np.arange(1, self.nr + 1) * step

    @property
    def q_tail_min(self) -> float:
        return self.qmax - self.q_tail_frac * (self.qmax - self.qmin)

    @property
    def r_fit_window(self) -> tuple[float, float]:
        lo = self.r_fit_min if self.r_fit_min is not None else self.r_cutoff + 0.2
        if self.r_fit_max is not None:
            hi = self.r_fit_max
        elif self.r0 is not None:
            hi = self.r0 - 0.25
        else:
            hi = lo + 1.0
        if hi <= lo:
            raise ValueError(f"empty low-r fit window [{lo}, {hi}]")
        return lo, hi


@dataclass(frozen=True)
class ScalingResult:
    """Auto-scaling outcome: fitted correction, arrays, and diagnostics.

    All arrays live on ``q`` / ``r``. ``gk``/``d_r``/``fk`` are the
    RMCProfile-ready functions; ``gk_enforced``/``d_r_enforced`` additionally
    apply the classic stog low-r replacement when ``enforce_cutoff`` is set.
    ``low_r_rms`` is the honest, *pre*-enforcement residual of ``g_filtered``
    against 0 over the fit window.
    """

    a: float
    b: float
    converged: bool
    iterations: int
    history: tuple[tuple[float, float, float], ...]
    q: np.ndarray
    sq_scaled: np.ndarray
    sq_filtered: np.ndarray
    sq_ft: np.ndarray
    r: np.ndarray
    g_filtered: np.ndarray
    gk: np.ndarray
    d_r: np.ndarray
    fk: np.ndarray
    gk_enforced: np.ndarray | None
    d_r_enforced: np.ndarray | None
    low_r_rms: float
    c1_tail_mean: float
    provenance: dict[str, Any] = field(repr=False)


def _bin_means(
    col_a: np.ndarray, col_b: np.ndarray, rhs: np.ndarray, n_bins: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Average LSQ rows over contiguous bins, weighted by sqrt(bin size)."""
    edges = np.linspace(0, col_a.size, min(n_bins, col_a.size) + 1).astype(int)
    out_a, out_b, out_rhs = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if hi <= lo:
            continue
        weight = np.sqrt(hi - lo)
        out_a.append(weight * col_a[lo:hi].mean())
        out_b.append(weight * col_b[lo:hi].mean())
        out_rhs.append(weight * rhs[lo:hi].mean())
    return np.asarray(out_a), np.asarray(out_b), np.asarray(out_rhs)


def _fit_windows(q: np.ndarray, r: np.ndarray, config: ScalingConfig):
    tail = q >= config.q_tail_min
    lo, hi = config.r_fit_window
    window = (r >= lo) & (r <= hi)
    if tail.sum() < 2 or window.sum() < 2:
        raise ValueError("fit windows contain fewer than 2 points")
    return tail, window


def _solve_affine(
    q: np.ndarray,
    sq: np.ndarray,
    delta_sq: np.ndarray,
    r: np.ndarray,
    tail: np.ndarray,
    window: np.ndarray,
    config: ScalingConfig,
) -> tuple[float, float]:
    """Closed-form LSQ for ``S_eff = a*sq + b - delta_sq`` against C1 and C2.

    C1 rows: ``a*sq + b - delta_sq - 1 = 0`` on the tail window.
    C2 rows: ``g_eff(r) = 0`` on the low-r window, where the transform of the
    affine ``S_eff`` decomposes into precomputed basis transforms.
    """
    rw = r[window]
    # Basis transforms on the fit window only (cheap: len(rw) outputs). The
    # omitted-low-Q correction is folded in analytically below — computing it
    # inside each basis transform would multiply-count its constant term.
    g_data = fq_to_gpdf(q, sq_to_fq(q, sq), rw, lorch=config.lorch)
    g_one = fq_to_gpdf(q, q, rw, lorch=config.lorch)  # transform of Q*1 (offset)
    g_delta = fq_to_gpdf(q, q * delta_sq, rw, lorch=config.lorch)  # filter (fixed)
    # G_PDF[a*sq + b - delta - 1] = a*(g_data + g_one) + b*g_one - g_delta - g_one
    # Low-Q correction (affine in S_eff(qmin) = a*sq[0] + b - delta[0]):
    #   + coef_s0*(a*sq[0] + b - delta[0]) - const
    if config.low_q_correction:
        coef_s0, const = low_q_correction_basis(q, rw, lorch=config.lorch)
    else:
        coef_s0 = const = np.zeros_like(rw)
    denom = 4.0 * np.pi * config.rho0 * rw
    w2 = np.sqrt(config.c2_weight)

    col_a_c1 = sq[tail]
    col_b_c1 = np.ones(int(tail.sum()))
    rhs_c1 = 1.0 + delta_sq[tail]

    col_a_c2 = w2 * (g_data + g_one + coef_s0 * sq[0]) / denom
    col_b_c2 = w2 * (g_one + coef_s0) / denom
    rhs_c2 = w2 * (
        (g_delta + g_one + coef_s0 * delta_sq[0] + const) / denom - 1.0
    )

    if config.c2_bins > 0:
        # Fit binned mean *levels* rather than every ripple point (the classic
        # "level matching" criterion). Pointwise residuals let the solver shrink
        # termination ripples by shrinking the scale — averaging removes that
        # bias while the bins keep enough shape to pin both a and b.
        col_a_c2, col_b_c2, rhs_c2 = _bin_means(
            col_a_c2, col_b_c2, rhs_c2, config.c2_bins
        )

    col_a = np.concatenate([col_a_c1, col_a_c2])
    col_b = np.concatenate([col_b_c1, col_b_c2])
    rhs = np.concatenate([rhs_c1, rhs_c2])

    if config.fit_offset:
        design = np.column_stack([col_a, col_b])
        solution, *_ = np.linalg.lstsq(design, rhs, rcond=None)
        return float(solution[0]), float(solution[1])
    solution = float(np.dot(col_a, rhs) / np.dot(col_a, col_a))
    return solution, 0.0


def _pipeline(
    q: np.ndarray,
    sq_scaled: np.ndarray,
    config: ScalingConfig,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fourier filter of the scaled data on the configured r grid."""
    return fourier_filter(
        q,
        sq_scaled,
        config.r_grid,
        rho0=config.rho0,
        cutoff=config.r_cutoff,
        lorch=config.lorch,
        low_q_correction=config.low_q_correction,
    )


def _low_r_rms(r: np.ndarray, g_filtered: np.ndarray, config: ScalingConfig) -> float:
    lo, hi = config.r_fit_window
    window = (r >= lo) & (r <= hi)
    return float(np.sqrt(np.mean(g_filtered[window] ** 2)))


def crop_sq(q: np.ndarray, sq: np.ndarray, config: ScalingConfig):
    """Crop to ``(0, qmax]`` ∩ ``[qmin, qmax]``, dropping non-finite rows.

    ``Q <= 0`` rows are dropped unconditionally: the S(Q) conversions divide
    by Q, and the analytic low-Q correction already models the omitted
    ``[0, Qmin]`` range, so a Q = 0 point carries no usable information.
    """
    q = np.asarray(q, dtype=float)
    sq = np.asarray(sq, dtype=float)
    keep = (
        np.isfinite(q)
        & np.isfinite(sq)
        & (q > 0)
        & (q >= config.qmin - 1e-12)
        & (q <= config.qmax + 1e-12)
    )
    if keep.sum() < 16:
        raise ValueError("fewer than 16 usable S(Q) points after cropping")
    return q[keep], sq[keep]


def scale_pipeline(
    q: np.ndarray,
    sq: np.ndarray,
    config: ScalingConfig,
    a: float,
    b: float,
    *,
    converged: bool = True,
    iterations: int = 0,
    history: tuple[tuple[float, float, float], ...] = (),
) -> ScalingResult:
    """Run the fixed-(a, b) pipeline: scale, filter, convert, (optionally) enforce.

    This is the manual-mode entry point (classic stog parity) and the final
    stage of :func:`autoscale`.
    """
    q, sq = crop_sq(q, sq, config)
    r = config.r_grid
    sq_scaled = a * sq + b
    sq_filtered, sq_ft, g_filtered = _pipeline(q, sq_scaled, config)

    gk = g_to_gk(g_filtered, config.b_avg_sq)
    d_r = gk_to_dr(r, gk, config.rho0)
    fk = sq_to_fk(sq_filtered, config.b_avg_sq)

    gk_enforced = d_r_enforced = None
    if config.enforce_cutoff is not None:
        gk_enforced = enforce_low_r(
            r, gk, cutoff=config.enforce_cutoff, b_avg_sq=config.b_avg_sq
        )
        d_r_enforced = gk_to_dr(r, gk_enforced, config.rho0)

    tail, _ = _fit_windows(q, r, config)
    provenance = {
        "model": "S_corr = a*S_meas + b",
        "mode": "manual",
        "a": a,
        "b": b,
        "config": {
            key: getattr(config, key)
            for key in (
                "qmin", "qmax", "rho0", "b_avg_sq", "b_sq_avg", "r_cutoff", "r0",
                "fit_offset", "q_tail_frac", "rmax", "nr", "lorch",
                "low_q_correction", "c2_weight", "c2_bins", "max_iter", "tol",
                "enforce_cutoff",
            )
        },
        "q_tail_window": [config.q_tail_min, config.qmax],
        "r_fit_window": list(config.r_fit_window),
        "n_q_points": int(q.size),
    }
    return ScalingResult(
        a=a,
        b=b,
        converged=converged,
        iterations=iterations,
        history=history,
        q=q,
        sq_scaled=sq_scaled,
        sq_filtered=sq_filtered,
        sq_ft=sq_ft,
        r=r,
        g_filtered=g_filtered,
        gk=gk,
        d_r=d_r,
        fk=fk,
        gk_enforced=gk_enforced,
        d_r_enforced=d_r_enforced,
        low_r_rms=_low_r_rms(r, g_filtered, config),
        c1_tail_mean=float(sq_filtered[tail].mean()),
        provenance=provenance,
    )


def autoscale(q: np.ndarray, sq: np.ndarray, config: ScalingConfig) -> ScalingResult:
    """Automatically determine ``(a, b)`` and run the full pipeline.

    Self-consistent loop: fit ``(a, b)`` with the Fourier-filter subtraction
    term held fixed, re-run the filter on the corrected data, repeat until the
    parameters converge (``tol``) or ``max_iter`` is reached.
    """
    q, sq = crop_sq(q, sq, config)
    r = config.r_grid
    tail, window = _fit_windows(q, r, config)

    delta_sq = np.zeros_like(q)
    a_prev = b_prev = np.inf
    a = 1.0
    b = 0.0
    history: list[tuple[float, float, float]] = []
    converged = False
    iterations = 0

    for iterations in range(1, config.max_iter + 1):
        a, b = _solve_affine(q, sq, delta_sq, r, tail, window, config)
        sq_scaled = a * sq + b
        _, sq_ft, g_filtered = _pipeline(q, sq_scaled, config)
        delta_sq = sq_ft - 1.0
        history.append((a, b, _low_r_rms(r, g_filtered, config)))
        if abs(a - a_prev) <= config.tol * max(1.0, abs(a)) and abs(
            b - b_prev
        ) <= config.tol * max(1.0, abs(b)):
            converged = True
            break
        a_prev, b_prev = a, b

    result = scale_pipeline(
        q,
        sq,
        config,
        a,
        b,
        converged=converged,
        iterations=iterations,
        history=tuple(history),
    )
    result.provenance["mode"] = "auto"
    return result


def diagnostics_summary(result: ScalingResult, config: ScalingConfig) -> dict[str, Any]:
    """Compact, human-readable fit-quality summary (JSON-friendly)."""
    lo, hi = config.r_fit_window
    window = (result.r >= lo) & (result.r <= hi)
    g_window_mean = float(result.g_filtered[window].mean())
    summary: dict[str, Any] = {
        "a": result.a,
        "b": result.b,
        "converged": result.converged,
        "iterations": result.iterations,
        "c1_tail_mean": result.c1_tail_mean,
        "low_r_rms_pre_enforcement": result.low_r_rms,
        "g_window_mean": g_window_mean,
        "gk_low_r_theory": -config.b_avg_sq,
        "d_r_low_r_slope_theory": -4.0 * np.pi * config.rho0 * config.b_avg_sq,
        # ONE-SIDED diagnostic. False proves the density limit could not be
        # satisfied by any (a, b) — the absolute scale is definitely not
        # recoverable from this data (e.g. severe missing-low-Q structure).
        # True only means the fit reached its C2 target; smooth low-Q
        # deficiencies are generically ABSORBED into a biased scale with all
        # residuals clean, so True does NOT certify the absolute scale —
        # validate externally (known composition/density, or fixed-a via
        # scale_pipeline) whenever Qmin is not small.
        "density_limit_satisfied": bool(abs(g_window_mean) < 0.1),
    }
    if config.b_sq_avg is not None:
        # Keen Eq. 14 diagnostic: FK(Q->0) -> -<b^2>. Report the lowest-Q value
        # actually measured for comparison (data rarely reach Q ~ 0).
        summary["fk_qmin"] = float(result.fk[0])
        summary["fk_q0_theory"] = -config.b_sq_avg
    return summary
