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

from dataclasses import dataclass, field, replace
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
    (8% on the synthetic benchmark vs 0.3% with it). The correction
    extrapolates S(Q) linearly to ``s0_target`` at Q = 0; ``s0_target=None``
    (default) resolves to the composition-derived Keen Eq. 21 limit
    ``1 - <b^2>/<b>^2`` when ``b_sq_avg`` is supplied, else to the pystog
    solid-state 0 — for negative-b compositions (Mn3Sn: S(0) ~ -12) the
    composition-aware target removes an O(1) low-r bias. Disable the
    correction only for strict classic-stog parity. ``lorch=True``
    additionally damps termination ripples and speeds loop convergence at the
    cost of real-space resolution.

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

    ``c1_mode`` selects the fit architecture. ``"sweep"`` (default) first
    determines the high-Q level by the criterion-driven two-sided window
    search (:func:`level_sweep`) and ties the offset to it (``b = 1 - a *
    level``), leaving the density limit a single amplitude dof — the
    "shift by the level, then scale" decomposition; ``fit_offset`` is
    superseded in this mode. ``"joint"`` is the original 2-dof fit (and the
    fallback when no statistically flat window exists).

    ``amplitude_criterion`` selects what pins the overall amplitude once the
    level anchors the offset. ``"density"`` (default) fits it against the
    low-r density limit (C2). ``"fz"`` instead sets it from the Q->0
    Faber-Ziman limit (Keen Eq. 21) via :func:`amplitude_from_fz_limit` — the
    "subtract the high-Q level, scale Q=0 onto S(0) = 1 - <b^2>/<b>^2, shift
    the level back to 1" procedure. It needs ``b_sq_avg``, requires
    ``c1_mode="sweep"`` (no level, no anchor), involves no self-consistent
    loop (the criterion is filter-independent), and inherits the caveat that
    a long Qmin -> 0 extrapolation owns the whole scale; the density-limit
    residuals in :func:`diagnostics_summary` then act as the independent
    cross-check.

    Robustness controls: ``robust`` (default on) runs a Huber IRLS loop over
    the joint fit so isolated outliers cannot drag the closed-form solution;
    per-point ``sigma`` (pass to :func:`autoscale`) 1/sigma-weights the high-Q
    rows. ``despike`` (OFF by default) drops narrow rolling-median outliers
    from the input before any transform — measured to restore clean recovery
    under detector-glitch spikes (which otherwise ring through the transform
    into the C2 window, a channel IRLS cannot reject) — but it also flags real
    Bragg maxima on crystalline data (12% of points on the 59438 benchmark),
    so enable it only for glitch-type contamination and check the reported
    ``n_despiked``. ``c1_slope_nuisance`` (experimental) adds a linear tail
    term absorbing f(Q)/Placzek-like drift in the C1 level estimate; a drift
    spanning the whole Q range also enters through the transform, which this
    does not correct.
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
    robust: bool = True
    c1_mode: str = "sweep"
    amplitude_criterion: str = "density"
    s0_target: float | None = None
    c1_slope_nuisance: bool = False
    despike: bool = False
    despike_window: int = 7
    despike_nsigma: float = 6.0
    max_iter: int = 50
    tol: float = 1.0e-6
    enforce_cutoff: float | None = None

    def __post_init__(self) -> None:
        if self.c1_mode not in ("sweep", "joint"):
            raise ValueError(f"c1_mode must be 'sweep' or 'joint', got {self.c1_mode!r}")
        if self.amplitude_criterion not in ("density", "fz"):
            raise ValueError(
                "amplitude_criterion must be 'density' or 'fz', "
                f"got {self.amplitude_criterion!r}"
            )
        if self.amplitude_criterion == "fz":
            if self.b_sq_avg is None:
                raise ValueError("amplitude_criterion='fz' requires b_sq_avg (<b^2>)")
            if self.c1_mode != "sweep":
                raise ValueError(
                    "amplitude_criterion='fz' requires c1_mode='sweep' "
                    "(the measured level anchors the offset)"
                )
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
    def effective_s0_target(self) -> float:
        """S(0) used by the omitted-low-Q extrapolation (see class docstring)."""
        if self.s0_target is not None:
            return float(self.s0_target)
        if self.b_sq_avg is not None:
            return 1.0 - float(self.b_sq_avg) / float(self.b_avg_sq)
        return 0.0

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
    sweep: LevelSweepResult | None
    a_fz: float | None
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


@dataclass(frozen=True)
class LevelSweepResult:
    """Outcome of the two-sided high-Q level sweep.

    ``level`` is the asymptote of the *measured* S(Q) over the optimal window
    ``(q_lo, q_hi)``; ``level_uncertainty`` is the spread of the level across
    every admissible window (the honest error bar the community question asks
    for). ``asymptote_found`` is False when no window has a statistically-zero
    slope — the data never reach a trustworthy asymptotic regime.
    """

    level: float
    level_uncertainty: float
    q_lo: float
    q_hi: float
    slope: float
    slope_sigma: float
    asymptote_found: bool
    n_admissible: int


def level_sweep(
    q: np.ndarray,
    sq: np.ndarray,
    *,
    min_width: float = 3.0,
    n_grid: int = 80,
    slope_nsigma: float = 2.0,
) -> LevelSweepResult:
    """Determine the high-Q level of S(Q) by a criterion-driven window search.

    Answers "what is high enough Q?" without hand-set tolerances: every
    candidate window ``[Q1, Q2]`` (both edges swept on an ``n_grid`` grid,
    width >= ``min_width``) gets an OLS line fit in O(1) via prefix sums. A
    window is **admissible** iff its slope is statistically zero given the
    fit's own residual noise (``|slope| < slope_nsigma * sigma_slope``) — the
    data define "flat", not a tolerance. Among admissible windows the one
    with the smallest predicted level variance wins (longest, quietest
    stretch); the level spread across all admissible windows is returned as
    the uncertainty. End-of-range artifacts (detector rolloff, dead tails)
    exclude themselves because any window touching them fails admissibility.
    """
    q = np.asarray(q, dtype=float)
    sq = np.asarray(sq, dtype=float)
    finite = np.isfinite(q) & np.isfinite(sq)
    q, sq = q[finite], sq[finite]
    if q.size < 32:
        raise ValueError("level_sweep needs at least 32 finite points")

    # Prefix sums for O(1) per-window OLS of y = c0 + c1*q.
    ones = np.concatenate([[0.0], np.cumsum(np.ones_like(q))])
    sum_q = np.concatenate([[0.0], np.cumsum(q)])
    sum_q2 = np.concatenate([[0.0], np.cumsum(q * q)])
    sum_y = np.concatenate([[0.0], np.cumsum(sq)])
    sum_qy = np.concatenate([[0.0], np.cumsum(q * sq)])
    sum_y2 = np.concatenate([[0.0], np.cumsum(sq * sq)])

    edges = np.unique(np.linspace(0, q.size - 1, n_grid).astype(int))
    best = None
    admissible_levels: list[float] = []

    for ii, i in enumerate(edges[:-1]):
        for j in edges[ii + 1 :]:
            if q[j] - q[i] < min_width:
                continue
            n = ones[j + 1] - ones[i]
            sq_sum = sum_q[j + 1] - sum_q[i]
            sq2 = sum_q2[j + 1] - sum_q2[i]
            sy = sum_y[j + 1] - sum_y[i]
            sqy = sum_qy[j + 1] - sum_qy[i]
            sy2 = sum_y2[j + 1] - sum_y2[i]
            det = n * sq2 - sq_sum * sq_sum
            if det <= 0 or n < 24:
                continue
            slope = (n * sqy - sq_sum * sy) / det
            intercept = (sq2 * sy - sq_sum * sqy) / det
            # Residual variance and standard errors.
            rss = sy2 - intercept * sy - slope * sqy
            dof = n - 2
            sigma2 = max(rss / dof, 0.0)
            var_slope = sigma2 * n / det
            qbar = sq_sum / n
            # Level = fitted value at the window centre; its variance.
            level = intercept + slope * qbar
            var_level = sigma2 / n
            slope_sigma = np.sqrt(max(var_slope, 1e-30))
            if abs(slope) < slope_nsigma * slope_sigma:
                admissible_levels.append(level)
                score = var_level
                if best is None or score < best[0]:
                    best = (score, level, q[i], q[j], slope, slope_sigma)

    if best is None:
        # No statistically flat window: report the least-sloped one honestly.
        return LevelSweepResult(
            level=float(np.median(sq[q >= q.max() - min_width])),
            level_uncertainty=float("nan"),
            q_lo=float(q.max() - min_width),
            q_hi=float(q.max()),
            slope=float("nan"),
            slope_sigma=float("nan"),
            asymptote_found=False,
            n_admissible=0,
        )
    _, level, q_lo, q_hi, slope, slope_sigma = best
    return LevelSweepResult(
        level=float(level),
        level_uncertainty=float(np.std(admissible_levels)),
        q_lo=float(q_lo),
        q_hi=float(q_hi),
        slope=float(slope),
        slope_sigma=float(slope_sigma),
        asymptote_found=True,
        n_admissible=len(admissible_levels),
    )


def detect_first_peak_onset(
    r: np.ndarray,
    g: np.ndarray,
    qmax: float,  # noqa: ARG001 - kept for signature stability / future ripple use
    *,
    search_min: float = 1.0,
    search_max: float = 6.0,
    fraction: float = 0.35,
    floor: float = 0.5,
) -> float | None:
    """Data-derived closest-approach r0: the left flank of the first shell.

    Finds the dominant |g| feature in ``[search_min, search_max]`` and walks
    left until |g| drops below ``max(floor, fraction * peak)``. Peak-relative
    (not absolute-threshold) because both the physical peak and the sub-r0
    truncation ripples scale with the fitted amplitude — on missing-low-Q
    data the ripples can reach O(peak/3), so no fixed threshold separates
    them, while the dominant coordination shell still towers above them. |g|
    is used because Faber-Ziman totals of negative-b compositions (e.g.
    Mn3Sn) can have an *inverted* first shell. Returns None when no feature
    exceeds ``floor`` or the flank never falls below the level inside the
    search range (feature not separable from the ripple field).
    """
    r = np.asarray(r, dtype=float)
    g = np.abs(np.asarray(g, dtype=float))
    selection = np.where((r >= search_min) & (r <= search_max))[0]
    if selection.size < 3:
        return None
    peak_index = int(selection[np.argmax(g[selection])])
    peak = g[peak_index]
    if peak < floor:
        return None
    level = max(floor, fraction * peak)
    index = peak_index
    while index > selection[0] and g[index] > level:
        index -= 1
    if g[index] > level:
        return None  # never dropped below the level inside the search range
    return float(r[index + 1])


_HUBER_C = 1.345  # 95% Gaussian efficiency


def _huber_weights(residuals: np.ndarray) -> np.ndarray:
    """Huber IRLS weights with a MAD scale (1 inside the core, down-weighted tails)."""
    med = np.median(residuals)
    scale = 1.4826 * np.median(np.abs(residuals - med))
    if scale <= 1e-14:
        return np.ones_like(residuals)
    z = np.abs(residuals) / scale
    return np.minimum(1.0, _HUBER_C / np.maximum(z, 1e-14))


def _solve_affine(
    q: np.ndarray,
    sq: np.ndarray,
    delta_sq: np.ndarray,
    r: np.ndarray,
    tail: np.ndarray,
    window: np.ndarray,
    config: ScalingConfig,
    sigma: np.ndarray | None = None,
    level: float | None = None,
) -> tuple[float, float]:
    """Closed-form (IRLS-)LSQ for ``S_eff = a*sq + b - delta_sq`` against C1 and C2.

    C1 rows: ``a*sq + b - delta_sq - 1 = 0`` on the tail window, optionally
    1/sigma-weighted, optionally with a slope-nuisance column ``m*(Q - Qbar)``
    that absorbs a linear tail drift (x-ray f(Q) residuals, Placzek leftovers)
    so it cannot bias the level.
    C2 rows: ``g_eff(r) = 0`` on the low-r window, where the transform of the
    affine ``S_eff`` decomposes into precomputed basis transforms.

    With ``config.robust`` (default), a short Huber IRLS loop re-weights rows
    per block (C1 and C2 scaled by their own MAD) so residual Bragg spikes or
    ripple bursts cannot drag the closed-form solution.
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
        coef_s0, const = low_q_correction_basis(
            q, rw, lorch=config.lorch, s0_target=config.effective_s0_target
        )
    else:
        coef_s0 = const = np.zeros_like(rw)
    denom = 4.0 * np.pi * config.rho0 * rw
    w2 = np.sqrt(config.c2_weight)

    col_a_c1 = sq[tail].copy()
    col_b_c1 = np.ones(int(tail.sum()))
    rhs_c1 = 1.0 + delta_sq[tail]
    q_tail = q[tail]
    col_m_c1 = q_tail - q_tail.mean()  # slope-nuisance basis (C1 rows only)
    if sigma is not None:
        w_sig = 1.0 / np.clip(sigma[tail], 1e-12, None)
        w_sig = w_sig / w_sig.mean()
        col_a_c1 = col_a_c1 * w_sig
        col_b_c1 = col_b_c1 * w_sig
        col_m_c1 = col_m_c1 * w_sig
        rhs_c1 = rhs_c1 * w_sig

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

    n1 = col_a_c1.size
    stacked_a = np.concatenate([col_a_c1, col_a_c2])
    stacked_b = np.concatenate([col_b_c1, col_b_c2])
    rhs = np.concatenate([rhs_c1, rhs_c2])

    if level is not None:
        # Sweep-anchored mode: the offset is tied to the measured asymptote,
        # b = 1 - a*level, leaving a single amplitude dof (the user's
        # "shift by the level, then scale" decomposition).
        columns = [stacked_a - level * stacked_b]
        rhs = rhs - stacked_b
    else:
        columns = [stacked_a]
        if config.fit_offset:
            columns.append(stacked_b)
    if config.c1_slope_nuisance:
        columns.append(np.concatenate([col_m_c1, np.zeros_like(col_a_c2)]))
    design = np.column_stack(columns)

    solution, *_ = np.linalg.lstsq(design, rhs, rcond=None)
    if config.robust:
        for _ in range(3):
            residuals = design @ solution - rhs
            weights = np.ones_like(rhs)
            weights[:n1] = _huber_weights(residuals[:n1])
            if rhs.size - n1 >= 4:
                weights[n1:] = _huber_weights(residuals[n1:])
            solution, *_ = np.linalg.lstsq(
                design * weights[:, np.newaxis], rhs * weights, rcond=None
            )
    a = float(solution[0])
    if level is not None:
        b = 1.0 - a * float(level)
    elif config.fit_offset:
        b = float(solution[1])
    else:
        b = 0.0
    return a, b


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
        s0_target=config.effective_s0_target,
    )


def _low_r_rms(r: np.ndarray, g_filtered: np.ndarray, config: ScalingConfig) -> float:
    lo, hi = config.r_fit_window
    window = (r >= lo) & (r <= hi)
    return float(np.sqrt(np.mean(g_filtered[window] ** 2)))


def amplitude_from_fz_limit(
    q: np.ndarray,
    sq: np.ndarray,
    level: float,
    config: ScalingConfig,
    *,
    fit_width: float = 1.0,
) -> float | None:
    """Independent amplitude estimate from the Q->0 Faber-Ziman limit.

    With the level-anchored model ``S_corr = a (S_meas - level) + 1``, Keen
    Eq. 21 fixes ``S_corr(0) = 1 - <b^2>/<b>^2``, so
    ``a_fz = (s0_target - 1) / (S_meas(0) - level)`` where ``S_meas(0)`` is a
    robust linear extrapolation of the first ``fit_width`` of measured data.
    Requires ``config.b_sq_avg``; returns None when unavailable or the
    extrapolation is degenerate. The caller should treat long extrapolations
    (Qmin >> fit_width) and Bragg-contaminated low-Q regions with suspicion —
    compare against the density-limit amplitude (``diagnostics_summary``'s
    concordance) rather than trusting either alone.
    """
    if config.b_sq_avg is None:
        return None
    s0_target = 1.0 - config.b_sq_avg / config.b_avg_sq
    head = q <= q[0] + fit_width
    if head.sum() < 8:
        return None
    qc = q[head] - q[head].mean()
    design = np.column_stack([np.ones_like(qc), qc])
    weights = np.ones_like(qc)
    solution = np.array([np.median(sq[head]), 0.0])
    for _ in range(4):
        solution, *_ = np.linalg.lstsq(
            design * weights[:, np.newaxis], sq[head] * weights, rcond=None
        )
        residuals = design @ solution - sq[head]
        weights = _huber_weights(residuals)
    s_meas_0 = float(solution[0] - solution[1] * q[head].mean())  # value at Q = 0
    denom = s_meas_0 - level
    if abs(denom) < 1e-9:
        return None
    return float((s0_target - 1.0) / denom)


def _despike_mask(sq: np.ndarray, window: int, nsigma: float) -> np.ndarray:
    """Keep-mask dropping narrow outliers vs a rolling median (glitch removal)."""
    from numpy.lib.stride_tricks import sliding_window_view

    pad = window // 2
    padded = np.pad(sq, pad, mode="edge")
    median = np.median(sliding_window_view(padded, window), axis=1)
    residual = sq - median
    mad = 1.4826 * np.median(np.abs(residual))
    return np.abs(residual) <= nsigma * max(mad, 1e-12)


def crop_sq(
    q: np.ndarray,
    sq: np.ndarray,
    config: ScalingConfig,
    sigma: np.ndarray | None = None,
):
    """Crop to ``(0, qmax]`` ∩ ``[qmin, qmax]``, dropping non-finite rows.

    ``Q <= 0`` rows are dropped unconditionally: the S(Q) conversions divide
    by Q, and the analytic low-Q correction already models the omitted
    ``[0, Qmin]`` range, so a Q = 0 point carries no usable information.
    Returns ``(q, sq, sigma)`` with ``sigma`` cropped alongside (or ``None``).
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
    q, sq = q[keep], sq[keep]
    if sigma is not None:
        sigma = np.asarray(sigma, dtype=float)[keep]
    if config.despike:
        keep2 = _despike_mask(sq, config.despike_window, config.despike_nsigma)
        q, sq = q[keep2], sq[keep2]
        if sigma is not None:
            sigma = sigma[keep2]
    return q, sq, sigma


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
    sweep: LevelSweepResult | None = None,
    a_fz: float | None = None,
) -> ScalingResult:
    """Run the fixed-(a, b) pipeline: scale, filter, convert, (optionally) enforce.

    This is the manual-mode entry point (classic stog parity) and the final
    stage of :func:`autoscale`.
    """
    n_despiked = 0
    if config.despike:
        q_ref, _, _ = crop_sq(q, sq, replace(config, despike=False))
        n_despiked = int(q_ref.size)
    q, sq, _ = crop_sq(q, sq, config)
    if config.despike:
        n_despiked -= int(q.size)
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
                "low_q_correction", "c2_weight", "c2_bins", "robust",
                "c1_mode", "amplitude_criterion", "s0_target",
                "c1_slope_nuisance", "despike", "despike_window",
                "despike_nsigma", "max_iter", "tol", "enforce_cutoff",
            )
        },
        "q_tail_window": [config.q_tail_min, config.qmax],
        "r_fit_window": list(config.r_fit_window),
        "n_q_points": int(q.size),
        "n_despiked": n_despiked,
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
        sweep=sweep,
        a_fz=a_fz,
        provenance=provenance,
    )


def autoscale(
    q: np.ndarray,
    sq: np.ndarray,
    config: ScalingConfig,
    sigma: np.ndarray | None = None,
) -> ScalingResult:
    """Automatically determine ``(a, b)`` and run the full pipeline.

    Self-consistent loop: fit ``(a, b)`` with the Fourier-filter subtraction
    term held fixed, re-run the filter on the corrected data, repeat until the
    parameters converge (``tol``) or ``max_iter`` is reached. ``sigma``
    (per-point uncertainties, e.g. the data file's third column) weights the
    high-Q C1 rows by 1/sigma.

    When the caller pins neither ``r0`` nor ``r_fit_max``, a second refinement
    pass runs: the first pass's filtered g(r) yields a data-derived closest
    approach (:func:`detect_first_peak_onset`), and if that moves the low-r
    window materially the fit is redone with the detected ``r0``. The result's
    provenance then carries ``r0_detected`` / ``window_refined`` — this is what
    lets composition + Q-window be the only required inputs.
    """
    result = _autoscale_pass(q, sq, config, sigma)
    onset = detect_first_peak_onset(
        result.r, result.g_filtered, config.qmax,
        search_min=config.r_cutoff + 0.3,
    )
    if onset is not None:
        result.provenance["r0_detected"] = float(onset)
    if (
        onset is not None
        and config.r0 is None
        and config.r_fit_max is None
        and onset - 0.25 > (config.r_fit_min if config.r_fit_min is not None else config.r_cutoff + 0.2)
        and abs((onset - 0.25) - config.r_fit_window[1]) > 0.05
    ):
        refined = replace(config, r0=float(onset))
        refined_result = _autoscale_pass(q, sq, refined, sigma)
        refined_result.provenance["r0_detected"] = float(onset)
        refined_result.provenance["window_refined"] = True
        return refined_result
    return result


def _autoscale_pass(
    q: np.ndarray,
    sq: np.ndarray,
    config: ScalingConfig,
    sigma: np.ndarray | None = None,
) -> ScalingResult:
    q, sq, sigma = crop_sq(q, sq, config, sigma)
    r = config.r_grid
    tail, window = _fit_windows(q, r, config)

    sweep: LevelSweepResult | None = None
    level: float | None = None
    if config.c1_mode == "sweep":
        sweep = level_sweep(q, sq)
        if sweep.asymptote_found:
            level = sweep.level

    if config.amplitude_criterion == "fz":
        # "Subtract the level, pin Q->0 on the Faber-Ziman limit, shift the
        # level back to 1": a single closed-form amplitude, no loop — the
        # criterion does not involve the Fourier filter.
        if level is None:
            raise ValueError(
                "amplitude_criterion='fz': level_sweep found no statistically "
                "flat high-Q window, so there is no measured level to anchor; "
                "inspect the tail or use the density-limit fit"
            )
        a_fz = amplitude_from_fz_limit(q, sq, level, config)
        if a_fz is None or not np.isfinite(a_fz) or a_fz <= 0:
            raise ValueError(
                "amplitude_criterion='fz': the Q->0 extrapolation is "
                f"degenerate (a_fz={a_fz}); the data head cannot support the "
                "Faber-Ziman limit"
            )
        result = scale_pipeline(
            q,
            sq,
            config,
            float(a_fz),
            float(1.0 - a_fz * level),
            converged=True,
            iterations=0,
            sweep=sweep,
            a_fz=float(a_fz),
        )
        result.provenance["mode"] = "auto"
        result.provenance["c1_mode_effective"] = "sweep"
        result.provenance["level_sweep"] = _sweep_provenance(sweep)
        return result

    delta_sq = np.zeros_like(q)
    a_prev = b_prev = np.inf
    a = 1.0
    b = 0.0
    history: list[tuple[float, float, float]] = []
    converged = False
    iterations = 0

    for iterations in range(1, config.max_iter + 1):
        a, b = _solve_affine(q, sq, delta_sq, r, tail, window, config, sigma, level)
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

    a_fz = None
    if level is not None:
        a_fz = amplitude_from_fz_limit(q, sq, level, config)

    result = scale_pipeline(
        q,
        sq,
        config,
        a,
        b,
        converged=converged,
        iterations=iterations,
        history=tuple(history),
        sweep=sweep,
        a_fz=a_fz,
    )
    result.provenance["mode"] = "auto"
    result.provenance["c1_mode_effective"] = "sweep" if level is not None else "joint"
    if sweep is not None:
        result.provenance["level_sweep"] = _sweep_provenance(sweep)
    return result


def _sweep_provenance(sweep: LevelSweepResult) -> dict[str, Any]:
    return {
        "level": sweep.level,
        "level_uncertainty": sweep.level_uncertainty,
        "q_window": [sweep.q_lo, sweep.q_hi],
        "asymptote_found": sweep.asymptote_found,
        "n_admissible": sweep.n_admissible,
    }


def diagnostics_summary(result: ScalingResult, config: ScalingConfig) -> dict[str, Any]:
    """Compact, human-readable fit-quality summary (JSON-friendly).

    Coefficients and the low-r window are read from the result's provenance
    when present, so a window-refined result (see :func:`autoscale`) reports
    against the window it was actually fitted with, not the caller's original
    config.
    """
    effective = result.provenance.get("config", {})
    b_avg_sq = float(effective.get("b_avg_sq", config.b_avg_sq))
    rho0 = float(effective.get("rho0", config.rho0))
    b_sq_avg = effective.get("b_sq_avg", config.b_sq_avg)
    amplitude_criterion = effective.get("amplitude_criterion", config.amplitude_criterion)
    lo, hi = result.provenance.get("r_fit_window", config.r_fit_window)
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
        "r_fit_window": [float(lo), float(hi)],
        "gk_low_r_theory": -b_avg_sq,
        "d_r_low_r_slope_theory": -4.0 * np.pi * rho0 * b_avg_sq,
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
    if result.sweep is not None:
        summary["level"] = result.sweep.level
        summary["level_uncertainty"] = result.sweep.level_uncertainty
        summary["level_window"] = [result.sweep.q_lo, result.sweep.q_hi]
        summary["asymptote_found"] = result.sweep.asymptote_found
    if "r0_detected" in result.provenance:
        summary["r0_detected"] = result.provenance["r0_detected"]
        summary["window_refined"] = bool(result.provenance.get("window_refined", False))
    if result.a_fz is not None:
        summary["a_fz"] = result.a_fz
        if amplitude_criterion != "fz":
            # Concordance of the two independent amplitude criteria: the
            # density-limit amplitude (result.a) vs the Q->0 Faber-Ziman-limit
            # amplitude. Agreement is evidence the absolute scale is
            # trustworthy; disagreement quantifies how much the data cannot
            # decide it. (In fz mode a IS a_fz, so the independent check is
            # the density-limit residual reported above instead.)
            summary["amplitude_concordance"] = float(result.a_fz / result.a)
            summary["amplitudes_concordant"] = bool(
                abs(result.a_fz / result.a - 1.0) < 0.1
            )
    if b_sq_avg is not None:
        # Keen Eq. 14 diagnostic: FK(Q->0) -> -<b^2>. Report the lowest-Q value
        # actually measured for comparison (data rarely reach Q ~ 0).
        summary["fk_qmin"] = float(result.fk[0])
        summary["fk_q0_theory"] = -float(b_sq_avg)
    return summary
