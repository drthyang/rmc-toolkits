# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Fourier termination ripples: what the Q range of a run's data puts into its G(r).

A G(r) is the sine transform of F(Q) measured from Q_min to Q_max. Cutting the
transform at Q_max convolves every peak with the termination kernel, a sinc of
period 2π/Q_max: the peak broadens and oscillating ripples spread on both sides
of it. RMC is a minimiser: fitted to a G(r) that carries those ripples, with a
model G(r) that does not, the box makes structure that reproduces them. This
module says, per dataset of a run, how large that risk is and what the ripples
are, so a small peak, a shoulder or a split in a g(r) is not read as structure
when it is a ripple:

* the fit space and how RMCProfile makes the model function it compares
  (established on RMCProfile's own files, docs/algorithms/run-dashboard.md,
  "Termination ripples"):

  - an X-ray F(Q) fitted in real space (``REAL_SPACE_FIT``): the data's and
    the model's G(r) are the same rectangle-rule transform of their F(Q) over
    the run's Q grid (the experiment to 1.5e-14, the model to 1.6e-5 of the
    largest |G|), so the model carries the data's ripples: ``matched``;
  - a neutron G(r) (``NEUTRON_REAL_SPACE_DATA``) without ``CONVOLVE``: the
    model is the configuration's histogram Σ w_ij [g_ij(r) − 1], to 6e-15 of
    the largest |G| on three runs — no truncation at all: ``high``;
  - the same with ``CONVOLVE :: Q_c``: the model is the histogram convolved
    with the sinc of Q_c in D(r) = r·G(r) space, to 1.8e-14, and set to the
    baseline below ``R_CUTOFF``. Matched only when Q_c is the data's Q_max
    and the data were transformed without a window; the run's files rarely
    say, so ``unknown`` unless they do;
  - an X-ray G(r) (``XRAY_REAL_SPACE_DATA``): the histogram transformed to Q
    on its QMIN–QMAX grid and back by the rectangle rule (to 8e-6 with
    QMIN = 0): matched only when its QMAX is the data's, ``unknown`` otherwise;
  - a dataset fitted in Q only: the box is not asked to make ripples,
    ``none`` — though a G(r) transformed from it still shows them.

* the transform: Q_min, Q_max, ΔQ, the window (Lorch where the files say so),
  the r-resolution — the Nyquist step π/Q_max and the kernel's full width at
  half maximum (3.79/Q_max for no window, computed on the run's own grid and
  window), the ripple period 2π/Q_max;
* the strong peaks, each a Gaussian A·[e^{−(r−r₀)²/2σ²} − e^{−(r+r₀)²/2σ²}] in
  D-space (``d = r·G`` for a G(r), the function itself for a D(r)-like one)
  fitted THROUGH the transform (its sine transform is A√(2π)σ e^{−Q²σ²/2}
  sin Q r₀, exactly, so the fitted model is the truncated peak), and its
  termination ripples
  Δ(r) = T[d](r) − T_{[Q_min,∞)}[d](r), the transform over the run's grid and
  window less the same transform with no upper cut: positions and amplitudes;
* the low-r level: below the first interatomic distance G(r) is its baseline
  (a constant, −(Σcb)² or −1 normalised, in G(r); the line −4πρ₀(Σcb)² r in
  D(r)); the rms about it measures ripple and noise, set against the ripples
  the strong peaks predict there. A low-r region set to the baseline (stog's
  and Auto StoG's RMC files) or Fourier-filtered (stog's filter, Auto StoG's
  r_cutoff) is detected: it no longer measures the ripples.
* features (:func:`classify_features`): a small peak, shoulder or split at r
  is ``ripple-consistent`` when a predicted ripple of the same sign lies within
  a quarter period and the feature is no larger than RIPPLE_SIZE_FACTOR times
  it; a split or a shoulder closer than the kernel's FWHM is below the
  resolution; across datasets with different Q_max a ripple moves and a real
  feature stays.

``web_app/frontend/src/ripples.js`` is a line-for-line port (parity golden
``ripples_fixture.json``, ``tests/generate_ripples_fixture.py``).
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Callable, Optional

import numpy as np

from .ft_uncertainty import FT_FILE_RE, RMC_FT_TOLERANCE, reproduce_transform, uniform_step

# ---------------------------------------------------------------------------
# Rules (each one named in the payload's `rules`)
# ---------------------------------------------------------------------------

#: A feature within this factor of the predicted ripple's amplitude there,
#: either way, is ripple-consistent (the strong peaks are Gaussians fitted
#: through the transform: an approximation of their true shape).
RIPPLE_SIZE_FACTOR = 1.5

#: A predicted ripple matches a feature's position within this fraction of
#: the ripple period 2π/Q_max (a quarter period: the ripple's own half-width).
RIPPLE_POSITION_FRACTION = 0.25

#: A predicted ripple is negligible at a feature below this fraction of the
#: feature's size.
NEGLIGIBLE_RIPPLE_FRACTION = 0.1

#: Strong-peak candidates: local extrema of d(r) whose prominence is at least
#: this fraction of the largest one; at most MAX_PEAK_CANDIDATES of them.
STRONG_PEAK_FRACTION = 0.1
MAX_PEAK_CANDIDATES = 24

#: Peaks kept for the prediction: at most MAX_STRONG_PEAKS, those whose
#: ripple scale |A|σ e^{−Q_max²σ²/2} is at least STRONG_RIPPLE_FRACTION of the
#: largest.
MAX_STRONG_PEAKS = 12
STRONG_RIPPLE_FRACTION = 0.02

#: A peak's fit window: the points around its extremum within this fraction of
#: its prominence (at least PEAK_WINDOW_MIN_POINTS on each side).
PEAK_WINDOW_FRACTION = 0.3
PEAK_WINDOW_MIN_POINTS = 2
PEAK_WINDOW_MAX_PERIODS = 1.0

#: A candidate is a peak only when its fitted height, the larger peaks'
#: truncated models taken off, is at least this fraction of its prominence:
#: otherwise those peaks' ripples explain it.
RIPPLE_EXPLAINED_FRACTION = 0.5

#: Backfitting sweeps after the first fits (each peak refitted with the
#: others' truncated models taken off).
BACKFIT_SWEEPS = 2

#: A peak's ripples are computed within this many ripple periods of it.
RIPPLE_REACH_PERIODS = 10

#: The ripples listed per peak on each side.
RIPPLES_LISTED = 2

#: The low-r window must span at least this many ripple periods.
LOW_R_MIN_PERIODS = 1.0

#: The low-r region is "set to the baseline" where it equals its baseline to
#: this fraction of the function's largest |value|.
ENFORCED_RTOL = 1e-9

#: The low-r level is called filtered when its rms is below this fraction of
#: the ripples the strong peaks predict there.
FILTERED_RATIO = 0.3

#: A feature is present in a dataset when a same-sign extremum of its
#: deviation lies within the position tolerance and is at least this fraction
#: of the feature's relative size in the dataset where it is largest.
PRESENT_FRACTION = 0.3

#: The background under a feature: the box average of the data less the
#: strong peaks' truncated models over this many ripple periods.
BACKGROUND_PERIODS = 2.0

#: Automatically found features: at most MAX_FEATURES, at least
#: FEATURE_MIN_FRACTION of the largest strong peak's height.
MAX_FEATURES = 8
FEATURE_MIN_FRACTION = 0.01

#: σ of a fitted peak is kept within [SIGMA_MIN_QMAX / Q_max, SIGMA_MAX] Å.
SIGMA_MIN_QMAX = 0.1
SIGMA_MAX = 1.0

#: Levenberg–Marquardt on (r₀, ln σ) to a step of FIT_STEP_TOL or a relative
#: cost decrease of FIT_COST_RTOL, then FIT_POLISH_STEPS Gauss–Newton steps to
#: the stationary point.
FIT_MAX_ITERATIONS = 30
FIT_STEP_TOL = 1e-7
FIT_COST_RTOL = 1e-10
FIT_POLISH_STEPS = 4

#: The grid step assumed for a transform whose own step the files do not give.
DEFAULT_DQ = 0.01

#: The input data file maps onto the fitted file (y = s·y_input + c) within
#: this fraction of the fitted file's largest |value|.
INPUT_MAP_TOLERANCE = 1e-5

RISKS = ("high", "matched", "unknown", "none")
RULES = {
    "rippleSizeFactor": RIPPLE_SIZE_FACTOR,
    "ripplePositionFraction": RIPPLE_POSITION_FRACTION,
    "strongPeakFraction": STRONG_PEAK_FRACTION,
    "maxStrongPeaks": MAX_STRONG_PEAKS,
    "lowRMinPeriods": LOW_R_MIN_PERIODS,
    "filteredRatio": FILTERED_RATIO,
    "presentFraction": PRESENT_FRACTION,
}

_SQRT_2PI = math.sqrt(2.0 * math.pi)


# ---------------------------------------------------------------------------
# The transform's grid and kernel
# ---------------------------------------------------------------------------


def lorch(q: np.ndarray, qmax: float) -> np.ndarray:
    """Lorch's window sin(πQ/Q_max)/(πQ/Q_max) (transforms.lorch_window)."""
    x = math.pi * np.asarray(q, dtype=float) / float(qmax)
    out = np.ones_like(x)
    nonzero = x != 0
    out[nonzero] = np.sin(x[nonzero]) / x[nonzero]
    return out


def q_grid(q_min: float, dq: float, n: int) -> np.ndarray:
    """The transform's Q points q_min + k·ΔQ, k = 0 … n − 1."""
    return float(q_min) + float(dq) * np.arange(int(n), dtype=float)


def transform_weights(q: np.ndarray, dq: float, rule: str, window: str) -> np.ndarray:
    """w_k of T[d](r) = (2/π) Σ_k w_k S(Q_k) sin(Q_k r): ΔQ (rectangle) or the trapezoid's, times the window."""
    weights = np.full(q.size, float(dq))
    if rule == "trapezoid" and q.size > 1:
        weights[0] *= 0.5
        weights[-1] *= 0.5
    if window == "lorch":
        weights = weights * lorch(q, float(q[-1]))
    return weights


def kernel(q: np.ndarray, weights: np.ndarray, x: float) -> float:
    """The termination kernel K(x) = Σ_k w_k cos(Q_k x): a δ peak's shape in r."""
    return float(np.sum(weights * np.cos(q * float(x))))


def _first_crossing(q: np.ndarray, weights: np.ndarray, level: float) -> float:
    """The first x > 0 where K(x) falls to ``level``, by stepping then bisection."""
    qmax = float(q[-1])
    step = 0.05 / qmax
    lo, hi = 0.0, step
    while kernel(q, weights, hi) > level and hi < 100.0 / qmax:
        lo, hi = hi, hi + step
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if kernel(q, weights, mid) > level:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def resolution(q: np.ndarray, weights: np.ndarray) -> dict:
    """The r-resolution of a transform: the Nyquist step, the kernel's FWHM and main lobe, the ripple period.

    The kernel is the run's own (its grid and window): the FWHM is twice the
    first x where K(x) = K(0)/2 (3.79/Q_max for a rectangle from 0), the main
    lobe its first zero (π/Q_max for a rectangle; wider with Lorch), where
    the ripples begin.
    """
    qmax = float(q[-1])
    return {
        "nyquistStep": math.pi / qmax,
        "fwhm": 2.0 * _first_crossing(q, weights, 0.5 * kernel(q, weights, 0.0)),
        "mainLobe": _first_crossing(q, weights, 0.0),
        "ripplePeriod": 2.0 * math.pi / qmax,
    }


# ---------------------------------------------------------------------------
# A Gaussian peak through the transform
# ---------------------------------------------------------------------------


def gaussian_d(r: np.ndarray, peak: dict) -> np.ndarray:
    """The peak in D-space: A[e^{−(r−r₀)²/2σ²} − e^{−(r+r₀)²/2σ²}]."""
    r = np.asarray(r, dtype=float)
    a, r0, s = peak["amplitude"], peak["r"], peak["sigma"]
    return a * (np.exp(-((r - r0) ** 2) / (2.0 * s * s)) - np.exp(-((r + r0) ** 2) / (2.0 * s * s)))


def gaussian_sine(q: np.ndarray, peak: dict) -> np.ndarray:
    """Its sine transform ∫₀^∞ d(r) sin(Qr) dr = A√(2π)σ e^{−Q²σ²/2} sin(Q r₀), exactly."""
    q = np.asarray(q, dtype=float)
    s = peak["sigma"]
    return peak["amplitude"] * _SQRT_2PI * s * np.exp(-0.5 * q * q * s * s) * np.sin(q * peak["r"])


def apply_transform(q: np.ndarray, weights: np.ndarray, spectrum: np.ndarray, r: np.ndarray) -> np.ndarray:
    """(2/π) Σ_k w_k S_k sin(Q_k r) at every r."""
    coefficients = (2.0 / math.pi) * weights * spectrum
    return np.sin(np.outer(np.asarray(r, dtype=float), q)) @ coefficients


def ripple_pattern(q: np.ndarray, weights: np.ndarray, dq: float, rule: str, peak: dict, r: np.ndarray) -> np.ndarray:
    """Δ(r) = T[d](r) − T_{[Q_min,∞)}[d](r): the peak's termination ripples (and window broadening).

    The transform over the whole grid Q_min + kΔQ, k ∈ ℤ, Q > 0, at weight ΔQ,
    is d itself (the rectangle rule of a smooth integrand that vanishes at
    Q = 0, to O(ΔQ³), below the aliasing distance 2π/ΔQ). So the same
    transform with no upper cut is d less its part below Q_min, at the fit's
    own weight at Q_min (half of ΔQ for the trapezoid): no sum to infinity.
    """
    r = np.asarray(r, dtype=float)
    truncated = apply_transform(q, weights, gaussian_sine(q, peak), r)
    below_count = int(math.floor(float(q[0]) / dq + 1e-9))
    below = float(q[0]) - dq * np.arange(1, below_count + 1, dtype=float)
    below = below[below > 0]
    low = apply_transform(below, np.full(below.size, dq), gaussian_sine(below, peak), r) if below.size else 0.0
    first_weight = dq * (0.5 if rule == "trapezoid" and q.size > 1 else 1.0)
    edge = (2.0 / math.pi) * (dq - first_weight) * gaussian_sine(q[:1], peak)[0] * np.sin(q[0] * r)
    return truncated - gaussian_d(r, peak) + low + edge


# ---------------------------------------------------------------------------
# Strong peaks, fitted through the transform
# ---------------------------------------------------------------------------


def local_extrema(y: np.ndarray) -> list[int]:
    """Interior indices i where y_i is a strict local maximum or minimum (plateaus to their first point)."""
    out = []
    n = y.size
    for i in range(1, n - 1):
        left, here, right = y[i - 1], y[i], y[i + 1]
        if (here > left and here >= right) or (here < left and here <= right):
            out.append(i)
    return out


def _prominences(y: np.ndarray, extrema: list[int]) -> list[tuple[int, float, float]]:
    """(index, prominence, base) of each extremum against the nearest opposite extrema (or the ends)."""
    out = []
    for position, i in enumerate(extrema):
        is_max = y[i] > y[i - 1]
        left = right = None
        for j in range(position - 1, -1, -1):
            k = extrema[j]
            if (y[k] < y[k - 1]) == is_max:  # opposite type
                left = k
                break
        for j in range(position + 1, len(extrema)):
            k = extrema[j]
            if (y[k] < y[k - 1]) == is_max:
                right = k
                break
        left_value = y[left] if left is not None else y[0]
        right_value = y[right] if right is not None else y[-1]
        base = 0.5 * (left_value + right_value)
        out.append((i, abs(float(y[i] - base)), float(base)))
    return out


def _fit_window(y: np.ndarray, centre: int, base: float) -> tuple[int, int]:
    """[lo, hi) around the extremum where |y − base| ≥ PEAK_WINDOW_FRACTION of its prominence."""
    height = y[centre] - base
    limit = PEAK_WINDOW_FRACTION * abs(height)
    lo = centre
    while lo - 1 >= 0 and (y[lo - 1] - base) * math.copysign(1.0, height) >= limit:
        lo -= 1
    hi = centre
    while hi + 1 < y.size and (y[hi + 1] - base) * math.copysign(1.0, height) >= limit:
        hi += 1
    lo = max(0, min(lo, centre - PEAK_WINDOW_MIN_POINTS))
    hi = min(y.size - 1, max(hi, centre + PEAK_WINDOW_MIN_POINTS))
    return lo, hi + 1


def _solve3(m: list[list[float]], v: list[float]) -> Optional[list[float]]:
    """The 3 × 3 system m x = v by Cramer's rule, or None when singular."""
    def det(a):
        return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
                - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
    d = det(m)
    if d == 0 or not math.isfinite(d):
        return None
    out = []
    for column in range(3):
        a = [[v[row] if k == column else m[row][k] for k in range(3)] for row in range(3)]
        out.append(det(a) / d)
    return out


def _projected(sines: np.ndarray, q: np.ndarray, weights: np.ndarray, target: np.ndarray, r: np.ndarray,
               r0: float, sigma: float):
    """For (r₀, σ): the unit peak through the transform φ, then A and a line b₀ + b₁(r − r₀) by least squares.

    Returns (A, b₀, b₁, residual) or None.
    """
    unit = {"amplitude": 1.0, "r": r0, "sigma": sigma}
    phi = sines @ ((2.0 / math.pi) * weights * gaussian_sine(q, unit))
    x = r - r0
    columns = (phi, np.ones_like(phi), x)
    m = [[float(np.dot(a, b)) for b in columns] for a in columns]
    v = [float(np.dot(a, target)) for a in columns]
    solution = _solve3(m, v)
    if solution is None:
        return None
    amplitude, b0, b1 = solution
    residual = target - amplitude * phi - b0 - b1 * x
    return amplitude, b0, b1, residual


def fit_peak(r: np.ndarray, target: np.ndarray, q: np.ndarray, weights: np.ndarray, r0: float, sigma: float,
             bounds: tuple[float, float]) -> Optional[dict]:
    """A Gaussian through the transform, fitted to ``target`` on ``r`` (a window): A, r₀, σ and a line.

    Levenberg–Marquardt on (r₀, ln σ) with A and the line b₀ + b₁(r − r₀)
    solved by least squares at every step (variable projection); central
    differences for the Jacobian. Returns ``{amplitude, r, sigma,
    background, slope, rms}`` or None.
    """
    qmax = float(q[-1])
    s_lo, s_hi = math.log(SIGMA_MIN_QMAX / qmax), math.log(SIGMA_MAX)
    sines = np.sin(np.outer(r, q))
    theta = [float(r0), min(max(math.log(sigma), s_lo), s_hi)]
    r_lo, r_hi = bounds

    def clamp(values):
        return [min(max(values[0], r_lo), r_hi), min(max(values[1], s_lo), s_hi)]

    def evaluate(values):
        return _projected(sines, q, weights, target, r, values[0], math.exp(values[1]))

    current = evaluate(theta)
    if current is None:
        return None
    cost = float(np.dot(current[3], current[3]))
    lam = 1e-3
    steps = (1e-5, 1e-5)
    for _ in range(FIT_MAX_ITERATIONS):
        columns = []
        for j in range(2):
            plus = list(theta)
            minus = list(theta)
            plus[j] += steps[j]
            minus[j] -= steps[j]
            fp, fm = evaluate(plus), evaluate(minus)
            if fp is None or fm is None:
                return None
            columns.append((fp[3] - fm[3]) / (2.0 * steps[j]))
        j0, j1 = columns
        f = current[3]
        h00, h01, h11 = float(np.dot(j0, j0)), float(np.dot(j0, j1)), float(np.dot(j1, j1))
        g0, g1 = float(np.dot(j0, f)), float(np.dot(j1, f))
        accepted = False
        while lam < 1e12:
            a00, a11 = h00 * (1.0 + lam), h11 * (1.0 + lam)
            det = a00 * a11 - h01 * h01
            if det <= 0 or not math.isfinite(det):
                lam *= 4.0
                continue
            d0 = -(a11 * g0 - h01 * g1) / det
            d1 = -(a00 * g1 - h01 * g0) / det
            trial = clamp([theta[0] + d0, theta[1] + d1])
            result = evaluate(trial)
            if result is not None:
                trial_cost = float(np.dot(result[3], result[3]))
                if trial_cost < cost:
                    moved = max(abs(trial[0] - theta[0]), abs(trial[1] - theta[1]))
                    theta, current, previous, cost = trial, result, cost, trial_cost
                    lam = max(lam / 3.0, 1e-12)
                    accepted = True
                    break
            lam *= 4.0
        if not accepted or moved < FIT_STEP_TOL or previous - cost <= FIT_COST_RTOL * previous:
            break
    # Polish: Gauss–Newton steps to the stationary point (no cost test, which
    # round-off decides near the minimum), so two engines that sum in another
    # order land on the same parameters, not merely on the same cost.
    for _ in range(FIT_POLISH_STEPS):
        columns = []
        for j in range(2):
            plus = list(theta)
            minus = list(theta)
            plus[j] += steps[j]
            minus[j] -= steps[j]
            fp, fm = evaluate(plus), evaluate(minus)
            if fp is None or fm is None:
                break
            columns.append((fp[3] - fm[3]) / (2.0 * steps[j]))
        if len(columns) < 2:
            break
        j0, j1 = columns
        f = current[3]
        h00, h01, h11 = float(np.dot(j0, j0)), float(np.dot(j0, j1)), float(np.dot(j1, j1))
        g0, g1 = float(np.dot(j0, f)), float(np.dot(j1, f))
        det = h00 * h11 - h01 * h01
        if det <= 0 or not math.isfinite(det):
            break
        trial = clamp([theta[0] - (h11 * g0 - h01 * g1) / det, theta[1] - (h00 * g1 - h01 * g0) / det])
        result = evaluate(trial)
        if result is None or float(np.dot(result[3], result[3])) > cost * (1.0 + 1e-9):
            break
        theta, current, cost = trial, result, float(np.dot(result[3], result[3]))
    amplitude, background, slope, residual = current
    return {
        "amplitude": float(amplitude),
        "r": float(theta[0]),
        "sigma": float(math.exp(theta[1])),
        "background": float(background),
        "slope": float(slope),
        "rms": float(math.sqrt(float(np.dot(residual, residual)) / residual.size)),
    }


def ripple_scale(peak: dict, qmax: float) -> float:
    """|A| σ e^{−Q_max²σ²/2}: the envelope of the peak's transform at Q_max, which sets its ripples."""
    s = peak["sigma"]
    return abs(peak["amplitude"]) * s * math.exp(-0.5 * qmax * qmax * s * s)


def _truncated_models(r: np.ndarray, q: np.ndarray, weights: np.ndarray, dq: float, rule: str,
                      peaks: list[dict]) -> np.ndarray:
    """Σ (d_j + Δ_j) on r: the peaks' models truncated at Q_max (not at Q_min)."""
    total = np.zeros(r.size)
    for peak in peaks:
        total += gaussian_d(r, peak) + ripple_pattern(q, weights, dq, rule, peak, r)
    return total


def _peak_like(peak: Optional[dict], sign: float, half_width: float) -> bool:
    """A fit is a peak when its amplitude has the extremum's sign and σ is within the window's half-width
    (a wider Gaussian fits the background under the window, not a peak) and below SIGMA_MAX."""
    return (peak is not None and peak["amplitude"] * sign > 0 and peak["sigma"] < SIGMA_MAX
            and peak["sigma"] <= half_width)


def strong_peaks(r: np.ndarray, d: np.ndarray, q: np.ndarray, weights: np.ndarray, dq: float, rule: str,
                 r_floor: Optional[float] = None, signs: str = "positive") -> list[dict]:
    """The strong peaks of d(r), each a Gaussian fitted through the transform, backfitted.

    Candidates: local maxima of d (``signs`` ``positive``) or extrema of either
    sign (``both``: a neutron dataset with a negative coefficient has negative
    peaks), at r ≥ r_floor, by prominence (≥ STRONG_PEAK_FRACTION of the
    largest, the MAX_PEAK_CANDIDATES largest). Taken largest first, each is
    fitted on its window (at most PEAK_WINDOW_MAX_PERIODS ripple periods each
    side) with the accepted peaks' truncated models (d + Δ) taken off, and
    accepted when it is peak-like (:func:`_peak_like`) and its fitted height is at least RIPPLE_EXPLAINED_FRACTION of its
    prominence — a ripple of a larger peak is explained by that peak's Δ and
    is not a peak. The accepted ones are fitted BACKFIT_SWEEPS times again
    with the others taken off. Returned: the MAX_STRONG_PEAKS with the largest
    ripple scale (≥ STRONG_RIPPLE_FRACTION of the largest), sorted by r.
    """
    qmax = float(q[-1])
    period = 2.0 * math.pi / qmax
    step = float(r[1] - r[0])
    scored = [entry for entry in _prominences(d, local_extrema(d))
              if (r_floor is None or r[entry[0]] >= r_floor)
              and (signs == "both" or d[entry[0]] > entry[2])]
    if not scored:
        return []
    largest = max(entry[1] for entry in scored)
    if not largest > 0:
        return []
    candidates = sorted((entry for entry in scored if entry[1] >= STRONG_PEAK_FRACTION * largest),
                        key=lambda entry: (-entry[1], entry[0]))[:MAX_PEAK_CANDIDATES]
    cap = max(PEAK_WINDOW_MIN_POINTS, int(round(PEAK_WINDOW_MAX_PERIODS * period / step)))
    accepted = []
    for index, prominence, base in candidates:
        lo, hi = _fit_window(d, index, base)
        lo, hi = max(lo, index - cap), min(hi, index + cap + 1)
        window_r = r[lo:hi]
        target = d[lo:hi] - _truncated_models(window_r, q, weights, dq, rule, [e["peak"] for e in accepted])
        half_width = max(step, 0.5 * (r[hi - 1] - r[lo]))
        sigma0 = max(half_width / 1.1774, 0.5 / qmax)
        peak = fit_peak(window_r, target, q, weights, float(r[index]), sigma0, (float(r[lo]), float(r[hi - 1])))
        sign = d[index] - base
        if not _peak_like(peak, sign, half_width):
            continue
        height = abs(peak["amplitude"] + 0.0) * abs(float(apply_transform(
            q, weights, gaussian_sine(q, {"amplitude": 1.0, "r": peak["r"], "sigma": peak["sigma"]}),
            np.array([peak["r"]]))[0]))
        if height < RIPPLE_EXPLAINED_FRACTION * prominence:
            continue
        accepted.append({"peak": peak, "window": (lo, hi), "sign": sign, "halfWidth": half_width})
    accepted.sort(key=lambda entry: entry["peak"]["r"])
    for _ in range(BACKFIT_SWEEPS):
        for j, entry in enumerate(accepted):
            lo, hi = entry["window"]
            window_r = r[lo:hi]
            others = [other["peak"] for i, other in enumerate(accepted) if i != j]
            target = d[lo:hi] - _truncated_models(window_r, q, weights, dq, rule, others)
            refit = fit_peak(window_r, target, q, weights, entry["peak"]["r"], entry["peak"]["sigma"],
                             (float(r[lo]), float(r[hi - 1])))
            if _peak_like(refit, entry["sign"], entry["halfWidth"]):
                entry["peak"] = refit
    kept = [{**entry["peak"], "rippleScale": ripple_scale(entry["peak"], qmax)} for entry in accepted]
    if not kept:
        return []
    top = max(peak["rippleScale"] for peak in kept)
    kept = [peak for peak in kept if peak["rippleScale"] >= STRONG_RIPPLE_FRACTION * top]
    kept.sort(key=lambda peak: (-peak["rippleScale"], peak["r"]))
    kept = kept[:MAX_STRONG_PEAKS]
    kept.sort(key=lambda peak: peak["r"])
    return kept


def pattern_on(r: np.ndarray, q: np.ndarray, weights: np.ndarray, dq: float, rule: str, peaks: list[dict]) -> np.ndarray:
    """Σ_j Δ_j(r): the strong peaks' predicted ripples at every r (a ripple decays only as 1/distance)."""
    r = np.asarray(r, dtype=float)
    total = np.zeros(r.size)
    for peak in peaks:
        total += ripple_pattern(q, weights, dq, rule, peak, r)
    return total


def peaks_on(r: np.ndarray, peaks: list[dict]) -> np.ndarray:
    """Σ_j d_j(r): the strong peaks with no truncation."""
    total = np.zeros(np.asarray(r).size)
    for peak in peaks:
        total += gaussian_d(r, peak)
    return total


def ripple_core(peak: dict, main_lobe: float) -> float:
    """Where a peak's ripples begin: beyond the kernel's main lobe and 2σ (within, Δ is the peak's broadening)."""
    return max(float(main_lobe), 2.0 * peak["sigma"])


def listed_ripples(r: np.ndarray, q: np.ndarray, weights: np.ndarray, dq: float, rule: str, peak: dict,
                   main_lobe: float, to_y: Callable[[np.ndarray, np.ndarray], np.ndarray]) -> list[dict]:
    """The first RIPPLES_LISTED extrema of one peak's Δ on each side, beyond its core (:func:`ripple_core`).

    ``{r, amplitude, relative}``: amplitude in the dataset's own units,
    relative to the peak's height there.
    """
    reach = RIPPLE_REACH_PERIODS * 2.0 * math.pi / float(q[-1])
    near = np.abs(r - peak["r"]) <= reach
    rr = r[near]
    if rr.size < 3:
        return []
    delta = to_y(rr, ripple_pattern(q, weights, dq, rule, peak, rr))
    height = float(to_y(np.array([peak["r"]]), np.array([peak["amplitude"]]))[0])
    core = ripple_core(peak, main_lobe)
    found = [i for i in local_extrema(delta) if abs(rr[i] - peak["r"]) > core]
    left = sorted((i for i in found if rr[i] < peak["r"]), key=lambda i: -rr[i])[:RIPPLES_LISTED]
    right = sorted((i for i in found if rr[i] > peak["r"]), key=lambda i: rr[i])[:RIPPLES_LISTED]
    out = []
    for i in sorted(left + right, key=lambda i: rr[i]):
        out.append({
            "r": float(rr[i]),
            "amplitude": float(delta[i]),
            "relative": float(delta[i] / height) if height != 0 else None,
        })
    return out


# ---------------------------------------------------------------------------
# The low-r level
# ---------------------------------------------------------------------------


def _baseline(r: np.ndarray, y: np.ndarray, form: str) -> float:
    """The baseline's one parameter by least squares: the constant (G-like) or the slope of a line through 0."""
    if form == "constant":
        return float(np.mean(y))
    denominator = float(np.dot(r, r))
    return float(np.dot(r, y)) / denominator if denominator > 0 else 0.0


def _baseline_values(r: np.ndarray, value: float, form: str) -> np.ndarray:
    return np.full(r.size, value) if form == "constant" else value * r


def enforced_below(r: np.ndarray, y: np.ndarray, form: str, scale: float) -> Optional[float]:
    """The largest r up to which y is exactly its baseline from the first point on (stog's RMC files), or None."""
    if r.size < 2 or not scale > 0:
        return None
    reference = y[0] if form == "constant" else (y[0] / r[0] if r[0] != 0 else None)
    if reference is None:
        return None
    last = None
    for i in range(r.size):
        expected = reference if form == "constant" else reference * r[i]
        if abs(y[i] - expected) <= ENFORCED_RTOL * scale:
            last = i
        else:
            break
    return float(r[last]) if last is not None and last >= 1 else None


def low_r_level(r: np.ndarray, y: np.ndarray, predicted: np.ndarray, *, form: str, r_first: float, fwhm: float,
                period: float, scale: float, model: Optional[np.ndarray] = None) -> dict:
    """The rms of y about its baseline below the first interatomic distance, against the predicted ripples.

    The window runs from the first r > 0 to ``r_first`` less the kernel's
    FWHM; ``predicted`` (Δ in y's units) gives the ripples the strong peaks
    put there. ``model``, where the model column is truncated as the data are
    (an X-ray F(Q) in real space), gives the model's own level. ``state``:
    ``measured``; ``enforced`` (y is its baseline exactly from the first point
    to ``enforcedBelow``, so the level cannot be read there); ``filtered``
    (the rms below FILTERED_RATIO of the predicted ripples: a low-r Fourier
    filter removed them); ``too-short`` (under LOW_R_MIN_PERIODS ripple
    periods); ``unavailable``.
    """
    out = {
        "state": "unavailable", "range": None, "baseline": None, "rms": None, "peak": None, "predictedRms": None,
        "residualRms": None, "modelRms": None, "dataMinusModelRms": None, "enforcedBelow": None,
    }
    r = np.asarray(r, dtype=float)
    hi = float(r_first) - float(fwhm)
    keep = (r > 0) & (r <= hi)
    if not np.any(keep):
        out["state"] = "too-short"
        return out
    rr, yy, pp = r[keep], y[keep], predicted[keep]
    enforced = enforced_below(rr, yy, form, scale)
    out["enforcedBelow"] = enforced
    if enforced is not None:
        rest = rr > enforced
        rr, yy, pp = rr[rest], yy[rest], pp[rest]
        if model is not None:
            model = np.asarray(model, dtype=float)[keep][rest]
        if rr.size < 2 or (rr[-1] - rr[0]) < LOW_R_MIN_PERIODS * period:
            out["state"] = "enforced"
            out["range"] = [float(r[keep][0]), float(r[keep][-1])]
            return out
    elif model is not None:
        model = np.asarray(model, dtype=float)[keep]
    if rr.size < 2 or (rr[-1] - rr[0]) < LOW_R_MIN_PERIODS * period:
        out["state"] = "too-short"
        out["range"] = [float(rr[0]), float(rr[-1])] if rr.size else None
        return out
    value = _baseline(rr, yy, form)
    about = yy - _baseline_values(rr, value, form)
    rms = float(math.sqrt(float(np.mean(about * about))))
    predicted_rms = float(math.sqrt(float(np.mean(pp * pp))))
    residual = about - pp
    out.update(
        range=[float(rr[0]), float(rr[-1])],
        baseline={"form": form, "value": value},
        rms=rms,
        peak=float(np.max(np.abs(about))),
        predictedRms=predicted_rms,
        residualRms=float(math.sqrt(float(np.mean(residual * residual)))),
        state="filtered" if predicted_rms > 0 and rms < FILTERED_RATIO * predicted_rms else "measured",
    )
    if model is not None and model.size == rr.size:
        model_value = _baseline(rr, model, form)
        model_about = model - _baseline_values(rr, model_value, form)
        out["modelRms"] = float(math.sqrt(float(np.mean(model_about * model_about))))
        out["baseline"]["model"] = model_value
        difference = yy - model
        out["dataMinusModelRms"] = float(math.sqrt(float(np.mean((difference - np.mean(difference)) ** 2))))
    return out


# ---------------------------------------------------------------------------
# One dataset: arrays → the assessment
# ---------------------------------------------------------------------------


def _d_space(r: np.ndarray, y: np.ndarray, space: str) -> np.ndarray:
    return np.asarray(y, dtype=float) * (np.asarray(r, dtype=float) if space == "G" else 1.0)


def _y_space(space: str) -> Callable[[np.ndarray, np.ndarray], np.ndarray]:
    def convert(r, d):
        r = np.asarray(r, dtype=float)
        d = np.asarray(d, dtype=float)
        if space != "G":
            return d
        out = np.zeros_like(d)
        nonzero = r != 0
        out[nonzero] = d[nonzero] / r[nonzero]
        return out
    return convert


def analyse_series(r, y, *, space: str, q_min: float, q_max: float, dq: float, rule: str, window: str,
                   r_first: Optional[float] = None, model=None, signs: str = "positive") -> dict:
    """The transform's numbers, the strong peaks and their ripples, the low-r level, for one series.

    ``r``, ``y``: the fitted function on a uniform grid (extended to low r where
    the run's files allow). ``space``: ``G`` (a G(r): D-space is r·y) or ``D``
    (a D(r)-like function, RMCProfile's X-ray transform included). The
    transform: Q_min + kΔQ up to Q_max, ``rule`` rectangle or trapezoid,
    ``window`` none or lorch. ``r_first``: the first interatomic distance (the
    configuration's closest approach), below which no peak is sought and the
    low-r level is read. ``model``: the model column on the same grid when it
    is truncated as the data are.
    """
    r = np.asarray(r, dtype=float)
    y = np.asarray(y, dtype=float)
    n = int(round((float(q_max) - float(q_min)) / float(dq))) + 1
    q = q_grid(q_min, dq, n)
    weights = transform_weights(q, dq, rule, window)
    res = resolution(q, weights)
    d = _d_space(r, y, space)
    to_y = _y_space(space)
    floor = None if r_first is None else float(r_first) - res["fwhm"]
    peaks = strong_peaks(r, d, q, weights, dq, rule, r_floor=floor, signs=signs)
    predicted_d = pattern_on(r, q, weights, dq, rule, peaks)
    predicted = to_y(r, predicted_d)
    listed = []
    for peak in peaks:
        height = float(to_y(np.array([peak["r"]]), np.array([peak["amplitude"]]))[0])
        listed.append({
            "r": peak["r"],
            "sigma": peak["sigma"],
            "amplitude": peak["amplitude"],
            "height": height,
            "rippleScale": peak["rippleScale"],
            "fitRms": float(peak["rms"]),
            "ripples": listed_ripples(r, q, weights, dq, rule, peak, res["mainLobe"], to_y),
        })
    scale = float(np.max(np.abs(y))) if y.size else 0.0
    low = None
    if r_first is not None:
        low = low_r_level(r, y, predicted, form="constant" if space == "G" else "line", r_first=float(r_first),
                          fwhm=res["fwhm"], period=res["ripplePeriod"], scale=scale, model=model)
    return {
        "space": space,
        "transform": {"qMin": float(q[0]), "qMax": float(q[-1]), "dQ": float(dq), "nQ": int(q.size), "rule": rule,
                      "window": window},
        "resolution": res,
        "peaks": listed,
        "lowR": low,
        "_arrays": {"r": r, "y": y, "d": d, "q": q, "weights": weights, "predicted": predicted,
                    "predictedD": predicted_d, "peaks": peaks},
    }


# ---------------------------------------------------------------------------
# Features
# ---------------------------------------------------------------------------


def _box_average(y: np.ndarray, half: int) -> np.ndarray:
    """The mean of y over i − half … i + half (fewer at the ends)."""
    n = y.size
    cumulative = np.concatenate([[0.0], np.cumsum(y)])
    out = np.empty(n)
    for i in range(n):
        lo, hi = max(0, i - half), min(n, i + half + 1)
        out[i] = (cumulative[hi] - cumulative[lo]) / (hi - lo)
    return out


def deviation(analysis: dict) -> dict:
    """A dataset's deviation from its strong peaks with no truncation, and the ripples predicted for it.

    In D-space h = d − P − B, P = Σ d_j and B the box average over
    BACKGROUND_PERIODS ripple periods of d − P − Δ (what the strong peaks and
    their ripples leave: weak peaks, the low-Q halo, the baseline), so h holds
    the ripples Δ and any feature; both returned in y's units.
    """
    arrays = analysis["_arrays"]
    r, d = arrays["r"], arrays["d"]
    step = float(r[1] - r[0])
    period = analysis["resolution"]["ripplePeriod"]
    peaks_d = peaks_on(r, arrays["peaks"])
    rest = d - peaks_d - arrays["predictedD"]
    half = max(1, int(round(0.5 * BACKGROUND_PERIODS * period / step)))
    background = _box_average(rest, half)
    h = d - peaks_d - background
    to_y = _y_space(analysis["space"])
    return {"r": r, "h": to_y(r, h), "predicted": arrays["predicted"]}


def _nearest_extremum(r: np.ndarray, values: np.ndarray, at: float, tolerance: float, sign: Optional[float] = None):
    """The local extremum of ``values`` nearest to ``at`` within ``tolerance`` (of the given sign), or None."""
    best = None
    for i in local_extrema(values):
        if abs(r[i] - at) > tolerance:
            continue
        if sign is not None and values[i] * sign <= 0:
            continue
        if best is None or abs(r[i] - at) < abs(r[best] - at):
            best = i
    return best


def _value_at(r: np.ndarray, values: np.ndarray, at: float) -> float:
    return float(np.interp(at, r, values))


def find_features(analysis: dict, r_min: float, r_max: float) -> list[float]:
    """Features in [r_min, r_max]: extrema of the deviation h beyond the strong peaks' cores (:func:`ripple_core`), by size.

    At least FEATURE_MIN_FRACTION of the largest strong peak's height; at most
    MAX_FEATURES; sorted by r.
    """
    dev = deviation(analysis)
    r, h = dev["r"], dev["h"]
    lobe = analysis["resolution"]["mainLobe"]
    peaks = analysis["peaks"]
    height = max((abs(peak["height"]) for peak in peaks), default=0.0)
    found = []
    for i in local_extrema(h):
        if not (r_min <= r[i] <= r_max):
            continue
        if any(abs(r[i] - peak["r"]) <= ripple_core(peak, lobe) for peak in peaks):
            continue
        if height > 0 and abs(h[i]) < FEATURE_MIN_FRACTION * height:
            continue
        found.append(i)
    found.sort(key=lambda i: (-abs(h[i]), r[i]))
    return sorted(float(r[i]) for i in found[:MAX_FEATURES])


def classify_feature(analysis: dict, at: float, noise: Optional[float] = None) -> dict:
    """One feature at ``at`` (Å) in one dataset: its size, the predicted ripple there, the reading.

    ``size``: the deviation h at the nearest extremum within the position
    tolerance (at ``at`` itself if none). ``ripple``: the nearest predicted
    extremum of the same sign within the tolerance. ``class``:
    ``ripple-consistent`` (|size| within a factor RIPPLE_SIZE_FACTOR of
    |ripple|, either way), ``larger-than-ripple`` (a ripple there, the feature
    larger: more than termination), ``smaller-than-ripple`` (the predicted
    ripple is not there in full: the data's transform differs from the one
    assumed — a window — or something of the other sign fills it),
    ``no-ripple-here`` (none of its sign within the tolerance, or one below
    NEGLIGIBLE_RIPPLE_FRACTION of it), ``within-noise`` (|size| < 2× the
    dataset's noise level, where known), or ``resolved-peak`` (a strong peak
    lies within the tolerance: fitted as a peak in its own right, which the
    larger peaks' ripples did not explain). ``fromPeak``: the nearest strong
    peak and the distance; ``belowResolution``: that distance under the
    kernel's FWHM.
    """
    dev = deviation(analysis)
    r, h, predicted = dev["r"], dev["h"], dev["predicted"]
    res = analysis["resolution"]
    tolerance = RIPPLE_POSITION_FRACTION * res["ripplePeriod"]
    nearest = None
    for peak in analysis["peaks"]:
        distance = abs(float(at) - peak["r"])
        if nearest is None or distance < nearest["distance"]:
            nearest = {"r": peak["r"], "distance": distance, "height": peak["height"], "sigma": peak["sigma"]}
    if nearest is not None and nearest["distance"] <= tolerance:
        # A strong peak itself: the larger peaks' ripples did not explain it (strong_peaks).
        return {
            "r": nearest["r"], "asked": float(at), "size": nearest["height"], "predictedHere": _value_at(r, predicted, nearest["r"]),
            "ripple": None, "class": "resolved-peak", "fromPeak": None, "belowResolution": False,
        }
    index = _nearest_extremum(r, h, at, tolerance)
    position = float(r[index]) if index is not None else float(at)
    size = float(h[index]) if index is not None else _value_at(r, h, at)
    sign = 1.0 if size >= 0 else -1.0
    ripple_index = _nearest_extremum(r, predicted, position, tolerance, sign)
    ripple = None if ripple_index is None else {"r": float(r[ripple_index]), "amplitude": float(predicted[ripple_index])}
    nearest = None
    for peak in analysis["peaks"]:
        distance = abs(position - peak["r"])
        if nearest is None or distance < nearest["distance"]:
            nearest = {"r": peak["r"], "distance": distance}
    if noise is not None and abs(size) < 2.0 * noise:
        kind = "within-noise"
    elif ripple is None or abs(ripple["amplitude"]) < NEGLIGIBLE_RIPPLE_FRACTION * abs(size):
        kind = "no-ripple-here"
    elif abs(size) > RIPPLE_SIZE_FACTOR * abs(ripple["amplitude"]):
        kind = "larger-than-ripple"
    elif abs(size) < abs(ripple["amplitude"]) / RIPPLE_SIZE_FACTOR:
        kind = "smaller-than-ripple"
    else:
        kind = "ripple-consistent"
    return {
        "r": position,
        "asked": float(at),
        "size": size,
        "predictedHere": _value_at(r, predicted, position),
        "ripple": ripple,
        "class": kind,
        "fromPeak": nearest,
        "belowResolution": bool(nearest is not None and nearest["distance"] < res["fwhm"]),
    }


def split_reading(analyses: list[dict], r1: float, r2: float) -> dict:
    """A split or shoulder of spacing |r₂ − r₁|: below each dataset's resolution (kernel FWHM) or not."""
    spacing = abs(float(r2) - float(r1))
    return {
        "spacing": spacing,
        "datasets": [{"fwhm": a["resolution"]["fwhm"], "belowResolution": bool(spacing < a["resolution"]["fwhm"])}
                     for a in analyses],
    }


def across_datasets(readings: list[dict], analyses: list[dict]) -> dict:
    """One feature read in several datasets: a ripple moves with Q_max, a real feature stays.

    ``present`` per dataset: an extremum of h within the tolerance, of the
    sign the feature has where its relative size (size over the dataset's
    largest strong-peak height) is largest, and at least PRESENT_FRACTION of
    that size (or a resolved peak there). ``verdict``:
    ``moves-with-qmax`` (ripple-consistent in one, absent in another, and the
    Q_max differ enough to tell: ``discriminating``, a ripple at this distance
    from its peak moving by more than the tolerance; ``not-discriminating``
    when they do not), ``in-every-dataset`` (present everywhere, a ripple in
    none), ``ripple-in-every-dataset``, or ``mixed``. Compare datasets that
    measure the same total function only: X-ray and neutron weight the pairs
    differently.
    """
    relative = []
    for reading, analysis in zip(readings, analyses):
        height = max((abs(peak["height"]) for peak in analysis["peaks"]), default=0.0)
        relative.append(abs(reading["size"]) / height if height > 0 else 0.0)
    largest = max(relative, default=0.0)
    reference = readings[relative.index(largest)] if readings else None
    sign = 1.0 if reference is None or reference["size"] >= 0 else -1.0
    present = [reading["class"] == "resolved-peak"
               or (bool(largest > 0 and rel >= PRESENT_FRACTION * largest) and reading["size"] * sign > 0
                   and reading["class"] != "within-noise")
               for rel, reading in zip(relative, readings)]
    ripple = [reading["class"] == "ripple-consistent" for reading in readings]
    qmaxes = [a["transform"]["qMax"] for a in analyses]
    distance = max((reading["fromPeak"]["distance"] for reading in readings if reading["fromPeak"]), default=0.0)
    tolerance = min(RIPPLE_POSITION_FRACTION * a["resolution"]["ripplePeriod"] for a in analyses)
    shift = distance * abs(1.0 - min(qmaxes) / max(qmaxes)) if qmaxes and max(qmaxes) > 0 else 0.0
    discriminating = bool(shift > tolerance)
    if any(ripple) and not all(present):
        verdict = "moves-with-qmax" if discriminating else "not-discriminating"
    elif all(present) and not any(ripple):
        verdict = "in-every-dataset"
    elif all(ripple):
        verdict = "ripple-in-every-dataset"
    else:
        verdict = "mixed"
    return {"present": present, "ripple": ripple, "verdict": verdict, "discriminating": discriminating}


def classify_features(analyses: list[dict], features: list[float], noise: Optional[list] = None,
                      totals: Optional[list] = None) -> list[dict]:
    """Each feature in every dataset, and across those that measure a total function when two or more do.

    ``totals``: per analysis, whether it is a total G(r) (a partial one is not
    compared across); all of them when None. ``across`` holds the indices it
    compared (``compared``).
    """
    compared = [k for k in range(len(analyses)) if totals is None or totals[k]]
    out = []
    for at in features:
        readings = [classify_feature(a, at, None if noise is None else noise[k]) for k, a in enumerate(analyses)]
        entry = {"asked": float(at), "datasets": readings}
        if len(compared) > 1:
            entry["across"] = {**across_datasets([readings[k] for k in compared], [analyses[k] for k in compared]),
                               "compared": compared}
        out.append(entry)
    return out


# ---------------------------------------------------------------------------
# The fit space, from the run-control file
# ---------------------------------------------------------------------------

_DAT_KEY_RE = re.compile(r"^([A-Z_\d]+)\s*::\s*(.*)$")
_DAT_SUB_RE = re.compile(r"^([A-Z_]+)\s*(?:::\s*(.*))?$")
_LINE_BREAK_RE = re.compile(r"\r\n|\r|\n")
_NUMBER_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")


def _numbers(value: str) -> list[float]:
    """The plain decimal numbers a value starts with."""
    out = []
    for token in value.split():
        if not _NUMBER_RE.match(token):
            break
        out.append(float(token))
    return out


def _first_number(value: str) -> Optional[float]:
    for token in value.split():
        return float(token) if _NUMBER_RE.match(token) else None
    return None


def ripple_settings(text: str) -> dict:
    """What the ``.dat`` says about each dataset's transform: its block, file, fit type and the keys that matter.

    ``{minimumDistances, datasets}``; each dataset ``{block, file, fitType,
    dataType, realSpaceFit, convolve, convolveQ, rCutoff, qMin, qMax, qStep,
    coefficients}`` in file order (``convolveQ``: the number after ``CONVOLVE ::``, a
    real-space block's stated Q range).
    """
    datasets: list[dict] = []
    block: Optional[dict] = None
    minimum = None
    for raw in _LINE_BREAK_RE.split(text):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if block is None:
                continue
            match = _DAT_SUB_RE.match(line[1:].strip())
            if not match:
                continue
            key, value = match.group(1), (match.group(2) or "").strip()
            if key == "FILENAME":
                block["file"] = value
            elif key == "FIT_TYPE":
                block["fitType"] = value
            elif key == "DATA_TYPE":
                block["dataType"] = value
            elif key == "REAL_SPACE_FIT":
                block["realSpaceFit"] = True
            elif key == "CONVOLVE":
                block["convolve"] = True
                block["convolveQ"] = _first_number(value)
            elif key == "R_CUTOFF":
                block["rCutoff"] = _first_number(value)
            elif key == "QMIN":
                block["qMin"] = _first_number(value)
            elif key == "QMAX":
                block["qMax"] = _first_number(value)
            elif key == "QSTEP":
                block["qStep"] = _first_number(value)
            elif key == "NEUTRON_COEFFICIENTS":
                block["coefficients"] = _numbers(value)
            continue
        match = _DAT_KEY_RE.match(line)
        block = None
        if not match:
            continue
        key, value = match.groups()
        if key == "END":
            break
        if key.endswith("_DATA"):
            block = {"block": key, "file": None, "fitType": None, "dataType": None, "realSpaceFit": False,
                     "convolve": False, "convolveQ": None, "rCutoff": None, "qMin": None, "qMax": None, "qStep": None,
                     "coefficients": None}
            datasets.append(block)
        elif key == "MINIMUM_DISTANCES":
            minimum = _numbers(value)
    return {"minimumDistances": minimum, "datasets": datasets}


def _nth_block(settings: Optional[dict], block: str, number: int) -> Optional[dict]:
    if not settings:
        return None
    blocks = [entry for entry in settings["datasets"] if entry["block"] == block]
    return blocks[number - 1] if 0 < number <= len(blocks) else None


def is_d_like(fit_type: Optional[str]) -> bool:
    """D(r), T(r) and PDFgui's G(r)P carry the factor r: their D-space is the function itself."""
    text = (fit_type or "").replace(" ", "").upper()
    return text.startswith("D(R)") or text.startswith("T(R)") or text.startswith("G(R)P")


def fit_space_risk(kind: str, dataset: Optional[dict], data_transform: Optional[dict]) -> dict:
    """How RMCProfile makes the model this dataset is compared with, and the risk that the box fits ripples.

    ``kind``: ``xpdf`` (the transform of an X-ray F(Q)), ``npdf`` (a neutron
    G(r) PDF<n> file), ``xray-real`` (an X-ray G(r) Xray_PDF<n> file),
    ``reciprocal`` (a fit in Q). ``data_transform``: the data's own transform
    as the files establish it ({qMax, window, source}) or None. Returns
    ``{fitSpace, model, modelQ, rCutoff, risk, because}`` — ``because`` a code
    the tool renders.
    """
    if kind == "reciprocal":
        return {"fitSpace": "reciprocal", "model": None, "modelQ": None, "rCutoff": None, "risk": "none",
                "because": "fitted-in-q"}
    if kind == "xpdf":
        return {"fitSpace": "real", "model": "same-transform", "modelQ": None, "rCutoff": None, "risk": "matched",
                "because": "same-transform"}
    known = data_transform is not None and data_transform.get("source") in ("provenance", "stog")
    if kind == "npdf":
        if not dataset or not dataset.get("convolve"):
            return {"fitSpace": "real", "model": "histogram", "modelQ": None, "rCutoff": None, "risk": "high",
                    "because": "untruncated-model"}
        model_q = dataset.get("convolveQ")
        r_cut = dataset.get("rCutoff")
        out = {"fitSpace": "real", "model": "sinc", "modelQ": model_q, "rCutoff": r_cut}
        if model_q is None:
            return {**out, "risk": "unknown", "because": "sinc-no-q"}
        if not known:
            return {**out, "risk": "unknown", "because": "sinc-data-unverified"}
        if data_transform.get("window") == "lorch":
            return {**out, "risk": "high", "because": "sinc-vs-lorch"}
        if abs(model_q - data_transform["qMax"]) > 0.02 * data_transform["qMax"]:
            return {**out, "risk": "high", "because": "sinc-q-differs"}
        return {**out, "risk": "matched", "because": "sinc-matches"}
    if kind == "xray-real":
        model_q = dataset.get("qMax") if dataset else None
        out = {"fitSpace": "real", "model": "xray-transform", "modelQ": model_q, "rCutoff": None}
        if model_q is None or not known:
            return {**out, "risk": "unknown", "because": "xray-real-unverified"}
        if data_transform.get("window") == "lorch":
            return {**out, "risk": "high", "because": "transform-vs-lorch"}
        if abs(model_q - data_transform["qMax"]) > 0.02 * data_transform["qMax"]:
            return {**out, "risk": "high", "because": "transform-q-differs"}
        return {**out, "risk": "matched", "because": "transform-matches"}
    return {"fitSpace": None, "model": None, "modelQ": None, "rCutoff": None, "risk": "unknown", "because": "unknown-kind"}


# ---------------------------------------------------------------------------
# What the data files say about their own transform
# ---------------------------------------------------------------------------


def provenance_transform(payload: dict, data_name: str) -> Optional[dict]:
    """An Auto StoG provenance JSON's transform, if one of its outputs is ``data_name`` (a basename)."""
    if not isinstance(payload, dict) or payload.get("tool") != "rmc-autoscale":
        return None
    outputs = payload.get("outputs") or {}
    names = {str(path).replace("\\", "/").split("/")[-1].lower() for path in outputs.values()}
    if data_name.lower() not in names:
        return None
    config = (payload.get("provenance") or {}).get("config") or {}
    try:
        q_min, q_max = float(config["qmin"]), float(config["qmax"])
    except (KeyError, TypeError, ValueError):
        return None
    n_points = (payload.get("provenance") or {}).get("n_q_points")
    dq = (q_max - q_min) / (float(n_points) - 1.0) if isinstance(n_points, (int, float)) and n_points > 1 else None
    enforcement = payload.get("enforcement") or None
    return {
        "source": "provenance",
        "qMin": q_min,
        "qMax": q_max,
        "dQ": dq,
        "window": "lorch" if config.get("lorch") else "none",
        "filterCutoff": float(config["r_cutoff"]) if config.get("r_cutoff") else None,
        "enforceCutoff": float(enforcement["cutoff"]) if isinstance(enforcement, dict) and enforcement.get("cutoff") is not None else None,
    }


def stog_transform(inp, data_name: str) -> Optional[dict]:
    """A classic stog input's transform (parsers.StogInput), if one of its outputs is ``data_name``."""
    names = {str(name).replace("\\", "/").split("/")[-1].lower()
             for name in (inp.out_rmc_gr, inp.out_rmc_dr, inp.out_rmc_fq, inp.out_ft_gr, inp.out_gr)}
    if data_name.lower() not in names:
        return None
    return {
        "source": "stog",
        "qMin": float(inp.qmin),
        "qMax": float(inp.qmax),
        "dQ": None,
        "window": "lorch" if inp.lorch else "none",
        "filterCutoff": float(inp.r_cutoff) if inp.use_filter and inp.r_cutoff > 0 else None,
        "enforceCutoff": float(inp.peak_cutoff) if inp.peak_cutoff > 0 else None,
    }


def auto_stog_title(text: str) -> bool:
    """Whether a data file's title line names Auto StoG (``rmc-autoscale <version>: …``)."""
    lines = [line for line in _LINE_BREAK_RE.split(text)[:3]]
    return any(line.strip().startswith("rmc-autoscale") for line in lines)


# ---------------------------------------------------------------------------
# One dataset of a run: its files → its assessment
# ---------------------------------------------------------------------------

_XRAY_PDF_RE = re.compile(r"^(.+)_Xray_PDF(\d+)\.csv$")
_PDF_RE = re.compile(r"^(.+)_PDF(\d+)\.csv$")
_FQ_RE = re.compile(r"^(.+)_FQ(\d+)\.csv$")
_SQ_RE = re.compile(r"^(.+)_SQ(\d+)\.csv$")

#: The .dat block each kind of file belongs to.
BLOCKS = {
    "xpdf": "XRAY_RECIPROCAL_SPACE_DATA",
    "npdf": "NEUTRON_REAL_SPACE_DATA",
    "xray-real": "XRAY_REAL_SPACE_DATA",
    "xray-reciprocal": "XRAY_RECIPROCAL_SPACE_DATA",
    "neutron-reciprocal": "NEUTRON_RECIPROCAL_SPACE_DATA",
}


def dataset_file(name: str) -> Optional[dict]:
    """What a run output file is: ``{kind, stem, number}`` for a fit this module reads, else None.

    ``<stem>_FT_XFQ<n>.csv`` the transform of an X-ray F(Q) (xpdf),
    ``<stem>_Xray_PDF<n>.csv`` an X-ray G(r), ``<stem>_PDF<n>.csv`` a neutron
    G(r), ``<stem>_FQ<n>.csv`` an X-ray F(Q), ``<stem>_SQ<n>.csv`` a neutron
    S(Q) or F(Q).
    """
    for pattern, kind in ((FT_FILE_RE, "xpdf"), (_XRAY_PDF_RE, "xray-real"), (_PDF_RE, "npdf"),
                          (_FQ_RE, "xray-reciprocal"), (_SQ_RE, "neutron-reciprocal")):
        match = pattern.match(name)
        if match:
            return {"kind": kind, "stem": match.group(1), "number": int(match.group(2))}
    return None


def closest_approach(columns: Optional[np.ndarray], minimum_distances: Optional[list]) -> Optional[dict]:
    """The first interatomic distance: the configuration's (the first r where any partial g(r) > 0), else the
    smallest positive MINIMUM_DISTANCES."""
    if columns is not None and columns.shape[0] >= 2:
        r = columns[0]
        partials = columns[1:]
        occupied = np.any(np.nan_to_num(partials) > 0, axis=0)
        if np.any(occupied):
            return {"value": float(r[int(np.argmax(occupied))]), "source": "partials"}
    positive = [value for value in (minimum_distances or []) if value > 0]
    if positive:
        return {"value": float(min(positive)), "source": "minimum-distances"}
    return None


def _uniform(r: np.ndarray) -> bool:
    if r.size < 3:
        return False
    step = float((r[-1] - r[0]) / (r.size - 1))
    return step > 0 and float(np.max(np.abs(r - (r[0] + step * np.arange(r.size))))) <= 1e-6 * step


def map_input(r_file: np.ndarray, y_file: np.ndarray, columns: np.ndarray) -> Optional[dict]:
    """The input data file in the fitted file's units: y_file ≈ s·y_input + c on the common r points.

    Returns ``{r, y, scale, offset, residual}`` over the input's r up to the
    fitted file's last r, or None when fewer than 10 points are common or the
    map misses by more than INPUT_MAP_TOLERANCE of the largest |y_file|.
    """
    r_in, y_in = np.asarray(columns[0], dtype=float), np.asarray(columns[1], dtype=float)
    finite = np.isfinite(r_in) & np.isfinite(y_in)
    r_in, y_in = r_in[finite], y_in[finite]
    if r_in.size < 10:
        return None
    step = float(r_file[1] - r_file[0]) if r_file.size > 1 else 0.0
    index = np.searchsorted(r_in, r_file)
    matched_file, matched_input = [], []
    for k, (r_value, i) in enumerate(zip(r_file, index)):
        for j in (i - 1, i):
            if 0 <= j < r_in.size and abs(r_in[j] - r_value) <= 1e-6 * max(step, 1e-12):
                matched_file.append(k)
                matched_input.append(j)
                break
    if len(matched_file) < 10:
        return None
    x = y_in[matched_input]
    target = y_file[matched_file]
    n = float(x.size)
    sx, sy, sxx, sxy = float(np.sum(x)), float(np.sum(target)), float(np.dot(x, x)), float(np.dot(x, target))
    det = n * sxx - sx * sx
    if det == 0:
        return None
    scale = (n * sxy - sx * sy) / det
    offset = (sxx * sy - sx * sxy) / det
    largest = float(np.max(np.abs(target)))
    residual = float(np.max(np.abs(target - scale * x - offset))) / largest if largest > 0 else math.inf
    if not residual <= INPUT_MAP_TOLERANCE or not scale > 0:
        return None
    keep = (r_in > 0) & (r_in <= float(r_file[-1]) + 1e-9)
    return {"r": r_in[keep], "y": scale * y_in[keep] + offset, "scale": float(scale), "offset": float(offset),
            "residual": residual}


def is_total(dataset: Optional[dict]) -> bool:
    """A dataset measures a total function unless its NEUTRON_COEFFICIENTS keep one pair only (a partial)."""
    coefficients = (dataset or {}).get("coefficients") or []
    return not (coefficients and sum(1 for value in coefficients if value != 0) == 1)


def expected_baseline(dataset: Optional[dict], space: str) -> Optional[float]:
    """RMCProfile's G(r) below the closest approach: −Σ w_ij of the block's coefficients (its model there)."""
    coefficients = (dataset or {}).get("coefficients")
    if space != "G" or not coefficients:
        return None
    return -float(sum(coefficients))


def _block_summary(dataset: Optional[dict]) -> Optional[dict]:
    if not dataset:
        return None
    return {key: dataset.get(key) for key in ("block", "file", "fitType", "realSpaceFit", "convolve", "convolveQ",
                                              "rCutoff", "qMin", "qMax")}


def assess_dataset(*, name: str, kind: str, number: int, settings: Optional[dict], x, calc, data,
                   fq: Optional[tuple] = None, input_columns: Optional[np.ndarray] = None,
                   data_transform: Optional[dict] = None, override: Optional[dict] = None,
                   r_first: Optional[dict] = None, processing: Optional[dict] = None) -> dict:
    """One fitted dataset of a run: how its model is made, its transform, ripples and low-r level.

    ``x``, ``calc``, ``data``: the file's columns (r or Q, model, experiment).
    ``fq``: an xpdf's run F(Q) (Q, model, experiment). ``input_columns``: a
    real-space block's input data file (read_stog_xy). ``data_transform``:
    the data's own transform where an Auto StoG provenance or a stog input
    names the file (:func:`provenance_transform`, :func:`stog_transform`).
    ``override``: the caller's ``{qMax, qMin, window}`` for a dataset whose
    files do not say. ``r_first``: :func:`closest_approach`. Returns the
    payload with the analysis kept under ``_analysis`` (removed by
    :func:`public_payload`).
    """
    dataset = _nth_block(settings, BLOCKS.get(kind, ""), number)
    real = kind in ("xpdf", "npdf", "xray-real")
    risk_kind = kind if real else "reciprocal"
    payload = {
        "status": "unavailable", "reason": None, "file": name, "kind": kind, "number": int(number),
        "total": is_total(dataset),
        "block": _block_summary(dataset), "fit": fit_space_risk(risk_kind, dataset, data_transform),
        "transform": None, "resolution": None, "peaks": [], "lowR": None, "processing": processing,
        "space": None, "rFirst": r_first, "_analysis": None,
    }
    if kind == "xpdf" and dataset is not None and not dataset.get("realSpaceFit"):
        payload["fit"] = {**payload["fit"], "fitSpace": "reciprocal", "risk": "none", "because": "fitted-in-q"}
    elif kind == "xpdf" and dataset is None:
        # No .dat to say whether the transform is fitted: if it is, the model is
        # transformed alike; if not, the fit is in Q. Either way no ripple is fitted.
        payload["fit"] = {**payload["fit"], "fitSpace": "unknown", "because": "same-transform-no-dat"}
    x = np.asarray(x, dtype=float)
    data = np.asarray(data, dtype=float)
    if not real:
        try:
            dq = uniform_step(x)
        except ValueError as error:
            payload.update(status="failed", reason=f"the Q grid: {error}")
            return payload
        q = x
        weights = transform_weights(q, dq, "rectangle", "none")
        payload["transform"] = {"qMin": float(q[0]), "qMax": float(q[-1]), "dQ": dq, "nQ": int(q.size),
                                "rule": "rectangle", "window": "none", "source": "q-file", "windowKnown": False}
        payload["resolution"] = resolution(q, weights)
        payload["status"] = "assessed"
        return payload

    model = None
    if kind == "xpdf":
        if fq is None:
            payload["reason"] = "no-fq-file"
            return payload
        q, f_calc, f_exp = (np.asarray(column, dtype=float) for column in fq)
        try:
            dq = uniform_step(q)
        except ValueError as error:
            payload.update(status="failed", reason=f"the Q grid: {error}")
            return payload
        reproduction = reproduce_transform(q, f_exp, x, data)
        if not reproduction["pass"]:
            payload.update(status="failed", reason="transform-not-reproduced",
                           transform={"reproduction": {**reproduction, "tolerance": RMC_FT_TOLERANCE}})
            return payload
        step = float(x[1] - x[0])
        count = int(math.floor(float(x[0]) / step - 1e-9))
        r_low = float(x[0]) - step * np.arange(count, 0, -1, dtype=float)
        r_low = r_low[r_low > 0]
        coefficients = (dq * q)
        low = np.sin(np.outer(r_low, q)) @ (coefficients * (f_exp + reproduction["offset"])) if r_low.size else np.array([])
        low_model = np.sin(np.outer(r_low, q)) @ (coefficients * f_calc) if r_low.size else np.array([])
        rr = np.concatenate([r_low, x])
        yy = np.concatenate([low, data])
        model = np.concatenate([low_model, np.asarray(calc, dtype=float)])
        q_min, q_max, rule, window, source, window_known = float(q[0]), float(q[-1]), "rectangle", "none", "run-fq", True
        transform_extra = {"reproduction": {"residual": reproduction["residual"], "offset": reproduction["offset"],
                                            "tolerance": RMC_FT_TOLERANCE}}
        space = "D"
    else:
        space = "D" if is_d_like(dataset.get("fitType") if dataset else None) else "G"
        override = override or {}
        if override.get("qMax") is not None:
            source = "argument"
            q_max = float(override["qMax"])
            q_min = override.get("qMin")
            window = override.get("window") or "none"
            window_known = override.get("window") is not None
            dq = DEFAULT_DQ
        elif data_transform is not None:
            source = data_transform["source"]
            q_max, q_min = data_transform["qMax"], data_transform["qMin"]
            window, window_known = data_transform["window"], True
            dq = data_transform.get("dQ") or DEFAULT_DQ
        elif kind == "npdf" and dataset and dataset.get("convolveQ"):
            source, q_max, q_min, window, window_known, dq = "convolve", float(dataset["convolveQ"]), None, "none", False, DEFAULT_DQ
        elif kind == "xray-real" and dataset and dataset.get("qMax"):
            source, q_max, q_min, window, window_known = "xray-block", float(dataset["qMax"]), dataset.get("qMin"), "none", False
            dq = dataset.get("qStep") or DEFAULT_DQ
        else:
            payload["reason"] = "qmax-unknown"
            return payload
        q_min = dq if q_min is None or q_min < dq else float(q_min)
        rule = "rectangle" if source == "xray-block" else "trapezoid"
        rr, yy = x, data
        transform_extra = {}
        if input_columns is not None:
            mapped = map_input(x, data, input_columns)
            if mapped is not None:
                rr, yy = mapped["r"], mapped["y"]
                transform_extra = {"input": {"scale": mapped["scale"], "offset": mapped["offset"],
                                             "residual": mapped["residual"]}}
    payload["space"] = space
    finite = np.isfinite(rr) & np.isfinite(yy)
    rr, yy = rr[finite], yy[finite]
    if model is not None:
        model = model[finite]
    if not _uniform(rr):
        payload.update(status="failed", reason="r-grid-not-uniform")
        return payload
    signs = "both" if dataset and any(value < 0 for value in (dataset.get("coefficients") or [])) else "positive"
    analysis = analyse_series(rr, yy, space=space, q_min=q_min, q_max=q_max, dq=dq, rule=rule, window=window,
                              r_first=None if r_first is None else r_first["value"], model=model, signs=signs)
    expected = expected_baseline(dataset, space)
    if expected is not None and analysis["lowR"] and analysis["lowR"].get("baseline"):
        analysis["lowR"]["baseline"]["expected"] = expected
    payload.update(
        status="assessed",
        transform={**analysis["transform"], "source": source, "windowKnown": window_known, **transform_extra},
        resolution=analysis["resolution"],
        peaks=analysis["peaks"],
        lowR=analysis["lowR"],
        range=[float(rr[0]), float(rr[-1])],
        _analysis=analysis,
    )
    return payload


def noise_level(payload: dict) -> Optional[float]:
    """A dataset's noise from its low-r level: the data less the model where both are truncated alike, else
    the rms left after the predicted ripples; None unless the level was measured."""
    low = payload.get("lowR")
    if not low or low.get("state") != "measured":
        return None
    if low.get("dataMinusModelRms") is not None:
        return float(low["dataMinusModelRms"])
    return float(low["residualRms"]) if low.get("residualRms") is not None else None


def pair_reading(r: np.ndarray, g: np.ndarray, features: list[float], payloads: list[dict]) -> list[dict]:
    """Features of a partial g(r) (the box's own histogram) against the real-space datasets' predicted ripples.

    Per feature: the partial's nearest local maximum within the position
    tolerance, the distance to the partial's nearest other peak (its maxima
    above 15 % of its largest) against the smallest kernel FWHM of the run,
    and, for every real-space dataset whose box may fit ripples (risk high or
    unknown), the nearest predicted ripple extremum within the tolerance with
    its size relative to that dataset's largest strong peak.
    """
    r = np.asarray(r, dtype=float)
    g = np.asarray(g, dtype=float)
    analysed = [p for p in payloads if p["status"] == "assessed" and p.get("_analysis") is not None]
    resolved = [p for p in payloads if p.get("resolution")]
    fwhm = min((p["resolution"]["fwhm"] for p in resolved), default=None)
    period = min((p["resolution"]["ripplePeriod"] for p in resolved), default=None)
    tolerance = RIPPLE_POSITION_FRACTION * period if period else 0.05
    largest = float(np.nanmax(g)) if g.size else 0.0
    maxima = [i for i in local_extrema(g) if g[i] > g[i - 1] and g[i] >= 0.15 * largest]
    out = []
    for at in features:
        near = [i for i in maxima if abs(r[i] - at) <= tolerance]
        position = float(r[min(near, key=lambda i: abs(r[i] - at))]) if near else float(at)
        others = [float(r[i]) for i in maxima if abs(r[i] - position) > tolerance]
        neighbour = min(others, key=lambda value: abs(value - position)) if others else None
        spacing = abs(neighbour - position) if neighbour is not None else None
        ripples = []
        for p in analysed:
            if p["fit"]["risk"] not in ("high", "unknown"):
                continue
            arrays = p["_analysis"]["_arrays"]
            predicted = arrays["predicted"]
            index = _nearest_extremum(arrays["r"], predicted, position, tolerance)
            if index is None:
                continue
            height = max((abs(peak["height"]) for peak in p["peaks"]), default=0.0)
            ripples.append({"file": p["file"], "r": float(arrays["r"][index]),
                            "relative": float(predicted[index] / height) if height > 0 else None})
        out.append({
            "asked": float(at), "r": position, "inPartial": bool(near), "neighbour": neighbour, "spacing": spacing,
            "belowResolution": bool(spacing is not None and fwhm is not None and spacing < fwhm),
            "atRipple": ripples,
        })
    return out


def assess_run(payloads: list[dict], *, features: Optional[list] = None, r_range: Optional[list] = None,
               split: Optional[list] = None, pair: Optional[dict] = None) -> dict:
    """The run's datasets, and the caller's features read in every real-space one (and across them).

    ``features``: positions (Å); else ``r_range`` [r_min, r_max] finds them
    (:func:`find_features` in each dataset, merged within the position
    tolerance, the MAX_FEATURES largest relative to their dataset's peaks).
    ``split``: [r₁, r₂]. ``pair``: ``{label, r, g}`` — the features are read in
    that partial g(r) of the box instead (:func:`pair_reading`).
    """
    analysed = [p for p in payloads if p["status"] == "assessed" and p.get("_analysis") is not None]
    analyses = [p["_analysis"] for p in analysed]
    found = None
    if features is None and r_range is not None and analyses:
        candidates = []
        for analysis in analyses:
            height = max((abs(peak["height"]) for peak in analysis["peaks"]), default=0.0)
            dev = None
            for at in find_features(analysis, float(r_range[0]), float(r_range[1])):
                if dev is None:
                    dev = deviation(analysis)
                size = abs(_value_at(dev["r"], dev["h"], at))
                candidates.append((size / height if height > 0 else 0.0, at))
        tolerance = min(RIPPLE_POSITION_FRACTION * a["resolution"]["ripplePeriod"] for a in analyses)
        merged: list[tuple[float, float]] = []
        for relative, at in sorted(candidates, key=lambda entry: (-entry[0], entry[1])):
            if all(abs(at - other) > tolerance for _, other in merged):
                merged.append((relative, at))
        found = sorted(at for _, at in merged[:MAX_FEATURES])
        features = found
    out = {
        "datasets": [public_payload(p) for p in payloads],
        "features": None,
        "found": found is not None,
        "split": None,
        "pair": None,
        "rules": RULES,
    }
    if features:
        if pair is not None:
            out["pair"] = {"label": pair["label"], "features": pair_reading(pair["r"], pair["g"], features, payloads)}
        elif analysed:
            readings = classify_features(analyses, [float(at) for at in features], [noise_level(p) for p in analysed],
                                         [bool(p.get("total", True)) for p in analysed])
            for entry in readings:
                for reading, p in zip(entry["datasets"], analysed):
                    reading["file"] = p["file"]
            out["features"] = readings
    if split is not None and len(split) == 2:
        resolved = [p for p in payloads if p.get("resolution")]
        reading = split_reading([{"resolution": p["resolution"]} for p in resolved], split[0], split[1])
        for entry, p in zip(reading["datasets"], resolved):
            entry["file"] = p["file"]
        out["split"] = reading
    return out


def public_payload(payload: dict) -> dict:
    """The payload without its internal analysis."""
    return {key: value for key, value in payload.items() if key != "_analysis"}


# ---------------------------------------------------------------------------
# The run folder
# ---------------------------------------------------------------------------


def _read_columns(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A fit CSV's x, model and experiment columns (parsers.rwp_columns' roles)."""
    from .parsers import read_rmc_csv, rwp_columns

    series = read_rmc_csv(path)
    roles = rwp_columns(series.labels, len(series.data))
    if roles is None:
        raise ValueError(f"{path.name} needs an experimental column beside the calculated one")
    return series.data[0], series.data[roles[0]], series.data[roles[1]]


def data_processing(folder: Path, data_name: str, input_text: Optional[str]) -> tuple[Optional[dict], Optional[dict]]:
    """The data file's own transform and processing, from the files beside it.

    An Auto StoG provenance (``*_provenance.json`` in the folder or its
    ``autoscale/``) or a classic stog input (``*.inp``, ``stog_input.dat``)
    whose outputs name ``data_name``; the data file's title line naming
    rmc-autoscale. Returns (transform, processing) — either None.
    """
    from .parsers import read_stog_inp

    transform = None
    for directory in (folder, folder / "autoscale"):
        if transform is not None or not directory.is_dir():
            continue
        for path in sorted(directory.glob("*_provenance.json")):
            try:
                found = provenance_transform(json.loads(path.read_text(encoding="utf-8")), data_name)
            except (OSError, ValueError):
                continue
            if found is not None:
                transform = {**found, "file": path.name}
                break
        if transform is None:
            for path in sorted([*directory.glob("*.inp"), *directory.glob("stog_input.dat")]):
                try:
                    found = stog_transform(read_stog_inp(path), data_name)
                except (OSError, ValueError, NotImplementedError, IndexError):
                    continue
                if found is not None:
                    transform = {**found, "file": path.name}
                    break
    auto = bool(input_text) and auto_stog_title(input_text)
    processing = None
    if transform is not None or auto:
        processing = {
            "autoStog": bool(auto or (transform is not None and transform["source"] == "provenance")),
            "filterCutoff": transform.get("filterCutoff") if transform else None,
            "enforceCutoff": transform.get("enforceCutoff") if transform else None,
            "source": transform.get("file") if transform else None,
        }
    return transform, processing


def ripples_for_run(paths, *, features=None, r_range=None, split=None, pair=None, qmax=None, qmin=None,
                    window=None, resolve: Optional[Callable[[Path], Path]] = None) -> dict:
    """:func:`assess_run` for a run's fit files (paths in one folder), read with the files beside them.

    Reads ``<stem>.dat`` (:func:`ripple_settings`), ``<stem>_PDFpartials.csv``
    (the closest approach; the pair's g(r)), for the transform of an X-ray
    F(Q) ``<stem>_XFQ<n>.csv`` else ``<stem>_FQ<n>.csv``, for a real-space
    block the input file it names (through ``resolve``, the caller's data-root
    guard) and the Auto StoG / stog files that made it. ``qmax``, ``qmin``,
    ``window``: the caller's word on a real-space dataset whose files do not
    give its transform. A file that cannot be read gives a ``failed`` dataset.
    """
    from .parsers import read_rmc_csv, read_stog_xy

    paths = [Path(path) for path in paths]
    override = {"qMax": qmax, "qMin": qmin, "window": window} if qmax is not None else None
    payloads = []
    cache: dict = {}

    def settings_for(folder: Path, stem: str):
        key = ("dat", folder, stem)
        if key not in cache:
            path = folder / f"{stem}.dat"
            cache[key] = ripple_settings(path.read_text(encoding="utf-8-sig", errors="replace")) if path.is_file() else None
        return cache[key]

    def partials_for(folder: Path, stem: str):
        key = ("partials", folder, stem)
        if key not in cache:
            path = folder / f"{stem}_PDFpartials.csv"
            cache[key] = read_rmc_csv(path) if path.is_file() else None
        return cache[key]

    pair_series = None
    for path in paths:
        info = dataset_file(path.name)
        if info is None:
            continue
        folder = path.parent
        try:
            settings = settings_for(folder, info["stem"])
            partials = partials_for(folder, info["stem"])
            r_first = closest_approach(partials.data if partials is not None else None,
                                       settings["minimumDistances"] if settings else None)
            if pair is not None and partials is not None and pair_series is None:
                for index, label in enumerate(partials.labels[1:], start=1):
                    parts = sorted(part.strip() for part in label.split("-"))
                    if parts == sorted(part.strip() for part in str(pair).split("-")):
                        pair_series = {"label": label.strip(), "r": partials.data[0], "g": partials.data[index]}
            x, calc, data = _read_columns(path)
            fq = None
            input_columns = None
            data_transform = None
            processing = None
            if info["kind"] == "xpdf":
                for candidate in (f"{info['stem']}_XFQ{info['number']}.csv", f"{info['stem']}_FQ{info['number']}.csv"):
                    if (folder / candidate).is_file():
                        fq = _read_columns(folder / candidate)
                        break
            elif info["kind"] in ("npdf", "xray-real"):
                dataset = _nth_block(settings, BLOCKS[info["kind"]], info["number"])
                if dataset and dataset.get("file"):
                    relative = dataset["file"].replace("\\", "/")
                    input_path: Optional[Path] = folder / relative
                    if resolve is not None:
                        try:
                            input_path = resolve(input_path)
                        except (PermissionError, ValueError):
                            input_path = None
                    text = None
                    if input_path is not None and input_path.is_file():
                        text = input_path.read_text(encoding="utf-8-sig", errors="replace")
                        input_columns = read_stog_xy(input_path)
                    data_transform, processing = data_processing(folder, relative.split("/")[-1], text)
            payloads.append(assess_dataset(
                name=path.name, kind=info["kind"], number=info["number"], settings=settings, x=x, calc=calc,
                data=data, fq=fq, input_columns=input_columns, data_transform=data_transform, override=override,
                r_first=r_first, processing=processing,
            ))
        except (OSError, ValueError) as error:
            message = str(error).replace(str(folder) + "/", "").replace(str(folder) + "\\", "")
            payloads.append({"status": "failed", "reason": f"could not read the run's files: {message}",
                             "file": path.name, "kind": info["kind"], "number": info["number"], "total": True, "block": None,
                             "fit": None, "transform": None, "resolution": None, "peaks": [], "lowR": None,
                             "processing": None, "space": None, "rFirst": None, "_analysis": None})
    return assess_run(payloads, features=features, r_range=r_range, split=split,
                      pair=pair_series if pair is not None else None)


__all__ = [
    "BACKFIT_SWEEPS",
    "DEFAULT_DQ",
    "FILTERED_RATIO",
    "MAX_STRONG_PEAKS",
    "PRESENT_FRACTION",
    "RIPPLE_POSITION_FRACTION",
    "RIPPLE_SIZE_FACTOR",
    "RISKS",
    "RULES",
    "across_datasets",
    "analyse_series",
    "auto_stog_title",
    "classify_feature",
    "classify_features",
    "deviation",
    "find_features",
    "fit_peak",
    "fit_space_risk",
    "gaussian_d",
    "gaussian_sine",
    "is_d_like",
    "kernel",
    "lorch",
    "low_r_level",
    "provenance_transform",
    "q_grid",
    "resolution",
    "ripple_pattern",
    "ripple_settings",
    "split_reading",
    "stog_transform",
    "strong_peaks",
    "transform_weights",
]
