# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""RMCProfile's Fourier transform of an X-ray F(Q), and its propagated uncertainty.

An RMCProfile run that fits X-ray F(Q) data also writes the data's Fourier
transform, ``<stem>_FT_XFQ<n>.csv`` (columns ``r(A), X_ray-calc,
X_ray_exp_renorm``), which the Run Dashboard draws as its G(r) panel. That file
carries no uncertainty. This module rebuilds the transform from the run's own
files, proves it reproduces the file, and propagates the measured F(Q)'s σ
through it: each r point's ±1σ and the full covariance between r points.

The transform (established on the committed demo run and two local runs, see
docs/algorithms/run-dashboard.md, "The G(r) error ribbon"):

    G(r) = Σ_k ΔQ · Q_k · F(Q_k) · sin(Q_k r)

over the Q points of the run's F(Q) file, ΔQ the uniform step — a rectangle
rule, not the trapezoid of stog and pystog (:mod:`rmc_toolkits.transforms`), no
2/π and no ρ₀: Keen's D(r) times π/2. No Lorch window, no omitted-low-Q term.
The F(Q) it transforms is the ``F(Q)_Expt`` column of ``<stem>_XFQ<n>.csv`` (or
``_FQ<n>.csv``) plus a constant, the run's fitted offset (``FITTED_OFFSET``),
recovered here by least squares.

That column is not the input data file the ``.dat`` names (``FILENAME``): it is
the input renormalised by the run's fitted scale s (``FITTED_SCALE``), shifted
by the offset, and, with ``CONVOLVE``, convolved with the box function of the
configuration's g(r) range R (the last r of ``<stem>_PDFpartials.csv``):

    (K F)(Q) = 1/(π Q) Σ_j ΔQ Q_j F(Q_j) [sin((Q − Q_j)R)/(Q − Q_j) − sin((Q + Q_j)R)/(Q + Q_j)]

so ``F(Q)_Expt = s · K F_input − c`` and the transformed function is
``s · K F_input``. Two reproduction tests establish all of it at once, each to
:data:`RMC_FT_TOLERANCE` of the function's largest value: the transform
reproduces the run's G(r) column from its F(Q) column, and the input file
reproduces the run's F(Q) column through s, c and K. Only then is σ propagated.

σ comes only from the input file's own error column (third column), and only
when it is clean (:func:`rmc_toolkits.scaling_cli.usable_sigma`). The transform
is linear, G = A·F with A_rk = ΔQ Q_k sin(Q_k r), so for independent σ_k

    Cov[G(r), G(r')] = Σ_k (ΔQ Q_k s σ_k)² sin(Q_k r) sin(Q_k r') = ½[P(r − r') − P(r + r')]

with P(x) = Σ_k (ΔQ Q_k s σ_k)² cos(Q_k x), exactly; the ribbon is its diagonal.
The box convolution is left out of this propagation: it removes only noise
components beyond r = R, so at r well inside R it changes σ(r) by a fraction of
a percent (≤ 0.2 % on the demo grids at r ≤ 20 Å). It is checked exactly at
the grid points nearest R (:func:`convolution_sigma_change`), and a run whose
r range reaches so close to R that it would change σ(r) by more than
:data:`CONVOLUTION_CHECK_LIMIT` gets no ribbon.

``web_app/frontend/src/ftUncertainty.js`` is a line-for-line port (parity
golden ``ft_uncertainty_fixture.json``, ``tests/generate_ft_uncertainty_fixture.py``).
"""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Callable, Optional

import numpy as np

#: A reproduction test passes when its largest deviation is at most this
#: fraction of the reproduced function's largest |value|. The demo's 15-digit
#: files reproduce to 1.5e-14; its 7-decimal F(Q) file to 4.5e-8. The smallest
#: convention error seen (the run's offset left out) is 3e-4.
RMC_FT_TOLERANCE = 1e-6

#: The box convolution may change σ(r) by at most this fraction at the checked
#: points before the ribbon is refused.
CONVOLUTION_CHECK_LIMIT = 0.01

#: A Q grid is uniform when every point lies within this fraction of the step
#: of ``q_0 + k·ΔQ``.
UNIFORM_GRID_RTOL = 1e-6

#: The checked points: r_max, and this many more spread over the last π/Q_max
#: (one period of the correction's oscillation), plus r_min.
CONVOLUTION_CHECK_POINTS = 8

#: The run's transform file and its dataset number.
FT_FILE_RE = re.compile(r"^(.+)_FT_XFQ(\d+)\.csv$")

_SIGMA_WORDS_RE = re.compile(r"sigma|σ|uncert|error")
_LABEL_SPLIT_RE = re.compile(r"[^a-z0-9σ]+")
_NUMBER_RE = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
_LINE_BREAK_RE = re.compile(r"\r\n|\r|\n")
_DAT_KEY_RE = re.compile(r"^([A-Z_\d]+)\s*::\s*(.*)$")
_DAT_SUB_RE = re.compile(r"^([A-Z_]+)\s*(?:::\s*(.*))?$")


# ---------------------------------------------------------------------------
# The transform and its uncertainty (pure numerics)
# ---------------------------------------------------------------------------


def uniform_step(q: np.ndarray, name: str = "Q") -> float:
    """The step ΔQ of a uniform, strictly increasing grid, or ``ValueError``.

    ΔQ = (q_last − q_0)/(n − 1); every point must lie within
    ``UNIFORM_GRID_RTOL · ΔQ`` of ``q_0 + k ΔQ``.
    """
    q = np.asarray(q, dtype=float)
    if q.ndim != 1 or q.size < 2:
        raise ValueError(f"the {name} grid needs at least two points")
    if not (np.all(np.isfinite(q)) and np.all(np.diff(q) > 0)):
        raise ValueError(f"the {name} grid is not strictly increasing")
    step = float((q[-1] - q[0]) / (q.size - 1))
    ideal = q[0] + step * np.arange(q.size)
    if float(np.max(np.abs(q - ideal))) > UNIFORM_GRID_RTOL * step:
        raise ValueError(f"the {name} grid is not uniform")
    return step


def rmc_ft_coefficients(q: np.ndarray) -> np.ndarray:
    """The operator's weights c_k = ΔQ · Q_k: G(r) = Σ_k c_k F_k sin(Q_k r)."""
    q = np.asarray(q, dtype=float)
    return uniform_step(q) * q


def _sines(r: np.ndarray, q: np.ndarray) -> np.ndarray:
    return np.sin(np.outer(np.asarray(r, dtype=float), np.asarray(q, dtype=float)))


def rmc_ft(q: np.ndarray, f: np.ndarray, r: np.ndarray) -> np.ndarray:
    """RMCProfile's transform of F(Q) onto ``r``: Σ_k ΔQ Q_k F_k sin(Q_k r)."""
    weights = rmc_ft_coefficients(q) * np.asarray(f, dtype=float)
    return _sines(r, q) @ weights


def rmc_ft_sigma(q: np.ndarray, sigma: np.ndarray, r: np.ndarray) -> np.ndarray:
    """σ of the transform at each r for independent σ_k of F(Q_k): √Σ_k (c_k σ_k sin Q_k r)²."""
    power = (rmc_ft_coefficients(q) * np.asarray(sigma, dtype=float)) ** 2
    return np.sqrt((_sines(r, q) ** 2) @ power)


def noise_autocorrelation(q: np.ndarray, sigma: np.ndarray, x: np.ndarray) -> np.ndarray:
    """P(x) = Σ_k (c_k σ_k)² cos(Q_k x), the covariance's Toeplitz kernel."""
    power = (rmc_ft_coefficients(q) * np.asarray(sigma, dtype=float)) ** 2
    return np.cos(np.outer(np.asarray(x, dtype=float), np.asarray(q, dtype=float))) @ power


def _uniform_r_step(r: np.ndarray) -> Optional[float]:
    """The step of a uniform r grid (within 1e-7 of a step), else None (MATERIA's rule)."""
    if r.size < 2:
        return None
    step = float((r[-1] - r[0]) / (r.size - 1))
    if not step > 0:
        return None
    ideal = r[0] + step * np.arange(r.size)
    return step if float(np.max(np.abs(r - ideal))) <= 1e-7 * step else None


def rmc_ft_covariance(q: np.ndarray, sigma: np.ndarray, r: np.ndarray) -> np.ndarray:
    """The full covariance Cov[G(r_i), G(r_j)] (n × n) for independent σ_k of F(Q_k).

    Uniform r grid: C_ij = ½[P((i − j)h) − P(2 r_0 + (i + j)h)] from O(n) values
    of P; otherwise the direct sum Σ_k (c_k σ_k)² sin(Q_k r_i) sin(Q_k r_j). The
    diagonal is :func:`rmc_ft_sigma` squared, the same numbers.
    """
    q = np.asarray(q, dtype=float)
    r = np.asarray(r, dtype=float)
    n = r.size
    step = _uniform_r_step(r)
    if step is not None:
        toeplitz = noise_autocorrelation(q, sigma, step * np.arange(n))
        hankel = noise_autocorrelation(q, sigma, 2.0 * r[0] + step * np.arange(2 * n - 1))
        index = np.arange(n)
        cov = 0.5 * (toeplitz[np.abs(index[:, None] - index[None, :])] - hankel[index[:, None] + index[None, :]])
    else:
        sines = _sines(r, q)
        power = (rmc_ft_coefficients(q) * np.asarray(sigma, dtype=float)) ** 2
        cov = (sines * power) @ sines.T
    np.fill_diagonal(cov, rmc_ft_sigma(q, sigma, r) ** 2)
    return cov


def correlation_summary(q: np.ndarray, sigma: np.ndarray, r: np.ndarray) -> dict:
    """How correlated the propagated errors of neighbouring r points are.

    ``neighbour`` and ``nyquist`` are the stationary correlations P(Δ)/P(0) at the
    grid step and at the Nyquist step π/Q_max; ``effectivePoints`` is the
    participation ratio n²/‖R‖²_F of the exact correlation matrix (Bretherton et
    al. 1999), ``nyquistPoints`` the sampling-theorem count (r_max − r_min)·Q_max/π
    (Farrow et al. 2011). Ported from MATERIA's grErrors.ts.
    """
    q = np.asarray(q, dtype=float)
    r = np.asarray(r, dtype=float)
    n = r.size
    grid_step = float((r[-1] - r[0]) / (n - 1)) if n > 1 else 0.0
    qmax = float(q[-1])
    nyquist_step = math.pi / qmax
    p0, p_grid, p_nyq = noise_autocorrelation(q, sigma, np.array([0.0, grid_step, nyquist_step]))
    cov = rmc_ft_covariance(q, sigma, r)
    sd = np.sqrt(np.diag(cov))
    keep = sd > 0
    m = int(np.count_nonzero(keep))
    if m:
        corr = cov[np.ix_(keep, keep)] / np.outer(sd[keep], sd[keep])
        frobenius = float(np.sum(corr * corr))
        effective = m * m / frobenius if frobenius > 0 else 0.0
    else:
        effective = 0.0
    return {
        "points": int(n),
        "gridStep": grid_step,
        "nyquistStep": nyquist_step,
        "neighbour": float(p_grid / p0) if p0 > 0 else 0.0,
        "nyquist": float(p_nyq / p0) if p0 > 0 else 0.0,
        "effectivePoints": float(effective),
        "nyquistPoints": float(max(0.0, r[-1] - r[0]) * qmax / math.pi) if n else 0.0,
    }


def _box_kernel_terms(q: np.ndarray, radius: float) -> tuple[np.ndarray, np.ndarray]:
    """T(m) = sin(m ΔQ R)/(m ΔQ) (R at m = 0) and H(s) = sin((2Q_0 + sΔQ)R)/(2Q_0 + sΔQ)."""
    q = np.asarray(q, dtype=float)
    step = uniform_step(q)
    if not q[0] > 0:
        raise ValueError("the box convolution needs Q > 0 (it divides by Q)")
    n = q.size
    lag = step * np.arange(n)
    toeplitz = np.empty(n)
    toeplitz[0] = float(radius)
    toeplitz[1:] = np.sin(lag[1:] * radius) / lag[1:]
    total = 2.0 * q[0] + step * np.arange(2 * n - 1)
    hankel = np.sin(total * radius) / total
    return toeplitz, hankel


def _box_kernel_apply(v: np.ndarray, toeplitz: np.ndarray, hankel: np.ndarray) -> np.ndarray:
    """Σ_k v_k [T(|k − j|) − H(k + j)] for every j (T symmetric, H Hankel)."""
    n = v.size
    full = np.concatenate([toeplitz[:0:-1], toeplitz])
    toeplitz_part = np.convolve(v, full)[n - 1 : 2 * n - 1]
    hankel_part = np.correlate(hankel, v, mode="valid")
    return toeplitz_part - hankel_part


def box_convolution(q: np.ndarray, f: np.ndarray, radius: float) -> np.ndarray:
    """RMCProfile's CONVOLVE: F(Q) convolved with the box function of radius R.

    (K F)(Q_k) = 1/(π Q_k) Σ_j ΔQ Q_j F_j [T(|k − j|) − H(k + j)]: the sine
    transform of F truncated at r = R, transformed back, on the same grid.
    """
    q = np.asarray(q, dtype=float)
    toeplitz, hankel = _box_kernel_terms(q, radius)
    weights = rmc_ft_coefficients(q) * np.asarray(f, dtype=float)
    return _box_kernel_apply(weights, toeplitz, hankel) / (math.pi * q)


def convolution_check_points(r: np.ndarray, qmax: float) -> np.ndarray:
    """r_min, and CONVOLUTION_CHECK_POINTS points over the last π/Q_max up to r_max."""
    r = np.asarray(r, dtype=float)
    r_min, r_max = float(r[0]), float(r[-1])
    span = min(math.pi / qmax, r_max - r_min)
    count = CONVOLUTION_CHECK_POINTS
    tail = [r_max - span * k / (count - 1) for k in range(count - 1, -1, -1)]
    return np.array([r_min] + tail)


def convolution_sigma_change(q: np.ndarray, sigma: np.ndarray, r_points: np.ndarray, radius: float) -> np.ndarray:
    """σ(r) with the box convolution over σ(r) without it, minus 1, at each point.

    With it, the operator from the input to G is B = A·K: B_rj = (ΔQ² Q_j/π)
    Σ_k sin(Q_k r)[T(|k − j|) − H(k + j)], and Var = Σ_j (B_rj σ_j)²; without
    it, Var = Σ_k (ΔQ Q_k σ_k sin Q_k r)².
    """
    q = np.asarray(q, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    step = uniform_step(q)
    toeplitz, hankel = _box_kernel_terms(q, radius)
    coefficients = step * q
    out = []
    for r in np.asarray(r_points, dtype=float):
        sines = np.sin(q * r)
        exact = (step * step * q / math.pi) * _box_kernel_apply(sines, toeplitz, hankel)
        with_box = float(np.sum((exact * sigma) ** 2))
        without = float(np.sum((coefficients * sines * sigma) ** 2))
        out.append(math.sqrt(with_box / without) - 1.0 if without > 0 else 0.0)
    return np.array(out)


# ---------------------------------------------------------------------------
# The two reproduction tests
# ---------------------------------------------------------------------------


def _deviation(target: np.ndarray, model: np.ndarray, tolerance: float) -> dict:
    residual = target - model
    scale = float(np.max(np.abs(target)))
    largest = float(np.max(np.abs(residual)))
    relative = largest / scale if scale > 0 else (0.0 if largest == 0 else math.inf)
    return {
        "residual": relative,
        "rms": float(math.sqrt(float(np.mean(residual * residual)))),
        "pass": bool(relative <= tolerance),
    }


def reproduce_transform(q, f, r, g, tolerance: float = RMC_FT_TOLERANCE) -> dict:
    """The run's G(r) column from its F(Q) column: G ≈ A(F + c), c by least squares.

    Returns ``{residual, rms, offset, pass}``: ``residual`` the largest |G − A(F + c)|
    over max |G|, ``offset`` the constant c (the run's fitted offset).
    """
    q = np.asarray(q, dtype=float)
    r = np.asarray(r, dtype=float)
    g = np.asarray(g, dtype=float)
    sines = _sines(r, q)
    coefficients = rmc_ft_coefficients(q)
    base = sines @ (coefficients * np.asarray(f, dtype=float))
    unit = sines @ coefficients
    norm = float(np.dot(unit, unit))
    offset = float(np.dot(g - base, unit)) / norm if norm > 0 else 0.0
    return {**_deviation(g, base + offset * unit, tolerance), "offset": offset}


def reproduce_input(q, f_run, f_input, *, radius: Optional[float] = None, tolerance: float = RMC_FT_TOLERANCE) -> dict:
    """The run's F(Q) column from the input file: F_run ≈ s·(K F_input) + c'.

    K is the box convolution when ``radius`` is given, else the identity; s and c'
    by least squares. Returns ``{residual, rms, scale, offset, pass}``.
    """
    q = np.asarray(q, dtype=float)
    f_run = np.asarray(f_run, dtype=float)
    base = box_convolution(q, f_input, radius) if radius is not None else np.asarray(f_input, dtype=float)
    n = float(base.size)
    sx = float(np.sum(base))
    sy = float(np.sum(f_run))
    sxx = float(np.dot(base, base))
    sxy = float(np.dot(base, f_run))
    det = n * sxx - sx * sx
    if det == 0:
        scale, offset = 0.0, sy / n
    else:
        scale = (n * sxy - sx * sy) / det
        offset = (sxx * sy - sx * sxy) / det
    return {**_deviation(f_run, scale * base + offset, tolerance), "scale": scale, "offset": offset}


# ---------------------------------------------------------------------------
# The run-control file
# ---------------------------------------------------------------------------


def _leading_numbers(value: str) -> list[float]:
    """The plain decimal numbers a value starts with (``1 2619 1 Angstrom`` → [1, 2619, 1])."""
    numbers: list[float] = []
    for token in value.split():
        if not _NUMBER_RE.match(token):
            break
        numbers.append(float(token))
    return numbers


def run_datasets(text: str) -> dict:
    """The ``.dat``'s data blocks, as the transform needs them.

    ``{rSpacing, numberDensity, datasets}``; each dataset ``{block, file,
    dataType, fitType, convolve, fittedScale, fittedOffset, reciprocalParameters,
    realParameters}`` in file order. ``> CONVOLVE`` (with or without ``::``) and
    ``> FITTED_SCALE`` / ``> FITTED_OFFSET`` are flags. Mirrors ``runDatasets``.
    """
    datasets: list[dict] = []
    block: Optional[dict] = None
    r_spacing = None
    number_density = None
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
            elif key == "DATA_TYPE":
                block["dataType"] = value
            elif key == "FIT_TYPE":
                block["fitType"] = value
            elif key == "CONVOLVE":
                block["convolve"] = True
            elif key == "FITTED_SCALE":
                block["fittedScale"] = True
            elif key == "FITTED_OFFSET":
                block["fittedOffset"] = True
            elif key == "RECIPROCAL_SPACE_PARAMETERS":
                block["reciprocalParameters"] = _leading_numbers(value)
            elif key == "REAL_SPACE_PARAMETERS":
                block["realParameters"] = _leading_numbers(value)
            continue
        match = _DAT_KEY_RE.match(line)
        block = None
        if not match:
            continue
        key, value = match.groups()
        if key == "END":
            break
        if key.endswith("_DATA"):
            block = {
                "block": key, "file": None, "dataType": None, "fitType": None, "convolve": False,
                "fittedScale": False, "fittedOffset": False, "reciprocalParameters": None, "realParameters": None,
            }
            datasets.append(block)
        elif key == "R_SPACING":
            numbers = _leading_numbers(value)
            r_spacing = numbers[0] if numbers else None
        elif key == "NUMBER_DENSITY":
            numbers = _leading_numbers(value)
            number_density = numbers[0] if numbers else None
    return {"rSpacing": r_spacing, "numberDensity": number_density, "datasets": datasets}


def xray_dataset(settings: dict, number: int) -> Optional[dict]:
    """The ``number``-th (1-based) ``XRAY_RECIPROCAL_SPACE_DATA`` block, or None."""
    blocks = [entry for entry in settings.get("datasets", []) if entry["block"] == "XRAY_RECIPROCAL_SPACE_DATA"]
    return blocks[number - 1] if 0 < number <= len(blocks) else None


def is_sigma_label(label) -> bool:
    """A header that names an uncertainty: sigma, σ, uncert…, error, or the word err / esd."""
    text = str(label).lower()
    if _SIGMA_WORDS_RE.search(text):
        return True
    return any(word in ("err", "esd") for word in _LABEL_SPLIT_RE.split(text))


def file_sigma_column(labels: list[str]) -> Optional[int]:
    """The index of a column whose header names an uncertainty (never the x column), or None."""
    for index, label in enumerate(labels[1:], start=1):
        if is_sigma_label(label):
            return index
    return None


# ---------------------------------------------------------------------------
# The assessment
# ---------------------------------------------------------------------------


def _usable_sigma(q: np.ndarray, f: np.ndarray, sigma: np.ndarray) -> Optional[np.ndarray]:
    """scaling_cli.usable_sigma's rule: σ finite and > 0 on every row with finite Q and F."""
    usable = np.isfinite(q) & np.isfinite(f)
    if not np.all(np.isfinite(sigma[usable])) or np.any(sigma[usable] <= 0):
        return None
    return sigma


def _matching_block(q_input: np.ndarray, q: np.ndarray, step: float) -> Optional[slice]:
    """The run's Q grid as a contiguous block of the input file's, or None."""
    tolerance = UNIFORM_GRID_RTOL * step
    hits = np.flatnonzero(np.abs(q_input - q[0]) <= tolerance)
    for start in hits:
        stop = int(start) + q.size
        if stop <= q_input.size and float(np.max(np.abs(q_input[start:stop] - q))) <= tolerance:
            return slice(int(start), stop)
    return None


def assess_ft_uncertainty(
    *,
    r,
    g,
    q=None,
    f=None,
    settings: Optional[dict] = None,
    dataset_number: int = 1,
    input_columns: Optional[np.ndarray] = None,
    radius: Optional[float] = None,
    names: Optional[dict] = None,
    tolerance: float = RMC_FT_TOLERANCE,
) -> dict:
    """Whether, and how, the G(r) file's error ribbon follows from the run's files.

    ``r``, ``g``: the transform file's r and experimental columns. ``q``, ``f``: the
    run's F(Q) file's Q and ``F(Q)_Expt`` columns (None when there is none).
    ``settings``: :func:`run_datasets` of ``<stem>.dat`` (None when absent).
    ``input_columns``: the input file's columns (Q, F[, σ, ...]) as
    :func:`rmc_toolkits.parsers.read_stog_xy` gives them (None when absent).
    ``radius``: the box radius R, the last r of ``<stem>_PDFpartials.csv``.
    ``names``: the files' names for the messages and ``sources`` (``transform``,
    ``fq``, ``settings``, ``input``, ``partials``; ``fqCandidates`` names the F(Q)
    files looked for, in a message only).

    Returns ``{status, reason, sigma, transform, reproduction,
    convolutionCheck, correlation, sources}``: ``status`` is ``propagated``
    (``sigma`` is σ(r) on the file's r grid), ``unavailable`` (an input the
    ribbon needs is not there: no σ, a missing file) or ``failed`` (a
    reproduction test or a check failed); ``reason`` says why there is no
    ribbon, without numbers (they are in ``reproduction`` and
    ``convolutionCheck``).
    """
    names = dict(names or {})
    transform_name = names.get("transform") or "the transform file"
    fq_name = names.get("fq")
    settings_name = names.get("settings") or "the run-control .dat"
    missing_settings = (
        f"no run-control file {names['settings']} beside {transform_name}" if names.get("settings")
        else f"no run-control .dat beside {transform_name}")
    input_name = names.get("input")
    partials_name = names.get("partials") or "the run's PDFpartials file"
    fq_candidates = names.get("fqCandidates") or "the run's F(Q) file"
    # The files actually read, for the tooltip.
    used = {
        "transform": True,
        "fq": q is not None,
        "settings": settings is not None,
        "input": input_columns is not None,
        "partials": radius is not None,
    }
    sources = {key: names[key] for key, present in used.items() if present and names.get(key)}

    r = np.asarray(r, dtype=float)
    g = np.asarray(g, dtype=float)
    payload = {
        "status": "unavailable",
        "reason": None,
        "sigma": None,
        "transform": None,
        "reproduction": {"tolerance": tolerance, "gr": None, "fq": None},
        "convolutionCheck": None,
        "correlation": None,
        "sources": sources,
    }

    # Why σ cannot be had at all, before anything is computed.
    dataset = xray_dataset(settings, dataset_number) if settings is not None else None
    sigma_input = None
    if settings is None:
        no_sigma = missing_settings
    elif dataset is None:
        no_sigma = f"{settings_name} names no X-ray reciprocal-space dataset {dataset_number}"
    elif not dataset.get("file"):
        no_sigma = f"{settings_name} names no data file for its X-ray dataset {dataset_number}"
    elif input_columns is None:
        no_sigma = f"F(Q) carries no uncertainties: {dataset['file']} (named in {settings_name}) is not in the run folder"
    elif input_columns.shape[0] < 3:
        no_sigma = f"F(Q) carries no uncertainties: {input_name} has no error column"
    else:
        sigma_input = _usable_sigma(input_columns[0], input_columns[1], input_columns[2])
        no_sigma = None if sigma_input is not None else (
            f"F(Q) carries no uncertainties: the error column of {input_name} is not usable "
            "(a zero, negative or non-finite σ)")
    convolve = bool(dataset and dataset.get("convolve"))

    problems: list[str] = []
    if not (r.size >= 2 and np.all(np.isfinite(r)) and np.all(np.diff(r) > 0)):
        problems.append(f"the r column of {transform_name} is not strictly increasing")
    if not np.all(np.isfinite(g)):
        problems.append(f"{transform_name} has non-finite values")

    step = None
    if q is not None:
        q = np.asarray(q, dtype=float)
        f = np.asarray(f, dtype=float)
        try:
            step = uniform_step(q)
        except ValueError as error:
            problems.append(f"{fq_name}: {error}")
        if not np.all(np.isfinite(f)):
            problems.append(f"{fq_name} has non-finite F(Q) values")
        if step is not None and not problems:
            payload["transform"] = {
                "qMin": float(q[0]), "qMax": float(q[-1]), "dQ": step, "nQ": int(q.size),
                "rMin": float(r[0]), "rMax": float(r[-1]), "nR": int(r.size),
                "convolve": convolve, "radius": float(radius) if (convolve and radius is not None) else None,
            }
            payload["reproduction"]["gr"] = reproduce_transform(q, f, r, g, tolerance)

    # The input file through the chain (also without σ: the tooltip says it).
    f_input = None
    if payload["reproduction"]["gr"] is not None and input_columns is not None and input_columns.shape[0] >= 2:
        block = _matching_block(np.asarray(input_columns[0], dtype=float), q, step)
        if block is None:
            problems.append(f"the Q grid of {input_name} does not hold the Q grid of {fq_name}")
        else:
            f_input = np.asarray(input_columns[1], dtype=float)[block]
            if not np.all(np.isfinite(f_input)):
                problems.append(f"{input_name} has non-finite F(Q) values on the run's Q grid")
                f_input = None
            elif not (convolve and radius is None):
                payload["reproduction"]["fq"] = {
                    **reproduce_input(q, f, f_input, radius=radius if convolve else None, tolerance=tolerance),
                    "convolve": convolve,
                }
            if sigma_input is not None:
                sigma_input = np.asarray(sigma_input, dtype=float)[block]

    gr = payload["reproduction"]["gr"]
    fq = payload["reproduction"]["fq"]
    if no_sigma is not None:
        payload["reason"] = no_sigma
        return payload
    if problems:
        payload.update(status="failed", reason=problems[0])
        return payload
    if q is None:
        payload["reason"] = f"no F(Q) file ({fq_candidates}) to reproduce the transform from"
        return payload
    if convolve and radius is None:
        payload["reason"] = (
            f"{partials_name} is not in the run folder: the radius of the box convolution RMCProfile "
            "applied is unknown")
        return payload
    if not gr["pass"]:
        payload.update(status="failed", reason=(
            f"RMCProfile's transform of {fq_name} does not reproduce {transform_name} within the tolerance"))
        return payload
    if not fq["pass"]:
        how = " and box convolution" if convolve else ""
        payload.update(status="failed", reason=(
            f"{input_name} does not reproduce {fq_name} through the run's renormalisation{how} within the tolerance"))
        return payload
    if not fq["scale"] > 0:
        payload.update(status="failed", reason=f"the renormalisation of {input_name} the run fitted is not positive")
        return payload

    sigma_q = fq["scale"] * sigma_input
    if convolve:
        points = convolution_check_points(r, float(q[-1]))
        changes = convolution_sigma_change(q, sigma_q, points, radius)
        worst = int(np.argmax(np.abs(changes)))
        payload["convolutionCheck"] = {
            "radius": float(radius),
            "limit": CONVOLUTION_CHECK_LIMIT,
            "maxChange": float(abs(changes[worst])),
            "at": float(points[worst]),
            "points": int(points.size),
        }
        if abs(changes[worst]) > CONVOLUTION_CHECK_LIMIT:
            payload.update(status="failed", reason=(
                f"the r range of {transform_name} reaches too close to the box-convolution radius: left out of the "
                "propagation, the convolution would change σ(r) beyond the limit"))
            return payload

    payload["sigma"] = rmc_ft_sigma(q, sigma_q, r).tolist()
    payload["correlation"] = correlation_summary(q, sigma_q, r)
    payload.update(status="propagated", reason=None)
    return payload


# ---------------------------------------------------------------------------
# The run folder
# ---------------------------------------------------------------------------


def ft_failure(message: str) -> dict:
    """The payload of a computation that raised: no ribbon, the message as the reason (``ftFailure``)."""
    return {
        "status": "failed",
        "reason": message,
        "sigma": None,
        "transform": None,
        "reproduction": {"tolerance": RMC_FT_TOLERANCE, "gr": None, "fq": None},
        "convolutionCheck": None,
        "correlation": None,
        "sources": {},
    }


def _experimental_column(labels: list[str], data: np.ndarray) -> np.ndarray:
    from .parsers import rwp_columns

    roles = rwp_columns(labels, len(data))
    if roles is None:
        raise ValueError("needs an experimental column beside the calculated one")
    return data[roles[1]]


def ft_uncertainty_for_run(
    ft_path: str | Path,
    *,
    resolve: Optional[Callable[[Path], Path]] = None,
    tolerance: float = RMC_FT_TOLERANCE,
) -> dict:
    """:func:`assess_ft_uncertainty` for ``<stem>_FT_XFQ<n>.csv`` from the files beside it.

    Reads ``<stem>.dat`` (the n-th X-ray dataset's ``FILENAME`` and ``CONVOLVE``),
    the input file it names (resolved against the folder, then through
    ``resolve``, the caller's data-root guard, which raises for a path it
    refuses), ``<stem>_XFQ<n>.csv`` (else ``<stem>_FQ<n>.csv``) and, for the box
    radius, ``<stem>_PDFpartials.csv``. A file that cannot be read or parsed
    raises (``ValueError``/``OSError``); the caller reports it.
    """
    from .parsers import read_rmc_csv, read_stog_xy

    ft_path = Path(ft_path)
    match = FT_FILE_RE.match(ft_path.name)
    if not match:
        raise ValueError(f"{ft_path.name} is not an RMCProfile transform file (<stem>_FT_XFQ<n>.csv)")
    stem, number = match.group(1), int(match.group(2))
    folder = ft_path.parent
    candidates = [f"{stem}_XFQ{number}.csv", f"{stem}_FQ{number}.csv"]
    names: dict = {
        "transform": ft_path.name,
        "settings": f"{stem}.dat",
        "fqCandidates": f"{candidates[0]} or {candidates[1]}",
    }

    transform = read_rmc_csv(ft_path)
    r = transform.data[0]
    g = _experimental_column(transform.labels, transform.data)

    q = f = None
    for candidate in candidates:
        path = folder / candidate
        if path.is_file():
            series = read_rmc_csv(path)
            q, f = series.data[0], _experimental_column(series.labels, series.data)
            names["fq"] = candidate
            break

    settings = None
    dat_path = folder / f"{stem}.dat"
    if dat_path.is_file():
        settings = run_datasets(dat_path.read_text(encoding="utf-8-sig", errors="replace"))
    dataset = xray_dataset(settings, number) if settings is not None else None

    input_columns = None
    radius = None
    if dataset is not None and dataset.get("file"):
        relative = dataset["file"].replace("\\", "/")
        names["input"] = relative.split("/")[-1]
        input_path: Optional[Path] = folder / relative
        if resolve is not None:
            try:
                input_path = resolve(input_path)
            except (PermissionError, ValueError):
                input_path = None  # outside what the caller may read: as if absent
        if input_path is not None and input_path.is_file():
            input_columns = read_stog_xy(input_path)
        if dataset.get("convolve"):
            partials_path = folder / f"{stem}_PDFpartials.csv"
            names["partials"] = partials_path.name
            if partials_path.is_file():
                radius = float(np.nanmax(read_rmc_csv(partials_path).data[0]))

    return assess_ft_uncertainty(
        r=r, g=g, q=q, f=f, settings=settings, dataset_number=number, input_columns=input_columns,
        radius=radius, names=names, tolerance=tolerance,
    )


__all__ = [
    "CONVOLUTION_CHECK_LIMIT",
    "CONVOLUTION_CHECK_POINTS",
    "FT_FILE_RE",
    "RMC_FT_TOLERANCE",
    "UNIFORM_GRID_RTOL",
    "box_convolution",
    "convolution_check_points",
    "convolution_sigma_change",
    "correlation_summary",
    "file_sigma_column",
    "ft_failure",
    "is_sigma_label",
    "assess_ft_uncertainty",
    "ft_uncertainty_for_run",
    "noise_autocorrelation",
    "reproduce_input",
    "reproduce_transform",
    "rmc_ft",
    "rmc_ft_coefficients",
    "rmc_ft_covariance",
    "rmc_ft_sigma",
    "run_datasets",
    "uniform_step",
    "xray_dataset",
]
