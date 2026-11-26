# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Fourier termination ripples: what a run's fit and the Q range of its data say about the ripples in its G(r).

A G(r) is the sine transform of F(Q) measured from Q_min to Q_max. Cutting the
transform at Q_max convolves every peak with the termination kernel, a sinc of
period 2π/Q_max: the peak broadens and oscillating ripples spread on both sides
of it. RMC is a minimiser: fitted to a G(r) that carries those ripples, with a
model G(r) that does not, the box makes structure that reproduces them. This
module says, per fitted dataset of a run, from the run's files and first
principles only:

* the fit space and how RMCProfile makes the model function it compares
  (established on RMCProfile's own files, docs/algorithms/run-dashboard.md,
  "Termination ripples"), so whether the box is asked to fit ripples:

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

* the transform: Q_min, Q_max, ΔQ, the window (Lorch where the files say so)
  and where they come from (the run's F(Q), an Auto StoG provenance or stog
  input naming the data file, CONVOLVE, the X-ray block, the caller's word);
* the r-resolution: the Nyquist step π/Q_max, the kernel's full width at half
  maximum (3.79/Q_max for no window, computed on the run's own grid and
  window), its main lobe and the ripple period 2π/Q_max.

It does not judge whether a given peak or shoulder is a ripple: that reading
(strong-peak fits, predicted ripples, the low-r level, feature classes) was
not reliable enough to ship and is kept on the branch ``wip/joint-fit``.

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
# Rules
# ---------------------------------------------------------------------------

#: The grid step assumed for a transform whose own step the files do not give.
DEFAULT_DQ = 0.01

#: The input data file maps onto the fitted file (y = s·y_input + c) within
#: this fraction of the fitted file's largest |value|.
INPUT_MAP_TOLERANCE = 1e-5

RISKS = ("high", "matched", "unknown", "none")


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



def _block_summary(dataset: Optional[dict]) -> Optional[dict]:
    if not dataset:
        return None
    return {key: dataset.get(key) for key in ("block", "file", "fitType", "realSpaceFit", "convolve", "convolveQ",
                                              "rCutoff", "qMin", "qMax")}


def assess_dataset(*, name: str, kind: str, number: int, settings: Optional[dict], x, data,
                   fq: Optional[tuple] = None, input_columns: Optional[np.ndarray] = None,
                   data_transform: Optional[dict] = None, override: Optional[dict] = None,
                   processing: Optional[dict] = None) -> dict:
    """One fitted dataset of a run: how its model is made, its transform and its r-resolution.

    ``x``, ``data``: the file's r or Q column and its experiment. ``fq``: an
    xpdf's run F(Q) (Q, model, experiment), from which RMCProfile's transform
    must reproduce the file. ``input_columns``: a real-space block's input
    data file (read_stog_xy), mapped onto the fitted file. ``data_transform``:
    the data's own transform where an Auto StoG provenance or a stog input
    names the file (:func:`provenance_transform`, :func:`stog_transform`).
    ``override``: the caller's ``{qMax, qMin, window}`` for a dataset whose
    files do not say.
    """
    dataset = _nth_block(settings, BLOCKS.get(kind, ""), number)
    real = kind in ("xpdf", "npdf", "xray-real")
    risk_kind = kind if real else "reciprocal"
    payload = {
        "status": "unavailable", "reason": None, "file": name, "kind": kind, "number": int(number),
        "block": _block_summary(dataset), "fit": fit_space_risk(risk_kind, dataset, data_transform),
        "transform": None, "resolution": None, "processing": processing,
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

    rr, yy = x, data
    if kind == "xpdf":
        if fq is None:
            payload["reason"] = "no-fq-file"
            return payload
        q, _, f_exp = (np.asarray(column, dtype=float) for column in fq)
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
        q_min, q_max, rule, window, source, window_known = float(q[0]), float(q[-1]), "rectangle", "none", "run-fq", True
        transform_extra = {"reproduction": {"residual": reproduction["residual"], "offset": reproduction["offset"],
                                            "tolerance": RMC_FT_TOLERANCE}}
    else:
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
        transform_extra = {}
        if input_columns is not None:
            mapped = map_input(x, data, input_columns)
            if mapped is not None:
                rr, yy = mapped["r"], mapped["y"]
                transform_extra = {"input": {"scale": mapped["scale"], "offset": mapped["offset"],
                                             "residual": mapped["residual"]}}
    if not _uniform(rr[np.isfinite(rr) & np.isfinite(yy)]):
        payload.update(status="failed", reason="r-grid-not-uniform")
        return payload
    n = int(round((float(q_max) - float(q_min)) / float(dq))) + 1
    q = q_grid(q_min, dq, n)
    weights = transform_weights(q, dq, rule, window)
    payload.update(
        status="assessed",
        transform={"qMin": float(q[0]), "qMax": float(q[-1]), "dQ": float(dq), "nQ": int(q.size), "rule": rule,
                   "window": window, "source": source, "windowKnown": window_known, **transform_extra},
        resolution=resolution(q, weights),
    )
    return payload


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


def ripples_for_run(paths, *, qmax=None, qmin=None, window=None,
                    resolve: Optional[Callable[[Path], Path]] = None) -> dict:
    """:func:`assess_dataset` for a run's fit files (paths in one folder), read with the files beside them.

    Reads ``<stem>.dat`` (:func:`ripple_settings`), for the transform of an
    X-ray F(Q) ``<stem>_XFQ<n>.csv`` else ``<stem>_FQ<n>.csv``, for a
    real-space block the input file it names (through ``resolve``, the
    caller's data-root guard) and the Auto StoG / stog files that made it.
    ``qmax``, ``qmin``, ``window``: the caller's word on a real-space dataset
    whose files do not give its transform. A file that cannot be read gives a
    ``failed`` dataset. Returns ``{datasets}``, one payload per fit file.
    """
    from .parsers import read_stog_xy

    paths = [Path(path) for path in paths]
    override = {"qMax": qmax, "qMin": qmin, "window": window} if qmax is not None else None
    payloads = []
    cache: dict = {}

    def settings_for(folder: Path, stem: str):
        key = (folder, stem)
        if key not in cache:
            path = folder / f"{stem}.dat"
            cache[key] = ripple_settings(path.read_text(encoding="utf-8-sig", errors="replace")) if path.is_file() else None
        return cache[key]

    for path in paths:
        info = dataset_file(path.name)
        if info is None:
            continue
        folder = path.parent
        try:
            settings = settings_for(folder, info["stem"])
            x, _, data = _read_columns(path)
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
                name=path.name, kind=info["kind"], number=info["number"], settings=settings, x=x, data=data, fq=fq,
                input_columns=input_columns, data_transform=data_transform, override=override, processing=processing,
            ))
        except (OSError, ValueError) as error:
            message = str(error).replace(str(folder) + "/", "").replace(str(folder) + "\\", "")
            payloads.append({"status": "failed", "reason": f"could not read the run's files: {message}",
                             "file": path.name, "kind": info["kind"], "number": info["number"], "block": None,
                             "fit": None, "transform": None, "resolution": None, "processing": None})
    return {"datasets": payloads}


__all__ = [
    "DEFAULT_DQ",
    "RISKS",
    "assess_dataset",
    "auto_stog_title",
    "data_processing",
    "dataset_file",
    "fit_space_risk",
    "kernel",
    "lorch",
    "map_input",
    "provenance_transform",
    "q_grid",
    "resolution",
    "ripple_settings",
    "ripples_for_run",
    "stog_transform",
    "transform_weights",
]
