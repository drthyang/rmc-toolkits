# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Plot builders for RMCProfile outputs and legacy STOG preprocessing files."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import io
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .parsers import fit_rwp, pdf_index, read_chi_log, read_exafs_csv, read_rmc_csv, read_stog, related_r_value_logs


@dataclass(frozen=True)
class PlotResult:
    figure: plt.Figure
    kind: str
    title: str
    # An entry may be None when the metric is not defined for the data (see rwp).
    metrics: dict[str, float | None] = field(default_factory=dict)


def detect_plot_kind(path: str | Path) -> str | None:
    name = Path(path).name
    if re.search(r"-EXAFS-.+_Q_OUTPUT\.csv$", name):
        return "exafs_q"
    if re.search(r"-EXAFS-.+_R_OUTPUT\.csv$", name):
        return "exafs_r"
    if re.search(r"_FT_XFQ\d+\.csv$", name):
        return "xpdf"
    if "PDF" in name and name.endswith(".csv"):
        return "pdf_partials" if "PDFpartials" in name else "npdf"
    # Reciprocal-space fits of any dataset number (_FQ2, _SQ3, …), like the
    # xPDF rule and the stem chooser; *_FQ1partials.csv is not one of them.
    if re.search(r"_FQ\d+\.csv$", name):
        return "xray_sq"
    if re.search(r"_SQ\d+\.csv$", name):
        return "neutron_sq"
    if re.search(r"_bragg(?:_.+)?\.csv$", name):
        return "bragg"
    if re.search(r"-\d{2,}\.log$", name):
        return "r_value"
    if name in {"scale_ft.gr", "scale_ft.sq", "scale_ft_rmc.fq"}:
        return "stog"
    return None


_FIT_FUNCTION_RE = re.compile(r"([A-Za-z])\(([QqRr])\)")
_RECIPROCAL_DATASET_RE = re.compile(r"_([FS])Q(\d+)\.csv$")


def fit_function_label(labels: list[str]) -> str | None:
    """The function a fit CSV holds, read from its own data-column headers.

    ``Q, F(Q)_RMC, F(Q)_Expt`` → ``F(Q)``; ``r, G(r)_RMC, …`` → ``G(r)``. ``None``
    when no data column names one. Mirrors ``fitFunctionLabel()`` in browserData.js.
    """
    for label in labels[1:]:
        match = _FIT_FUNCTION_RE.search(str(label))
        if match:
            return f"{match.group(1)}({match.group(2)})"
    return None


def series_titles(kind: str, name: str, labels: list[str]) -> tuple[str, str]:
    """``(title, y label)`` of a dashboard fit/partials CSV, from what the file holds.

    * ``*_FQn.csv`` / ``*_SQn.csv``: RMCProfile writes F(Q) (``F(Q)_RMC, F(Q)_Expt``,
      → −⟨b⟩² at low Q and 0 at high Q), not S(Q) (→ 1). The function is read from
      the column headers, else from the file name (``FQ`` → F(Q), ``SQ`` → S(Q)); the
      title is that function, with ``#n`` for dataset n > 1, and claims no radiation.
    * ``*_PDFpartials.csv``: the partial pair distribution functions g_ij(r) (0
      below the closest approach, → 1 at large r), not G(r) (→ 0).
    * ``xpdf`` / ``npdf``: the header's function if it names one, else ``G(r)``.

    Mirrors ``seriesTitles()`` in browserData.js; ``/api/plot/data`` and the
    matplotlib figures use it too.
    """
    if kind in ("xray_sq", "neutron_sq"):
        match = _RECIPROCAL_DATASET_RE.search(name)
        letter, index = (match.group(1), int(match.group(2))) if match else ("F" if kind == "xray_sq" else "S", 1)
        function = fit_function_label(labels) or f"{letter}(Q)"
        return (function if index == 1 else f"{function} #{index}"), function
    if kind == "pdf_partials":
        return "Partial g(r)", "g(r)"
    if kind == "xpdf":
        return "xPDF", fit_function_label(labels) or "G(r)"
    if kind == "npdf":
        return Path(name).stem.split("_")[-1], fit_function_label(labels) or "G(r)"
    if kind == "exafs_q":
        return "EXAFS Q-space", "χ(k) k²"
    if kind == "exafs_r":
        return "EXAFS R-space", "FT[χ(k) k²]"
    if kind == "bragg":
        return "BRAGG", "Intensity"
    return name, "data"


def stog_function_label(name: str) -> str:
    """Name default for a STOG file: ``.gr`` G(r) or g(r), ``.fq`` F(Q), else S(Q).

    Classic stog ``scale.gr`` and the filtered ``*_ft.gr`` (``scale_ft.gr``,
    rmc-autoscale's ``<stem>_ft.gr``) hold g(r), → 1 at large r; ``*_rmc.gr``
    is Keen's G_K(r) and any other ``.gr`` keeps the generic G(r).
    ``scale_ft_rmc.fq`` is Keen's F(Q) (→ 0 at high Q); only ``.sq`` holds S(Q).
    The browser (``stogFunctionLabel``, same rule) prefers the run-control
    ``FIT_TYPE`` when it knows it.
    """
    lower = name.lower()
    if lower.endswith(".gr"):
        return "g(r)" if lower == "scale.gr" or lower.endswith("_ft.gr") else "G(r)"
    if lower.endswith(".fq"):
        return "F(Q)"
    return "S(Q)"


def bragg_is_tof(header: str | None) -> bool:
    """True when a Bragg dataset's column header is time-of-flight.

    RMCProfile time-of-flight Bragg CSVs label the first column ``Flight time (us)``;
    those are shown as ToF in microseconds (the conventional neutron unit), so the
    raw x-values are used as-is. Anything else (e.g. ``Q or theta``) stays a Q axis.
    """
    return bool(re.search(r"tof|flight|time", (header or "").lower()))


def _series_plot(
    path: Path,
    title: str,
    xlabel: str,
    ylabel: str = "data",
    calculate_rwp: bool = False,
    reader=read_rmc_csv,
) -> PlotResult:
    series = reader(path)
    if len(series.data) < 2:
        raise ValueError(f"{path} needs at least two numeric columns")

    metrics: dict[str, float | None] = {}
    if calculate_rwp and len(series.data) >= 3:
        # Normalized by the experiment: RMCProfile writes (x, calculated,
        # experimental), and a header that names the roles overrides that order.
        metrics["rwp"] = fit_rwp(series.labels, series.data)

    fig = plt.figure(figsize=(6.75, 4.05))
    ax = fig.add_subplot(111)
    for idx, label in enumerate(series.labels[1:], start=1):
        if idx < len(series.data):
            ax.plot(series.data[0], series.data[idx], label=label.strip(), lw=1.0, alpha=0.65)
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.legend(loc=1, fontsize=9, frameon=False)
    fig.suptitle(title, fontsize=14)
    return PlotResult(fig, detect_plot_kind(path) or "series", title, metrics)


def chi_history_ln(chi: np.ndarray) -> np.ndarray:
    """``ln(max(chi, 1e-12))`` per log row; a non-finite chi^2 stays ``NaN`` (a gap).

    The same clamp as ``plotDataFromText`` in browserData.js, so the PNG, the
    JSON series and the browser chart plot identical numbers.
    """
    chi = np.asarray(chi, dtype=float)
    out = np.full(chi.shape, np.nan)
    finite = np.isfinite(chi)
    out[finite] = np.log(np.maximum(chi[finite], 1e-12))
    return out


CHI_HISTORY_Y_LABEL = "ln(χ²)"
UNNAMED_CHI_COLUMN = "last log column"


def chi_history_labels(column: str | None) -> tuple[str, str]:
    """``(title, series label)`` of the chi^2 history of one log column.

    The plotted series is the LAST column of the RMCProfile ``.log`` — the chi^2
    of one fit term (``X_ray_(R)1``: the X-ray real-space fit), not a total or an
    R-factor — so it is named by its header, never as "R-value". Mirrors
    ``chiHistoryLabels()`` in browserData.js.
    """
    name = column or UNNAMED_CHI_COLUMN
    return f"χ² history: {name}", name


def _chi_plot(path: Path) -> PlotResult:
    log_paths = related_r_value_logs(path)
    log = read_chi_log(log_paths)
    chi_r = log.chi_r
    if len(chi_r) == 0:
        raise ValueError(f"{path} does not contain chi values")
    title, label = chi_history_labels(log.column)

    fig = plt.figure(figsize=(6.75, 4.05))
    ax = fig.add_subplot(111)
    ax.plot(chi_history_ln(chi_r), label=label, lw=1.0, alpha=0.65)
    ax.set_xlabel("Time steps", fontsize=11)
    ax.set_ylabel(r"ln($\chi^2$)", fontsize=11)
    ax.legend(loc=1, fontsize=9, frameon=False)
    fig.suptitle(title, fontsize=14)
    return PlotResult(fig, "r_value", title, {"final_chi_r": float(chi_r[-1])})


def _stog_plot(path: Path) -> PlotResult:
    data = read_stog(path)
    title = path.name
    xlabel = r"r ($\mathrm{\AA}$)" if path.name.endswith(".gr") else r"Q ($\mathrm{\AA^{-1}}$)"
    ylabel = stog_function_label(path.name)

    fig = plt.figure(figsize=(6.75, 4.725))
    ax = fig.add_subplot(111)
    ax.plot(data[0], data[1], label=path.name, lw=1.0, alpha=1.0, color="r")
    ax.hlines(0 if path.name.endswith(".fq") else 1, data[0][0], data[0][-1], ls="--", lw=0.5, color="black")
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.legend(loc=1, fontsize=9, frameon=False)
    return PlotResult(fig, "stog", title)


def make_plot(path: str | Path) -> PlotResult:
    path = Path(path)
    kind = detect_plot_kind(path)
    if kind is None:
        raise ValueError(f"Unsupported plot file type: {path.name}")

    if kind == "r_value":
        return _chi_plot(path)
    if kind == "stog":
        return _stog_plot(path)

    reader = read_exafs_csv if kind in ("exafs_q", "exafs_r") else read_rmc_csv
    labels = reader(path).labels
    title, y_label = series_titles(kind, path.name, labels)
    if kind == "exafs_q":
        return _series_plot(path, title, r"k ($\mathrm{\AA^{-1}}$)", r"$\chi(k) k^2$", reader=read_exafs_csv)
    if kind == "exafs_r":
        return _series_plot(path, title, r"r ($\mathrm{\AA}$)", r"FT[$\chi(k) k^2$]", reader=read_exafs_csv)
    if kind == "xpdf":
        return _series_plot(path, title, r"r ($\mathrm{\AA}$)", y_label, calculate_rwp=True)
    if kind == "npdf":
        result = _series_plot(path, title, r"r ($\mathrm{\AA}$)", y_label, calculate_rwp=True)
        metrics = dict(result.metrics)
        metrics["pdf_index"] = float(pdf_index(path))
        return PlotResult(result.figure, kind, title, metrics)
    if kind == "pdf_partials":
        return _series_plot(path, title, r"r ($\mathrm{\AA}$)", y_label, calculate_rwp=False)
    if kind in ("xray_sq", "neutron_sq"):
        return _series_plot(path, title, r"Q ($\mathrm{\AA^{-1}}$)", y_label, calculate_rwp=True)
    if kind == "bragg":
        x_label = r"ToF ($\mu$s)" if bragg_is_tof(labels[0] if labels else None) else r"Q ($\mathrm{\AA^{-1}}$)"
        return _series_plot(path, title, x_label, y_label, calculate_rwp=True)
    raise ValueError(f"Unsupported plot file type: {path.name}")


def plot_to_png(result: PlotResult, dpi: int = 150) -> bytes:
    image = io.BytesIO()
    result.figure.savefig(image, format="png", bbox_inches="tight", dpi=dpi)
    plt.close(result.figure)
    return image.getvalue()


def close_plot(result: PlotResult) -> None:
    plt.close(result.figure)
