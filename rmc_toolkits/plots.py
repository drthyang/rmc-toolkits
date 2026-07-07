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

from .parsers import pdf_index, read_chi, read_exafs_csv, read_rmc_csv, read_stog, related_r_value_logs, rwp


@dataclass(frozen=True)
class PlotResult:
    figure: plt.Figure
    kind: str
    title: str
    metrics: dict[str, float] = field(default_factory=dict)


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
    if name.endswith("_FQ1.csv"):
        return "xray_sq"
    if name.endswith("_SQ1.csv"):
        return "neutron_sq"
    if re.search(r"_bragg(?:_.+)?\.csv$", name):
        return "bragg"
    if re.search(r"-\d{2,}\.log$", name):
        return "r_value"
    if name in {"scale_ft.gr", "scale_ft.sq", "scale_ft_rmc.fq"}:
        return "stog"
    return None


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

    metrics: dict[str, float] = {}
    if calculate_rwp and len(series.data) >= 3:
        metrics["rwp"] = rwp(series.data[0], series.data[1], series.data[2])

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


def _chi_plot(path: Path) -> PlotResult:
    log_paths = related_r_value_logs(path)
    _, chi_r = read_chi(log_paths)
    if len(chi_r) == 0:
        raise ValueError(f"{path} does not contain chi values")

    fig = plt.figure(figsize=(6.75, 4.05))
    ax = fig.add_subplot(111)
    ax.plot(np.log(chi_r), label="R", lw=1.0, alpha=0.65)
    ax.set_xlabel("Time steps", fontsize=11)
    ax.set_ylabel(r"log($\chi$)", fontsize=11)
    ax.legend(loc=1, fontsize=9, frameon=False)
    fig.suptitle("R-value", fontsize=14)
    return PlotResult(fig, "r_value", "R-value", {"final_chi_r": float(chi_r[-1])})


def _stog_plot(path: Path) -> PlotResult:
    data = read_stog(path)
    title = path.name
    xlabel = r"r ($\mathrm{\AA}$)" if path.name.endswith(".gr") else r"Q ($\mathrm{\AA^{-1}}$)"
    ylabel = "G(r)" if path.name.endswith(".gr") else "S(Q)"

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

    if kind == "exafs_q":
        return _series_plot(
            path,
            "EXAFS Q-space",
            r"k ($\mathrm{\AA^{-1}}$)",
            r"$\chi(k) k^2$",
            reader=read_exafs_csv,
        )
    if kind == "exafs_r":
        return _series_plot(
            path,
            "EXAFS R-space",
            r"r ($\mathrm{\AA}$)",
            r"FT[$\chi(k) k^2$]",
            reader=read_exafs_csv,
        )
    if kind == "xpdf":
        return _series_plot(path, "xPDF", r"r ($\mathrm{\AA}$)", calculate_rwp=True)
    if kind == "npdf":
        neutron_index = pdf_index(path)
        title = path.stem.split("_")[-1]
        result = _series_plot(path, title, r"r ($\mathrm{\AA}$)", calculate_rwp=True)
        metrics = dict(result.metrics)
        metrics["pdf_index"] = float(neutron_index)
        return PlotResult(result.figure, kind, title, metrics)
    if kind == "pdf_partials":
        return _series_plot(path, path.stem.split("_")[-1], r"r ($\mathrm{\AA}$)", calculate_rwp=False)
    if kind == "xray_sq":
        labels = read_rmc_csv(path).labels
        return _series_plot(path, "S(Q) (x-ray)", labels[0] if labels else r"Q ($\mathrm{\AA^{-1}}$)", calculate_rwp=True)
    if kind == "neutron_sq":
        labels = read_rmc_csv(path).labels
        return _series_plot(path, "S(Q) (neutron)", labels[0] if labels else r"Q ($\mathrm{\AA^{-1}}$)", calculate_rwp=True)
    if kind == "bragg":
        labels = read_rmc_csv(path).labels
        x_label = r"ToF ($\mu$s)" if bragg_is_tof(labels[0] if labels else None) else r"Q ($\mathrm{\AA^{-1}}$)"
        return _series_plot(path, "BRAGG", x_label, calculate_rwp=True)
    if kind == "r_value":
        return _chi_plot(path)
    return _stog_plot(path)


def plot_to_png(result: PlotResult, dpi: int = 150) -> bytes:
    image = io.BytesIO()
    result.figure.savefig(image, format="png", bbox_inches="tight", dpi=dpi)
    plt.close(result.figure)
    return image.getvalue()


def close_plot(result: PlotResult) -> None:
    plt.close(result.figure)
