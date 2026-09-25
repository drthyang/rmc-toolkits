# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
import io
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import threading

import numpy as np
from flask import Flask, jsonify, request, send_file, send_from_directory
from flask.json.provider import DefaultJSONProvider
from flask_cors import CORS

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
FRONTEND_DIST = PROJECT_ROOT / "web_app" / "frontend" / "dist"

from rmc_toolkits.kde import UnitCellPositions, load_unit_cell_positions, oriented_kde_slice
from rmc_toolkits.orientation import site_orientation_histogram
from rmc_toolkits.pca_kde import (
    SiteDisplacements,
    load_site_displacements,
    site_ellipsoids,
    site_pca_kde,
)
from rmc_toolkits.triplets import APP_MAX_ANGLES, cached_bond_angle_summary
from rmc_toolkits.parsers import (
    parse_rmc6f_atoms,
    read_cell_vectors,
    read_moves_metadata,
    read_chi_log,
    rmc6f_problem,
    read_dat_header,
    read_exafs_csv,
    read_rmc_csv,
    read_stog,
    read_stog_inp,
    read_stog_xy,
    related_r_value_logs,
    write_frac_from_rmc6f,
)
from rmc_toolkits.plots import (
    CHI_HISTORY_Y_LABEL,
    bragg_is_tof,
    chi_history_labels,
    chi_history_ln,
    close_plot,
    detect_plot_kind,
    make_plot,
    plot_to_png,
    series_titles,
    stog_function_label,
)
from rmc_toolkits.scaling import (
    ScalingConfig,
    autoscale,
    diagnostics_summary,
    scale_pipeline,
)
from rmc_toolkits.scaling_cli import (  # shared writer keeps CLI/API outputs identical
    CliError,
    _json_safe,
    _write_outputs,
    _resolve_targets as _resolve_scaling_targets,
    refuse_failed_fit,
    resolve_coefficients,
    stog_inp_closest_approach,
    usable_sigma,
)
from rmc_toolkits.scaling import auto_enforcement_cutoff, detect_first_peak_onset
from rmc_toolkits.scattering import faber_ziman, number_density_from_mass_density
from rmc_toolkits.transforms import first_peak_zero, g_to_gk, gk_to_dr


def _finite_json(value):
    """``value`` with every non-finite float (NaN, +/-Inf) replaced by ``None``.

    Walks dicts, lists/tuples, NumPy arrays and NumPy scalars, so a masked region
    (NaN in an RMCProfile CSV or log) reaches the browser as JSON ``null`` — a gap
    the chart skips — instead of the bare token ``NaN`` that ``JSON.parse`` rejects.
    """
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {key: _finite_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_finite_json(item) for item in value]
    if isinstance(value, np.ndarray):
        return _finite_json(value.tolist())
    if isinstance(value, np.generic):
        return _finite_json(value.item())
    return value


class StrictJSONProvider(DefaultJSONProvider):
    """App-wide JSON provider that only ever emits strict (RFC 8259) JSON.

    Non-finite floats become ``null``; ``allow_nan=False`` is the guard that makes
    any path that slips past the sanitizer fail loudly instead of emitting ``NaN``.
    The common all-finite payload is serialized in one pass; only a payload that
    trips the guard is walked by :func:`_finite_json` and serialized again.
    """

    @staticmethod
    def default(o):
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, np.generic):
            return o.item()
        return DefaultJSONProvider.default(o)

    def dumps(self, obj, **kwargs):
        kwargs.setdefault("default", self.default)
        kwargs.setdefault("ensure_ascii", self.ensure_ascii)
        kwargs.setdefault("sort_keys", self.sort_keys)
        kwargs["allow_nan"] = False
        try:
            return json.dumps(obj, **kwargs)
        except ValueError:
            return json.dumps(_finite_json(obj), **kwargs)


app = Flask(__name__, static_folder=str(FRONTEND_DIST), static_url_path="")
app.json = StrictJSONProvider(app)
CORS(app)

DATA_ROOT = Path(os.environ.get("RMC_TOOLKITS_DATA_ROOT", PROJECT_ROOT)).expanduser().resolve()
SELECTED_DATA_ROOTS: set[Path] = set()
SUPPORTED_PATTERNS = (
    "*.csv", "*.log", "*.rmc6f", "Frac*.txt", "scale_ft.*", "stog_input.dat",
    "*.inp", "*.sq", "*.dat",
)
MAX_STRUCTURE_POINTS = 1_000_000


def _is_inside_root(candidate: Path, root: Path) -> bool:
    return candidate == root or root in candidate.parents


def _resolve_inside_root(raw_path: str | None) -> Path:
    candidate = Path(raw_path or ".").expanduser()
    if not candidate.is_absolute():
        candidate = DATA_ROOT / candidate
    resolved = candidate.resolve()
    allowed_roots = (DATA_ROOT, *SELECTED_DATA_ROOTS)
    if not any(_is_inside_root(resolved, root) for root in allowed_roots):
        allowed = ", ".join(str(root) for root in allowed_roots)
        raise PermissionError(f"Path is outside configured data roots: {allowed}")
    return resolved


def _choose_folder(initial_dir: Path) -> Path | None:
    if platform.system() == "Darwin":
        osascript = shutil.which("osascript") or "/usr/bin/osascript"
        if not Path(osascript).exists():
            raise FileNotFoundError("macOS folder picker is unavailable: osascript was not found")

        escaped_initial_dir = str(initial_dir).replace('"', '\\"')
        script = (
            'POSIX path of (choose folder with prompt "Select an RMCProfile run folder" '
            f'default location POSIX file "{escaped_initial_dir}")'
        )
        result = subprocess.run(
            [osascript, "-e", script],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            return None
        return Path(result.stdout.strip()).expanduser().resolve()

    import tkinter as tk
    from tkinter import filedialog

    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    try:
        selected = filedialog.askdirectory(
            initialdir=str(initial_dir),
            title="Select an RMCProfile run folder",
            mustexist=True,
        )
    finally:
        root.destroy()

    return Path(selected).expanduser().resolve() if selected else None


def _nearest_existing_directory(path: Path) -> Path:
    current = path if path.is_dir() else path.parent
    while not current.exists() or not current.is_dir():
        if current == current.parent:
            return DATA_ROOT
        current = current.parent
    return current


def _file_payload(path: Path, kind: str = "file") -> dict[str, object]:
    payload = {
        "name": path.name,
        "path": str(path),
        "type": kind,
        "plotKind": detect_plot_kind(path) if kind == "file" else None,
    }
    try:
        stat = path.stat()
        payload["modified"] = stat.st_mtime
        payload["size"] = stat.st_size if kind == "file" else None
    except OSError:
        payload["modified"] = None
        payload["size"] = None
    return payload


def _run_stem_from_output_name(name: str) -> tuple[int, str] | None:
    patterns = (
        (0, r"^(.+)-\d{2,}\.log$"),
        (1, r"^(.+)-EXAFS-.+_[QR]_OUTPUT\.csv$"),
        (1, r"^(.+)_FT_XFQ\d+\.csv$"),
        (1, r"^(.+)_[FS]Q\d+\.csv$"),
        (1, r"^(.+)_bragg(?:_.+)?\.csv$"),
        (1, r"^(.+)_PDF(?:partials|\d+)?\.csv$"),
        (2, r"^Frac_coord_(.+)\.txt$"),
    )
    for priority, pattern in patterns:
        match = re.match(pattern, name)
        if match:
            return priority, match.group(1)
    return None


def _find_rmc6f(directory: Path) -> Path:
    if directory.is_file() and directory.suffix == ".rmc6f":
        return directory
    all_rmc6f = sorted(directory.glob("*.rmc6f"))
    if not all_rmc6f:
        raise FileNotFoundError(f"No .rmc6f file found in {directory}")
    # A 0-byte or marker-less configuration (e.g. left by a killed run) must not
    # hide the valid ones beside it: skip it and fall through to the next match.
    # With no usable candidate at all the stem rule still names one, and the
    # caller's reader reports what is wrong with it (see _require_usable_rmc6f).
    rmc6f_files = [path for path in all_rmc6f if rmc6f_problem(path) is None] or all_rmc6f
    rmc6f_by_stem = {path.stem: path for path in rmc6f_files}
    output_stems: list[tuple[int, str, str]] = []
    for item in sorted(directory.iterdir(), key=lambda path: path.name.lower()):
        if not item.is_file():
            continue
        stem_match = _run_stem_from_output_name(item.name)
        if stem_match:
            priority, stem = stem_match
            output_stems.append((priority, item.name.lower(), stem))
    for _, _, stem in sorted(output_stems):
        if stem in rmc6f_by_stem:
            return rmc6f_by_stem[stem]
    return rmc6f_files[0]


def _require_usable_rmc6f(rmc6f_path: Path, target: Path) -> None:
    """Raise a FileNotFoundError that lists every candidate when none is usable."""
    if rmc6f_problem(rmc6f_path) is None:
        return
    directory = target if target.is_dir() else rmc6f_path.parent
    candidates = sorted(directory.glob("*.rmc6f")) if target.is_dir() else [rmc6f_path]
    listing = ", ".join(f"{path.name} ({rmc6f_problem(path)})" for path in candidates)
    raise FileNotFoundError(f"No usable .rmc6f file in {directory}: {listing}")


def _sample_atoms_by_site(atoms: list[dict], max_points: int) -> tuple[list[dict], int]:
    if len(atoms) <= max_points:
        return atoms, 1

    # Legacy coords-only atoms carry no reference number; they form one group
    # (key None), sorted after the numbered sites.
    by_reference: dict[int | None, list[dict]] = {}
    for atom in atoms:
        by_reference.setdefault(atom["reference_number"], []).append(atom)

    quota = max(1, max_points // len(by_reference))
    sampled: list[dict] = []
    for reference_number in sorted(by_reference, key=lambda ref: (ref is None, ref or 0)):
        group = by_reference[reference_number]
        stride = max(1, len(group) // quota)
        sampled.extend(group[::stride][:quota])

    return sampled[:max_points], max(1, len(atoms) // max_points)


def _clean_axis_label(label: str) -> str:
    normalized = label.strip()
    if normalized == "Q":
        return "Q (Å^{-1})"
    if normalized in ("r", "R"):
        return "r (Å)"
    return (
        normalized.replace("(A^-1)", "(Å^{-1})")
        .replace("(A^{-1})", "(Å^{-1})")
        .replace("(A)", "(Å)")
    )


# --- Request parameter parsing ------------------------------------------------
# Every numeric query-string or JSON parameter goes through _number(). It
# rejects anything that is not a finite number (text, lists, objects, booleans,
# NaN, +-inf), non-integral values for integer parameters, and values outside
# the parameter's documented range, by raising ValueError naming the parameter.
# Every route maps ValueError to HTTP 400, so a bad value never reaches an
# engine that would answer 200 with bare NaN/Infinity tokens (invalid JSON for
# a browser) or with a silently empty map. Grid sizes are clamped (not
# rejected) to the same limits the engines apply. The accepted ranges are
# listed with the endpoints in docs/REFERENCE.md.

KDE_GRID_CLAMP = (16, 400)  # kde.kde_slice clamps to the same range
KDE_MAX_LEVELS = 64
PCA_GRID_CLAMP = (8, 128)  # pca_kde.pca_kde_volume clamps to the same range
MAX_ORIENTATION_SMOOTHING = 64  # neighbour-diffusion passes (the UI offers 0-12)


def _number(
    raw,
    name: str,
    *,
    integer: bool = False,
    gt: float | None = None,
    ge: float | None = None,
    lt: float | None = None,
    le: float | None = None,
    clamp: tuple[int, int] | None = None,
) -> float | int:
    """Parse one request value as a finite number and check its range."""
    if isinstance(raw, bool) or not isinstance(raw, (int, float, str)):
        raise ValueError(f"{name} must be a number, got {raw!r}")
    try:
        value = float(raw)
    except (ValueError, OverflowError):
        raise ValueError(f"{name} must be a number, got {raw!r}") from None
    if not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number, got {raw!r}")
    if integer:
        if not value.is_integer():
            raise ValueError(f"{name} must be an integer, got {raw!r}")
        value = int(value)
    if gt is not None and not value > gt:
        raise ValueError(f"{name} must be > {gt:g}, got {value:g}")
    if ge is not None and not value >= ge:
        raise ValueError(f"{name} must be >= {ge:g}, got {value:g}")
    if lt is not None and not value < lt:
        raise ValueError(f"{name} must be < {lt:g}, got {value:g}")
    if le is not None and not value <= le:
        raise ValueError(f"{name} must be <= {le:g}, got {value:g}")
    if clamp is not None:
        value = min(max(value, clamp[0]), clamp[1])
    return value


def _query_number(name: str, default, **rules):
    """A numeric query-string argument; missing or blank means ``default``."""
    raw = request.args.get(name)
    if raw is None or not raw.strip():
        return default
    return _number(raw, name, **rules)


def _strict_result_response(payload):
    """JSON response for a computed numeric result; NaN/Infinity is a ValueError (400).

    Finite but extreme parameters (a bandwidth of 1e-200, an extent of 1e300)
    pass the range checks yet make the float64 kernels underflow or overflow.
    Such a result is an error, not data: never a 200 with bare NaN tokens
    (invalid JSON for a browser) nor with nulls (a silently empty map).

    The check calls the stdlib encoder with ``allow_nan=False`` itself instead
    of going through ``app.json``, so it holds whatever JSON provider the app
    installs -- including one that writes non-finite floats as ``null`` (the
    right policy for a masked *data series*, the wrong one for a computed
    density). The provider still supplies ``default`` for non-JSON types.
    """
    provider = app.json
    try:
        body = json.dumps(
            payload,
            allow_nan=False,
            default=getattr(provider, "default", None),
            ensure_ascii=getattr(provider, "ensure_ascii", True),
            sort_keys=getattr(provider, "sort_keys", True),
        )
    except ValueError:
        raise ValueError(
            "the result contains NaN or Infinity for these parameters (an extreme "
            "bandwidth, scale or extent?); use less extreme values"
        ) from None
    return app.response_class(f"{body}\n", mimetype=getattr(provider, "mimetype", "application/json"))


# --- Parsed-file caches ----------------------------------------------------------
# Parsing a 50k-atom .rmc6f takes about a second, so the analysis routes keep
# small LRU caches of parsed files. Every cache is keyed on the file signature
# below, never on st_mtime alone: sshfs/SFTP mounts, `scp -p` and rsync from a
# coarse filesystem report whole-second mtimes, so a half-written file and the
# finished one can share an mtime -- and a parse of the half-written file would
# then be served until the server restarted.


def _file_signature(path: str | Path) -> tuple[int, int, int, int]:
    """Freshness key of a file: ``(st_mtime_ns, st_ctime_ns, st_size, st_ino)``.

    The size catches a completed write, the inode an atomic replace (write a
    temporary file, then rename), and the ctime -- set by the kernel on every
    write or utime, never by the writer -- a same-size rewrite whose mtime was
    reset (``scp -p``, ``rsync -t``) on a filesystem with sub-second ctimes.
    """
    stat = os.stat(path)
    return (stat.st_mtime_ns, stat.st_ctime_ns, stat.st_size, stat.st_ino)


class SourceChangedError(RuntimeError):
    """The source file kept changing while it was being read (HTTP 409)."""


class _FileCache:
    """Thread-safe LRU of values parsed from one file, keyed on its signature.

    ``get`` re-takes the signature after loading and stores the value only if
    the file did not change during the read. A file that changed is read once
    more under its new signature; if it changes again (a writer still busy),
    ``SourceChangedError`` is raised -- a torn read is never cached or served.
    Storing a fresh entry drops the entries of older signatures of that path.
    """

    def __init__(self, maxsize: int):
        self.maxsize = maxsize
        self._entries: OrderedDict = OrderedDict()
        self._lock = threading.Lock()

    def get(self, path: str | Path, params: tuple, load):
        path = str(path)
        for _attempt in range(2):
            signature = _file_signature(path)
            key = (path, signature, params)
            with self._lock:
                if key in self._entries:
                    self._entries.move_to_end(key)
                    return self._entries[key]
            value = load()
            if _file_signature(path) != signature:
                continue
            with self._lock:
                for stale in [k for k in self._entries if k[0] == path and k[1] != signature]:
                    del self._entries[stale]
                self._entries[key] = value
                while len(self._entries) > self.maxsize:
                    self._entries.popitem(last=False)
            return value
        raise SourceChangedError(
            f"{Path(path).name} changed while it was being read (it is probably still "
            "being written); retry in a moment"
        )

    def clear(self) -> None:
        with self._lock:
            self._entries.clear()


_POSITIONS_CACHE = _FileCache(16)  # /api/kde/slice: per (file, element)
_SITES_CACHE = _FileCache(8)  # /api/pca/sites, /api/pca/kde, /api/pca/orientation
_TRIPLETS_CACHE = _FileCache(16)  # /api/triplets: per (file, every parameter)
_SCALING_CACHE = _FileCache(8)  # /api/scaling/*: per (data file, config, mode, a, b, sigma)


@app.route("/api/health", methods=["GET"])
def health():
    return jsonify({"status": "ok", "dataRoot": str(DATA_ROOT)})


@app.route("/api/files", methods=["GET"])
def list_files():
    try:
        directory = _resolve_inside_root(request.args.get("dir", "."))
        if not directory.exists() or not directory.is_dir():
            return jsonify({"error": "Directory not found"}), 404

        paths: dict[Path, dict[str, object]] = {}
        for item in sorted(directory.iterdir(), key=lambda path: path.name.lower()):
            if item.is_dir() and not item.name.startswith("."):
                paths[item] = _file_payload(item, "directory")

        for pattern in SUPPORTED_PATTERNS:
            for item in directory.glob(pattern):
                if item.is_file():
                    paths[item] = _file_payload(item)

        files = sorted(paths.values(), key=lambda item: (item["type"] != "directory", item["name"].lower()))
        return jsonify({"root": str(DATA_ROOT), "currentPath": str(directory), "files": files})
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/dialog/folder", methods=["POST"])
def choose_folder():
    try:
        payload = request.get_json(silent=True) or {}
        initial_dir = _resolve_inside_root(payload.get("dir", "."))
        if initial_dir.is_file():
            initial_dir = initial_dir.parent
        initial_dir = _nearest_existing_directory(initial_dir)

        selected = _choose_folder(initial_dir)
        if selected is None:
            return jsonify({"error": "Folder selection cancelled"}), 400
        if not selected.exists() or not selected.is_dir():
            return jsonify({"error": "Selected path is not a folder"}), 400

        SELECTED_DATA_ROOTS.add(selected)
        return jsonify({"path": str(selected), "name": selected.name})
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/plot", methods=["GET"])
def plot_file():
    try:
        path = _resolve_inside_root(request.args.get("path"))
        if not path.exists() or not path.is_file():
            return jsonify({"error": "File not found"}), 404

        result = make_plot(path)
        image = io.BytesIO(plot_to_png(result))
        return send_file(
            image,
            mimetype="image/png",
            download_name=f"{path.stem}.png",
        )
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/plot/metadata", methods=["GET"])
def plot_metadata():
    try:
        path = _resolve_inside_root(request.args.get("path"))
        result = make_plot(path)
        metadata = {"kind": result.kind, "title": result.title, "metrics": result.metrics}
        close_plot(result)
        return jsonify(metadata)
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/plot/data", methods=["GET"])
def plot_data():
    try:
        path = _resolve_inside_root(request.args.get("path"))
        if not path.exists() or not path.is_file():
            return jsonify({"error": "File not found"}), 404

        kind = detect_plot_kind(path)
        if kind is None:
            return jsonify({"error": f"Unsupported plot file type: {path.name}"}), 400

        metadata_result = make_plot(path)
        metadata = {"kind": metadata_result.kind, "title": metadata_result.title, "metrics": metadata_result.metrics}
        close_plot(metadata_result)

        if kind == "r_value":
            log = read_chi_log(related_r_value_logs(path))
            chi_r = log.chi_r
            _, series_label = chi_history_labels(log.column)
            return jsonify(
                {
                    **metadata,
                    "xLabel": "Time steps",
                    "yLabel": CHI_HISTORY_Y_LABEL,
                    "chiColumn": log.column,
                    "series": [
                        {
                            "label": series_label,
                            "x": list(range(len(chi_r))),
                            # Non-finite chi^2 rows stay in the series (null in JSON).
                            "y": chi_history_ln(chi_r).tolist(),
                        }
                    ],
                }
            )

        if kind == "stog":
            data = read_stog(path)
            return jsonify(
                {
                    **metadata,
                    "xLabel": "r (Å)" if path.name.endswith(".gr") else "Q (Å^{-1})",
                    "yLabel": stog_function_label(path.name),
                    "series": [{"label": path.name, "x": data[0].tolist(), "y": data[1].tolist()}],
                }
            )

        series = read_exafs_csv(path) if kind in ("exafs_q", "exafs_r") else read_rmc_csv(path)
        x_values = series.data[0].tolist()
        payload_series = []
        for idx, label in enumerate(series.labels[1:], start=1):
            if idx < len(series.data):
                payload_series.append({"label": label.strip() or f"Series {idx}", "x": x_values, "y": series.data[idx].tolist()})

        # One label source for Flask, the PNGs and (mirrored) the browser:
        # F(Q) for *_FQn.csv, partial g(r) for PDFpartials, from the file's headers.
        _, y_label = series_titles(kind, path.name, series.labels)
        x_label = series.labels[0] if series.labels else "x"
        if kind == "exafs_q":
            x_label = "k (Å^{-1})"
        elif kind in ("exafs_r", "xpdf", "npdf", "pdf_partials"):
            x_label = "r (Å)"
        elif kind in ("xray_sq", "neutron_sq"):
            x_label = "Q (Å^{-1})"
        elif kind == "bragg":
            x_label = "ToF (µs)" if bragg_is_tof(series.labels[0] if series.labels else None) else "Q (Å^{-1})"
        else:
            x_label = _clean_axis_label(x_label)

        return jsonify({**metadata, "xLabel": x_label, "yLabel": y_label, "series": payload_series})
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/convert/frac", methods=["POST"])
def convert_frac():
    try:
        payload = request.get_json(silent=True) or {}
        source = _resolve_inside_root(payload.get("path"))
        if source.suffix != ".rmc6f":
            return jsonify({"error": "Expected a .rmc6f file"}), 400

        output_raw = payload.get("outputPath")
        output = _resolve_inside_root(output_raw) if output_raw else None
        out_path = write_frac_from_rmc6f(
            source,
            output_path=output,
            overwrite=bool(payload.get("overwrite", False)),
        )
        return jsonify({"path": str(out_path), "name": out_path.name})
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileExistsError as exc:
        return jsonify({"error": str(exc)}), 409
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/structure", methods=["GET"])
def structure():
    try:
        target = _resolve_inside_root(request.args.get("dir", "."))
        max_points = _query_number(
            "maxPoints", MAX_STRUCTURE_POINTS, integer=True, clamp=(100, MAX_STRUCTURE_POINTS)
        )
        rmc6f_path = _find_rmc6f(target)
        _require_usable_rmc6f(rmc6f_path, target)
        lattice_vectors, supercell = read_cell_vectors(rmc6f_path)
        moves = read_moves_metadata(rmc6f_path)

        # Same line grammar as the browser parser (rmc6f.js): full-layout and
        # legacy coords-only atoms both count; non-finite / unparsed lines are
        # reported, and zero parsed atoms is an error naming what was found.
        atoms, parse_report = parse_rmc6f_atoms(rmc6f_path, include_coords_only=True)
        index_sets: dict[str, set[int]] = {}
        for atom in atoms:
            if atom["reference_number"] is not None:
                index_sets.setdefault(atom["element"], set()).add(int(atom["reference_number"]))
        atom_indices = {element: sorted(indices) for element, indices in index_sets.items()}
        sampled, stride = _sample_atoms_by_site(atoms, max_points)
        points = []
        counts: dict[str, int] = {}
        for atom in atoms:
            counts[atom["element"]] = counts.get(atom["element"], 0) + 1
        for atom in sampled:
            # Fold the box coordinate into one unit cell; subtracting the cell
            # index first only removes an integer, so coords-only atoms fold too.
            unit_cell = (atom["coords"] * supercell) % 1.0
            points.append(
                {
                    "element": atom["element"],
                    "referenceNumber": atom["reference_number"],
                    "boxX": float(atom["coords"][0]),
                    "boxY": float(atom["coords"][1]),
                    "boxZ": float(atom["coords"][2]),
                    "x": float(unit_cell[0]),
                    "y": float(unit_cell[1]),
                    "z": float(unit_cell[2]),
                }
            )

        return jsonify(
            {
                "source": str(rmc6f_path),
                "totalAtoms": len(atoms),
                "sampledAtoms": len(points),
                "sampleStride": stride,
                "elements": sorted(counts.keys()),
                "elementCounts": counts,
                "atomIndices": atom_indices,
                "supercell": supercell.tolist(),
                "latticeVectors": lattice_vectors.tolist(),
                "moves": moves,
                "parseReport": parse_report.to_dict(),
                "parseWarning": parse_report.warning(),
                "points": points,
            }
        )
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


def _cached_positions(rmc6f_path: Path, element: str | None) -> UnitCellPositions:
    return _POSITIONS_CACHE.get(
        rmc6f_path, (element,), lambda: load_unit_cell_positions(str(rmc6f_path), element=element)
    )


SLICE_ORIENTATIONS = {
    "a": {"normal": (1.0, 0.0, 0.0), "u": (0.0, 1.0, 0.0), "v": (0.0, 0.0, 1.0)},
    "b": {"normal": (0.0, 1.0, 0.0), "u": (1.0, 0.0, 0.0), "v": (0.0, 0.0, 1.0)},
    "c": {"normal": (0.0, 0.0, 1.0), "u": (1.0, 0.0, 0.0), "v": (0.0, 1.0, 0.0)},
}


def _bw_argument(raw: str | None, default: str = "scott") -> str | float:
    """Bandwidth query arg: the names 'scott'/'silverman', else a positive finite number."""
    if raw is None or not str(raw).strip():
        return default
    name = str(raw).strip().lower()
    if name in ("scott", "silverman"):
        return name
    try:
        return _number(raw, "bw", gt=0.0)
    except ValueError:
        raise ValueError(
            f"bw must be 'scott', 'silverman' or a positive finite number, got {raw!r}"
        ) from None


def _slice_orientation_from_request():
    orientation = (request.args.get("orientation") or "c").lower()
    if orientation in SLICE_ORIENTATIONS:
        config = SLICE_ORIENTATIONS[orientation]
        return orientation, config["normal"], config["u"], config["v"]

    # A zero normal is rejected by kde._plane_basis (ValueError -> 400).
    normal = (
        _query_number("nx", 0.0),
        _query_number("ny", 0.0),
        _query_number("nz", 1.0),
    )
    return "custom", normal, None, None


@app.route("/api/kde/slice", methods=["GET"])
def kde_slice_endpoint():
    try:
        target = _resolve_inside_root(request.args.get("dir", "."))
        rmc6f_path = _find_rmc6f(target)

        element = request.args.get("element") or None
        if element in ("", "all"):
            element = None

        orientation, normal, u_axis, v_axis = _slice_orientation_from_request()

        # z and dz arrive as fractions of the projection range along the slice
        # normal. Keep the KDE slice in fractional coordinates so non-orthogonal
        # cells can be projected through the actual cell basis in the frontend.
        # z is clamped to [0, 1] by the engine (echoed back as "center").
        z_frac = _query_number("z", 0.5)
        dz_frac = _query_number("dz", 0.08, gt=0.0, le=1.0)
        bw = _query_number("bw", 0.03, gt=0.0)
        grid = _query_number("grid", 120, integer=True, clamp=KDE_GRID_CLAMP)
        levels = _query_number("levels", 8, integer=True, ge=0, le=KDE_MAX_LEVELS)
        log = request.args.get("log", "false").lower() in ("1", "true", "yes")

        positions = _cached_positions(rmc6f_path, element)
        cell_lengths = positions.cell_lengths

        result = oriented_kde_slice(
            positions.fractional_positions,
            center=z_frac,
            thickness=dz_frac,
            normal=normal,
            u_axis=u_axis,
            v_axis=v_axis,
            bw=bw,
            grid=grid,
            log=log,
            n_levels=levels,
        )
        result["cellLengths"] = cell_lengths.tolist()
        result["unitVectors"] = positions.unit_vectors.tolist()
        result["orientation"] = orientation
        result["source"] = str(rmc6f_path)
        result["element"] = element or "all"
        return _strict_result_response(result)
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except ValueError as exc:  # includes numpy.linalg.LinAlgError
        return jsonify({"error": str(exc)}), 400
    except SourceChangedError as exc:
        return jsonify({"error": str(exc)}), 409
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


def _cached_site_displacements(rmc6f_path: Path) -> SiteDisplacements:
    # Looked up at call time (not bound here) so tests can intercept the loader.
    return _SITES_CACHE.get(rmc6f_path, (), lambda: load_site_displacements(str(rmc6f_path)))


@app.route("/api/pca/sites", methods=["GET"])
def pca_sites_endpoint():
    """Anisotropic displacement tensor + thermal ellipsoid for every site.

    One cheap batched pass over the whole configuration; the frontend uses it to
    populate the site picker and the per-site ellipsoid summary table.
    """
    try:
        target = _resolve_inside_root(request.args.get("dir", "."))
        rmc6f_path = _find_rmc6f(target)
        probability = _query_number("probability", 0.5, gt=0.0, lt=1.0)
        sites = _cached_site_displacements(rmc6f_path)
        ellipsoids = site_ellipsoids(sites, probability=probability)
        return jsonify(
            {
                "source": str(rmc6f_path),
                "referenceNumbers": sites.reference_numbers.tolist(),
                # Every species present, including the minority species of a
                # mixed-occupancy site (site labels carry only the majority).
                "elements": sites.species,
                "totalAtoms": int(sites.counts.sum()),
                "latticeVectors": sites.lattice_vectors.tolist(),
                "supercell": sites.supercell.tolist(),
                "probability": probability,
                "sites": ellipsoids,
                # Atom lines the shared .rmc6f grammar skipped (non-finite
                # coordinates, unparsed lines, a header count mismatch), or
                # null for a clean file -- a dropped atom is never silent.
                "parseWarning": sites.parse_warning,
            }
        )
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except ValueError as exc:
        # Bad input (no parseable atom, probability outside (0, 1)): a clear
        # 400 naming the problem, not a 500.
        return jsonify({"error": str(exc)}), 400
    except SourceChangedError as exc:
        return jsonify({"error": str(exc)}), 409
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/pca/kde", methods=["GET"])
def pca_kde_endpoint():
    """PCA-aligned 3D KDE volume for one site (or one element's pooled sites)."""
    try:
        target = _resolve_inside_root(request.args.get("dir", "."))
        rmc6f_path = _find_rmc6f(target)

        reference_number = _query_number("referenceNumber", None, integer=True)
        element = request.args.get("element") or None
        if element in ("", "all"):
            element = None
        bw = _bw_argument(request.args.get("bw"))
        bw_scale = _query_number("bwScale", 1.0, gt=0.0)
        grid = _query_number("grid", 48, integer=True, clamp=PCA_GRID_CLAMP)
        extent = _query_number("extent", 3.0, gt=0.0)
        probability = _query_number("probability", 0.5, gt=0.0, lt=1.0)

        sites = _cached_site_displacements(rmc6f_path)
        result = site_pca_kde(
            sites,
            reference_number=reference_number,
            element=element,
            bw=bw,
            bw_scale=bw_scale,
            grid=grid,
            extent=extent,
            cubic_box=request.args.get("cubicBox", "false").lower() in ("1", "true", "yes"),
            probability=probability,
            projections=request.args.get("projections", "true").lower() in ("1", "true", "yes"),
        )
        result["source"] = str(rmc6f_path)
        return _strict_result_response(result)
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except SourceChangedError as exc:
        return jsonify({"error": str(exc)}), 409
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


# Work budget for one /api/triplets request (same constant as the browser
# worker's): an over-budget spec is a 400 before any angle is formed.
TRIPLETS_MAX_ANGLES = APP_MAX_ANGLES


@app.route("/api/triplets", methods=["GET"])
def triplets_endpoint():
    """Bond-angle (triplet) summary of the run's configuration.

    Params: end1/apex/end2 (elements, apex central), r12Min/r12Max and the
    optional r23Min/r23Max windows (angstrom, inclusive), binWidth (degrees).
    Payload shape is defined by rmc_toolkits.triplets.bond_angle_summary and
    mirrored by the browser worker's 'triplets' request. Both boundaries apply
    the same caps (the engine itself is unrestricted for library/CLI use):
    rmax <= 15 A bounds the neighbour search, binWidth >= 0.05 deg the
    response size, and TRIPLETS_MAX_ANGLES -- the engine's APP_MAX_ANGLES,
    checked against the exact angle count before any angle is formed -- the
    pairing work, which grows ~rmax^6 and which the rmax cap does not bound.
    """
    try:
        target = _resolve_inside_root(request.args.get("dir", "."))
        rmc6f_path = _find_rmc6f(target)

        def window(min_key, max_key, fallback=None):
            raw_min, raw_max = request.args.get(min_key), request.args.get(max_key)
            # Blank (empty or whitespace) is missing, exactly as in the browser
            # worker -- never a 0 A bound.
            missing = [
                key
                for key, raw in ((min_key, raw_min), (max_key, raw_max))
                if raw is None or not raw.strip()
            ]
            if len(missing) == 2 and fallback is not None:
                return fallback
            if missing:
                raise ValueError(f"{min_key}/{max_key} are required together; missing {missing[0]}")
            bounds = _number(raw_min, min_key), _number(raw_max, max_key)
            if bounds[1] > 15.0:
                raise ValueError(f"{max_key} is capped at 15 A for API requests, got {bounds[1]}")
            return bounds

        window12 = window("r12Min", "r12Max")
        window23 = window("r23Min", "r23Max", fallback=window12)
        bin_width = _query_number("binWidth", 1.0)
        if bin_width < 0.05:
            raise ValueError(f"binWidth is capped at >= 0.05 deg for API requests, got {bin_width}")
        params = (
            # Normalized here so 'se' and 'Se' share one cache entry.
            request.args.get("end1", "").strip().capitalize(),
            request.args.get("apex", "").strip().capitalize(),
            request.args.get("end2", "").strip().capitalize(),
            *window12,
            *window23,
            bin_width,
            TRIPLETS_MAX_ANGLES,
        )
        # The library function's own lru_cache is keyed on the caller's mtime;
        # call the uncached body (__wrapped__, which ignores that key) under
        # the file-signature cache instead. `params` is both that cache key and
        # the engine's argument list, so every engine argument -- a work budget
        # included -- belongs in it (TripletsWorkBudgetTests guards this).
        result = dict(
            _TRIPLETS_CACHE.get(
                rmc6f_path,
                params,
                lambda: cached_bond_angle_summary.__wrapped__(str(rmc6f_path), None, *params),
            )
        )
        result["source"] = str(rmc6f_path)
        return jsonify(result)
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except SourceChangedError as exc:
        return jsonify({"error": str(exc)}), 409
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/pca/orientation", methods=["GET"])
def pca_orientation_endpoint():
    """Hex-binned solid-angle histogram of one site's displacement directions.

    Same site selection as /api/pca/kde (referenceNumber, or element for pooled
    sites); the binning and output are documented in rmc_toolkits.orientation.
    """
    try:
        target = _resolve_inside_root(request.args.get("dir", "."))
        rmc6f_path = _find_rmc6f(target)

        reference_number = _query_number("referenceNumber", None, integer=True)
        element = request.args.get("element") or None
        if element in ("", "all"):
            element = None
        # frequency's [1, 64] range is enforced by orientation.goldberg_tiling.
        frequency = _query_number("frequency", None, integer=True)
        min_amplitude = _query_number("minAmplitude", 0.0, ge=0.0)
        min_amplitude_quantile = _query_number("minAmplitudeQuantile", 0.0, ge=0.0, lt=1.0)
        smoothing = _query_number(
            "smoothing", 0, integer=True, ge=0, le=MAX_ORIENTATION_SMOOTHING
        )

        sites = _cached_site_displacements(rmc6f_path)
        result = site_orientation_histogram(
            sites,
            reference_number=reference_number,
            element=element,
            frequency=frequency,
            weight=request.args.get("weight", "count"),
            min_amplitude=min_amplitude,
            min_amplitude_quantile=min_amplitude_quantile,
            smoothing=smoothing,
            frame=request.args.get("frame", "cartesian"),
            geometry=request.args.get("geometry", "true").lower() in ("1", "true", "yes"),
        )
        result["source"] = str(rmc6f_path)
        result["parseWarning"] = sites.parse_warning
        return _strict_result_response(result)
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except SourceChangedError as exc:
        return jsonify({"error": str(exc)}), 409
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


# --- Auto StoG scaling API ---------------------------------------------------
# Thin HTTP face over rmc_toolkits.scaling; file writing is shared with the
# rmc-autoscale CLI so API outputs are identical to a CLI/classic-stog session.


def _payload_float(payload: dict, key: str, **rules) -> float | None:
    """A numeric JSON-body field (see ``_number``); missing, null or blank is None."""
    value = payload.get(key)
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    return _number(value, key, **rules)


def _payload_bool(payload: dict, key: str, default: bool) -> bool:
    value = payload.get(key)
    if value is None:
        return default
    if isinstance(value, str):
        return value.lower() in ("1", "true", "yes")
    return bool(value)


def _resolve_scaling_source(payload: dict):
    """Resolve the request source into (inp, inp_path, data_path, header)."""
    source = _resolve_inside_root(payload.get("path"))
    if not source.exists() or not source.is_file():
        raise FileNotFoundError(f"Source file not found: {source}")
    kind = (payload.get("kind") or "auto").lower()
    inp = None
    inp_path = None
    looks_like_inp = source.suffix == ".inp" or "input" in source.name.lower()
    if kind == "inp" or (kind == "auto" and looks_like_inp):
        try:
            inp = read_stog_inp(source)
            inp_path = source
        except (ValueError, NotImplementedError):
            if kind == "inp":
                raise
    data_path = (source.parent / inp.data_file).resolve() if inp is not None else source
    if inp is not None and not any(
        _is_inside_root(data_path, root) for root in (DATA_ROOT, *SELECTED_DATA_ROOTS)
    ):
        raise PermissionError("stog input's data file lies outside the configured data roots")
    if not data_path.exists():
        raise FileNotFoundError(f"Data file not found: {data_path}")
    header: dict = {}
    try:
        header = read_dat_header(data_path)
    except OSError:
        pass
    return inp, inp_path, data_path, header


def _resolve_scaling_config(payload: dict, inp, header: dict) -> ScalingConfig:
    def pick(key: str, fallback, **rules):
        value = _payload_float(payload, key, **rules)
        return fallback if value is None else value

    if inp is not None:
        qmin = pick("qmin", inp.qmin)
        qmax = pick("qmax", inp.qmax)
        rho0 = pick("rho0", inp.rho0)
        b_avg_sq = pick("bAvgSq", inp.b_avg_sq)
        r_cutoff = pick("rCutoff", inp.r_cutoff)
        rmax = pick("rmax", inp.rmax)
        nr = int(pick("nr", inp.nr, integer=True))
        lorch = _payload_bool(payload, "lorch", inp.lorch)
    else:
        qmin = _payload_float(payload, "qmin")
        qmax = _payload_float(payload, "qmax")
        if qmin is None or qmax is None:
            raise CliError("data mode requires qmin and qmax")
        rho0 = _payload_float(payload, "rho0")
        if rho0 is None:
            rho0 = header.get("number_density")
        if rho0 is None:
            mass_density = _payload_float(payload, "massDensity")
            formula_raw = (payload.get("formula") or "").strip()
            if mass_density is not None and formula_raw:
                rho0 = number_density_from_mass_density(formula_raw, mass_density)
        if rho0 is None:
            raise CliError(
                "number density unknown: set rho0, or massDensity with formula, "
                "or use a data file with a NUMBER_DENSITY :: header"
            )
        b_avg_sq = _payload_float(payload, "bAvgSq")
        r_cutoff = pick("rCutoff", 1.0)
        rmax = pick("rmax", 50.0)
        nr = int(pick("nr", 5000, integer=True))
        lorch = _payload_bool(payload, "lorch", False)

    # One consistent source for <b>^2 and <b^2> (CLI parity): a formula's <b^2>
    # is never paired with a <b>^2 from another radiation/normalization.
    if _payload_float(payload, "bAvgSq") is not None:
        b_avg_sq_source = "bAvgSq"
    elif inp is not None:
        b_avg_sq_source = "stog.inp"
    else:
        b_avg_sq_source = None
    resolved = resolve_coefficients(
        b_avg_sq=b_avg_sq, b_avg_sq_source=b_avg_sq_source,
        b_sq_avg=_payload_float(payload, "bSqAvg"),
        formula=payload.get("formula"),
    )
    b_avg_sq, b_sq_avg = resolved["b_avg_sq"], resolved["b_sq_avg"]
    if b_avg_sq is None:
        raise CliError("data mode requires <b>^2: set bAvgSq or formula")

    r0 = _payload_float(payload, "r0")
    if r0 is None and "min_distance" in header:
        r0 = float(header["min_distance"])
    if r0 is None and inp is not None:
        r0 = stog_inp_closest_approach(inp, float(r_cutoff))

    config = ScalingConfig(
        qmin=float(qmin),
        qmax=float(qmax),
        rho0=float(rho0),
        b_avg_sq=float(b_avg_sq),
        b_sq_avg=None if b_sq_avg is None else float(b_sq_avg),
        r_cutoff=float(r_cutoff),
        r0=r0,
        r_fit_min=_payload_float(payload, "rFitMin"),
        r_fit_max=_payload_float(payload, "rFitMax"),
        rmax=float(rmax),
        nr=nr,
        lorch=lorch,
        low_q_correction=_payload_bool(payload, "lowQCorrection", True),
        robust=_payload_bool(payload, "robust", True),
        c1_mode=(payload.get("c1Mode") or "sweep").lower(),
        amplitude_criterion=(payload.get("amplitude") or "density").lower(),
        despike=_payload_bool(payload, "despike", False),
    )
    config.r_fit_window  # validate eagerly with a clean 400
    return config


def _scaling_enforce_flag(payload: dict) -> bool | None:
    """The ``enforce`` flag, parsed ONCE as a tri-state (CLI ``--enforce`` parity).

    None (absent / empty: enforcement on by default), True, or False — with
    the same string forms as every other boolean of this API ("false", "0"
    and 0 are False).
    """
    if payload.get("enforce") in (None, ""):
        return None
    return _payload_bool(payload, "enforce", True)


def _resolve_scaling_enforcement(
    payload: dict, inp, enforce_flag: bool | None
) -> tuple[float, float, float] | None:
    """Explicit enforcement triple, mirroring ``scaling_cli._resolve_enforcement``.

    Returns None when enforcement is off (``enforce_flag is False``) or when no
    explicit cutoff exists (then the caller applies the automatic first-shell
    cutoff unless ``enforce_flag is False``).
    """
    window = payload.get("peakWindow")
    if window is not None and not (isinstance(window, (list, tuple)) and len(window) == 2):
        raise CliError("peakWindow must be a two-element list [rmin, rmax]")
    explicit_cutoff = _payload_float(payload, "enforceCutoff")
    if enforce_flag is False:
        if explicit_cutoff is not None or window is not None:
            raise CliError("enforce=false contradicts enforceCutoff/peakWindow")
        return None
    cutoff = explicit_cutoff
    if cutoff is None and inp is not None:
        cutoff = inp.peak_cutoff
    if cutoff is None:
        if window is not None:
            raise CliError("peakWindow requires enforceCutoff (or a stog input)")
        return None  # resolved post-run: the automatic first-shell cutoff
    if window is not None:
        peak_rmin = _number(window[0], "peakWindow[0]")
        peak_rmax = _number(window[1], "peakWindow[1]")
    elif inp is not None and explicit_cutoff is None:
        peak_rmin, peak_rmax = inp.peak_rmin, inp.peak_rmax
    else:
        peak_rmin = peak_rmax = cutoff
    return float(cutoff), float(peak_rmin), float(peak_rmax)


def _resolve_scaling_mode(payload: dict, inp) -> tuple[str, float, float]:
    mode = (payload.get("mode") or "auto").lower()
    if mode not in ("auto", "manual"):
        raise CliError(f"mode must be 'auto' or 'manual', got {mode!r}")
    a = _payload_float(payload, "a")
    b = _payload_float(payload, "b")
    if mode == "manual":
        if a is None and inp is not None:
            a = inp.a
            if b is None:
                b = inp.b
        if a is None:
            raise CliError("manual mode requires a scale 'a' (or a stog input to take it from)")
        if not math.isfinite(float(a)) or float(a) == 0.0:
            # S_corr = a*S_meas + b: a = 0 discards the data (and the raw-S(Q)
            # series (S_corr - b)/a would be 0/0).
            raise ValueError(f"manual mode requires a finite, non-zero scale 'a', got {a}")
        if b is None:
            b = 0.0
        return mode, float(a), float(b)
    return mode, 0.0, 0.0


def _cached_scaling(data_path: Path, config: ScalingConfig, mode: str, a: float, b: float, use_sigma: bool):
    return _SCALING_CACHE.get(
        data_path,
        (config, mode, a, b, use_sigma),
        lambda: _compute_scaling(str(data_path), config, mode, a, b, use_sigma),
    )


def _compute_scaling(path_str: str, config: ScalingConfig, mode: str, a: float, b: float, use_sigma: bool):
    data = read_stog_xy(path_str)
    q, sq = data[0], data[1]
    sigma = None
    if use_sigma and data.shape[0] >= 3:
        sigma = usable_sigma(q, sq, data[2])  # shared CLI/page guard
    if mode == "manual":
        return scale_pipeline(q, sq, config, a, b)
    return autoscale(q, sq, config, sigma=sigma)


def _scaling_request(payload: dict):
    inp, inp_path, data_path, header = _resolve_scaling_source(payload)
    config = _resolve_scaling_config(payload, inp, header)
    enforce_flag = _scaling_enforce_flag(payload)
    enforcement = _resolve_scaling_enforcement(payload, inp, enforce_flag)
    mode, a, b = _resolve_scaling_mode(payload, inp)
    use_sigma = _payload_bool(payload, "useSigma", True)
    result = _cached_scaling(data_path, config, mode, a, b, use_sigma)
    # The cached ScalingResult is shared by every identical request (and
    # request thread): annotate a per-request copy, never the cached object.
    from dataclasses import replace as _replace_result

    result = _replace_result(result, provenance=dict(result.provenance))
    # No explicit cutoff and enforcement not refused: enforce automatically
    # at the foot of the first shell (CLI-mirroring auto default).
    if enforcement is None and enforce_flag is not False:
        r0_detected = result.provenance.get("r0_detected")
        if r0_detected is None:
            r0_detected = detect_first_peak_onset(
                result.r, result.g_filtered, config.qmax,
                search_min=config.r_cutoff + 0.3,
            )
            if r0_detected is not None:
                result.provenance["r0_detected"] = float(r0_detected)
        cutoff = auto_enforcement_cutoff(
            result.r, result.g_filtered, config, onset=r0_detected
        )
        if cutoff is not None:
            enforcement = (cutoff,) * 3
    return inp, inp_path, data_path, header, config, enforcement, mode, result


def _header_payload(header: dict) -> dict:
    return {
        "title": header.get("title"),
        "numberDensity": header.get("number_density"),
        "minDistance": header.get("min_distance"),
    }


def _inp_payload(inp) -> dict | None:
    if inp is None:
        return None
    return {
        "a": inp.a,
        "b": inp.b,
        "qmin": inp.qmin,
        "qmax": inp.qmax,
        "rho0": inp.rho0,
        "bAvgSq": inp.b_avg_sq,
        "rCutoff": inp.r_cutoff,
        "rmax": inp.rmax,
        "nr": inp.nr,
        "lorch": inp.lorch,
        "dataFile": inp.data_file,
        "peakCutoff": inp.peak_cutoff,
        "peakRmin": inp.peak_rmin,
        "peakRmax": inp.peak_rmax,
    }


@app.route("/api/scaling/preview", methods=["POST"])
def scaling_preview():
    try:
        payload = request.get_json(silent=True) or {}
        if payload.get("inspect"):
            inp, inp_path, data_path, header = _resolve_scaling_source(payload)
            return jsonify(
                {
                    "source": str(inp_path or data_path),
                    "kind": "inp" if inp is not None else "data",
                    "dataFile": str(data_path),
                    "inp": _inp_payload(inp),
                    "header": _header_payload(header),
                }
            )

        inp, inp_path, data_path, header, config, enforcement, mode, result = _scaling_request(payload)
        summary = diagnostics_summary(result, config)

        gk_enforced = dr_enforced = None
        if enforcement is not None:
            cutoff, peak_rmin, peak_rmax = enforcement
            g_final = first_peak_zero(
                result.r, result.g_filtered,
                cutoff=cutoff, peak_rmin=peak_rmin, peak_rmax=peak_rmax,
            )
            gk_array = g_to_gk(g_final, config.b_avg_sq)
            gk_enforced = gk_array.tolist()
            dr_enforced = gk_to_dr(result.r, gk_array, config.rho0).tolist()

        response = {
            "source": str(inp_path or data_path),
            "dataFile": str(data_path),
            "kind": "inp" if inp is not None else "data",
            "mode": mode,
            "inp": _inp_payload(inp),
            "header": _header_payload(header),
            "result": {
                "a": result.a,
                "b": result.b,
                "converged": result.converged,
                "iterations": result.iterations,
                "lowRRms": result.low_r_rms,
                "c1TailMean": result.c1_tail_mean,
                "history": [list(item) for item in result.history],
            },
            "diagnostics": _json_safe(summary),
            "provenance": _json_safe(result.provenance),
            "enforcement": None
            if enforcement is None
            else dict(zip(("cutoff", "peakRmin", "peakRmax"), enforcement)),
            "guides": _json_safe(
                {
                    "asymptote": 1.0,
                    "gkLowR": -config.b_avg_sq,
                    "drSlope": -4.0 * math.pi * config.rho0 * config.b_avg_sq,
                    "s0Target": None
                    if config.b_sq_avg is None
                    else 1.0 - config.b_sq_avg / config.b_avg_sq,
                    "qTailMin": config.q_tail_min,
                    "rFitWindow": summary.get("r_fit_window", list(config.r_fit_window)),
                    "r0Detected": summary.get("r0_detected"),
                    "level": None if result.sweep is None else result.sweep.level,
                    "levelWindow": None
                    if result.sweep is None
                    else [result.sweep.q_lo, result.sweep.q_hi],
                }
            ),
            "series": {
                "q": result.q.tolist(),
                "sqRaw": ((result.sq_scaled - result.b) / result.a).tolist(),
                "sqScaled": result.sq_scaled.tolist(),
                "sqFiltered": result.sq_filtered.tolist(),
                "r": result.r.tolist(),
                "gk": result.gk.tolist(),
                "dr": result.d_r.tolist(),
                "gkEnforced": gk_enforced,
                "drEnforced": dr_enforced,
            },
        }
        return jsonify(response)
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except (CliError, ValueError, NotImplementedError) as exc:
        return jsonify({"error": str(exc)}), 400
    except SourceChangedError as exc:
        return jsonify({"error": str(exc)}), 409
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/scaling/run", methods=["POST"])
def scaling_run():
    try:
        payload = request.get_json(silent=True) or {}
        inp, inp_path, data_path, header, config, enforcement, mode, result = _scaling_request(payload)
        if mode == "auto":
            refuse_failed_fit(result)  # a <= 0 never becomes RMCProfile input (CLI parity)
        summary = diagnostics_summary(result, config)

        from types import SimpleNamespace

        out_dir = payload.get("outDir")
        args = SimpleNamespace(
            out_dir=str(_resolve_inside_root(out_dir)) if out_dir else None,
            out_stem=(payload.get("outStem") or "").strip() or None,
            force=_payload_bool(payload, "force", False),
        )
        targets = _resolve_scaling_targets(args, inp, inp_path, data_path)
        for target in targets.values():
            if not any(
                _is_inside_root(target.resolve().parent, root) or target.resolve().parent == root
                for root in (DATA_ROOT, *SELECTED_DATA_ROOTS)
            ):
                raise PermissionError("Output directory is outside configured data roots")

        provenance_payload = {
            "tool": "rmc-autoscale (web API)",
            "source": str(inp_path or data_path),
            "data_file": str(data_path),
            "stog_inp_reference": None if inp is None else {"a": inp.a, "b": inp.b},
            "enforcement": None
            if enforcement is None
            else dict(zip(("cutoff", "peak_rmin", "peak_rmax"), enforcement)),
            "outputs": {key: str(path) for key, path in targets.items()},
            "diagnostics": summary,
            "provenance": result.provenance,
        }
        targets["provenance"].parent.mkdir(parents=True, exist_ok=True)
        _write_outputs(result, config, targets, enforcement, provenance_payload)
        return jsonify(
            {
                "a": result.a,
                "b": result.b,
                "mode": mode,
                "outputs": {key: str(path) for key, path in targets.items()},
                "outDir": str(targets["provenance"].parent),
                "diagnostics": _json_safe(summary),
            }
        )
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except CliError as exc:
        status = 409 if "refusing to overwrite" in str(exc) else 400
        return jsonify({"error": str(exc)}), status
    except (ValueError, NotImplementedError) as exc:
        return jsonify({"error": str(exc)}), 400
    except SourceChangedError as exc:
        return jsonify({"error": str(exc)}), 409
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/", defaults={"path": ""})
@app.route("/<path:path>")
def serve_frontend(path: str):
    if path.startswith("api/"):
        return jsonify({"error": "API endpoint not found"}), 404

    static_folder = Path(app.static_folder or FRONTEND_DIST)
    requested = static_folder / path
    if path and requested.is_file():
        return send_from_directory(static_folder, path)

    index = static_folder / "index.html"
    if index.exists():
        return send_from_directory(static_folder, "index.html")

    return jsonify(
        {
            "error": "Frontend build not found",
            "hint": "Run `npm run build` in web_app/frontend or use the Dockerfile.",
        }
    ), 404


if __name__ == "__main__":
    port = int(os.environ.get("PORT", os.environ.get("RMC_TOOLKITS_PORT", 5000)))
    app.run(debug=True, host="0.0.0.0", port=port)
