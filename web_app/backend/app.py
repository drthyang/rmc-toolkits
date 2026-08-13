# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import io
import math
import os
import platform
import re
import shutil
import subprocess
import sys

import numpy as np
from flask import Flask, jsonify, request, send_file, send_from_directory
from flask_cors import CORS

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
FRONTEND_DIST = PROJECT_ROOT / "web_app" / "frontend" / "dist"

from rmc_toolkits.kde import UnitCellPositions, load_unit_cell_positions, oriented_kde_slice
from rmc_toolkits.orientation import site_orientation_histogram
from rmc_toolkits.pca_kde import (
    SiteDisplacements,
    cached_site_displacements,
    site_ellipsoids,
    site_pca_kde,
)
from rmc_toolkits.triplets import cached_bond_angle_summary
from rmc_toolkits.parsers import (
    iter_rmc6f_atoms,
    read_atom_indices,
    read_cell_vectors,
    read_moves_metadata,
    read_chi,
    read_dat_header,
    read_exafs_csv,
    read_rmc_csv,
    read_stog,
    read_stog_inp,
    read_stog_xy,
    related_r_value_logs,
    write_frac_from_rmc6f,
)
from rmc_toolkits.plots import bragg_is_tof, close_plot, detect_plot_kind, make_plot, plot_to_png
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
)
from rmc_toolkits.scaling import detect_first_peak_onset
from rmc_toolkits.scattering import faber_ziman, number_density_from_mass_density
from rmc_toolkits.transforms import first_peak_zero, g_to_gk, gk_to_dr


app = Flask(__name__, static_folder=str(FRONTEND_DIST), static_url_path="")
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
    rmc6f_files = sorted(directory.glob("*.rmc6f"))
    if not rmc6f_files:
        raise FileNotFoundError(f"No .rmc6f file found in {directory}")
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


def _sample_atoms_by_site(atoms: list[dict], max_points: int) -> tuple[list[dict], int]:
    if len(atoms) <= max_points:
        return atoms, 1

    by_reference: dict[int, list[dict]] = {}
    for atom in atoms:
        by_reference.setdefault(atom["reference_number"], []).append(atom)

    quota = max(1, max_points // len(by_reference))
    sampled: list[dict] = []
    for reference_number in sorted(by_reference):
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
            _, chi_r = read_chi(related_r_value_logs(path))
            y_values = [float(value) for value in chi_r]
            return jsonify(
                {
                    **metadata,
                    "xLabel": "Time steps",
                    "yLabel": "log(χ)",
                    "series": [
                        {
                            "label": "R",
                            "x": list(range(len(y_values))),
                            "y": [float(math.log(max(value, 1e-12))) for value in y_values],
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
                    "yLabel": "G(r)" if path.name.endswith(".gr") else "S(Q)",
                    "series": [{"label": path.name, "x": data[0].tolist(), "y": data[1].tolist()}],
                }
            )

        series = read_exafs_csv(path) if kind in ("exafs_q", "exafs_r") else read_rmc_csv(path)
        x_values = series.data[0].tolist()
        payload_series = []
        for idx, label in enumerate(series.labels[1:], start=1):
            if idx < len(series.data):
                payload_series.append({"label": label.strip() or f"Series {idx}", "x": x_values, "y": series.data[idx].tolist()})

        x_label = series.labels[0] if series.labels else "x"
        if kind == "exafs_q":
            x_label = "k (Å^{-1})"
            y_label = "χ(k) k²"
        elif kind == "exafs_r":
            x_label = "r (Å)"
            y_label = "FT[χ(k) k²]"
        elif kind in ("xpdf", "npdf", "pdf_partials"):
            x_label = "r (Å)"
            y_label = "G(r)"
        elif kind in ("xray_sq", "neutron_sq"):
            x_label = "Q (Å^{-1})"
            y_label = "S(Q)"
        elif kind == "bragg":
            x_label = "ToF (µs)" if bragg_is_tof(series.labels[0] if series.labels else None) else "Q (Å^{-1})"
            y_label = "Intensity"
        else:
            x_label = _clean_axis_label(x_label)
            y_label = "data"

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
        max_points = max(100, min(int(request.args.get("maxPoints", MAX_STRUCTURE_POINTS)), MAX_STRUCTURE_POINTS))
        rmc6f_path = _find_rmc6f(target)
        lattice_vectors, supercell = read_cell_vectors(rmc6f_path)
        atom_indices = read_atom_indices(rmc6f_path)
        moves = read_moves_metadata(rmc6f_path)

        atoms = list(iter_rmc6f_atoms(rmc6f_path))
        sampled, stride = _sample_atoms_by_site(atoms, max_points)
        points = []
        counts: dict[str, int] = {}
        for atom in atoms:
            counts[atom["element"]] = counts.get(atom["element"], 0) + 1
        for atom in sampled:
            reduced = atom["coords"] - (atom["cell_indices"] / supercell)
            unit_cell = (reduced * supercell) % 1.0
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
                "points": points,
            }
        )
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@lru_cache(maxsize=16)
def _cached_positions(path_str: str, mtime: float, element: str | None) -> UnitCellPositions:
    return load_unit_cell_positions(path_str, element=element)


SLICE_ORIENTATIONS = {
    "a": {"normal": (1.0, 0.0, 0.0), "u": (0.0, 1.0, 0.0), "v": (0.0, 0.0, 1.0)},
    "b": {"normal": (0.0, 1.0, 0.0), "u": (1.0, 0.0, 0.0), "v": (0.0, 0.0, 1.0)},
    "c": {"normal": (0.0, 0.0, 1.0), "u": (1.0, 0.0, 0.0), "v": (0.0, 1.0, 0.0)},
}


def _bw_argument(raw: str) -> str | float:
    """Bandwidth query arg: the names 'scott'/'silverman', else a positive float."""
    name = str(raw).strip().lower()
    if name in ("scott", "silverman"):
        return name
    return float(raw)


def _slice_orientation_from_request():
    orientation = (request.args.get("orientation") or "c").lower()
    if orientation in SLICE_ORIENTATIONS:
        config = SLICE_ORIENTATIONS[orientation]
        return orientation, config["normal"], config["u"], config["v"]

    normal = (
        float(request.args.get("nx", 0.0)),
        float(request.args.get("ny", 0.0)),
        float(request.args.get("nz", 1.0)),
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
        positions = _cached_positions(str(rmc6f_path), rmc6f_path.stat().st_mtime, element)
        cell_lengths = positions.cell_lengths

        orientation, normal, u_axis, v_axis = _slice_orientation_from_request()

        # z and dz arrive as fractions of the projection range along the slice
        # normal. Keep the KDE slice in fractional coordinates so non-orthogonal
        # cells can be projected through the actual cell basis in the frontend.
        z_frac = float(request.args.get("z", 0.5))
        dz_frac = float(request.args.get("dz", 0.08))
        bw = float(request.args.get("bw", 0.03))
        grid = int(request.args.get("grid", 120))
        levels = int(request.args.get("levels", 8))
        log = request.args.get("log", "false").lower() in ("1", "true", "yes")

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
        return jsonify(result)
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


def _cached_site_displacements(rmc6f_path: Path) -> SiteDisplacements:
    return cached_site_displacements(str(rmc6f_path), rmc6f_path.stat().st_mtime)


@app.route("/api/pca/sites", methods=["GET"])
def pca_sites_endpoint():
    """Anisotropic displacement tensor + thermal ellipsoid for every site.

    One cheap batched pass over the whole configuration; the frontend uses it to
    populate the site picker and the per-site ellipsoid summary table.
    """
    try:
        target = _resolve_inside_root(request.args.get("dir", "."))
        rmc6f_path = _find_rmc6f(target)
        probability = float(request.args.get("probability", 0.5))
        sites = _cached_site_displacements(rmc6f_path)
        ellipsoids = site_ellipsoids(sites, probability=probability)
        return jsonify(
            {
                "source": str(rmc6f_path),
                "referenceNumbers": sites.reference_numbers.tolist(),
                "elements": sorted(set(sites.elements)),
                "totalAtoms": int(sites.counts.sum()),
                "latticeVectors": sites.lattice_vectors.tolist(),
                "supercell": sites.supercell.tolist(),
                "probability": probability,
                "sites": ellipsoids,
            }
        )
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/pca/kde", methods=["GET"])
def pca_kde_endpoint():
    """PCA-aligned 3D KDE volume for one site (or one element's pooled sites)."""
    try:
        target = _resolve_inside_root(request.args.get("dir", "."))
        rmc6f_path = _find_rmc6f(target)
        sites = _cached_site_displacements(rmc6f_path)

        reference_raw = request.args.get("referenceNumber")
        reference_number = int(reference_raw) if reference_raw not in (None, "") else None
        element = request.args.get("element") or None
        if element in ("", "all"):
            element = None

        result = site_pca_kde(
            sites,
            reference_number=reference_number,
            element=element,
            bw=_bw_argument(request.args.get("bw", "scott")),
            bw_scale=float(request.args.get("bwScale", 1.0)),
            grid=int(request.args.get("grid", 48)),
            extent=float(request.args.get("extent", 3.0)),
            cubic_box=request.args.get("cubicBox", "false").lower() in ("1", "true", "yes"),
            probability=float(request.args.get("probability", 0.5)),
            projections=request.args.get("projections", "true").lower() in ("1", "true", "yes"),
        )
        result["source"] = str(rmc6f_path)
        return jsonify(result)
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/triplets", methods=["GET"])
def triplets_endpoint():
    """Bond-angle (triplet) summary of the run's configuration.

    Params: end1/apex/end2 (elements, apex central), r12Min/r12Max and the
    optional r23Min/r23Max windows (angstrom, inclusive), binWidth (degrees).
    Payload shape is defined by rmc_toolkits.triplets.bond_angle_summary and
    mirrored by the browser worker's 'triplets' request; both boundaries cap
    rmax at 15 A and binWidth at >= 0.05 deg (the engine itself is
    unrestricted for library/CLI use) so one request cannot blow up the
    stencil volume or the response size.
    """
    try:
        target = _resolve_inside_root(request.args.get("dir", "."))
        rmc6f_path = _find_rmc6f(target)

        def window(min_key, max_key, fallback=None):
            raw_min, raw_max = request.args.get(min_key), request.args.get(max_key)
            missing = [key for key, raw in ((min_key, raw_min), (max_key, raw_max)) if raw in (None, "")]
            if len(missing) == 2 and fallback is not None:
                return fallback
            if missing:
                raise ValueError(f"{min_key}/{max_key} are required together; missing {missing[0]}")
            bounds = float(raw_min), float(raw_max)
            if bounds[1] > 15.0:
                raise ValueError(f"{max_key} is capped at 15 A for API requests, got {bounds[1]}")
            return bounds

        window12 = window("r12Min", "r12Max")
        window23 = window("r23Min", "r23Max", fallback=window12)
        bin_width = float(request.args.get("binWidth", 1.0))
        if bin_width < 0.05:
            raise ValueError(f"binWidth is capped at >= 0.05 deg for API requests, got {bin_width}")
        result = dict(
            cached_bond_angle_summary(
                str(rmc6f_path),
                rmc6f_path.stat().st_mtime,
                # Normalized here so 'se' and 'Se' share one cache entry.
                request.args.get("end1", "").strip().capitalize(),
                request.args.get("apex", "").strip().capitalize(),
                request.args.get("end2", "").strip().capitalize(),
                *window12,
                *window23,
                bin_width,
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
        sites = _cached_site_displacements(rmc6f_path)

        reference_raw = request.args.get("referenceNumber")
        reference_number = int(reference_raw) if reference_raw not in (None, "") else None
        element = request.args.get("element") or None
        if element in ("", "all"):
            element = None

        frequency_raw = request.args.get("frequency")
        result = site_orientation_histogram(
            sites,
            reference_number=reference_number,
            element=element,
            frequency=int(frequency_raw) if frequency_raw not in (None, "") else None,
            weight=request.args.get("weight", "count"),
            min_amplitude=float(request.args.get("minAmplitude", 0.0)),
            min_amplitude_quantile=float(request.args.get("minAmplitudeQuantile", 0.0)),
            smoothing=int(request.args.get("smoothing", 0)),
            frame=request.args.get("frame", "cartesian"),
            geometry=request.args.get("geometry", "true").lower() in ("1", "true", "yes"),
        )
        result["source"] = str(rmc6f_path)
        return jsonify(result)
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


# --- Auto StoG scaling API ---------------------------------------------------
# Thin HTTP face over rmc_toolkits.scaling; file writing is shared with the
# rmc-autoscale CLI so API outputs are identical to a CLI/classic-stog session.


def _payload_float(payload: dict, key: str) -> float | None:
    value = payload.get(key)
    if value in (None, ""):
        return None
    return float(value)


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
    def pick(key: str, fallback):
        value = _payload_float(payload, key)
        return fallback if value is None else value

    if inp is not None:
        qmin = pick("qmin", inp.qmin)
        qmax = pick("qmax", inp.qmax)
        rho0 = pick("rho0", inp.rho0)
        b_avg_sq = pick("bAvgSq", inp.b_avg_sq)
        r_cutoff = pick("rCutoff", inp.r_cutoff)
        rmax = pick("rmax", inp.rmax)
        nr = int(pick("nr", inp.nr))
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
        nr = int(pick("nr", 5000))
        lorch = _payload_bool(payload, "lorch", False)

    b_sq_avg = _payload_float(payload, "bSqAvg")
    formula = (payload.get("formula") or "").strip()
    if formula:
        coefficients = faber_ziman(formula)
        if b_sq_avg is None:
            b_sq_avg = coefficients.b_sq_avg_barn
        if b_avg_sq is None:
            b_avg_sq = coefficients.b_avg_sq_barn
    if b_avg_sq is None:
        raise CliError("data mode requires <b>^2: set bAvgSq or formula")

    r0 = _payload_float(payload, "r0")
    if r0 is None and "min_distance" in header:
        r0 = float(header["min_distance"])
    if r0 is None and inp is not None:
        candidate = max(inp.peak_cutoff, inp.peak_rmin)
        if candidate - 0.25 > float(r_cutoff) + 0.2:
            r0 = candidate

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


def _resolve_scaling_enforcement(payload: dict, inp) -> tuple[float, float, float] | None:
    if _payload_bool(payload, "enforce", inp is not None) is False:
        return None
    cutoff = _payload_float(payload, "enforceCutoff")
    if cutoff is None and inp is not None:
        cutoff = inp.peak_cutoff
    if cutoff is None:
        return None  # data mode: resolved post-run from the detected r0
    window = payload.get("peakWindow")
    if isinstance(window, (list, tuple)) and len(window) == 2:
        peak_rmin, peak_rmax = float(window[0]), float(window[1])
    elif inp is not None and _payload_float(payload, "enforceCutoff") is None:
        peak_rmin, peak_rmax = inp.peak_rmin, inp.peak_rmax
    else:
        peak_rmin = peak_rmax = cutoff
    return float(cutoff), peak_rmin, peak_rmax


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
        if b is None:
            b = 0.0
        return mode, float(a), float(b)
    return mode, 0.0, 0.0


@lru_cache(maxsize=8)
def _cached_scaling(path_str: str, mtime: float, config: ScalingConfig, mode: str, a: float, b: float, use_sigma: bool):
    data = read_stog_xy(path_str)
    q, sq = data[0], data[1]
    sigma = None
    if use_sigma and data.shape[0] >= 3:
        sigma = data[2]
        usable = np.isfinite(q) & np.isfinite(sq)
        if not np.all(np.isfinite(sigma[usable])) or np.any(sigma[usable] <= 0):
            sigma = None
    if mode == "manual":
        return scale_pipeline(q, sq, config, a, b)
    return autoscale(q, sq, config, sigma=sigma)


def _scaling_request(payload: dict):
    inp, inp_path, data_path, header = _resolve_scaling_source(payload)
    config = _resolve_scaling_config(payload, inp, header)
    enforcement = _resolve_scaling_enforcement(payload, inp)
    mode, a, b = _resolve_scaling_mode(payload, inp)
    use_sigma = _payload_bool(payload, "useSigma", True)
    result = _cached_scaling(
        str(data_path), data_path.stat().st_mtime, config, mode, a, b, use_sigma
    )
    # No explicit cutoff and enforcement not refused: enforce at the
    # data-derived closest approach (CLI-mirroring auto default).
    if enforcement is None and payload.get("enforce") is not False:
        r0_detected = result.provenance.get("r0_detected")
        if r0_detected is None:
            r0_detected = detect_first_peak_onset(
                result.r, result.g_filtered, config.qmax,
                search_min=config.r_cutoff + 0.3,
            )
            if r0_detected is not None:
                result.provenance["r0_detected"] = float(r0_detected)
        if r0_detected is not None:
            enforcement = (float(r0_detected),) * 3
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
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


@app.route("/api/scaling/run", methods=["POST"])
def scaling_run():
    try:
        payload = request.get_json(silent=True) or {}
        inp, inp_path, data_path, header, config, enforcement, mode, result = _scaling_request(payload)
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
