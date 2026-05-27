from __future__ import annotations

from functools import lru_cache
from pathlib import Path
import io
import math
import os
import sys

from flask import Flask, jsonify, request, send_file
from flask_cors import CORS

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from rmc_toolkits.kde import UnitCellPositions, kde_slice, load_unit_cell_positions
from rmc_toolkits.parsers import (
    iter_rmc6f_atoms,
    read_atom_indices,
    read_cell_vectors,
    read_chi,
    read_rmc_csv,
    read_stog,
    write_frac_from_rmc6f,
)
from rmc_toolkits.plots import close_plot, detect_plot_kind, make_plot, plot_to_png


app = Flask(__name__)
CORS(app)

DATA_ROOT = Path(os.environ.get("RMC_TOOLKITS_DATA_ROOT", PROJECT_ROOT)).expanduser().resolve()
SUPPORTED_PATTERNS = ("*.csv", "*.log", "*.rmc6f", "Frac*.txt", "scale_ft.*", "stog_input.dat")


def _resolve_inside_root(raw_path: str | None) -> Path:
    candidate = Path(raw_path or ".").expanduser()
    if not candidate.is_absolute():
        candidate = DATA_ROOT / candidate
    resolved = candidate.resolve()
    if resolved != DATA_ROOT and DATA_ROOT not in resolved.parents:
        raise PermissionError(f"Path is outside configured data root: {DATA_ROOT}")
    return resolved


def _file_payload(path: Path, kind: str = "file") -> dict[str, str]:
    return {
        "name": path.name,
        "path": str(path),
        "type": kind,
        "plotKind": detect_plot_kind(path) if kind == "file" else None,
    }


def _find_rmc6f(directory: Path) -> Path:
    if directory.is_file() and directory.suffix == ".rmc6f":
        return directory
    rmc6f_files = sorted(directory.glob("*.rmc6f"))
    if not rmc6f_files:
        raise FileNotFoundError(f"No .rmc6f file found in {directory}")
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

        paths: dict[Path, dict[str, str]] = {}
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
            _, chi_r = read_chi([path])
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

        series = read_rmc_csv(path)
        x_values = series.data[0].tolist()
        payload_series = []
        for idx, label in enumerate(series.labels[1:], start=1):
            if idx < len(series.data):
                payload_series.append({"label": label.strip() or f"Series {idx}", "x": x_values, "y": series.data[idx].tolist()})

        x_label = series.labels[0] if series.labels else "x"
        if kind in ("xpdf", "npdf", "pdf_partials"):
            x_label = "r (Å)"
        elif kind == "bragg":
            x_label = "Q (Å^{-1})"
        else:
            x_label = _clean_axis_label(x_label)

        return jsonify({**metadata, "xLabel": x_label, "yLabel": "data", "series": payload_series})
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
        max_points = max(100, min(int(request.args.get("maxPoints", 12000)), 75000))
        rmc6f_path = _find_rmc6f(target)
        lattice_vectors, supercell = read_cell_vectors(rmc6f_path)
        atom_indices = read_atom_indices(rmc6f_path)

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

        # z and dz arrive as fractions of the cell edge (matching the slider
        # semantics in the frontend); convert to Angstrom for the KDE math.
        z_frac = float(request.args.get("z", 0.5))
        dz_frac = float(request.args.get("dz", 0.08))
        bw = float(request.args.get("bw", 0.03))
        grid = int(request.args.get("grid", 120))
        levels = int(request.args.get("levels", 8))
        log = request.args.get("log", "false").lower() in ("1", "true", "yes")

        result = kde_slice(
            positions.positions,
            z_frac * cell_lengths[2],
            dz_frac * cell_lengths[2],
            xlim=(0.0, float(cell_lengths[0])),
            ylim=(0.0, float(cell_lengths[1])),
            bw=bw,
            grid=grid,
            log=log,
            n_levels=levels,
        )
        result["cellLengths"] = cell_lengths.tolist()
        result["source"] = str(rmc6f_path)
        result["element"] = element or "all"
        return jsonify(result)
    except PermissionError as exc:
        return jsonify({"error": str(exc)}), 403
    except FileNotFoundError as exc:
        return jsonify({"error": str(exc)}), 404
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500


if __name__ == "__main__":
    app.run(debug=True, port=int(os.environ.get("RMC_TOOLKITS_PORT", 5000)))
