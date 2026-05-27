"""File parsers used by the CLI scripts and web application."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re

import numpy as np


@dataclass(frozen=True)
class CsvSeries:
    labels: list[str]
    data: np.ndarray


@dataclass(frozen=True)
class RmcStructure:
    atom_indices: dict[str, list[int]]
    lattice_vectors: np.ndarray
    supercell: np.ndarray
    atom_types: list[str]
    positions: np.ndarray


def read_rmc_csv(path: str | Path) -> CsvSeries:
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        lines = handle.readlines()

    if not lines:
        raise ValueError(f"{path} is empty")

    labels = [label.strip() for label in lines[0].split(",")]
    rows: list[list[float]] = []
    for line in lines[1:]:
        values = [value.strip() for value in line.split(",") if value.strip()]
        if values:
            rows.append([float(value) for value in values])

    if not rows:
        raise ValueError(f"{path} does not contain numeric rows")

    return CsvSeries(labels=labels, data=np.asarray(rows, dtype=float).T)


def read_chi(paths: list[str | Path]) -> tuple[np.ndarray, np.ndarray]:
    chi_q: list[float] = []
    chi_r: list[float] = []
    for path in paths:
        with Path(path).open("r", encoding="utf-8") as handle:
            lines = handle.readlines()
        for line in lines[2:]:
            parts = line.split()
            if len(parts) >= 2:
                chi_q.append(float(parts[-2]))
                chi_r.append(float(parts[-1]))
    return np.asarray(chi_q, dtype=float), np.asarray(chi_r, dtype=float)


def read_stog(path: str | Path) -> np.ndarray:
    rows: list[list[float]] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle.readlines()[2:]:
            parts = line.split()
            if parts:
                rows.append([float(value) for value in parts])
    if not rows:
        raise ValueError(f"{path} does not contain STOG numeric rows")
    return np.asarray(rows, dtype=float).T


def pdf_index(path: str | Path) -> int:
    match = re.search(r"PDF(\d+)\.csv$", str(path))
    return int(match.group(1)) if match else 0


def rwp(x: np.ndarray, observed: np.ndarray, fitted: np.ndarray) -> float:
    denom = float(np.sum(observed * observed))
    if denom == 0:
        return 0.0
    residual = fitted - observed
    return float(np.sqrt(np.sum(residual * residual) / denom))


def read_atom_indices(rmc6f_path: str | Path) -> dict[str, list[int]]:
    lines = Path(rmc6f_path).read_text(encoding="utf-8", errors="replace").splitlines()
    start = next((idx for idx, line in enumerate(lines) if line.split()[:1] == ["Atoms:"]), None)
    if start is None:
        raise ValueError(f"{rmc6f_path} does not contain an Atoms section")

    atom_indices: dict[str, set[int]] = {}
    for line in lines[start + 1 :]:
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            atom_index = int(parts[-4])
        except ValueError:
            continue
        atom_indices.setdefault(parts[1], set()).add(atom_index)

    return {atom: sorted(indices) for atom, indices in atom_indices.items()}


def read_cell_vectors(rmc6f_path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    lines = Path(rmc6f_path).read_text(encoding="utf-8", errors="replace").splitlines()
    lattice_vectors: np.ndarray | None = None
    supercell: np.ndarray | None = None

    for idx, line in enumerate(lines):
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "Supercell":
            supercell = np.asarray(parts[-3:], dtype=float)
        elif parts[0] == "Lattice":
            lattice_vectors = np.asarray(
                [lines[idx + 1].split(), lines[idx + 2].split(), lines[idx + 3].split()],
                dtype=float,
            )

    if lattice_vectors is None or supercell is None:
        raise ValueError(f"{rmc6f_path} is missing lattice or supercell metadata")

    return lattice_vectors, supercell


def iter_rmc6f_atoms(rmc6f_path: str | Path):
    """Yield atom records from an RMCProfile `.rmc6f` file."""
    in_atoms = False
    with Path(rmc6f_path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "Atoms:":
                in_atoms = True
                continue
            if not in_atoms:
                continue
            if len(parts) < 10:
                continue
            try:
                yield {
                    "atom_number": int(parts[0]),
                    "element": parts[1],
                    "type_label": parts[2],
                    "coords": np.asarray(parts[3:6], dtype=float),
                    "reference_number": int(parts[6]),
                    "cell_indices": np.asarray(parts[7:10], dtype=int),
                }
            except ValueError:
                continue


def frac_lines_from_rmc6f(rmc6f_path: str | Path) -> list[str]:
    """Build `Frac_coord*.txt` content from an RMCProfile `.rmc6f` file."""
    _, supercell = read_cell_vectors(rmc6f_path)
    lines = [
        " RN - reference number (a column in rmc6f file indicating an atom type\n",
        " in the unit cell)\n",
        " XYZ - fractional coordinates of the atom reduced to unit cell\n",
        " Nx,Ny,Nz - unit cell indices in the box\n",
        " RN    X    Y     Z    Nx    Ny    Nz\n",
    ]

    for atom in iter_rmc6f_atoms(rmc6f_path):
        reduced = atom["coords"] - (atom["cell_indices"] / supercell)
        rn = atom["reference_number"]
        nx, ny, nz = atom["cell_indices"]
        lines.append(
            f"{rn:3d}    {reduced[0]:.5f}    {reduced[1]:.5f}    {reduced[2]:.5f}  "
            f"{nx:d}  {ny:d}  {nz:d}\n"
        )
    return lines


def write_frac_from_rmc6f(
    rmc6f_path: str | Path,
    output_path: str | Path | None = None,
    overwrite: bool = False,
) -> Path:
    """Write a `Frac_coord*.txt` file derived from an RMCProfile `.rmc6f` file."""
    rmc6f_path = Path(rmc6f_path)
    if output_path is None:
        output_path = rmc6f_path.with_name(f"Frac_coord_{rmc6f_path.stem}.txt")
    output_path = Path(output_path)
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"{output_path} already exists; pass overwrite=True to replace it")

    lines = frac_lines_from_rmc6f(rmc6f_path)
    output_path.write_text("".join(lines), encoding="utf-8")
    return output_path


def read_structure(
    directory: str | Path,
    element: str | int | None = None,
    mode: str = "cartesian",
) -> RmcStructure:
    directory = Path(directory)
    frac_path = next(iter(sorted(directory.glob("Frac*.txt"))), None)
    rmc6f_path = next(iter(sorted(directory.glob("*.rmc6f"))), None)
    if frac_path is None:
        raise FileNotFoundError(f"No Frac*.txt file found in {directory}")
    if rmc6f_path is None:
        raise FileNotFoundError(f"No .rmc6f file found in {directory}")

    atom_indices = read_atom_indices(rmc6f_path)
    lattice_vectors, supercell = read_cell_vectors(rmc6f_path)
    unit_vectors = lattice_vectors / supercell[:, None]
    selected_indices = None if element in (None, 0, "0", "all") else set(atom_indices[str(element)])

    atom_types: list[str] = []
    positions: list[np.ndarray] = []
    with frac_path.open("r", encoding="utf-8") as handle:
        lines = handle.readlines()[5:]

    for line in lines:
        parts = line.split()
        if len(parts) < 4:
            continue
        atom_id = int(parts[0])
        if selected_indices is not None and atom_id not in selected_indices:
            continue
        frac = np.asarray(parts[1:4], dtype=float) * supercell
        folded = np.mod(frac, 1.0)
        atom_types.append(parts[0])
        if mode == "fractional":
            positions.append(folded)
        else:
            positions.append(folded[0] * unit_vectors[0] + folded[1] * unit_vectors[1] + folded[2] * unit_vectors[2])

    return RmcStructure(
        atom_indices=atom_indices,
        lattice_vectors=lattice_vectors,
        supercell=supercell,
        atom_types=atom_types,
        positions=np.asarray(positions, dtype=float),
    )
