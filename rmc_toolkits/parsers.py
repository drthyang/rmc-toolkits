# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""File parsers used by the CLI scripts and web application."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Iterator, TypedDict

import numpy as np

R_VALUE_LOG_RE = re.compile(r"^(.+)-(\d{2,})\.log$")


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


class Rmc6fAtom(TypedDict):
    atom_number: int
    element: str
    type_label: str
    coords: np.ndarray
    reference_number: int
    cell_indices: np.ndarray


def read_rmc_csv(path: str | Path) -> CsvSeries:
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        lines = handle.readlines()

    if not lines:
        raise ValueError(f"{path} is empty")

    labels = [label.strip() for label in lines[0].split(",")]
    rows: list[list[float]] = []
    expected_columns = len(labels)
    for line_number, line in enumerate(lines[1:], start=2):
        values = [value.strip() for value in line.split(",") if value.strip()]
        if values:
            if len(values) != expected_columns:
                raise ValueError(
                    f"{path} line {line_number} has {len(values)} values; expected {expected_columns}"
                )
            rows.append([float(value) for value in values])

    if not rows:
        raise ValueError(f"{path} does not contain numeric rows")

    return CsvSeries(labels=labels, data=np.asarray(rows, dtype=float).T)


def _csv_values(line: str) -> list[str]:
    return [value.strip() for value in line.split(",") if value.strip()]


def _numeric_csv_values(line: str) -> list[float] | None:
    values = _csv_values(line)
    if not values:
        return None
    try:
        return [float(value) for value in values]
    except ValueError:
        return None


def read_exafs_csv(path: str | Path) -> CsvSeries:
    """Read RMCProfile EXAFS Q/R output CSV files.

    Q output files include a descriptive title row before the column header,
    while R output files start directly with the column header. The data rows
    are detected by scanning for the first fully numeric CSV row.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"{path} is empty")

    data_start = None
    for idx, line in enumerate(lines):
        if _numeric_csv_values(line) is not None:
            data_start = idx
            break
    if data_start is None or data_start == 0:
        raise ValueError(f"{path} does not contain an EXAFS column header and numeric rows")

    labels = _csv_values(lines[data_start - 1])
    rows: list[list[float]] = []
    expected_columns = len(labels)
    for line_number, line in enumerate(lines[data_start:], start=data_start + 1):
        values = _numeric_csv_values(line)
        if values is None:
            continue
        if len(values) != expected_columns:
            raise ValueError(
                f"{path} line {line_number} has {len(values)} values; expected {expected_columns}"
            )
        rows.append(values)

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
                try:
                    chi_q.append(float(parts[-2]))
                    chi_r.append(float(parts[-1]))
                except ValueError:
                    continue
    return np.asarray(chi_q, dtype=float), np.asarray(chi_r, dtype=float)


def r_value_log_parts(path: str | Path) -> tuple[str, int] | None:
    match = R_VALUE_LOG_RE.match(Path(path).name)
    if not match:
        return None
    return match.group(1), int(match.group(2))


def sort_r_value_logs(paths: list[str | Path]) -> list[Path]:
    def sort_key(path: str | Path) -> tuple[str, int, str]:
        parsed = r_value_log_parts(path)
        resolved = Path(path)
        if parsed:
            stem, sequence = parsed
            return stem.lower(), sequence, resolved.name.lower()
        return resolved.stem.lower(), -1, resolved.name.lower()

    return sorted((Path(path) for path in paths), key=sort_key)


def related_r_value_logs(path: str | Path) -> list[Path]:
    path = Path(path)
    parsed = r_value_log_parts(path)
    if not parsed or not path.parent.exists():
        return [path]

    stem, _ = parsed
    matches: list[Path] = []
    for candidate in path.parent.iterdir():
        candidate_parts = r_value_log_parts(candidate)
        if candidate.is_file() and candidate_parts and candidate_parts[0] == stem:
            matches.append(candidate)
    return sort_r_value_logs(matches) or [path]


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


@dataclass(frozen=True)
class StogInput:
    """Parsed classic stog/stog_new input file (the 23-line ``stog.inp``).

    ``yoffset``/``yscale`` follow the Fortran convention
    ``S_scaled = S_raw / yscale + yoffset``; the equivalent multiply-convention
    correction used by :mod:`rmc_toolkits.scaling` is exposed as ``a``/``b``.
    """

    n_files: int
    data_file: str
    qmin: float
    qmax: float
    yoffset: float
    yscale: float
    qoffset: float
    out_sq: str
    out_gr: str
    rmax: float
    nr: int
    lorch: bool
    rho0: float
    yoffset2: float
    try_again: bool
    use_filter: bool
    r_cutoff: float
    out_ft_sq: str
    out_ft_gr: str
    b_avg_sq: float
    out_rmc_fq: str
    out_rmc_gr: str
    out_rmc_dr: str
    peak_cutoff: float
    peak_rmin: float
    peak_rmax: float

    @property
    def a(self) -> float:
        return 1.0 / self.yscale

    @property
    def b(self) -> float:
        return self.yoffset


def _stog_flag(token: str) -> bool:
    return token.strip().upper().startswith("Y")


def read_stog_inp(path: str | Path) -> StogInput:
    """Parse a classic stog input file (single-dataset, filter-on layout).

    The interactive Fortran program's recorded input layout depends on the
    answers given; only the canonical layout exercised by the validation
    example is supported. Variants (multiple files, nonzero Q offset, the
    "try again" rescale loop, filter disabled) raise ``NotImplementedError``
    so silent misparses are impossible.
    """
    path = Path(path)
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    lines = [line for line in lines if line]
    if len(lines) < 22:
        raise ValueError(f"{path} has {len(lines)} non-empty lines; expected >= 22")

    n_files = int(lines[0].split()[0])
    if n_files != 1:
        raise NotImplementedError(f"{path}: only single-dataset inputs supported")
    qmin, qmax = (float(value) for value in lines[2].split()[:2])
    yoffset, yscale = (float(value) for value in lines[3].split()[:2])
    if yscale == 0 or not np.isfinite(yscale) or not np.isfinite(yoffset):
        raise ValueError(
            f"{path}: invalid yoffset/yscale line {lines[3]!r}; the Fortran "
            "convention divides by yscale, so it must be finite and nonzero"
        )
    qoffset = float(lines[4].split()[0])
    if qoffset != 0:
        raise NotImplementedError(f"{path}: nonzero Q offset not supported")
    yoffset2 = float(lines[11].split()[0])
    if yoffset2 != 0:
        raise NotImplementedError(f"{path}: nonzero second y offset not supported")
    try_again = _stog_flag(lines[12])
    if try_again:
        raise NotImplementedError(f"{path}: interactive 'try again' loops not supported")
    use_filter = _stog_flag(lines[13])
    if not use_filter:
        raise NotImplementedError(f"{path}: only filter-enabled inputs supported")
    peak_cutoff, peak_rmin, peak_rmax = (float(value) for value in lines[21].split()[:3])

    return StogInput(
        n_files=n_files,
        data_file=lines[1],
        qmin=qmin,
        qmax=qmax,
        yoffset=yoffset,
        yscale=yscale,
        qoffset=qoffset,
        out_sq=lines[5],
        out_gr=lines[6],
        rmax=float(lines[7].split()[0]),
        nr=int(lines[8].split()[0]),
        lorch=_stog_flag(lines[9]),
        rho0=float(lines[10].split()[0]),
        yoffset2=yoffset2,
        try_again=try_again,
        use_filter=use_filter,
        r_cutoff=float(lines[14].split()[0]),
        out_ft_sq=lines[15],
        out_ft_gr=lines[16],
        b_avg_sq=float(lines[17].split()[0]),
        out_rmc_fq=lines[18],
        out_rmc_gr=lines[19],
        out_rmc_dr=lines[20],
        peak_cutoff=peak_cutoff,
        peak_rmin=peak_rmin,
        peak_rmax=peak_rmax,
    )


def read_stog_xy(path: str | Path) -> np.ndarray:
    """Robustly read a whitespace-separated STOG-style x/y(/err) data file.

    Skips count headers, stray scalar lines, and text titles; keeps rows whose
    tokens all parse as floats with at least two columns (``NaN`` tokens are
    kept, so rebinned files retain their padding rows for the caller to mask).
    Returns the columns transposed, matching :func:`read_stog`.
    """
    groups: dict[int, list[list[float]]] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                values = [float(value) for value in parts]
            except ValueError:
                continue
            groups.setdefault(len(values), []).append(values)
    if not groups:
        raise ValueError(f"{path} does not contain STOG numeric rows")
    # Keep the modal column-count group so a stray numeric header line (e.g.
    # "count  qmin" on one line) cannot become the column template and
    # silently discard every real data row.
    rows = max(groups.values(), key=len)
    return np.asarray(rows, dtype=float).T


def read_dat_header(path: str | Path) -> dict[str, object]:
    """Parse ``KEY :: value`` metadata from an RMCProfile/STOG ``.dat`` header.

    Returns the raw string fields plus parsed conveniences:
    ``number_density`` (float, from ``NUMBER_DENSITY``) and ``min_distance``
    (float, the smallest ``MINIMUM_DISTANCES`` entry) when present.
    """
    raw: dict[str, str] = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "::" not in line:
                continue
            key, _, value = line.partition("::")
            raw[key.strip().upper()] = value.strip()

    result: dict[str, object] = {"raw": raw}
    if "TITLE" in raw:
        result["title"] = raw["TITLE"]
    density = raw.get("NUMBER_DENSITY")
    if density:
        for token in density.split():
            try:
                result["number_density"] = float(token)
                break
            except ValueError:
                continue
    distances = raw.get("MINIMUM_DISTANCES")
    if distances:
        values = []
        for token in distances.split():
            try:
                values.append(float(token))
            except ValueError:
                continue
        if values:
            result["min_distance"] = min(values)
    return result


def write_stog_xy(
    path: str | Path,
    x: np.ndarray,
    y: np.ndarray,
    *,
    title: str = "",
    extra: np.ndarray | None = None,
) -> Path:
    """Write x/y(/extra) columns in the classic STOG layout (count, title, rows).

    ``extra`` adds a third value column (e.g. the D(r) column of the classic
    ``scale_ft.gr``).
    """
    path = Path(path)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape:
        raise ValueError(f"x and y shapes differ: {x.shape} vs {y.shape}")
    if extra is not None:
        extra = np.asarray(extra, dtype=float)
        if extra.shape != x.shape:
            raise ValueError(f"extra column shape differs: {extra.shape} vs {x.shape}")
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"{x.size:>12d}\n")
        handle.write(f"{title}\n")
        for index, (xi, yi) in enumerate(zip(x, y)):
            row = f"  {xi:.16E}  {yi:.16E}"
            if extra is not None:
                row += f"  {extra[index]:.16E}"
            handle.write(row + "\n")
    return path


def pdf_index(path: str | Path) -> int:
    match = re.search(r"PDF(\d+)\.csv$", str(path))
    return int(match.group(1)) if match else 0


def rwp(x: np.ndarray, observed: np.ndarray, fitted: np.ndarray) -> float | None:
    """Weighted profile residual of ``fitted`` against ``observed``.

    Returns ``None`` — not ``0.0`` — for the two degenerate cases, so a caller
    can report the metric as unavailable instead of as a perfect fit: no finite
    ``(observed, fitted)`` pair (e.g. an all-NaN column), and a zero denominator,
    which offers no scale to normalize the residual against. Mirrors ``rwp()`` in
    ``web_app/frontend/src/browserData.js``.
    """
    observed = np.asarray(observed, dtype=float)
    fitted = np.asarray(fitted, dtype=float)
    paired = np.isfinite(observed) & np.isfinite(fitted)
    if not paired.any():
        return None

    obs = observed[paired]
    denom = float(np.dot(obs, obs))
    if denom == 0.0:
        return None

    residual = fitted[paired] - obs
    return float(np.sqrt(float(np.dot(residual, residual)) / denom))


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


_MOVE_COUNTERS = {
    "generated": re.compile(r"Number of moves generated:\s*([\d.]+)"),
    "tried": re.compile(r"Number of moves tried:\s*([\d.]+)"),
    "accepted": re.compile(r"Number of moves accepted:\s*([\d.]+)"),
    "accumulatedTimeS": re.compile(r"Accumulated time \(s\)[^:]*:\s*([\d.]+)"),
}


def read_moves_metadata(rmc6f_path: str | Path) -> dict[str, float] | None:
    """Run-history counters from a `.rmc6f` header.

    How many moves the run generated / tried / accepted, plus the accumulated running
    time. Taken together with the atom count these gauge sampling sufficiency — the
    raw totals mean little without the box size.

    Only the header is scanned, so this stays cheap on multi-megabyte configurations.
    Returns ``None`` when the file carries none of the counters; individual missing
    counters are simply absent from the mapping. Keys are camelCase to match the JSON
    the frontend already receives from the browser-side parser.
    """
    text = Path(rmc6f_path).read_text(encoding="utf-8", errors="replace")
    marker = text.find("Atoms:")
    header = text[: marker if marker > 0 else 4000]

    moves: dict[str, float] = {}
    for key, pattern in _MOVE_COUNTERS.items():
        match = pattern.search(header)
        if match:
            moves[key] = float(match.group(1))
    return moves or None


def iter_rmc6f_atoms(rmc6f_path: str | Path) -> Iterator[Rmc6fAtom]:
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
            # Index from the END so both the current format (with a bracketed type
            # label between the element and coordinates) and older files that omit
            # it parse: the reference number and three cell indices are always the
            # last four fields, the fractional coordinates the three before them.
            #   current:  id element [type] x y z ref cellx celly cellz  (10 fields)
            #   older:    id element        x y z ref cellx celly cellz  ( 9 fields)
            n = len(parts)
            if n < 9:
                continue
            try:
                yield {
                    "atom_number": int(parts[0]),
                    "element": parts[1].capitalize(),
                    "type_label": " ".join(parts[2 : n - 7]),
                    "coords": np.asarray(parts[n - 7 : n - 4], dtype=float),
                    "reference_number": int(parts[n - 4]),
                    "cell_indices": np.asarray(parts[n - 3 : n], dtype=int),
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
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("".join(lines), encoding="utf-8")
    return output_path


def read_structure(
    directory: str | Path,
    element: str | int | None = None,
    mode: str = "cartesian",
) -> RmcStructure:
    if mode not in {"cartesian", "fractional"}:
        raise ValueError("mode must be either 'cartesian' or 'fractional'")

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
    if element in (None, 0, "0", "all"):
        selected_indices = None
    else:
        element_key = str(element)
        if element_key not in atom_indices:
            available = ", ".join(sorted(atom_indices)) or "none"
            raise ValueError(
                f"Unknown element/reference label {element_key!r}; available labels: {available}"
            )
        selected_indices = set(atom_indices[element_key])

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
            positions.append(
                folded[0] * unit_vectors[0]
                + folded[1] * unit_vectors[1]
                + folded[2] * unit_vectors[2]
            )

    return RmcStructure(
        atom_indices=atom_indices,
        lattice_vectors=lattice_vectors,
        supercell=supercell,
        atom_types=atom_types,
        positions=np.asarray(positions, dtype=float),
    )
