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
    # None only for legacy coords-only lines (``iter_rmc6f_atoms(...,
    # include_coords_only=True)``); always set for the full layout.
    reference_number: int | None
    cell_indices: np.ndarray | None


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


_LINE_BREAK_RE = re.compile(r"\r\n|\r|\n")


@dataclass(frozen=True)
class ChiLog:
    """The chi^2 columns of one or more RMCProfile ``-NN.log`` files, concatenated.

    ``chi_r`` is the LAST log column and ``chi_q`` the second-to-last, one entry
    per complete data row (``NaN`` where the token is non-finite — ``NaN``,
    ``Inf``, a Fortran ``****`` overflow — or not a number). ``column`` is the
    header name of the last column (e.g. ``X_ray_(R)1``: the chi^2 of one fit
    term, not a total), ``None`` when the log has no column-name header.
    ``skipped_rows`` counts data lines dropped for a token count that differs
    from the header's, plus unterminated final lines.
    """

    chi_q: np.ndarray
    chi_r: np.ndarray
    column: str | None
    skipped_rows: int


def read_chi_log(paths: list[str | Path]) -> ChiLog:
    """Read RMCProfile ``.log`` chi^2 history, robust to a file still being written.

    Per file: line 1 names the columns (``Time moves_acc moves_gen F(Q)_1 …
    X_ray_(R)1``), line 2 carries the ``WEIGHT PARAMETERS`` and is skipped.
    A data row is kept only when it has exactly as many tokens as line 1 names
    (when line 1 names fewer than two columns — a log with no column-name header
    — the first data row sets the count), and a final line without a newline is
    dropped: in Live Data the log is re-read while RMCProfile appends to it, and
    a half-written last line would otherwise become the "final" chi^2 (a move
    counter, ``0.``, or a truncated mantissa). Rows are never dropped for their
    VALUE: a non-finite chi^2 stays in the series as ``NaN`` so a blown-up run
    shows as such. Mirrors ``readChi()`` in browserData.js.
    """
    chi_q: list[float] = []
    chi_r: list[float] = []
    column: str | None = None
    skipped = 0
    for path in paths:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
        lines = _LINE_BREAK_RE.split(text)
        # Every complete line ends with a break, so the last element is either ''
        # (a terminated file) or a line RMCProfile is still writing.
        unterminated = lines.pop()
        if unterminated.strip() and len(lines) >= 2:
            skipped += 1
        header = lines[0].split() if lines else []
        expected = len(header) if len(header) >= 2 else None
        if expected is not None and column is None:
            column = header[-1]
        for line in lines[2:]:
            parts = line.split()
            if not parts:
                continue
            if expected is None:
                expected = len(parts)
            if len(parts) != expected or len(parts) < 2:
                skipped += 1
                continue
            values = [parse_fortran_number(token) for token in parts[-2:]]
            chi_q.append(float("nan") if values[0] is None else values[0])
            chi_r.append(float("nan") if values[1] is None else values[1])
    return ChiLog(
        chi_q=np.asarray(chi_q, dtype=float),
        chi_r=np.asarray(chi_r, dtype=float),
        column=column,
        skipped_rows=skipped,
    )


def read_chi(paths: list[str | Path]) -> tuple[np.ndarray, np.ndarray]:
    """``(second-to-last, last)`` log columns — see :func:`read_chi_log`.

    The names are historical: in current RMCProfile logs the last column is the
    chi^2 of the last fitted term (``X_ray_(R)1`` in the demo run) and the
    second-to-last is often a constraint term, not a reciprocal-space chi^2.
    """
    log = read_chi_log(paths)
    return log.chi_q, log.chi_r


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


# Column-role vocabulary of RMCProfile fit CSV headers: ``F(Q)_Expt``,
# ``X_ray_exp_renorm``, ``observed`` name the measurement; ``F(Q)_RMC``,
# ``X_ray-calc``, ``calculated``, ``fitted`` the model curve. Mirrors
# ``EXPERIMENTAL_LABEL`` / ``CALCULATED_LABEL`` in browserData.js.
_EXPERIMENTAL_LABEL = re.compile(r"exp|obs", re.IGNORECASE)
_CALCULATED_LABEL = re.compile(r"calc|rmc|fit", re.IGNORECASE)


def rwp_columns(labels: list[str], n_columns: int | None = None) -> tuple[int, int] | None:
    """``(calculated, experimental)`` column indices for the R-factor of a fit CSV.

    RMCProfile writes its fit files as ``(x, calculated, experimental)`` —
    ``Q, F(Q)_RMC, F(Q)_Expt`` and ``r(A), X_ray-calc, X_ray_exp_renorm`` — so
    that positional layout is the default. When the header names both roles
    explicitly (a column matching ``exp``/``obs`` and another matching
    ``calc``/``rmc``/``fit``), the header wins, so a file written in another
    order is still normalized by its measurement. Returns ``None`` when there
    are fewer than three columns. Mirrors ``rwpColumns()`` in browserData.js.
    """
    count = len(labels) if n_columns is None else n_columns
    if count < 3:
        return None
    names = [str(label) for label in labels[1:count]]
    experimental = [idx for idx, name in enumerate(names, start=1) if _EXPERIMENTAL_LABEL.search(name)]
    calculated = [
        idx
        for idx, name in enumerate(names, start=1)
        if _CALCULATED_LABEL.search(name) and not _EXPERIMENTAL_LABEL.search(name)
    ]
    if experimental and calculated:
        return calculated[0], experimental[0]
    return 1, 2


def fit_rwp(labels: list[str], data: np.ndarray) -> float | None:
    """The dashboard R-factor of a parsed fit CSV: ``rwp`` normalized by the experiment.

    ``data`` is the transposed column array of :class:`CsvSeries`. The column
    roles come from :func:`rwp_columns`; ``None`` when fewer than three columns
    exist or when :func:`rwp` itself is undefined.
    """
    roles = rwp_columns(labels, len(data))
    if roles is None:
        return None
    calculated, experimental = roles
    return rwp(data[0], observed=data[experimental], fitted=data[calculated])


# --- .rmc6f atom-line grammar --------------------------------------------------
#
# One grammar, shared verbatim with ``web_app/frontend/src/rmc6f.js`` (keep the
# two in sync). An atom line is ``id element [label] <data>`` where
#
#   id       a non-negative integer (the atom number);
#   element  a token starting with a letter (normalized ``str.capitalize()``);
#   label    optional: a bracket group (``[1]``, or split as ``[ 1]``) or one
#            non-numeric token;
#   data     exactly 7 tokens  ``x y z ref cx cy cz``   (full layout), or
#            exactly 3 tokens  ``x y z``                (legacy coords-only).
#
# Numbers accept Fortran ``D`` exponents. ``ref`` must be a positive integer and
# each cell index an integer in [0, N_i) (N from the ``Supercell`` header line,
# when known). A line whose layout is valid but whose coordinates are non-finite
# (NaN, Inf, or Fortran ``****`` overflow) is skipped and counted separately;
# every other line after the ``Atoms`` marker that fits no layout is counted as
# unparsed. Nothing is guessed from the end of the line any more: an extra
# trailing field used to shift every column silently in the browser.

_RMC6F_NUMBER_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?$")
_RMC6F_NON_FINITE_RE = re.compile(r"^(?:[+-]?(?:nan|inf|infinity)|\*+)$", re.IGNORECASE)
_RMC6F_INTEGER_RE = re.compile(r"^[+-]?\d+$")
_RMC6F_ELEMENT_RE = re.compile(r"^[A-Za-z]")
_RMC6F_ATOMS_MARKER_RE = re.compile(r"^\s*atoms\b", re.IGNORECASE)
_RMC6F_DECLARED_ATOMS_RE = re.compile(r"^\s*Number of atoms\s*:\s*(\d+)", re.IGNORECASE)


def parse_fortran_number(token: str) -> float | None:
    """A numeric token of an RMCProfile text file as a float, Fortran-aware.

    Accepts plain/``E``/``D`` exponent forms (``0.117D-03``); returns ``nan`` for
    an explicit non-finite token (``NaN``, ``Inf``, ``Infinity``, or the all-``*``
    field Fortran prints on overflow) and ``None`` for anything that is not a
    number. Mirrors ``parseFortranNumber()`` in ``rmc6f.js``.
    """
    if _RMC6F_NUMBER_RE.match(token):
        return float(token.replace("D", "E").replace("d", "e"))
    if _RMC6F_NON_FINITE_RE.match(token):
        return float("nan")
    return None


def is_rmc6f_atoms_marker(line: str) -> bool:
    """True for the line that opens the atom list: ``Atoms:``, ``Atoms :``,
    ``atoms:``, ``Atoms (fractional coordinates):`` … (case-insensitive)."""
    return bool(_RMC6F_ATOMS_MARKER_RE.match(line))


@dataclass
class Rmc6fParseReport:
    """What an ``.rmc6f`` atom section held, line by line.

    ``declared_atoms`` is the header's ``Number of atoms:`` (``None`` if absent);
    ``atom_lines`` counts the non-blank lines after the ``Atoms`` marker. Each of
    them is exactly one of: a full-layout atom (``parsed_atoms``), a legacy
    coords-only atom (``coords_only_atoms``), a line skipped for non-finite
    coordinates (``non_finite_lines``) or an unparsed line (``invalid_lines``).
    Mirrors the ``report`` of ``parseRmc6fAtoms()`` in ``rmc6f.js``.
    """

    has_atoms_section: bool = False
    declared_atoms: int | None = None
    atom_lines: int = 0
    parsed_atoms: int = 0
    coords_only_atoms: int = 0
    non_finite_lines: int = 0
    invalid_lines: int = 0
    first_invalid_line: str | None = None
    first_non_finite_line: str | None = None

    @property
    def accepted_atoms(self) -> int:
        return self.parsed_atoms + self.coords_only_atoms

    def warning(self) -> str | None:
        """The human-readable problem list, or ``None`` for a clean atom section."""
        problems: list[str] = []
        if self.declared_atoms is not None and self.accepted_atoms != self.declared_atoms:
            problems.append(
                f"parsed {self.accepted_atoms} of {self.declared_atoms} atoms declared in the header"
            )
        if self.invalid_lines:
            problems.append(
                f"{self.invalid_lines} of {self.atom_lines} atom lines unparsed "
                f"(first: '{self.first_invalid_line}')"
            )
        if self.non_finite_lines:
            problems.append(
                f"{self.non_finite_lines} atom lines skipped for non-finite coordinates "
                f"(first: '{self.first_non_finite_line}')"
            )
        return "; ".join(problems) or None

    def to_dict(self) -> dict[str, object]:
        """camelCase mapping, the same keys the browser parser reports."""
        return {
            "declaredAtoms": self.declared_atoms,
            "atomLines": self.atom_lines,
            "parsedAtoms": self.parsed_atoms,
            "coordsOnlyAtoms": self.coords_only_atoms,
            "nonFiniteLines": self.non_finite_lines,
            "invalidLines": self.invalid_lines,
            "firstInvalidLine": self.first_invalid_line,
            "firstNonFiniteLine": self.first_non_finite_line,
        }


def _rmc6f_integer(token: str) -> int | None:
    return int(token) if _RMC6F_INTEGER_RE.match(token) else None


def classify_rmc6f_atom_line(
    parts: list[str],
    supercell: np.ndarray | None = None,
) -> tuple[str, Rmc6fAtom | None]:
    """Classify one whitespace-split atom line: ``(kind, atom)``.

    ``kind`` is ``"atom"`` (full layout), ``"coords"`` (legacy coords-only; the
    record's ``reference_number``/``cell_indices`` are ``None``), ``"non_finite"``
    (a valid layout with a NaN/Inf/``****`` coordinate; ``atom`` is ``None``) or
    ``"invalid"``. See the grammar comment above; mirrors ``classifyAtomLine()``
    in ``rmc6f.js``.
    """
    if len(parts) < 5:
        return "invalid", None
    atom_number = _rmc6f_integer(parts[0])
    if atom_number is None or atom_number < 0 or not _RMC6F_ELEMENT_RE.match(parts[1]):
        return "invalid", None
    element = parts[1].capitalize()

    index = 2
    label_tokens: list[str] = []
    if parts[index].startswith("["):
        while index < len(parts):
            label_tokens.append(parts[index])
            index += 1
            if label_tokens[-1].endswith("]"):
                break
        else:
            return "invalid", None
    elif parse_fortran_number(parts[index]) is None:
        label_tokens.append(parts[index])
        index += 1
    data = parts[index:]
    if len(data) not in (3, 7):
        return "invalid", None

    coords = [parse_fortran_number(token) for token in data[:3]]
    if any(value is None for value in coords):
        return "invalid", None
    record: Rmc6fAtom = {
        "atom_number": atom_number,
        "element": element,
        "type_label": " ".join(label_tokens),
        "coords": np.asarray(coords, dtype=float),
        "reference_number": None,
        "cell_indices": None,
    }
    if len(data) == 7:
        reference = _rmc6f_integer(data[3])
        cells = [_rmc6f_integer(token) for token in data[4:7]]
        if reference is None or reference < 1 or any(cell is None for cell in cells):
            return "invalid", None
        for axis, cell in enumerate(cells):
            limit = supercell[axis] if supercell is not None else None
            if cell < 0 or (limit is not None and np.isfinite(limit) and limit >= 1 and cell >= limit):
                return "invalid", None
        record["reference_number"] = reference
        record["cell_indices"] = np.asarray(cells, dtype=int)
    if not np.all(np.isfinite(record["coords"])):
        return "non_finite", None
    return ("atom" if len(data) == 7 else "coords"), record


def read_atom_indices(rmc6f_path: str | Path) -> dict[str, list[int]]:
    """Distinct reference numbers per element over the full-layout atom lines.

    Built from :func:`iter_rmc6f_atoms` so the site table and the atom list can
    never disagree (it used to read ``parts[-4]`` of any line, and reported cell
    indices as "sites" when a line carried an extra field).
    """
    report = Rmc6fParseReport()
    atom_indices: dict[str, set[int]] = {}
    for atom in iter_rmc6f_atoms(rmc6f_path, report=report):
        atom_indices.setdefault(atom["element"], set()).add(int(atom["reference_number"]))
    if not report.has_atoms_section:
        raise ValueError(f"{rmc6f_path} does not contain an Atoms section")
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
    marker = re.search(r"^[ \t]*atoms\b", text, re.IGNORECASE | re.MULTILINE)
    header = text[: marker.start() if marker and marker.start() > 0 else 4000]

    moves: dict[str, float] = {}
    for key, pattern in _MOVE_COUNTERS.items():
        match = pattern.search(header)
        if match:
            moves[key] = float(match.group(1))
    return moves or None


def iter_rmc6f_atoms(
    rmc6f_path: str | Path,
    *,
    include_coords_only: bool = False,
    report: Rmc6fParseReport | None = None,
) -> Iterator[Rmc6fAtom]:
    """Yield atom records from an RMCProfile ``.rmc6f`` file.

    Lines follow the grammar documented above :func:`classify_rmc6f_atom_line`
    (identical to the browser parser). By default only full-layout atoms are
    yielded, so ``reference_number`` and ``cell_indices`` are always set;
    ``include_coords_only=True`` also yields legacy ``id element [label] x y z``
    atoms, with those two fields ``None`` (enough for position-only consumers such
    as bond angles). Lines with non-finite coordinates are never yielded. Pass a
    :class:`Rmc6fParseReport` as ``report`` to learn what was skipped and whether
    the count matches the header's ``Number of atoms:``.
    """
    if report is None:
        report = Rmc6fParseReport()
    supercell: np.ndarray | None = None
    in_atoms = False
    with Path(rmc6f_path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.split()
            if not parts:
                continue
            if not in_atoms:
                if is_rmc6f_atoms_marker(line):
                    in_atoms = True
                    report.has_atoms_section = True
                    continue
                declared = _RMC6F_DECLARED_ATOMS_RE.match(line)
                if declared:
                    report.declared_atoms = int(declared.group(1))
                elif parts[0] == "Supercell" and len(parts) >= 3:
                    values = [parse_fortran_number(token) for token in parts[-3:]]
                    if all(value is not None for value in values):
                        supercell = np.asarray(values, dtype=float)
                continue
            report.atom_lines += 1
            kind, atom = classify_rmc6f_atom_line(parts, supercell)
            if kind == "atom":
                report.parsed_atoms += 1
                yield atom
            elif kind == "coords":
                report.coords_only_atoms += 1
                if include_coords_only:
                    yield atom
            elif kind == "non_finite":
                report.non_finite_lines += 1
                if report.first_non_finite_line is None:
                    report.first_non_finite_line = line.strip()
            else:
                report.invalid_lines += 1
                if report.first_invalid_line is None:
                    report.first_invalid_line = line.strip()


def parse_rmc6f_atoms(
    rmc6f_path: str | Path,
    *,
    include_coords_only: bool = True,
) -> tuple[list[Rmc6fAtom], Rmc6fParseReport]:
    """All atoms of an ``.rmc6f`` file plus the :class:`Rmc6fParseReport`.

    Raises ``ValueError`` when the file has no ``Atoms`` section, or when not a
    single atom line could be parsed (naming what was found instead of letting a
    caller report an empty model).
    """
    report = Rmc6fParseReport()
    atoms = list(iter_rmc6f_atoms(rmc6f_path, include_coords_only=include_coords_only, report=report))
    if not report.has_atoms_section:
        raise ValueError(f"{rmc6f_path} does not contain an Atoms section")
    if not atoms:
        detail = report.warning() or "the Atoms section is empty"
        raise ValueError(f"{rmc6f_path}: no atoms could be parsed — {detail}")
    return atoms, report


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


_RMC6F_HEAD_BYTES = 65536
_FRAC_STEM_RE = re.compile(r"^Frac_coord_(.+)\.txt$")


def rmc6f_problem(path: str | Path) -> str | None:
    """Why ``path`` cannot be a run's configuration, or ``None`` when it can.

    A killed run can leave a 0-byte (or header-only) ``.rmc6f`` beside valid
    ones; picking it hid every usable model in the folder. A candidate must be
    non-empty and show the ``Atoms`` marker (see :func:`is_rmc6f_atoms_marker`)
    within its first 64 KiB. Mirrors ``structureFileProblem()`` in browserData.js.
    """
    path = Path(path)
    try:
        size = path.stat().st_size
    except OSError as exc:
        return f"unreadable ({exc.strerror or exc})"
    if size == 0:
        return "empty (0 bytes)"
    with path.open("rb") as handle:
        head = handle.read(_RMC6F_HEAD_BYTES).decode("utf-8", errors="replace")
    if not any(is_rmc6f_atoms_marker(line) for line in head.splitlines()):
        return "no Atoms section in its first 64 KiB"
    return None


def _structure_pair(directory: Path, frac_path, rmc6f_path) -> tuple[Path, Path]:
    """The (Frac*.txt, .rmc6f) pair of ONE configuration for :func:`read_structure`."""
    frac_files = sorted(directory.glob("Frac*.txt"))
    rmc6f_files = sorted(directory.glob("*.rmc6f"))
    usable = [path for path in rmc6f_files if rmc6f_problem(path) is None]

    def rmc6f_for(frac: Path) -> Path | None:
        match = _FRAC_STEM_RE.match(frac.name)
        candidate = frac.with_name(f"{match.group(1)}.rmc6f") if match else None
        return candidate if candidate is not None and candidate in usable else None

    def frac_for(rmc6f: Path) -> Path | None:
        candidate = rmc6f.with_name(f"Frac_coord_{rmc6f.stem}.txt")
        return candidate if candidate.exists() else None

    if frac_path is not None and rmc6f_path is not None:
        return Path(frac_path), Path(rmc6f_path)
    if rmc6f_path is not None:
        rmc6f_path = Path(rmc6f_path)
        frac = frac_for(rmc6f_path) or (frac_files[0] if len(frac_files) == 1 else None)
        if frac is None:
            raise FileNotFoundError(
                f"No Frac_coord_{rmc6f_path.stem}.txt beside {rmc6f_path}; pass frac_path="
            )
        return frac, rmc6f_path
    if frac_path is not None:
        frac_path = Path(frac_path)
        rmc6f = rmc6f_for(frac_path) or (usable[0] if len(usable) == 1 else None)
        if rmc6f is None:
            raise FileNotFoundError(f"No usable .rmc6f pairs with {frac_path.name}; pass rmc6f_path=")
        return frac_path, rmc6f

    if not frac_files:
        raise FileNotFoundError(f"No Frac*.txt file found in {directory}")
    if not rmc6f_files:
        raise FileNotFoundError(f"No .rmc6f file found in {directory}")
    for frac in frac_files:
        rmc6f = rmc6f_for(frac)
        if rmc6f is not None:
            return frac, rmc6f
    if len(frac_files) == 1 and len(usable) == 1:
        # A single-configuration folder whose files do not share a stem.
        return frac_files[0], usable[0]
    skipped = [f"{path.name} ({rmc6f_problem(path)})" for path in rmc6f_files if path not in usable]
    raise ValueError(
        f"Cannot tell which configuration the Frac*.txt files in {directory} belong to: "
        f"no Frac_coord_<stem>.txt pairs with a usable <stem>.rmc6f "
        f"(Frac: {', '.join(path.name for path in frac_files)}; "
        f".rmc6f: {', '.join(path.name for path in usable) or 'none usable'}"
        f"{'; skipped ' + ', '.join(skipped) if skipped else ''}). "
        "Pass frac_path= and rmc6f_path= explicitly."
    )


def read_structure(
    directory: str | Path,
    element: str | int | None = None,
    mode: str = "cartesian",
    *,
    frac_path: str | Path | None = None,
    rmc6f_path: str | Path | None = None,
) -> RmcStructure:
    """Folded unit-cell positions from a ``Frac_coord_<stem>.txt`` and its ``<stem>.rmc6f``.

    The two files must describe the SAME configuration — the ``.rmc6f`` supplies
    the supercell used for the fold and the element → reference-number map — so
    they are paired by stem (``Frac_coord_<stem>.txt`` ⟷ ``<stem>.rmc6f``),
    skipping empty or marker-less ``.rmc6f`` candidates. A folder with exactly
    one Frac file and one usable ``.rmc6f`` pairs them regardless of name; any
    other ambiguity raises instead of pairing files from different runs. Pass
    ``frac_path`` and/or ``rmc6f_path`` to choose explicitly. The pair is
    cross-checked: every Frac cell index must lie inside the ``.rmc6f``
    supercell and every Frac reference number must be one of its sites.
    """
    if mode not in {"cartesian", "fractional"}:
        raise ValueError("mode must be either 'cartesian' or 'fractional'")

    directory = Path(directory)
    frac_path, rmc6f_path = _structure_pair(directory, frac_path, rmc6f_path)

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

    known_references = {index for indices in atom_indices.values() for index in indices}
    max_cells = np.full(3, -1)
    frac_references: set[int] = set()
    atom_types: list[str] = []
    positions: list[np.ndarray] = []
    with frac_path.open("r", encoding="utf-8") as handle:
        lines = handle.readlines()[5:]

    for line in lines:
        parts = line.split()
        if len(parts) < 4:
            continue
        atom_id = int(parts[0])
        frac_references.add(atom_id)
        if len(parts) >= 7:
            max_cells = np.maximum(max_cells, [int(value) for value in parts[4:7]])
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

    outside = np.nonzero(max_cells >= supercell)[0]
    unknown = sorted(frac_references - known_references)
    if outside.size or unknown:
        details = []
        if outside.size:
            details.append(
                f"cell indices up to {max_cells.tolist()} exceed the supercell {supercell.tolist()}"
            )
        if unknown:
            details.append(f"reference numbers {unknown[:5]}{'…' if len(unknown) > 5 else ''} are not sites of it")
        raise ValueError(
            f"{frac_path.name} does not belong to {rmc6f_path.name}: {'; '.join(details)}. "
            "Pass the matching frac_path= / rmc6f_path=."
        )

    return RmcStructure(
        atom_indices=atom_indices,
        lattice_vectors=lattice_vectors,
        supercell=supercell,
        atom_types=atom_types,
        positions=np.asarray(positions, dtype=float),
    )
