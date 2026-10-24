# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Polyhedral distortion: how far each coordination polyhedron is from regular.

For a centre element C and a ligand element L, the L atoms within an inclusive
bond window [r_min, r_max] of a C atom are the vertices of its coordination
polyhedron. The neighbours come from the triplets engine's search
(``triplets._neighbor_bonds``: the linked-cell periodic search with explicit
image shifts, triclinic-safe, the same ``WINDOW_TOL`` on both bounds), never a
second one. The polyhedra measured are those of the most common coordination n
among the centres with at least one ligand (the lowest count on a tie); a
centre with any other count is left out and counted (``skipped``, and the full
``coordination`` histogram).

Three standard indices, each 0, 1 and 0 for a regular polyhedron:

* Baur's bond-length distortion index (W. H. Baur, Acta Cryst. B30 (1974)
  1195-1215), for any n::

      DI = (1/n) sum_i |l_i - <l>| / <l>,   <l> = (1/n) sum_i l_i

* Robinson's quadratic elongation (K. Robinson, G. V. Gibbs and P. H. Ribbe,
  Science 172 (1971) 567-570), for octahedra (n = 6) and tetrahedra (n = 4)::

      <lambda> = (1/n) sum_i (l_i / l_0)^2

  with l_0 the centre-to-vertex distance of the regular polyhedron of the same
  volume V: V = (4/3) l_0^3 for the octahedron, V = (8 sqrt(3) / 27) l_0^3 for
  the tetrahedron. <lambda> >= 1, and 1 only for the regular polyhedron with
  the centre at its centre: an off-centred atom in a regular octahedron gives
  1 + (delta / l_0)^2.

* Robinson's bond-angle variance (the same paper), in degrees squared::

      sigma^2_theta = sum_{i=1}^{m} (theta_i - theta_0)^2 / (m - 1)

  over the polyhedron's m cis angles L-C-L: the angles the polyhedron's edges
  subtend at the centre. For a tetrahedron that is all 6 angles, theta_0 =
  arccos(-1/3) = 109.4712 deg (Robinson's 109.47 deg; the exact value, so a
  regular tetrahedron gives 0). For an octahedron it is the 12 angles that are
  90 deg when it is regular (theta_0 = 90 deg): its 15 angles less the three
  trans ones.

The three trans pairs of an octahedron are the perfect matching of its six
vertices (one of the 15) whose three angles have the largest sum; a matching
within ``TRANS_TIE_DEG`` of that sum goes to the first in a fixed order
(``OCTAHEDRON_MATCHINGS``). The octahedron's volume is that of the polyhedron
bounded by the 8 faces holding one vertex of each trans pair: summed over the
8 signed face tetrahedra, multilinearity of the determinant gives::

      V = |d_1 . (d_2 x d_3)| / 6,   d_k = v_k - v_k'  (the k-th trans pair)

and the tetrahedron's is |(v_2 - v_1) . ((v_3 - v_1) x (v_4 - v_1))| / 6. A
polyhedron of zero volume (flat) has no <lambda>: it is left out of that
summary, whose ``count`` then falls short of ``polyhedra``. The mean of each
octahedron's three trans angles is reported too (``transAngle``: 180 deg for
a regular octahedron), so a polyhedron that is not octahedral shows.

Each statistic is summarised over the polyhedra as ``{count, mean, sd, p10,
median, p90}`` (sd the sample sd, null for one polyhedron; quantiles by linear
interpolation between order statistics).

The same indices are measured on the AVERAGE structure, so a static distortion
separates from the spread thermal motion adds: the average-structure basis
(``structure_basis``: one site per reference number at its mean position
``mean`` in the unit cell, lattice / supercell), each centre site's polyhedron
from the same window and the same search. A site counts as the species that
holds most of its copies (ties to the code-point first symbol, the PCA
engine's site label; the basis's ``el``, its first atom's species, is what the
symmetry finder reads). Every index is non-negative and is 0 (``<lambda>``: 1)
only for a regular polyhedron, so thermal motion and any disorder the average
folds away raise the configuration's values but not the average's: a
distortion the average polyhedron carries too is static, one only the
configuration's polyhedra show is local.

The browser port is ``workers/distortion.js``; keep the two in sync (the parity
fixture ``tests/generate_distortion_fixture.py`` pins them).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Sequence

import numpy as np

from .parsers import Rmc6fAtom, parse_rmc6f_atoms, read_cell_vectors
from .structure_basis import structure_basis
from .triplets import WINDOW_TOL, _Bonds, _neighbor_bonds, _sort_by_center, coordination_summary

#: Work budget of one app request: the bonds the window holds, counted by a
#: search that stores nothing before any is stored. The Flask route and the
#: browser worker refuse a request above it (the engine itself is unrestricted,
#: ``max_bonds=None``). Mirrored as APP_MAX_BONDS in workers/distortion.js.
APP_MAX_BONDS = 1_000_000

#: Robinson's ideal angles, degrees: an octahedron's cis angle and the
#: tetrahedral angle arccos(-1/3).
OCTAHEDRAL_ANGLE_DEG = 90.0
TETRAHEDRAL_ANGLE_DEG = math.degrees(math.acos(-1.0 / 3.0))

#: l_0^3 per unit volume of the regular polyhedron: V = (4/3) l_0^3 for the
#: octahedron, V = (8 sqrt(3) / 27) l_0^3 for the tetrahedron.
OCTAHEDRON_RADIUS_CUBED_PER_VOLUME = 3.0 / 4.0
TETRAHEDRON_RADIUS_CUBED_PER_VOLUME = 9.0 * math.sqrt(3.0) / 8.0

#: Two trans matchings whose angle sums lie within this many degrees are a tie,
#: which goes to the first in OCTAHEDRON_MATCHINGS: far above the float noise of
#: a sum of angles (~1e-13 deg), far below any real difference.
TRANS_TIE_DEG = 1e-9

#: The quantiles each summary gives (p10, median, p90).
QUANTILES = (0.1, 0.5, 0.9)

#: The ideal polyhedra: coordination -> (shape, ideal cis angle).
IDEAL_POLYHEDRA = {6: ("octahedron", OCTAHEDRAL_ANGLE_DEG), 4: ("tetrahedron", TETRAHEDRAL_ANGLE_DEG)}

#: Vertex pairs (i < j) of a polyhedron of n vertices, in this order.
PAIRS = {n: tuple(combinations(range(n), 2)) for n in IDEAL_POLYHEDRA}


def _perfect_matchings(items: list[int]) -> list[list[tuple[int, int]]]:
    """Every perfect matching of ``items``: the first pairs with each later one
    in turn, the rest recursively (workers/distortion.js builds the same order)."""
    if not items:
        return [[]]
    first, rest = items[0], items[1:]
    matchings = []
    for k, partner in enumerate(rest):
        for tail in _perfect_matchings(rest[:k] + rest[k + 1:]):
            matchings.append([(first, partner)] + tail)
    return matchings


#: The 15 ways to choose an octahedron's three trans pairs, as indices into
#: PAIRS[6], in the order a tie is broken.
OCTAHEDRON_MATCHINGS = tuple(
    tuple(PAIRS[6].index(pair) for pair in matching) for matching in _perfect_matchings(list(range(6)))
)


@dataclass(frozen=True)
class DistortionConfiguration:
    """What the engine reads of a configuration.

    ``fractional`` (N, 3) supercell fractions of every atom (coords-only lines
    included, as the bond-angle engine reads them) with their ``elements``;
    ``lattice_vectors`` the (3, 3) supercell rows in angstrom; ``supercell``
    its tiling; ``basis`` the average-structure basis (``structure_basis``
    rows; empty for a coords-only file); ``parse_warning`` what the ``.rmc6f``
    grammar skipped, or None.
    """

    fractional: np.ndarray
    elements: tuple[str, ...]
    lattice_vectors: np.ndarray
    supercell: np.ndarray
    basis: tuple[dict, ...]
    parse_warning: str | None = None


def distortion_configuration(
    atoms: Sequence[Rmc6fAtom], lattice_vectors, supercell, parse_warning: str | None = None
) -> DistortionConfiguration:
    """The engine's input from parsed atoms (``parse_rmc6f_atoms``, both layouts)."""
    atoms = list(atoms)
    return DistortionConfiguration(
        fractional=np.array([atom["coords"] for atom in atoms], dtype=float).reshape(-1, 3),
        elements=tuple(str(atom["element"]) for atom in atoms),
        lattice_vectors=np.asarray(lattice_vectors, dtype=float).reshape(3, 3),
        supercell=np.asarray(supercell, dtype=float).reshape(3),
        basis=tuple(structure_basis(atoms, lattice_vectors, supercell)),
        parse_warning=parse_warning,
    )


def read_distortion_configuration(path: str | Path) -> DistortionConfiguration:
    """:func:`distortion_configuration` of an ``.rmc6f`` file (one read)."""
    lattice_vectors, supercell = read_cell_vectors(path)
    atoms, report = parse_rmc6f_atoms(path, include_coords_only=True)
    return distortion_configuration(atoms, lattice_vectors, supercell, report.warning())


def _symbol(value) -> str:
    return str(value).strip().capitalize()


def _validate_window(r_min, r_max) -> tuple[float, float]:
    try:
        lo, hi = float(r_min), float(r_max)
    except (TypeError, ValueError) as error:
        raise ValueError(f"window must be two numbers 0 <= r_min < r_max, got [{r_min}, {r_max}]") from error
    if not (math.isfinite(lo) and math.isfinite(hi)) or lo < 0 or hi <= lo:
        raise ValueError(f"window must be 0 <= r_min < r_max, got [{r_min}, {r_max}]")
    return lo, hi


def site_species(element_counts: dict) -> str:
    """The species holding most of a site's copies, ties to the code-point first symbol."""
    return min(element_counts.items(), key=lambda item: (-item[1], item[0]))[0]


def _quantile(ordered: np.ndarray, q: float) -> float:
    """Linear interpolation between the order statistics of ``ordered`` (sorted)."""
    position = (ordered.size - 1) * q
    lower = int(math.floor(position))
    upper = min(lower + 1, ordered.size - 1)
    return float(ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower))


def summarize(values) -> dict | None:
    """``{count, mean, sd, p10, median, p90}`` of the finite values; None when there is none."""
    values = np.asarray(values, dtype=float).reshape(-1)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return None
    ordered = np.sort(values)
    mean = float(values.sum() / values.size)
    deviation = values - mean
    return {
        "count": int(values.size),
        "mean": mean,
        "sd": float(math.sqrt(float((deviation * deviation).sum()) / (values.size - 1))) if values.size > 1 else None,
        "p10": _quantile(ordered, QUANTILES[0]),
        "median": _quantile(ordered, QUANTILES[1]),
        "p90": _quantile(ordered, QUANTILES[2]),
    }


def _triple(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    """a . (b x c), rowwise, written out in the JS port's evaluation order."""
    cross_x = b[:, 1] * c[:, 2] - b[:, 2] * c[:, 1]
    cross_y = b[:, 2] * c[:, 0] - b[:, 0] * c[:, 2]
    cross_z = b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0]
    return (a[:, 0] * cross_x + a[:, 1] * cross_y) + a[:, 2] * cross_z


def polyhedron_indices(vectors: np.ndarray, lengths: np.ndarray) -> dict:
    """DI, <lambda>, sigma^2_theta and the mean trans angle of polyhedra of one coordination.

    ``vectors`` (M, n, 3) are the centre-to-ligand vectors in angstrom and
    ``lengths`` (M, n) their lengths. Returns arrays of M values:
    ``distortionIndex`` always; ``quadraticElongation`` (NaN for a flat
    polyhedron) and ``angleVariance`` only for n = 4 or 6, and ``transAngle``
    only for n = 6 (None otherwise).
    """
    vectors = np.asarray(vectors, dtype=float)
    lengths = np.asarray(lengths, dtype=float)
    count, n = lengths.shape
    mean_length = lengths.sum(axis=1) / n
    result = {
        "distortionIndex": np.abs(lengths - mean_length[:, None]).sum(axis=1) / n / mean_length,
        "quadraticElongation": None,
        "angleVariance": None,
        "transAngle": None,
    }
    if n not in IDEAL_POLYHEDRA or count == 0:
        return result
    pairs = np.asarray(PAIRS[n])
    first, second = pairs[:, 0], pairs[:, 1]
    one, two = vectors[:, first, :], vectors[:, second, :]
    # The bond-angle engine's evaluation order (bitwise the JS port's).
    dot = (one[..., 0] * two[..., 0] + one[..., 1] * two[..., 1]) + one[..., 2] * two[..., 2]
    cosine = dot / (lengths[:, first] * lengths[:, second])
    angles = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))  # (M, pairs)
    rows = np.arange(count)
    cis = np.ones(angles.shape, dtype=bool)
    if n == 6:
        matchings = np.asarray(OCTAHEDRON_MATCHINGS)
        sums = (angles[:, matchings[:, 0]] + angles[:, matchings[:, 1]]) + angles[:, matchings[:, 2]]
        best = sums.max(axis=1)
        trans = matchings[np.argmax(sums >= (best - TRANS_TIE_DEG)[:, None], axis=1)]  # (M, 3) pair indices
        cis[rows[:, None], trans] = False
        trans_angles = angles[rows[:, None], trans]
        result["transAngle"] = ((trans_angles[:, 0] + trans_angles[:, 1]) + trans_angles[:, 2]) / 3.0
        tips = pairs[trans]  # (M, 3, 2) vertex indices of the trans pairs
        diagonals = vectors[rows[:, None], tips[..., 0], :] - vectors[rows[:, None], tips[..., 1], :]
        volume = np.abs(_triple(diagonals[:, 0], diagonals[:, 1], diagonals[:, 2])) / 6.0
        radius_cubed = OCTAHEDRON_RADIUS_CUBED_PER_VOLUME * volume
    else:
        edges = vectors[:, 1:, :] - vectors[:, :1, :]
        volume = np.abs(_triple(edges[:, 0], edges[:, 1], edges[:, 2])) / 6.0
        radius_cubed = TETRAHEDRON_RADIUS_CUBED_PER_VOLUME * volume
    deviation = angles - IDEAL_POLYHEDRA[n][1]
    squared = np.where(cis, deviation * deviation, 0.0)
    result["angleVariance"] = squared.sum(axis=1) / (int(cis[0].sum()) - 1)
    radius = np.cbrt(radius_cubed)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        ratio = lengths / radius[:, None]
        elongation = (ratio * ratio).sum(axis=1) / n
    result["quadraticElongation"] = np.where((radius > 0) & np.isfinite(elongation), elongation, np.nan)
    return result


def _gather(grouped: _Bonds, starts: np.ndarray, centres: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    """(M, n, 3) vectors and (M, n) lengths of the given centres, each with n
    bonds (``grouped``: the bonds sorted by centre, the search's order within a
    centre; ``starts``: each centre's first bond)."""
    index = starts[np.asarray(centres, dtype=int)][:, None] + np.arange(n)[None, :]
    return grouped.vectors[index], grouped.lengths[index]


def _block(vectors: np.ndarray, lengths: np.ndarray) -> dict:
    """The summaries of polyhedra of one coordination (``polyhedra`` of them)."""
    count = int(lengths.shape[0])
    if count == 0:
        return {"polyhedra": 0, "bondLength": None, "distortionIndex": None,
                "quadraticElongation": None, "angleVariance": None, "transAngle": None}
    indices = polyhedron_indices(vectors, lengths)
    return {
        "polyhedra": count,
        "bondLength": summarize(lengths),
        **{key: (None if value is None else summarize(value)) for key, value in indices.items()},
    }


def _finite_or_none(value) -> float | None:
    return float(value) if value is not None and math.isfinite(value) else None


def _average_structure(config: DistortionConfiguration, centre: str, ligand: str, window, n: int) -> dict:
    """The average structure's polyhedra: per centre site, and summarised over
    the sites with the configuration's coordination n (``reason`` says why
    there is no summary)."""
    empty = {**_block(np.empty((0, n, 3)), np.empty((0, n))), "reason": None, "sites": []}
    basis = list(config.basis)
    if not basis:
        return {**empty, "reason": "no_basis"}
    labels = [site_species(row["elementCounts"]) for row in basis]
    centre_sites = np.asarray([k for k, label in enumerate(labels) if label == centre], dtype=int)
    ligand_sites = np.asarray([k for k, label in enumerate(labels) if label == ligand], dtype=int)
    if centre_sites.size == 0:
        return {**empty, "reason": "no_centre_sites"}
    cells = np.asarray(config.supercell, dtype=float).reshape(3)
    unit = np.asarray(config.lattice_vectors, dtype=float).reshape(3, 3) / np.maximum(cells, 1.0)[:, None]
    positions = np.asarray([row["mean"] for row in basis], dtype=float).reshape(-1, 3)
    wrapped = positions - np.floor(positions)
    bonds = _neighbor_bonds(wrapped[centre_sites], centre_sites, wrapped[ligand_sites], ligand_sites, unit, [window])
    grouped, sizes, starts = _sort_by_center(bonds, centre_sites.size)

    sites = [
        {
            "referenceNumber": int(basis[site]["referenceNumber"]),
            "element": labels[site],
            "mixed": len(basis[site]["elementCounts"]) > 1,
            "copies": int(basis[site]["count"]),
            "coordination": int(sizes[local]),
            "bondLengths": [],
            "distortionIndex": None,
            "quadraticElongation": None,
            "angleVariance": None,
            "transAngle": None,
        }
        for local, site in enumerate(centre_sites.tolist())
    ]
    # Each site's own polyhedron, whatever its coordination: one pass per count.
    for size in sorted(set(sizes.tolist()) - {0}):
        members = np.flatnonzero(sizes == size)
        vectors, lengths = _gather(grouped, starts, members, size)
        indices = polyhedron_indices(vectors, lengths)
        for row, local in enumerate(members.tolist()):
            sites[local]["bondLengths"] = sorted(float(value) for value in lengths[row])
            for key, values in indices.items():
                sites[local][key] = None if values is None else _finite_or_none(values[row])
    matching = np.flatnonzero(sizes == n)
    block = _block(*_gather(grouped, starts, matching, n))
    return {**block, "reason": None if matching.size else "coordination_differs", "sites": sites}


def polyhedral_distortion(
    config: DistortionConfiguration,
    centre: str,
    ligand: str,
    r_min: float,
    r_max: float,
    *,
    max_bonds: int | None = None,
) -> dict:
    """Baur's DI, Robinson's <lambda> and sigma^2_theta of the C-centred L polyhedra.

    ``centre`` and ``ligand`` are element symbols (any case); ``r_min`` and
    ``r_max`` the inclusive bond window in angstrom. Raises ``ValueError`` for
    a malformed window, an element the configuration lacks, no polyhedra (no
    centre has a ligand in the window, or the most common count is one), and,
    with ``max_bonds`` set (the app boundaries pass ``APP_MAX_BONDS``), a
    window that holds more bonds than that. Returns the JSON-safe payload of
    ``/api/distortion`` and the browser worker's ``distortion`` request (see the
    module doc); ``workers/distortion.js`` is its port.
    """
    lo, hi = _validate_window(r_min, r_max)
    centre, ligand = _symbol(centre), _symbol(ligand)
    fractional = np.asarray(config.fractional, dtype=float).reshape(-1, 3)
    if not np.all(np.isfinite(fractional)):
        raise ValueError("fractional coordinates contain non-finite values")
    lattice = np.asarray(config.lattice_vectors, dtype=float).reshape(3, 3)
    symbols = np.asarray([_symbol(element) for element in config.elements])
    if symbols.size != fractional.shape[0]:
        raise ValueError(f"{symbols.size} elements for {fractional.shape[0]} coordinates")
    available = ", ".join(sorted(set(symbols.tolist())))
    for role, symbol in (("centre", centre), ("ligand", ligand)):
        if not np.any(symbols == symbol):
            raise ValueError(f"unknown element {symbol!r} for the {role}: the configuration holds {available}")

    # Fold every coordinate into [0, 1); the image bookkeeping restores the
    # true relative geometry (as in the bond-angle engine).
    wrapped = fractional - np.floor(fractional)
    centre_rows = np.flatnonzero(symbols == centre)
    ligand_rows = np.flatnonzero(symbols == ligand)

    def search(count_only: bool = False):
        return _neighbor_bonds(
            wrapped[centre_rows], centre_rows, wrapped[ligand_rows], ligand_rows, lattice, [(lo, hi)],
            count_only=count_only,
        )

    if max_bonds is not None:
        counted = int(search(count_only=True)[:, 1].sum())
        if counted > max_bonds:
            raise ValueError(
                f"{centre}-{ligand} within {lo:g}-{hi:g} A holds {counted:,} bonds, over the limit of "
                f"{int(max_bonds):,} for one request; narrow the window"
            )

    n_centres = int(centre_rows.size)
    grouped, sizes, starts = _sort_by_center(search(), n_centres)
    coordination = np.bincount(sizes, minlength=1)
    bonded = coordination[1:]
    if not bonded.any():
        raise ValueError(f"no polyhedra: no {centre} atom has a {ligand} neighbour within [{lo:g}, {hi:g}] A")
    n = 1 + int(np.argmax(bonded))  # the most common count among bonded centres, the lowest on a tie
    if n < 2:
        raise ValueError(
            f"no polyhedra: the most common {centre} coordination within [{lo:g}, {hi:g}] A is one {ligand} "
            "(a polyhedron needs at least 2 ligands)"
        )
    chosen = np.flatnonzero(sizes == n)
    shape, ideal = IDEAL_POLYHEDRA.get(n, (None, None))
    return {
        "centre": centre,
        "ligand": ligand,
        "window": [lo, hi],
        "centres": n_centres,
        # coordination[k] = centre atoms with exactly k ligands in the window.
        "coordination": coordination.tolist(),
        "coordinationSummary": coordination_summary(coordination),
        "polyhedron": {
            "coordination": n,
            "shape": shape,
            "idealAngle": ideal,
            "cisAngles": None if shape is None else len(PAIRS[n]) - (3 if n == 6 else 0),
        },
        "skipped": n_centres - int(chosen.size),
        "configuration": _block(*_gather(grouped, starts, chosen, n)),
        "average": _average_structure(config, centre, ligand, (lo, hi), n),
        "rules": {
            "windowTol": WINDOW_TOL,
            "transTieDeg": TRANS_TIE_DEG,
            "tetrahedralAngle": TETRAHEDRAL_ANGLE_DEG,
            "quantiles": list(QUANTILES),
        },
    }


def polyhedral_distortion_from_file(
    path: str | Path,
    centre: str,
    ligand: str,
    r_min: float,
    r_max: float,
    max_bonds: int | None = None,
) -> dict:
    """:func:`polyhedral_distortion` of an ``.rmc6f`` file, read afresh, plus
    ``parseWarning`` (what the shared grammar skipped, or None). Callers that
    cache it key the cache on the file's content signature, as ``app.py``'s
    ``_DISTORTION_CACHE`` does."""
    config = read_distortion_configuration(path)
    result = polyhedral_distortion(config, centre, ligand, r_min, r_max, max_bonds=max_bonds)
    result["parseWarning"] = config.parse_warning
    return result
