# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Bond-angle (triplet) distribution engine for RMC configurations.

Given a triplet of atom types ``(A, B, C)`` -- with **B the central atom**, as
in RMCProfile's ``triplets`` utility -- and a bond-length window for each of
the A--B and B--C bonds, this module finds every A--B--C triplet in a periodic
configuration and histograms the angle at B. The workflow mirrors
``triplets_new_bonds_sinth`` from the RMCProfile tool set: select three atom
types, bound the two bond lengths, and read off the angle distribution with an
optional sin(theta) geometric correction.

Bond search
-----------
An ``.rmc6f`` configuration is periodic in the supercell, so neighbours are
found with a linked-cell search that carries *explicit* periodic-image shifts
instead of assuming the minimum-image convention. The box is divided into
``n_i`` cells per lattice direction, sized so each cell's perpendicular
thickness is at least ``rmax`` where possible, and every candidate cell within
``k_i`` layers of a central atom's cell is visited; an out-of-range cell index
wraps and records the whole-box shift it wrapped by. Because atoms in cells
``q`` layers apart are strictly more than ``(q - 1)`` cell thicknesses apart,
``k_i = ceil(rmax * n_i / w_i)`` layers (``w_i`` = perpendicular width of the
box along direction ``i``) cover every pair within ``rmax`` -- for any cell
shape, triclinic included, and even when the box is *smaller* than ``rmax``,
in which case multiple images of the same atom are genuine distinct
neighbours. Bond windows are inclusive at both ends: ``rmin <= r <= rmax`` --
except that a pair at exactly zero length (bitwise-coincident atoms under
``rmin = 0``) is never a bond, since a zero vector subtends no angle.

Angle counting
--------------
For each central B atom the A-bond list and the C-bond list are combined:

- If the two ends are equivalent (same element, same window) each *unordered*
  pair of distinct bonds contributes one angle, so an octahedron's six bonds
  give the expected ``C(6, 2) = 15`` angles (12 x 90 deg + 3 x 180 deg).
- Otherwise every (A-bond, C-bond) combination contributes one angle, minus
  the combinations where both bonds reach the *same atom in the same periodic
  image* -- the degenerate zero-degree "angle" of a bond with itself, which
  arises when A and C name the same element with overlapping windows.

Histogram conventions
---------------------
Angles are binned uniformly over [0, 180] degrees. Alongside the raw counts
the result carries:

- ``density`` -- counts / (total * bin_width), a probability density per
  degree with unit integral over [0, 180].
- ``sin_corrected`` -- the count fraction divided by the *exact* isotropic
  reference fraction per bin, ``(cos(edge_lo) - cos(edge_hi)) / 2``. For bonds
  pointing in independent uniformly-random directions the angle density is
  proportional to sin(theta); dividing by the bin-integrated reference makes
  that case flat at 1.0, which is the "sinth" normalization -- computed from
  the bin integral rather than ``1 / sin(theta_center)`` so the 0 and 180
  degree bins stay finite.

Conventions
-----------
Positions are fractions of the supercell (the ``.rmc6f`` storage convention)
and map to Cartesian angstrom through the row-vector product
``x_cart = frac @ lattice_vectors``, matching ``pca_kde.py``. Element symbols
are matched after ``str.capitalize()``, the same normalization
``parsers.iter_rmc6f_atoms`` applies.

Cross-engine caveat: ``workers/triplets.js`` reproduces this module's
histograms exactly for every tested configuration, but transcendental libm
functions (``acos``) may differ by 1 ulp between platforms, so a
*bitwise-ideal* geometry (an undisplaced average configuration whose cosine
lands exactly on a bin edge, e.g. cos = 0.5) can shift one count into the
neighbouring bin between engines. Real RMC configurations never hit this.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from pathlib import Path
from typing import Sequence

import numpy as np

from .parsers import iter_rmc6f_atoms, read_cell_vectors

# Cap on linked-list cells per lattice direction. Beyond this the per-cell
# occupancy for any realistic RMC box is far below one atom and finer cells
# only cost memory; capping keeps the cell-count arrays small for tiny rmax.
MAX_CELLS_PER_AXIS = 64

# Relative headroom on the neighbour-layer reach. The strict-inequality
# argument in the module docstring makes ceil() exact in real arithmetic; the
# headroom absorbs float rounding at exact-integer ratios at the cost of one
# spurious (empty) extra layer in those rare cases.
REACH_HEADROOM = 1e-9


@dataclass(frozen=True)
class BondAngleDistribution:
    """Bond-angle histogram for one triplet spec, plus its bond statistics.

    ``triplet`` is ``(A, B, C)`` with B central; ``bond12`` / ``bond23`` are
    the inclusive (rmin, rmax) windows in angstrom for the A--B and B--C
    bonds. Bin arrays are in degrees over [0, 180]. ``angles`` holds the raw
    angle list (degrees) only when the engine was asked to collect it.
    """

    triplet: tuple[str, str, str]
    bond12: tuple[float, float]
    bond23: tuple[float, float]
    bin_edges: np.ndarray  # (nbins + 1,) degrees
    bin_centers: np.ndarray  # (nbins,) degrees
    counts: np.ndarray  # (nbins,) int
    density: np.ndarray  # (nbins,) per-degree probability density
    sin_corrected: np.ndarray  # (nbins,) isotropic-reference enhancement
    angle_count: int
    mean_angle: float | None  # degrees; None when no angles
    std_angle: float | None
    apex_count: int  # atoms of the central element
    bond12_count: int
    bond23_count: int
    mean_length12: float | None  # angstrom; None when no bonds
    mean_length23: float | None
    angles: np.ndarray | None  # (angle_count,) degrees, optional


def _validate_window(name: str, window: Sequence[float]) -> tuple[float, float]:
    try:
        rmin, rmax = (float(value) for value in window)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a (rmin, rmax) pair") from error
    if not (np.isfinite(rmin) and np.isfinite(rmax)):
        raise ValueError(f"{name} bounds must be finite, got ({rmin}, {rmax})")
    if rmin < 0 or rmax <= rmin:
        raise ValueError(f"{name} needs 0 <= rmin < rmax, got ({rmin}, {rmax})")
    return rmin, rmax


def _bin_count(bin_width: float) -> int:
    """Bin count for a requested width: round half UP, not half to even.

    ``floor(x + 0.5)`` matches JavaScript's ``Math.round`` exactly, so the
    browser port builds the same number of bins for every width -- Python's
    banker's ``round()`` would disagree at exact .5 ratios (e.g. 8-degree
    bins: 180 / 8 = 22.5 -> 23 bins in both engines, not 22 vs 23).
    """
    return max(1, int(np.floor(180.0 / bin_width + 0.5)))


def _perpendicular_widths(lattice_vectors: np.ndarray) -> np.ndarray:
    """Perpendicular width of the box along each lattice direction.

    ``w_i = V / |a_j x a_k|`` is the distance between the two box faces
    spanned by the other two vectors -- the length scale that decides how many
    cells fit along direction ``i`` and how far a periodic image can reach.
    """
    a, b, c = lattice_vectors
    cross = np.stack([np.cross(b, c), np.cross(c, a), np.cross(a, b)])
    volume = abs(float(np.dot(a, np.cross(b, c))))
    areas = np.linalg.norm(cross, axis=1)
    if volume <= 0 or not np.all(areas > 0):
        raise ValueError("lattice_vectors are singular (zero cell volume)")
    return volume / areas


def _ragged_ranks(counts: np.ndarray) -> np.ndarray:
    """0..k-1 rank within consecutive groups of the given sizes, flattened."""
    total = int(counts.sum())
    starts = np.concatenate(([0], np.cumsum(counts)[:-1]))
    return np.arange(total) - np.repeat(starts, counts)


@dataclass(frozen=True)
class _Bonds:
    """All (center, candidate image) pairs whose length falls in one window.

    ``center_pos`` indexes into the *selection* array the search was given
    (not the full configuration) so angles can be grouped by central atom;
    ``candidate_row`` is the candidate's row in the full configuration, and
    with ``image`` -- the integer whole-box shift applied to it -- identifies
    a physical atom image uniquely across differently-selected bond lists.
    """

    center_pos: np.ndarray  # (P,) int
    candidate_row: np.ndarray  # (P,) int, global atom row
    image: np.ndarray  # (P, 3) int
    vectors: np.ndarray  # (P, 3) Cartesian angstrom, center -> candidate
    lengths: np.ndarray  # (P,)


def _neighbor_bonds(
    frac_centers: np.ndarray,
    center_rows: np.ndarray,
    frac_candidates: np.ndarray,
    candidate_rows: np.ndarray,
    lattice_vectors: np.ndarray,
    rmin: float,
    rmax: float,
) -> _Bonds:
    """Linked-cell periodic neighbour search with explicit image shifts."""
    widths = _perpendicular_widths(lattice_vectors)
    cells = np.minimum(
        np.maximum(1, np.floor(widths / rmax).astype(int)), MAX_CELLS_PER_AXIS
    )
    reach = np.ceil(rmax * cells / widths + REACH_HEADROOM).astype(int)

    def cell_of(frac: np.ndarray) -> np.ndarray:
        return np.minimum(np.floor(frac * cells).astype(int), cells - 1)

    candidate_cells = cell_of(frac_candidates)
    flat_candidates = (
        candidate_cells[:, 0] * cells[1] + candidate_cells[:, 1]
    ) * cells[2] + candidate_cells[:, 2]
    n_cells = int(np.prod(cells))
    order = np.argsort(flat_candidates, kind="stable")
    per_cell = np.bincount(flat_candidates, minlength=n_cells)
    cell_starts = np.concatenate(([0], np.cumsum(per_cell)[:-1]))

    center_cells = cell_of(frac_centers)
    n_centers = frac_centers.shape[0]
    rmin_sq, rmax_sq = rmin * rmin, rmax * rmax

    found: list[tuple[np.ndarray, ...]] = []
    offsets = product(*(range(-int(k), int(k) + 1) for k in reach))
    for offset in offsets:
        shifted = center_cells + np.asarray(offset)
        wrapped = np.mod(shifted, cells)
        # Whole-box image shift of the visited cell relative to its wrapped
        # copy; applied to the candidate to place its image beside the center.
        image = (shifted - wrapped) // cells
        flat = (wrapped[:, 0] * cells[1] + wrapped[:, 1]) * cells[2] + wrapped[:, 2]
        counts = per_cell[flat]
        if not counts.any():
            continue
        center_rep = np.repeat(np.arange(n_centers), counts)
        slots = cell_starts[flat][center_rep] + _ragged_ranks(counts)
        candidate_pos = order[slots]

        delta = (
            frac_candidates[candidate_pos]
            + image[center_rep]
            - frac_centers[center_rep]
        )
        # einsum, not `@`: numpy on Apple's Accelerate BLAS emits spurious
        # divide/overflow/invalid warnings for large (P, 3) @ (3, 3) matmuls
        # of finite values (results are correct; einsum skips that path).
        vectors = np.einsum("ij,jk->ik", delta, lattice_vectors)
        dist_sq = np.einsum("ij,ij->i", vectors, vectors)
        # The lower bound is exclusive at exactly zero even when rmin == 0: a
        # zero-length pair (bitwise-coincident atoms) has no direction, and
        # admitting it would put a 0/0 NaN in every angle it joins.
        keep = (dist_sq >= rmin_sq) & (dist_sq <= rmax_sq) & (dist_sq > 0)
        # A center is never its own neighbour in the unshifted image; other
        # images of the same atom are genuine neighbours and stay.
        keep &= ~(
            (center_rows[center_rep] == candidate_rows[candidate_pos])
            & np.all(image[center_rep] == 0, axis=1)
        )
        if not keep.any():
            continue
        found.append(
            (
                center_rep[keep],
                candidate_rows[candidate_pos[keep]],
                image[center_rep[keep]],
                vectors[keep],
                np.sqrt(dist_sq[keep]),
            )
        )

    if not found:
        empty = np.empty(0, dtype=int)
        return _Bonds(
            center_pos=empty,
            candidate_row=empty,
            image=np.empty((0, 3), dtype=int),
            vectors=np.empty((0, 3)),
            lengths=np.empty(0),
        )
    columns = [np.concatenate(parts) for parts in zip(*found)]
    return _Bonds(*columns)


def _sort_by_center(bonds: _Bonds, n_centers: int) -> tuple[_Bonds, np.ndarray, np.ndarray]:
    """Bonds reordered by central atom, with per-center group sizes/starts."""
    order = np.argsort(bonds.center_pos, kind="stable")
    sorted_bonds = _Bonds(
        center_pos=bonds.center_pos[order],
        candidate_row=bonds.candidate_row[order],
        image=bonds.image[order],
        vectors=bonds.vectors[order],
        lengths=bonds.lengths[order],
    )
    sizes = np.bincount(bonds.center_pos, minlength=n_centers)
    starts = np.concatenate(([0], np.cumsum(sizes)[:-1]))
    return sorted_bonds, sizes, starts


def _pair_angles(
    bonds12: _Bonds, bonds23: _Bonds, n_centers: int, shared_ends: bool
) -> np.ndarray:
    """Angles (degrees) at each center between its window-1 and window-2 bonds."""
    first, sizes1, starts1 = _sort_by_center(bonds12, n_centers)
    if shared_ends:
        second, sizes2, starts2 = first, sizes1, starts1
    else:
        second, sizes2, starts2 = _sort_by_center(bonds23, n_centers)

    pair_counts = sizes1 * sizes2
    total = int(pair_counts.sum())
    if total == 0:
        return np.empty(0)
    group = np.repeat(np.arange(n_centers), pair_counts)
    rank = _ragged_ranks(pair_counts)
    i = starts1[group] + rank // sizes2[group]
    j = starts2[group] + rank % sizes2[group]

    if shared_ends:
        # One list paired with itself: keep the strict upper triangle so each
        # unordered pair of distinct bonds counts once.
        keep = i < j
    else:
        # Two lists: drop the self-angle a bond makes with itself when both
        # windows catch the same atom in the same periodic image.
        keep = ~(
            (first.candidate_row[i] == second.candidate_row[j])
            & np.all(first.image[i] == second.image[j], axis=1)
        )
    i, j = i[keep], j[keep]
    if i.size == 0:
        return np.empty(0)

    cosine = np.einsum("ij,ij->i", first.vectors[i], second.vectors[j])
    cosine /= first.lengths[i] * second.lengths[j]
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


@dataclass(frozen=True)
class _TripletCore:
    """Shared mid-stage state between the public result builders."""

    triplet: tuple[str, str, str]
    window12: tuple[float, float]
    window23: tuple[float, float]
    shared_ends: bool
    apex_count: int
    bonds12: _Bonds
    bonds23: _Bonds
    angles: np.ndarray


def _triplet_core(
    fractional: np.ndarray,
    elements: Sequence[str],
    lattice_vectors: np.ndarray,
    triplet: Sequence[str],
    bond12: Sequence[float],
    bond23: Sequence[float] | None,
) -> _TripletCore:
    """Validate inputs, find both bond sets, and form every angle."""
    fractional = np.asarray(fractional, dtype=float)
    if fractional.ndim != 2 or fractional.shape[1] != 3:
        raise ValueError(f"fractional must be (N, 3), got {fractional.shape}")
    if not np.all(np.isfinite(fractional)):
        raise ValueError("fractional coordinates contain non-finite values")
    lattice_vectors = np.asarray(lattice_vectors, dtype=float)
    if lattice_vectors.shape != (3, 3) or not np.all(np.isfinite(lattice_vectors)):
        raise ValueError("lattice_vectors must be a finite (3, 3) matrix")
    symbols = [str(symbol).strip().capitalize() for symbol in elements]
    if len(symbols) != fractional.shape[0]:
        raise ValueError(
            f"{len(symbols)} elements for {fractional.shape[0]} coordinates"
        )
    if len(tuple(triplet)) != 3:
        raise ValueError("triplet must name three atom types (A, B, C)")
    end1, apex, end2 = (str(symbol).strip().capitalize() for symbol in triplet)
    window12 = _validate_window("bond12", bond12)
    window23 = window12 if bond23 is None else _validate_window("bond23", bond23)

    symbol_array = np.asarray(symbols)
    available = sorted(set(symbols))
    selections = {}
    for symbol in {end1, apex, end2}:
        rows = np.flatnonzero(symbol_array == symbol)
        if rows.size == 0:
            raise ValueError(
                f"No {symbol!r} atoms in the configuration; available: "
                + ", ".join(available)
            )
        selections[symbol] = rows

    # Fold every coordinate into [0, 1); the image bookkeeping restores the
    # true relative geometry, so pre-wrapped and drifted inputs agree.
    wrapped = fractional - np.floor(fractional)

    apex_rows = selections[apex]
    shared_ends = end1 == end2 and window12 == window23
    bonds12 = _neighbor_bonds(
        wrapped[apex_rows],
        apex_rows,
        wrapped[selections[end1]],
        selections[end1],
        lattice_vectors,
        *window12,
    )
    if shared_ends:
        bonds23 = bonds12
    else:
        bonds23 = _neighbor_bonds(
            wrapped[apex_rows],
            apex_rows,
            wrapped[selections[end2]],
            selections[end2],
            lattice_vectors,
            *window23,
        )
    angles = _pair_angles(bonds12, bonds23, apex_rows.size, shared_ends)

    return _TripletCore(
        triplet=(end1, apex, end2),
        window12=window12,
        window23=window23,
        shared_ends=shared_ends,
        apex_count=int(apex_rows.size),
        bonds12=bonds12,
        bonds23=bonds23,
        angles=angles,
    )


def bond_angle_distribution(
    fractional: np.ndarray,
    elements: Sequence[str],
    lattice_vectors: np.ndarray,
    *,
    triplet: Sequence[str],
    bond12: Sequence[float],
    bond23: Sequence[float] | None = None,
    bin_width: float = 1.0,
    collect_angles: bool = False,
) -> BondAngleDistribution:
    """Histogram the A--B--C bond angles of a periodic configuration.

    ``fractional`` are supercell-fraction coordinates (N, 3); ``elements`` the
    matching symbols; ``lattice_vectors`` the (3, 3) supercell rows in
    angstrom. ``triplet`` is ``(A, B, C)`` with **B the central atom**;
    ``bond12`` bounds the A--B length and ``bond23`` the B--C length
    (inclusive; ``None`` reuses ``bond12``). ``bin_width`` is in degrees --
    the bin count is ``round(180 / bin_width)``, so widths that do not divide
    180 are adjusted to the nearest exact tiling.
    """
    bin_width = float(bin_width)
    if not (np.isfinite(bin_width) and 0 < bin_width <= 180):
        raise ValueError(f"bin_width must be in (0, 180], got {bin_width}")
    core = _triplet_core(fractional, elements, lattice_vectors, triplet, bond12, bond23)
    angles, bonds12, bonds23 = core.angles, core.bonds12, core.bonds23

    nbins = _bin_count(bin_width)
    width = 180.0 / nbins
    edges = np.linspace(0.0, 180.0, nbins + 1)
    counts, _ = np.histogram(angles, bins=nbins, range=(0.0, 180.0))
    total = int(counts.sum())

    edges_rad = np.radians(edges)
    # Exact isotropic reference fraction per bin: integral of sin(theta)/2.
    isotropic = (np.cos(edges_rad[:-1]) - np.cos(edges_rad[1:])) / 2.0
    if total:
        density = counts / (total * width)
        sin_corrected = (counts / total) / isotropic
    else:
        density = np.zeros(nbins)
        sin_corrected = np.zeros(nbins)

    def stats(values: np.ndarray) -> tuple[float | None, float | None]:
        if values.size == 0:
            return None, None
        return float(np.mean(values)), float(np.std(values))

    mean_angle, std_angle = stats(angles)
    mean12 = float(np.mean(bonds12.lengths)) if bonds12.lengths.size else None
    mean23 = float(np.mean(bonds23.lengths)) if bonds23.lengths.size else None

    return BondAngleDistribution(
        triplet=core.triplet,
        bond12=core.window12,
        bond23=core.window23,
        bin_edges=edges,
        bin_centers=(edges[:-1] + edges[1:]) / 2.0,
        counts=counts,
        density=density,
        sin_corrected=sin_corrected,
        angle_count=int(angles.size),
        mean_angle=mean_angle,
        std_angle=std_angle,
        apex_count=core.apex_count,
        bond12_count=int(bonds12.lengths.size),
        bond23_count=int(bonds23.lengths.size),
        mean_length12=mean12,
        mean_length23=mean23,
        angles=angles if collect_angles else None,
    )


# Bond-length histogram resolution inside each window, used by the summary
# payload. Fixed count (not width) so any window renders at the same detail.
LENGTH_BINS = 40


def bond_angle_summary(
    fractional: np.ndarray,
    elements: Sequence[str],
    lattice_vectors: np.ndarray,
    *,
    triplet: Sequence[str],
    bond12: Sequence[float],
    bond23: Sequence[float] | None = None,
    bin_width: float = 1.0,
) -> dict:
    """JSON-safe payload for the app: angles + bond lengths + coordination.

    One engine pass feeding every panel of the Local Geometry page. The dict
    (camelCase keys, plain lists/scalars) is the payload contract shared by
    the Flask ``/api/triplets`` route and the browser worker's ``triplets``
    request -- keep ``workers/triplets.js`` in sync with any change here.
    """
    bin_width = float(bin_width)
    if not (np.isfinite(bin_width) and 0 < bin_width <= 180):
        raise ValueError(f"bin_width must be in (0, 180], got {bin_width}")
    core = _triplet_core(fractional, elements, lattice_vectors, triplet, bond12, bond23)

    nbins = _bin_count(bin_width)
    width = 180.0 / nbins
    edges = np.linspace(0.0, 180.0, nbins + 1)
    counts, _ = np.histogram(core.angles, bins=nbins, range=(0.0, 180.0))
    total = int(counts.sum())
    edges_rad = np.radians(edges)
    isotropic = (np.cos(edges_rad[:-1]) - np.cos(edges_rad[1:])) / 2.0
    density = counts / (total * width) if total else np.zeros(nbins)
    sin_corrected = (counts / total) / isotropic if total else np.zeros(nbins)

    def length_histogram(bonds: _Bonds, window: tuple[float, float]) -> dict:
        length_counts, length_edges = np.histogram(
            bonds.lengths, bins=LENGTH_BINS, range=window
        )
        return {
            "binCenters": ((length_edges[:-1] + length_edges[1:]) / 2.0).tolist(),
            "counts": length_counts.tolist(),
            "count": int(bonds.lengths.size),
            "meanLength": float(np.mean(bonds.lengths)) if bonds.lengths.size else None,
        }

    coordination = np.bincount(
        np.bincount(core.bonds12.center_pos, minlength=core.apex_count)
    )

    angles = core.angles
    return {
        "triplet": list(core.triplet),
        "bond12": list(core.window12),
        "bond23": list(core.window23),
        "sharedEnds": core.shared_ends,
        "binWidth": width,
        "binCenters": ((edges[:-1] + edges[1:]) / 2.0).tolist(),
        "counts": counts.tolist(),
        "density": density.tolist(),
        "sinCorrected": sin_corrected.tolist(),
        "angleCount": int(angles.size),
        "meanAngle": float(np.mean(angles)) if angles.size else None,
        "stdAngle": float(np.std(angles)) if angles.size else None,
        "apexCount": core.apex_count,
        "lengths12": length_histogram(core.bonds12, core.window12),
        # Shared ends reuse the window-1 bonds, so the page shows one histogram.
        "lengths23": None if core.shared_ends else length_histogram(core.bonds23, core.window23),
        # coordination[n] = how many central atoms have exactly n window-1 bonds.
        "coordination": coordination.tolist(),
    }


@lru_cache(maxsize=16)
def cached_bond_angle_summary(
    path: str,
    mtime: float,
    end1: str,
    apex: str,
    end2: str,
    r12_min: float,
    r12_max: float,
    r23_min: float,
    r23_max: float,
    bin_width: float,
) -> dict:
    """``bond_angle_summary`` of an ``.rmc6f`` file, memoized for API callers.

    Keyed on (path, mtime) plus every parameter, mirroring
    ``pca_kde.cached_site_displacements``. Windows arrive resolved (bond23
    defaults applied by the caller) so equal windows hit one cache entry.
    """
    lattice_vectors, _ = read_cell_vectors(path)
    coords: list[np.ndarray] = []
    elements: list[str] = []
    for atom in iter_rmc6f_atoms(path):
        coords.append(atom["coords"])
        elements.append(atom["element"])
    if not coords:
        raise ValueError(f"{path} does not contain any atoms")
    return bond_angle_summary(
        np.asarray(coords, dtype=float),
        elements,
        lattice_vectors,
        triplet=(end1, apex, end2),
        bond12=(r12_min, r12_max),
        bond23=(r23_min, r23_max),
        bin_width=bin_width,
    )


def bond_angles_from_rmc6f(
    rmc6f_path: str | Path,
    *,
    triplet: Sequence[str],
    bond12: Sequence[float],
    bond23: Sequence[float] | None = None,
    bin_width: float = 1.0,
    collect_angles: bool = False,
) -> BondAngleDistribution:
    """Run ``bond_angle_distribution`` on an ``.rmc6f`` configuration file."""
    rmc6f_path = Path(rmc6f_path)
    lattice_vectors, _ = read_cell_vectors(rmc6f_path)
    coords: list[np.ndarray] = []
    elements: list[str] = []
    for atom in iter_rmc6f_atoms(rmc6f_path):
        coords.append(atom["coords"])
        elements.append(atom["element"])
    if not coords:
        raise ValueError(f"{rmc6f_path} does not contain any atoms")
    return bond_angle_distribution(
        np.asarray(coords, dtype=float),
        elements,
        lattice_vectors,
        triplet=triplet,
        bond12=bond12,
        bond23=bond23,
        bin_width=bin_width,
        collect_angles=collect_angles,
    )
