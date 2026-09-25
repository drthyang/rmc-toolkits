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
"Inclusive" holds for ideal geometries too: a distance within
``WINDOW_TOL`` (1e-9 A) of a bound counts as on it, so a bound typed exactly
at an ideal shell distance keeps the whole shell instead of whichever bonds
float rounding happens to leave inside.

A bond vector is ``((f_cand - f_center) + m) @ L``: the fractional difference
is taken *before* the integer image shift ``m`` is added, so the vector from
the other end, ``((f_center - f_cand) - m) @ L``, is its exact negative
(IEEE rounding is symmetric under negation). A bond therefore has the same
length from both of its ends and sits inside or outside a window for both
-- which is what makes the undirected bond count below exact.

Bond counts
-----------
Bonds are found from each central atom, so ``bond12_count`` counts
B-centred bond vectors. When the end element is the central element (A = B,
e.g. Nb-Nb-Nb) every bond is found from both of its ends, and the directed
count is exactly twice the number of physical bonds; ``unique_bonds12`` /
``unique_bonds23`` report the physical (undirected) count -- half the
directed one when the end element is the central one, the directed count
otherwise. Per-central-atom coordination is the directed count per B.

Angle counting
--------------
For each central B atom the A-bond list and the C-bond list are combined:

- If A and C name the same element, every angle is one physical triplet
  ``{x, B, y}``: each *unordered* pair of distinct bonds (distinct atom
  images) contributes one angle when one bond lies in the A--B window and the
  other in the B--C window, under either assignment. With equal windows that
  is every unordered pair, so an octahedron's six bonds give the expected
  ``C(6, 2) = 15`` angles (12 x 90 deg + 3 x 180 deg). The rule is continuous
  in the windows: moving a bound changes the count only by the triplets whose
  bonds cross it (an earlier ordered rule counted every overlap triplet twice
  as soon as the two windows differed by any amount).
- Otherwise every (A-bond, C-bond) combination contributes one angle; A and C
  atoms are then always distinct atoms.

Work and memory
---------------
The exact angle count follows from the per-center bond counts before any
angle exists (``_angles_per_center``). Angles are then streamed: centers are
paired in runs of at most ``PAIR_CHUNK`` bond pairs, each run is binned and
folded into running moments, and nothing outlives its run unless
``collect_angles`` asks for the raw list -- so memory does not grow with the
angle count, which scales as ~rmax^6. With ``max_angles`` set (the app
boundaries pass ``APP_MAX_ANGLES``) a count-only search that stores nothing
refuses an oversized spec before either the bond lists or any angle is built.

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
  that case flat at exactly 1.0. Since ``cos(c - w/2) - cos(c + w/2) =
  2 sin(c) sin(w/2)``, the reference *is* ``sin(c) sin(w/2)``: the curve is
  the bin-centre ``1 / sin(c)`` correction times the constant
  ``1 / sin(w/2)`` -- the same shape, scaled so random reads 1. Neither form
  diverges (bin centres lie in [w/2, 180 - w/2]); only a per-angle
  ``1 / sin(theta_i)`` weight would.
- Relation to RMCProfile's TRIPLETS output: its ``norm`` column is
  ``density`` and its ``norm/sin(theta)`` column is ``density / sin(c)`` =
  ``sin_corrected * sin(w/2) / w_deg`` (``w/2`` in radians), ~``pi / 360``
  -- a constant factor, so the same shape but not the same numbers.

Conventions
-----------
Positions are fractions of the supercell (the ``.rmc6f`` storage convention)
and map to Cartesian angstrom through the row-vector product
``x_cart = frac @ lattice_vectors``, matching ``pca_kde.py``. Element symbols
are matched after ``str.capitalize()``, the same normalization
``parsers.iter_rmc6f_atoms`` applies.

Bin edges and ideal geometries: bins are half-open ``[edge_k, edge_k+1)``
(the last one closed at 180 deg), numpy's convention. Edges are multiples of
``180 / nbins``, so every symmetry angle of an undisplaced configuration
(60, 90, 120 deg ... -- an RMCProfile start configuration built from a CIF)
sits exactly on one, and float noise puts its computed value a few ulp
either side. An angle within ``EDGE_SNAP_DEG`` (1e-9 deg) of an edge is
therefore binned as if exactly on it -- into the bin that edge starts -- so
a symmetry class lands whole in one bin, whatever the rounding, a rigid
shift of the configuration, or the platform's ``acos`` (numpy's and V8's
differ by 1 ulp on ~17% of inputs). ``workers/triplets.js`` applies the same
rule, and the two engines then agree bin for bin on ideal geometries too.
Only the binning snaps; means, standard deviations and raw angles keep the
computed values.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from itertools import product
from pathlib import Path
from typing import Sequence

import numpy as np

from .parsers import Rmc6fParseReport, iter_rmc6f_atoms, read_cell_vectors

# Cap on linked-list cells per lattice direction. Beyond this the per-cell
# occupancy for any realistic RMC box is far below one atom and finer cells
# only cost memory; capping keeps the cell-count arrays small for tiny rmax.
MAX_CELLS_PER_AXIS = 64

# Relative headroom on the neighbour-layer reach. The strict-inequality
# argument in the module docstring makes ceil() exact in real arithmetic; the
# headroom absorbs float rounding at exact-integer ratios at the cost of one
# spurious (empty) extra layer in those rare cases.
REACH_HEADROOM = 1e-9

# Work budget for one app request: the exact number of angles a spec forms is
# counted from the bond lists *before* any angle exists, and the Flask route
# and the browser worker both refuse a spec above this (400 / thrown Error).
# The rmax <= 15 A request cap bounds the neighbour search, not this: the
# angle count grows ~rmax^6 (1.3e9 Se-Nb-Se angles at 15 A on the 52 000-atom
# sample). Library and CLI callers are unrestricted (``max_angles=None``).
# Mirrored as APP_MAX_ANGLES in workers/triplets.js -- keep the two equal.
APP_MAX_ANGLES = 50_000_000

# Distances within this many angstrom of a window bound count as on it
# (inside: windows are inclusive). Far above float noise in a length
# (~1e-13 A) and far below any real displacement. Mirrored as WINDOW_TOL in
# workers/triplets.js.
WINDOW_TOL = 1e-9

# Angles closer than this to a bin edge are binned as exactly on it (see the
# module docstring). Far above float noise in an angle (~1e-13 deg) and far
# below any meaningful bin width or real displacement. Mirrored as
# EDGE_SNAP_DEG in workers/triplets.js.
EDGE_SNAP_DEG = 1e-9

# Streaming chunk: the bond pairs formed at once while histogramming. Peak
# pairing memory is ~100 B per pair, so ~25 MB here, whatever the angle count.
PAIR_CHUNK = 1 << 18

# Neighbour-search block: candidate pairs examined at once per stencil offset
# (~100 B each), so the search's transient memory is ~100 MB at most.
SEARCH_CHUNK = 1 << 20


@dataclass(frozen=True)
class BondAngleDistribution:
    """Bond-angle histogram for one triplet spec, plus its bond statistics.

    ``triplet`` is ``(A, B, C)`` with B central; ``bond12`` / ``bond23`` are
    the inclusive (rmin, rmax) windows in angstrom for the A--B and B--C
    bonds. Bin arrays are in degrees over [0, 180]. ``angles`` holds the raw
    angle list (degrees) only when the engine was asked to collect it.
    ``bond12_count`` / ``bond23_count`` count B-centred bond vectors (an
    A--B bond with A = B is found from both ends and counted twice);
    ``unique_bonds12`` / ``unique_bonds23`` count physical bonds once.
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
    # Physical (undirected) bonds: bond12_count / 2 when A = B (each bond is
    # found from both ends), else bond12_count. Same for window 2 with C.
    unique_bonds12: int
    unique_bonds23: int
    # What the .rmc6f atom section skipped (Rmc6fParseReport.warning()), set
    # only by the file loaders; None for a clean file or in-memory input.
    parse_warning: str | None = None


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


def _angle_bins(angles: np.ndarray, nbins: int) -> np.ndarray:
    """Bin index of each angle (degrees) over [0, 180] in ``nbins`` bins.

    Exactly numpy.histogram's uniform-bin assignment -- a truncated first
    guess corrected against the ``linspace`` edge values, half-open bins, the
    last one closed at 180 -- written out as ``histogramIndex`` in the JS
    port, plus the edge snap: an angle within ``EDGE_SNAP_DEG`` of an edge
    goes to the bin that edge starts (180 to the last bin).
    """
    width = 180.0 / nbins
    edges = np.linspace(0.0, 180.0, nbins + 1)
    index = np.floor(angles * (nbins / 180.0)).astype(np.int64)
    np.clip(index, 0, nbins - 1, out=index)
    index -= angles < edges[index]
    index += (angles >= edges[index + 1]) & (index != nbins - 1)
    nearest = np.floor(angles / width + 0.5)
    snap = np.abs(angles - nearest * width) < EDGE_SNAP_DEG
    index[snap] = np.minimum(nearest[snap].astype(np.int64), nbins - 1)
    return index


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
    member: np.ndarray  # (P, W) bool: which of the searched windows hold it

    def subset(self, mask: np.ndarray) -> "_Bonds":
        return _Bonds(
            center_pos=self.center_pos[mask],
            candidate_row=self.candidate_row[mask],
            image=self.image[mask],
            vectors=self.vectors[mask],
            lengths=self.lengths[mask],
            member=self.member[mask],
        )


def _neighbor_bonds(
    frac_centers: np.ndarray,
    center_rows: np.ndarray,
    frac_candidates: np.ndarray,
    candidate_rows: np.ndarray,
    lattice_vectors: np.ndarray,
    windows: Sequence[tuple[float, float]],
    *,
    count_only: bool = False,
) -> _Bonds | np.ndarray:
    """Linked-cell periodic neighbour search with explicit image shifts.

    One pass serves every window in ``windows``: a pair is kept when it falls
    in at least one of them, and ``member[:, w]`` records window ``w``.

    Centers are processed in blocks sized so one stencil offset examines at
    most ~``SEARCH_CHUNK`` candidate pairs, bounding the transient memory
    whatever rmax and the box. With ``count_only`` nothing is stored: the
    result is an ``(n_centers, 2**W)`` table counting each center's bonds by
    window-membership pattern (bit ``w`` = window ``w``) -- exactly the pairs
    the storing search would keep, since both run the same arithmetic.
    """
    rmax = max(window[1] for window in windows)
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
    # Inclusive bounds, widened by WINDOW_TOL against float noise.
    bounds_sq = [
        (max(lo - WINDOW_TOL, 0.0) ** 2, (hi + WINDOW_TOL) ** 2) for lo, hi in windows
    ]
    n_patterns = 1 << len(windows)
    pattern_weights = 1 << np.arange(len(windows))
    patterns = np.zeros((n_centers, n_patterns), dtype=np.int64) if count_only else None
    block = max(1, SEARCH_CHUNK // max(1, int(per_cell.max(initial=0))))
    offsets = [np.asarray(offset) for offset in product(*(range(-int(k), int(k) + 1) for k in reach))]

    found: list[tuple[np.ndarray, ...]] = []
    for block_start in range(0, n_centers, block):
        block_stop = min(block_start + block, n_centers)
        block_cells = center_cells[block_start:block_stop]
        block_frac = frac_centers[block_start:block_stop]
        block_rows = center_rows[block_start:block_stop]
        for offset in offsets:
            shifted = block_cells + offset
            wrapped = np.mod(shifted, cells)
            # Whole-box image shift of the visited cell relative to its wrapped
            # copy; applied to the candidate to place its image beside the center.
            image = (shifted - wrapped) // cells
            flat = (wrapped[:, 0] * cells[1] + wrapped[:, 1]) * cells[2] + wrapped[:, 2]
            counts = per_cell[flat]
            if not counts.any():
                continue
            local = np.repeat(np.arange(block_stop - block_start), counts)
            slots = cell_starts[flat][local] + _ragged_ranks(counts)
            candidate_pos = order[slots]

            # Difference first, image shift second: the vector from the other
            # end of the same bond is then the exact negative of this one.
            delta = (frac_candidates[candidate_pos] - block_frac[local]) + image[local]
            # Row-vector product written out term by term -- the evaluation
            # order of workers/triplets.js, so both engines produce bitwise
            # identical vectors (and no BLAS matmul, whose Accelerate build
            # emits spurious warnings on large finite inputs).
            vectors = (
                delta[:, 0:1] * lattice_vectors[0] + delta[:, 1:2] * lattice_vectors[1]
            ) + delta[:, 2:3] * lattice_vectors[2]
            dist_sq = (
                vectors[:, 0] * vectors[:, 0] + vectors[:, 1] * vectors[:, 1]
            ) + vectors[:, 2] * vectors[:, 2]
            member = np.stack(
                [(dist_sq >= lo_sq) & (dist_sq <= hi_sq) for lo_sq, hi_sq in bounds_sq],
                axis=1,
            )
            # The lower bound is exclusive at exactly zero even when rmin == 0: a
            # zero-length pair (bitwise-coincident atoms) has no direction, and
            # admitting it would put a 0/0 NaN in every angle it joins.
            keep = member.any(axis=1) & (dist_sq > 0)
            # A center is never its own neighbour in the unshifted image; other
            # images of the same atom are genuine neighbours and stay.
            keep &= ~(
                (block_rows[local] == candidate_rows[candidate_pos])
                & np.all(image[local] == 0, axis=1)
            )
            if not keep.any():
                continue
            if count_only:
                pattern = member[keep] @ pattern_weights
                patterns[block_start:block_stop] += np.bincount(
                    local[keep] * n_patterns + pattern,
                    minlength=(block_stop - block_start) * n_patterns,
                ).reshape(-1, n_patterns)
                continue
            found.append(
                (
                    local[keep] + block_start,
                    candidate_rows[candidate_pos[keep]],
                    image[local[keep]],
                    vectors[keep],
                    np.sqrt(dist_sq[keep]),
                    member[keep],
                )
            )

    if count_only:
        return patterns
    if not found:
        empty = np.empty(0, dtype=int)
        return _Bonds(
            center_pos=empty,
            candidate_row=empty,
            image=np.empty((0, 3), dtype=int),
            vectors=np.empty((0, 3)),
            lengths=np.empty(0),
            member=np.empty((0, len(windows)), dtype=bool),
        )
    columns = [np.concatenate(parts) for parts in zip(*found)]
    return _Bonds(*columns)


def _membership_patterns(bonds: _Bonds, n_centers: int) -> np.ndarray:
    """Per-center bond counts by window-membership pattern, from stored bonds."""
    n_windows = bonds.member.shape[1]
    pattern = bonds.member @ (1 << np.arange(n_windows))
    return np.bincount(
        bonds.center_pos * (1 << n_windows) + pattern,
        minlength=n_centers * (1 << n_windows),
    ).reshape(n_centers, 1 << n_windows)


def _angles_per_center(
    same_end: bool, patterns1: np.ndarray, patterns2: np.ndarray | None = None
) -> np.ndarray:
    """Exact number of angles at each center, from membership-pattern counts.

    Different ends (``patterns1`` for A--B, ``patterns2`` for B--C, one window
    each): every (A-bond, C-bond) combination, ``n1 * n2``. Same end (one
    table over both windows; first window A--B, last B--C): the combined
    list's unordered pairs minus those with no bond in one of the windows,
    ``C(a+b+c, 2) - C(a, 2) - C(b, 2)`` with ``a`` bonds only in A--B, ``b``
    only in B--C and ``c`` in both.
    """
    if not same_end:
        return patterns1[:, 1].astype(np.int64) * patterns2[:, 1].astype(np.int64)
    n_patterns = patterns1.shape[1]
    last_bit = n_patterns.bit_length() - 2  # window count - 1
    pattern = np.arange(n_patterns)
    in12 = (pattern & 1) != 0
    in23 = ((pattern >> last_bit) & 1) != 0
    total = patterns1.sum(axis=1)
    only12 = patterns1[:, in12 & ~in23].sum(axis=1)
    only23 = patterns1[:, in23 & ~in12].sum(axis=1)
    return _pairs2(total) - _pairs2(only12) - _pairs2(only23)


def _sort_by_center(bonds: _Bonds, n_centers: int) -> tuple[_Bonds, np.ndarray, np.ndarray]:
    """Bonds reordered by central atom, with per-center group sizes/starts."""
    order = np.argsort(bonds.center_pos, kind="stable")
    sorted_bonds = bonds.subset(order)
    sizes = np.bincount(bonds.center_pos, minlength=n_centers)
    starts = np.concatenate(([0], np.cumsum(sizes)[:-1]))
    return sorted_bonds, sizes, starts


def _pairs2(n: np.ndarray) -> np.ndarray:
    """``n (n - 1) / 2``: unordered pairs among ``n`` items, as int64."""
    n = n.astype(np.int64)
    return n * (n - 1) // 2


@dataclass(frozen=True)
class _Pairing:
    """Per-center bond lists ready to be paired into angles, chunk by chunk.

    ``same_end``: A and C are one element and ``first`` is the single search
    over both windows (``second`` is the same object); its ``member`` columns
    say which window each bond is in -- first column A--B, last column B--C.
    Otherwise ``first`` holds the A--B bonds and ``second`` the B--C bonds.
    """

    first: _Bonds
    sizes1: np.ndarray
    starts1: np.ndarray
    second: _Bonds
    sizes2: np.ndarray
    starts2: np.ndarray
    same_end: bool

    def angle_counts(self) -> np.ndarray:
        """Exact number of angles at each center (see ``_angles_per_center``)."""
        n_centers = self.sizes1.size
        if self.same_end:
            return _angles_per_center(True, _membership_patterns(self.first, n_centers))
        return _angles_per_center(
            False,
            _membership_patterns(self.first, n_centers),
            _membership_patterns(self.second, n_centers),
        )

    def chunk_angles(self, start: int, stop: int) -> np.ndarray:
        """Angles (degrees) at centers ``start:stop``."""
        sizes1 = self.sizes1[start:stop]
        sizes2 = self.sizes2[start:stop]
        pair_counts = sizes1 * sizes2
        if not pair_counts.any():
            return np.empty(0)
        local = np.repeat(np.arange(stop - start), pair_counts)
        rank = _ragged_ranks(pair_counts)
        i = self.starts1[start:stop][local] + rank // sizes2[local]
        j = self.starts2[start:stop][local] + rank % sizes2[local]
        del local, rank
        first, second = self.first, self.second
        if self.same_end:
            # One list paired with itself: the strict upper triangle makes each
            # unordered pair of distinct bond images count once, and the pair
            # is a triplet when either assignment puts one bond in each window.
            in12, in23 = first.member[:, 0], first.member[:, -1]
            keep = (i < j) & ((in12[i] & in23[j]) | (in23[i] & in12[j]))
            i, j = i[keep], j[keep]
        if i.size == 0:
            return np.empty(0)
        one, two = first.vectors[i], second.vectors[j]
        # Written out in the JS port's evaluation order (bitwise parity).
        cosine = (one[:, 0] * two[:, 0] + one[:, 1] * two[:, 1]) + one[:, 2] * two[:, 2]
        cosine /= first.lengths[i] * second.lengths[j]
        return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


@dataclass(frozen=True)
class _AngleHistogram:
    counts: np.ndarray  # (nbins,) int64
    angle_count: int
    mean: float | None
    std: float | None
    angles: np.ndarray | None  # only when collected


def _stream_angles(pairing: _Pairing, nbins: int, collect: bool) -> _AngleHistogram:
    """Histogram every angle chunk by chunk; memory is O(PAIR_CHUNK).

    Centers are taken in consecutive runs whose combined bond-pair count
    stays within ``PAIR_CHUNK`` (a single center above it forms a chunk of
    its own). Counts accumulate per chunk; mean and variance combine per
    chunk with Chan's parallel update, so no angle outlives its chunk unless
    ``collect`` asks for the raw list.
    """
    pair_counts = pairing.sizes1.astype(np.int64) * pairing.sizes2.astype(np.int64)
    cumulative = np.cumsum(pair_counts)
    n_centers = pair_counts.size
    counts = np.zeros(nbins, dtype=np.int64)
    total, mean, m2 = 0, 0.0, 0.0
    parts: list[np.ndarray] = []
    start = 0
    while start < n_centers:
        base = int(cumulative[start - 1]) if start else 0
        stop = int(np.searchsorted(cumulative, base + PAIR_CHUNK, side="right"))
        stop = min(max(stop, start + 1), n_centers)
        angles = pairing.chunk_angles(start, stop)
        start = stop
        if angles.size == 0:
            continue
        counts += np.bincount(_angle_bins(angles, nbins), minlength=nbins)
        size = int(angles.size)
        chunk_mean = float(np.mean(angles))
        chunk_m2 = float(np.sum((angles - chunk_mean) ** 2))
        combined = total + size
        delta = chunk_mean - mean
        mean += delta * size / combined
        m2 += chunk_m2 + delta * delta * total * size / combined
        total = combined
        if collect:
            parts.append(angles)
    return _AngleHistogram(
        counts=counts,
        angle_count=total,
        mean=mean if total else None,
        std=float(np.sqrt(m2 / total)) if total else None,
        angles=(np.concatenate(parts) if parts else np.empty(0)) if collect else None,
    )


def _unique_bonds(bonds: _Bonds, end: str, apex: str) -> int:
    """Physical (undirected) bond count of a B-centred bond list.

    With the end element equal to the central element every bond is found
    from both of its ends -- exactly, since a bond's two vectors are exact
    negatives (see ``_neighbor_bonds``) -- so the directed count is even and
    halves; otherwise each bond is found once, from its B end.
    """
    directed = int(bonds.lengths.size)
    if end != apex:
        return directed
    if directed % 2:
        raise RuntimeError(
            f"internal error: {directed} directed {end}-{apex} bonds do not pair up"
        )
    return directed // 2


@dataclass(frozen=True)
class _TripletCore:
    """Shared mid-stage state between the public result builders."""

    triplet: tuple[str, str, str]
    window12: tuple[float, float]
    window23: tuple[float, float]
    shared_ends: bool
    apex_count: int
    bonds12: _Bonds  # the A--B window's bonds
    bonds23: _Bonds  # the B--C window's bonds
    unique_bonds12: int  # physical (undirected) A--B bonds
    unique_bonds23: int
    pairing: _Pairing
    angle_count: int  # exact, counted before any angle is formed


def _triplet_core(
    fractional: np.ndarray,
    elements: Sequence[str],
    lattice_vectors: np.ndarray,
    triplet: Sequence[str],
    bond12: Sequence[float],
    bond23: Sequence[float] | None,
    max_angles: int | None = None,
) -> _TripletCore:
    """Validate inputs, find both bond sets, and count the angles exactly.

    With ``max_angles`` set, a spec that would form more angles than that is
    refused with ``ValueError`` *before* any angle is formed.
    """
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
    n_centers = int(apex_rows.size)
    same_end = end1 == end2
    shared_ends = same_end and window12 == window23

    def search(end: str, windows: list[tuple[float, float]], count_only: bool = False):
        return _neighbor_bonds(
            wrapped[apex_rows],
            apex_rows,
            wrapped[selections[end]],
            selections[end],
            lattice_vectors,
            windows,
            count_only=count_only,
        )

    windows_same_end = [window12] if shared_ends else [window12, window23]

    def refuse(angle_count: int) -> None:
        raise ValueError(
            f"{end1}-{apex}-{end2} with bond windows {window12[0]:g}-{window12[1]:g} / "
            f"{window23[0]:g}-{window23[1]:g} A would form {angle_count:,} angles, over "
            f"the limit of {int(max_angles):,} for one request; narrow the bond windows "
            "(the angle count grows roughly as rmax^6)"
        )

    if max_angles is not None:
        # Budgeted (app) request: an exact count from a search that stores
        # nothing, so an oversized spec is refused before either the bond
        # lists or a single angle take up memory.
        if same_end:
            per_center = _angles_per_center(True, search(end1, windows_same_end, True))
        else:
            per_center = _angles_per_center(
                False, search(end1, [window12], True), search(end2, [window23], True)
            )
        counted = int(per_center.sum())
        if counted > max_angles:
            refuse(counted)

    if same_end:
        # One search over both windows: each bond image appears once, tagged
        # with the window(s) it falls in, so pairing can count each physical
        # triplet once however the two windows overlap.
        both = search(end1, windows_same_end)
        first, sizes1, starts1 = _sort_by_center(both, n_centers)
        pairing = _Pairing(first, sizes1, starts1, first, sizes1, starts1, True)
        bonds12 = both.subset(both.member[:, 0])
        bonds23 = bonds12 if shared_ends else both.subset(both.member[:, -1])
    else:
        bonds12 = search(end1, [window12])
        bonds23 = search(end2, [window23])
        pairing = _Pairing(
            *_sort_by_center(bonds12, n_centers),
            *_sort_by_center(bonds23, n_centers),
            False,
        )

    angle_count = int(pairing.angle_counts().sum())

    return _TripletCore(
        triplet=(end1, apex, end2),
        window12=window12,
        window23=window23,
        shared_ends=shared_ends,
        apex_count=n_centers,
        bonds12=bonds12,
        bonds23=bonds23,
        unique_bonds12=_unique_bonds(bonds12, end1, apex),
        unique_bonds23=_unique_bonds(bonds23, end2, apex),
        pairing=pairing,
        angle_count=angle_count,
    )


def _validate_bin_width(bin_width: float) -> float:
    bin_width = float(bin_width)
    if not (np.isfinite(bin_width) and 0 < bin_width <= 180):
        raise ValueError(f"bin_width must be in (0, 180], got {bin_width}")
    return bin_width


def _normalized(counts: np.ndarray, nbins: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bin edges, per-degree density and sin-corrected curve for ``counts``."""
    width = 180.0 / nbins
    edges = np.linspace(0.0, 180.0, nbins + 1)
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
    return edges, density, sin_corrected


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
    max_angles: int | None = None,
) -> BondAngleDistribution:
    """Histogram the A--B--C bond angles of a periodic configuration.

    ``fractional`` are supercell-fraction coordinates (N, 3); ``elements`` the
    matching symbols; ``lattice_vectors`` the (3, 3) supercell rows in
    angstrom. ``triplet`` is ``(A, B, C)`` with **B the central atom**;
    ``bond12`` bounds the A--B length and ``bond23`` the B--C length
    (inclusive; ``None`` reuses ``bond12``). ``bin_width`` is in degrees --
    the bin count is ``round(180 / bin_width)``, so widths that do not divide
    180 are adjusted to the nearest exact tiling.

    Angles are streamed into the histogram, so memory does not grow with the
    angle count unless ``collect_angles`` keeps the raw list. ``max_angles``
    (default unlimited) refuses, before pairing, a spec whose exact angle
    count exceeds it -- the app boundaries pass ``APP_MAX_ANGLES``.
    """
    bin_width = _validate_bin_width(bin_width)
    core = _triplet_core(
        fractional, elements, lattice_vectors, triplet, bond12, bond23, max_angles
    )
    nbins = _bin_count(bin_width)
    histogram = _stream_angles(core.pairing, nbins, collect_angles)
    edges, density, sin_corrected = _normalized(histogram.counts, nbins)
    bonds12, bonds23 = core.bonds12, core.bonds23
    mean12 = float(np.mean(bonds12.lengths)) if bonds12.lengths.size else None
    mean23 = float(np.mean(bonds23.lengths)) if bonds23.lengths.size else None

    return BondAngleDistribution(
        triplet=core.triplet,
        bond12=core.window12,
        bond23=core.window23,
        bin_edges=edges,
        bin_centers=(edges[:-1] + edges[1:]) / 2.0,
        counts=histogram.counts,
        density=density,
        sin_corrected=sin_corrected,
        angle_count=histogram.angle_count,
        mean_angle=histogram.mean,
        std_angle=histogram.std,
        apex_count=core.apex_count,
        bond12_count=int(bonds12.lengths.size),
        bond23_count=int(bonds23.lengths.size),
        mean_length12=mean12,
        mean_length23=mean23,
        angles=histogram.angles,
        unique_bonds12=core.unique_bonds12,
        unique_bonds23=core.unique_bonds23,
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
    max_angles: int | None = None,
) -> dict:
    """JSON-safe payload for the app: angles + bond lengths + coordination.

    One engine pass feeding every panel of the Bond Geometry page. The dict
    (camelCase keys, plain lists/scalars) is the payload contract shared by
    the Flask ``/api/triplets`` route and the browser worker's ``triplets``
    request -- keep ``workers/triplets.js`` in sync with any change here.
    ``max_angles`` as in :func:`bond_angle_distribution`.
    """
    bin_width = _validate_bin_width(bin_width)
    core = _triplet_core(
        fractional, elements, lattice_vectors, triplet, bond12, bond23, max_angles
    )
    nbins = _bin_count(bin_width)
    histogram = _stream_angles(core.pairing, nbins, False)
    edges, density, sin_corrected = _normalized(histogram.counts, nbins)

    def length_histogram(bonds: _Bonds, window: tuple[float, float], unique: int) -> dict:
        # Clipped into the window: a bond admitted by WINDOW_TOL just outside
        # a bound belongs to the edge bin, so the histogram totals ``count``.
        length_counts, length_edges = np.histogram(
            np.clip(bonds.lengths, *window), bins=LENGTH_BINS, range=window
        )
        return {
            "binCenters": ((length_edges[:-1] + length_edges[1:]) / 2.0).tolist(),
            "counts": length_counts.tolist(),
            # B-centred bond vectors (the histogram total): an A-B bond with
            # A = B is found from both ends and appears twice here.
            "count": int(bonds.lengths.size),
            # Physical bonds, each once.
            "uniqueBonds": unique,
            "meanLength": float(np.mean(bonds.lengths)) if bonds.lengths.size else None,
        }

    coordination = np.bincount(
        np.bincount(core.bonds12.center_pos, minlength=core.apex_count)
    )

    return {
        "triplet": list(core.triplet),
        "bond12": list(core.window12),
        "bond23": list(core.window23),
        "sharedEnds": core.shared_ends,
        "binWidth": 180.0 / nbins,
        "binCenters": ((edges[:-1] + edges[1:]) / 2.0).tolist(),
        "counts": histogram.counts.tolist(),
        "density": density.tolist(),
        "sinCorrected": sin_corrected.tolist(),
        "angleCount": histogram.angle_count,
        "meanAngle": histogram.mean,
        "stdAngle": histogram.std,
        "apexCount": core.apex_count,
        "lengths12": length_histogram(core.bonds12, core.window12, core.unique_bonds12),
        # Shared ends reuse the window-1 bonds, so the page shows one histogram.
        "lengths23": (
            None
            if core.shared_ends
            else length_histogram(core.bonds23, core.window23, core.unique_bonds23)
        ),
        # coordination[n] = how many central atoms have exactly n window-1 bonds.
        "coordination": coordination.tolist(),
    }


def _read_configuration(
    path: str | Path,
) -> tuple[np.ndarray, list[str], np.ndarray, str | None]:
    """Every atom's element and supercell-fraction position, plus the parse warning.

    Bond angles need only those two, so legacy coords-only lines count too --
    the same atom set as the browser worker (``parseRmc6fAtoms`` keeps both
    layouts). Non-finite and unparsed lines are skipped by the shared grammar
    and reported through the returned ``Rmc6fParseReport.warning()`` (``None``
    for a clean atom section); no atom at all is a ``ValueError`` naming what
    was found.
    """
    lattice_vectors, _ = read_cell_vectors(path)
    coords: list[np.ndarray] = []
    elements: list[str] = []
    report = Rmc6fParseReport()
    for atom in iter_rmc6f_atoms(path, include_coords_only=True, report=report):
        coords.append(atom["coords"])
        elements.append(atom["element"])
    if not coords:
        detail = report.warning() or (
            "the Atoms section is empty" if report.has_atoms_section else "there is no Atoms section"
        )
        raise ValueError(f"{path}: no atoms could be parsed — {detail}")
    return np.asarray(coords, dtype=float), elements, lattice_vectors, report.warning()


def bond_angle_summary_from_file(
    path: str | Path,
    end1: str,
    apex: str,
    end2: str,
    r12_min: float,
    r12_max: float,
    r23_min: float,
    r23_max: float,
    bin_width: float,
    max_angles: int | None = None,
) -> dict:
    """``bond_angle_summary`` of an ``.rmc6f`` file, read afresh on every call.

    The payload of ``/api/triplets`` and of the browser worker's ``triplets``
    request, plus ``parseWarning``: the atom lines the shared ``.rmc6f``
    grammar skipped (non-finite coordinates, unparsed lines, a count short of
    the header's), or ``None``. Windows arrive resolved (bond23 defaults
    applied by the caller). Callers that cache it must key the cache on the
    file's content signature, as ``app.py``'s ``_TRIPLETS_CACHE`` does.
    """
    coords, elements, lattice_vectors, parse_warning = _read_configuration(path)
    summary = bond_angle_summary(
        coords,
        elements,
        lattice_vectors,
        triplet=(end1, apex, end2),
        bond12=(r12_min, r12_max),
        bond23=(r23_min, r23_max),
        bin_width=bin_width,
        max_angles=max_angles,
    )
    summary["parseWarning"] = parse_warning
    return summary


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
    max_angles: int | None = None,
) -> dict:
    """:func:`bond_angle_summary_from_file`, memoized for library callers.

    Keyed on (path, ``mtime``) plus every parameter, mirroring
    ``pca_kde.cached_site_displacements``. ``mtime`` is whatever the CALLER
    passes: the file is not re-stat'ed here, so a caller that passes a stale
    value (or a coarse mtime that a quick rewrite does not change) gets the
    previous configuration's result. For files that change while they are
    served (Live Data), call the uncached :func:`bond_angle_summary_from_file`
    under a cache keyed on the full file signature, as ``app.py`` does. A
    refused (over-budget) request raises and is therefore never cached.
    """
    return bond_angle_summary_from_file(
        path, end1, apex, end2, r12_min, r12_max, r23_min, r23_max, bin_width, max_angles
    )


def bond_angles_from_rmc6f(
    rmc6f_path: str | Path,
    *,
    triplet: Sequence[str],
    bond12: Sequence[float],
    bond23: Sequence[float] | None = None,
    bin_width: float = 1.0,
    collect_angles: bool = False,
    max_angles: int | None = None,
) -> BondAngleDistribution:
    """Run ``bond_angle_distribution`` on an ``.rmc6f`` configuration file.

    ``parse_warning`` of the result names the atom lines the shared grammar
    skipped (``None`` for a clean file).
    """
    rmc6f_path = Path(rmc6f_path)
    coords, elements, lattice_vectors, parse_warning = _read_configuration(rmc6f_path)
    result = bond_angle_distribution(
        coords,
        elements,
        lattice_vectors,
        triplet=triplet,
        bond12=bond12,
        bond23=bond23,
        bin_width=bin_width,
        collect_angles=collect_angles,
        max_angles=max_angles,
    )
    return replace(result, parse_warning=parse_warning)
