# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Server-side KDE slice computation for RMC structures.

Ports the XY ``gaussian_kde`` slab math from ``src/RMC_KDE.py`` into a reusable
function that returns plain arrays/segments so a web frontend can render the
density with its own colormap and contour styling.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path

import numpy as np
from scipy.linalg import cholesky, solve_triangular
from scipy.stats import gaussian_kde

from .parsers import iter_rmc6f_atoms, read_cell_vectors

# Cap on the number of slab atoms fed to gaussian_kde. The density estimate is
# stable well below the full population, and the eval cost scales with the
# number of fit points, so subsampling keeps slider interaction responsive.
MAX_KDE_FIT_POINTS = 6000

# Why a slab produced no density. The browser worker
# (web_app/frontend/src/workers/localKdeWorker.js, KDE_MESSAGES) returns the
# same strings for the same conditions, checked in the same order, so both
# runtimes either draw the same kernel or decline with the same reason.
KDE_MESSAGES = {
    "empty": "No atoms in this slab.",
    "bandwidth": "The bandwidth must be a positive finite number.",
    "too_few": "Fewer than 5 slab rows: too few atoms for a 2D KDE.",
    "few_unique": (
        "The slab atoms occupy fewer than 3 distinct in-plane positions, "
        "so their covariance (the KDE bandwidth) is undefined."
    ),
    "collinear": (
        "The slab atoms are collinear in the slice plane, "
        "so their covariance (the KDE bandwidth) is singular."
    ),
    "singular": (
        "The slab covariance is singular to within round-off, "
        "so the KDE bandwidth is undefined."
    ),
}

# Slab membership is |d - z_c| <= dz/2 + SLAB_FACE_TOLERANCE, with d the atom's
# depth normalised to [0, 1] across the unit cube's projection range (the
# slider's units). Atoms of an ideal or unrelaxed configuration sit exactly on
# slider-reachable faces (z = 0.125 against z_c = 0.165, dz = 0.08), where two
# roundings of the same inequality disagree and a whole site flips in or out.
# The worker (localKdeWorker.js) and the Slab-In-Cell highlight
# (StructurePage.jsx) use the same expression and constant (slabSelection.js).
SLAB_FACE_TOLERANCE = 1e-9

# A covariance whose in-plane correlation coefficient rho satisfies
# 1 - rho^2 <= this limit is declined as numerically singular. Below it the
# Cholesky pivot of C sits at the level of summation round-off (~N*eps), so
# whether it comes out positive depends on the summation order, and the two
# runtimes could disagree on whether (and how) to draw. The corresponding
# kernel would be a needle with an aspect ratio above ~2e5, far below any
# grid spacing. localKdeWorker.js applies the same test.
COVARIANCE_CONDITION_LIMIT = 1e-10


def _well_conditioned(covariance: np.ndarray) -> bool:
    """True when ``covariance`` is safely positive definite (see COVARIANCE_CONDITION_LIMIT)."""
    c00, c01, c11 = float(covariance[0, 0]), float(covariance[0, 1]), float(covariance[1, 1])
    if not (c00 > 0.0 and c11 > 0.0):
        return False
    return 1.0 - (c01 * c01) / (c00 * c11) > COVARIANCE_CONDITION_LIMIT


# Warnings attached to a drawn map (``warnings`` in the payload, as
# {"code", "message"}); the worker returns the same codes and strings.
KERNEL_SUBGRID_RATIO = 0.5
KDE_WARNINGS = {
    "subgrid": (
        "The kernel is narrower than half a grid step along its minor axis, so the map is "
        "aliased: peak values, contours and the integrated density depend on the grid size. "
        "Raise the bandwidth or the grid."
    ),
}


def _kernel_warnings(kernel: dict, grid_step: float) -> list[dict]:
    """``subgrid`` when the minor kernel sigma is below ``KERNEL_SUBGRID_RATIO`` grid steps.

    A Gaussian sampled at spacing h keeps its integral to ~1 % while
    sigma >= h/2; below that most atoms fall between nodes and the few that land
    on one spike, so the drawn map is set by the grid, not the atoms.
    """
    if kernel["sigmaMinor"] < KERNEL_SUBGRID_RATIO * grid_step:
        return [{"code": "subgrid", "message": KDE_WARNINGS["subgrid"]}]
    return []


def _valid_bandwidth(bw) -> bool:
    """A usable bandwidth factor: a finite number > 0 (booleans excluded)."""
    if isinstance(bw, bool) or not isinstance(bw, (int, float, np.integer, np.floating)):
        return False
    return math.isfinite(float(bw)) and float(bw) > 0.0


def _kernel_summary(covariance: np.ndarray, cholesky_factor: np.ndarray) -> dict:
    """Kernel matrix H and its principal standard deviations (in-plane units).

    ``cholesky_factor`` is the lower Cholesky factor of ``covariance``. The minor
    eigenvalue is taken as det(H)/lambda_max with det(H) from the Cholesky
    diagonal, which stays accurate for needle-shaped kernels where the
    closed-form ``mean - radius`` would cancel. The worker computes the same
    closed form (``kernelSummary`` in localKdeWorker.js).
    """
    h00, h01, h11 = float(covariance[0, 0]), float(covariance[0, 1]), float(covariance[1, 1])
    root_det = float(cholesky_factor[0, 0]) * float(cholesky_factor[1, 1])
    lambda_major = 0.5 * (h00 + h11) + math.hypot(0.5 * (h00 - h11), h01)
    lambda_minor = root_det * root_det / lambda_major if lambda_major > 0 else 0.0
    return {
        "covariance": [[h00, h01], [h01, h11]],
        "sigmaMinor": math.sqrt(max(lambda_minor, 0.0)),
        "sigmaMajor": math.sqrt(max(lambda_major, 0.0)),
    }


_CUBE_CORNERS = np.asarray(
    [[float(x), float(y), float(z)] for x in (0, 1) for y in (0, 1) for z in (0, 1)],
    dtype=float,
)
_CUBE_EDGES = [
    (start, end)
    for start, end in combinations(range(len(_CUBE_CORNERS)), 2)
    if np.count_nonzero(_CUBE_CORNERS[start] != _CUBE_CORNERS[end]) == 1
]


@dataclass(frozen=True)
class UnitCellPositions:
    """Cartesian (Angstrom) atom positions folded into a single unit cell."""

    positions: np.ndarray  # (N, 3)
    fractional_positions: np.ndarray  # (N, 3)
    unit_vectors: np.ndarray  # (3, 3)
    cell_lengths: np.ndarray  # (3,) unit-cell edge lengths


def load_unit_cell_positions(
    rmc6f_path: str | Path,
    element: str | None = None,
) -> UnitCellPositions:
    """Load atom positions from an ``.rmc6f`` file folded into one unit cell.

    Coordinates are returned in Angstrom (cartesian), matching the desktop
    ``RMC_KDE.py`` convention so axis limits and aspect ratios line up.
    """
    rmc6f_path = Path(rmc6f_path)
    lattice_vectors, supercell = read_cell_vectors(rmc6f_path)
    unit_vectors = lattice_vectors / supercell[:, None]

    select = element if element not in (None, "", "all") else None
    folded: list[np.ndarray] = []
    fractional: list[np.ndarray] = []
    for atom in iter_rmc6f_atoms(rmc6f_path):
        if select is not None and atom["element"] != select:
            continue
        unit_frac = (atom["coords"] * supercell) % 1.0
        fractional.append(unit_frac)
        cartesian = (
            unit_frac[0] * unit_vectors[0]
            + unit_frac[1] * unit_vectors[1]
            + unit_frac[2] * unit_vectors[2]
        )
        folded.append(cartesian)

    positions = np.asarray(folded, dtype=float) if folded else np.empty((0, 3))
    fractional_positions = np.asarray(fractional, dtype=float) if fractional else np.empty((0, 3))
    cell_lengths = np.linalg.norm(unit_vectors, axis=1)
    return UnitCellPositions(
        positions=positions,
        fractional_positions=fractional_positions,
        unit_vectors=unit_vectors,
        cell_lengths=cell_lengths,
    )


def _contour_segments(
    grid_x: np.ndarray,
    grid_y: np.ndarray,
    density: np.ndarray,
    n_levels: int,
) -> list[dict]:
    """Extract contour polylines for a density field without rendering a figure.

    ``density`` may be log-transformed, so no sign test here: whether there is
    anything to contour is decided on the linear density by the caller.
    """
    if n_levels <= 0 or not np.isfinite(density).any():
        return []

    # contourpy ships with matplotlib; use it directly to avoid pyplot state.
    from contourpy import contour_generator

    finite_max = float(np.nanmax(density))
    finite_min = float(np.nanmin(density))
    if finite_max <= finite_min:
        return []

    levels = np.linspace(finite_min, finite_max, n_levels + 2)[1:-1]
    generator = contour_generator(grid_x, grid_y, density)
    segments: list[dict] = []
    for level in levels:
        lines = generator.lines(float(level))
        polylines = [line.tolist() for line in lines if len(line) >= 2]
        if polylines:
            segments.append({"level": float(level), "lines": polylines})
    return segments


def _image_rank(offset: tuple[float, float, float]) -> int:
    """Order of preference among an atom's periodic images: fewest, then earliest shifts.

    ``|offset|^2 * 27`` plus the offset's position in ``product((-1, 0, 1), repeat=3)``
    (the loop order of both runtimes), so the unshifted atom ranks first among
    its images and every rank is unique. ``imageRank`` in localKdeWorker.js.
    """
    ox, oy, oz = (int(value) for value in offset)
    return (ox * ox + oy * oy + oz * oz) * 27 + (ox + 1) * 9 + (oy + 1) * 3 + (oz + 1)


def _augment_periodic_images(
    positions: np.ndarray,
    margin: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tile fractional positions into neighbor cells within ``margin`` of the cube.

    Folding atoms into one unit cell drops their periodic neighbors, so a KDE
    evaluated near a cell face misses the density that should wrap around from
    the opposite face. Adding the images restores those contributions. Returns
    the augmented positions plus, for every row, the index of its source atom
    (so callers can count unique atoms and normalize out the image duplicates)
    and its image rank (``_image_rank``; lower is a smaller shift), which picks
    one representative row per source atom for the bandwidth.
    """
    positions = np.asarray(positions, dtype=float)
    n = positions.shape[0]
    source_index = np.arange(n)
    identity_rank = np.full(n, _image_rank((0.0, 0.0, 0.0)))
    if n == 0 or margin <= 0:
        return positions, source_index, identity_rank

    augmented = [positions]
    sources = [source_index]
    ranks = [identity_rank]
    for offset in product((-1.0, 0.0, 1.0), repeat=3):
        if offset == (0.0, 0.0, 0.0):
            continue
        shifted = positions + np.asarray(offset)
        inside = np.all((shifted >= -margin) & (shifted <= 1.0 + margin), axis=1)
        if inside.any():
            augmented.append(shifted[inside])
            sources.append(source_index[inside])
            ranks.append(np.full(int(inside.sum()), _image_rank(offset)))
    return np.concatenate(augmented), np.concatenate(sources), np.concatenate(ranks)


def _dot3(points: np.ndarray, vector: np.ndarray) -> np.ndarray:
    """Row-wise ``points . vector`` summed left to right, like the worker's ``dot()``.

    An explicit ``(x*v0 + y*v1) + z*v2`` rather than a BLAS matmul, so depths
    and in-plane coordinates round exactly as in the browser.
    """
    points = np.atleast_2d(points)
    return points[:, 0] * vector[0] + points[:, 1] * vector[1] + points[:, 2] * vector[2]


def _normalize_vector(vector: np.ndarray, name: str) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(vector))
    if vector.shape != (3,) or norm <= 1e-12:
        raise ValueError(f"{name} must be a non-zero 3D vector")
    return vector / norm


def _orthogonal_axis(vector: np.ndarray) -> np.ndarray:
    axis = np.eye(3)[int(np.argmin(np.abs(vector)))]
    candidate = axis - np.dot(axis, vector) * vector
    return _normalize_vector(candidate, "plane axis")


def _plane_basis(
    normal: np.ndarray,
    u_axis: np.ndarray | None = None,
    v_axis: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    normal = _normalize_vector(normal, "normal")
    if u_axis is None:
        u_axis = _orthogonal_axis(normal)
    else:
        u_axis = np.asarray(u_axis, dtype=float)
        u_axis = u_axis - np.dot(u_axis, normal) * normal
        u_axis = _normalize_vector(u_axis, "u_axis")

    if v_axis is None:
        v_axis = np.cross(normal, u_axis)
    else:
        v_axis = np.asarray(v_axis, dtype=float)
        v_axis = v_axis - np.dot(v_axis, normal) * normal - np.dot(v_axis, u_axis) * u_axis
    v_axis = _normalize_vector(v_axis, "v_axis")
    return normal, u_axis, v_axis


def _plane_section_vertices(normal: np.ndarray, offset: float) -> list[list[float]]:
    vertices: list[np.ndarray] = []
    for start, end in _CUBE_EDGES:
        p0 = _CUBE_CORNERS[start]
        p1 = _CUBE_CORNERS[end]
        d0 = float(np.dot(p0, normal) - offset)
        d1 = float(np.dot(p1, normal) - offset)
        if abs(d0) <= 1e-9:
            vertices.append(p0)
        if abs(d1) <= 1e-9:
            vertices.append(p1)
        if d0 * d1 < 0:
            t = d0 / (d0 - d1)
            vertices.append(p0 + t * (p1 - p0))

    unique: list[np.ndarray] = []
    for vertex in vertices:
        if not any(np.linalg.norm(vertex - existing) <= 1e-8 for existing in unique):
            unique.append(vertex)
    if len(unique) < 3:
        return []

    _, u_axis, v_axis = _plane_basis(normal)
    center = np.mean(unique, axis=0)
    ordered = sorted(
        unique,
        key=lambda vertex: np.arctan2(np.dot(vertex - center, v_axis), np.dot(vertex - center, u_axis)),
    )
    return [vertex.tolist() for vertex in ordered]


def oriented_kde_slice(
    positions: np.ndarray,
    center: float,
    thickness: float,
    *,
    normal: np.ndarray,
    u_axis: np.ndarray | None = None,
    v_axis: np.ndarray | None = None,
    bw: float = 0.03,
    grid: int = 120,
    log: bool = False,
    n_levels: int = 8,
    rng_seed: int = 0,
) -> dict:
    """Compute a KDE slice through fractional coordinates along any direction.

    ``center`` and ``thickness`` are fractions of the unit-cube projection range
    along ``normal`` (they equal cell-edge fractions only for the a/b/c axis
    normals; nothing is converted to Angstrom). For example, normal
    ``[0, 0, 1]`` matches the original c-axis slice semantics. An atom is in
    the slab when its normalised depth d satisfies
    ``|d - center| <= thickness / 2 + SLAB_FACE_TOLERANCE``; the payload's
    ``z``/``dz`` echo ``center``/``thickness`` and ``depth``/``depthThickness``
    give the same slab in absolute depth units.

    Positions are treated as periodic: images from neighbor cells within a
    margin of the unit cube join the slab selection and the density sum, so the
    density is correct at cell faces, edges, and corners instead of decaying
    toward the boundary. The kernel itself is fitted to the slab's source
    atoms (one row each), so the margin does not change it (see kde_slice).
    """
    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions must be a numeric array with shape (N, 3)")

    # Margin must cover the kernel reach (sigma scales with bw times the data
    # spread, which is O(1) in fractional units) and the slab depth, so both
    # the in-plane density and the depth selection wrap correctly. An invalid
    # bandwidth contributes nothing here; kde_slice declines it.
    bw_reach = float(bw) if _valid_bandwidth(bw) else 0.0
    margin = min(0.5, max(0.1, 2.0 * bw_reach, thickness))
    positions, source_index, image_rank = _augment_periodic_images(positions, margin)

    normal, u_axis, v_axis = _plane_basis(normal, u_axis, v_axis)
    corner_depths = _dot3(_CUBE_CORNERS, normal)
    depth_min = float(np.min(corner_depths))
    depth_max = float(np.max(corner_depths))
    depth_span = max(depth_max - depth_min, 1e-12)
    center = max(0.0, min(float(center), 1.0))
    thickness = max(float(thickness), 1e-12)
    center_depth = depth_min + center * depth_span
    thickness_depth = thickness * depth_span

    projected_corners = np.column_stack([_dot3(_CUBE_CORNERS, u_axis), _dot3(_CUBE_CORNERS, v_axis)])
    xlim = (float(np.min(projected_corners[:, 0])), float(np.max(projected_corners[:, 0])))
    ylim = (float(np.min(projected_corners[:, 1])), float(np.max(projected_corners[:, 1])))

    # The slab test runs on the depth normalised across the cube's projection
    # range -- the slider's own units, so z/dz in the payload echo the inputs --
    # with the expression the worker and the Slab-In-Cell highlight use.
    normalized_depth = (_dot3(positions, normal) - depth_min) / depth_span
    projected_positions = np.column_stack(
        [_dot3(positions, u_axis), _dot3(positions, v_axis), normalized_depth]
    )
    result = kde_slice(
        projected_positions,
        z_center=center,
        dz=thickness,
        xlim=xlim,
        ylim=ylim,
        bw=bw,
        grid=grid,
        log=log,
        n_levels=n_levels,
        rng_seed=rng_seed,
        source_index=source_index,
        image_rank=image_rank,
    )

    slab_start = max(depth_min, center_depth - thickness_depth / 2)
    slab_end = min(depth_max, center_depth + thickness_depth / 2)
    plane_vertices = _plane_section_vertices(normal, center_depth)
    result.update(
        {
            "center": center,
            "thickness": thickness,
            "depth": float(center_depth),
            "depthThickness": float(thickness_depth),
            "depthRange": [depth_min, depth_max],
            "normal": normal.tolist(),
            "uVector": u_axis.tolist(),
            "vVector": v_axis.tolist(),
            "planeVertices": plane_vertices,
            "planePolygon": [
                [float(np.dot(vertex, u_axis)), float(np.dot(vertex, v_axis))]
                for vertex in np.asarray(plane_vertices, dtype=float)
            ]
            if plane_vertices
            else [],
            "slabVertices": [
                _plane_section_vertices(normal, slab_start),
                _plane_section_vertices(normal, slab_end),
            ],
        }
    )
    return result


def _source_atom_rows(
    slab: np.ndarray,
    mask: np.ndarray,
    source_index: np.ndarray | None,
    image_rank: np.ndarray | None,
) -> np.ndarray:
    """One in-plane row per source atom in the slab: the bandwidth's data.

    Among an atom's slab rows the one with the lowest image rank (its smallest
    periodic shift) represents it; ties cannot occur because ranks are unique
    per atom. Rows come out ordered by source index, which is also the order
    the worker's makeSlab() collects them in.
    """
    if source_index is None:
        return slab
    sources = np.asarray(source_index)[mask]
    if sources.size == 0:
        return slab
    ranks = np.asarray(image_rank)[mask] if image_rank is not None else np.flatnonzero(mask)
    order = np.lexsort((ranks, sources))
    ordered_sources = sources[order]
    first = np.ones(order.size, dtype=bool)
    first[1:] = ordered_sources[1:] != ordered_sources[:-1]
    return slab[order[first]]


class _FixedCovarianceKDE(gaussian_kde):
    """``scipy.stats.gaussian_kde`` with a supplied data covariance.

    gaussian_kde estimates ``C`` from its own dataset. Here the dataset is the
    slab rows -- periodic images included, possibly subsampled -- while ``C``
    must come from the source atoms alone, so that the kernel ``bw**2 * C``
    depends on neither the image margin nor the subsample. Only the covariance
    estimate is replaced: ``_compute_covariance`` sets the attributes scipy's
    own ``_compute_covariance`` sets (``factor``, ``covariance``, ``cho_cov``,
    ``log_det``), and the density is still scipy's compiled Gaussian sum.
    Because that relies on scipy's internals, construction checks one value
    against the direct formula and raises rather than return a different
    kernel on a scipy version that evaluates differently.
    """

    def __init__(self, dataset: np.ndarray, data_covariance: np.ndarray, bw: float):
        self._fixed_covariance = np.atleast_2d(np.asarray(data_covariance, dtype=float))
        super().__init__(dataset, bw_method=bw)
        self._check_scipy_honours_the_covariance()

    def _compute_covariance(self):
        self.factor = self.covariance_factor()
        self._data_covariance = self._fixed_covariance
        self._data_cho_cov = cholesky(self._data_covariance, lower=True)
        self.covariance = self._data_covariance * self.factor**2
        self.cho_cov = (self._data_cho_cov * self.factor).astype(np.float64)
        self.log_det = 2 * np.log(np.diag(self.cho_cov * np.sqrt(2 * np.pi))).sum()

    def _check_scipy_honours_the_covariance(self) -> None:
        point = self.dataset[:, :1]
        whitened = solve_triangular(self.cho_cov, self.dataset - point, lower=True)
        expected = float(
            np.exp(-0.5 * np.sum(whitened * whitened, axis=0)).sum()
            / (self.n * 2 * np.pi * self.cho_cov[0, 0] * self.cho_cov[1, 1])
        )
        actual = float(self.evaluate(point)[0])
        # scipy whitens the point and the data separately, so needle kernels
        # differ from this direct form at ~1e-10; a scipy that ignored the
        # supplied covariance would differ at O(1).
        if not abs(actual - expected) <= 1e-6 * expected:
            raise RuntimeError(
                "scipy.stats.gaussian_kde no longer evaluates a supplied covariance "
                f"(scipy internals changed: {actual!r} != {expected!r}); update rmc_toolkits.kde"
            )


def kde_slice(
    positions: np.ndarray,
    z_center: float,
    dz: float,
    *,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    bw: float = 0.03,
    grid: int = 120,
    log: bool = False,
    n_levels: int = 8,
    rng_seed: int = 0,
    source_index: np.ndarray | None = None,
    image_rank: np.ndarray | None = None,
) -> dict:
    """Compute an XY ``gaussian_kde`` density for a z-slab of a structure.

    Returns a JSON-serializable dict with the density grid, plot extent,
    contour polylines, the slab atom count, the kernel, and (when no density
    was drawn) a ``message`` saying why. A row is in the slab when
    ``|z - z_center| <= dz / 2 + SLAB_FACE_TOLERANCE`` (an absolute tolerance
    in the units of ``z``; oriented_kde_slice passes normalised depths).

    The kernel is scipy's: ``H = bw**2 * C``. ``C`` is the covariance of the
    slab's *source atoms*, one row per atom, and the density sums that fixed
    kernel over every slab row (periodic images included, subsampled to
    ``MAX_KDE_FIT_POINTS``), so neither the image margin nor the subsample
    changes the kernel.

    ``source_index`` maps each position row to its source atom when the input
    contains periodic images. The reported ``slabCount`` is then the number of
    unique atoms in the slab, and the density is rescaled to per-atom
    normalization (``gaussian_kde`` divides by every fit point, images included).
    ``image_rank`` orders a source atom's rows (lower = smaller periodic shift,
    see ``_augment_periodic_images``); the lowest-ranked slab row represents
    the atom in ``C``. Without it the first slab row of each atom does; without
    ``source_index`` every row is its own atom.
    """
    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions must be a numeric array with shape (N, 3)")

    grid = int(max(16, min(grid, 400)))
    grid_x = np.linspace(xlim[0], xlim[1], grid)
    grid_y = np.linspace(ylim[0], ylim[1], grid)
    mesh_x, mesh_y = np.meshgrid(grid_x, grid_y)
    extent = [float(xlim[0]), float(xlim[1]), float(ylim[0]), float(ylim[1])]

    density = np.zeros_like(mesh_x)
    slab_count = 0
    fit_count = 0
    kernel = None
    warnings: list[dict] = []
    message = KDE_MESSAGES["empty"]
    if positions.shape[0]:
        x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
        half = 0.5 * max(dz, 1e-12)
        mask = np.abs(z - z_center) <= half + SLAB_FACE_TOLERANCE
        slab = np.column_stack([x[mask], y[mask]])
        slab_total = int(slab.shape[0])
        atoms = _source_atom_rows(slab, mask, source_index, image_rank)
        slab_count = int(atoms.shape[0])

        # Decline reasons, in the order the browser worker checks them.
        if slab_total == 0:
            message = KDE_MESSAGES["empty"]
        elif not _valid_bandwidth(bw):
            message = KDE_MESSAGES["bandwidth"]
        elif slab_total < 5:
            message = KDE_MESSAGES["too_few"]
        elif np.unique(atoms, axis=0).shape[0] < 3:
            message = KDE_MESSAGES["few_unique"]
        elif np.linalg.matrix_rank(atoms - atoms.mean(axis=0)) < 2:
            message = KDE_MESSAGES["collinear"]
        else:
            covariance = np.cov(atoms, rowvar=False)
            if slab_total > MAX_KDE_FIT_POINTS:
                rng = np.random.default_rng(rng_seed)
                choice = rng.choice(slab_total, MAX_KDE_FIT_POINTS, replace=False)
                slab = slab[choice]
            try:
                if not _well_conditioned(covariance):
                    raise np.linalg.LinAlgError("slab covariance is numerically singular")
                kde = _FixedCovarianceKDE(slab.T, covariance, float(bw))
            except (np.linalg.LinAlgError, ValueError):
                # The source atoms' covariance is not (safely) positive
                # definite even though they passed the rank test.
                message = KDE_MESSAGES["singular"]
            else:
                sample = np.vstack([mesh_x.ravel(), mesh_y.ravel()])
                density = kde(sample).reshape(mesh_x.shape)
                if slab_total > slab_count > 0:
                    # gaussian_kde divides by every fit point, periodic images
                    # included; rescale to per-source-atom normalization so the
                    # amplitude matches the cell interior.
                    density *= slab_total / slab_count
                fit_count = int(slab.shape[0])
                kernel = _kernel_summary(kde.covariance, kde.cho_cov)
                grid_step = max(
                    (float(xlim[1]) - float(xlim[0])) / (grid - 1),
                    (float(ylim[1]) - float(ylim[0])) / (grid - 1),
                )
                warnings = _kernel_warnings(kernel, grid_step)
                message = None

    # Contour only a map with positive density, tested on the linear values:
    # after log10 a smooth field whose peak is <= 1 per unit fractional area
    # (common for oblique slices at bw >= 0.1) has a negative maximum and must
    # still be contoured. The worker's guard is equivalent (its vmax > vmin is
    # false only for the all-zero grid of a declined slab).
    has_density = bool(np.nanmax(density) > 0) if density.size else False
    if log:
        density = np.log10(density + 1e-12)

    contours = _contour_segments(grid_x, grid_y, density, n_levels) if has_density else []

    return {
        "density": density.tolist(),
        "extent": extent,
        "grid": grid,
        "z": float(z_center),
        "dz": float(dz),
        "bw": float(bw) if _valid_bandwidth(bw) else None,
        "log": bool(log),
        "slabCount": slab_count,
        "fitCount": fit_count,
        # H = bw^2 * Cov (in-plane units of `positions`) and its principal
        # sigmas; None when the slab was declined, and then `message` says why.
        "kernel": kernel,
        "message": message,
        "warnings": warnings,
        "vmin": float(np.nanmin(density)) if density.size else 0.0,
        "vmax": float(np.nanmax(density)) if density.size else 0.0,
        "contours": contours,
    }
