"""Server-side KDE slice computation for RMC structures.

Ports the XY ``gaussian_kde`` slab math from ``src/RMC_KDE.py`` into a reusable
function that returns plain arrays/segments so a web frontend can render the
density with its own colormap and contour styling.
"""

from __future__ import annotations

from contextlib import suppress
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import gaussian_kde

from .parsers import iter_rmc6f_atoms, read_cell_vectors

# Cap on the number of slab atoms fed to gaussian_kde. The density estimate is
# stable well below the full population, and the eval cost scales with the
# number of fit points, so subsampling keeps slider interaction responsive.
MAX_KDE_FIT_POINTS = 6000
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
    """Extract contour polylines for a density field without rendering a figure."""
    if n_levels <= 0 or not np.isfinite(density).any() or float(density.max()) <= 0:
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
    along ``normal``. For example, normal ``[0, 0, 1]`` matches the original
    c-axis slice semantics.
    """
    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3:
        raise ValueError("positions must be a numeric array with shape (N, 3)")

    normal, u_axis, v_axis = _plane_basis(normal, u_axis, v_axis)
    corner_depths = _CUBE_CORNERS @ normal
    depth_min = float(np.min(corner_depths))
    depth_max = float(np.max(corner_depths))
    depth_span = max(depth_max - depth_min, 1e-12)
    center = max(0.0, min(float(center), 1.0))
    thickness = max(float(thickness), 1e-12)
    center_depth = depth_min + center * depth_span
    thickness_depth = thickness * depth_span

    projected_corners = np.column_stack([_CUBE_CORNERS @ u_axis, _CUBE_CORNERS @ v_axis])
    xlim = (float(np.min(projected_corners[:, 0])), float(np.max(projected_corners[:, 0])))
    ylim = (float(np.min(projected_corners[:, 1])), float(np.max(projected_corners[:, 1])))

    projected_positions = np.column_stack([positions @ u_axis, positions @ v_axis, positions @ normal])
    result = kde_slice(
        projected_positions,
        z_center=center_depth,
        dz=thickness_depth,
        xlim=xlim,
        ylim=ylim,
        bw=bw,
        grid=grid,
        log=log,
        n_levels=n_levels,
        rng_seed=rng_seed,
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
) -> dict:
    """Compute an XY ``gaussian_kde`` density for a z-slab of a structure.

    Returns a JSON-serializable dict with the density grid, plot extent,
    contour polylines, and the slab atom count.
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
    if positions.shape[0]:
        x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
        half = 0.5 * max(dz, 1e-12)
        mask = (z >= z_center - half) & (z <= z_center + half)
        slab = np.column_stack([x[mask], y[mask]])
        slab_count = int(slab.shape[0])

        if slab_count >= 5:
            has_enough_unique_points = np.unique(slab, axis=0).shape[0] >= 3
            centered_slab = slab - slab.mean(axis=0)
            has_two_dimensional_spread = np.linalg.matrix_rank(centered_slab) >= 2
        else:
            has_enough_unique_points = False
            has_two_dimensional_spread = False

        if has_enough_unique_points and has_two_dimensional_spread:
            if slab_count > MAX_KDE_FIT_POINTS:
                rng = np.random.default_rng(rng_seed)
                choice = rng.choice(slab_count, MAX_KDE_FIT_POINTS, replace=False)
                slab = slab[choice]
            with suppress(np.linalg.LinAlgError, ValueError):
                kde = gaussian_kde(slab.T, bw_method=bw)
                sample = np.vstack([mesh_x.ravel(), mesh_y.ravel()])
                density = kde(sample).reshape(mesh_x.shape)
                fit_count = int(slab.shape[0])

    if log:
        density = np.log10(density + 1e-12)

    contours = _contour_segments(grid_x, grid_y, density, n_levels)

    return {
        "density": density.tolist(),
        "extent": extent,
        "grid": grid,
        "z": float(z_center),
        "dz": float(dz),
        "bw": float(bw),
        "log": bool(log),
        "slabCount": slab_count,
        "fitCount": fit_count,
        "vmin": float(np.nanmin(density)) if density.size else 0.0,
        "vmax": float(np.nanmax(density)) if density.size else 0.0,
        "contours": contours,
    }
