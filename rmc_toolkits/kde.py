"""Server-side KDE slice computation for RMC structures.

Ports the XY ``gaussian_kde`` slab math from ``src/RMC_KDE.py`` into a reusable
function that returns plain arrays/segments so a web frontend can render the
density with its own colormap and contour styling.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.stats import gaussian_kde

from .parsers import iter_rmc6f_atoms, read_cell_vectors

# Cap on the number of slab atoms fed to gaussian_kde. The density estimate is
# stable well below the full population, and the eval cost scales with the
# number of fit points, so subsampling keeps slider interaction responsive.
MAX_KDE_FIT_POINTS = 6000


@dataclass(frozen=True)
class UnitCellPositions:
    """Cartesian (Angstrom) atom positions folded into a single unit cell."""

    positions: np.ndarray  # (N, 3)
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
    for atom in iter_rmc6f_atoms(rmc6f_path):
        if select is not None and atom["element"] != select:
            continue
        unit_frac = (atom["coords"] * supercell) % 1.0
        cartesian = (
            unit_frac[0] * unit_vectors[0]
            + unit_frac[1] * unit_vectors[1]
            + unit_frac[2] * unit_vectors[2]
        )
        folded.append(cartesian)

    positions = np.asarray(folded, dtype=float) if folded else np.empty((0, 3))
    cell_lengths = np.linalg.norm(unit_vectors, axis=1)
    return UnitCellPositions(positions=positions, cell_lengths=cell_lengths)


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
    grid = int(max(16, min(grid, 400)))
    grid_x = np.linspace(xlim[0], xlim[1], grid)
    grid_y = np.linspace(ylim[0], ylim[1], grid)
    mesh_x, mesh_y = np.meshgrid(grid_x, grid_y)
    extent = [float(xlim[0]), float(xlim[1]), float(ylim[0]), float(ylim[1])]

    density = np.zeros_like(mesh_x)
    slab_count = 0
    if positions.shape[0]:
        x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]
        half = 0.5 * max(dz, 1e-12)
        mask = (z >= z_center - half) & (z <= z_center + half)
        slab = np.column_stack([x[mask], y[mask]])
        slab_count = int(slab.shape[0])

        if slab_count >= 5:
            if slab_count > MAX_KDE_FIT_POINTS:
                rng = np.random.default_rng(rng_seed)
                choice = rng.choice(slab_count, MAX_KDE_FIT_POINTS, replace=False)
                slab = slab[choice]
            kde = gaussian_kde(slab.T, bw_method=bw)
            sample = np.vstack([mesh_x.ravel(), mesh_y.ravel()])
            density = kde(sample).reshape(mesh_x.shape)

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
        "fitCount": int(min(slab_count, MAX_KDE_FIT_POINTS)) if slab_count >= 5 else 0,
        "vmin": float(np.nanmin(density)) if density.size else 0.0,
        "vmax": float(np.nanmax(density)) if density.size else 0.0,
        "contours": contours,
    }
