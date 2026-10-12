# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Synthetic ``.rmc6f`` configurations with known cluster distortions.

Shared by ``tests/test_clusters.py`` and ``tests/generate_clusters_fixture.py``
(the JS parity golden). One cubic 8 A cell holds an Nb4 tetrahedron (edge
3.0 A) centred at (1/2, 1/2, 1/2), a Ga2 dimer (bond 2.5 A along x) and a
lone Se, repeated over the supercell. Every atom gets independent Gaussian
noise; ``mode`` plants a distortion on the tetrahedra:

* ``independent`` -- nothing more: the shuffle null's own picture;
* ``trigonal`` -- in every cluster one vertex, picked at random, moves out
  along its own 3-fold axis by ``delta``: a strong local distortion whose
  average over the supercell keeps the average tetrahedron regular;
* ``static`` -- the same vertex (0) moves out in every cluster: the average
  tetrahedron is distorted and the clusters are not, beyond it;
* ``doubled`` -- every Nb moves by +delta/3 along z in even-z cells and by
  -delta/3 in odd ones: order with wavevector (0, 0, 1/2) that doubles the
  cell along c and averages away in the folded structure.
"""

from __future__ import annotations

import numpy as np

CELL_A = 8.0
TETRA_EDGE_A = 3.0
DIMER_BOND_A = 2.5
_VERTICES = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float)


def cluster_sites() -> list[tuple[str, np.ndarray]]:
    """(element, unit-cell fraction) of every reference site, in reference order."""
    half = TETRA_EDGE_A / (2.0 * np.sqrt(2.0)) / CELL_A
    sites = [("Nb", np.array([0.5, 0.5, 0.5]) + half * v) for v in _VERTICES]
    sites += [
        ("Ga", np.array([0.25 - DIMER_BOND_A / 2 / CELL_A, 0.75, 0.25])),
        ("Ga", np.array([0.25 + DIMER_BOND_A / 2 / CELL_A, 0.75, 0.25])),
        ("Se", np.array([0.0, 0.0, 0.0])),
    ]
    return sites


def cluster_config_text(mode: str = "independent", *, supercell=(4, 4, 4), sigma=0.06, delta=0.3, seed=7) -> str:
    """The ``.rmc6f`` text of one synthetic configuration (see the module doc)."""
    rng = np.random.default_rng(seed)
    supercell = np.asarray(supercell, dtype=int)
    lattice = np.diag(CELL_A * supercell.astype(float))
    sites = cluster_sites()
    lines = []
    atom = 0
    for ix in range(supercell[0]):
        for iy in range(supercell[1]):
            for iz in range(supercell[2]):
                cell = np.array([ix, iy, iz])
                if mode == "trigonal":
                    pushed = int(rng.integers(4))
                elif mode == "static":
                    pushed = 0
                else:
                    pushed = None
                for ref, (element, frac) in enumerate(sites, start=1):
                    shift = rng.normal(0.0, sigma, size=3)
                    if element == "Nb" and pushed is not None and ref - 1 == pushed:
                        shift = shift + delta * _VERTICES[pushed] / np.sqrt(3.0)
                    if element == "Nb" and mode == "doubled":
                        shift = shift + np.array([0.0, 0.0, delta / 3.0 * (1 if iz % 2 == 0 else -1)])
                    coord = (frac + cell + shift / CELL_A) / supercell
                    atom += 1
                    lines.append(
                        f"{atom:6d} {element:>2} [1] {coord[0]:.12f} {coord[1]:.12f} {coord[2]:.12f}"
                        f" {ref:4d} {ix} {iy} {iz}"
                    )
    header = [
        "(Version 6f format configuration file)",
        f"Number of atoms: {atom}",
        f"Supercell dimensions: {supercell[0]} {supercell[1]} {supercell[2]}",
        "Lattice vectors (Ang):",
        *(" ".join(f"{value:.6f}" for value in row) for row in lattice),
        "Atoms:",
    ]
    return "\n".join(header + lines) + "\n"


ORIGIN_CELL_A = 5.0
ORIGIN_SITES = [("V", np.array([0.0, 0.5, 0.5])), ("V", np.array([0.45, 0.5, 0.5]))]


def origin_dimer_config_text(bias: float, *, supercell=(4, 4, 4), sigma=0.08, seed=1) -> str:
    """A rigid V2 dimer whose first site sits at x = 0: both atoms of a dimer
    share one random shift (sd ``sigma``) plus 0.01 A of their own, and the
    first site's mean is set to ``bias`` A along x, so a negative bias puts it
    a hair below the cell edge. Coordinates are wrapped into [0, 1)."""
    rng = np.random.default_rng(seed)
    supercell = np.asarray(supercell, dtype=int)
    cells = [(ix, iy, iz) for ix in range(supercell[0]) for iy in range(supercell[1]) for iz in range(supercell[2])]
    shifts = {}
    for cell in cells:
        common = rng.normal(0.0, sigma, 3)
        for ref in (1, 2):
            shifts[(cell, ref)] = common + rng.normal(0.0, 0.01, 3)
    first = np.array([shifts[(cell, 1)] for cell in cells]).mean(axis=0)
    for cell in cells:
        shifts[(cell, 1)] = shifts[(cell, 1)] - first + np.array([bias, 0.0, 0.0])
    lines = []
    atom = 0
    for cell in cells:
        for ref, (element, frac) in enumerate(ORIGIN_SITES, start=1):
            coord = np.mod((frac + np.array(cell) + shifts[(cell, ref)] / ORIGIN_CELL_A) / supercell, 1.0)
            atom += 1
            lines.append(
                f"{atom:6d} {element:>2} [1] {coord[0]:.12f} {coord[1]:.12f} {coord[2]:.12f}"
                f" {ref:4d} {cell[0]} {cell[1]} {cell[2]}"
            )
    header = [
        "(Version 6f format configuration file)",
        f"Number of atoms: {atom}",
        f"Supercell dimensions: {supercell[0]} {supercell[1]} {supercell[2]}",
        "Lattice vectors (Ang):",
        *(" ".join(f"{value:.6f}" for value in row) for row in np.diag(ORIGIN_CELL_A * supercell.astype(float))),
        "Atoms:",
    ]
    return "\n".join(header + lines) + "\n"


CHAINS_CELL_A = 10.0
# Chain A's middle atom is its second site, chain B's its first: sorted by
# site, chain A is bonded 0-1, 1-2 and chain B 0-1, 0-2.
CHAINS_SITES = [
    ("Mo", [0.1, 0.2, 0.2]), ("Mo", [0.3, 0.2, 0.2]), ("Mo", [0.5, 0.2, 0.2]),
    ("Mo", [0.6, 0.7, 0.7]), ("Mo", [0.4, 0.7, 0.7]), ("Mo", [0.8, 0.7, 0.7]),
]


def chains_config_text(*, supercell=(4, 4, 4), sigma=0.03, seed=3) -> str:
    """Two linear Mo3 chains per cell (bonds 2.0 A along x), independent
    Gaussian noise of ``sigma`` A on every atom (the 2026-10-12 review's repro)."""
    rng = np.random.default_rng(seed)
    supercell = np.asarray(supercell, dtype=int)
    lines = []
    atom = 0
    for ix in range(supercell[0]):
        for iy in range(supercell[1]):
            for iz in range(supercell[2]):
                for ref, (element, frac) in enumerate(CHAINS_SITES, start=1):
                    shift = rng.normal(0.0, sigma, 3) / CHAINS_CELL_A
                    coord = np.mod((np.array(frac) + [ix, iy, iz] + shift) / supercell, 1.0)
                    atom += 1
                    lines.append(
                        f"{atom:6d} {element:>2} [1] {coord[0]:.12f} {coord[1]:.12f} {coord[2]:.12f}"
                        f" {ref:4d} {ix} {iy} {iz}"
                    )
    header = [
        "(Version 6f format configuration file)",
        f"Number of atoms: {atom}",
        f"Supercell dimensions: {supercell[0]} {supercell[1]} {supercell[2]}",
        "Lattice vectors (Ang):",
        *(" ".join(f"{value:.6f}" for value in row) for row in np.diag(CHAINS_CELL_A * supercell.astype(float))),
        "Atoms:",
    ]
    return "\n".join(header + lines) + "\n"


LONG_WAVE_CELL_A = 5.0
LONG_WAVE_SITES = [("Ti", [0.5, 0.5, 0.5]), ("O", [0.0, 0.5, 0.5])]


def long_wave_config_text(*, supercell=(6, 6, 6), zone_amp=0.06, long_amp=0.15, sigma=0.08, seed=5) -> str:
    """Ti displaced by +-``zone_amp`` A along z in alternate c layers (order at
    q = (0, 0, 1/2)) and by ``long_amp`` cos(2 pi x / N1) along y (a long wave
    at q = (1/N1, 0, 0), stronger), with independent noise of ``sigma`` on
    every atom (the 2026-10-12 review's repro)."""
    rng = np.random.default_rng(seed)
    supercell = np.asarray(supercell, dtype=int)
    lines = []
    atom = 0
    for ix in range(supercell[0]):
        for iy in range(supercell[1]):
            for iz in range(supercell[2]):
                for ref, (element, frac) in enumerate(LONG_WAVE_SITES, start=1):
                    shift = rng.normal(0.0, sigma, 3)
                    if element == "Ti":
                        shift[2] += zone_amp * (1 if iz % 2 == 0 else -1)
                        shift[1] += long_amp * np.cos(2 * np.pi * ix / supercell[0])
                    coord = np.mod((np.array(frac) + [ix, iy, iz] + shift / LONG_WAVE_CELL_A) / supercell, 1.0)
                    atom += 1
                    lines.append(
                        f"{atom:6d} {element:>2} [1] {coord[0]:.12f} {coord[1]:.12f} {coord[2]:.12f}"
                        f" {ref:4d} {ix} {iy} {iz}"
                    )
    header = [
        "(Version 6f format configuration file)",
        f"Number of atoms: {atom}",
        f"Supercell dimensions: {supercell[0]} {supercell[1]} {supercell[2]}",
        "Lattice vectors (Ang):",
        *(" ".join(f"{value:.6f}" for value in row) for row in np.diag(LONG_WAVE_CELL_A * supercell.astype(float))),
        "Atoms:",
    ]
    return "\n".join(header + lines) + "\n"


def with_vacancies(text: str, element: str, fraction: float, *, seed: int = 0) -> str:
    """The configuration with about ``fraction`` of ``element``'s atom lines
    dropped at random (the 2026-10-12 review's repro of missing atoms)."""
    rng = np.random.default_rng(seed)
    kept = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 10 and parts[1] == element and rng.random() < fraction:
            continue
        kept.append(line)
    return "\n".join(kept) + "\n"


TIE_SITES = [
    ("Mo", [0.1, 0.2, 0.2]), ("Mo", [0.3, 0.2, 0.2]), ("Mo", [0.2, 0.37, 0.2]),  # a triangle
    ("Mo", [0.4, 0.7, 0.7]), ("Mo", [0.6, 0.7, 0.7]),  # a dimer
]


def tie_config_text(*, supercell=(4, 4, 4), sigma=0.03, seed=3) -> str:
    """One Mo3 triangle and one Mo2 dimer per 10 A cell (bonds ~2.0 A), so
    the two kinds hold as many clusters (the 2026-10-12 review's repro)."""
    rng = np.random.default_rng(seed)
    supercell = np.asarray(supercell, dtype=int)
    lines = []
    atom = 0
    for ix in range(supercell[0]):
        for iy in range(supercell[1]):
            for iz in range(supercell[2]):
                for ref, (element, frac) in enumerate(TIE_SITES, start=1):
                    shift = rng.normal(0.0, sigma, 3) / CHAINS_CELL_A
                    coord = np.mod((np.array(frac) + [ix, iy, iz] + shift) / supercell, 1.0)
                    atom += 1
                    lines.append(
                        f"{atom:6d} {element:>2} [1] {coord[0]:.12f} {coord[1]:.12f} {coord[2]:.12f}"
                        f" {ref:4d} {ix} {iy} {iz}"
                    )
    header = [
        "(Version 6f format configuration file)",
        f"Number of atoms: {atom}",
        f"Supercell dimensions: {supercell[0]} {supercell[1]} {supercell[2]}",
        "Lattice vectors (Ang):",
        *(" ".join(f"{value:.6f}" for value in row) for row in np.diag(CHAINS_CELL_A * supercell.astype(float))),
        "Atoms:",
    ]
    return "\n".join(header + lines) + "\n"
