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
