# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Synthetic configurations for the polyhedral-distortion engine's tests and
parity fixture (tests/test_distortion.py, tests/generate_distortion_fixture.py).

Each builder returns the atoms (``Rmc6fAtom``-like records, as
``parsers.parse_rmc6f_atoms`` yields them) with the supercell's lattice
vectors and tiling, so a test can run the engine on them directly
(``distortion_configuration``) or write the ``.rmc6f`` text both engines parse
(:func:`rmc6f_text`).
"""

from __future__ import annotations

import math

import numpy as np

AXES = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float)
TETRAHEDRON = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]], dtype=float) / math.sqrt(3.0)


def supercell_atoms(unit_cell, sites, supercell, displace=None, seed=0, centred=False):
    """Atoms of a supercell of ``unit_cell`` (3 x 3 rows, angstrom).

    ``sites``: ``[(element, fractional position in the unit cell)]``, the
    reference numbers 1, 2, ... in order; an element may be a callable
    ``(rng, cell) -> symbol`` (a mixed site). ``displace(rng, cell, site,
    element) -> Cartesian displacement`` (angstrom) moves each copy; with
    ``centred`` each site's displacements are shifted to a zero mean over its
    copies, so the average structure is the undisplaced one exactly.
    Returns ``(atoms, lattice_vectors, supercell)``.
    """
    unit_cell = np.asarray(unit_cell, dtype=float)
    tiling = np.asarray(supercell, dtype=int)
    lattice = unit_cell * tiling[:, None]
    inverse = np.linalg.inv(lattice)
    rng = np.random.default_rng(seed)
    cells = [(i, j, k) for i in range(tiling[0]) for j in range(tiling[1]) for k in range(tiling[2])]
    records = []
    for cell in cells:
        for index, (element, frac) in enumerate(sites):
            symbol = element(rng, cell) if callable(element) else element
            shift = np.zeros(3) if displace is None else np.asarray(displace(rng, cell, index, symbol), dtype=float)
            records.append([cell, index, symbol, np.asarray(frac, dtype=float), shift])
    if centred:
        for index in range(len(sites)):
            members = [record for record in records if record[1] == index]
            mean = np.mean([record[4] for record in members], axis=0)
            for record in members:
                record[4] = record[4] - mean
    atoms = []
    for number, (cell, index, symbol, frac, shift) in enumerate(records, start=1):
        box = (np.asarray(cell, dtype=float) + frac) / tiling + shift @ inverse
        atoms.append({
            "atom_number": number,
            "element": symbol,
            "type_label": symbol,
            "coords": np.mod(box, 1.0),
            "reference_number": index + 1,
            "cell_indices": np.asarray(cell, dtype=int),
        })
    return atoms, lattice, tiling


def rmc6f_text(atoms, lattice, supercell, *, coords_only=False, digits=12) -> str:
    """The ``.rmc6f`` text of ``atoms`` (the full layout, or legacy coords-only lines)."""
    lines = [
        "(Version 6f format configuration file)",
        f"Number of atoms:                    {len(atoms)}",
        f"Supercell dimensions:               {supercell[0]} {supercell[1]} {supercell[2]}",
        "Lattice vectors (Ang):",
        *("  ".join(f"{value:.12f}" for value in row) for row in np.asarray(lattice, dtype=float)),
        "Atoms:",
    ]
    for atom in atoms:
        x, y, z = (f"{value:.{digits}f}" for value in atom["coords"])
        if coords_only:
            lines.append(f"{atom['atom_number']} {atom['element']} {x} {y} {z}")
        else:
            cx, cy, cz = (int(value) for value in atom["cell_indices"])
            lines.append(
                f"{atom['atom_number']} {atom['element']} [{atom['reference_number']}] {x} {y} {z} "
                f"{atom['reference_number']} {cx} {cy} {cz}"
            )
    return "\n".join(lines) + "\n"


def isolated_polyhedron(vectors, *, a=10.0, supercell=(2, 2, 2), centre="Cu", ligand="F", offset=(0.0, 0.0, 0.0)):
    """One polyhedron per cube of edge ``a``: the centre at the cell's middle,
    a ligand at each Cartesian vector from it (angstrom). ``offset`` moves the
    centre off it (the ligands stay)."""
    middle = np.full(3, 0.5)
    sites = [(centre, middle + np.asarray(offset, dtype=float) / a)]
    sites += [(ligand, middle + np.asarray(vector, dtype=float) / a) for vector in vectors]
    return supercell_atoms(np.eye(3) * a, sites, supercell)


def octahedron(lengths=(2.0,) * 6):
    """Bond vectors along ±x, ±y, ±z with the given lengths (in that order)."""
    return [AXES[k] * lengths[k] for k in range(6)]


def trigonal_octahedron(psi_deg: float, radius: float = 2.0):
    """Six bonds of one length at polar angle psi from [0 0 1], azimuths 0, 120,
    240 deg and their inverses: a regular octahedron at psi = arccos(1/sqrt(3))
    = 54.7356 deg, compressed along its threefold axis above it."""
    psi = math.radians(psi_deg)
    top = [radius * np.array([math.sin(psi) * math.cos(phi), math.sin(psi) * math.sin(phi), math.cos(psi)])
           for phi in (0.0, 2 * math.pi / 3, 4 * math.pi / 3)]
    return top + [-vector for vector in top]


def flattened_tetrahedron(cos_alpha: float, radius: float = 2.0):
    """Four bonds of one length, (±s, 0, c) and (0, ±s, -c) (D2d): regular at
    c = 1/sqrt(3), flattened along [0 0 1] below it."""
    c = cos_alpha
    s = math.sqrt(1.0 - c * c)
    return [radius * np.array(v) for v in ((s, 0, c), (-s, 0, c), (0, s, -c), (0, -s, -c))]


def jahn_teller_network(short=1.93, long=2.25, a=8.0, supercell=(3, 3, 3), sigma=0.0, seed=7):
    """Isolated CuF6 octahedra, each 4 x short + 2 x long along the axis
    (i + j + k) mod 3 of its cell, so every F site is long in exactly a third
    of the cells: the average octahedron is regular at (long + 2 short) / 3.
    ``sigma`` adds isotropic Gaussian noise (angstrom per axis) to every atom."""
    mean = (long + 2 * short) / 3

    def displace(rng, cell, site, element):
        noise = rng.normal(scale=sigma, size=3) if sigma else np.zeros(3)
        if site == 0:
            return noise
        axis = AXES[site - 1]
        bond = long if axis[sum(cell) % 3] != 0 else short
        return noise + axis * (bond - mean)

    sites = [("Cu", (0.5, 0.5, 0.5))] + [("F", 0.5 + AXES[k] * mean / a) for k in range(6)]
    return supercell_atoms(np.eye(3) * a, sites, supercell, displace=displace, seed=seed)


PEROVSKITE_SITES = [
    ("Sr", (0.0, 0.0, 0.0)), ("Ti", (0.5, 0.5, 0.5)),
    ("O", (0.5, 0.5, 0.0)), ("O", (0.5, 0.0, 0.5)), ("O", (0.0, 0.5, 0.5)),
]


def perovskite(a=3.9, supercell=(4, 4, 4), sigma_o=0.0, sigma_cation=0.0, seed=11, centred=False, unit_cell=None, b_site=None):
    """Cubic ABO3 (Sr at the corner, Ti at the centre, O at the face centres):
    corner-sharing TiO6 octahedra. ``sigma_o`` / ``sigma_cation``: isotropic
    Gaussian displacements (angstrom per axis) of the O and of the cations;
    ``unit_cell`` replaces the cubic cell (a triclinic one shears the
    octahedra); ``b_site`` replaces the B site's element (a callable: a
    mixed site)."""
    def displace(rng, cell, site, element):
        sigma = sigma_o if site >= 2 else sigma_cation
        return rng.normal(scale=sigma, size=3) if sigma else np.zeros(3)

    sites = list(PEROVSKITE_SITES)
    if b_site is not None:
        sites[1] = (b_site, sites[1][1])
    cell = np.eye(3) * a if unit_cell is None else np.asarray(unit_cell, dtype=float)
    return supercell_atoms(cell, sites, supercell, displace=displace, seed=seed, centred=centred)


def zincblende(a=5.4, supercell=(3, 3, 3), sigma=0.0, seed=13):
    """Zincblende ZnS: corner-sharing ZnS4 and SZn4 tetrahedra."""
    fcc = [(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)]
    sites = [("Zn", f) for f in fcc] + [("S", tuple(np.add(f, 0.25))) for f in fcc]

    def displace(rng, cell, site, element):
        return rng.normal(scale=sigma, size=3) if sigma else np.zeros(3)

    return supercell_atoms(np.eye(3) * a, sites, supercell, displace=displace, seed=seed)
