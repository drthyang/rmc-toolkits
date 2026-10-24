# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Synthetic configurations for the modulation engine's null and diffuse test.

Shared by ``tests/test_modulation_calibration.py``, ``tests/test_modulations.py``
and ``tests/generate_modulations_fixture.py``. Each family is the physics of
one of the validation's DEV generator kinds (web_app/frontend/validation/
cases.js), with numpy's generator and seeds of its own:

* ``thermal`` -- an ABX3 perovskite with independent Gaussian displacements, X
  flattened across its bond (the ``thermal_null`` kind);
* ``jahn_teller`` -- isolated BX6 octahedra elongated along x, y or z at random
  per cell: the six X of a cell share one variable, the long axis
  (``jahn_teller``);
* ``jt_off_centre`` -- the same with every B also +-0.15 A off the centre
  along z (``jahn_teller`` with ``offCentre``, the two-feature judgment kind);
* ``octahedral`` -- M6 octahedral clusters elongated along a random axis per
  cell (``elongated_clusters``);
* ``breathing`` -- M3 triangles scaled about their centre by 1 + eps per cell
  (``breathing``);
* ``rigid_dimer`` -- M2 dimers moving as one body (``rigid_dimer``);

and two with order between cells: ``diffuse`` -- the perovskite's B sites
displaced by a field correlated between neighbouring cells along one axis,
u(c) = e(c) + s kappa (e(c + a) + e(c - a)), e independent per cell: short-range
order with no single dominant q (s = -1: antiphase, the power on the whole
sheet q_a = 1/2; s = +1: in phase); ``zone_mode`` -- one ordered zone-boundary
mode of B (a cell doubling) over thermal motion.

A configuration is built as displacements (cells x sites x 3, Angstrom) of
sites in an orthogonal cell: :func:`site_displacements` makes the engine's
input in memory (the calibration's thousands of seeds), :func:`rmc6f_text` the
``.rmc6f`` file the parity fixture carries.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from rmc_toolkits.pca_kde import SiteDisplacements

AXES = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], dtype=float)
PEROVSKITE_FRACTIONS = [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
# A and B isotropic; X flattened across its B-X bond (the axis it lies on).
THERMAL = {"A": 0.08, "B": 0.06, "X_PAR": 0.05, "X_PERP": 0.09}


@dataclass
class Configuration:
    """Sites (element, unit-cell fraction) of an orthogonal cell, the supercell
    and every atom's displacement (cells in C order x sites x 3, Angstrom)."""

    elements: list[str]
    fractions: np.ndarray
    cell: np.ndarray
    supercell: tuple[int, int, int]
    displacements: np.ndarray


def _cells(supercell) -> np.ndarray:
    return np.array([(ix, iy, iz) for ix in range(supercell[0]) for iy in range(supercell[1]) for iz in range(supercell[2])])


def _perovskite_thermal(rng, cells: int) -> np.ndarray:
    d = np.zeros((cells, 5, 3))
    d[:, 0] = rng.normal(0.0, THERMAL["A"], (cells, 3))
    d[:, 1] = rng.normal(0.0, THERMAL["B"], (cells, 3))
    for site, bond in ((2, 2), (3, 1), (4, 0)):
        x = rng.normal(0.0, THERMAL["X_PERP"], (cells, 3))
        x[:, bond] = rng.normal(0.0, THERMAL["X_PAR"], cells)
        d[:, site] = x
    return d


def thermal(seed: int, *, supercell=(8, 8, 8), elements=("Ba", "Ti", "O"), a=4.0) -> Configuration:
    rng = np.random.default_rng(seed)
    d = _perovskite_thermal(rng, int(np.prod(supercell)))
    names = [elements[0], elements[1]] + [elements[2]] * 3
    return Configuration(names, np.array(PEROVSKITE_FRACTIONS, float), np.full(3, a), tuple(supercell), d)


def jahn_teller(seed: int, *, supercell=(6, 6, 6), elements=("K", "Cu", "F"), a=8.0, short=1.93, long=2.25,
                sigma=0.05, off=None) -> Configuration:
    """BX6 octahedra (B at the body centre, A at the origin), each elongated
    along x, y or z at random; with ``off``, B also +-off along z at random."""
    rng = np.random.default_rng(seed)
    cells = int(np.prod(supercell))
    mean = (long + 2 * short) / 3.0
    d = rng.normal(0.0, sigma, (cells, 8, 3))
    axis = rng.integers(0, 3, cells)
    for k, v in enumerate(AXES):
        bond = np.where(np.abs(v[axis]) > 0, long, short)
        d[:, 2 + k] += (bond - mean)[:, None] * v[None, :]
    if off is not None:
        d[:, 1, 2] += np.where(rng.random(cells) < 0.5, off, -off)
    fractions = np.array([[0, 0, 0], [0.5, 0.5, 0.5]] + [list(0.5 + mean * v / a) for v in AXES])
    return Configuration([elements[0], elements[1]] + [elements[2]] * 6, fractions, np.full(3, a), tuple(supercell), d)


def jt_off_centre(seed: int, **options) -> Configuration:
    return jahn_teller(seed, elements=("Cs", "Ag", "F"), a=8.5, short=2.05, long=2.4, off=0.15, **options)


def octahedral(seed: int, *, supercell=(8, 8, 8), elements=("W", "Br"), a=8.5, edge=2.65, elongation=0.5,
               sigma=0.05) -> Configuration:
    """M6 octahedra about the body centre (an anion at the origin), each
    elongated along a random axis: its two apical atoms out by ``elongation``."""
    rng = np.random.default_rng(seed)
    cells = int(np.prod(supercell))
    r0 = edge / np.sqrt(2.0)
    mean = r0 + elongation / 3.0
    d = rng.normal(0.0, sigma, (cells, 7, 3))
    axis = rng.integers(0, 3, cells)
    for k, v in enumerate(AXES):
        radius = np.where(np.abs(v[axis]) > 0, r0 + elongation, r0)
        d[:, 1 + k] += (radius - mean)[:, None] * v[None, :]
    fractions = np.array([[0, 0, 0]] + [list(0.5 + mean * v / a) for v in AXES])
    return Configuration([elements[1]] + [elements[0]] * 6, fractions, np.full(3, a), tuple(supercell), d)


def breathing(seed: int, *, supercell=(8, 8, 8), elements=("Mo", "O"), a=7.0, edge=2.8, breathe=0.04,
              sigma=0.05) -> Configuration:
    """M3 triangles in the (111) plane about the body centre, each scaled by 1 + eps."""
    rng = np.random.default_rng(seed)
    cells = int(np.prod(supercell))
    directions = np.array([[2, -1, -1], [-1, 2, -1], [-1, -1, 2]], dtype=float)
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    radius = edge / np.sqrt(3.0)
    d = rng.normal(0.0, sigma, (cells, 4, 3))
    eps = rng.normal(0.0, breathe, cells)
    for j, direction in enumerate(directions):
        d[:, 1 + j] += (radius * eps)[:, None] * direction[None, :]
    fractions = np.array([[0, 0, 0]] + [list(0.5 + radius * v / a) for v in directions])
    return Configuration([elements[1]] + [elements[0]] * 3, fractions, np.full(3, a), tuple(supercell), d)


def rigid_dimer(seed: int, *, supercell=(8, 8, 8), element="V", a=5.0, bond=2.25, common=0.08, own=0.01) -> Configuration:
    rng = np.random.default_rng(seed)
    cells = int(np.prod(supercell))
    shared = rng.normal(0.0, common, (cells, 3))
    d = np.stack([shared + rng.normal(0.0, own, (cells, 3)), shared + rng.normal(0.0, own, (cells, 3))], axis=1)
    fractions = np.array([[0.25, 0.5, 0.5], [0.25 + bond / a, 0.5, 0.5]])
    return Configuration([element, element], fractions, np.full(3, a), tuple(supercell), d)


def diffuse(seed: int, *, supercell=(8, 8, 8), elements=("Ag", "Nb", "O"), a=3.95, kappa=0.3, amp=0.06, axis=0,
            sign=-1.0, anions=True) -> Configuration:
    """Short-range order between cells: the B sites' field e(c) + sign kappa
    (e(c + 1) + e(c - 1)) along ``axis`` (e ~ N(0, amp^2) per component and
    cell) over the perovskite's thermal motion. Neighbouring cells' B
    displacements correlate by 2 sign kappa amp^2 / (amp^2 (1 + 2 kappa^2) +
    B thermal variance) (:func:`diffuse_correlation`): -0.27 for the
    defaults. ``anions=False`` keeps only the A and B sites (the parity
    fixture's small file)."""
    rng = np.random.default_rng(seed)
    shape = tuple(int(n) for n in supercell)
    cells = int(np.prod(shape))
    d = _perovskite_thermal(rng, cells)
    e = rng.normal(0.0, amp, shape + (3,))
    field = e + sign * kappa * (np.roll(e, 1, axis=axis) + np.roll(e, -1, axis=axis))
    d[:, 1] += field.reshape(cells, 3)
    names = [elements[0], elements[1]] + [elements[2]] * 3
    fractions = np.array(PEROVSKITE_FRACTIONS, float)
    if not anions:
        return Configuration(names[:2], fractions[:2], np.full(3, a), shape, d[:, :2].copy())
    return Configuration(names, fractions, np.full(3, a), shape, d)


def diffuse_correlation(kappa=0.3, amp=0.06, sign=-1.0) -> float:
    """The neighbour-cell correlation of :func:`diffuse`'s B displacements
    (each component: the field's covariance at lag 1 over its variance plus
    the thermal one)."""
    return 2.0 * sign * kappa * amp ** 2 / (amp ** 2 * (1.0 + 2.0 * kappa ** 2) + THERMAL["B"] ** 2)


def zone_mode(seed: int, *, supercell=(8, 8, 8), elements=("Sr", "Ti", "O"), a=3.9, delta=0.1, q=(0, 0, 0.5)) -> Configuration:
    """B displaced delta cos(2 pi q.n) along z over the perovskite's thermal motion."""
    rng = np.random.default_rng(seed)
    cells = _cells(supercell)
    d = _perovskite_thermal(rng, len(cells))
    d[:, 1, 2] += delta * np.cos(2.0 * np.pi * cells @ np.asarray(q, dtype=float))
    names = [elements[0], elements[1]] + [elements[2]] * 3
    return Configuration(names, np.array(PEROVSKITE_FRACTIONS, float), np.full(3, a), tuple(supercell), d)


#: The families without order between cells (the calibration's nulls).
NULL_FAMILIES = {
    "thermal": thermal,
    "jahn_teller": jahn_teller,
    "jt_off_centre": jt_off_centre,
    "octahedral": octahedral,
    "breathing": breathing,
    "rigid_dimer": rigid_dimer,
}


def site_displacements(config: Configuration) -> SiteDisplacements:
    """The engine's input, built in memory: every site's displacements
    centred on their mean (as load_site_displacements centres them)."""
    supercell = np.asarray(config.supercell, dtype=int)
    cells = _cells(supercell)
    count = len(cells)
    sites = len(config.elements)
    rows = [config.displacements[:, s, :] - config.displacements[:, s, :].mean(axis=0) for s in range(sites)]
    unit = np.diag(np.asarray(config.cell, dtype=float))
    return SiteDisplacements(
        reference_numbers=np.arange(1, sites + 1),
        elements=list(config.elements),
        counts=np.full(sites, count),
        displacements=np.concatenate(rows),
        site_index=np.repeat(np.arange(sites), count),
        site_fractional=np.asarray(config.fractions, dtype=float) % 1.0,
        lattice_vectors=unit * supercell[:, None],
        supercell=supercell,
        unit_vectors=unit,
        atom_elements=np.repeat(np.array(config.elements), count),
        element_counts=tuple({name: count} for name in config.elements),
        parse_warning=None,
        cell_indices=np.tile(cells, (sites, 1)),
    )


def rmc6f_text(config: Configuration) -> str:
    """The configuration as an ``.rmc6f`` file (coordinates wrapped into [0, 1))."""
    supercell = np.asarray(config.supercell, dtype=int)
    cell = np.asarray(config.cell, dtype=float)
    lines = []
    for c, n in enumerate(_cells(supercell)):
        for s, (element, fraction) in enumerate(zip(config.elements, config.fractions)):
            coord = np.mod((np.asarray(fraction) + n + config.displacements[c, s] / cell) / supercell, 1.0)
            lines.append(f"{len(lines) + 1:6d} {element:>2} [1] {coord[0]:.12f} {coord[1]:.12f} {coord[2]:.12f}"
                         f" {s + 1:4d} {n[0]} {n[1]} {n[2]}")
    header = [
        "(Version 6f format configuration file)",
        f"Number of atoms: {len(lines)}",
        f"Supercell dimensions: {supercell[0]} {supercell[1]} {supercell[2]}",
        "Lattice vectors (Ang):",
        *(" ".join(f"{value:.6f}" for value in row) for row in np.diag(cell * supercell)),
        "Atoms:",
    ]
    return "\n".join(header + lines) + "\n"
