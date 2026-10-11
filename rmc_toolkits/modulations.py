# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Displacement modulations: is there order that enlarges the cell?

The average structure folds every unit cell of the supercell into one, so it
cannot show order whose period is longer than one cell -- the cell doubling of
GaNb4Se8 below 29 K (P2_12_12_1, c' = 2c) or of GaTa4Se8 below 50 K, octahedral
tilts, charge-density waves. The configuration itself can: an ordered
displacement pattern with wavevector q puts power into the Fourier component
of every site's displacement field at q.

For each reference site s with copies in the cells c of an N1 x N2 x N3
supercell, the displacements u_s(c) (from the site mean, Angstrom) are Fourier
transformed over the cell grid,

    U_s(q) = sum_c u_s(c) exp(-2 pi i q.c),   q = (h/N1, k/N2, l/N3),

and the power P_s(q) = |U_s(q)|^2 / (n_s tr C_s) is normalised by the copy
count and the site's total displacement variance, so independent (white)
displacements give 1 on average at every q != 0. The powers of an element's
sites are averaged. Two kinds of q are reported against the empirical
distribution over all q != 0 (the configuration's own diffuse background,
which carries RMC's long-wavelength correlations that a white-noise null
ignores):

* the zone-boundary points (1/2,0,0) ... (1/2,1/2,1/2) of the cell's
  reciprocal lattice -- order doubling the cell along those directions -- when
  the supercell is even along each halved axis;
* the strongest modes overall (+q and -q counted once).

Each with its ratio to the median power and its percentile (the share of
q != 0 with lower power). The browser port is ``workers/modulations.js``
(a separable DFT, the same numbers to floating-point noise).
"""

from __future__ import annotations

from itertools import product

import numpy as np

from rmc_toolkits.pca_kde import SiteDisplacements

#: How many of the strongest modes to report per element.
TOP_MODES = 5
#: Relative tolerance on "lower power" for the percentile: Friedel pairs
#: P(q) = P(-q) differ only by rounding, and must count alike in both engines.
PERCENTILE_RTOL = 1e-9


def zone_boundary_points(supercell: np.ndarray) -> list[tuple[int, int, int]]:
    """The half-integer points (as 0/1 flags per axis) the supercell can hold."""
    points = []
    for flags in product((0, 1), repeat=3):
        if not any(flags):
            continue
        if all(supercell[axis] % 2 == 0 for axis in range(3) if flags[axis]):
            points.append(flags)
    return points


def _folded(index: int, size: int) -> float:
    """Grid index -> fractional wavevector component in (-1/2, 1/2]."""
    value = int(index) / int(size)
    return value - 1.0 if value > 0.5 else value


def displacement_modulations(sites: SiteDisplacements, element: str | None = None, *, top: int = TOP_MODES) -> dict:
    """Per-element displacement power at every wavevector of the supercell grid.

    ``element`` restricts the analysis to one element; None analyses each
    element present. Raises ``ValueError`` for an unknown element or a file
    without cell indices.
    """
    if sites.cell_indices is None:
        raise ValueError("modulation analysis needs the reference-site and cell columns of the .rmc6f file")
    supercell = np.asarray(sites.supercell, dtype=int)
    n_q = int(np.prod(supercell))
    if n_q < 8:
        raise ValueError(f"a {supercell.tolist()} supercell has too few wavevectors to compare")
    labels = list(sites.elements)
    if element is not None:
        element = element.strip().capitalize()
        if element not in labels:
            raise ValueError(f"no site holds {element}; elements: {', '.join(sites.species)}")
        elements = [element]
    else:
        elements = sorted(set(labels))
    cells = np.mod(np.asarray(sites.cell_indices, dtype=int), supercell)
    boundary = zone_boundary_points(supercell)
    result = {
        "supercell": supercell.tolist(),
        "wavevectors": n_q - 1,
        "untestable": [list(flags) for flags in product((0, 1), repeat=3)
                       if any(flags) and list(flags) not in [list(b) for b in boundary]],
        "elements": [],
        "parseWarning": sites.parse_warning,
    }
    for name in elements:
        members = [s for s, label in enumerate(labels) if label == name]
        power = np.zeros(tuple(supercell))
        missing = 0
        used = 0
        for s in members:
            rows = np.flatnonzero(sites.site_index == s)
            u = sites.displacements[rows]
            variance = float(((u - u.mean(axis=0)) ** 2).mean(axis=0).sum())
            if not variance > 0:
                continue
            grid = np.zeros(tuple(supercell) + (3,))
            grid[cells[rows, 0], cells[rows, 1], cells[rows, 2]] = u
            missing += n_q - rows.size
            transform = np.fft.fftn(grid, axes=(0, 1, 2))
            power += (np.abs(transform) ** 2).sum(axis=-1) / (rows.size * variance)
            used += 1
        if not used:
            result["elements"].append({"element": name, "sites": 0})
            continue
        power /= used
        flat = power.ravel()
        others = np.delete(flat, 0)
        median = float(np.median(others))

        def describe(index: tuple[int, int, int]) -> dict:
            value = float(power[index])
            return {
                "q": [float(_folded(index[axis], supercell[axis])) for axis in range(3)],
                "ratio": value / median if median > 0 else None,
                "percentile": float((others < value * (1.0 - PERCENTILE_RTOL)).mean() * 100.0),
            }

        zone = []
        for flags in boundary:
            index = tuple(supercell[axis] // 2 if flags[axis] else 0 for axis in range(3))
            entry = describe(index)
            entry["q"] = [0.5 * flag for flag in flags]
            zone.append(entry)
        order = np.argsort(-flat, kind="stable")
        strongest = []
        seen = set()
        for linear in order:
            if linear == 0:
                continue
            index = tuple(int(v) for v in np.unravel_index(int(linear), tuple(supercell)))
            partner = tuple(int((-index[axis]) % supercell[axis]) for axis in range(3))
            # +q and -q carry the same power (to rounding): report the pair
            # once, as its member with the smaller grid index, in both engines.
            key = min(index, partner)
            if key in seen:
                continue
            seen.add(key)
            strongest.append(describe(key))
            if len(strongest) == top:
                break
        result["elements"].append({
            "element": name,
            "sites": used,
            "missingCopies": int(missing),
            "medianPower": median,
            "zoneBoundary": zone,
            "strongest": strongest,
        })
    return result
