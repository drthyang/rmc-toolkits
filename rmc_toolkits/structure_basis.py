# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The average-structure basis of an RMC configuration: one site per reference number.

Every atom of an RMCProfile ``.rmc6f`` carries a reference number naming the site
of the unit cell it is a copy of. Folding the copies into one cell and averaging
them site by site gives the *average structure* the space-group finder reads: the
Detected SG card and its symmetry-averaged CIF (``averageStructure.js``), the AI
Copilot's ``detect_symmetry`` and ``measure_symmetry_lowering``, and the run
context's per-orbit displacements. This module is the source of truth for the
server path -- ``/api/structure`` returns :func:`structure_basis` as ``basis`` --
and ``structureFromRmc6f()`` (``web_app/frontend/src/browserData.js``) computes the
same rows for a browser-loaded run; the golden that
``tests/generate_structure_basis_fixture.py`` writes holds the two to round-off.

Every copy of a site (each atom carrying its reference number, whatever its
element) contributes, in file order, with ``x`` its fractional box coordinate and
``N`` the supercell:

* ``w = ((x N) mod 1 + 1) mod 1`` -- its fraction in one unit cell (``mod`` is C's
  ``fmod``, JavaScript's ``%``);
* ``frac`` -- per axis, the circular mean ``atan2(sum sin 2 pi w, sum cos 2 pi w) / 2 pi``
  wrapped into [0, 1) as ``a - floor(a)``, so a site straddling the cell edge
  (w ~ 0 = 1) lands on its true position rather than on 0.5;
* ``d = (w - frac) - round(w - frac)`` -- the copy's offset from that centre, wrapped
  to the nearest image (a tie rounds towards +infinity, as ``Math.round``);
* ``mean = frac + <d>`` -- the arithmetic mean of the copies unwrapped about the
  circular mean (which leans towards the mode of a skewed or split site); it may
  leave [0, 1) by a fraction of the spread;
* ``covFrac = <d d^T> - <d><d>^T`` -- the population covariance (divided by n, so a
  one-copy site has a zero covariance, not NaN) in unit-cell fractions squared,
  its diagonal floored at 0;
* ``dispA = sqrt(max(<|dr|^2> - |<dr>|^2, 0))`` with ``dr = sum_i d_i a_i`` through the
  unit-cell vectors ``a_i = L_i / N_i`` -- the rms displacement in Angstrom with
  the full metric (the square root of the Cartesian covariance's trace,
  sqrt(3 U_iso) in the PCA page's terms);
* ``count`` -- the copies of every species; ``elementCounts`` -- the copies per
  species, in the order the species first appear;
* ``el`` -- the element of the FIRST atom in the file that carries the reference
  number. A mixed or swapped site is handed to the finder as that one species
  (``elementCounts`` has them all; docs/algorithms/run-dashboard.md, Part A Step 5).

Sites come sorted by reference number, as
``{el, referenceNumber, frac, dispA, count, elementCounts, mean, covFrac}``.
Legacy coords-only atoms carry no reference number and are left out, so a file of
only those has an empty basis (and no symmetry analysis).

The sums run in file order, one copy after another, and every per-copy step is
the browser's own floating-point operation in the browser's order, so the two
runtimes differ only where their ``cos``/``sin``/``atan2`` differ (an ulp).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

import numpy as np

from .parsers import Rmc6fAtom, parse_rmc6f_atoms, read_cell_vectors

TWO_PI = 2.0 * math.pi

# The covariance's upper triangle, in the order browserData.js accumulates it.
_PAIRS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


def _round_half_up(values: np.ndarray) -> np.ndarray:
    """JavaScript's ``Math.round``: the nearest integer, a tie going towards +infinity.

    ``np.round`` sends a tie to the even neighbour, and ``floor(x + 0.5)`` is one
    too high at 0.49999999999999994 (the sum rounds up to 1); ``x - floor(x)`` is
    exact, so comparing it with 0.5 is the specification's rule.
    """
    floor = np.floor(values)
    return floor + (values - floor >= 0.5)


def _at_least_zero(value: float) -> float:
    """``Math.max(value, 0)``: NaN stays NaN; -0 and every negative value give +0."""
    return value if value > 0.0 or value != value else 0.0


def structure_basis(atoms: Iterable[Rmc6fAtom], lattice_vectors, supercell) -> list[dict]:
    """The average-structure basis of a parsed configuration (see the module doc).

    ``atoms`` are the records :func:`rmc_toolkits.parsers.parse_rmc6f_atoms` yields
    (coords-only records included or not: they are skipped either way);
    ``lattice_vectors`` the box's rows in Angstrom and ``supercell`` its tiling, as
    :func:`rmc_toolkits.parsers.read_cell_vectors` returns them. Returns one dict
    per reference number, sorted by it, with plain Python numbers (JSON-ready);
    ``[]`` when no atom carries a reference number.
    """
    referenced = [atom for atom in atoms if atom["reference_number"] is not None]
    if not referenced:
        return []
    lattice = np.asarray(lattice_vectors, dtype=float).reshape(3, 3)
    cells = np.asarray(supercell, dtype=float).reshape(3)
    references = np.fromiter(
        (atom["reference_number"] for atom in referenced), dtype=np.int64, count=len(referenced)
    )
    coords = np.array([atom["coords"] for atom in referenced], dtype=float).reshape(-1, 3)
    site_numbers, site_index = np.unique(references, return_inverse=True)
    site_index = site_index.reshape(-1)
    n_sites = len(site_numbers)

    def per_site(weights: np.ndarray) -> np.ndarray:
        # bincount adds the weights in input (file) order: the browser's running sums.
        return np.bincount(site_index, weights=weights, minlength=n_sites)

    within = np.fmod(np.fmod(coords * cells, 1.0) + 1.0, 1.0)

    frac = np.empty((n_sites, 3))
    for axis in range(3):
        angle = TWO_PI * within[:, axis]
        turn = np.arctan2(per_site(np.sin(angle)), per_site(np.cos(angle))) / TWO_PI
        frac[:, axis] = turn - np.floor(turn)

    offset = within - frac[site_index]
    offset -= _round_half_up(offset)
    unit = lattice / np.maximum(cells, 1.0)[:, None]
    # dr_k = (d0 a0k + d1 a1k) + d2 a2k, summed axis by axis as browserData.js does.
    dr = offset[:, 0:1] * unit[0]
    dr = dr + offset[:, 1:2] * unit[1]
    dr = dr + offset[:, 2:3] * unit[2]
    squared = dr[:, 0] * dr[:, 0] + dr[:, 1] * dr[:, 1] + dr[:, 2] * dr[:, 2]

    counts = np.bincount(site_index, minlength=n_sites)
    sum_dr = [per_site(dr[:, k]) for k in range(3)]
    sum_squared = per_site(squared)
    sum_d = [per_site(offset[:, i]) for i in range(3)]
    sum_dd = [per_site(offset[:, i] * offset[:, j]) for i, j in _PAIRS]

    first_element: list[str | None] = [None] * n_sites
    element_counts: list[dict[str, int]] = [{} for _ in range(n_sites)]
    for atom, site in zip(referenced, site_index.tolist()):
        element = atom["element"]
        if first_element[site] is None:
            first_element[site] = element
        tally = element_counts[site]
        tally[element] = tally.get(element, 0) + 1

    basis = []
    for site in range(n_sites):
        n = int(counts[site])
        s0, s1, s2 = (float(values[site]) for values in sum_dr)
        mean_squared = (s0 * s0 + s1 * s1 + s2 * s2) / (n * n)
        m = [float(values[site]) / n for values in sum_d]
        c = [float(sum_dd[k][site]) / n - m[i] * m[j] for k, (i, j) in enumerate(_PAIRS)]
        centre = [float(value) for value in frac[site]]
        basis.append(
            {
                "el": first_element[site],
                "referenceNumber": int(site_numbers[site]),
                "frac": centre,
                "dispA": math.sqrt(_at_least_zero(float(sum_squared[site]) / n - mean_squared)),
                "count": n,
                "elementCounts": element_counts[site],
                "mean": [centre[i] + m[i] for i in range(3)],
                "covFrac": [
                    [_at_least_zero(c[0]), c[3], c[4]],
                    [c[3], _at_least_zero(c[1]), c[5]],
                    [c[4], c[5], _at_least_zero(c[2])],
                ],
            }
        )
    return basis


def read_structure_basis(rmc6f_path: str | Path) -> list[dict]:
    """:func:`structure_basis` of an ``.rmc6f`` file.

    Reads the header with :func:`rmc_toolkits.parsers.read_cell_vectors` and the
    atoms with :func:`rmc_toolkits.parsers.parse_rmc6f_atoms` (both layouts), so
    a bad header or a file with no parseable atom raises their ``ValueError``.
    """
    lattice_vectors, supercell = read_cell_vectors(rmc6f_path)
    atoms, _ = parse_rmc6f_atoms(rmc6f_path, include_coords_only=True)
    return structure_basis(atoms, lattice_vectors, supercell)
