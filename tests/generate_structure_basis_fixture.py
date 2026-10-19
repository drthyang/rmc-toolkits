# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the JS parity golden for the average-structure basis.

Writes ``web_app/frontend/src/__tests__/fixtures/structure_basis_fixture.json``:
``.rmc6f`` configurations and the basis ``rmc_toolkits.structure_basis`` makes
of them -- the ``basis`` of ``/api/structure``, which the Detected SG card, its
CIF export and the AI Copilot's ``detect_symmetry`` and
``measure_symmetry_lowering`` read for a server folder.
``structureFromRmc6f()`` (``browserData.js``) is held to it
(``src/__tests__/structureBasisParity.test.js``), and the frontend's Flask-mode
tests use the demo run's rows as the server's payload.
``tests/test_structure_basis_fixture.py`` fails when the committed file no
longer matches the Python engine.

Cases:

* every site of the committed demo run (``web_app/frontend/public/demo``,
  GaTa4Se8 250 K, 52 sites of 1000 copies);
* a synthetic configuration in an oblique (triclinic) cell, 4 x 4 x 2 cells, its
  atom lines shuffled so no site's copies are contiguous: a site at x ~ 0 and
  z ~ 1 whose copies straddle the cell edge, a 20/12 Nb/Ta mixed site whose
  first line in the file is the minority Ta (``el`` is the first atom's
  element), a split site with unequal wells (its arithmetic ``mean`` leaves
  the circular-mean ``frac``), a correlated anisotropic site (off-diagonal
  ``covFrac``), a site whose copies all sit exactly on it, a one-copy site
  (zero covariance), reference numbers out of order, and three coords-only
  atoms that the basis leaves out;
* a file of coords-only atoms only: an empty basis.

The configurations' text is committed in the fixture, so the staleness test
checks the engine, not NumPy's random stream. Regenerate from the checkout
under test whenever the basis's conventions change (without ``PYTHONPATH=.``
another checkout's ``rmc_toolkits`` may be imported):

    PYTHONPATH=. python tests/generate_structure_basis_fixture.py
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np

from rmc_toolkits.structure_basis import read_structure_basis

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "structure_basis_fixture.json"
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"
DEMO_FILE = "GTS_250K.rmc6f"

SUPERCELL = (4, 4, 2)
# Unit-cell rows (Angstrom): an oblique cell, so dispA needs the full metric.
UNIT_CELL = ((6.2, 0.0, 0.0), (-1.9, 5.8, 0.0), (0.7, 1.2, 7.1))


def _gaussian(sd):
    return lambda rng, n: rng.normal(size=(n, 3)) * sd


def _still(rng, n):
    return np.zeros((n, 3))


def _split(rng, n):
    # Two wells 0.06 apart along x, 22 copies in one and 10 in the other.
    wells = np.where(np.arange(n) < 22, -0.03, 0.03)
    return np.column_stack([wells, np.zeros(n), np.zeros(n)]) + rng.normal(size=(n, 3)) * 0.004


def _correlated(rng, n):
    mixing = np.array([[0.011, 0.0, 0.0], [0.006, 0.008, 0.0], [-0.004, 0.003, 0.005]])
    return rng.normal(size=(n, 3)) @ mixing.T


# (reference number, species and copies, unit-cell fraction, displacement model).
# Reference numbers are deliberately not in file order and not contiguous.
SITES = [
    (3, [("Ga", 32)], (0.004, 0.37, 0.996), _gaussian(0.012)),
    (1, [("Nb", 20), ("Ta", 12)], (0.31, 0.62, 0.18), _gaussian(0.008)),
    (7, [("Se", 32)], (0.61, 0.52, 0.28), _split),
    (4, [("Se", 32)], (0.82, 0.21, 0.57), _correlated),
    (12, [("Ga", 32)], (0.5, 0.25, 0.75), _still),
    (9, [("Te", 1)], (0.83, 0.14, 0.66), _gaussian(0.01)),
]
COORDS_ONLY = [("O", (0.11, 0.42, 0.93)), ("O", (0.71, 0.05, 0.38)), ("O", (0.47, 0.88, 0.02))]


def _header(count: int) -> list[str]:
    lattice = np.asarray(UNIT_CELL) * np.asarray(SUPERCELL, dtype=float)[:, None]
    return [
        "(Version 6f format configuration file)",
        f"Number of atoms: {count}",
        f"Supercell dimensions: {SUPERCELL[0]} {SUPERCELL[1]} {SUPERCELL[2]}",
        "Lattice vectors (Ang):",
        *(" ".join(f"{value:.6f}" for value in row) for row in lattice),
        "Atoms:",
    ]


def synthetic_config_text(seed: int = 2026) -> str:
    """The ``.rmc6f`` text of the synthetic configuration (see the module doc)."""
    rng = np.random.default_rng(seed)
    supercell = np.asarray(SUPERCELL, dtype=float)
    cells = [(ix, iy, iz) for ix in range(SUPERCELL[0]) for iy in range(SUPERCELL[1]) for iz in range(SUPERCELL[2])]
    records = []  # (reference, element, coord, cell)
    for reference, species, fraction, model in SITES:
        elements = [element for element, copies in species for _ in range(copies)]
        chosen = rng.permutation(len(cells))[: len(elements)]
        shifts = model(rng, len(elements))
        for slot, index in enumerate(chosen):
            cell = np.asarray(cells[index])
            coord = np.mod((np.asarray(fraction) + shifts[slot] + cell) / supercell, 1.0)
            records.append((reference, elements[slot], coord, cell))
    for element, fraction in COORDS_ONLY:
        records.append((None, element, np.asarray(fraction) / supercell, None))
    order = rng.permutation(len(records))
    records = [records[index] for index in order]
    # The mixed site's first line in the file is its minority species.
    mixed = [index for index, record in enumerate(records) if record[0] == 1]
    first_ta = next(index for index in mixed if records[index][1] == "Ta")
    records[mixed[0]], records[first_ta] = records[first_ta], records[mixed[0]]
    lines = []
    for number, (reference, element, coord, cell) in enumerate(records, start=1):
        position = f"{coord[0]:.15f} {coord[1]:.15f} {coord[2]:.15f}"
        if reference is None:
            lines.append(f"{number:6d} {element:>2} {position}")
        else:
            lines.append(
                f"{number:6d} {element:>2} [{reference}] {position} {reference:4d} {cell[0]} {cell[1]} {cell[2]}"
            )
    return "\n".join(_header(len(records)) + lines) + "\n"


def coords_only_config_text() -> str:
    """Eight coords-only atoms (no reference numbers): no basis at all."""
    supercell = np.asarray(SUPERCELL, dtype=float)
    lines = []
    for number in range(1, 9):
        coord = np.mod((np.array([0.1, 0.3, 0.7]) * number) / supercell, 1.0)
        lines.append(f"{number:6d} Ga {coord[0]:.15f} {coord[1]:.15f} {coord[2]:.15f}")
    return "\n".join(_header(len(lines)) + lines) + "\n"


def _basis_of_text(text: str) -> list[dict]:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "config.rmc6f"
        path.write_text(text, encoding="utf-8")
        return read_structure_basis(path)


def build_fixture(synthetic_text: str | None = None, coords_only_text: str | None = None) -> dict:
    """The whole fixture as a JSON-ready dict (``main`` writes it).

    The two texts replace :func:`synthetic_config_text` and
    :func:`coords_only_config_text` -- the staleness test passes the committed
    ones, so it checks the engine, not NumPy's random stream.
    """
    if synthetic_text is None:
        synthetic_text = synthetic_config_text()
    if coords_only_text is None:
        coords_only_text = coords_only_config_text()
    return {
        "synthetic": synthetic_text,
        "coordsOnly": coords_only_text,
        "cases": [
            {
                "name": "demo run (GaTa4Se8, 250 K)",
                "source": {"demo": DEMO_FILE},
                "basis": read_structure_basis(DEMO / DEMO_FILE),
            },
            {
                "name": "synthetic: edge-straddling, mixed, split, correlated, still and one-copy sites",
                "source": {"synthetic": True},
                "basis": _basis_of_text(synthetic_text),
            },
            {
                "name": "coords-only file: no reference numbers, no basis",
                "source": {"coordsOnly": True},
                "basis": _basis_of_text(coords_only_text),
            },
        ],
    }


def main() -> None:
    text = json.dumps(build_fixture(), separators=(",", ":"), allow_nan=False)
    OUT.write_text(text + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    main()
