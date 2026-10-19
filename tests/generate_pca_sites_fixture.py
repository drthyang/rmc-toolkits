# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the JS parity golden for the per-site displacement statistics.

Writes ``web_app/frontend/src/__tests__/fixtures/pca_sites_fixture.json``:
``.rmc6f`` configurations and the rows ``rmc_toolkits.pca_kde.site_ellipsoids``
makes of them (the table behind ``/api/pca/sites`` and the PCA Ellipsoid and
Displacement Directions pages) -- the one written interface that
``workers/pcaKde.js`` (``siteEllipsoids``) is held to
(``src/__tests__/pcaSitesParity.test.js``). ``tests/test_pca_sites_fixture.py``
fails when the committed file no longer matches the Python engine.

Cases: every site of the committed demo run (``web_app/frontend/public/demo``,
52 sites of 1000 copies), and a synthetic configuration with five sites of
256 copies each -- isotropic, split (two wells, negative kurtosis), zero spread
(every copy exactly on its site), a 60/40 mixed-occupancy prolate site, and a
50/50 mixed site with heavy tails (positive kurtosis, labelled by the
alphabetically first species) -- at two ellipsoid probabilities. Every cloud
stays below ``MAX_PCA_FIT_POINTS`` (20 000), above which the engines subsample
differently.

Regenerate from the checkout under test whenever the engine's conventions
change (without ``PYTHONPATH=.`` another checkout's ``rmc_toolkits`` may be
imported):

    PYTHONPATH=. python tests/generate_pca_sites_fixture.py
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np

from rmc_toolkits.pca_kde import MAX_PCA_FIT_POINTS, load_site_displacements, site_ellipsoids

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "pca_sites_fixture.json"
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"
DEMO_FILE = "GTS_250K.rmc6f"

CELL_A = 6.0
# Powers of two: a copy exactly on its site, (n + s) / N with s in 1/4ths, is
# exact in binary and in the file, so the zero-spread site has no round-off.
SUPERCELL = (8, 8, 4)


def _rotation(angles) -> np.ndarray:
    """A fixed rotation (z-y-x Euler angles, radians): oblique principal axes."""
    a, b, c = angles
    rz = np.array([[np.cos(a), -np.sin(a), 0.0], [np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])
    ry = np.array([[np.cos(b), 0.0, np.sin(b)], [0.0, 1.0, 0.0], [-np.sin(b), 0.0, np.cos(b)]])
    rx = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(c), -np.sin(c)], [0.0, np.sin(c), np.cos(c)]])
    return rz @ ry @ rx


def _isotropic(rng, n):
    return rng.normal(size=(n, 3)) * 0.08


def _split(rng, n):
    # Two wells 0.7 A apart along x: a bimodal axis, excess kurtosis near -2.
    wells = np.where(rng.random(n) < 0.5, -0.35, 0.35)
    return np.column_stack([wells, np.zeros(n), np.zeros(n)]) + rng.normal(size=(n, 3)) * [0.03, 0.04, 0.05]


def _zero(rng, n):
    return np.zeros((n, 3))


def _prolate(rng, n):
    return (rng.normal(size=(n, 3)) * [0.10, 0.06, 0.04]) @ _rotation((0.4, -0.7, 1.1)).T


def _heavy_tailed(rng, n):
    return (rng.laplace(size=(n, 3)) * [0.05, 0.04, 0.035]) @ _rotation((-0.9, 0.3, 0.5)).T


# (reference number, species and copies in cell order, unit-cell fraction, model)
SITES = [
    (1, [("Ta", 256)], (0.13, 0.27, 0.41), _isotropic),
    (2, [("Se", 256)], (0.62, 0.55, 0.31), _split),
    (3, [("Ga", 256)], (0.25, 0.75, 0.5), _zero),
    (4, [("Nb", 154), ("Ta", 102)], (0.81, 0.18, 0.57), _prolate),
    # Near the box edge along c (s = 0.06): copies wrap past 0 and are unwrapped.
    (5, [("Te", 128), ("Se", 128)], (0.44, 0.93, 0.06), _heavy_tailed),
]


def synthetic_config_text(seed: int = 2026) -> str:
    """The ``.rmc6f`` text of the synthetic configuration (see the module doc)."""
    rng = np.random.default_rng(seed)
    supercell = np.asarray(SUPERCELL, dtype=int)
    lattice = np.diag(CELL_A * supercell.astype(float))
    cells = [(ix, iy, iz) for ix in range(supercell[0]) for iy in range(supercell[1]) for iz in range(supercell[2])]
    lines = []
    atom = 0
    for reference, species, fraction, model in SITES:
        elements = [element for element, copies in species for _ in range(copies)]
        order = rng.permutation(len(cells))  # which cells get which species
        shifts = model(rng, len(cells))
        for slot, index in enumerate(order):
            cell = np.asarray(cells[index])
            coord = np.mod((np.asarray(fraction) + cell + shifts[slot] / CELL_A) / supercell, 1.0)
            atom += 1
            lines.append(
                f"{atom:6d} {elements[slot]:>2} [{reference}] {coord[0]:.12f} {coord[1]:.12f} {coord[2]:.12f}"
                f" {reference:4d} {cell[0]} {cell[1]} {cell[2]}"
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


def _case(name: str, source: dict, sites, probability: float) -> dict:
    if int(sites.counts.max()) >= MAX_PCA_FIT_POINTS:
        raise ValueError(f"{name}: a site holds {int(sites.counts.max())} copies, not below {MAX_PCA_FIT_POINTS}")
    return {
        "name": name,
        "source": source,
        "probability": probability,
        "referenceNumbers": sites.reference_numbers.tolist(),
        "species": sites.species,
        "parseWarning": sites.parse_warning,
        "expected": site_ellipsoids(sites, probability=probability),
    }


def build_fixture(synthetic_text: str | None = None) -> dict:
    """The whole fixture as a JSON-ready dict (``main`` writes it).

    ``synthetic_text`` replaces :func:`synthetic_config_text` -- the staleness
    test passes the committed text, so it checks the engine, not NumPy's
    random stream.
    """
    if synthetic_text is None:
        synthetic_text = synthetic_config_text()
    demo = load_site_displacements(DEMO / DEMO_FILE)
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "synthetic.rmc6f"
        path.write_text(synthetic_text, encoding="utf-8")
        synthetic = load_site_displacements(path)
    return {
        "synthetic": synthetic_text,
        "cases": [
            _case("demo run (GaTa4Se8, 250 K)", {"demo": DEMO_FILE}, demo, 0.5),
            _case("synthetic sites, p = 0.5", {"synthetic": True}, synthetic, 0.5),
            _case("synthetic sites, p = 0.95", {"synthetic": True}, synthetic, 0.95),
        ],
    }


def main() -> None:
    text = json.dumps(build_fixture(), separators=(",", ":"), allow_nan=False)
    OUT.write_text(text + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    main()
