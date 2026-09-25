# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the JS parity golden for the triplets engine.

Writes ``web_app/frontend/src/__tests__/fixtures/triplets_fixture.json``:
deterministic input configurations plus the reference-grade Python results
(`rmc_toolkits.triplets.bond_angle_distribution`). The vitest suite runs the
`workers/triplets.js` port on the same inputs and asserts exact integer
histograms and float agreement, so regenerate this file whenever the engine's
conventions change:

    python tests/generate_triplets_fixture.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from rmc_toolkits.triplets import (
    APP_MAX_ANGLES,
    EDGE_SNAP_DEG,
    bond_angle_distribution,
    bond_angle_summary,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "triplets_fixture.json"


def case_random_triclinic() -> dict:
    rng = np.random.default_rng(11)
    count = 48
    positions = rng.uniform(size=(count, 3))
    elements = ["Se" if value < 0.6 else "Nb" for value in rng.uniform(size=count)]
    lattice = [[6.0, 0.0, 0.0], [3.0, 5.0, 0.0], [1.0, 1.0, 7.0]]
    return {
        "name": "random-triclinic",
        "fractional": np.round(positions, 12).tolist(),
        "elements": elements,
        "lattice": lattice,
        "specs": [
            {"triplet": ["Se", "Nb", "Se"], "bond12": [1.0, 3.4], "bond23": None, "binWidth": 5.0},
            {"triplet": ["Se", "Nb", "Nb"], "bond12": [1.0, 3.0], "bond23": [1.5, 3.4], "binWidth": 5.0},
            {"triplet": ["Se", "Se", "Se"], "bond12": [1.0, 3.2], "bond23": None, "binWidth": 2.0},
            # Same end element, overlapping distinct windows: each physical
            # triplet counts once under either window assignment.
            {"triplet": ["Se", "Nb", "Se"], "bond12": [1.0, 3.0], "bond23": [1.5, 3.4], "binWidth": 5.0},
            {"triplet": ["Se", "Se", "Se"], "bond12": [1.0, 3.0], "bond23": [2.0, 3.4], "binWidth": 3.0},
        ],
    }


def case_small_box_images() -> dict:
    return {
        "name": "small-box-images",
        "fractional": [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0], [0.5, 0.55, 0.1]],
        "elements": ["Nb", "Se", "Se"],
        "lattice": [[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]],
        "specs": [
            {"triplet": ["Se", "Nb", "Se"], "bond12": [0.5, 3.5], "bond23": None, "binWidth": 1.0}
        ],
    }


def case_ideal_perovskite() -> dict:
    """Undisplaced cubic SrTiO3 3x3x3 (fractions (i + x) / 3, as a CIF-built
    start configuration stores them): TiO6 octahedra and SrO12 cages whose
    60/90/120/180 deg angles all sit on bin edges up to float noise -- the
    engines must bin each symmetry class whole and identically."""
    cells, a = 3, 3.905
    basis = [("Sr", (0, 0, 0)), ("Ti", (0.5, 0.5, 0.5)), ("O", (0.5, 0.5, 0)),
             ("O", (0.5, 0, 0.5)), ("O", (0, 0.5, 0.5))]
    fractional, elements = [], []
    for i in range(cells):
        for j in range(cells):
            for k in range(cells):
                for element, (x, y, z) in basis:
                    fractional.append([(i + x) / cells, (j + y) / cells, (k + z) / cells])
                    elements.append(element)
    return {
        "name": "ideal-perovskite",
        "fractional": fractional,
        "elements": elements,
        "lattice": [[a * cells, 0.0, 0.0], [0.0, a * cells, 0.0], [0.0, 0.0, a * cells]],
        "specs": [
            {"triplet": ["O", "Ti", "O"], "bond12": [1.5, 2.5], "bond23": None, "binWidth": 1.0},
            {"triplet": ["O", "Sr", "O"], "bond12": [2.0, 3.0], "bond23": None, "binWidth": 1.0},
            {"triplet": ["O", "O", "O"], "bond12": [2.0, 3.0], "bond23": None, "binWidth": 0.5},
            {"triplet": ["Sr", "Ti", "O"], "bond12": [3.0, 3.6], "bond23": [1.5, 2.5], "binWidth": 5.0},
        ],
    }


def case_ideal_fcc() -> dict:
    """Undisplaced fcc Cu 3x3x3 conventional cells: 60/90/120/180 deg only.
    The case where libm and V8 acos used to put every 60 deg angle in a
    different bin (59 vs 60)."""
    cells, a = 3, 3.61
    basis = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
    fractional = [
        [(i + x) / cells, (j + y) / cells, (k + z) / cells]
        for i in range(cells) for j in range(cells) for k in range(cells)
        for x, y, z in basis
    ]
    return {
        "name": "ideal-fcc",
        "fractional": fractional,
        "elements": ["Cu"] * len(fractional),
        "lattice": [[a * cells, 0.0, 0.0], [0.0, a * cells, 0.0], [0.0, 0.0, a * cells]],
        "specs": [
            {"triplet": ["Cu", "Cu", "Cu"], "bond12": [2.3, 2.8], "bond23": None, "binWidth": 1.0},
            {"triplet": ["Cu", "Cu", "Cu"], "bond12": [2.3, 2.8], "bond23": None, "binWidth": 3.0},
        ],
    }


def angle_samples(sorted_angles: np.ndarray, limit: int = 200, ends: int = 50) -> dict:
    values = [round(float(v), 6) for v in sorted_angles]
    if len(values) <= limit:
        return {"sortedAngles": values}
    return {"sortedAnglesHead": values[:ends], "sortedAnglesTail": values[-ends:]}


def evaluate(case: dict) -> dict:
    expectations = []
    for spec in case["specs"]:
        result = bond_angle_distribution(
            np.asarray(case["fractional"], dtype=float),
            case["elements"],
            np.asarray(case["lattice"], dtype=float),
            triplet=spec["triplet"],
            bond12=spec["bond12"],
            bond23=spec["bond23"],
            bin_width=spec["binWidth"],
            collect_angles=True,
        )
        def rounded(values, digits):
            return [None if v is None else round(float(v), digits) for v in values]

        expectations.append(
            {
                "counts": result.counts.tolist(),
                "density": rounded(result.density, 12),
                "sinCorrected": rounded(result.sin_corrected, 10),
                "angleCount": result.angle_count,
                "meanAngle": rounded([result.mean_angle], 8)[0],
                "stdAngle": rounded([result.std_angle], 8)[0],
                "apexCount": result.apex_count,
                "bond12Count": result.bond12_count,
                "bond23Count": result.bond23_count,
                "meanLength12": rounded([result.mean_length12], 8)[0],
                "meanLength23": rounded([result.mean_length23], 8)[0],
                # 1e-6 deg rounding: far coarser than float noise, far finer
                # than the vitest tolerance, and it keeps the fixture small.
                # Exact bin counts pin the histogram, so big cases only need
                # spot checks at the sorted list's ends (plus mean/std above).
                **angle_samples(np.sort(result.angles)),
                # The app payload: lengths + coordination parity for the port.
                "summary": bond_angle_summary(
                    np.asarray(case["fractional"], dtype=float),
                    case["elements"],
                    np.asarray(case["lattice"], dtype=float),
                    triplet=spec["triplet"],
                    bond12=spec["bond12"],
                    bond23=spec["bond23"],
                    bin_width=spec["binWidth"],
                ),
            }
        )
    return {**case, "expected": expectations}


def main() -> None:
    payload = {
        "note": "Generated by tests/generate_triplets_fixture.py -- do not edit by hand.",
        # The app-boundary work budget; the JS port's APP_MAX_ANGLES must match.
        "appMaxAngles": APP_MAX_ANGLES,
        # The bin-edge snap; the JS port's EDGE_SNAP_DEG must match.
        "edgeSnapDeg": EDGE_SNAP_DEG,
        "cases": [
            evaluate(case_random_triclinic()),
            evaluate(case_small_box_images()),
            evaluate(case_ideal_perovskite()),
            evaluate(case_ideal_fcc()),
        ],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload) + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
