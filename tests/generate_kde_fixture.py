# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the JS parity golden for the Structure-page KDE slice.

Writes ``web_app/frontend/src/__tests__/fixtures/kde_parity_fixture.json``: the
reference-grade Python densities (`rmc_toolkits.kde.oriented_kde_slice`) for a
set of slices on real configurations and on small synthetic point sets. The
vitest suite (``workers/__tests__/kdeParity.test.js``) runs the browser worker
(`computeKde`) on the same points and asserts agreement to 1e-6 of the peak,
identical slab/fit counts, identical kernels and identical decline messages
(and their ``messageCode``).
``tests/test_kde_parity_fixture.py`` checks that the committed golden is still
what the Python engine produces. Regenerate whenever the engine changes:

    python tests/generate_kde_fixture.py

Real-data cases point at repo-relative ``.rmc6f`` files. The bundled demo run
(``web_app/frontend/public/demo/GTS_250K.rmc6f``) is committed; the GaNb4Se8
``data/5K_try1`` run is not (``data/`` is gitignored), so its cases are marked
``requiresData`` and both test suites skip them when the file is absent. Every
case keeps the slab below the 6000-point fit cap, so both runtimes sum the same
rows (the two engines draw different subsamples above it).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from rmc_toolkits.kde import MAX_KDE_FIT_POINTS, load_unit_cell_positions, oriented_kde_slice

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "kde_parity_fixture.json"
DEMO = "web_app/frontend/public/demo/GTS_250K.rmc6f"
DATA_5K = "data/5K_try1/GaNb4Se8_5K.rmc6f"
GRID = 32

PRESETS = {
    "a": ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]),
    "b": ([0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
    "c": ([0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
}


def browser_plane_basis(normal: list[float]) -> tuple[list[float], list[float], list[float]]:
    """The frame StructurePage.makeSliceConfig() hands the worker for a custom normal.

    The library's default frame (``kde._orthogonal_axis``) differs for some
    normals, so the page sends this frame to /api/kde/slice (ux..vz) and the
    parity cases pin it: both runtimes grid the same plane the same way.
    """
    n = np.asarray(normal, dtype=float)
    n = n / np.sqrt(float(np.dot(n, n)))
    reference = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.85 else np.array([0.0, 1.0, 0.0])
    u = reference - n * float(np.dot(reference, n))
    u = u / np.sqrt(float(np.dot(u, u)))
    v = np.cross(n, u)
    v = v / np.sqrt(float(np.dot(v, v)))
    return n.tolist(), u.tolist(), v.tolist()


def _real(name, path, element, orientation, z, dz, bw, *, requires_data=False, custom=None):
    if custom is not None:
        normal, u, v = browser_plane_basis(custom)
    else:
        normal, u, v = PRESETS[orientation]
    return {
        "name": name,
        "rmc6f": path,
        "element": element,
        "requiresData": requires_data,
        "normal": normal,
        "u": u,
        "v": v,
        "z": z,
        "dz": dz,
        "bw": bw,
        "grid": GRID,
    }


def _synthetic(name, points, z, dz, bw):
    normal, u, v = PRESETS["c"]
    return {
        "name": name,
        "points": np.round(np.asarray(points, dtype=float), 10).tolist(),
        "normal": normal,
        "u": u,
        "v": v,
        "z": z,
        "dz": dz,
        "bw": bw,
        "grid": GRID,
    }


def synthetic_cases() -> list[dict]:
    rng = np.random.default_rng(20260924)
    cases = []
    # One isotropic site: a full-rank but compact slab (the old browser path
    # swapped in an inflated isotropic kernel here).
    site = np.array([0.5, 0.5, 0.5]) + 0.008 * rng.standard_normal((300, 3))
    cases.append(_synthetic("single-site", site, 0.5, 0.08, 0.03))
    # Two sites on the cell anti-diagonal: a needle-shaped (but full-rank) kernel.
    a = np.array([0.25, 0.75, 0.5]) + 0.006 * rng.standard_normal((150, 3))
    b = np.array([0.75, 0.25, 0.5]) + 0.006 * rng.standard_normal((150, 3))
    cases.append(_synthetic("two-site-needle", np.vstack([a, b]), 0.5, 0.08, 0.01))
    # An exactly collinear line (dyadic coordinates, so no rounding breaks it),
    # far from the faces so no periodic image leaves the line: rank 1.
    t = np.arange(16, 49) / 64.0
    line = np.column_stack([t, 0.25 + 0.5 * t, np.full(t.size, 0.5)])
    cases.append(_synthetic("collinear", line, 0.5, 0.08, 0.03))
    # The same line through decimal coordinates: rank 2 to numpy's tolerance,
    # but collinear to within round-off (1 - rho^2 below the conditioning limit).
    decimal = np.linspace(0.2, 0.8, 40)
    rounded = np.column_stack([decimal, 0.3 + 0.7 * decimal, np.full(40, 0.5)])
    cases.append(_synthetic("collinear-to-round-off", rounded, 0.5, 0.08, 0.03))
    # Nearly collinear but genuinely two-dimensional: drawn, as a needle.
    near = np.column_stack([decimal, 0.3 + 0.7 * decimal + 1e-4 * rng.standard_normal(40), np.full(40, 0.5)])
    cases.append(_synthetic("near-collinear", near, 0.5, 0.08, 0.03))
    # Twenty atoms on two positions: fewer than 3 distinct points.
    pair = np.array([[0.3, 0.4, 0.5]] * 10 + [[0.6, 0.7, 0.5]] * 10)
    cases.append(_synthetic("two-positions", pair, 0.5, 0.08, 0.03))
    cases.append(_synthetic("three-atoms", [[0.3, 0.4, 0.5], [0.6, 0.7, 0.5], [0.4, 0.5, 0.5]], 0.5, 0.08, 0.03))
    cases.append(_synthetic("zero-bandwidth", site[:20], 0.5, 0.08, 0.0))
    # Two sites on y = 0.5, midway between node rows: the needle kernel misses
    # every node, the map is tails only (grid mass ~7e-12), flagged `unresolved`.
    between = [(x0 + dx, 0.5 + dy, 0.5) for x0 in (0.2, 0.8) for dx in (-0.01, 0.01) for dy in (-0.01, 0.01)]
    cases.append(_synthetic("unresolved-needle", between, 0.5, 0.08, 0.076))
    return cases


def real_cases() -> list[dict]:
    return [
        _real("demo Ga c bw=0.01", DEMO, "Ga", "c", 0.25, 0.08, 0.01),
        _real("demo Ga c bw=0.03", DEMO, "Ga", "c", 0.25, 0.08, 0.03),
        _real("demo Ta a bw=0.005", DEMO, "Ta", "a", 0.13, 0.08, 0.005),
        _real("demo Se c bw=0.03", DEMO, "Se", "c", 0.37, 0.02, 0.03),
        _real("demo all (111) bw=0.03", DEMO, "all", None, 0.37, 0.01, 0.03, custom=[1, 1, 1]),
        # The two custom planes whose browser frame differs from the Flask
        # route's default one; tests/test_kde_custom_frame.py checks the route
        # returns exactly these maps when the page sends its u/v.
        _real("demo Ta (110) bw=0.06", DEMO, "Ta", None, 0.37, 0.02, 0.06, custom=[1, 1, 0]),
        _real("demo Se (101) bw=0.05", DEMO, "Se", None, 0.37, 0.01, 0.05, custom=[1, 0, 1]),
        _real("5K Ga c bw=0.01", DATA_5K, "Ga", "c", 0.25, 0.08, 0.01, requires_data=True),
        _real("5K Ga c bw=0.015", DATA_5K, "Ga", "c", 0.25, 0.08, 0.015, requires_data=True),
        _real("5K Ga c bw=0.03", DATA_5K, "Ga", "c", 0.25, 0.08, 0.03, requires_data=True),
        _real("5K Nb c bw=0.005", DATA_5K, "Nb", "c", 0.15, 0.08, 0.005, requires_data=True),
        _real("5K Se (110) bw=0.03", DATA_5K, "Se", None, 0.49, 0.02, 0.03, requires_data=True, custom=[1, 1, 0]),
    ]


def case_positions(case: dict, root: Path = ROOT) -> np.ndarray | None:
    """Fractional positions for a case, or None when its data file is absent."""
    if "points" in case:
        return np.asarray(case["points"], dtype=float)
    path = root / case["rmc6f"]
    if not path.exists():
        return None
    element = None if case["element"] == "all" else case["element"]
    return load_unit_cell_positions(path, element=element).fractional_positions


def compute_case(case: dict, positions: np.ndarray) -> dict:
    result = oriented_kde_slice(
        positions,
        center=case["z"],
        thickness=case["dz"],
        normal=np.asarray(case["normal"], dtype=float),
        u_axis=np.asarray(case["u"], dtype=float),
        v_axis=np.asarray(case["v"], dtype=float),
        bw=case["bw"],
        grid=case["grid"],
        log=False,
        n_levels=0,
    )
    return {
        "slabCount": result["slabCount"],
        "fitCount": result["fitCount"],
        "message": result["message"],
        "messageCode": result["messageCode"],
        "warnings": [warning["code"] for warning in result["warnings"]],
        "kernel": result["kernel"],
        "vmax": result["vmax"],
        # The grid relative to its peak, rounded to 1e-8: far below the 1e-6
        # parity tolerance, and it keeps the golden small.
        "densityOverPeak": [
            [round(value / result["vmax"], 8) if result["vmax"] > 0 else 0.0 for value in row]
            for row in result["density"]
        ],
    }


def main() -> None:
    cases = []
    for case in synthetic_cases() + real_cases():
        positions = case_positions(case)
        if positions is None:
            raise SystemExit(f"{case['rmc6f']} is missing: run the generator where data/ is present")
        expected = compute_case(case, positions)
        if expected["fitCount"] >= MAX_KDE_FIT_POINTS:
            raise SystemExit(f"{case['name']}: slab reaches the fit cap; the runtimes would subsample differently")
        cases.append({**case, "expected": expected})
    payload = {"note": "Generated by tests/generate_kde_fixture.py -- do not edit by hand.", "cases": cases}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload) + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
