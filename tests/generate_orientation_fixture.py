# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the JS parity golden for the displacement-direction engine.

Writes ``web_app/frontend/src/__tests__/fixtures/orientation_fixture.json``:
deterministic inputs and what ``rmc_toolkits.orientation`` makes of them --
the one written interface that ``workers/orientation.js`` is held to
(``src/__tests__/orientationParity.test.js``), in place of values copied by
hand into both suites. ``tests/test_orientation_fixture.py`` fails when the
committed file no longer matches the Python engine. Contents:

* ``histograms`` -- full ``orientation_histogram`` / ``site_orientation_histogram``
  payloads: synthetic clouds (a preferred direction, isotropic, antipodally
  asymmetric, sparse, singly occupied pairs, 8 wells) stored in ``inputs``, and
  sites of the committed demo run (``web_app/frontend/public/demo``), at
  several resolutions, smoothing passes, weights and frames. ``keys`` is each
  payload's whole key set; the per-cell geometry (``centers``, ``areas``,
  ``sizes``, ``antipode``) depends only on the frequency and sits once per
  frequency in ``tilings``.
* ``assignments`` -- ``assign_cells`` of the exact tie families (the axes, the
  26 <100>/<110>/<111> directions).
* ``tilingChecks`` -- ``neighbors`` checksums (the sum over cells c and slots k
  of (k + 1)(neighbors[c][k] + 1)(c % 97 + 1)) and the first polygon vertex of
  the two nu = 38 cells whose start once differed between the engines.
* ``recommendedFrequency`` -- the boundary table, where 12 (10 nu^2 + 2) == N.
* ``pearsonMoments`` -- the exact null mean, variance and third moment of
  Pearson's X^2 for unequal cells.

Regenerate from the checkout under test whenever the engine's conventions
change (without ``PYTHONPATH=.`` another checkout's ``rmc_toolkits`` may be
imported):

    PYTHONPATH=. python tests/generate_orientation_fixture.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from rmc_toolkits.orientation import (
    _pearson_null_moments,
    assign_cells,
    goldberg_tiling,
    orientation_histogram,
    recommended_frequency,
    site_orientation_histogram,
)
from rmc_toolkits.pca_kde import displacement_cloud, load_site_displacements

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "orientation_fixture.json"
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"
DEMO_FILE = "GTS_250K.rmc6f"

# Per-cell geometry that depends on the frequency alone (stored in `tilings`).
TILING_KEYS = ("centers", "areas", "sizes", "antipode")
# Per-cell values stored only for the cases marked `cells` (the integer
# `counts`, which pin where every point went, are stored for every case).
PER_CELL_KEYS = ("mass", "density", "enhancement", "expected", "zScore", "cellMeanAmplitude")

# The JS options and their Python keyword names.
OPTION_NAMES = {
    "frequency": "frequency",
    "weight": "weight",
    "minAmplitude": "min_amplitude",
    "minAmplitudeQuantile": "min_amplitude_quantile",
    "smoothing": "smoothing",
    "frame": "frame",
    "geometry": "geometry",
    "targetPerCell": "target_per_cell",
}


# --- inputs --------------------------------------------------------------------


def golden_sphere(n_sphere: int = 900) -> np.ndarray:
    """A Fibonacci-spiral sphere, radius modulated by index (RNG-free)."""
    i = np.arange(n_sphere, dtype=float)
    z = 1.0 - (2.0 * i + 1.0) / n_sphere
    r = np.sqrt(1.0 - z * z)
    phi = i * 2.399963229728653
    radius = 0.05 + 0.03 * ((i * 7) % 11) / 11.0
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), z]) * radius[:, None]


def golden_lobe(n_lobe: int, lobe=(0.3, -0.5, 0.81)) -> np.ndarray:
    """A tight lobe of ``n_lobe`` points around ``lobe`` (RNG-free)."""
    j = np.arange(n_lobe, dtype=float)
    unit = np.asarray(lobe) / np.linalg.norm(lobe)
    jitter = 0.02 * np.column_stack([np.cos(j * 1.3), np.sin(j * 1.7), np.cos(j * 0.9)])
    return (unit[None, :] + jitter) * 0.12


def golden_cloud(n_sphere: int = 900, n_lobe: int = 60, lobe=(0.3, -0.5, 0.81)) -> np.ndarray:
    """The golden sphere plus ``n_lobe`` lobe points: a cross-engine input."""
    return np.vstack([golden_sphere(n_sphere), golden_lobe(n_lobe, lobe)])


def special_directions() -> np.ndarray:
    """<111>, <100> and <110> (26 directions), in a fixed order."""
    signs = (1.0, -1.0)
    directions = [[a, b, c] for a in signs for b in signs for c in signs]
    for axis in range(3):
        for sign in signs:
            v = [0.0, 0.0, 0.0]
            v[axis] = sign
            directions.append(v)
    for i, j in ((0, 1), (0, 2), (1, 2)):
        for a in signs:
            for b in signs:
                v = [0.0, 0.0, 0.0]
                v[i], v[j] = a, b
                directions.append(v)
    return np.asarray(directions)


def well_model(length: float = 0.1) -> np.ndarray:
    """Atoms on the eight <111> wells: 30 at +[111], 30 at -[111], 20 at each other."""
    rows = []
    for direction in special_directions()[:8]:
        weight = 30 if abs(direction.sum()) == 3 else 20
        rows.append(np.repeat(direction[None, :] * length / np.sqrt(3.0), weight, axis=0))
    return np.vstack(rows)


def neighbor_checksum(neighbors) -> int:
    """Sum over cells c and slots k of (k + 1) * (neighbors[c][k] + 1) * (c % 97 + 1)."""
    neighbors = np.asarray(neighbors)
    cells = np.arange(neighbors.shape[0])
    weights = (np.arange(neighbors.shape[1]) + 1)[None, :] * ((cells % 97) + 1)[:, None]
    return int((weights * (neighbors + 1)).sum())


def pinned_probabilities(cells: int) -> np.ndarray:
    """Deterministic unequal cell probabilities."""
    raw = 1.0 + 0.3 * np.sin(np.arange(1, cells + 1, dtype=float))
    return raw / raw.sum()


def _rounded(vectors: np.ndarray) -> np.ndarray:
    # Random inputs are stored to 1e-6 A; both engines read the rounded values.
    return np.round(vectors, 6)


def _isotropic(n: int, seed: int) -> np.ndarray:
    return _rounded(np.random.default_rng(seed).normal(size=(n, 3)) * 0.1)


def _hemisphere(n: int, seed: int) -> np.ndarray:
    vectors = np.random.default_rng(seed).normal(size=(n, 3)) * 0.1
    vectors[:, 0] = np.abs(vectors[:, 0])
    return _rounded(vectors)


def _prolate(n: int, seed: int) -> np.ndarray:
    return _rounded(np.random.default_rng(seed).normal(size=(n, 3)) * np.array([0.125, 0.1, 0.1]))


def _singly_occupied_pairs() -> np.ndarray:
    # One atom in each of 20 cells whose antipodes are empty: every pair holds
    # one atom, so the asymmetry equals its null exactly and has no z.
    tiling = goldberg_tiling(4)
    cells = np.flatnonzero(np.arange(tiling.cell_count) < tiling.antipode)[:20]
    return tiling.centers[cells] * 0.1


def inputs() -> dict[str, np.ndarray]:
    """The named synthetic clouds (Cartesian Angstrom), as stored in the fixture."""
    return {
        "golden-sphere": golden_sphere(),
        "golden-lobe": golden_lobe(400),
        "well-model": well_model(),
        "isotropic": _isotropic(600, 101),
        "hemisphere": _hemisphere(600, 102),
        "prolate": _prolate(600, 103),
        "sparse": _isotropic(8, 104),
        "singly-occupied-pairs": _singly_occupied_pairs(),
    }


def golden(n_lobe: int) -> dict:
    """Input spec of golden_cloud(n_lobe=n_lobe): the sphere, then the lobe's first rows."""
    parts = [["golden-sphere", 900]]
    if n_lobe:
        parts.append(["golden-lobe", n_lobe])
    return {"vectors": parts}


def named(name: str) -> dict:
    return {"vectors": [[name, None]]}


def demo_site(reference_number: int) -> dict:
    return {"demo": DEMO_FILE, "referenceNumber": reference_number}


def demo_element(element: str) -> dict:
    return {"demo": DEMO_FILE, "element": element}


# (name, input spec, JS options, cells). Every case sets geometry false unless
# it names it: the polygons are pinned once, at nu = 3. `cells` stores the
# per-cell values (PER_CELL_KEYS, and the frequency's tiling) beside the
# integer counts: once for each weight, frame, smoothing and amplitude cut.
HISTOGRAMS = [
    # The golden cloud: a lobe toward (0.3, -0.5, 0.81) on a near-uniform sphere.
    ("golden lobe, nu 6, 1 pass", golden(60), {"frequency": 6, "smoothing": 1}, True),
    ("golden lobe, nu 10, 2 passes", golden(60), {"frequency": 10, "smoothing": 2}, False),
    ("golden lobe, recommended nu", golden(60), {"frequency": None, "smoothing": 0}, True),
    ("golden lobe, amplitude2 weight", golden(60), {"frequency": 6, "smoothing": 1, "weight": "amplitude2"}, False),
    # A Fibonacci sphere is more uniform than counting noise.
    ("golden sphere, nu 10, 2 passes", golden(0), {"frequency": 10, "smoothing": 2}, False),
    ("golden sphere with 6 lobe points, nu 2", golden(6), {"frequency": 2, "smoothing": 0}, False),
    # A one-sided lobe: antipodally asymmetric.
    ("golden sphere with 400 lobe points, nu 4", golden(400), {"frequency": 4, "smoothing": 0}, False),
    ("golden lobe, nu 3, 2 passes, geometry", golden(60), {"frequency": 3, "smoothing": 2, "geometry": True}, True),
    # Exact ties: the eight <111> wells sit on Voronoi edges and vertices.
    ("8-well <111> model, nu 2", named("well-model"), {"frequency": 2}, True),
    ("8-well <111> model, nu 4", named("well-model"), {"frequency": 4}, False),
    ("isotropic Gaussian, recommended nu", named("isotropic"), {}, False),
    ("isotropic Gaussian, nu 6, 1 pass", named("isotropic"), {"frequency": 6, "smoothing": 1}, False),
    ("hemisphere (one-sided), nu 10, 2 passes", named("hemisphere"), {"frequency": 10, "smoothing": 2}, False),
    (
        "prolate along x, PCA frame, amplitude weight, 10% cut",
        named("prolate"),
        {"frequency": 4, "frame": "pca", "weight": "amplitude", "minAmplitudeQuantile": 0.1},
        True,
    ),
    ("prolate along x, nu 10, 2 passes", named("prolate"), {"frequency": 10, "smoothing": 2}, False),
    # Fewer than MAP_TEST_MIN_PAIRS coincident pairs: no map p-value.
    ("8 atoms at nu 10 (map test withheld)", named("sparse"), {"frequency": 10}, False),
    ("singly occupied pairs (no asymmetry z)", named("singly-occupied-pairs"), {"frequency": 4}, False),
    # The committed demo run (GaTa4Se8, 250 K; 1000 copies per site).
    ("demo Ta site 21, recommended nu", demo_site(21), {}, False),
    ("demo Ta site 21, nu 6, 1 pass", demo_site(21), {"frequency": 6, "smoothing": 1}, False),
    (
        "demo Se site 18, nu 4, PCA frame, 2 passes",
        demo_site(18),
        {"frequency": 4, "frame": "pca", "smoothing": 2},
        True,
    ),
    ("demo Ga pooled, nu 10, 2 passes", demo_element("Ga"), {"frequency": 10, "smoothing": 2}, False),
]

ASSIGNMENTS = [
    ("+x, +y, +z, -x, -y, -z", np.vstack([np.eye(3), -np.eye(3)]), (3,)),
    ("8 x <111>, 6 x <100>, 12 x <110>", special_directions(), (1, 2, 3, 4, 10, 11)),
]

NEIGHBOR_CHECKSUM_FREQUENCIES = (4, 6, 10, 38)
POLYGON_STARTS = {38: (2052, 4275)}

# (N, target per cell, max frequency): the boundaries sit where 12 (10 nu^2 + 2) == N.
RECOMMENDED_FREQUENCY = [
    (0, 12, 24), (294, 12, 24), (300, 12, 24), (503, 12, 24), (504, 12, 24), (774, 12, 24),
    (1000, 12, 24), (1103, 12, 24), (1104, 12, 24), (12000, 12, 24), (12023, 12, 24),
    (12024, 12, 24), (10_000_000, 12, 24), (1000, 5, 24), (5000, 12, 3),
]

# (N, C) for pinned_probabilities(C).
PEARSON_MOMENTS = [(3, 4), (6, 5), (9, 4), (12, 3), (2, 7)]


# --- the fixture -----------------------------------------------------------------


def _vectors(spec: dict, named_inputs: dict[str, np.ndarray]) -> np.ndarray:
    return np.vstack([
        named_inputs[name] if count is None else named_inputs[name][:count]
        for name, count in spec["vectors"]
    ])


def _python_options(options: dict) -> dict:
    merged = {"geometry": False, **options}
    return {OPTION_NAMES[key]: value for key, value in merged.items()}


def _pca_spectrum(vectors: np.ndarray) -> list[float]:
    """Descending eigenvalues behind ``pcaAxes``: the covariance of every vector
    (before any amplitude cut), or the identity below four vectors."""
    if vectors.shape[0] < 4:
        return [1.0, 1.0, 1.0]
    centered = vectors - vectors.mean(axis=0)
    eigenvalues = np.linalg.eigvalsh(np.cov(centered, rowvar=False, bias=False))[::-1]
    return np.maximum(eigenvalues, 0.0).tolist()


def build_fixture(named_inputs: dict[str, np.ndarray] | None = None) -> dict:
    """The whole fixture as a JSON-ready dict (``main`` writes it).

    ``named_inputs`` replaces the synthetic clouds of :func:`inputs` -- the
    staleness test passes the committed ones, so it checks the engine, not
    NumPy's random stream.
    """
    if named_inputs is None:
        named_inputs = inputs()
    named_inputs = {name: np.asarray(vectors, dtype=float) for name, vectors in named_inputs.items()}
    demo_sites = load_site_displacements(DEMO / DEMO_FILE)
    tilings: dict[str, dict] = {}
    histograms = []
    for name, spec, options, cells in HISTOGRAMS:
        kwargs = _python_options(options)
        if "demo" in spec:
            selection = {"reference_number": spec.get("referenceNumber"), "element": spec.get("element")}
            vectors = displacement_cloud(demo_sites, **selection)
            result = site_orientation_histogram(demo_sites, **selection, **kwargs)
        else:
            vectors = _vectors(spec, named_inputs)
            result = orientation_histogram(vectors, **kwargs)
        keys = sorted(result)
        geometry = {key: result.pop(key) for key in TILING_KEYS}
        if cells:
            tilings.setdefault(str(result["frequency"]), geometry)
        else:
            for key in PER_CELL_KEYS:
                del result[key]
        histograms.append({
            "name": name,
            "input": spec,
            "options": options,
            "cells": cells,
            "keys": keys,
            # The eigenvalues each axis set belongs to, so a comparison can
            # tell a fixed axis from one inside a degenerate eigenspace.
            "spectra": {
                "pcaAxes": _pca_spectrum(vectors),
                "orientationAxes": result["orientationEigenvalues"],
            },
            "expected": result,
        })

    assignments = [
        {
            "name": name,
            "directions": directions.tolist(),
            "cells": {str(nu): assign_cells(goldberg_tiling(nu), directions).tolist() for nu in frequencies},
        }
        for name, directions, frequencies in ASSIGNMENTS
    ]
    tiling_checks = {
        "neighborChecksums": {
            str(nu): neighbor_checksum(goldberg_tiling(nu).neighbors) for nu in NEIGHBOR_CHECKSUM_FREQUENCIES
        },
        "polygonStarts": {
            str(nu): {str(cell): goldberg_tiling(nu).polygons[cell, 0].tolist() for cell in cells}
            for nu, cells in POLYGON_STARTS.items()
        },
    }
    frequencies = [
        {
            "n": n,
            "targetPerCell": target,
            "maxFrequency": cap,
            "frequency": recommended_frequency(n, target_per_cell=target, max_frequency=cap),
        }
        for n, target, cap in RECOMMENDED_FREQUENCY
    ]
    moments = [
        {
            "used": used,
            "probabilities": pinned_probabilities(cells).tolist(),
            "moments": [float(value) for value in _pearson_null_moments(used, pinned_probabilities(cells))],
        }
        for used, cells in PEARSON_MOMENTS
    ]
    return {
        "inputs": {name: vectors.tolist() for name, vectors in named_inputs.items()},
        "histograms": histograms,
        "tilings": tilings,
        "assignments": assignments,
        "tilingChecks": tiling_checks,
        "recommendedFrequency": frequencies,
        "pearsonMoments": moments,
    }


def main() -> None:
    text = json.dumps(build_fixture(), separators=(",", ":"), allow_nan=False)
    OUT.write_text(text + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    main()
