# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the JS parity golden for the PCA-KDE fit subsample and volume.

Writes ``web_app/frontend/src/__tests__/fixtures/pca_kde_fixture.json``, which
``src/workers/__tests__/pcaKdeParity.test.js`` holds ``workers/pcaKde.js`` to and
``tests/test_pca_kde_fixture.py`` checks against the Python engine. Above
``MAX_PCA_FIT_POINTS`` (20 000) both engines fit the KDE on the rows
``fit_subsample_indices`` / ``fitSubsampleIndices`` draw -- mulberry32 and a
partial Fisher-Yates -- so the golden pins the draw index for index and then the
volume that follows from it:

* ``stream``: the first outputs of mulberry32 for several seeds, as uint32
  (negative and > 2**32 seeds pin the ``seed >>> 0`` derivation);
* ``draws``: drawn indices for (count, limit, seed) cases -- in full for the
  small ones; for the 28^3-copy site (21 952 of 20 000) as the complement, which
  with the ascending order fixes every index; for larger clouds as checksums;
* ``volumes``: ``pca_kde_volume`` / ``site_pca_kde`` above the cap: a 21 952-copy
  synthetic site (a 28^3 box, the smallest cubic box above the cap) whose
  points both runtimes build from the same mulberry32 stream with exactly
  rounded arithmetic only, so no point cloud is stored; and the committed demo
  run's pooled Se cloud (32 000 atoms written cell by cell, so it also pins the
  site-major order of a pooled cloud).

Regenerate from the checkout under test whenever the draw or the KDE changes
(without ``PYTHONPATH=.`` another checkout's ``rmc_toolkits`` may be imported):

    PYTHONPATH=. python tests/generate_pca_kde_fixture.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from rmc_toolkits.pca_kde import (
    MAX_PCA_FIT_POINTS,
    displacement_cloud,
    fit_subsample_indices,
    load_site_displacements,
    mulberry32_uint32,
    pca_kde_volume,
    site_pca_kde,
)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "pca_kde_fixture.json"
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"
DEMO_FILE = "GTS_250K.rmc6f"

STREAM_SEEDS = [0, 1, 2026, 4294967295, -1, 2**32 + 5]
STREAM_LENGTH = 8

# (count, limit, seed, how the indices are stored)
DRAWS = [
    (50, 50, 0, "indices"),  # at the cap: nothing drawn
    (10, 4, 0, "indices"),
    (100, 37, 7, "indices"),
    (1000, 999, 3, "indices"),
    (21952, MAX_PCA_FIT_POINTS, 0, "excluded"),  # one site of a 28^3 box
    (27000, MAX_PCA_FIT_POINTS, 1, "checksum"),
    (32000, MAX_PCA_FIT_POINTS, 0, "checksum"),  # the demo's pooled Se cloud
]

# The synthetic site: point i is M^T e_i with e_ik = (u + u' + u'') - 1.5 from
# nine consecutive mulberry32 outputs (Irwin-Hall, excess kurtosis -0.4), in Å.
SYNTHETIC = {
    "seed": 2028,
    "count": 21952,
    "transform": [[0.16, 0.03, -0.02], [0.0, 0.11, 0.04], [0.05, 0.0, 0.08]],
}

VOLUMES = [
    {
        "name": "synthetic site, 28^3 = 21 952 copies",
        "source": {"synthetic": True},
        "options": {"grid": 12, "extent": 3.0, "projections": True},
    },
    {
        "name": "demo run, Se pooled (32 000 atoms)",
        "source": {"demo": DEMO_FILE, "element": "Se"},
        "options": {"grid": 10, "extent": 4.0, "cubicBox": True, "probability": 0.5, "projections": True},
    },
]

PYTHON_NAMES = {"cubicBox": "cubic_box", "bwScale": "bw_scale"}


def synthetic_cloud(spec: dict = SYNTHETIC) -> np.ndarray:
    """The synthetic site's points. Only +, - and * (no fused, reordered or
    library operation), in the order pcaKdeParity.test.js uses, so the two
    runtimes build the same doubles."""
    count = int(spec["count"])
    u = mulberry32_uint32(spec["seed"], 9 * count).astype(float) / 4294967296.0
    u = u.reshape(count, 3, 3)
    e = u[:, :, 0] + u[:, :, 1] + u[:, :, 2] - 1.5
    m = spec["transform"]
    return np.column_stack(
        [e[:, 0] * m[0][j] + e[:, 1] * m[1][j] + e[:, 2] * m[2][j] for j in range(3)]
    )


def index_checksum(indices: np.ndarray) -> int:
    """sum((i + 1) * index_i): exact in a double for these sizes, and order-sensitive."""
    indices = np.asarray(indices, dtype=np.int64)
    return int(((np.arange(indices.size, dtype=np.int64) + 1) * indices).sum())


def cloud_checksum(cloud: np.ndarray) -> float:
    """sum(((i mod 1000) + 1) * (x_i + 2 y_i + 3 z_i)): which row sits where."""
    weights = (np.arange(cloud.shape[0]) % 1000 + 1).astype(float)
    return float(weights @ (cloud[:, 0] + 2.0 * cloud[:, 1] + 3.0 * cloud[:, 2]))


def _draw(count: int, limit: int, seed: int, stored: str) -> dict:
    indices = fit_subsample_indices(count, limit, seed)
    entry = {"count": count, "limit": limit, "seed": seed, "size": int(indices.size),
             "checksum": index_checksum(indices)}
    if stored == "indices":
        entry["indices"] = indices.tolist()
    elif stored == "excluded":
        entry["excluded"] = np.setdiff1d(np.arange(count), indices).tolist()
    return entry


def _python_options(options: dict) -> dict:
    return {PYTHON_NAMES.get(name, name): value for name, value in options.items()}


def _volume(case: dict) -> dict:
    source = case["source"]
    kwargs = _python_options(case["options"])
    if source.get("synthetic"):
        cloud = synthetic_cloud()
        result = pca_kde_volume(cloud, **kwargs)
    else:
        sites = load_site_displacements(DEMO / source["demo"])
        selection = {"element": source.get("element"), "reference_number": source.get("referenceNumber")}
        cloud = displacement_cloud(sites, **selection)
        result = site_pca_kde(sites, **selection, **kwargs)
    if cloud.shape[0] <= MAX_PCA_FIT_POINTS:
        raise SystemExit(f"{case['name']}: {cloud.shape[0]} points do not exceed the fit cap")
    return {
        "name": case["name"],
        "source": source,
        "options": case["options"],
        "cloud": {
            "count": int(cloud.shape[0]),
            "head": cloud[:3].tolist(),
            "tail": cloud[-3:].tolist(),
            "checksum": cloud_checksum(cloud),
        },
        "fitChecksum": index_checksum(fit_subsample_indices(cloud.shape[0], MAX_PCA_FIT_POINTS, 0)),
        "expected": result,
    }


def build_fixture() -> dict:
    """The whole fixture as a JSON-ready dict (``main`` writes it)."""
    return {
        "note": "Generated by tests/generate_pca_kde_fixture.py -- do not edit by hand.",
        "maxFitPoints": MAX_PCA_FIT_POINTS,
        "stream": [
            {"seed": seed, "uint32": mulberry32_uint32(seed, STREAM_LENGTH).tolist()} for seed in STREAM_SEEDS
        ],
        "draws": [_draw(*draw) for draw in DRAWS],
        "synthetic": SYNTHETIC,
        "volumes": [_volume(case) for case in VOLUMES],
    }


def main() -> None:
    text = json.dumps(build_fixture(), separators=(",", ":"), allow_nan=False)
    OUT.write_text(text + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    main()
