# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the JS parity golden for the cluster-distortion engine.

Writes ``web_app/frontend/src/__tests__/fixtures/clusters_fixture.json``:
synthetic ``.rmc6f`` configurations (tests/cluster_configs.py) and the
reference results of ``rmc_toolkits.clusters.cluster_distortion``. The vitest
suite parses the same text with ``siteDisplacementsFromRmc6f``, runs
``workers/clusters.js`` and asserts exact integers (templates, counts, axis
fractions, the shuffles' hash) and floats to 1e-9. Regenerate whenever the
engine's conventions change:

    python tests/generate_clusters_fixture.py
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cluster_configs import (  # noqa: E402
    chains_config_text,
    cluster_config_text,
    origin_dimer_config_text,
    tie_config_text,
    with_vacancies,
)
from rmc_toolkits.clusters import SHUFFLE_SEED, cluster_distortion, mix32, shuffle_permutation  # noqa: E402
from rmc_toolkits.pca_kde import load_site_displacements  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "clusters_fixture.json"

CASES = [
    ("trigonal", lambda: cluster_config_text("trigonal", supercell=(4, 4, 4)),
     [("Nb", [2.5, 3.5]), ("Ga", [2.0, 3.0]), ("Se", [7.5, 8.0]), ("Nb", [0.5, 1.0])]),
    ("static", lambda: cluster_config_text("static", supercell=(4, 4, 4)), [("Nb", [2.5, 3.5])]),
    # A site a hair below fraction 0 (review 2026-10-12 #1).
    ("origin-dimer", lambda: origin_dimer_config_text(-1e-4), [("V", [2.0, 2.5])]),
    # Two chain templates with differently listed edges (review 2026-10-12 #2).
    ("chains", chains_config_text, [("Mo", [1.8, 2.2])]),
    # Clusters missing an atom (review 2026-10-12 #16).
    ("vacancies", lambda: with_vacancies(cluster_config_text("independent", supercell=(4, 4, 4)), "Nb", 0.1), [("Nb", [2.5, 3.5])]),
    # As many triangles as dimers (review 2026-10-12 #17).
    ("tie", tie_config_text, [("Mo", [1.8, 2.2])]),
]


def main() -> None:
    cases = []
    with tempfile.TemporaryDirectory() as directory:
        for mode, make_text, queries in CASES:
            text = make_text()
            path = Path(directory) / f"{mode}.rmc6f"
            path.write_text(text, encoding="utf-8")
            sites = load_site_displacements(path)
            expected = [cluster_distortion(sites, element, tuple(window)) for element, window in queries]
            cases.append({
                "name": mode,
                "text": text,
                "queries": [{"element": element, "window": window} for element, window in queries],
                "expected": expected,
            })
    hashes = {
        "mix32": {str(value): int(mix32(np.uint32(value))) for value in (0, 1, 12345, 0xFFFFFFFF)},
        "permutation": {"n": 20, "shuffle": 3, "template": 1, "member": 2, "seed": SHUFFLE_SEED,
                        "order": shuffle_permutation(20, 3, 1, 2).tolist()},
    }
    OUT.write_text(json.dumps({"cases": cases, "hashes": hashes}, indent=None, separators=(",", ":")), encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    main()
