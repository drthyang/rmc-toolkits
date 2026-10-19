# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the JS parity golden for the cluster-distortion engine.

Writes ``web_app/frontend/src/__tests__/fixtures/clusters_fixture.json``:
synthetic ``.rmc6f`` configurations (tests/cluster_configs.py) and the
reference results of ``rmc_toolkits.clusters.cluster_distortion``, verdicts
included, with the verdicts' closed vocabularies (``enums``). The vitest suite
parses the same text with ``siteDisplacementsFromRmc6f``, runs
``workers/clusters.js`` and asserts exact integers, strings and booleans
(templates, counts, axis fractions, the shuffles' hash, every answer and
``because``) and floats to 1e-9. ``tests/test_clusters_fixture.py`` fails
when the committed file no longer matches the Python engine.

Every statistic a verdict decides on sits more than ``MARGIN`` from its
threshold in every case (:func:`assert_margins`; the generator refuses to
write otherwise), so the two engines -- which agree to ~1e-12 -- cannot
disagree on a verdict here. Regenerate from the checkout under test whenever
the engine's conventions change:

    PYTHONPATH=. python tests/generate_clusters_fixture.py
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
    octahedra_config_text,
    origin_dimer_config_text,
    rigid_dimer_config_text,
    tie_config_text,
    with_vacancies,
)
from rmc_toolkits.clusters import (  # noqa: E402
    CHANGE_KINDS,
    CLUSTER_BECAUSE,
    CLUSTER_RULES,
    EVIDENCE_AGAINST,
    SHUFFLE_SEED,
    cluster_distortion,
    mix32,
    shuffle_permutation,
)
from rmc_toolkits.pca_kde import load_site_displacements  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "clusters_fixture.json"

#: The least distance (in its own units) between a decided statistic and its threshold.
MARGIN = 1e-6

CASES = [
    ("trigonal", lambda: cluster_config_text("trigonal", supercell=(4, 4, 4)),
     [("Nb", [2.5, 3.5]), ("Ga", [2.0, 3.0]), ("Se", [7.5, 8.0]), ("Nb", [0.5, 1.0])]),
    ("static", lambda: cluster_config_text("static", supercell=(4, 4, 4)), [("Nb", [2.5, 3.5])]),
    # Collective breathing: size, not shape (review 2026-10-14 L2).
    ("breathing", lambda: cluster_config_text("breathing", supercell=(4, 4, 4)), [("Nb", [2.5, 3.5])]),
    # A site a hair below fraction 0 (review 2026-10-12 #1).
    ("origin-dimer", lambda: origin_dimer_config_text(-1e-4), [("V", [2.0, 2.5])]),
    # Two chain templates with differently listed edges (review 2026-10-12 #2).
    ("chains", chains_config_text, [("Mo", [1.8, 2.2])]),
    # Clusters missing an atom (review 2026-10-12 #16).
    ("vacancies", lambda: with_vacancies(cluster_config_text("independent", supercell=(4, 4, 4)), "Nb", 0.1), [("Nb", [2.5, 3.5])]),
    # As many triangles as dimers (review 2026-10-12 #17): mixed sizes.
    ("tie", tie_config_text, [("Mo", [1.8, 2.2])]),
    # The verdicts' known answers (2026-10-19): one edge pinched (a shape
    # distortion beyond independent atoms, one edge against its opposite),
    # rigid tetrahedra (less distorted, rigid), Jahn-Teller F6 octahedra with
    # every F site split (split member sites, the pattern), rigid dimers.
    ("pinched", lambda: cluster_config_text("pinched", supercell=(4, 4, 4)), [("Nb", [2.5, 3.5])]),
    ("rigid", lambda: cluster_config_text("rigid", supercell=(4, 4, 4)), [("Nb", [2.5, 3.5])]),
    ("split-octahedra", octahedra_config_text, [("F", [2.5, 3.25])]),
    ("rigid-dimer", rigid_dimer_config_text, [("V", [2.0, 2.5])]),
]


def _records(verdict: dict):
    """Every Evidence record of a verdict, with the name of its finding."""
    for name, finding in verdict.items():
        if name == "rules":
            continue
        items = finding.items() if name == "window" else [(name, finding)]
        for key, item in items:
            for record in item["evidence"]:
                yield key, record


def assert_margins(verdict: dict, where: str) -> None:
    """Every statistic a verdict compares with a threshold sits more than MARGIN from it."""
    for name, record in _records(verdict):
        side = record.get("side")
        for field, bar in record["bar"].items():
            value = record[field]
            if value is None:
                continue
            thresholds = [bar] if side == "above" else [-bar] if side == "below" else [bar, -bar]
            for threshold in thresholds:
                if abs(value - threshold) <= MARGIN:
                    raise AssertionError(
                        f"{where}: {name} {record['statistic']} {field} = {value!r} lies within {MARGIN} of {threshold}"
                    )


def build_fixture(texts: dict | None = None) -> dict:
    """The whole fixture as a JSON-ready dict (``main`` writes it).

    ``texts`` (case name -> ``.rmc6f`` text) replaces the generated
    configurations: the staleness test passes the committed ones, so it checks
    the engine, not NumPy's random stream.
    """
    cases = []
    with tempfile.TemporaryDirectory() as directory:
        for name, make_text, queries in CASES:
            text = texts[name] if texts is not None else make_text()
            path = Path(directory) / f"{name}.rmc6f"
            path.write_text(text, encoding="utf-8")
            sites = load_site_displacements(path)
            expected = []
            for element, window in queries:
                result = cluster_distortion(sites, element, tuple(window))
                assert_margins(result["verdict"], f"{name} {element} {window}")
                expected.append(result)
            cases.append({
                "name": name,
                "text": text,
                "queries": [{"element": element, "window": window} for element, window in queries],
                "expected": expected,
            })
    hashes = {
        "mix32": {str(value): int(mix32(np.uint32(value))) for value in (0, 1, 12345, 0xFFFFFFFF)},
        "permutation": {"n": 20, "shuffle": 3, "template": 1, "member": 2, "seed": SHUFFLE_SEED,
                        "order": shuffle_permutation(20, 3, 1, 2).tolist()},
    }
    enums = {
        "because": {name: list(codes) for name, codes in CLUSTER_BECAUSE.items()},
        "against": list(EVIDENCE_AGAINST),
        "changeKinds": list(CHANGE_KINDS),
        "rules": dict(CLUSTER_RULES),
    }
    return {"cases": cases, "hashes": hashes, "enums": enums}


def main() -> None:
    text = json.dumps(build_fixture(), indent=None, separators=(",", ":"), allow_nan=False)
    OUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    main()
