# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the JS parity golden for the displacement-modulation engine.

Writes ``web_app/frontend/src/__tests__/fixtures/modulations_fixture.json``:
synthetic ``.rmc6f`` configurations (tests/cluster_configs.py) and the
reference results of ``rmc_toolkits.modulations.displacement_modulations``,
verdicts included, with the verdicts' closed vocabularies (``enums``). The
vitest suite runs ``workers/modulations.js`` on the same text and asserts
integers, strings and booleans exactly and floats to 1e-9.
``tests/test_modulations_fixture.py`` fails when the committed file no longer
matches the Python engine.

Every ratio and significance a verdict decides on sits more than ``MARGIN``
from its threshold in every case (:func:`assert_margins`; the generator
refuses to write otherwise), so the two engines cannot disagree on a verdict
here. Regenerate from the checkout under test whenever the engine's
conventions change:

    PYTHONPATH=. python tests/generate_modulations_fixture.py
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cluster_configs import cluster_config_text, long_wave_config_text, wave_config_text  # noqa: E402
from rmc_toolkits.modulations import (  # noqa: E402
    LONG_WAVE_RATIO,
    MODE_KINDS,
    MODULATION_BECAUSE,
    MODULATION_RULES,
    displacement_modulations,
)
from rmc_toolkits.pca_kde import load_site_displacements  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "modulations_fixture.json"

#: The least distance between a decided ratio or sigma and its threshold.
MARGIN = 1e-6

CASES = [
    ("doubled (4, 4, 4)", lambda: cluster_config_text("doubled", supercell=(4, 4, 4)), [None, "Nb"]),
    ("independent (3, 4, 4)", lambda: cluster_config_text("independent", supercell=(3, 4, 4)), [None]),
    # Order under a stronger long wave (review 2026-10-12 #6).
    ("long-wave (6, 6, 6)", lambda: long_wave_config_text(supercell=(6, 6, 6)), ["Ti"]),
    # The verdicts' known answers (2026-10-19): order at a general wavevector
    # (3/8, 0, 0), and the same order under a long wave (1/8, 0, 0) that is the
    # strongest mode -- each element alone, and both.
    ("general-q (8, 4, 4)", wave_config_text, [None]),
    ("long wave and general order (8, 4, 4)", lambda: wave_config_text(long_amp=0.15, order_amp=0.12), [None, "Ti"]),
]


def assert_margins(result: dict, where: str) -> None:
    """Every ratio and sigma a verdict compares with a threshold sits more than MARGIN from it."""
    rules = result["rules"]
    for entry in result["elements"]:
        if not entry["sites"]:
            continue
        checks = []
        for point in entry["zoneBoundary"] + entry["strongest"]:
            if point["bar"] is None:
                continue
            checks += [(point["q"], "ratio", point["ratio"], point["bar"]["ratio"]),
                       (point["q"], "sigma", point["sigma"], point["bar"]["sigma"])]
        first = entry["strongest"][0]
        if first["kind"] == "long_wave":
            checks.append((first["q"], "long-wave ratio", first["ratio"], rules["longWaveRatio"]))
        for q, what, value, threshold in checks:
            if value is not None and abs(value - threshold) <= MARGIN:
                raise AssertionError(f"{where}: {entry['element']} {q} {what} = {value!r} lies within {MARGIN} of {threshold}")


def build_fixture(texts: dict | None = None) -> dict:
    """The whole fixture as a JSON-ready dict (``main`` writes it).

    ``texts`` (case name -> ``.rmc6f`` text) replaces the generated
    configurations: the staleness test passes the committed ones.
    """
    cases = []
    with tempfile.TemporaryDirectory() as directory:
        for index, (name, make_text, elements) in enumerate(CASES):
            text = texts[name] if texts is not None else make_text()
            path = Path(directory) / f"case{index}.rmc6f"
            path.write_text(text, encoding="utf-8")
            sites = load_site_displacements(path)
            expected = []
            for element in elements:
                result = displacement_modulations(sites, element)
                assert_margins(result, f"{name} {element}")
                expected.append(result)
            cases.append({"name": name, "text": text, "elements": elements, "expected": expected})
    enums = {
        "because": {name: list(codes) for name, codes in MODULATION_BECAUSE.items()},
        "kinds": list(MODE_KINDS),
        "rules": dict(MODULATION_RULES),
        "longWaveRatio": LONG_WAVE_RATIO,
    }
    return {"cases": cases, "enums": enums}


def main() -> None:
    OUT.write_text(json.dumps(build_fixture(), separators=(",", ":"), allow_nan=False), encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    main()
