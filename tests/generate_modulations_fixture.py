# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the JS parity golden for the displacement-modulation engine.

Writes ``web_app/frontend/src/__tests__/fixtures/modulations_fixture.json``:
synthetic ``.rmc6f`` configurations (tests/cluster_configs.py) and the
reference results of ``rmc_toolkits.modulations.displacement_modulations``.
The vitest suite runs ``workers/modulations.js`` on the same text and asserts
integers exactly and floats to 1e-9. Regenerate whenever the engine's
conventions change:

    python tests/generate_modulations_fixture.py
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from cluster_configs import cluster_config_text, long_wave_config_text  # noqa: E402
from rmc_toolkits.modulations import displacement_modulations  # noqa: E402
from rmc_toolkits.pca_kde import load_site_displacements  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "modulations_fixture.json"

CASES = [
    ("doubled", (4, 4, 4), [None, "Nb"]),
    ("independent", (3, 4, 4), [None]),
    # Order under a stronger long wave (review 2026-10-12 #6).
    ("long-wave", (6, 6, 6), ["Ti"]),
]


def main() -> None:
    cases = []
    with tempfile.TemporaryDirectory() as directory:
        for mode, supercell, elements in CASES:
            text = (
                long_wave_config_text(supercell=supercell) if mode == "long-wave"
                else cluster_config_text(mode, supercell=supercell)
            )
            path = Path(directory) / f"{mode}.rmc6f"
            path.write_text(text, encoding="utf-8")
            sites = load_site_displacements(path)
            cases.append({
                "name": f"{mode} {supercell}",
                "text": text,
                "elements": elements,
                "expected": [displacement_modulations(sites, element) for element in elements],
            })
    OUT.write_text(json.dumps({"cases": cases}, separators=(",", ":")), encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    main()
