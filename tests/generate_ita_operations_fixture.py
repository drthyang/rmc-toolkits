# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate spglib's space-group operations for the CIF export's ITA table.

Writes ``web_app/frontend/src/__tests__/fixtures/ita_operations_fixture.json``: for each of
the 230 space groups, every symmetry operation (centring included) of its ITA standard
setting as spglib's database gives it, as coordinate triplets. The setting is the one the
app's tables use (``itaOperations.js``, ``wyckoffTable.js``): origin choice 2 where ITA gives
two origins, hexagonal axes for the R groups, and spglib's first (ITA standard) setting
otherwise — unique axis b, cell choice 1 for the monoclinic groups.

``itaOperations.test.js`` closes the app's generators and requires the same sets, so the
origin the CIF export moves a structure to is ITA's own, checked against an independent
source. spglib is a development dependency only (``pip install spglib``); the app does not
use it:

    python tests/generate_ita_operations_fixture.py
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

import spglib

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "ita_operations_fixture.json"


def standard_hall_numbers() -> dict[int, tuple[int, str]]:
    """IT number -> (Hall number, spglib setting choice) of the setting the app uses."""
    settings: dict[int, list[tuple[int, str]]] = {}
    for hall in range(1, 531):
        kind = spglib.get_spacegroup_type(hall)
        settings.setdefault(kind.number, []).append((hall, kind.choice))
    chosen = {}
    for number, options in settings.items():
        by_choice = {choice: hall for hall, choice in options}
        if "2" in by_choice:
            chosen[number] = (by_choice["2"], "2")
        elif "H" in by_choice:
            chosen[number] = (by_choice["H"], "H")
        else:
            chosen[number] = options[0]
    return chosen


def triplet(rotation, translation) -> str:
    """'-y+1/2,x-y,z+1/4': the operation as ITA prints it, translation in [0, 1)."""
    parts = []
    for row, shift in zip(rotation, translation):
        text = ""
        for coefficient, axis in zip(row, "xyz"):
            if coefficient:
                sign = "-" if coefficient < 0 else "+"
                magnitude = "" if abs(coefficient) == 1 else str(abs(int(coefficient)))
                text += f"{sign}{magnitude}{axis}"
        fraction = Fraction(float(shift) % 1.0).limit_denominator(48) % 1
        if fraction:
            text += f"+{fraction.numerator}/{fraction.denominator}"
        parts.append(text.lstrip("+") or "0")
    return ",".join(parts)


def main() -> None:
    groups = {}
    for number, (hall, choice) in sorted(standard_hall_numbers().items()):
        symmetry = spglib.get_symmetry_from_database(hall)
        operations = sorted({
            triplet(rotation, translation)
            for rotation, translation in zip(symmetry["rotations"], symmetry["translations"])
        })
        groups[str(number)] = {"hall": hall, "choice": choice, "operations": operations}
    payload = {
        "description": "spglib ITA standard-setting operations (origin choice 2, hexagonal R) per IT number",
        "spglib": spglib.__version__,
        "groups": groups,
    }
    OUT.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {OUT} ({len(groups)} groups)")


if __name__ == "__main__":
    main()
