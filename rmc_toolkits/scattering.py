# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Neutron scattering lengths and Faber-Ziman coefficient calculator.

``COHERENT_B_FM`` holds bound coherent neutron scattering lengths for natural
elements in femtometres, from the NIST NCNR table
(https://www.ncnr.nist.gov/resources/n-lengths/, after Sears, Neutron News 3,
26 (1992)). For the strong absorbers with complex b (B, Cd, Dy, Eu, Gd, In,
Sm) the REAL part is stored, the standard choice for Faber-Ziman weighting.

The classic stog input's "Faber-Ziman coefficient" is ``<b>^2 = (sum_i c_i
b_i)^2`` in **barns** (1 barn = 100 fm^2); Keen 2001 Eq. 14's Q->0 limit uses
the different number ``<b^2> = sum_i c_i b_i^2``. Both are computed here, in
both unit systems, because the ecosystem is inconsistent (e.g. pystog's argon
example config quotes ``<b_coh>^2`` in fm^2).
"""

from __future__ import annotations

import re
from dataclasses import dataclass

COMPLEX_B_ELEMENTS = ("B", "Cd", "Dy", "Eu", "Gd", "In", "Sm")

COHERENT_B_FM: dict[str, float] = {
    "Ag": 5.922, "Al": 3.449, "Am": 8.3, "Ar": 1.909, "As": 6.58, "Au": 7.63,
    "B": 5.3, "Ba": 5.07, "Be": 7.79, "Bi": 8.532, "Br": 6.795, "C": 6.646,
    "Ca": 4.7, "Cd": 4.87, "Ce": 4.84, "Cl": 9.577, "Co": 2.49, "Cr": 3.635,
    "Cs": 5.42, "Cu": 7.718, "Dy": 16.9, "Er": 7.79, "Eu": 7.22, "F": 5.654,
    "Fe": 9.45, "Ga": 7.288, "Gd": 6.5, "Ge": 8.185, "H": -3.739, "He": 3.26,
    "Hf": 7.7, "Hg": 12.692, "Ho": 8.01, "I": 5.28, "In": 4.065, "Ir": 10.6,
    "K": 3.67, "Kr": 7.81, "La": 8.24, "Li": -1.9, "Lu": 7.21, "Mg": 5.375,
    "Mn": -3.73, "Mo": 6.715, "N": 9.36, "Na": 3.63, "Nb": 7.054, "Nd": 7.69,
    "Ne": 4.566, "Ni": 10.3, "Np": 10.55, "O": 5.803, "Os": 10.7, "P": 5.13,
    "Pa": 9.1, "Pb": 9.405, "Pd": 5.91, "Pm": 12.6, "Pr": 4.58, "Pt": 9.6,
    "Ra": 10.0, "Rb": 7.09, "Re": 9.2, "Rh": 5.88, "Ru": 7.03, "S": 2.847,
    "Sb": 5.57, "Sc": 12.29, "Se": 7.97, "Si": 4.1491, "Sm": 0.8, "Sn": 6.225,
    "Sr": 7.02, "Ta": 6.91, "Tb": 7.38, "Tc": 6.8, "Te": 5.8, "Th": 10.31,
    "Ti": -3.438, "Tl": 8.776, "Tm": 7.07, "U": 8.417, "V": -0.3824, "W": 4.86,
    "Xe": 4.92, "Y": 7.75, "Yb": 12.43, "Zn": 5.68, "Zr": 7.16,
}

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)|(\()|(\))(\d*\.?\d*)")


def parse_formula(formula: str) -> dict[str, float]:
    """Parse a chemical formula into element counts.

    Supports decimal stoichiometries and parentheses:
    ``"GaTa4Se8"``, ``"Sr0.5Ba0.5TiO3"``, ``"Ca(OH)2"``, ``"Al2(SO4)3"``.
    """
    compact = formula.replace(" ", "")
    if not compact:
        raise ValueError("empty formula")

    def parse_group(i: int) -> tuple[dict[str, float], int]:
        counts: dict[str, float] = {}
        while i < len(compact):
            char = compact[i]
            if char == ")":
                return counts, i
            if char == "(":
                inner, j = parse_group(i + 1)
                if j >= len(compact) or compact[j] != ")":
                    raise ValueError(f"unbalanced parentheses in {formula!r}")
                j += 1
                m = re.match(r"\d*\.?\d*", compact[j:])
                mult = float(m.group()) if m.group() else 1.0
                j += len(m.group())
                for el, n in inner.items():
                    counts[el] = counts.get(el, 0.0) + n * mult
                i = j
                continue
            m = re.match(r"([A-Z][a-z]?)(\d*\.?\d*)", compact[i:])
            if not m:
                raise ValueError(f"cannot parse formula {formula!r} at {compact[i:]!r}")
            el = m.group(1)
            count = float(m.group(2)) if m.group(2) else 1.0
            counts[el] = counts.get(el, 0.0) + count
            i += m.end()
        return counts, i

    counts, i = parse_group(0)
    if i != len(compact):
        raise ValueError(f"unbalanced parentheses in {formula!r}")
    if not counts:
        raise ValueError(f"no elements found in {formula!r}")
    return counts


@dataclass(frozen=True)
class FaberZiman:
    """Composition-derived Faber-Ziman coefficients.

    ``b_avg_sq_barn`` is the classic stog "Faber-Ziman coefficient" input line
    and :class:`rmc_toolkits.scaling.ScalingConfig`'s ``b_avg_sq``;
    ``b_sq_avg_barn`` is its ``b_sq_avg`` (Keen Eq. 14 diagnostic).
    """

    fractions: dict[str, float]
    b_coh_fm: dict[str, float]
    b_avg_sq_fm2: float
    b_sq_avg_fm2: float
    weights: dict[tuple[str, str], float]

    @property
    def b_avg_sq_barn(self) -> float:
        return self.b_avg_sq_fm2 / 100.0

    @property
    def b_sq_avg_barn(self) -> float:
        return self.b_sq_avg_fm2 / 100.0


def faber_ziman(
    composition: str | dict[str, float],
    b_overrides_fm: dict[str, float] | None = None,
) -> FaberZiman:
    """Compute Faber-Ziman coefficients for a formula or counts dict.

    ``b_overrides_fm`` substitutes per-element scattering lengths (fm), e.g.
    for isotopically enriched samples.
    """
    counts = parse_formula(composition) if isinstance(composition, str) else dict(composition)
    if not counts:
        raise ValueError("empty composition")
    total = float(sum(counts.values()))
    if total <= 0:
        raise ValueError("composition counts must sum to a positive number")

    b_values: dict[str, float] = {}
    for el in counts:
        if b_overrides_fm and el in b_overrides_fm:
            b_values[el] = float(b_overrides_fm[el])
        elif el in COHERENT_B_FM:
            b_values[el] = COHERENT_B_FM[el]
        else:
            raise ValueError(
                f"no coherent scattering length for element {el!r}; "
                "supply it via b_overrides_fm"
            )

    fractions = {el: counts[el] / total for el in counts}
    b_avg = sum(fractions[el] * b_values[el] for el in counts)
    b_avg_sq = b_avg * b_avg
    b_sq_avg = sum(fractions[el] * b_values[el] ** 2 for el in counts)
    if b_avg_sq < 1e-4 * b_sq_avg:
        raise ValueError(
            "near-null-matrix composition: <b>^2 is negligible relative to "
            "<b^2>, so Faber-Ziman normalization is undefined (e.g. Ti/Zr "
            "null matrices); adjust the composition or b_overrides_fm"
        )
    weights = {
        (i, j): fractions[i] * fractions[j] * b_values[i] * b_values[j] / b_avg_sq
        for i in counts
        for j in counts
    }
    return FaberZiman(
        fractions=fractions,
        b_coh_fm=b_values,
        b_avg_sq_fm2=b_avg_sq,
        b_sq_avg_fm2=b_sq_avg,
        weights=weights,
    )
