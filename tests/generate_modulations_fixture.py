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

Every ratio, significance and correlation a verdict decides on sits more than
``MARGIN`` from its threshold in every case, every polarisation is well
conditioned (its two largest eigenvalues more than ``POLARISATION_GAP`` of the
trace apart, and its sense not decided within ``AXIS_TIE_RTOL`` x 10 of a
tie), and every discrete choice of the diffuse test is clear of its tolerance
(rank ties, the lags the grid resolves, the strongest lag) and of the null's
(equal weights) (:func:`assert_margins`; the generator refuses to write
otherwise), so the two engines cannot disagree on a verdict or a direction
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

import numpy as np  # noqa: E402

import modulation_configs as configs  # noqa: E402
from cluster_configs import cluster_config_text, long_wave_config_text, polarised_config_text, wave_config_text  # noqa: E402
from rmc_toolkits.modulations import (  # noqa: E402
    AXIS_TIE_RTOL,
    DIFFUSE_LAGS,
    DIFFUSE_PHASES,
    DIFFUSE_RANK_RTOL,
    LONG_WAVE_RATIO,
    MODE_KINDS,
    MODULATION_BECAUSE,
    MODULATION_RULES,
    NULL_EQUAL_RTOL,
    PERCENTILE_RTOL,
    _folded,
    displacement_modulations,
    element_spectrum,
    independent_cell_null,
    is_long_wave,
)
from rmc_toolkits.pca_kde import SiteDisplacements, load_site_displacements  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "modulations_fixture.json"

#: The least distance between a decided ratio or sigma and its threshold.
MARGIN = 1e-6
#: The least gap between a polarisation's two largest eigenvalues, relative to
#: the trace: a rounding change of 1e-15 then moves its eigenvector by at most
#: ~1e-11, far inside the 1e-9 the engines and NumPy releases are held to.
POLARISATION_GAP = 1e-4

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
    # Plan item 4.3's known answer: both Zr sites displaced along the
    # direct-lattice [1 0 1] of a monoclinic cell at q = (1/2, 0, 0).
    ("polarised [1 0 1], monoclinic (8, 4, 4)", polarised_config_text, [None]),
    # F3's known answer: a weak Nb doubling whose (0, 0, 1/2) point sits at
    # 2.99 sigma over every element's points (against the independent-cell
    # null, F1) -- each element, and Nb named.
    ("weak doubling at the bar (4, 4, 4)", lambda: cluster_config_text("doubled", supercell=(4, 4, 4), delta=0.06, seed=33),
     [None, "Nb"]),
    # F1's known answer: Jahn-Teller octahedra whose six F share the long axis
    # per cell. An F mode five times the median read 3.66 sigma against the
    # null of independent sites; against the independent-cell null, 1.4.
    ("Jahn-Teller octahedra, a shared axis per cell (6, 6, 6)", lambda: configs.rmc6f_text(configs.jahn_teller(11)), [None]),
    # Plan item 4.4's known answers: B displacements correlated between
    # neighbouring cells, in antiphase along a and in phase along c, with no
    # single mode standing out (the A and B sites only, a small file).
    ("diffuse order in antiphase along a (6, 6, 6)",
     lambda: configs.rmc6f_text(configs.diffuse(1, supercell=(6, 6, 6), kappa=0.35, amp=0.08, anions=False)), [None, "Nb"]),
    ("diffuse order in phase along c (6, 6, 6)",
     lambda: configs.rmc6f_text(configs.diffuse(1, supercell=(6, 6, 6), kappa=0.35, amp=0.08, anions=False, sign=1.0, axis=2)),
     [None]),
]


def _diffuse_margins(sites: SiteDisplacements, result: dict, where: str) -> None:
    """The diffuse test's discrete choices sit clear of their tolerances in
    every element: the rank ties (each relative gap between sorted powers
    either an exact tie, below 1e-12, or above 1e-8, against PERCENTILE_RTOL
    = 1e-9), the lags the grid resolves (each template's residual after the
    earlier ones either below 1e-11 or above 1e-7 of its norm, against
    DIFFUSE_RANK_RTOL = 1e-9) and the strongest lag (the two largest
    |correlation| apart by more than ten times PERCENTILE_RTOL); and the null
    is clear of NULL_EQUAL_RTOL (its weights' spread either 0 or above 1e-9)."""
    supercell = np.asarray(result["supercell"], dtype=int)
    n_q = int(np.prod(supercell))
    grid = np.stack(np.unravel_index(np.arange(n_q), tuple(supercell)), axis=1)
    partner = np.ravel_multi_index(tuple(np.mod(-grid, supercell).T), tuple(supercell))
    q_frac = np.array([[_folded(index[axis], supercell[axis]) for axis in range(3)] for index in grid])
    compared = (np.arange(n_q) <= partner) & (np.arange(n_q) != 0) & ~np.array([is_long_wave(q) for q in q_frac])
    boundary = np.all((grid == 0) | (2 * grid == supercell), axis=1)
    for entry in result["elements"]:
        if not entry["sites"]:
            continue
        name = entry["element"]
        weights = independent_cell_null(sites, name)
        spread = float(weights.max() - weights.min()) / float(weights.max())
        if NULL_EQUAL_RTOL / 100 < spread < 1e-9:
            raise AssertionError(f"{where}: {name} null weights' spread {spread!r} sits at NULL_EQUAL_RTOL")
        power = element_spectrum(sites, name)[0].ravel()
        found = {tuple(round(value * size) % size for value, size in zip(mode["q"], supercell))
                 for mode in entry["strongest"] if mode["kind"] == "general" and mode["ordered"] is True}
        modes = [index for index in np.flatnonzero(compared & ~boundary) if tuple(grid[index]) not in found]
        values = np.sort(power[modes])
        gaps = np.diff(values) / np.maximum(np.abs(values[:-1]), 1e-300)
        if np.any((gaps >= 1e-12) & (gaps <= 1e-8)):
            raise AssertionError(f"{where}: {name} powers tie near PERCENTILE_RTOL ({gaps[(gaps >= 1e-12) & (gaps <= 1e-8)]!r})")
        q = q_frac[modes]
        basis = []
        for column in [np.ones(len(modes))] + [np.cos(2.0 * np.pi * (q @ np.asarray(lag, float))) for lag in DIFFUSE_LAGS]:
            norm = float(np.linalg.norm(column))
            for unit in basis:
                column = column - float(unit @ column) * unit
            rest = float(np.linalg.norm(column))
            if 1e-11 * norm <= rest <= 1e-7 * norm:
                raise AssertionError(f"{where}: {name} a lag's template is resolved near DIFFUSE_RANK_RTOL ({rest / norm!r})")
            if rest > DIFFUSE_RANK_RTOL * norm:
                basis.append(column / rest)
        magnitudes = sorted((abs(c) for c in entry["diffuse"]["correlations"] if c is not None), reverse=True)
        if len(magnitudes) > 1 and 1.0 - magnitudes[1] / magnitudes[0] <= 10 * PERCENTILE_RTOL:
            raise AssertionError(f"{where}: {name} the strongest lag is decided at the tie tolerance")


def _polarisation_margins(sites: SiteDisplacements, result: dict, where: str) -> None:
    """Every polarisation is well conditioned: the leading eigenvalue of the
    engine's own tensor at q (element_spectrum) apart from the next by more
    than POLARISATION_GAP of the trace, so its eigenvector is defined to far
    better than the parity tolerance; and its sense decided clear of the tie
    tolerance (the two largest |cos| to a, b, c either tie well inside
    AXIS_TIE_RTOL or differ by more than ten times it)."""
    cell = np.asarray(sites.unit_vectors, dtype=float)
    lengths = np.linalg.norm(cell, axis=1)
    supercell = np.asarray(result["supercell"], dtype=int)
    for entry in result["elements"]:
        if not entry["sites"]:
            continue
        tensor = element_spectrum(sites, entry["element"])[1]
        for point in entry["zoneBoundary"] + entry["strongest"]:
            pol = point["polarisation"]
            if pol is None:
                continue
            index = tuple(int(round(value * size)) % size for value, size in zip(point["q"], supercell))
            values = np.linalg.eigvalsh(tensor[index])
            if (values[-1] - values[-2]) <= POLARISATION_GAP * values.sum():
                raise AssertionError(f"{where}: {entry['element']} {point['q']} polarisation has no eigenvalue gap ({values!r})")
            magnitude = np.sort(np.abs(cell @ np.asarray(pol["direction"]) / lengths))[::-1]
            gap = 1.0 - magnitude[1] / magnitude[0]
            if AXIS_TIE_RTOL / 10.0 <= gap <= AXIS_TIE_RTOL * 10.0:
                raise AssertionError(f"{where}: {entry['element']} {point['q']} polarisation sense sits at the tie tolerance ({gap!r})")


def assert_margins(result: dict, where: str, sites: SiteDisplacements | None = None) -> None:
    """Every ratio and sigma a verdict compares with a threshold sits more than
    MARGIN from it; with ``sites``, every polarisation is well conditioned."""
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
        diffuse = entry["diffuse"]
        checks.append(("diffuse", "sigma", diffuse["sigma"], rules["diffuseSigma"]))
        if diffuse["correlation"] is not None:
            checks.append(("diffuse", "|correlation|", abs(diffuse["correlation"]), rules["diffuseCorrelation"]))
        for q, what, value, threshold in checks:
            if value is not None and abs(value - threshold) <= MARGIN:
                raise AssertionError(f"{where}: {entry['element']} {q} {what} = {value!r} lies within {MARGIN} of {threshold}")
    if sites is not None:
        _polarisation_margins(sites, result, where)
        _diffuse_margins(sites, result, where)


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
                assert_margins(result, f"{name} {element}", sites)
                expected.append(result)
            cases.append({"name": name, "text": text, "elements": elements, "expected": expected})
    enums = {
        "because": {name: list(codes) for name, codes in MODULATION_BECAUSE.items()},
        "kinds": list(MODE_KINDS),
        "rules": dict(MODULATION_RULES),
        "longWaveRatio": LONG_WAVE_RATIO,
        "axisTieRtol": AXIS_TIE_RTOL,
        "diffusePhases": list(DIFFUSE_PHASES),
        "diffuseRankRtol": DIFFUSE_RANK_RTOL,
        "nullEqualRtol": NULL_EQUAL_RTOL,
    }
    return {"cases": cases, "enums": enums}


def main() -> None:
    OUT.write_text(json.dumps(build_fixture(), separators=(",", ":"), allow_nan=False), encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size // 1024} KiB)")


if __name__ == "__main__":
    main()
