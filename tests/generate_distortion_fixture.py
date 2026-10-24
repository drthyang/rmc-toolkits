# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the JS parity golden for the polyhedral-distortion engine.

Writes ``web_app/frontend/src/__tests__/fixtures/distortion_fixture.json``:
synthetic ``.rmc6f`` configurations (tests/distortion_configs.py), the queries
run on each, and the reference results of
``rmc_toolkits.distortion.polyhedral_distortion`` (or the refusal's message),
with the engine's constants. The vitest suite (``distortion.test.js``) parses
the same text as the worker does (``siteDisplacementsFromRmc6f`` for the
atoms, ``structureFromRmc6f`` for the average-structure basis), runs
``workers/distortion.js`` and compares every value by the golden rule
(goldenCompare.js: integers, strings, booleans and nulls exactly, other
numbers to 1e-9 relative + 1e-12); ``tests/test_distortion_fixture.py`` fails
when the committed file no longer matches the Python engine.

The two engines take every decision on the same numbers only where those
numbers agree bitwise (the configuration's bond vectors) or sit far from the
decision's threshold, so the generator refuses to write a case in which a bond
length lies within ``MARGIN_A`` of a window bound (the configuration's or the
average structure's, whose site means agree to ~1e-16 between the runtimes),
or an octahedron's largest trans-matching sum lies within ``MARGIN_DEG`` of the
next (:func:`assert_margins`). Regenerate from the checkout under test whenever
the engine's conventions change:

    PYTHONPATH=. python tests/generate_distortion_fixture.py
"""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import distortion_configs as configs  # noqa: E402
from rmc_toolkits.distortion import (  # noqa: E402
    APP_MAX_BONDS,
    OCTAHEDRON_MATCHINGS,
    PAIRS,
    QUANTILES,
    TETRAHEDRAL_ANGLE_DEG,
    TRANS_TIE_DEG,
    _gather,
    polyhedral_distortion,
    read_distortion_configuration,
    site_species,
)
from rmc_toolkits.triplets import WINDOW_TOL, _neighbor_bonds, _sort_by_center  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "distortion_fixture.json"

#: The least distance between a bond length and a window bound (angstrom).
MARGIN_A = 1e-7
#: The least gap between an octahedron's largest trans-matching angle sum and the next (degrees).
MARGIN_DEG = 1e-6

SHEARED = [[3.9, 0.0, 0.0], [0.25, 3.9, 0.0], [0.15, 0.2, 3.95]]


def _mixed_b_site(rng, cell):
    return "Zr" if (cell[0] + 2 * cell[1] + 3 * cell[2]) % 4 == 0 else "Ti"


def _tie():
    """Two Cu, one with 4 F and one with 2: a tie, which goes to the lower count."""
    sites = [("Cu", (0.25, 0.5, 0.5)), ("Cu", (0.75, 0.5, 0.5))]
    sites += [("F", tuple(np.array([0.25, 0.5, 0.5]) + 0.2 * v)) for v in configs.AXES[:4]]
    sites += [("F", tuple(np.array([0.75, 0.5, 0.5]) + 0.2 * v)) for v in configs.AXES[:2]]
    return configs.supercell_atoms(np.eye(3) * 10.0, sites, (1, 1, 1))


def _random_triclinic():
    """48 atoms at random in a triclinic box: every coordination at once."""
    rng = np.random.default_rng(11)
    lattice = np.array([[6.0, 0.0, 0.0], [3.0, 5.0, 0.0], [1.0, 1.0, 7.0]])
    atoms = []
    for number in range(1, 49):
        element = "Se" if rng.uniform() < 0.6 else "Nb"
        atoms.append({
            "atom_number": number, "element": element, "type_label": element,
            "coords": rng.uniform(size=3), "reference_number": number, "cell_indices": np.zeros(3, dtype=int),
        })
    return atoms, lattice, np.array([1, 1, 1])


def _text(built, **options):
    return configs.rmc6f_text(*built, **options)


def query(centre, ligand, r_min, r_max, max_bonds=None):
    entry = {"centre": centre, "ligand": ligand, "rMin": r_min, "rMax": r_max}
    if max_bonds is not None:
        entry["maxBonds"] = max_bonds
    return entry


CASES = [
    # Orientationally disordered Jahn-Teller octahedra with thermal noise: a
    # 4 + 2 configuration whose average octahedron is regular; the first peak
    # alone (4 F, near-flat: a large <lambda>), the reverse pair (one Cu per F:
    # no polyhedra) and the work budget.
    ("jahn-teller", lambda: _text(configs.jahn_teller_network(sigma=0.03, seed=7)),
     [query("Cu", "F", 1.6, 2.6), query("Cu", "F", 1.6, 2.1), query("f", "cu", 1.6, 2.6),
      query("Cu", "F", 1.6, 2.6, max_bonds=100), query("Cu", "F", 0.5, 1)]),
    # Corner-sharing TiO6 with thermal motion of every atom: octahedra, the
    # SrO12 cuboctahedra (DI only), linear Ti-O-Ti (n = 2), an unknown element
    # and an inverted window.
    ("perovskite-thermal", lambda: _text(configs.perovskite(sigma_o=0.05, sigma_cation=0.03, seed=11)),
     [query("Ti", "O", 1.5, 2.5), query("Sr", "O", 2.2, 3.4), query("O", "Ti", 1.5, 2.5),
      query("Ti", "Xx", 1.5, 2.5), query("Ti", "O", 2.5, 1.5)]),
    # A sheared (triclinic) cell: a static distortion the average carries too.
    ("perovskite-triclinic", lambda: _text(configs.perovskite(unit_cell=SHEARED, sigma_o=0.02, seed=3)),
     [query("Ti", "O", 1.5, 2.5)]),
    # Corner-sharing tetrahedra both ways.
    ("zincblende", lambda: _text(configs.zincblende(sigma=0.06)),
     [query("Zn", "S", 1.9, 2.9), query("S", "Zn", 1.9, 2.9)]),
    # A mixed B site (Ti majority): Zr-centred polyhedra have no average site.
    ("mixed-site", lambda: _text(configs.perovskite(b_site=_mixed_b_site, sigma_o=0.03, seed=2)),
     [query("Ti", "O", 1.5, 2.5), query("Zr", "O", 1.5, 2.5)]),
    # Legacy coords-only lines: no reference numbers, no average structure.
    ("coords-only", lambda: _text(configs.perovskite(sigma_o=0.03, seed=4), coords_only=True),
     [query("Ti", "O", 1.5, 2.5)]),
    # One cell smaller than the window: each TiO6 through periodic images.
    ("small-box", lambda: _text(configs.perovskite(supercell=(1, 1, 1), sigma_o=0.05, seed=6)),
     [query("Ti", "O", 1.5, 2.5)]),
    ("coordination-tie", lambda: _text(_tie()), [query("Cu", "F", 1.0, 2.5)]),
    # Irregular polyhedra of every kind: tetrahedra and octahedra far from
    # regular (the trans matching on arbitrary angles), and 18 neighbours.
    ("random-triclinic", lambda: _text(_random_triclinic()),
     [query("Nb", "Se", 0.5, 2.0), query("Se", "Se", 0.5, 2.4), query("Se", "Nb", 0.5, 2.5), query("Nb", "Se", 1.0, 3.4)]),
]


def _octahedron_gaps(vectors: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Each octahedron's largest trans-matching angle sum minus the next."""
    pairs = np.asarray(PAIRS[6])
    one, two = vectors[:, pairs[:, 0]], vectors[:, pairs[:, 1]]
    cosine = ((one[..., 0] * two[..., 0] + one[..., 1] * two[..., 1]) + one[..., 2] * two[..., 2]) / (
        lengths[:, pairs[:, 0]] * lengths[:, pairs[:, 1]]
    )
    angles = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    matchings = np.asarray(OCTAHEDRON_MATCHINGS)
    sums = np.sort((angles[:, matchings[:, 0]] + angles[:, matchings[:, 1]]) + angles[:, matchings[:, 2]], axis=1)
    return sums[:, -1] - sums[:, -2]


def _searches(config, centre: str, ligand: str):
    """(centre positions, centre rows, ligand positions, ligand rows, lattice) of
    the configuration and of the average structure (None without one)."""
    symbols = np.asarray([str(element).strip().capitalize() for element in config.elements])
    fractional = np.asarray(config.fractional, dtype=float)
    wrapped = fractional - np.floor(fractional)
    centres, ligands = np.flatnonzero(symbols == centre), np.flatnonzero(symbols == ligand)
    searches = [(wrapped[centres], centres, wrapped[ligands], ligands, np.asarray(config.lattice_vectors, dtype=float))]
    if config.basis:
        labels = [site_species(row["elementCounts"]) for row in config.basis]
        sites_c = np.asarray([k for k, label in enumerate(labels) if label == centre], dtype=int)
        sites_l = np.asarray([k for k, label in enumerate(labels) if label == ligand], dtype=int)
        positions = np.asarray([row["mean"] for row in config.basis], dtype=float)
        means = positions - np.floor(positions)
        unit = np.asarray(config.lattice_vectors, dtype=float) / np.maximum(np.asarray(config.supercell, dtype=float), 1.0)[:, None]
        if sites_c.size:
            searches.append((means[sites_c], sites_c, means[sites_l], sites_l, unit))
    return searches


def assert_margins(name: str, config, spec: dict) -> None:
    """Refuse a query whose decisions sit at their thresholds (see the module doc)."""
    centre, ligand = spec["centre"].strip().capitalize(), spec["ligand"].strip().capitalize()
    lo, hi = float(spec["rMin"]), float(spec["rMax"])
    if hi <= lo or centre not in set(config.elements) or ligand not in set(config.elements):
        return  # a refusal: the engines refuse before measuring anything
    for centres_at, centres, ligands_at, ligands, lattice in _searches(config, centre, ligand):
        def count(window):
            table = _neighbor_bonds(centres_at, centres, ligands_at, ligands, lattice, [window], count_only=True)
            return table[:, 1]
        wide = count((max(lo - MARGIN_A, 0.0), hi + MARGIN_A))
        narrow = count((lo + MARGIN_A, hi - MARGIN_A))
        if not np.array_equal(wide, narrow):
            raise SystemExit(f"{name} {spec}: a bond length lies within {MARGIN_A} A of a window bound (+/- {WINDOW_TOL})")
        bonds = _neighbor_bonds(centres_at, centres, ligands_at, ligands, lattice, [(lo, hi)])
        grouped, sizes, starts = _sort_by_center(bonds, centres.size)
        six = np.flatnonzero(sizes == 6)
        if six.size:
            gaps = _octahedron_gaps(*_gather(grouped, starts, six, 6))
            if gaps.min() < MARGIN_DEG:
                raise SystemExit(f"{name} {spec}: an octahedron's trans matchings tie within {gaps.min():.3g} deg")


def evaluate(name: str, text: str, queries: list[dict]) -> list[dict]:
    with tempfile.TemporaryDirectory() as folder:
        path = Path(folder) / f"{name}.rmc6f"
        path.write_text(text, encoding="utf-8")
        config = read_distortion_configuration(path)
    expected = []
    for spec in queries:
        assert_margins(name, config, spec)
        try:
            result = polyhedral_distortion(
                config, spec["centre"], spec["ligand"], spec["rMin"], spec["rMax"], max_bonds=spec.get("maxBonds")
            )
        except ValueError as error:
            result = {"error": str(error)}
        expected.append(json.loads(json.dumps(result, allow_nan=False)))
    return expected


def constants() -> dict:
    return {
        "APP_MAX_BONDS": APP_MAX_BONDS,
        "WINDOW_TOL": WINDOW_TOL,
        "TRANS_TIE_DEG": TRANS_TIE_DEG,
        "TETRAHEDRAL_ANGLE_DEG": TETRAHEDRAL_ANGLE_DEG,
        "QUANTILES": list(QUANTILES),
        "PAIRS": {str(n): [list(pair) for pair in pairs] for n, pairs in PAIRS.items()},
        "OCTAHEDRON_MATCHINGS": [list(matching) for matching in OCTAHEDRON_MATCHINGS],
    }


def build_fixture(texts: dict[str, str] | None = None) -> dict:
    """The fixture; ``texts`` (case name -> .rmc6f text) reuses committed
    configurations instead of building them (the staleness test)."""
    cases = []
    for name, build, queries in CASES:
        text = texts[name] if texts is not None else build()
        cases.append({"name": name, "text": text, "queries": queries, "expected": evaluate(name, text, queries)})
    return {
        "about": [
            "Golden for workers/distortion.js: rmc_toolkits.distortion.polyhedral_distortion on each case's",
            ".rmc6f text (or its refusal's message). Written by tests/generate_distortion_fixture.py;",
            "compared by goldenCompare.js (distortion.test.js) and tests/golden_compare.py",
            "(tests/test_distortion_fixture.py). Do not edit by hand.",
        ],
        "constants": constants(),
        "cases": cases,
    }


def main() -> None:
    fixture = build_fixture()
    OUT.write_text(json.dumps(fixture, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    queries = sum(len(case["queries"]) for case in fixture["cases"])
    print(f"wrote {OUT.relative_to(ROOT)}: {len(fixture['cases'])} cases, {queries} queries")


if __name__ == "__main__":
    main()
