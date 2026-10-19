# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""The average-structure basis (``rmc_toolkits/structure_basis.py``) and ``/api/structure``'s ``basis``.

Each convention ``structureFromRmc6f()`` (``browserData.js``) follows is pinned
here on a configuration small enough to work by hand -- the circular-mean
centre of a site that straddles the cell edge, the arithmetic mean unwrapped
about it, the population (divide by n) covariance, the rms displacement through
the full metric of an oblique cell, the first atom's element for a mixed site,
coords-only atoms left out, sites sorted by reference number, and
``Math.round``'s ties -- so the parity golden
(``tests/generate_structure_basis_fixture.py``) holds two runtimes that agree
on purpose, not by luck. The route tests check that ``/api/structure`` returns
the basis of the same read, an empty one for a coords-only file, a 400 for a
basis that is not finite (a computed result, never nulled), and above the
space-group finder's site cap (``MAX_BASIS_SITES`` = ``MAX_SYMMETRY_SITES``,
pinned here against symmetryModel.js) only ``basisOmitted: {sites, cap}``.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
import unittest

import numpy as np

from rmc_toolkits import MAX_BASIS_SITES, basis_site_count, read_structure_basis, structure_basis
from rmc_toolkits.parsers import parse_rmc6f_atoms
from rmc_toolkits.structure_basis import _round_half_up

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "structure_basis_fixture.json"
SYMMETRY_MODEL = ROOT / "web_app" / "frontend" / "src" / "symmetryModel.js"
KEYS = ["count", "covFrac", "dispA", "el", "elementCounts", "frac", "mean", "referenceNumber"]


def _one_copy_box(sites: int) -> str:
    """A 1x1x1 box of `sites` atoms, each its own reference site (a glass, or an
    imported configuration): the case the finder's site cap is for."""
    lines = []
    for number in range(1, sites + 1):
        x, y, z = ((number * 0.6180339887) % 1.0, (number * 0.4142135623) % 1.0, (number * 0.7320508075) % 1.0)
        lines.append(f"{number} {'Si' if number % 3 else 'O'} [{number}] {x:.10f} {y:.10f} {z:.10f} {number} 0 0 0")
    return _config(lines, supercell=(1, 1, 1), lattice=((40.0, 0.0, 0.0), (0.0, 40.0, 0.0), (0.0, 0.0, 40.0)))


def _config(lines, supercell=(2, 1, 1), lattice=((20.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0))):
    header = [
        "(Version 6f format configuration file)",
        f"Number of atoms: {len(lines)}",
        "Supercell dimensions: {} {} {}".format(*supercell),
        "Lattice vectors (Ang):",
        *(" ".join(f"{value:.10f}" for value in row) for row in lattice),
        "Atoms:",
    ]
    return "\n".join(header + list(lines)) + "\n"


def _basis(text: str) -> list[dict]:
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "config.rmc6f"
        path.write_text(text, encoding="utf-8")
        return read_structure_basis(path)


def _line(number, element, reference, x, y, z, cell=(0, 0, 0)):
    return f"{number} {element} [{reference}] {x!r} {y!r} {z!r} {reference} {cell[0]} {cell[1]} {cell[2]}"


class StructureBasisConventionTests(unittest.TestCase):
    def test_two_copies_give_the_population_moments(self):
        # A 2x1x1 supercell of a 10 A cubic cell: copies at within-cell 0.5 +/- 0.02.
        text = _config([
            _line(1, "Ga", 1, (0.5 + 0.02) / 2, 0.25, 0.75),
            _line(2, "Ga", 1, (0.5 - 0.02 + 1) / 2, 0.25, 0.75, (1, 0, 0)),
        ])
        (site,) = _basis(text)
        self.assertEqual(sorted(site), KEYS)
        self.assertEqual((site["el"], site["referenceNumber"], site["count"]), ("Ga", 1, 2))
        self.assertEqual(site["elementCounts"], {"Ga": 2})
        np.testing.assert_allclose(site["frac"], [0.5, 0.25, 0.75], atol=1e-12)
        np.testing.assert_allclose(site["mean"], [0.5, 0.25, 0.75], atol=1e-12)
        # Divided by n (0.02^2), not n - 1 (twice that): a one-copy site is 0, not NaN.
        self.assertAlmostEqual(site["covFrac"][0][0], 0.02**2, delta=1e-15)
        self.assertEqual(site["covFrac"][1][1], 0.0)
        # +/-0.02 of a 10 A edge: 0.2 A about the mean.
        self.assertAlmostEqual(site["dispA"], 0.2, delta=1e-12)

    def test_a_site_straddling_the_cell_edge_is_centred_on_it(self):
        # Within-cell 0.99 and 0.01: the circular mean is 0 (= 1), not 0.5.
        (site,) = _basis(_config([
            _line(1, "Ga", 1, 0.495, 0.25, 0.75),
            _line(2, "Ga", 1, 0.505, 0.25, 0.75, (1, 0, 0)),
        ]))
        self.assertLess(min(site["frac"][0], 1.0 - site["frac"][0]), 1e-12)
        self.assertLess(abs(site["mean"][0] - round(site["mean"][0])), 1e-12)
        self.assertAlmostEqual(site["dispA"], 0.1, delta=1e-12)

    def test_the_mean_is_arithmetic_about_the_circular_centre(self):
        # Within-cell 0.98, 0.99 and 0.06 unwrap to -0.02, -0.01 and 0.06: mean 0.01,
        # while the circular mean of the three leans elsewhere.
        (site,) = _basis(_config(
            [_line(1, "Se", 1, w, 0.5, 0.5) for w in (0.98, 0.99, 0.06)], supercell=(1, 1, 1),
            lattice=((10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0)),
        ))
        self.assertAlmostEqual(site["mean"][0], 0.01, delta=1e-12)
        self.assertGreater(abs(site["frac"][0] - 0.01), 1e-6)
        offsets = np.array([-0.02, -0.01, 0.06])
        self.assertAlmostEqual(site["covFrac"][0][0], offsets.var(), delta=1e-15)

    def test_the_rms_displacement_uses_the_full_metric(self):
        # Hexagonal cell (gamma = 120 deg): copies at +/-(0.01, 0.01, 0) move by
        # 0.01 (a + b), |a + b| = 10 A, so 0.1 A -- not 0.01 * sqrt(2) * 10.
        root3 = math.sqrt(3.0)
        lattice = ((10.0, 0.0, 0.0), (-5.0, 5.0 * root3, 0.0), (0.0, 0.0, 12.0))
        (site,) = _basis(_config(
            [_line(1, "Ga", 1, 0.31, 0.42, 0.5), _line(2, "Ga", 1, 0.33, 0.44, 0.5)],
            supercell=(1, 1, 1), lattice=lattice,
        ))
        self.assertAlmostEqual(site["dispA"], 0.1, delta=1e-12)
        np.testing.assert_allclose(site["covFrac"], [[1e-4, 1e-4, 0.0], [1e-4, 1e-4, 0.0], [0.0, 0.0, 0.0]], atol=1e-15)

    def test_a_mixed_site_is_its_first_atom_species_with_every_count(self):
        lines = [
            _line(1, "Ta", 5, 0.1, 0.2, 0.3),
            _line(2, "Nb", 5, 0.1, 0.2, 0.3),
            _line(3, "Nb", 5, 0.1, 0.2, 0.3),
        ]
        (site,) = _basis(_config(lines, supercell=(1, 1, 1)))
        self.assertEqual(site["el"], "Ta")
        self.assertEqual(site["count"], 3)
        self.assertEqual(list(site["elementCounts"].items()), [("Ta", 1), ("Nb", 2)])

    def test_sites_are_sorted_by_reference_number_whatever_the_file_order(self):
        lines = [
            _line(1, "Se", 7, 0.7, 0.1, 0.1),
            _line(2, "Ga", 2, 0.2, 0.2, 0.2),
            _line(3, "Se", 7, 0.7, 0.1, 0.1),
            _line(4, "Ta", 4, 0.4, 0.4, 0.4),
        ]
        basis = _basis(_config(lines, supercell=(1, 1, 1)))
        self.assertEqual([site["referenceNumber"] for site in basis], [2, 4, 7])
        self.assertEqual([site["count"] for site in basis], [1, 1, 2])

    def test_a_one_copy_site_has_no_spread(self):
        (site,) = _basis(_config([_line(1, "Te", 3, 0.123, 0.456, 0.789)], supercell=(1, 1, 1)))
        self.assertEqual(site["covFrac"], [[0.0] * 3] * 3)
        self.assertEqual(site["dispA"], 0.0)
        self.assertEqual(site["mean"], site["frac"])

    def test_coords_only_atoms_are_left_out(self):
        full = _line(1, "Ga", 1, 0.25, 0.25, 0.25)
        coords = ["2 O 0.5 0.5 0.5", "3 O 0.75 0.1 0.2"]
        basis = _basis(_config([full, *coords], supercell=(1, 1, 1)))
        self.assertEqual([(site["referenceNumber"], site["count"]) for site in basis], [(1, 1)])
        self.assertEqual(_basis(_config(coords, supercell=(1, 1, 1))), [])
        self.assertEqual(structure_basis([], np.eye(3), np.ones(3)), [])

    def test_the_parsed_atoms_and_the_file_give_one_basis(self):
        text = _config([_line(1, "Ga", 1, 0.495, 0.25, 0.75), _line(2, "Ga", 1, 0.505, 0.25, 0.75, (1, 0, 0))])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.rmc6f"
            path.write_text(text, encoding="utf-8")
            atoms, _ = parse_rmc6f_atoms(path)
            lattice = np.array([[20.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]])
            self.assertEqual(structure_basis(atoms, lattice, np.array([2.0, 1.0, 1.0])), read_structure_basis(path))

    def test_the_site_count_is_the_basis_length_without_the_rows(self):
        fixture = json.loads(FIXTURE.read_text(encoding="utf-8"))
        for text in (fixture["synthetic"], fixture["coordsOnly"], _one_copy_box(25)):
            with tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "config.rmc6f"
                path.write_text(text, encoding="utf-8")
                atoms, _ = parse_rmc6f_atoms(path)
                self.assertEqual(basis_site_count(atoms), len(read_structure_basis(path)))

    def test_the_cap_is_the_symmetry_finders(self):
        # MAX_SYMMETRY_SITES in symmetryModel.js: above it the finder declines,
        # so /api/structure sends no basis. (The JS parity test holds the JS
        # constant to the golden's maxBasisSites too.)
        match = re.search(r"^export const MAX_SYMMETRY_SITES = (\d+);$", SYMMETRY_MODEL.read_text(encoding="utf-8"), re.MULTILINE)
        self.assertIsNotNone(match)
        self.assertEqual(int(match.group(1)), MAX_BASIS_SITES)
        self.assertEqual(MAX_BASIS_SITES, 2000)

    def test_offsets_round_like_math_round(self):
        # Ties towards +infinity, and 0.49999999999999994 (where floor(x + 0.5) is 1) to 0.
        values = np.array([-1.5, -0.5, -0.49999999999999994, 0.49999999999999994, 0.5, 1.5, 2.5])
        np.testing.assert_array_equal(_round_half_up(values), [-1.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0])


# --- /api/structure ------------------------------------------------------------

if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402


class StructureBasisEndpointTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        cls.fixture = json.loads(FIXTURE.read_text(encoding="utf-8"))

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(dir=ROOT)
        self.directory = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def _structure(self, text: str, **query):
        path = self.directory / "run.rmc6f"
        path.write_text(text, encoding="utf-8")
        return path, self.client.get("/api/structure", query_string={"dir": str(self.directory), **query})

    def test_the_route_returns_the_basis_of_its_own_read(self):
        path, response = self._structure(self.fixture["synthetic"], maxPoints=100)
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["basis"], read_structure_basis(path))
        self.assertIsNone(payload["basisOmitted"])
        self.assertEqual([sorted(site) for site in payload["basis"]], [KEYS] * len(payload["basis"]))
        # The basis sees every copy, whatever maxPoints samples for the points.
        self.assertEqual(sum(site["count"] for site in payload["basis"]), payload["totalAtoms"] - 3)
        self.assertLess(payload["sampledAtoms"], payload["totalAtoms"])

    def test_the_demo_run_gives_its_52_sites(self):
        response = self.client.get(
            "/api/structure", query_string={"dir": "web_app/frontend/public/demo", "maxPoints": 100}
        )
        self.assertEqual(response.status_code, 200)
        basis = response.get_json()["basis"]
        self.assertEqual(len(basis), 52)
        self.assertEqual(basis, read_structure_basis(ROOT / "web_app" / "frontend" / "public" / "demo" / "GTS_250K.rmc6f"))
        self.assertEqual({site["count"] for site in basis}, {1000})

    def test_a_coords_only_file_has_an_empty_basis(self):
        _, response = self._structure(self.fixture["coordsOnly"])
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertEqual(payload["basis"], [])
        self.assertIsNone(payload["basisOmitted"])
        self.assertEqual(payload["totalAtoms"], 8)

    def test_above_the_finders_cap_only_the_site_count_is_sent(self):
        # 2001 one-copy sites: nothing in the app can use the rows (the card
        # says "not analysed"), so the route sends their number instead.
        _, response = self._structure(_one_copy_box(MAX_BASIS_SITES + 1), maxPoints=100)
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertIsNone(payload["basis"])
        self.assertEqual(payload["basisOmitted"], {"sites": MAX_BASIS_SITES + 1, "cap": MAX_BASIS_SITES})
        # The rest of the model is there as before.
        self.assertEqual(payload["totalAtoms"], MAX_BASIS_SITES + 1)
        self.assertEqual(len(payload["points"]), payload["sampledAtoms"])
        self.assertEqual(sum(len(refs) for refs in payload["atomIndices"].values()), MAX_BASIS_SITES + 1)

    def test_at_the_cap_the_basis_is_sent(self):
        path, response = self._structure(_one_copy_box(MAX_BASIS_SITES), maxPoints=100)
        self.assertEqual(response.status_code, 200)
        payload = response.get_json()
        self.assertIsNone(payload["basisOmitted"])
        self.assertEqual(len(payload["basis"]), MAX_BASIS_SITES)
        self.assertEqual(payload["basis"], read_structure_basis(path))

    def test_a_basis_that_is_not_finite_is_a_400(self):
        # Finite in the file, but x * N overflows: NaN positions are refused, not nulled.
        text = _config([
            _line(1, "Ga", 1, 0.25, 0.25, 0.25),
            _line(2, "Ga", 1, 1.7e308, 0.25, 0.25, (1, 0, 0)),
        ])
        with np.errstate(all="ignore"):
            _, response = self._structure(text)
        self.assertEqual(response.status_code, 400)
        self.assertIn("basis is not finite", response.get_json()["error"])


if __name__ == "__main__":
    unittest.main()
