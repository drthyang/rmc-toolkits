# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

import unittest

from rmc_toolkits.scattering import (
    COHERENT_B_FM,
    FaberZiman,
    faber_ziman,
    parse_formula,
)


class TableTests(unittest.TestCase):
    def test_known_scattering_lengths(self):
        # Spot checks against Sears (1992) / NIST NCNR values (fm).
        self.assertAlmostEqual(COHERENT_B_FM["O"], 5.803, places=3)
        self.assertAlmostEqual(COHERENT_B_FM["Si"], 4.1491, places=4)
        self.assertAlmostEqual(COHERENT_B_FM["H"], -3.739, places=3)
        self.assertAlmostEqual(COHERENT_B_FM["Ti"], -3.438, places=3)
        self.assertAlmostEqual(COHERENT_B_FM["Mn"], -3.73, places=2)
        self.assertAlmostEqual(COHERENT_B_FM["Nb"], 7.054, places=3)

    def test_table_covers_common_elements(self):
        for el in ("C", "N", "F", "Na", "Al", "S", "Cl", "Fe", "Cu", "Ba", "Pb", "U"):
            self.assertIn(el, COHERENT_B_FM)


class ParseFormulaTests(unittest.TestCase):
    def test_simple(self):
        self.assertEqual(parse_formula("GaTa4Se8"), {"Ga": 1.0, "Ta": 4.0, "Se": 8.0})

    def test_decimal_stoichiometry(self):
        counts = parse_formula("Sr0.5Ba0.5TiO3")
        self.assertEqual(counts, {"Sr": 0.5, "Ba": 0.5, "Ti": 1.0, "O": 3.0})

    def test_parentheses(self):
        self.assertEqual(parse_formula("Ca(OH)2"), {"Ca": 1.0, "O": 2.0, "H": 2.0})
        self.assertEqual(
            parse_formula("Al2(SO4)3"), {"Al": 2.0, "S": 3.0, "O": 12.0}
        )

    def test_repeated_element_accumulates(self):
        self.assertEqual(parse_formula("CH3COOH"), {"C": 2.0, "H": 4.0, "O": 2.0})

    def test_errors(self):
        for bad in ("", "(", "Ca(OH", "Ca)2", "123", "xyz"):
            with self.assertRaises(ValueError):
                parse_formula(bad)


class FaberZimanTests(unittest.TestCase):
    def test_argon_matches_pystog_example(self):
        # pystog's bundled argon config: "<b_coh>^2": 3.644, "<b_tot^2>": 5.435
        # — quoted in fm^2 (b_Ar = 1.909 fm).
        fz = faber_ziman("Ar")
        self.assertAlmostEqual(fz.b_avg_sq_fm2, 3.644, places=2)
        self.assertAlmostEqual(fz.b_avg_sq_barn, 0.03644, places=4)

    def test_monatomic_identity(self):
        # For one element <b^2> = <b>^2 exactly (Laue term vanishes).
        fz = faber_ziman("Ni")
        self.assertAlmostEqual(fz.b_avg_sq_fm2, fz.b_sq_avg_fm2, places=10)

    def test_polyatomic_laue_inequality_and_weights(self):
        fz = faber_ziman("GaTa4Se8")
        # <b^2> >= <b>^2 always (Cauchy-Schwarz); strict for unequal b.
        self.assertGreater(fz.b_sq_avg_fm2, fz.b_avg_sq_fm2)
        self.assertAlmostEqual(sum(fz.weights.values()), 1.0, places=10)
        self.assertIsInstance(fz, FaberZiman)

    def test_dict_composition_and_overrides(self):
        by_formula = faber_ziman("SiO2")
        by_dict = faber_ziman({"Si": 1, "O": 2})
        self.assertAlmostEqual(by_formula.b_avg_sq_fm2, by_dict.b_avg_sq_fm2, places=12)
        # Isotopic override: null-scattering "Si" makes only O contribute.
        fz = faber_ziman({"Si": 1, "O": 2}, b_overrides_fm={"Si": 0.0})
        expected_avg = (2.0 / 3.0) * COHERENT_B_FM["O"]
        self.assertAlmostEqual(fz.b_avg_sq_fm2, expected_avg**2, places=10)

    def test_unknown_element_requires_override(self):
        with self.assertRaises(ValueError):
            faber_ziman("Po")
        fz = faber_ziman("Po", b_overrides_fm={"Po": 5.0})
        self.assertAlmostEqual(fz.b_avg_sq_fm2, 25.0, places=10)

    def test_null_matrix_rejected(self):
        # Ti/Zr null matrix: <b> ~ 0 makes Faber-Ziman weights undefined.
        with self.assertRaises(ValueError):
            faber_ziman({"Ti": 0.676, "Zr": 0.324})


if __name__ == "__main__":
    unittest.main()
