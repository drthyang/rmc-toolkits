# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

from pathlib import Path
import unittest

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib

import rmc_toolkits


ROOT = Path(__file__).resolve().parents[1]


class PackageApiTests(unittest.TestCase):
    def test_top_level_exports_current_package_helpers(self):
        expected = {
            "detect_plot_kind",
            "kde_slice",
            "load_unit_cell_positions",
            "read_exafs_csv",
            "read_rmc_csv",
            "read_structure",
            "write_frac_from_rmc6f",
        }

        self.assertTrue(expected.issubset(set(rmc_toolkits.__all__)))
        for name in expected:
            self.assertTrue(hasattr(rmc_toolkits, name), name)

    def test_every_exported_name_resolves_once(self):
        self.assertEqual(len(rmc_toolkits.__all__), len(set(rmc_toolkits.__all__)))
        for name in rmc_toolkits.__all__:
            self.assertTrue(hasattr(rmc_toolkits, name), name)

    def test_top_level_exports_the_1_0_engine_api(self):
        # Public functions and constants the CLIs, the API and the docs use,
        # re-exported beside the rest of their modules' API.
        from rmc_toolkits import scaling, transforms, triplets

        expected = {
            scaling: (
                "auto_enforcement_cutoff", "first_shell_foot", "fz_limit_fit", "alias_limit",
                "RHO0_PHYSICAL_RANGE", "RHO0_SEED", "FZ_REL_SE_MAX",
            ),
            transforms: ("gpdf_slope_at_zero",),
            triplets: ("bond_angle_summary_from_file",),
        }
        for module, names in expected.items():
            for name in names:
                with self.subTest(name=name):
                    self.assertIn(name, rmc_toolkits.__all__)
                    self.assertIs(getattr(rmc_toolkits, name), getattr(module, name))

    def test_package_version_matches_pyproject(self):
        metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))

        self.assertEqual(rmc_toolkits.__version__, metadata["project"]["version"])


if __name__ == "__main__":
    unittest.main()
