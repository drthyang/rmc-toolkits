import unittest

import rmc_toolkits


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


if __name__ == "__main__":
    unittest.main()
