# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""estimate_rho0's ``extrapolated`` flag follows the data, not config.qmin (0.6.0 audit, stog-b).

The flag marks estimates whose Q->0 extrapolation is longer than the ~1 A^-1 head
it rests on. It was computed as ``config.qmin > 1.0``: a NaN-padded rebinned file
(or a stog.inp qmin below the first finite point) turned the warning off although
the estimate was identical — a 22 %-biased FeCoSn density was labelled a
measurement. It now uses the first cropped Q.
"""

from pathlib import Path
import sys
import unittest

import numpy as np

from rmc_toolkits.scaling import FZ_FIT_WIDTH, estimate_rho0

sys.path.insert(0, str(Path(__file__).resolve().parent))
from test_stog_b_rho0 import synthetic_config, synthetic_sq  # noqa: E402


class ExtrapolatedFlagTests(unittest.TestCase):
    def test_flag_uses_the_first_measured_q(self):
        q, sq, b_sq_avg = synthetic_sq()
        padded = sq.copy()
        padded[q < 1.3] = np.nan  # rebin-style NaN padding below 1.3 A^-1
        low = estimate_rho0(q, padded, synthetic_config(b_sq_avg=b_sq_avg, qmin=0.6))
        at_data = estimate_rho0(q, padded, synthetic_config(b_sq_avg=b_sq_avg, qmin=1.3))
        self.assertAlmostEqual(low["q_first"], 1.32, places=9)
        # The same estimate (qmin only moves the C1 tail window edge) ...
        self.assertLess(abs(low["rho0"] - at_data["rho0"]) / at_data["rho0"], 1e-6)
        self.assertTrue(low["extrapolated"])  # ... flagged the same way
        self.assertTrue(at_data["extrapolated"])

    def test_data_starting_below_the_fit_width_is_not_extrapolated(self):
        q, sq, b_sq_avg = synthetic_sq()
        est = estimate_rho0(q, sq, synthetic_config(b_sq_avg=b_sq_avg, qmin=0.3))
        self.assertAlmostEqual(est["q_first"], 0.6, places=9)
        self.assertLess(est["q_first"], FZ_FIT_WIDTH)
        self.assertFalse(est["extrapolated"])


if __name__ == "__main__":
    unittest.main()
