# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Atomic Density (hkl) slices: Flask draws the map in the page's in-plane frame.

For a custom normal the Flask route used to pick its own in-plane axes
(``kde._orthogonal_axis``: argmin |n|) while the browser worker, the Slab In
Cell panel and the panel's aspect use StructurePage's ``makeFreePlaneBasis``
(reference a unless |n_a| >= 0.85). For the default (1 1 0) plane the Flask map
was therefore rotated 90 degrees against the browser and against the page's own
Slab In Cell panel. The page now sends its u/v with a custom slice and the
route passes them to ``oriented_kde_slice``; the densities are then the ones
the committed browser-parity golden pins (``kdeParity.test.js``).
"""

from __future__ import annotations

from pathlib import Path
import json
import os
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))
os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import app as backend_app  # noqa: E402
from generate_kde_fixture import OUT, browser_plane_basis  # noqa: E402
from rmc_toolkits.kde import load_unit_cell_positions, oriented_kde_slice  # noqa: E402

DEMO_DIR = "web_app/frontend/public/demo"
DEMO = ROOT / DEMO_DIR / "GTS_250K.rmc6f"
CUBE_CORNERS = np.array([[i, j, k] for k in (0, 1) for j in (0, 1) for i in (0, 1)], dtype=float)
SQRT_HALF = np.sqrt(0.5)

# (fixture case name, element, (h k l), z, dz, bw) -- the same slices the
# browser-parity golden holds, so route == golden == browser worker.
CASES = (
    ("demo Ta (110) bw=0.06", "Ta", (1, 1, 0), 0.37, 0.02, 0.06),
    ("demo Se (101) bw=0.05", "Se", (1, 0, 1), 0.37, 0.01, 0.05),
)


def frame_query(u, v) -> dict:
    return dict(zip(("ux", "uy", "uz", "vx", "vy", "vz"), (*u, *v)))


class CustomSliceFrameTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()

    def slice(self, element, hkl, z, dz, bw, **extra):
        normal = np.asarray(hkl, dtype=float) / np.linalg.norm(hkl)
        params = {
            "dir": DEMO_DIR, "element": element, "orientation": "custom",
            "nx": normal[0], "ny": normal[1], "nz": normal[2],
            "z": z, "dz": dz, "bw": bw, "grid": 32, "levels": 0, **extra,
        }
        return self.client.get("/api/kde/slice", query_string=params)

    def test_the_browser_frames_of_the_default_planes(self):
        # The literal frames the page (and the JS test customSliceFrame.test.js) use.
        _, u, v = browser_plane_basis([1, 1, 0])
        np.testing.assert_allclose(u, [SQRT_HALF, -SQRT_HALF, 0.0], atol=1e-15)
        np.testing.assert_allclose(v, [0.0, 0.0, -1.0], atol=1e-15)
        _, u, v = browser_plane_basis([1, 0, 1])
        np.testing.assert_allclose(u, [SQRT_HALF, 0.0, -SQRT_HALF], atol=1e-15)
        np.testing.assert_allclose(v, [0.0, 1.0, 0.0], atol=1e-15)

    def test_route_draws_in_the_frame_the_page_sends(self):
        fixture = {case["name"]: case for case in json.loads(OUT.read_text(encoding="utf-8"))["cases"]}
        for name, element, hkl, z, dz, bw in CASES:
            with self.subTest(name):
                normal, u, v = browser_plane_basis(list(hkl))
                response = self.slice(element, hkl, z, dz, bw, **frame_query(u, v))
                self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
                payload = response.get_json()
                np.testing.assert_allclose(payload["uVector"], u, atol=1e-12)
                np.testing.assert_allclose(payload["vVector"], v, atol=1e-12)

                # The extent the browser worker computes from the same frame.
                us, vs = CUBE_CORNERS @ np.asarray(u), CUBE_CORNERS @ np.asarray(v)
                np.testing.assert_allclose(
                    payload["extent"], [us.min(), us.max(), vs.min(), vs.max()], atol=1e-12
                )

                # Identical to the engine called with that frame ...
                positions = load_unit_cell_positions(DEMO, element=element).fractional_positions
                direct = oriented_kde_slice(
                    positions, center=z, thickness=dz, normal=normal,
                    u_axis=np.asarray(u), v_axis=np.asarray(v), bw=bw, grid=32, n_levels=0,
                )
                np.testing.assert_array_equal(np.asarray(payload["density"]), np.asarray(direct["density"]))

                # ... and to the golden the browser worker is held to.
                expected = fixture[name]["expected"]
                self.assertEqual(payload["slabCount"], expected["slabCount"])
                self.assertEqual(payload["fitCount"], expected["fitCount"])
                density = np.asarray(payload["density"]) / payload["vmax"]
                np.testing.assert_allclose(density, expected["densityOverPeak"], rtol=0, atol=1e-8)

    def test_a_frame_that_is_not_in_the_plane_is_a_400(self):
        normal, u, v = browser_plane_basis([1, 1, 0])
        bad = (
            (frame_query([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]), "u must be orthogonal to the slice normal"),
            (frame_query(u, [SQRT_HALF, SQRT_HALF, 0.0]), "v must be orthogonal to the slice normal"),
            (frame_query(u, [SQRT_HALF, -SQRT_HALF, 0.0]), "u and v must be orthogonal"),
            (frame_query([0.0, 0.0, 0.0], v), "u must be a non-zero 3D vector"),
            ({"ux": u[0], "uy": u[1]}, "ux/uy/uz/vx/vy/vz are required together"),
            ({**frame_query(u, v), "vz": "nan"}, "vz must be a finite number"),
        )
        for query, message in bad:
            with self.subTest(message):
                response = self.slice("Ga", (1, 1, 0), 0.5, 0.02, 0.03, **query)
                self.assertEqual(response.status_code, 400, response.get_data(as_text=True))
                self.assertIn(message, response.get_json()["error"])

    def test_presets_keep_their_fixed_frames(self):
        response = self.client.get(
            "/api/kde/slice",
            query_string={"dir": DEMO_DIR, "element": "Ga", "orientation": "c", "z": 0.25, "grid": 32, "levels": 0},
        )
        self.assertEqual(response.status_code, 200, response.get_data(as_text=True))
        payload = response.get_json()
        self.assertEqual(payload["uVector"], [1.0, 0.0, 0.0])
        self.assertEqual(payload["vVector"], [0.0, 1.0, 0.0])


if __name__ == "__main__":
    unittest.main()
