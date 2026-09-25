# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Request-parameter validation at the Flask boundary.

Every numeric query/JSON parameter must be a finite number inside its
documented range; anything else is an HTTP 400 with a message naming the
parameter -- never a 200 whose body carries bare NaN/Infinity tokens (invalid
JSON for browsers), never a silently empty map, and never a 500.

Runs on synthetic fixtures under results/ (the backend's data root is the repo
root, as in tests/test_backend_api.py), so it needs no sample data.
"""

from pathlib import Path
import json
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "web_app" / "backend") not in sys.path:
    sys.path.insert(0, str(ROOT / "web_app" / "backend"))

os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "rmc_toolkits_cache"))
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["XDG_CACHE_HOME"]).mkdir(parents=True, exist_ok=True)

import numpy as np  # noqa: E402
from flask.json.provider import DefaultJSONProvider  # noqa: E402

import app as backend_app  # noqa: E402


def strict_json(body: str):
    """Parse like a browser's JSON.parse: NaN/Infinity tokens are an error."""

    def reject(token):
        raise ValueError(f"non-standard JSON token {token}")

    return json.loads(body, parse_constant=reject)


def write_synthetic_rmc6f(path: Path, *, supercell=(6, 6, 6), seed=3) -> int:
    """Two-element cubic run (Nb at the origin, Se at the body centre), 8 A cell."""
    rng = np.random.default_rng(seed)
    n1, n2, n3 = supercell
    lines = [
        f"Supercell dimensions {n1} {n2} {n3}",
        "Lattice vectors (Ang):",
        f"{8.0 * n1} 0.0 0.0",
        f"0.0 {8.0 * n2} 0.0",
        f"0.0 0.0 {8.0 * n3}",
        "Atoms:",
    ]
    atom = 0
    for ix in range(n1):
        for iy in range(n2):
            for iz in range(n3):
                cell = np.array([ix, iy, iz], dtype=float)
                for reference, element, basis, sigma in (
                    (1, "Nb", (0.0, 0.0, 0.0), (0.004, 0.002, 0.002)),
                    (2, "Se", (0.5, 0.5, 0.5), (0.003, 0.003, 0.003)),
                ):
                    atom += 1
                    coord = (cell + np.asarray(basis)) / np.asarray(supercell)
                    coord = (coord + rng.normal(size=3) * np.asarray(sigma)) % 1.0
                    lines.append(
                        f"{atom} {element} [{reference}] {coord[0]:.10f} {coord[1]:.10f} "
                        f"{coord[2]:.10f} {reference} {ix} {iy} {iz}"
                    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return atom


class _ValidationCase(unittest.TestCase):
    RUN = "results/backend_validation_test"

    @classmethod
    def setUpClass(cls):
        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        cls.run_dir = ROOT / cls.RUN
        cls.run_dir.mkdir(parents=True, exist_ok=True)
        write_synthetic_rmc6f(cls.run_dir / "synthetic.rmc6f")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.run_dir, ignore_errors=True)

    def get(self, url, **params):
        return self.client.get(url, query_string={"dir": self.RUN, **params})

    def assertOk(self, response):
        body = response.get_data(as_text=True)
        self.assertEqual(response.status_code, 200, body[:300])
        return strict_json(body)

    def assertBadRequest(self, response, *fragments):
        body = response.get_data(as_text=True)
        self.assertEqual(response.status_code, 400, body[:300])
        message = strict_json(body)["error"]
        for fragment in fragments:
            self.assertIn(fragment, message)
        return message


class KdeSliceValidationTests(_ValidationCase):
    BASE = {"element": "Nb", "z": 0.0, "dz": 0.1, "grid": 32, "levels": 2}

    def slice(self, **overrides):
        return self.get("/api/kde/slice", **{**self.BASE, **overrides})

    def test_baseline_is_valid_json(self):
        payload = self.assertOk(self.slice())
        self.assertGreater(payload["slabCount"], 0)
        self.assertGreater(payload["vmax"], 0)

    def test_non_finite_values_are_rejected(self):
        for key in ("bw", "dz", "z", "nx", "ny", "nz"):
            for raw in ("nan", "inf", "-inf"):
                with self.subTest(key=key, raw=raw):
                    self.assertBadRequest(self.slice(orientation="custom", **{key: raw}), key)

    def test_bandwidth_must_be_positive(self):
        # bw=0 used to answer 200 with an all-zero map while reporting the slab
        # count; bw<0 was silently treated as |bw| and echoed back negative.
        for raw in ("0", "-0.03"):
            with self.subTest(bw=raw):
                self.assertBadRequest(self.slice(bw=raw), "bw")

    def test_thickness_must_lie_in_unit_interval(self):
        for raw in ("0", "-0.08", "1.5"):
            with self.subTest(dz=raw):
                self.assertBadRequest(self.slice(dz=raw), "dz")

    def test_integer_parameters_reject_text_and_fractions(self):
        for key, raw in (("grid", "abc"), ("grid", "nan"), ("grid", "32.5"), ("levels", "x")):
            with self.subTest(key=key, raw=raw):
                self.assertBadRequest(self.slice(**{key: raw}), key)

    def test_levels_must_be_in_range(self):
        for raw in ("-3", "1000"):
            with self.subTest(levels=raw):
                self.assertBadRequest(self.slice(levels=raw), "levels")

    def test_grid_is_clamped_like_the_engine(self):
        self.assertEqual(self.assertOk(self.slice(grid=4))["grid"], 16)

    def test_zero_custom_normal_is_bad_request(self):
        self.assertBadRequest(self.slice(orientation="custom", nx=0, ny=0, nz=0), "normal")

    def test_underflowing_bandwidth_is_a_flagged_zero_map_not_a_nan_map(self):
        # Finite and positive, so it passes the range check, but the kernel is
        # far below double precision (det H underflows, the normaliser would
        # overflow). The contract: never a silent empty or NaN map -- a finite
        # zero map that carries the `unresolved` warning (the page draws the
        # warning, not the map), not NaN, not a decline that blames SciPy.
        assert_flagged_zero_map(self, self.assertOk(self.slice(bw="1e-200")))

    def test_a_bandwidth_below_the_coordinate_resolution_is_flagged_too(self):
        # bw = 1e-30: det H is representable, but scipy's whitening of the
        # absolute coordinates cannot resolve the kernel (its self-check used
        # to fail and the slice declined with the SciPy "engine" message).
        assert_flagged_zero_map(self, self.assertOk(self.slice(bw="1e-30")))


def assert_flagged_zero_map(case, payload):
    """An unresolvable kernel's slice: finite zeros, the kernel, `unresolved`."""
    density = payload["density"]
    case.assertTrue(
        all(isinstance(value, float) and np.isfinite(value) for row in density for value in row),
        "every node must be a finite number (no NaN, no null)",
    )
    case.assertEqual(payload["vmax"], 0.0)
    case.assertIsNone(payload["message"])
    case.assertIsNotNone(payload["kernel"])
    case.assertEqual([warning["code"] for warning in payload["warnings"]], ["subgrid", "unresolved"])
    case.assertEqual(payload["contours"], [])


def _nulls(value):
    """``value`` with every non-finite float replaced by None."""
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, dict):
        return {key: _nulls(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_nulls(item) for item in value]
    return value


class _NullingProvider(DefaultJSONProvider):
    """Writes non-finite floats as ``null`` whatever the caller asks for.

    The shape of an app-wide policy for masked data series (a NaN region of an
    RMCProfile CSV must reach the chart as a gap): it forces ``allow_nan=False``
    and, when that raises, serializes a NaN-free copy -- so asking this provider
    for ``allow_nan=False`` never raises.
    """

    def dumps(self, obj, **kwargs):
        kwargs.setdefault("default", self.default)
        kwargs["allow_nan"] = False
        try:
            return json.dumps(obj, **kwargs)
        except ValueError:
            return json.dumps(_nulls(obj), **kwargs)


class _NanTokenProvider(DefaultJSONProvider):
    """Always emits bare NaN/Infinity tokens (ignores ``allow_nan=False``)."""

    def dumps(self, obj, **kwargs):
        kwargs.setdefault("default", self.default)
        kwargs["allow_nan"] = True
        return json.dumps(obj, **kwargs)


class NonFiniteResultUnderAnyProviderTests(_ValidationCase):
    """A NaN/Infinity *computed result* is a 400 whatever JSON provider is installed.

    The guard must not depend on ``app.json`` raising for ``allow_nan=False``:
    a provider that writes NaN as null would otherwise turn an underflowed KDE
    into a 200 with an all-null map -- the silently empty result the guard
    exists to prevent.
    """

    PROVIDERS = (_NullingProvider, _NanTokenProvider)

    def with_provider(self, provider_class):
        original = backend_app.app.json
        backend_app.app.json = provider_class(backend_app.app)
        self.addCleanup(setattr, backend_app.app, "json", original)

    def test_underflowing_kde_slice_is_a_flagged_finite_map(self):
        # A NaN map would come out as nulls under _NullingProvider and as bare
        # NaN tokens (rejected by strict_json) under _NanTokenProvider; the
        # engine returns finite zeros with the `unresolved` warning instead.
        for provider_class in self.PROVIDERS:
            with self.subTest(provider=provider_class.__name__):
                self.with_provider(provider_class)
                assert_flagged_zero_map(
                    self,
                    self.assertOk(self.get("/api/kde/slice", **{**KdeSliceValidationTests.BASE, "bw": "1e-200"})),
                )

    def test_a_nan_kde_slice_is_a_bad_request(self):
        # The guard itself: were the engine to return a NaN map, the route
        # answers 400 whatever the provider, never a 200 with nulls or tokens.
        real = backend_app.oriented_kde_slice

        def nan_slice(*args, **kwargs):
            result = real(*args, **kwargs)
            result["density"][0][0] = float("nan")
            return result

        for provider_class in self.PROVIDERS:
            with self.subTest(provider=provider_class.__name__):
                self.with_provider(provider_class)
                with mock.patch.object(backend_app, "oriented_kde_slice", nan_slice):
                    self.assertBadRequest(
                        self.get("/api/kde/slice", **{**KdeSliceValidationTests.BASE, "bw": "0.07"}),
                        "NaN or Infinity",
                    )

    def test_extreme_pca_kde_is_a_bad_request(self):
        base = {"referenceNumber": 1, "grid": 12, "projections": "false"}
        for provider_class in self.PROVIDERS:
            for key, raw in (("bw", "1e-300"), ("extent", "1e300")):
                with self.subTest(provider=provider_class.__name__, key=key):
                    self.with_provider(provider_class)
                    self.assertBadRequest(
                        self.get("/api/pca/kde", **{**base, key: raw}), "NaN or Infinity"
                    )

    def test_finite_results_are_unaffected(self):
        for provider_class in self.PROVIDERS:
            with self.subTest(provider=provider_class.__name__):
                self.with_provider(provider_class)
                payload = self.assertOk(self.get("/api/kde/slice", **KdeSliceValidationTests.BASE))
                self.assertGreater(payload["vmax"], 0)
                self.assertOk(
                    self.get("/api/pca/orientation", referenceNumber=1, frequency=4, geometry="false")
                )


class StructureValidationTests(_ValidationCase):
    def test_max_points_must_be_an_integer(self):
        for raw in ("abc", "nan", "inf", "12.5"):
            with self.subTest(maxPoints=raw):
                self.assertBadRequest(self.get("/api/structure", maxPoints=raw), "maxPoints")

    def test_max_points_is_clamped(self):
        payload = self.assertOk(self.get("/api/structure", maxPoints=5))
        self.assertLessEqual(payload["sampledAtoms"], 100)


class PcaValidationTests(_ValidationCase):
    def test_sites_probability_errors_are_bad_requests(self):
        # Used to be a 500 here while /api/pca/kde answered 400 for the same value.
        for raw in ("nan", "1.5", "0", "abc"):
            with self.subTest(probability=raw):
                self.assertBadRequest(self.get("/api/pca/sites", probability=raw), "probability")

    def test_kde_baseline_is_valid_json(self):
        payload = self.assertOk(
            self.get("/api/pca/kde", referenceNumber=1, grid=12, projections="false")
        )
        self.assertTrue(np.isfinite(payload["mass"]))

    def test_kde_non_finite_and_non_positive_parameters_are_rejected(self):
        base = {"referenceNumber": 1, "grid": 12, "projections": "false"}
        cases = [
            ("extent", "nan"), ("extent", "inf"), ("extent", "0"), ("extent", "-1"),
            ("bwScale", "nan"), ("bwScale", "inf"), ("bwScale", "0"),
            ("bw", "nan"), ("bw", "inf"), ("bw", "0"), ("bw", "-0.5"), ("bw", "wide"),
            ("probability", "nan"), ("grid", "abc"), ("referenceNumber", "1.5"),
        ]
        for key, raw in cases:
            with self.subTest(key=key, raw=raw):
                self.assertBadRequest(self.get("/api/pca/kde", **{**base, key: raw}), key)

    def test_kde_extreme_finite_parameters_are_an_error_not_a_nan_volume(self):
        base = {"referenceNumber": 1, "grid": 12, "projections": "false"}
        for key, raw in (("bw", "1e-300"), ("bwScale", "1e300"), ("extent", "1e300")):
            with self.subTest(key=key, raw=raw):
                self.assertBadRequest(
                    self.get("/api/pca/kde", **{**base, key: raw}), "NaN or Infinity"
                )

    def test_kde_accepts_named_and_numeric_bandwidths(self):
        base = {"referenceNumber": 1, "grid": 12, "projections": "false"}
        for raw in ("scott", "Silverman", "0.4"):
            with self.subTest(bw=raw):
                self.assertOk(self.get("/api/pca/kde", **{**base, "bw": raw}))


class OrientationValidationTests(_ValidationCase):
    BASE = {"referenceNumber": 1, "frequency": 4, "geometry": "false"}

    def orientation(self, **overrides):
        return self.get("/api/pca/orientation", **{**self.BASE, **overrides})

    def test_baseline_is_valid_json(self):
        self.assertOk(self.orientation())

    def test_bad_values_are_rejected_with_the_parameter_name(self):
        cases = [
            ("smoothing", "-5"), ("smoothing", "1000"), ("smoothing", "1.5"),
            ("minAmplitude", "nan"), ("minAmplitude", "-0.1"),
            ("minAmplitudeQuantile", "nan"), ("frequency", "2.5"),
        ]
        for key, raw in cases:
            with self.subTest(key=key, raw=raw):
                self.assertBadRequest(self.orientation(**{key: raw}), key)


class TripletsValidationTests(_ValidationCase):
    BASE = {"end1": "Se", "apex": "Nb", "end2": "Se", "r12Min": 6.0, "r12Max": 7.5}

    def triplets(self, **overrides):
        return self.get("/api/triplets", **{**self.BASE, **overrides})

    def test_baseline_is_valid_json(self):
        self.assertGreater(self.assertOk(self.triplets())["angleCount"], 0)

    def test_non_finite_window_and_bin_width_are_rejected_by_name(self):
        for key in ("r12Min", "r12Max", "binWidth"):
            for raw in ("nan", "inf", "abc"):
                with self.subTest(key=key, raw=raw):
                    self.assertBadRequest(self.triplets(**{key: raw}), key)


@unittest.skipUnless(
    hasattr(backend_app, "TRIPLETS_MAX_ANGLES"),
    "this backend has no /api/triplets work budget (TRIPLETS_MAX_ANGLES)",
)
class TripletsWorkBudgetTests(_ValidationCase):
    """Integration guard: once the backend has a triplets work budget, the route
    must forward it to the engine AND key its cache on it.

    /api/triplets calls the uncached bond_angle_summary_from_file under the file-signature cache
    with one ``params`` tuple that is both the cache key and the argument list;
    a budget left out of that tuple would silently fall back to the engine's
    unrestricted default (``max_angles=None``).
    """

    PARAMS = {"end1": "Se", "apex": "Nb", "end2": "Se", "r12Min": 6.0, "r12Max": 7.5}

    def test_budget_is_forwarded_to_the_engine_and_part_of_the_cache_key(self):
        backend_app._TRIPLETS_CACHE.clear()
        # Under the real budget: computed and cached.
        self.assertGreater(self.assertOk(self.get("/api/triplets", **self.PARAMS))["angleCount"], 1)
        original = backend_app.TRIPLETS_MAX_ANGLES
        backend_app.TRIPLETS_MAX_ANGLES = 1
        self.addCleanup(setattr, backend_app, "TRIPLETS_MAX_ANGLES", original)
        # Same request under a budget of one angle: refused, not served from
        # the entry computed under the larger budget.
        self.assertBadRequest(self.get("/api/triplets", **self.PARAMS), "limit")


class ScalingValidationTests(unittest.TestCase):
    RUN = "results/backend_validation_scaling"
    RHO0 = 0.05
    B2 = 0.02

    @classmethod
    def setUpClass(cls):
        from rmc_toolkits.parsers import write_stog_xy
        from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

        backend_app.app.config.update(TESTING=True)
        cls.client = backend_app.app.test_client()
        cls.run_dir = ROOT / cls.RUN
        cls.run_dir.mkdir(parents=True, exist_ok=True)
        q = np.arange(20, 981) * 0.03
        r = np.arange(1, 12001) * 0.005
        onset = 0.5 * (1.0 + np.tanh((r - 2.65) / 0.07))
        peak = 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
        gpdf = g_to_gpdf(r, onset + peak, cls.RHO0)
        sq_true = fq_to_sq(q, gpdf_to_fq(r, gpdf, q))
        write_stog_xy(cls.run_dir / "synth.dat", q, (sq_true + 9.0) / 10.0, title="synthetic")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.run_dir, ignore_errors=True)

    def manual(self, **overrides):
        body = {
            "path": f"{self.RUN}/synth.dat",
            "qmin": 0.6, "qmax": 30, "rho0": self.RHO0, "bAvgSq": self.B2, "r0": 2.65,
            "mode": "manual", "a": 10.0, "b": -9.0, "enforce": False,
        }
        body.update(overrides)
        return self.client.post(
            "/api/scaling/preview",
            data=json.dumps(body, allow_nan=True),
            content_type="application/json",
        )

    def assertBadRequest(self, response, fragment):
        body = response.get_data(as_text=True)
        self.assertEqual(response.status_code, 400, body[:300])
        self.assertIn(fragment, strict_json(body)["error"])

    def test_manual_baseline_is_valid_json(self):
        response = self.manual()
        self.assertEqual(response.status_code, 200)
        strict_json(response.get_data(as_text=True))

    def test_zero_manual_scale_is_rejected(self):
        # a=0 used to answer 200 with sqRaw=(S-b)/0 serialized as bare NaN.
        self.assertBadRequest(self.manual(a=0), "non-zero scale 'a'")

    def test_non_numeric_and_non_finite_payload_numbers_are_rejected(self):
        cases = [
            ("qmin", [1]), ("qmin", True), ("qmin", "abc"), ("qmax", float("inf")),
            ("rho0", float("nan")), ("b", "nan"), ("a", float("inf")), ("r0", {"x": 1}),
            ("nr", 2.5), ("rmax", "-inf"),
        ]
        for key, raw in cases:
            with self.subTest(key=key, raw=raw):
                self.assertBadRequest(self.manual(**{key: raw}), key)


if __name__ == "__main__":
    unittest.main()
