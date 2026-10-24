# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the Python ⟷ browser golden of the G(r) error ribbon.

Writes ``web_app/frontend/src/__tests__/fixtures/ft_uncertainty_fixture.json``:
run folders and the payload ``rmc_toolkits.ft_uncertainty.ft_uncertainty_for_run``
makes of each one's ``<stem>_FT_XFQ1.csv`` -- the ``uncertainty`` that Flask's
``/api/plot/data`` adds to an xPDF plot. The browser's port
(``src/ftUncertainty.js``, through ``src/ftUncertaintyRun.js`` on the same
files) is held to it by ``src/__tests__/ftUncertaintyParity.test.js``, and
``tests/test_ft_uncertainty_fixture.py`` fails when the committed file no
longer matches the Python engine.

Cases:

* the committed demo run (``web_app/frontend/public/demo``), read in place:
  its input file ``scale_ft_rmc.fq`` is not in the folder, so no ribbon, but
  the transform reproduces its G(r);
* synthetic runs built the way RMCProfile writes them: an input F(Q) (a toy
  S(Q) − 1 with σ = 0.002 (1 + 0.15 Q)), the run's F(Q)_Expt = s · K F − c
  (K the box convolution of radius R when the .dat says CONVOLVE), and the
  G(r) column Σ ΔQ Q (F_Expt + c) sin(Q r): with and without CONVOLVE, the
  7-decimal ``_FQ1.csv`` alone, no error column, a zero σ, each file missing
  in turn, a G(r) made with the trapezoid rule, an input that does not give
  the run's F(Q), and an r range so close to R that the box convolution
  would change σ(r) beyond the limit.

The files' text is committed in the fixture, so the staleness test checks the
engine, not this construction. Regenerate from the checkout under test:

    PYTHONPATH=. python tests/generate_ft_uncertainty_fixture.py
"""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import numpy as np

from rmc_toolkits.ft_uncertainty import box_convolution, ft_uncertainty_for_run, rmc_ft_coefficients

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "ft_uncertainty_fixture.json"
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"
DEMO_FILES = ["GTS_250K.dat", "GTS_250K_FT_XFQ1.csv", "GTS_250K_XFQ1.csv", "GTS_250K_FQ1.csv", "GTS_250K_PDFpartials.csv"]

STEM = "rmc_run"
INPUT = "scale_ft_rmc.fq"
SHELLS = ((2.5, 0.6, 0.006), (3.6, 0.45, 0.009), (4.4, 0.35, 0.011))


def toy_f(q: np.ndarray) -> np.ndarray:
    """S(Q) − 1 of damped Debye shells (MATERIA's toySofQ − 1)."""
    s = np.ones_like(q)
    for d, n, u in SHELLS:
        s = s + (n * np.sin(q * d)) / (q * d) * np.exp(-0.5 * u * q * q * d)
    return s - 1.0


def _row(values, sep: str = " , ") -> str:
    return "   " + sep.join(repr(float(value)) for value in values) + "   "


def _csv(header: str, columns, fmt=None) -> str:
    rows = []
    for values in zip(*columns):
        if fmt is None:
            rows.append(_row(values))
        else:
            rows.append(" , ".join(fmt.format(float(value)).rjust(15) for value in values))
    return "\n".join([header, *rows]) + "\n"


def _dat(*, file: str = INPUT, convolve: bool = True) -> str:
    lines = [
        "TITLE :: synthetic run",
        "",
        "NUMBER_DENSITY :: 0.046690 Angstrom^(-3)",
        "R_SPACING :: 0.0400 Angstrom",
        "",
        "FLAGS ::",
        "  > NO_MOVEOUT",
        "",
        "XRAY_RECIPROCAL_SPACE_DATA ::",
        f"  > FILENAME :: {file}",
        "  > DATA_TYPE :: F(Q)",
        "  > FIT_TYPE ::  F(Q)",
        "  > RECIPROCAL_SPACE_PARAMETERS :: 1 300 1",
        "  > REAL_SPACE_PARAMETERS :: 38 250 1",
        "  > FITTED_OFFSET",
        "  > FITTED_SCALE",
    ]
    if convolve:
        lines.append("  > CONVOLVE ::")
    lines += ["", "END ::", ""]
    return "\n".join(lines)


def synthetic_run(
    *,
    convolve: bool = True,
    radius: float = 25.0,
    build_radius: float | None = None,
    sigma: str = "clean",
    rule: str = "rmc",
    r_max: float = 10.0,
    scale: float = 1.15,
    offset: float = 2.7e-4,
    fq_name: str = f"{STEM}_XFQ1.csv",
    fq_fmt: str | None = None,
) -> dict:
    """A run folder's files as {name: text} (every file; cases leave some out)."""
    q = 0.8 + 0.05 * np.arange(300)
    f_input = toy_f(q)
    base = box_convolution(q, f_input, build_radius or radius) if convolve else f_input
    f_run = scale * base - offset
    r = 1.5 + 0.04 * np.arange(int(round((r_max - 1.5) / 0.04)) + 1)
    weights = rmc_ft_coefficients(q)
    if rule == "trapezoid":
        weights = weights.copy()
        weights[[0, -1]] *= 0.5
    g = np.sin(np.outer(r, q)) @ (weights * (f_run + offset))
    sig = 0.002 * (1.0 + 0.15 * q)
    if sigma == "zero":
        sig = sig.copy()
        sig[100] = 0.0
    input_columns = [q, f_input] if sigma == "none" else [q, f_input, sig]
    input_rows = ["  ".join(repr(float(value)) for value in values) for values in zip(*input_columns)]
    partials_r = radius * np.array([0.25, 0.5, 0.75, 1.0])
    return {
        f"{STEM}.dat": _dat(convolve=convolve),
        INPUT: "\n".join([f"{q.size:>12d}", "rmc S(Q)#L", *input_rows]) + "\n",
        fq_name: _csv("Q , F(Q)_RMC, F(Q)_Expt", [q, 0.97 * f_run, f_run], fq_fmt),
        f"{STEM}_FT_XFQ1.csv": _csv("r(A) ,   X_ray-calc  ,   X_ray_exp_renorm", [r, 0.97 * g, g]),
        f"{STEM}_PDFpartials.csv": "\n".join(
            ["   r (Ang),     A-A,     A-B"]
            + [f"   {float(value)!r} ,   1.0 ,   1.0 , " for value in partials_r]
        ) + "\n",
    }


def _without(files: dict, *names: str) -> dict:
    return {name: text for name, text in files.items() if name not in names}


def synthetic_cases() -> list[dict]:
    full = synthetic_run()
    return [
        {"name": "propagated, CONVOLVE", "files": full},
        {"name": "propagated, no CONVOLVE", "files": synthetic_run(convolve=False, scale=0.9, offset=-1e-3)},
        {"name": "the 7-decimal _FQ1.csv alone", "files": synthetic_run(fq_name=f"{STEM}_FQ1.csv", fq_fmt="{:.7f}")},
        {"name": "no error column", "files": synthetic_run(sigma="none")},
        {"name": "a zero sigma", "files": synthetic_run(sigma="zero")},
        {"name": "input file missing", "files": _without(full, INPUT)},
        {"name": "run-control file missing", "files": _without(full, f"{STEM}.dat")},
        {"name": "F(Q) file missing", "files": _without(full, f"{STEM}_XFQ1.csv")},
        {"name": "PDFpartials missing", "files": _without(full, f"{STEM}_PDFpartials.csv")},
        {"name": "G(r) by the trapezoid rule", "files": synthetic_run(rule="trapezoid")},
        {"name": "input built with another radius", "files": synthetic_run(build_radius=20.0)},
        {"name": "r range at the box radius", "files": synthetic_run(radius=10.15)},
    ]


def payload_for(files: dict) -> dict:
    with tempfile.TemporaryDirectory() as folder:
        for name, text in files.items():
            (Path(folder) / name).write_text(text, encoding="utf-8")
        return ft_uncertainty_for_run(Path(folder) / f"{STEM}_FT_XFQ1.csv")


def build_fixture(cases: list[dict] | None = None) -> dict:
    """The golden: the demo case and the synthetic cases, each with its payload."""
    if cases is None:
        cases = synthetic_cases()
    texts: dict[str, str] = {}
    keyed = []
    for case in cases:
        refs = {}
        for name, text in case["files"].items():
            key = next((k for k, v in texts.items() if v == text), None)
            if key is None:
                key = f"t{len(texts)}"
                texts[key] = text
            refs[name] = key
        keyed.append({"name": case["name"], "files": refs, "expected": payload_for(case["files"])})
    demo = ft_uncertainty_for_run(DEMO / "GTS_250K_FT_XFQ1.csv")
    return {
        "demo": {"folder": "public/demo", "files": DEMO_FILES, "transform": "GTS_250K_FT_XFQ1.csv", "expected": demo},
        "stem": STEM,
        "texts": texts,
        "cases": keyed,
    }


def committed_cases(fixture: dict) -> list[dict]:
    """The synthetic cases back from a fixture's texts (the staleness test re-runs these)."""
    return [
        {"name": case["name"], "files": {name: fixture["texts"][key] for name, key in case["files"].items()}}
        for case in fixture["cases"]
    ]


def main() -> None:
    fixture = build_fixture()
    OUT.write_text(json.dumps(fixture, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    statuses = [(case["name"], case["expected"]["status"]) for case in fixture["cases"]]
    print(f"wrote {OUT}")
    for name, status in statuses:
        print(f"  {status:12s} {name}")
    print(f"  {fixture['demo']['expected']['status']:12s} demo")
    if not math.isfinite(fixture["demo"]["expected"]["reproduction"]["gr"]["residual"]):
        raise SystemExit("the demo's reproduction residual is not finite")


if __name__ == "__main__":
    main()
