# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the Python ⟷ browser golden of the termination-ripple assessment.

Writes ``web_app/frontend/src/__tests__/fixtures/ripples_fixture.json``: run
folders, the calls made on them, and what ``rmc_toolkits.ripples.ripples_for_run``
returns for each — what Flask's ``/api/ripples`` sends. The browser's port
(``src/ripples.js`` through ``src/ripplesRun.js`` on the same files) is held to
it by ``src/__tests__/ripplesParity.test.js``, and
``tests/test_ripples_fixture.py`` fails when the committed file no longer
matches the Python engine.

Cases: the committed Demo (an X-ray F(Q) fitted in real space), read in
place, with an r range, named features, a split and a pair; and synthetic run
folders built the way RMCProfile and stog write them — a neutron G(r) from a
stog input (Lorch, its Fourier filter), from Auto StoG (its provenance in
autoscale/, the low r set to the baseline), with CONVOLVE and nothing else, with
no word on Q_max (and then with the caller's), an X-ray G(r) with its
QMIN/QMAX, a partial beside two totals of different Q_max (features read
across), and the failures (an unreadable file, no F(Q) beside a transform, a
transform the F(Q) does not reproduce).

The files' text is committed in the fixture, so the staleness test checks the
engine, not this construction. Regenerate from the checkout under test:

    PYTHONPATH=. python tests/generate_ripples_fixture.py
"""

from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import numpy as np

from rmc_toolkits import ripples as R

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "ripples_fixture.json"
DEMO = ROOT / "web_app" / "frontend" / "public" / "demo"
DEMO_FILES = ["GTS_250K.dat", "GTS_250K_FT_XFQ1.csv", "GTS_250K_XFQ1.csv", "GTS_250K_FQ1.csv", "GTS_250K_PDFpartials.csv"]
DEMO_FITS = ["GTS_250K_FQ1.csv", "GTS_250K_FT_XFQ1.csv"]

STEM = "rmc_run"
DQ = 0.01
R_GRID = np.round(0.01 * np.arange(1, 701), 10)
PEAKS = (
    {"amplitude": 2.4, "r": 2.5, "sigma": 0.06},
    {"amplitude": 1.6, "r": 3.6, "sigma": 0.08},
    {"amplitude": 2.8, "r": 4.4, "sigma": 0.1},
    {"amplitude": 0.25, "r": 5.3, "sigma": 0.07},
)


def d_space(qmax: float, *, rule: str = "trapezoid", window: str = "none", q_min: float = DQ) -> np.ndarray:
    """The peaks' D-space function through a transform cut at Q_max."""
    n = int(round((qmax - q_min) / DQ)) + 1
    q = R.q_grid(q_min, DQ, n)
    weights = R.transform_weights(q, DQ, rule, window)
    return R.apply_transform(q, weights, sum(R.gaussian_sine(q, peak) for peak in PEAKS), R_GRID)


def g_of(d: np.ndarray, baseline: float) -> np.ndarray:
    return d / R_GRID + baseline


def gr_text(title: str, y: np.ndarray) -> str:
    return f"{R_GRID.size:>10d}\n{title}\n" + "".join(f"  {x:.4f}  {v:.10f}\n" for x, v in zip(R_GRID, y))


def pdf_csv(y: np.ndarray, calc: np.ndarray | None = None, start: float = 1.5) -> str:
    keep = R_GRID >= start - 1e-9
    calc = y if calc is None else calc
    rows = [f"   {x:.4f} ,  {c:.10e} ,  {v:.10e}" for x, c, v in zip(R_GRID[keep], calc[keep], y[keep])]
    return "\n".join(["r (Ang), PDF (RMC), PDF (Expt)", *rows]) + "\n"


def partials_csv(first: float = 2.2) -> str:
    r = np.round(0.01 * np.arange(1, 801), 10)
    rows = [f"   {x:.5f} ,   {0.0 if x < first else 1.0:.1f} ,   {0.0 if x < first + 0.3 else 1.0:.1f} , " for x in r]
    return "\n".join(["   r (Ang),     A-A,     A-B", *rows]) + "\n"


def dat(blocks: list[str], minimum: str = "2.2 2.2") -> str:
    return "\n".join([
        "TITLE :: synthetic run", "", f"MINIMUM_DISTANCES :: {minimum} Angstrom", "R_SPACING :: 0.0100 Angstrom", "",
        "ATOMS :: A B", "", *blocks, "END ::", "",
    ])


def neutron_block(file: str, extra: list[str] = (), coefficients: str = "0.03 0.02 0.01") -> str:
    return "\n".join([
        "NEUTRON_REAL_SPACE_DATA :: 1", f"  > FILENAME :: {file}", "  > DATA_TYPE :: G(r)", "  > FIT_TYPE :: G(r)",
        "  > START_POINT :: 150", "  > END_POINT :: 700", f"  > NEUTRON_COEFFICIENTS :: {coefficients}", *extra, "",
    ])


def stog_inp(out_rmc_gr: str, *, qmax: float, lorch: bool, r_cutoff: float) -> str:
    lines = ["1", "data.sq", f"0.01 {qmax}", "0 1", "0", "scale.fq", "scale.gr", "7", "700", "Y" if lorch else "N",
             "0.09", "0", "N", "Y", f"{r_cutoff}", "scale_ft.sq", "scale_ft.gr", "1.0", "scale_ft_rmc.fq", out_rmc_gr,
             "scale_ft_rmc.dr", "1.0 0 0"]
    return "\n".join(lines) + "\n"


def xray_ft_files(qmax: float, *, reproduce: bool = True, with_fq: bool = True) -> dict:
    """An X-ray F(Q) fitted in real space: <stem>_FT_XFQ1.csv = Σ ΔQ Q F sin(Qr) of <stem>_XFQ1.csv."""
    dq = 0.02
    q = R.q_grid(0.5, dq, int(round((qmax - 0.5) / dq)) + 1)
    weights = R.transform_weights(q, dq, "rectangle", "none")
    # Q F(Q) = (2/π)·S(Q) of the peaks: then Σ ΔQ Q F sin(Qr) = (2/π) Σ ΔQ S sin(Qr).
    spectrum = sum(R.gaussian_sine(q, peak) for peak in PEAKS)
    f = np.array([float(f"{v:.12g}") for v in (2.0 / math.pi) * spectrum / q])
    r = np.round(1.5 + 0.01 * np.arange(551), 10)
    g = np.sin(np.outer(r, q)) @ (dq * q * f)
    if not reproduce:
        g = g * (1.0 + 0.01 * np.sin(3.0 * r))
    files = {
        f"{STEM}_FT_XFQ1.csv": "\n".join(["r(A) ,   X_ray-calc  ,   X_ray_exp_renorm",
                                          *(f"   {x:.4f} , {0.98 * v:.12g} , {v:.12g}" for x, v in zip(map(float, r), map(float, g)))]) + "\n",
    }
    if with_fq:
        files[f"{STEM}_XFQ1.csv"] = "\n".join(["Q , F(Q)_RMC, F(Q)_Expt",
                                              *(f"  {x:.4f} , {0.98 * v:.12g} , {v:.12g}" for x, v in zip(map(float, q), map(float, f)))]) + "\n"
    return files


def synthetic_cases() -> list[dict]:
    lorch_d = d_space(25.0, window="lorch")
    rect_d = d_space(25.0)
    low_d = d_space(18.0)
    g_lorch = g_of(lorch_d, -0.06)
    g_rect = g_of(rect_d, -0.06)
    enforced = g_rect.copy()
    enforced[R_GRID <= 2.0] = -0.06
    provenance = {
        "tool": "rmc-autoscale", "outputs": {"rmc_gr": "/data/autoscale/scale_ft_rmc.gr"},
        "enforcement": {"cutoff": 2.0, "peak_rmin": 2.0, "peak_rmax": 2.0, "source": "auto (first-shell foot)"},
        "provenance": {"n_q_points": 2500, "config": {"qmin": 0.01, "qmax": 25.0, "lorch": False, "r_cutoff": 0.8}},
    }
    xray_real = R.apply_transform(R.q_grid(DQ, DQ, 2000), R.transform_weights(R.q_grid(DQ, DQ, 2000), DQ, "rectangle", "none"),
                                  sum(R.gaussian_sine(R.q_grid(DQ, DQ, 2000), p) for p in PEAKS), R_GRID)
    partial = g_of(rect_d * 0.4, -0.02)
    return [
        {
            "name": "a neutron G(r) from a stog input (Lorch, its filter)",
            "files": {
                f"{STEM}.dat": dat([neutron_block("scale_ft_rmc.gr")]),
                "scale_ft_rmc.gr": gr_text("stog G(r)", g_lorch),
                "stog.inp": stog_inp("scale_ft_rmc.gr", qmax=25.0, lorch=True, r_cutoff=1.0),
                f"{STEM}_PDF1.csv": pdf_csv(g_lorch),
                f"{STEM}_PDFpartials.csv": partials_csv(),
            },
            "fits": [f"{STEM}_PDF1.csv"],
            "calls": [{"rRange": [2.0, 5.0]}, {"features": [2.26, 5.3], "split": [3.55, 3.66]}],
        },
        {
            "name": "a neutron G(r) from Auto StoG (provenance in autoscale/, enforced below 2 Å)",
            "files": {
                f"{STEM}.dat": dat([neutron_block("scale_ft_rmc.gr", ["  > FITTED_SCALE"])]),
                "scale_ft_rmc.gr": gr_text("rmc-autoscale 0.8.0: a=1.1 b=-0.2", enforced),
                "autoscale/stog_input_provenance.json": json.dumps(provenance, indent=1) + "\n",
                f"{STEM}_PDF1.csv": pdf_csv(0.8 * enforced + 0.001),
                f"{STEM}_PDFpartials.csv": partials_csv(),
            },
            "fits": [f"{STEM}_PDF1.csv"],
            "calls": [{"features": [2.26, 2.74]}],
        },
        {
            "name": "a neutron G(r) with CONVOLVE and no word on its data's transform",
            "files": {
                f"{STEM}.dat": dat([neutron_block("missing.gr", ["  > CONVOLVE :: 25.0", "  > R_CUTOFF :: 1.9"])]),
                f"{STEM}_PDF1.csv": pdf_csv(g_rect, start=0.01),
                f"{STEM}_PDFpartials.csv": partials_csv(),
            },
            "fits": [f"{STEM}_PDF1.csv"],
            "calls": [{"rRange": [2.0, 4.0]}],
        },
        {
            "name": "a neutron G(r) with no word on Q_max",
            "files": {
                f"{STEM}.dat": dat([neutron_block("missing.gr", [], "0.06")], minimum="2.2"),
                f"{STEM}_PDF1.csv": pdf_csv(g_rect, start=0.01),
            },
            "fits": [f"{STEM}_PDF1.csv"],
            "calls": [{}, {"qmax": 25.0, "features": [2.26]}, {"qmax": 25.0, "qmin": 0.5, "window": "lorch", "features": [2.26]}],
        },
        {
            "name": "an X-ray G(r) with its QMIN-QMAX",
            "files": {
                f"{STEM}.dat": dat(["\n".join(["XRAY_REAL_SPACE_DATA ::", "  > FILENAME :: x.gr", "  > DATA_TYPE :: G(r)",
                                               "  > FIT_TYPE :: D(r)", "  > QMIN :: 0.0", "  > QMAX :: 20.0",
                                               "  > QSTEP :: 0.01", ""])]),
                f"{STEM}_Xray_PDF1.csv": pdf_csv(xray_real, start=0.01),
                f"{STEM}_PDFpartials.csv": partials_csv(),
            },
            "fits": [f"{STEM}_Xray_PDF1.csv"],
            "calls": [{"features": [2.25]}],
        },
        {
            "name": "two totals of different Q_max and a partial",
            "files": {
                f"{STEM}.dat": dat([
                    "\n".join(["XRAY_RECIPROCAL_SPACE_DATA ::", "  > FILENAME :: scale_ft_rmc.fq", "  > FIT_TYPE :: F(Q)",
                               "  > REAL_SPACE_FIT :: 1 5000 1", ""]),
                    neutron_block("low.gr"),
                    neutron_block("partial.gr", ["  > CONVOLVE :: 18.0"], "0.0 0.02 0.0").replace(":: 1", ":: 2", 1),
                ]),
                **xray_ft_files(25.0),
                "low.gr": gr_text("rmc G(r)", g_of(low_d, -0.06)),
                "partial.gr": gr_text("partial", partial),
                "stog.inp": stog_inp("low.gr", qmax=18.0, lorch=False, r_cutoff=0.0),
                f"{STEM}_PDF1.csv": pdf_csv(g_of(low_d, -0.06)),
                f"{STEM}_PDF2.csv": pdf_csv(partial),
                f"{STEM}_PDFpartials.csv": partials_csv(),
            },
            "fits": [f"{STEM}_FT_XFQ1.csv", f"{STEM}_PDF1.csv", f"{STEM}_PDF2.csv"],
            "calls": [{"features": [2.26, 2.72, 5.3]}, {"pair": "A-B", "features": [2.5, 2.76]}],
        },
        {
            "name": "the failures",
            "files": {
                f"{STEM}.dat": dat(["\n".join(["XRAY_RECIPROCAL_SPACE_DATA ::", "  > FILENAME :: scale_ft_rmc.fq",
                                               "  > REAL_SPACE_FIT :: 1 5000 1", ""])]),
                **xray_ft_files(25.0, reproduce=False),
                f"{STEM}_SQ1.csv": "Q , S(Q)_RMC, S(Q)_Expt\n 0.5 , 1 , x\n",
                f"other_FT_XFQ1.csv": xray_ft_files(25.0, with_fq=False)[f"{STEM}_FT_XFQ1.csv"],
            },
            "fits": [f"{STEM}_FT_XFQ1.csv", f"{STEM}_SQ1.csv", "other_FT_XFQ1.csv"],
            "calls": [{}],
        },
    ]


def call_kwargs(call: dict) -> dict:
    return {
        "features": call.get("features"), "r_range": call.get("rRange"), "split": call.get("split"),
        "pair": call.get("pair"), "qmax": call.get("qmax"), "qmin": call.get("qmin"), "window": call.get("window"),
    }


def result_for(files: dict, fits: list[str], call: dict) -> dict:
    with tempfile.TemporaryDirectory() as folder:
        for name, text in files.items():
            path = Path(folder) / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8")
        return R.ripples_for_run([Path(folder) / name for name in fits], **call_kwargs(call))


def _plain(value):
    """JSON-safe: tuples to lists, numpy scalars to floats."""
    return json.loads(json.dumps(value, allow_nan=False, default=float))


DEMO_CALLS = [
    {"rRange": [2.0, 6.0], "split": [2.45, 2.6]},
    {"features": [2.32, 2.8, 3.18]},
    {"pair": "Ta-Se", "features": [2.55, 2.75]},
]


def build_fixture(cases: list[dict] | None = None) -> dict:
    """The golden: the Demo's calls and the synthetic cases' calls, each with its result."""
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
        keyed.append({
            "name": case["name"], "files": refs, "fits": case["fits"],
            "calls": [{"args": call, "expected": _plain(result_for(case["files"], case["fits"], call))}
                      for call in case["calls"]],
        })
    demo = [{"args": call, "expected": _plain(R.ripples_for_run([DEMO / name for name in DEMO_FITS], **call_kwargs(call)))}
            for call in DEMO_CALLS]
    return {
        "demo": {"folder": "public/demo", "files": DEMO_FILES, "fits": DEMO_FITS, "calls": demo},
        "stem": STEM,
        "texts": texts,
        "cases": keyed,
    }


def committed_cases(fixture: dict) -> list[dict]:
    """The synthetic cases back from a fixture's texts (the staleness test re-runs these)."""
    return [
        {"name": case["name"], "fits": case["fits"], "calls": [call["args"] for call in case["calls"]],
         "files": {name: fixture["texts"][key] for name, key in case["files"].items()}}
        for case in fixture["cases"]
    ]


#: A reading this close to a rule's threshold (relative) is refused: the two
#: engines' fits agree to ~1e-9, and the golden must not depend on which side.
THRESHOLD_MARGIN = 1e-6


def near_thresholds(fixture: dict) -> list[str]:
    """Readings within THRESHOLD_MARGIN of a rule's threshold (the classes would be round-off)."""
    found = []
    calls = [("demo", call) for call in fixture["demo"]["calls"]]
    calls += [(case["name"], call) for case in fixture["cases"] for call in case["calls"]]
    for name, call in calls:
        for feature in call["expected"]["features"] or []:
            for reading in feature["datasets"]:
                ripple = reading["ripple"]
                if ripple and ripple["amplitude"]:
                    ratio = abs(reading["size"]) / abs(ripple["amplitude"])
                    for bound in (R.RIPPLE_SIZE_FACTOR, 1.0 / R.RIPPLE_SIZE_FACTOR, 1.0 / R.NEGLIGIBLE_RIPPLE_FRACTION):
                        if abs(ratio / bound - 1.0) < THRESHOLD_MARGIN:
                            found.append(f"{name}: feature at {reading['r']} has |size|/|ripple| at {bound}")
    return found


def main() -> None:
    fixture = build_fixture()
    close = near_thresholds(fixture)
    if close:
        raise SystemExit("refusing a golden on a threshold:\n  " + "\n  ".join(close))
    OUT.write_text(json.dumps(fixture, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")
    for case in [{"name": "demo", "calls": fixture["demo"]["calls"]}, *fixture["cases"]]:
        for call in case["calls"]:
            statuses = [(d["file"], d["status"], (d["fit"] or {}).get("risk")) for d in call["expected"]["datasets"]]
            print(f"  {case['name'][:60]:60s} {json.dumps(call['args'])[:50]:50s} {statuses}")


if __name__ == "__main__":
    main()
