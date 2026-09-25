# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Regenerate the Flask ⟷ browser parity golden for RMCProfile output files.

Writes ``web_app/frontend/src/__tests__/fixtures/plot_parity_fixture.json``:
synthetic files in the RMCProfile v6 output layouts that no real run on hand
covers — neutron ``*_PDFn.csv`` / ``*_SQn.csv``, ``*_bragg*.csv`` (time-of-flight
and Q), EXAFS ``_Q_OUTPUT`` / ``_R_OUTPUT`` — with NaN-masked regions, CRLF line
endings, trailing commas, a leading blank line and E-notation, each paired with
the payload the FLASK path produces for it (``/api/plot/data`` and
``/api/plot/metadata`` through the test client). ``plotParity.test.js`` runs
``browserData.plotDataFromText`` on the same text and must reproduce it;
``tests/test_parsers_plot_payload.py`` regenerates it in memory and fails if the
committed file is stale. Each case also records, by construction, which column
is the calculation and which the experiment (``truth``), so the tests can check
the R-factor's column roles independently of either implementation.

    python tests/generate_plot_parity_fixture.py
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "plot_parity_fixture.json"


def _fmt(value: float) -> str:
    return "NaN" if not math.isfinite(value) else f"{value:.7f}"


def _table(header: str, columns: list[np.ndarray], *, newline: str = "\n", trailing_comma: bool = False,
           leading_blank: bool = False, title: str | None = None, exponent: bool = False) -> str:
    lines = [""] if leading_blank else []
    if title is not None:
        lines.append(title)
    lines.append(header)
    for row in zip(*columns):
        cells = [(f"{value:.15E}" if exponent and math.isfinite(value) else _fmt(value)) for value in row]
        lines.append(" , ".join(cells) + (" , " if trailing_comma else ""))
    return newline.join(lines) + newline


def build_cases() -> list[dict]:
    r = np.round(np.arange(0.5, 6.0, 0.25), 6)
    g = np.sin(2.0 * r) * np.exp(-0.2 * r)
    q = np.round(np.arange(0.4, 8.0, 0.4), 6)
    f = np.cos(q) * np.exp(-0.3 * q) - 0.1
    s = 1.0 + f
    tof = np.round(np.arange(1000.0, 2000.0, 50.0), 3)
    peak = 5.0 + 40.0 * np.exp(-(((tof - 1500.0) / 60.0) ** 2))
    k = np.round(np.arange(3.0, 6.0, 0.25), 6)
    chi = np.sin(3.0 * k) * k * 0.1

    masked = g.copy()
    masked[3:6] = np.nan
    chi_masked = chi.copy()
    chi_masked[4] = np.nan

    return [
        {   # neutron G(r), roles named, CRLF line endings
            "name": "NPDF_PDF1.csv",
            "text": _table("r(A), G(r)_RMC, G(r)_Expt", [r, 0.9 * g + 0.01, g], newline="\r\n"),
            "truth": {"calculated": 1, "experimental": 2},
        },
        {   # second neutron bank, NaN-masked experiment
            "name": "NPDF_PDF2.csv",
            "text": _table("r, calc, expt", [r, 0.95 * g, masked]),
            "truth": {"calculated": 1, "experimental": 2},
        },
        {   # neutron S(Q), trailing commas on the data rows
            "name": "NSQ_SQ1.csv",
            "text": _table("Q, S(Q)_RMC, S(Q)_Expt", [q, 0.97 * s, s], trailing_comma=True),
            "truth": {"calculated": 1, "experimental": 2},
        },
        {   # second reciprocal dataset written as F(Q), experiment first
            "name": "NSQ_SQ2.csv",
            "text": _table("Q, F(Q)_Expt, F(Q)_RMC", [q, f, 0.8 * f]),
            "truth": {"calculated": 2, "experimental": 1},
        },
        {   # time-of-flight Bragg profile with extra background / difference columns
            "name": "TOF_bragg.csv",
            "text": _table("Flight time (us), Experiment, RMC, Background, Difference",
                           [tof, peak, 0.9 * peak + 0.5, np.full_like(tof, 5.0), peak - (0.9 * peak + 0.5)]),
            "truth": {"calculated": 2, "experimental": 1},
        },
        {   # Q-axis Bragg bank, observed/calculated naming
            "name": "Q_bragg_1.csv",
            "text": _table("Q or theta, Iobs, Icalc", [q, 3.0 + s, 3.1 + 0.9 * s]),
            "truth": {"calculated": 2, "experimental": 1},
        },
        {   # x-ray F(Q), leading blank line and E-notation
            "name": "XR_FQ1.csv",
            "text": _table("Q, F(Q)_RMC, F(Q)_Expt", [q, 0.7 * f, f], leading_blank=True, exponent=True),
            "truth": {"calculated": 1, "experimental": 2},
        },
        {   # EXAFS k-space output: descriptive title row, one NaN-masked row
            "name": "Nb-EXAFS-1_Q_OUTPUT.csv",
            "text": _table("      k    ,  calculated  ,  experiment", [k, 0.9 * chi, chi_masked],
                           title=" EXAFS #1,   chi(k)*k^2"),
            "truth": None,
        },
        {   # EXAFS R-space output: header first
            "name": "Nb-EXAFS-1_R_OUTPUT.csv",
            "text": _table("r, Re_Calc, Im_Calc, Mod_Calc, Re_Ex, Im_Ex, Mod_Ex",
                           [r, g, 0.5 * g, np.abs(g), 0.9 * g, 0.4 * g, 0.9 * np.abs(g)]),
            "truth": None,
        },
        {   # partial g(r), trailing commas
            "name": "RUN_PDFpartials.csv",
            "text": _table("r (Ang), A-A, A-B, B-B", [r, 1.0 + g, 1.0 - 0.5 * g, 1.0 + 0.2 * g], trailing_comma=True),
            "truth": None,
        },
    ]


def _finite(values: list) -> list[float]:
    return [value for value in values if value is not None and math.isfinite(value)]


def summarize(data: dict, metadata: dict) -> dict:
    """The comparable part of a plot payload (floats as JSON numbers, NaN as null)."""
    series = []
    for entry in data["series"]:
        finite = _finite(entry["y"])
        series.append({
            "label": entry["label"],
            "n": len(entry["y"]),
            "finite": len(finite),
            "first": finite[0] if finite else None,
            "last": finite[-1] if finite else None,
        })
    return {
        "kind": data["kind"],
        "title": data["title"],
        "metadataTitle": metadata["title"],
        "xLabel": data["xLabel"],
        "yLabel": data["yLabel"],
        "metrics": data["metrics"],
        "series": series,
    }


def flask_payloads() -> list[dict]:
    """Every case through the real Flask endpoints."""
    backend = ROOT / "web_app" / "backend"
    if str(backend) not in sys.path:
        sys.path.insert(0, str(backend))
    os.environ.setdefault("RMC_TOOLKITS_DATA_ROOT", str(ROOT))
    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rmc_toolkits_matplotlib"))
    import app as backend_app  # noqa: E402

    client = backend_app.app.test_client()
    cases = []
    with tempfile.TemporaryDirectory(dir=ROOT) as tmpdir:
        for case in build_cases():
            path = Path(tmpdir) / case["name"]
            path.write_bytes(case["text"].encode("utf-8"))
            data = client.get("/api/plot/data", query_string={"path": str(path)})
            metadata = client.get("/api/plot/metadata", query_string={"path": str(path)})
            if data.status_code != 200 or metadata.status_code != 200:
                raise RuntimeError(f"{case['name']}: {data.get_data(as_text=True)}")
            cases.append({**case, "expected": summarize(data.get_json(), metadata.get_json())})
    return cases


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({"cases": flask_payloads()}, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
