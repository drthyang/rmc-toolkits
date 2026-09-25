# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Golden fixture for the static-mode JS engine parity tests.

Runs the Python Auto StoG engine on the shared synthetic model and writes
web_app/frontend/src/__tests__/fixtures/autoscale_fixture.json. Regenerate
whenever the engine's math changes:

    .venv/bin/python tests/generate_autoscale_fixture.py
"""

import json
from pathlib import Path
import re

import numpy as np

from rmc_toolkits.scaling import (
    ScalingConfig,
    _place_low_r_window,
    auto_enforcement_cutoff,
    autoscale,
    detect_first_peak_onset,
    first_shell_candidates,
    first_shell_foot,
    estimate_rho0,
    level_sweep,
    scale_pipeline,
)
from rmc_toolkits.transforms import fq_to_sq, g_to_gpdf, gpdf_to_fq

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "web_app" / "frontend" / "src" / "__tests__" / "fixtures" / "autoscale_fixture.json"

RHO0, B2 = 0.05, 0.02
A_TRUE, B_TRUE = 10.0, -9.0


def synthetic_g(r):
    """Same synthetic model as tests/test_scaling.py (kept in sync)."""
    onset = 0.5 * (1.0 + np.tanh((r - 2.65) / 0.07))
    peak = 1.6 * np.exp(-0.5 * ((r - 2.8) / 0.15) ** 2)
    return onset + peak


def _gauss(r, centre, sigma):
    return np.exp(-0.5 * ((r - centre) / sigma) ** 2)


def _continuum(r, start):
    return 0.5 * (1.0 + np.tanh((r - start) / 0.08))


def detector_cases() -> dict:
    """First-shell detector parity cases (models of tests/test_stog_a_detection.py and
    tests/test_stog_a_placement.py)."""
    from test_stog_a_placement import mn3sn_lobe_g

    r = np.arange(1, 801) * 0.01
    ripple = 0.9 * np.sin(2 * np.pi * (r - 1.3) / 0.24) * ((r > 1.3) & (r < 2.55))
    profiles = {
        "invertedFirst": -2.8 * _gauss(r, 1.95, 0.07) + 10.0 * _gauss(r, 2.76, 0.09)
        + _continuum(r, 3.4),
        "weakFirst": 0.6 * _gauss(r, 2.1, 0.08) + 5.6 * _gauss(r, 2.9, 0.09)
        + _continuum(r, 3.6),
        "rippleField": ripple - 3.7 * _gauss(r, 2.84, 0.06) + 1.8 * _gauss(r, 4.0, 0.08)
        + _continuum(r, 4.4),
        "shellAtSearchStart": 6.0 * _gauss(r, 1.37, 0.045) + 3.0 * _gauss(r, 2.4, 0.08)
        + _continuum(r, 3.0),
        # Review follow-up: the real Mn3Sn sub-shell lobe (2.9x, 35 %).
        "mn3snLobe": mn3sn_lobe_g(r),
        "twoShell": 1.5 * _gauss(r, 1.95, 0.06) + 4.0 * _gauss(r, 2.76, 0.08)
        + _continuum(r, 3.4),
    }
    profiles = {name: np.round(g, 12) for name, g in profiles.items()}
    cases = []
    for name, g in profiles.items():
        for qmax in (28.0, 0.0):
            onset = detect_first_peak_onset(r, g, qmax, search_min=1.3)
            candidates = first_shell_candidates(r, g, qmax, search_min=1.3)
            cases.append({
                "name": name, "qmax": qmax, "onset": onset, "candidates": candidates,
            })
    return {
        "r": r.tolist(),
        "searchMin": 1.3,
        "profiles": {name: g.tolist() for name, g in profiles.items()},
        "cases": cases,
    }


def window_cases() -> dict:
    """Low-r window placement parity cases (models from tests/test_stog_a_window.py).

    Composition-free configs (b_sq_avg unset) so the engines' S(0) handling in
    the self-consistent loop cannot differ; coarse grids keep the JS test fast.
    """
    from test_stog_a_window import B2O3_LIKE, PEROVSKITE, crystal_sq, shell_sq

    q = np.arange(17, 1001) * 0.03  # 0.51 .. 30.0
    cases = []
    sq, values = crystal_sq("SrTiO3", PEROVSKITE, 3.905, q=q)
    formula, rho0, shells, r_continuum = B2O3_LIKE
    sq_short, values_short = shell_sq(formula, rho0, shells, r_continuum, q=q)
    for name, data, values in (("srtio3", sq, values), ("shortBond", sq_short, values_short)):
        config = {
            "qmin": values["qmin"], "qmax": values["qmax"], "rho0": values["rho0"],
            "bAvgSq": values["b_avg_sq"], "rmax": 25.0, "nr": 1000,
        }
        py_config = ScalingConfig(
            qmin=config["qmin"], qmax=config["qmax"], rho0=config["rho0"],
            b_avg_sq=config["bAvgSq"], rmax=25.0, nr=1000,
        )
        case = {"name": name, "config": config, "sqMeas": data.tolist()}
        try:
            result = autoscale(q, data, py_config)
            case["expected"] = {
                "a": result.a,
                "b": result.b,
                "r0Detected": result.provenance["r0_detected"],
                "rFitWindow": list(result.provenance["r_fit_window"]),
            }
        except ValueError as exc:
            case["error"] = str(exc)
        cases.append(case)
    return {"q": q.tolist(), "aTrue": 10.0, "cases": cases}


#: Numbers in an error message (not the digit of "r0"): the Python and JS messages
#: word units differently (A / Å, r0 / r₀) but carry the same numbers.
NUMBER = re.compile(r"(?<![A-Za-z_\d.])-?\d+(?:\.\d+)?")


def placement_cases() -> dict:
    """Low-r window placement loop on scripted passes (tests/test_stog_a_placement.py).

    The JS test replays each scenario with the same scripted-pass rules; both
    engines detect the candidates on the same (rounded) profiles.
    """
    from test_stog_a_placement import (
        ERROR_KINDS, PROFILES, R, SCENARIOS, scenario_config, scripted_pass,
    )

    profiles = {name: np.round(build(R), 12) for name, build in PROFILES.items()}
    cases = []
    for name, scenario in SCENARIOS.items():
        run, calls = scripted_pass(scenario, profiles)
        case = {"name": name, **{key: scenario[key] for key in ("qmax", "trials", "refits", "default")}}
        try:
            result = _place_low_r_window(run, scenario_config(scenario))
            case["expected"] = {
                "a": result.a,
                "r0Detected": result.provenance["r0_detected"],
                "rFitWindow": list(result.provenance["r_fit_window"]),
            }
        except ValueError as exc:
            message = str(exc)
            case["error"] = {
                "kind": next(kind for kind in ERROR_KINDS if kind in message),
                "numbers": [float(value) for value in NUMBER.findall(message)],
                "message": message,
            }
        case["refitOnsets"] = [call[0] for call in calls if call[0] is not None]
        cases.append(case)
    return {
        "r": R.tolist(),
        "profiles": {name: g.tolist() for name, g in profiles.items()},
        "config": {"qmin": 0.5, "rho0": 0.05, "bAvgSq": 1.0},
        "cases": cases,
    }


def enforcement_cases() -> dict:
    """Automatic enforcement cutoff parity (models from tests/test_stog_a_enforcement.py)."""
    from test_stog_a_enforcement import exact_sq

    cases = []
    for sigma, qmax, lorch, r0 in (
        (0.10, 26.0, False, None), (0.15, 26.0, True, None), (0.08, 40.0, False, None),
        (0.10, 26.0, False, 2.2),  # a pinned r0 below the detected onset caps the cutoff
    ):
        q, sq = exact_sq(sigma, qmax)
        config = ScalingConfig(
            qmin=0.01, qmax=qmax, rho0=RHO0, b_avg_sq=1.0, lorch=lorch, rmax=20.0, nr=2000,
            r0=r0,
        )
        result = scale_pipeline(q, sq, config, 1.0, 0.0)
        keep = result.r <= 6.5
        r, g = result.r[keep], result.g_filtered[keep]
        onset = detect_first_peak_onset(r, g, qmax, search_min=config.r_cutoff + 0.3)
        cases.append({
            "config": {"qmax": qmax, "rCutoff": config.r_cutoff, "r0": r0},
            "r": r.tolist(),
            "g": g.tolist(),
            "onset": onset,
            "foot": first_shell_foot(r, g, onset),
            "cutoff": auto_enforcement_cutoff(r, g, config),
        })
    return {"cases": cases}


def main() -> None:
    q = np.arange(20, 981) * 0.03  # 0.60 .. 29.40
    r = np.arange(1, 12001) * 0.005
    gpdf = g_to_gpdf(r, synthetic_g(r), RHO0)
    sq_true = fq_to_sq(q, gpdf_to_fq(r, gpdf, q))
    sq_meas = (sq_true - B_TRUE) / A_TRUE

    base = dict(
        qmin=0.6, qmax=30.0, rho0=RHO0, b_avg_sq=B2, r_cutoff=1.0, r0=2.5,
        r_fit_min=1.2, r_fit_max=2.25, rmax=25.0, nr=1000,
    )
    config = ScalingConfig(**base)
    auto = autoscale(q, sq_meas, config)
    sweep = level_sweep(q, sq_meas)

    head = q <= q[0] + 1.0
    _, s_true_0 = np.polyfit(q[head], sq_true[head], 1)
    fz_config = ScalingConfig(
        **base, b_sq_avg=float(B2 * (1.0 - s_true_0)), amplitude_criterion="fz"
    )
    fz = autoscale(q, sq_meas, fz_config)

    detect_config = ScalingConfig(**{**base, "r0": None, "r_fit_max": None})
    detected = autoscale(q, sq_meas, detect_config)

    # Density self-consistency, seeded away from the truth (rho0 = 0.02).
    est_config = ScalingConfig(
        **{**base, "rho0": 0.02}, b_sq_avg=float(B2 * (1.0 - s_true_0))
    )
    estimate = estimate_rho0(q, sq_meas, est_config)

    manual = scale_pipeline(q, sq_meas, config, A_TRUE, B_TRUE)
    r_sample_idx = [50, 200, 500, 999]   # on the r grid (nr = 1000)
    q_sample_idx = [50, 200, 500, 950]   # on the cropped q grid (961 pts)

    payload = {
        "note": "generated by tests/generate_autoscale_fixture.py — do not edit",
        "aTrue": A_TRUE,
        "bTrue": B_TRUE,
        "config": {
            "qmin": 0.6, "qmax": 30.0, "rho0": RHO0, "bAvgSq": B2,
            "rCutoff": 1.0, "r0": 2.5, "rFitMin": 1.2, "rFitMax": 2.25,
            "rmax": 25.0, "nr": 1000,
        },
        "fzBSqAvg": float(B2 * (1.0 - s_true_0)),
        "q": q.tolist(),
        "sqMeas": sq_meas.tolist(),
        "expected": {
            "auto": {
                "a": auto.a,
                "b": auto.b,
                "iterations": auto.iterations,
                "converged": bool(auto.converged),
                "lowRRms": auto.low_r_rms,
                "c1TailMean": auto.c1_tail_mean,
            },
            "sweep": {
                "level": sweep.level,
                "levelUncertainty": sweep.level_uncertainty,
                "qLo": sweep.q_lo,
                "qHi": sweep.q_hi,
                "nAdmissible": sweep.n_admissible,
            },
            "fz": {"a": fz.a, "b": fz.b},
            "rho0Estimate": {
                "rho0": estimate["rho0"],
                "converged": estimate["converged"],
                "iterations": estimate["iterations"],
                "concordance": estimate["concordance"],
            },
            "autoDetect": {
                "a": detected.a,
                "r0Detected": detected.provenance.get("r0_detected"),
                "windowRefined": bool(detected.provenance.get("window_refined", False)),
                "rFitWindow": list(detected.provenance["r_fit_window"]),
            },
            "detector": detector_cases(),
            "window": window_cases(),
            "enforcement": enforcement_cases(),
            "placement": placement_cases(),
            "manual": {
                "lowRRms": manual.low_r_rms,
                "c1TailMean": manual.c1_tail_mean,
                "rSampleIdx": r_sample_idx,
                "qSampleIdx": q_sample_idx,
                "gkSamples": [manual.gk[i] for i in r_sample_idx],
                "sqFilteredSamples": [manual.sq_filtered[i] for i in q_sample_idx],
            },
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload))
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")
    print(f"auto a={auto.a} b={auto.b} iters={auto.iterations}; fz a={fz.a}; "
          f"level={sweep.level} n_adm={sweep.n_admissible}")


if __name__ == "__main__":
    main()
