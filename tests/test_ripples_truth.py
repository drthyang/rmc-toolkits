# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Known-truth acceptance test for check_ripples (rmc_toolkits/ripples.py), on synthetic data only.

The data are built here with numpy alone, independent of the engine's own helpers:

* the structure is fcc Ni (a = 3.524 Å, public), its coordination shells to 16 Å
  ("crowded": 22 shells in the 0–12 Å grid, above the 12 strong peaks the engine
  fits) or to 6 Å ("uncrowded": 5 shells);
* D(r) = Σ_j (N_j / r_j) [g(r − r_j) − g(r + r_j)], g a unit-area Gaussian of
  width σ; its F(Q) = Σ_j (N_j / r_j) e^{−Q²σ²/2} sin(Q r_j) is exact;
* the data are (2/π) Σ_k ΔQ F(Q_k) sin(Q_k r) over Q_k = ΔQ … Q_max (the
  rectangle rule, as RMCProfile), and the true termination artefact is the data
  less the untruncated D(r).

The engine is run as ``assess_run`` runs it: ``analyse_series`` per dataset, then
``find_features`` over 1.5–8 Å or ``classify_features`` at given positions.

Three groups, one test per case:

* A, the ripples are found: asked at every true ripple of at least 1 % of the
  largest peak (outside the peaks' cores, 1.5–8 Å), at least 90 % read
  "ripple-consistent"; every feature found by the engine itself lies at a true
  ripple (a quarter period) — or is called something other than a ripple;
* B, quiet without ripples: no feature is called "ripple-consistent" when the
  data have no ripples, and nothing smaller than 1 % of the largest peak
  (the engine's own FEATURE_MIN_FRACTION) is called one;
* C, real features are not called ripples: a real peak on or near a ripple is
  not "ripple-consistent" at one Q_max when it is at least 2× the ripple, and
  never "moves-with-qmax" across Q_max 20 and 28 Å⁻¹; bare ripples (controls)
  are still read as ripples.

The cases that fail on the engine as first proposed (PR #144) are marked
``expectedFailure`` (KNOWN_FAILURES); a fix that makes one pass turns it into an
unexpected success, and the mark is then removed.
"""

from __future__ import annotations

import itertools
import math
import unittest

import numpy as np

from rmc_toolkits import ripples as R

A_NI = 3.524
DQ = 0.01
R_GRID = np.round(np.arange(1, 1201) * 0.01, 10)  # 0.01 … 12 Å
R_FIRST = A_NI / math.sqrt(2)
WINDOW = (1.5, 8.0)
FEATURE_SIGMA = 0.07
MIN_FRACTION = 0.01  # the engine's FEATURE_MIN_FRACTION, the size below which nothing is a feature
RECALL = 0.9


def shells(rmax):
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    n = int(rmax / A_NI) + 2
    points = np.array([(np.array(cell) + b) * A_NI
                       for cell in itertools.product(range(-n, n + 1), repeat=3) for b in basis])
    d = np.linalg.norm(points, axis=1)
    d = d[(d > 1e-6) & (d <= rmax)]
    rs, counts = np.unique(np.round(d, 6), return_counts=True)
    return list(zip(rs.tolist(), counts.tolist()))


STRUCTURES = {"uncrowded": shells(6.0), "crowded": shells(16.0)}


def peaks(structure, sigma, extra=()):
    return [(r0, n / r0, sigma) for r0, n in STRUCTURES[structure]] + list(extra)


def d_true(pk):
    y = np.zeros_like(R_GRID)
    for r0, amp, s in pk:
        c = amp / (s * math.sqrt(2 * math.pi))
        y += c * (np.exp(-(R_GRID - r0) ** 2 / (2 * s * s)) - np.exp(-(R_GRID + r0) ** 2 / (2 * s * s)))
    return y


def d_measured(pk, qmax):
    q = np.arange(1, int(round(qmax / DQ)) + 1) * DQ
    f = np.zeros_like(q)
    for r0, amp, s in pk:
        f += amp * np.exp(-q * q * s * s / 2) * np.sin(q * r0)
    return (2 / math.pi) * DQ * (np.sin(np.outer(R_GRID, q)) @ f)


def extrema(y):
    i = np.arange(1, y.size - 1)
    return i[((y[i] > y[i - 1]) & (y[i] >= y[i + 1])) | ((y[i] < y[i - 1]) & (y[i] <= y[i + 1]))]


def true_ripples(pk, qmax):
    """(r, amplitude) of the artefact's extrema outside the cores, in WINDOW, ≥ MIN_FRACTION of the largest peak."""
    clean = d_true(pk)
    artefact = d_measured(pk, qmax) - clean
    top = float(np.max(np.abs(clean)))
    out = []
    for i in extrema(artefact):
        r = float(R_GRID[i])
        if WINDOW[0] <= r <= WINDOW[1] and all(abs(r - r0) >= 2.5 * s for r0, _, s in pk) \
                and abs(artefact[i]) >= MIN_FRACTION * top:
            out.append((r, float(artefact[i])))
    return out


def analyse(y, qmax):
    return R.analyse_series(R_GRID, y, space="D", q_min=DQ, q_max=qmax, dq=DQ, rule="rectangle", window="none",
                            r_first=R_FIRST)


def found_features(analysis):
    """``assess_run``'s r_range path for one dataset."""
    features = R.find_features(analysis, *WINDOW)
    return R.classify_features([analysis], features) if features else []


def tolerance(qmax):
    return R.RIPPLE_POSITION_FRACTION * 2 * math.pi / qmax


# ---------------------------------------------------------------------------
# The cases
# ---------------------------------------------------------------------------

A_CASES = [(s, qmax, sigma) for s in ("uncrowded", "crowded")
           for qmax, sigma in ((20.0, 0.05), (20.0, 0.07), (20.0, 0.09), (15.0, 0.07), (25.0, 0.05))]

C_CASES = (  # name, which true ripple (the largest + or −, at Q_max 20, σ 0.07), offset in periods, height × |ripple|
    ("on_plus_3x", "+", 0.0, 3.0),
    ("on_plus_1x", "+", 0.0, 1.0),
    ("on_minus_1x", "-", 0.0, 1.0),
    ("quarter_off_plus_2x", "+", 0.25, 2.0),
)

KNOWN_FAILURES = {
    # A: the crowded G(r): the global cap of 12 strong peaks skips shells (4.98 Å at Q_max 20), and some
    # fits reach the σ floor; recall 9/14 – 17/19.
    "test_A_crowded_q20_s0.05", "test_A_crowded_q20_s0.07", "test_A_crowded_q20_s0.09",
    "test_A_crowded_q15_s0.07", "test_A_crowded_q25_s0.05",
    # B: wiggles of 1e-5 of the largest peak called ripple-consistent (no floor when the noise is unknown).
    "test_B_uncrowded_q60_asked_at_q20_ripples", "test_B_crowded_q60_asked_at_q20_ripples",
    # C: real features called ripples, at one Q_max or across Q_max 20 and 28.
    "test_C_uncrowded_on_plus_3x", "test_C_uncrowded_on_plus_1x", "test_C_uncrowded_quarter_off_plus_2x",
    "test_C_crowded_on_plus_1x", "test_C_crowded_on_minus_1x",
}


class KnownTruth(unittest.TestCase):
    """Termination ripples of known Q_max in synthetic fcc Ni; see the module docstring."""

    maxDiff = None


def _mark(name, fn):
    fn.__name__ = name
    setattr(KnownTruth, name, unittest.expectedFailure(fn) if name in KNOWN_FAILURES else fn)


def _a_case(structure, qmax, sigma):
    def test(self):
        pk = peaks(structure, sigma)
        analysis = analyse(d_measured(pk, qmax), qmax)
        truth = true_ripples(pk, qmax)
        self.assertGreater(len(truth), 10)
        classes = [R.classify_feature(analysis, r)["class"] for r, _ in truth]
        recall = classes.count("ripple-consistent") / len(classes)
        self.assertGreaterEqual(recall, RECALL, f"ripple-consistent at {classes.count('ripple-consistent')} of "
                                                f"{len(classes)} true ripples: {classes}")
        off = [(round(e["asked"], 3), e["datasets"][0]["class"]) for e in found_features(analysis)
               if e["datasets"][0]["class"] == "ripple-consistent"
               and min(abs(e["asked"] - r) for r, _ in truth) > tolerance(qmax)]
        self.assertEqual(off, [], "found and called a ripple where there is none")
        far = [(round(e["asked"], 3), e["datasets"][0]["class"]) for e in found_features(analysis)
               if min(abs(e["asked"] - r) for r, _ in truth) > tolerance(qmax)
               and all(abs(e["asked"] - r0) > 2.5 * s for r0, _, s in pk)]
        self.assertEqual(far, [], "found a feature at neither a true ripple nor a true peak")
    return test


def _b_cases(structure):
    pk = peaks(structure, 0.07)
    tests = {}

    def removed(self):
        analysis = analyse(d_true(pk), 20.0)
        called = [round(e["asked"], 3) for e in found_features(analysis)
                  if e["datasets"][0]["class"] == "ripple-consistent"]
        self.assertEqual(called, [])
    tests[f"test_B_{structure}_ripples_removed"] = removed

    def high_qmax(self):
        analysis = analyse(d_measured(pk, 60.0), 60.0)
        called = [round(e["asked"], 3) for e in found_features(analysis)
                  if e["datasets"][0]["class"] == "ripple-consistent"]
        self.assertEqual(called, [])
    tests[f"test_B_{structure}_q60_found"] = high_qmax

    def asked(self):
        analysis = analyse(d_measured(pk, 60.0), 60.0)
        top = float(np.max(np.abs(d_true(pk))))
        tiny = [(r, R.classify_feature(analysis, r)) for r, _ in true_ripples(pk, 20.0)]
        called = [round(r, 3) for r, reading in tiny
                  if reading["class"] == "ripple-consistent" and abs(reading["size"]) < MIN_FRACTION * top]
        self.assertEqual(called, [], "below 1 % of the largest peak, called a ripple")
    tests[f"test_B_{structure}_q60_asked_at_q20_ripples"] = asked
    return tests


def _c_case(structure, which, offset, multiple):
    def test(self):
        pk0 = peaks(structure, 0.07)
        truth = true_ripples(pk0, 20.0)
        r_at, amplitude = (max if which == "+" else min)(truth, key=lambda t: t[1])
        r_at += offset * 2 * math.pi / 20.0
        height = multiple * abs(amplitude)
        pk = peaks(structure, 0.07, [(r_at, height * FEATURE_SIGMA * math.sqrt(2 * math.pi), FEATURE_SIGMA)])
        low, high = analyse(d_measured(pk, 20.0), 20.0), analyse(d_measured(pk, 28.0), 28.0)
        both = R.classify_features([low, high], [r_at])[0]
        self.assertNotEqual(both["across"]["verdict"], "moves-with-qmax",
                            f"a real feature at {r_at:.3f} Å: {[d['class'] for d in both['datasets']]}")
        if multiple >= 2.0:
            self.assertNotEqual(both["datasets"][0]["class"], "ripple-consistent",
                                f"a real feature {multiple}× the ripple at {r_at:.3f} Å, at Q_max 20")
    return test


def _c_control(structure, which):
    def test(self):
        pk = peaks(structure, 0.07)
        r_at, _ = (max if which == "+" else min)(true_ripples(pk, 20.0), key=lambda t: t[1])
        low, high = analyse(d_measured(pk, 20.0), 20.0), analyse(d_measured(pk, 28.0), 28.0)
        both = R.classify_features([low, high], [r_at])[0]
        self.assertEqual(both["datasets"][0]["class"], "ripple-consistent")
        self.assertIn(both["across"]["verdict"], ("moves-with-qmax", "ripple-in-every-dataset", "mixed"))
    return test


for _s, _q, _sg in A_CASES:
    _mark(f"test_A_{_s}_q{_q:g}_s{_sg}", _a_case(_s, _q, _sg))
for _s in STRUCTURES:
    for _name, _fn in _b_cases(_s).items():
        _mark(_name, _fn)
    for _name, _which, _offset, _multiple in C_CASES:
        _mark(f"test_C_{_s}_{_name}", _c_case(_s, _which, _offset, _multiple))
    for _which, _label in (("+", "plus"), ("-", "minus")):
        _mark(f"test_C_{_s}_control_bare_{_label}", _c_control(_s, _which))


if __name__ == "__main__":
    unittest.main()
