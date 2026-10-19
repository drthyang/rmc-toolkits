# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Compare a regenerated Python golden with the committed fixture.

The rule ``web_app/frontend/src/__tests__/goldenCompare.js`` applies to the JS
engines, so a fixture that passes here is the same interface there:

* strings, booleans and None exactly;
* an int on both sides exactly (counts, cell indices, frequencies);
* any other number to ``RTOL`` relative plus ``ATOL`` absolute, or to ``RTOL``
  relative alone for a ``relative`` key (p-values reach 1e-142);
* an eigenvector set (a ``bases`` key, given the descending eigenvalues of
  its rows) eigenspace by eigenspace to ``AXIS_TOL``: a separated axis up to
  its sign, a group whose eigenvalues lie within ``DEGENERATE_GAP`` of the
  largest by the projector onto its span (LAPACK builds pick different bases
  of a degenerate eigenspace);
* dicts by their whole key set.

The tolerances absorb what differs between the supported NumPy/SciPy
releases (summation order, LAPACK, the incomplete gamma's tails); a change of
the engine's result moves values far more.
"""

from __future__ import annotations

import math

import numpy as np

RTOL = 1e-9
ATOL = 1e-12
AXIS_TOL = 1e-9
DEGENERATE_GAP = 1e-6


def eigenspaces(spectrum) -> list[list[int]]:
    """Groups of row indices whose (descending) eigenvalues are degenerate."""
    scale = max(abs(float(spectrum[0])), 1e-300)
    groups = [[0]]
    for index in range(1, len(spectrum)):
        if float(spectrum[index - 1]) - float(spectrum[index]) > DEGENERATE_GAP * scale:
            groups.append([index])
        else:
            groups[-1].append(index)
    return groups


def _basis_misfit(got, want, spectrum) -> float:
    got = np.asarray(got, dtype=float)
    want = np.asarray(want, dtype=float)
    misfit = 0.0
    for group in eigenspaces(spectrum):
        if len(group) == 1:
            row = group[0]
            misfit = max(misfit, min(np.abs(got[row] - want[row]).max(), np.abs(got[row] + want[row]).max()))
        else:
            difference = got[group].T @ got[group] - want[group].T @ want[group]
            misfit = max(misfit, float(np.abs(difference).max()))
    return misfit


def assert_same(case, got, want, path: str = "$", *, relative=frozenset(), bases=None, key=None,
                message: str = "") -> None:
    """Assert ``got`` (regenerated) matches ``want`` (committed) under the rule above."""
    bases = bases or {}
    where = f"{message} ({path})" if message else path
    if want is None or isinstance(want, (bool, str)):
        case.assertEqual(got, want, where)
        return
    if isinstance(want, (int, float)):
        case.assertIsInstance(got, (int, float), where)
        case.assertNotIsInstance(got, bool, where)
        if isinstance(want, int) and isinstance(got, int):
            case.assertEqual(got, want, where)
            return
        bound = RTOL * abs(want) + (1e-300 if key in relative else ATOL)
        case.assertTrue(math.isfinite(got), where)
        case.assertLessEqual(abs(got - want), bound, f"{where}: {got!r} vs {want!r}")
        return
    if isinstance(want, list):
        case.assertIsInstance(got, list, where)
        case.assertEqual(len(got), len(want), f"{where} length")
        if key is not None and key in bases:
            case.assertLessEqual(_basis_misfit(got, want, bases[key]), AXIS_TOL, f"{where} eigenbasis")
            return
        for index, (got_item, want_item) in enumerate(zip(got, want)):
            assert_same(case, got_item, want_item, f"{path}[{index}]", relative=relative, key=key, message=message)
        return
    if isinstance(want, dict):
        case.assertIsInstance(got, dict, where)
        case.assertEqual(sorted(got), sorted(want), f"{where} keys")
        for name in want:
            assert_same(case, got[name], want[name], f"{path}.{name}", relative=relative, bases=bases,
                        key=name, message=message)
        return
    case.fail(f"{where}: unexpected value {want!r}")
