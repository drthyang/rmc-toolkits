# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Displacement modulations: is there order that enlarges the cell?

The average structure folds every unit cell of the supercell into one, so it
cannot show order whose period is longer than one cell -- the cell doubling of
GaNb4Se8 below 29 K (P2_12_12_1, c' = 2c) or of GaTa4Se8 below 50 K, octahedral
tilts, charge-density waves. The configuration itself can: an ordered
displacement pattern with wavevector q puts power into the Fourier component
of every site's displacement field at q.

For each reference site s with copies in the cells c of an N1 x N2 x N3
supercell, the displacements u_s(c) (from the site mean, Angstrom) are Fourier
transformed over the cell grid,

    U_s(q) = sum_c u_s(c) exp(-2 pi i q.c),   q = (h/N1, k/N2, l/N3),

and the power P_s(q) = |U_s(q)|^2 / (n_s tr C_s) is normalised by the copy
count and the site's total displacement variance, so independent (white)
displacements give 1 on average at every q != 0. The powers of an element's
sites are averaged. Two kinds of q are reported against the empirical
distribution over all q != 0 (the configuration's own diffuse background,
which carries RMC's long-wavelength correlations that a white-noise null
ignores):

* the zone-boundary points (1/2,0,0) ... (1/2,1/2,1/2) of the cell's
  reciprocal lattice -- order doubling the cell along those directions -- when
  the supercell is even along each halved axis;
* the strongest modes overall (+q and -q counted once).

Each with its ratio to the median power over all q != 0, and its percentile
and rank among the *compared* modes: one per +q/-q pair (they carry the same
power), leaving out the long-wavelength modes (every |q_i| <= ``LONG_WAVE_Q``),
the box-scale correlations RMC leaves. Counting both members of a pair, and
the long waves, pushed a strong zone-boundary order below the 99th
percentile (the 2026-10-12 review: 55x the median at percentile 98.6, under
a long wave and its partner).

Every reported point also carries its ``polarisation`` (plan item 4.3): the
direction the element's atoms move in that mode, the leading axis of the
element's power tensor at q, with the share of the power along it
(:func:`mode_polarisation`).

The engine also decides what these show (the verdict): every zone-boundary
point and every strongest mode carries its calibrated significance ``sigma``
against independent displacements (family-wise, Sidak-corrected over the
points or modes of every element the configuration has, whichever the call
names), its ``kind`` ('zone_boundary', 'long_wave', 'general'), whether it is
``ordered`` (ratio >= ``MODULATION_RATIO`` and sigma >= ``MODULATION_SIGMA``),
``because`` and its ``bar``; each element a ``verdict`` -- ``cellDoubling``,
``generalOrder``, ``longWave`` -- with the q that decided it and the nearest
miss; the payload the ``family`` the corrections used (its sizes and
elements) and the ``rules``. The browser port is ``workers/modulations.js``
(a separable DFT, the same numbers to floating-point noise; the parity
fixture pins the verdicts and polarisations too).
"""

from __future__ import annotations

import math
from itertools import product

import numpy as np
from scipy.special import gammaincc, ndtri

from rmc_toolkits.pca_kde import SiteDisplacements

#: How many of the strongest modes to report per element.
TOP_MODES = 5
#: Relative tolerance on "lower power" for the percentile: Friedel pairs
#: P(q) = P(-q) differ only by rounding, and must count alike in both engines.
PERCENTILE_RTOL = 1e-9
#: A mode with every component within this of the zone centre (fractions of
#: the cell's reciprocal axes) is a long-wavelength wave, left out of the
#: comparison.
LONG_WAVE_Q = 0.25
#: The long-wave test's tolerance, |q_i| <= LONG_WAVE_Q + LONG_WAVE_TOL. A grid
#: component k/N is 1/4 exactly or at least 1/(4N) from it, so on any grid the
#: tolerance changes nothing (tests/test_modulations.py checks every grid
#: point up to N = 256): it only keeps the test robust to rounding.
LONG_WAVE_TOL = 1e-12
#: A polarisation is an axis (its sign is arbitrary): it is given in the sense
#: that makes the closest crystal axis acute, as the PCA tables give an axis
#: (pcaCrystalFrame.js crystalOrientationRows). Axes whose |cos| lies within
#: this relative tolerance of the largest tie, and the first of them decides,
#: so both engines pick the same sense for a direction equally close to two.
AXIS_TIE_RTOL = 1e-9

# --- the verdict's rules --------------------------------------------------------
#
# A zone-boundary point is ordered displacement when its power is at least
# MODULATION_RATIO times the median over all wavevectors (the effect size) AND
# it is significant at MODULATION_SIGMA after a family-wise (Sidak) correction
# over every zone-boundary point of every element the configuration has --
# whichever elements the call names: choosing the element to ask about after
# seeing the data is selective inference, so the family never shrinks with the
# request (finding F3, 2026-10-19; until then it was the elements of the call,
# and naming one element of three took a 2.98 sigma point to 3.30 sigma). The null: for independent
# displacements a site's normalised power at a half-integer q (a real Fourier
# sum) is sum_c w_c chi2_1 with sum w_c = 1, so an element's mean over m sites
# has at least m degrees of freedom: P >= chi2_m / m is the heaviest tail any
# site anisotropy allows (conservative). The same holds at a general q (not a
# zone-boundary point, not a long wave) with a complex sum: P >= chi2_2m / 2m,
# corrected over the compared modes. The earlier rule (top 1 % of the modes, or
# rank 1) had no correction: the validation cases found "order" in 6 of 12
# configurations without any. These were the Copilot tool's rules
# (copilotTools.js summarizeModulations) until 2026-10-19, ported unchanged.

#: Ratio to the median power an ordered point needs (the effect size).
MODULATION_RATIO = 3
#: Family-wise one-sided significance (sigma) an ordered point needs.
MODULATION_SIGMA = 3
#: An element's strongest mode this many times the median, and a long wave, is
#: a box-scale displacement wave (common in RMC), not a superstructure.
LONG_WAVE_RATIO = 5
#: sigma is clipped here: an upper tail that underflows to 0 (sigma = +inf) is
#: reported as SIGMA_MAX, which no finite double reaches (p = 5e-324 gives
#: about 38.5). A sigma of -inf (power far below the null mean) is None.
SIGMA_MAX = 40
#: The p-value is capped just below 1 before the Sidak step (as the tool did).
_P_CAP = 1.0 - 1e-16

MODULATION_RULES = {
    "ratio": MODULATION_RATIO,
    "sigma": MODULATION_SIGMA,
    "longWaveRatio": LONG_WAVE_RATIO,
    "longWaveQ": LONG_WAVE_Q,
    "sigmaMax": SIGMA_MAX,
}
#: Every ``because`` a point and each element finding can give (closed sets;
#: the parity fixture pins the strings).
MODULATION_BECAUSE = {
    "point": ("ordered", "below_ratio", "below_sigma", "undefined_statistic", "zone_boundary", "long_wave"),
    "cellDoubling": ("ordered_point", "no_ordered_point", "undefined_statistic", "untestable", "no_displacements"),
    "generalOrder": ("ordered_mode", "no_ordered_mode", "undefined_statistic", "no_general_mode", "no_displacements"),
    "longWave": ("long_wave_dominates", "below_ratio", "strongest_not_long_wave", "undefined_statistic",
                 "no_displacements"),
}
MODE_KINDS = ("zone_boundary", "long_wave", "general")


def is_long_wave(q) -> bool:
    """A long-wavelength mode: every component within LONG_WAVE_Q of the zone centre."""
    return all(abs(float(value)) <= LONG_WAVE_Q + LONG_WAVE_TOL for value in q)


def _sidak_sigma(p: float, family: float) -> float:
    """One-sided sigma of 1 - (1 - p)^family, computed without cancellation
    for small p; +inf where it underflows to 0, -inf where it rounds to 1."""
    family_p = -math.expm1(family * math.log1p(-min(p, _P_CAP)))
    if family_p <= 0:
        return math.inf
    if family_p >= 1:
        return -math.inf
    return float(-ndtri(family_p))


def _reported(sigma: float) -> float | None:
    """sigma as reported: -inf as None, +inf as SIGMA_MAX."""
    return None if sigma == -math.inf else min(sigma, SIGMA_MAX)


def zone_boundary_sigma(power, sites: int, family: int) -> float | None:
    """The family-wise significance (one-sided sigma) of an element's
    zone-boundary power ``power`` (normalised: 1 = independent displacements)
    over ``sites`` sites, among ``family`` points tested. None where undefined
    (power <= 0) or -inf; +inf is SIGMA_MAX."""
    if power is None or not (power > 0) or not (sites > 0):
        return None
    p = float(gammaincc(sites / 2.0, sites * power / 2.0))
    return _reported(_sidak_sigma(p, family))


def general_mode_sigma(power, sites: int, family: int) -> float | None:
    """The same at a general wavevector, where the Fourier sum is complex
    (P >= chi2_2m / 2m); ``family`` = the compared modes. Order there is
    correlation between cells that does not double the cell: incommensurate,
    or short-range (the 2026-10-14 review's L4)."""
    if power is None or not (power > 0) or not (sites > 0) or not (family >= 1):
        return None
    p = float(gammaincc(sites, sites * power))
    return _reported(_sidak_sigma(p, family))


def mode_polarisation(tensor, unit_vectors) -> dict | None:
    """The displacement direction of a mode, from its power tensor at q.

    ``tensor`` is the real part of the element's cross-spectral tensor at q,
    (1/m) sum_s Re[U_s(q) U_s(q)^H] / (n_s tr C_s) over its m sites, whose
    trace is the element's power there. The real part is the mode's
    mean-square displacement tensor (an elliptical mode a + ib moves along an
    ellipse with second moment aa^T + bb^T); the imaginary part is
    antisymmetric, the ellipse's circulation, and carries no axis. It is the
    same at -q.

    Returns ``direction``, the leading eigenvector (a Cartesian unit vector);
    ``uvw``, the same direction as a direct-lattice [u v w] (M^-T d, M the
    rows a, b, c of ``unit_vectors``), scaled so its largest component is 1 --
    the form of the app's other [u v w] readings (list_sites' pc1_uvw), never
    a plane normal; and ``share``, lambda_1 / tr: 1 for a linearly polarised
    mode, 1/3 when no direction is preferred. An axis has no sign: both are
    given in the sense that makes the closest crystal axis acute
    (``AXIS_TIE_RTOL``). None when the tensor holds no power.
    """
    tensor = np.asarray(tensor, dtype=float)
    cell = np.asarray(unit_vectors, dtype=float)
    trace = float(tensor[0, 0] + tensor[1, 1] + tensor[2, 2])
    if not (math.isfinite(trace) and trace > 0):
        return None
    values, vectors = np.linalg.eigh(tensor)
    direction = vectors[:, -1] / np.linalg.norm(vectors[:, -1])
    cosines = cell @ direction / np.linalg.norm(cell, axis=1)
    magnitude = np.abs(cosines)
    dominant = int(np.flatnonzero(magnitude >= magnitude.max() * (1.0 - AXIS_TIE_RTOL))[0])
    if cosines[dominant] < 0:
        direction = -direction
    uvw = np.linalg.solve(cell.T, direction)
    uvw = uvw / np.abs(uvw).max()
    # + 0.0: a -0.0 component reads 0 in both engines' JSON.
    return {
        "direction": [float(value) + 0.0 for value in direction],
        "uvw": [float(value) + 0.0 for value in uvw],
        "share": min(float(values[-1]) / trace, 1.0),
    }


def _site_variance(u: np.ndarray) -> float:
    """A site's total displacement variance tr C (population); 0 for a site that does not move."""
    return float(((u - u.mean(axis=0)) ** 2).mean(axis=0).sum())


def zone_boundary_points(supercell: np.ndarray) -> list[tuple[int, int, int]]:
    """The half-integer points (as 0/1 flags per axis) the supercell can hold."""
    points = []
    for flags in product((0, 1), repeat=3):
        if not any(flags):
            continue
        if all(supercell[axis] % 2 == 0 for axis in range(3) if flags[axis]):
            points.append(flags)
    return points


def _folded(index: int, size: int) -> float:
    """Grid index -> fractional wavevector component in (-1/2, 1/2]."""
    value = int(index) / int(size)
    return value - 1.0 if value > 0.5 else value


def element_spectrum(sites: SiteDisplacements, name: str) -> tuple[np.ndarray, np.ndarray, int, int]:
    """One element's displacement spectrum over the supercell's cell grid.

    Returns ``(power, tensor, missing, used)``: the power P(q) averaged over
    the element's sites that move (shape N1 x N2 x N3), the real part of their
    cross-spectral tensor normalised the same way (N1 x N2 x N3 x 3 x 3; its
    trace is the power, mode_polarisation reads it), the copies missing from
    the grid and the number of sites used. ``used`` is 0 (and the arrays
    zero) when no site of ``name`` moves.
    """
    supercell = np.asarray(sites.supercell, dtype=int)
    n_q = int(np.prod(supercell))
    cells = np.mod(np.asarray(sites.cell_indices, dtype=int), supercell)
    power = np.zeros(tuple(supercell))
    tensor = np.zeros(tuple(supercell) + (3, 3))
    missing = 0
    used = 0
    for s, label in enumerate(sites.elements):
        if label != name:
            continue
        rows = np.flatnonzero(sites.site_index == s)
        u = sites.displacements[rows]
        variance = _site_variance(u)
        if not variance > 0:
            continue
        grid = np.zeros(tuple(supercell) + (3,))
        grid[cells[rows, 0], cells[rows, 1], cells[rows, 2]] = u
        missing += n_q - rows.size
        transform = np.fft.fftn(grid, axes=(0, 1, 2))
        power += (np.abs(transform) ** 2).sum(axis=-1) / (rows.size * variance)
        tensor += np.einsum("...j,...k->...jk", transform, transform.conj()).real / (rows.size * variance)
        used += 1
    if used:
        power /= used
        tensor /= used
    return power, tensor, int(missing), used


def displacement_modulations(sites: SiteDisplacements, element: str | None = None, *, top: int = TOP_MODES) -> dict:
    """Per-element displacement power at every wavevector of the supercell grid, and the verdict.

    ``element`` restricts the report to one element; None reports each
    element present. The significance family is every element of the
    configuration that has displacements, whichever ``element`` names (F3:
    naming one element must not shrink the correction), so a point's sigma
    and verdict are the same either way. Raises ``ValueError`` for an unknown
    element or a file without cell indices.
    """
    if sites.cell_indices is None:
        raise ValueError("modulation analysis needs the reference-site and cell columns of the .rmc6f file")
    supercell = np.asarray(sites.supercell, dtype=int)
    n_q = int(np.prod(supercell))
    if n_q < 8:
        raise ValueError(f"a {supercell.tolist()} supercell has too few wavevectors to compare")
    labels = list(sites.elements)
    every = sorted(set(labels))
    if element is not None:
        element = element.strip().capitalize()
        if element not in labels:
            raise ValueError(f"no site holds {element}; elements: {', '.join(sites.species)}")
        elements = [element]
    else:
        elements = every
    # The family (F3): every element with a site that moves, as a call
    # without ``element`` tests them.
    moving = {name: False for name in every}
    for s, label in enumerate(labels):
        if not moving[label]:
            moving[label] = _site_variance(sites.displacements[np.flatnonzero(sites.site_index == s)]) > 0
    family_elements = [name for name in every if moving[name]]
    boundary = zone_boundary_points(supercell)
    # The compared modes: one per +q/-q pair (the smaller grid index), q != 0,
    # no long wave.
    grid_index = np.stack(np.unravel_index(np.arange(n_q), tuple(supercell)), axis=1)
    partner = np.mod(-grid_index, supercell)
    partner_linear = np.ravel_multi_index(tuple(partner.T), tuple(supercell))
    q_frac = np.array([[_folded(index[axis], supercell[axis]) for axis in range(3)] for index in grid_index])
    long_wave = np.array([is_long_wave(q) for q in q_frac], dtype=bool)
    compared = (np.arange(n_q) <= partner_linear) & (np.arange(n_q) != 0) & ~long_wave
    result = {
        "supercell": supercell.tolist(),
        "wavevectors": n_q - 1,
        "untestable": [list(flags) for flags in product((0, 1), repeat=3)
                       if any(flags) and list(flags) not in [list(b) for b in boundary]],
        "elements": [],
        "parseWarning": sites.parse_warning,
    }
    for name in elements:
        power, tensor, missing, used = element_spectrum(sites, name)
        if not used:
            result["elements"].append({"element": name, "sites": 0, "verdict": _untested_verdict()})
            continue
        flat = power.ravel()
        others = np.delete(flat, 0)
        median = float(np.median(others))
        reference = flat[compared]

        def describe(index: tuple[int, int, int]) -> dict:
            value = float(power[index])
            return {
                "q": [float(_folded(index[axis], supercell[axis])) for axis in range(3)],
                "ratio": value / median if median > 0 else None,
                "percentile": (
                    float((reference < value * (1.0 - PERCENTILE_RTOL)).mean() * 100.0) if reference.size else None
                ),
                "rank": int((reference > value * (1.0 + PERCENTILE_RTOL)).sum()) + 1 if reference.size else None,
                "polarisation": mode_polarisation(tensor[index], sites.unit_vectors),
            }

        zone = []
        for flags in boundary:
            index = tuple(supercell[axis] // 2 if flags[axis] else 0 for axis in range(3))
            entry = describe(index)
            entry["q"] = [0.5 * flag for flag in flags]
            zone.append(entry)
        order = np.argsort(-flat, kind="stable")
        strongest = []
        seen = set()
        for linear in order:
            if linear == 0:
                continue
            index = tuple(int(v) for v in np.unravel_index(int(linear), tuple(supercell)))
            partner = tuple(int((-index[axis]) % supercell[axis]) for axis in range(3))
            # +q and -q carry the same power (to rounding): report the pair
            # once, as its member with the smaller grid index, in both engines.
            key = min(index, partner)
            if key in seen:
                continue
            seen.add(key)
            mode = describe(key)
            # A zone-boundary point (every index 0 or N/2) is tested in
            # zoneBoundary; a long wave is not tested; the rest are general.
            boundary_mode = all(key[axis] == 0 or 2 * key[axis] == supercell[axis] for axis in range(3))
            mode["kind"] = "zone_boundary" if boundary_mode else ("long_wave" if is_long_wave(mode["q"]) else "general")
            strongest.append(mode)
            if len(strongest) == top:
                break
        result["elements"].append({
            "element": name,
            "sites": used,
            "missingCopies": missing,
            "medianPower": median,
            "comparedModes": int(reference.size),
            "zoneBoundary": zone,
            "strongest": strongest,
        })
    _decide(result, family_elements, len(boundary), int(compared.sum()))
    return result


def _decide(result: dict, family_elements: list[str], points: int, modes: int) -> None:
    """Annotate every tested point with its sigma, ordered, because and bar,
    and give each element its verdict, over the family: the ``points``
    zone-boundary points and ``modes`` compared modes of each of
    ``family_elements`` -- every element the configuration has with
    displacements, whichever the call reports (F3). The sizes are the same
    for every element, so this is what a call without ``element`` counted."""
    tested = [entry for entry in result["elements"] if entry["sites"]]
    family = {
        "zoneBoundaryPoints": points * len(family_elements),
        "comparedModes": modes * len(family_elements),
        "elements": list(family_elements),
    }
    for entry in tested:
        for point in entry["zoneBoundary"]:
            point["kind"] = "zone_boundary"
            _test(point, "zone", entry, family["zoneBoundaryPoints"])
        for mode in entry["strongest"]:
            if mode["kind"] == "general":
                _test(mode, "general", entry, family["comparedModes"])
            else:
                mode.update({"sigma": None, "ordered": None, "because": mode["kind"], "bar": None})
        entry["verdict"] = _element_verdict(entry)
    result["family"] = family
    result["rules"] = dict(MODULATION_RULES)


def _test(point: dict, test: str, entry: dict, family: int) -> None:
    """One point against its bar. The power is ratio x median, as the tool
    formed it: sigma far below the background is ill-conditioned (1 - (1 - p)^F
    rounds), and this keeps it bit for bit the value the tool reported."""
    ratio = point["ratio"]
    power = ratio * entry["medianPower"] if ratio is not None else None
    sites = entry["sites"]
    defined = power is not None and power > 0 and sites > 0 and (test == "zone" or family >= 1)
    sigma = None
    if not defined:
        ordered, because = None, "undefined_statistic"
    else:
        sigma = zone_boundary_sigma(power, sites, family) if test == "zone" else general_mode_sigma(power, sites, family)
        if ratio < MODULATION_RATIO:
            ordered, because = False, "below_ratio"
        elif sigma is None or sigma < MODULATION_SIGMA:
            ordered, because = False, "below_sigma"
        else:
            ordered, because = True, "ordered"
    point.update({"sigma": sigma, "ordered": ordered, "because": because,
                  "bar": {"ratio": MODULATION_RATIO, "sigma": MODULATION_SIGMA}})


def _nearest(points: list[dict]) -> dict | None:
    """The point with the largest sigma (None as -inf), the first on a tie: the nearest miss."""
    best, best_key = None, -math.inf
    for point in points:
        key = point["sigma"] if point["sigma"] is not None else -math.inf
        if best is None or key > best_key:
            best, best_key = point, key
    if best is None:
        return None
    return {"q": list(best["q"]), "ratio": best["ratio"], "sigma": best["sigma"], "because": best["because"]}


def _modulation_finding(answer, because: str, decided=(), nearest=None) -> dict:
    return {"answer": answer, "because": because, "decidedBy": [list(point["q"]) for point in decided],
            "nearest": nearest}


def _element_verdict(entry: dict) -> dict:
    zone = entry["zoneBoundary"]
    if not zone:
        cell_doubling = _modulation_finding(None, "untestable")
    else:
        ordered = [point for point in zone if point["ordered"] is True]
        if ordered:
            cell_doubling = _modulation_finding(True, "ordered_point", ordered)
        elif all(point["because"] == "undefined_statistic" for point in zone):
            cell_doubling = _modulation_finding(None, "undefined_statistic", nearest=_nearest(zone))
        else:
            cell_doubling = _modulation_finding(False, "no_ordered_point", nearest=_nearest(zone))

    general = [mode for mode in entry["strongest"] if mode["kind"] == "general"]
    ordered = [mode for mode in general if mode["ordered"] is True]
    if not general:
        general_order = _modulation_finding(None, "no_general_mode")
    elif ordered:
        general_order = _modulation_finding(True, "ordered_mode", ordered)
    elif all(mode["because"] == "undefined_statistic" for mode in general):
        general_order = _modulation_finding(None, "undefined_statistic", nearest=_nearest(general))
    else:
        general_order = _modulation_finding(False, "no_ordered_mode", nearest=_nearest(general))

    first = entry["strongest"][0] if entry["strongest"] else None
    if first is None or first["ratio"] is None:
        long_wave = _modulation_finding(None, "undefined_statistic")
    elif first["kind"] != "long_wave":
        long_wave = _modulation_finding(False, "strongest_not_long_wave")
    elif first["ratio"] >= LONG_WAVE_RATIO:
        long_wave = _modulation_finding(True, "long_wave_dominates", [first])
    else:
        long_wave = _modulation_finding(False, "below_ratio")
    return {"cellDoubling": cell_doubling, "generalOrder": general_order, "longWave": long_wave}


def _untested_verdict() -> dict:
    """An element with no displacement to analyse (zero spread everywhere)."""
    return {name: _modulation_finding(None, "no_displacements") for name in ("cellDoubling", "generalOrder", "longWave")}
