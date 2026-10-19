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
against independent cells (family-wise, Sidak-corrected over the points or
modes of every element the configuration has, whichever the call names), its
``kind`` ('zone_boundary', 'long_wave', 'general'), whether it is ``ordered``
(ratio >= ``MODULATION_RATIO`` and sigma >= ``MODULATION_SIGMA``),
``because`` and its ``bar``; each element a ``verdict`` -- ``cellDoubling``,
``generalOrder``, ``longWave``, ``diffuseOrder`` -- with the q that decided it
and the nearest miss (``diffuseOrder``: its evidence); the payload the
``family`` the corrections used (its sizes and elements) and the ``rules``.
The browser port is ``workers/modulations.js`` (a separable DFT, the same
numbers to floating-point noise; the parity fixture pins the verdicts and
polarisations too).

The null is the configuration's own cells placed at random (finding F1,
2026-10-19; :func:`independent_cell_null`): every correlation inside a cell is
kept -- the ligands of a Jahn-Teller octahedron, the vertices of a cluster,
which share one variable per cell -- and every correlation between cells is
destroyed. Its cross-spectral matrix is white, E[U(q) U(q)^H] =
N/(N - 1) sum_c y(c) y(c)^T at every q != 0, with y(c) the displacements of
all the element's sites in cell c, so the element's power is a weighted sum of
chi^2_1 (a zone-boundary point, a real Fourier sum) or of chi^2_2 / 2 (a
general q) with the eigenvalues of that matrix, normalised as the power, as
weights; its tail is a saddlepoint approximation (:func:`null_tail`). Until
2026-10-19 the null took the element's m sites as independent (P >= chi^2_m /
m): sites that share a cell-level variable carry fewer independent components,
and it found order in 3 of 12 seeds of a two-feature configuration, 2 of 16
Jahn-Teller and 1 of 10 octahedral ones; for independent sites it was
conservative (one component per site where an isotropic site has three).

``diffuseOrder`` (plan item 4.4; :func:`diffuse_order`) asks whether the power
is spread over many q in the pattern that correlation between neighbouring
cells makes -- short-range order, which no single mode carries: the normal
scores of the ranks of the power over the compared modes (no long wave, no
zone-boundary point, no mode already found ordered) regressed on the 13
neighbour-cell lags cos(2 pi q.R), against the null of independent cells,
under which the power at every such q is identically distributed. Parseval
fixes the mean power over q != 0, so the test is about how the power is spread,
never its total.
"""

from __future__ import annotations

import math
from itertools import product

import numpy as np
from scipy.special import gammaincc, ndtri

from rmc_toolkits.clusters import evidence_clears
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
# and naming one element of three took a 2.98 sigma point to 3.30 sigma). The
# null is the independent-cell null (finding F1, 2026-10-19): the element's
# power at q is U^H W U with U's covariance N/(N - 1) sum_c y(c) y(c)^T (its
# cells placed at random), a weighted sum of chi2_1 at a half-integer q (a
# real Fourier sum) and of chi2_2 / 2 at a general q (a complex one), whose
# weights are the eigenvalues of that covariance normalised as the power
# (independent_cell_null); the tail is the Lugannani-Rice saddlepoint
# (null_tail), within +1 to +7 % of the exact tail from p = 1e-2 to 1e-9 on
# every weight pattern checked, never below it. Until F1 the null took the
# sites as independent, P >= chi2_m / m (and chi2_2m / 2m at a general q):
# anti-conservative where sites share a cell-level variable (a Jahn-Teller
# axis, a cluster's breathing), conservative where they do not. The general
# modes are corrected over the compared modes. The earlier rule (top 1 % of
# the modes, or rank 1) had no correction: the validation cases found "order"
# in 6 of 12 configurations without any.

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

# --- the independent-cell null (F1) ---------------------------------------------

#: A null weight (an eigenvalue of the normalised cell covariance) below this
#: fraction of the largest is a direction no site moves in (or rounding): left out.
NULL_WEIGHT_RTOL = 1e-12
#: Weights within this relative spread of each other are one weight: the null
#: is then exactly a gamma (as for one isotropic site), no saddlepoint needed.
NULL_EQUAL_RTOL = 1e-12
#: Bisection steps for the saddlepoint equation K'(t) = x (the bracket shrinks
#: to adjacent doubles long before; both engines run the same steps).
SADDLEPOINT_STEPS = 200
#: Within this of the null mean (|w|, the signed root deviance) the
#: Lugannani-Rice formula cancels: its limit there is used.
SADDLEPOINT_W_MIN = 1e-6

# --- diffuse order between cells (plan item 4.4) ---------------------------------

#: The lags R (in cells) whose displacement correlation the diffuse test reads:
#: the 26 cells around a cell (face, edge and corner neighbours), one of each
#: +R/-R pair. The power of correlation at lag R is cos(2 pi q.R) across q.
DIFFUSE_LAGS = (
    (1, 0, 0), (0, 1, 0), (0, 0, 1),
    (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1),
    (1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1),
)
#: Family-wise one-sided significance (sigma) diffuse order needs (Sidak over
#: every element of the configuration, as the other tests).
DIFFUSE_SIGMA = 3
#: The effect size diffuse order needs: the strongest neighbour-cell
#: displacement correlation |Gamma(R)| (a +-2 |Gamma| modulation of the power
#: across q). With thousands of cells a correlation of 0.02 is significant and
#: still means nothing physically.
DIFFUSE_CORRELATION = 0.1
#: The test needs at least this many compared modes per term of the
#: regression (the constant and every lag the grid resolves): fewer, and the
#: chi^2 reference is not reached (a 4 x 4 x 4 grid has 15 such modes).
DIFFUSE_MODES_PER_TERM = 2
#: A lag whose template is, within this relative norm, a combination of the
#: constant and the lags before it adds nothing on this grid (an axis of two
#: cells makes +R and -R one lag): left out, and the degrees of freedom with it.
DIFFUSE_RANK_RTOL = 1e-9

MODULATION_RULES = {
    "ratio": MODULATION_RATIO,
    "sigma": MODULATION_SIGMA,
    "longWaveRatio": LONG_WAVE_RATIO,
    "longWaveQ": LONG_WAVE_Q,
    "sigmaMax": SIGMA_MAX,
    "diffuseSigma": DIFFUSE_SIGMA,
    "diffuseCorrelation": DIFFUSE_CORRELATION,
    "diffuseLags": [list(lag) for lag in DIFFUSE_LAGS],
    "diffuseModesPerTerm": DIFFUSE_MODES_PER_TERM,
}
#: Every ``because`` a point and each element finding can give (closed sets;
#: the parity fixture pins the strings).
MODULATION_BECAUSE = {
    "point": ("ordered", "below_ratio", "below_sigma", "undefined_statistic", "zone_boundary", "long_wave"),
    "cellDoubling": ("ordered_point", "no_ordered_point", "undefined_statistic", "untestable", "no_displacements"),
    "generalOrder": ("ordered_mode", "no_ordered_mode", "undefined_statistic", "no_general_mode", "no_displacements"),
    "longWave": ("long_wave_dominates", "below_ratio", "strongest_not_long_wave", "undefined_statistic",
                 "no_displacements"),
    "diffuseOrder": ("correlated_in_phase", "correlated_antiphase", "below_correlation", "below_sigma",
                     "undefined_statistic", "too_few_modes", "no_displacements"),
}
MODE_KINDS = ("zone_boundary", "long_wave", "general")
#: The phase of the strongest neighbour-cell correlation: in phase (Gamma > 0,
#: the power toward the zone centre) or in antiphase (Gamma < 0, toward the
#: zone boundary).
DIFFUSE_PHASES = ("in_phase", "antiphase")


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


def _null_weights(weights) -> np.ndarray:
    """The null's weights as an array: an int m stands for m equal weights 1/m
    (the bound P >= chi2_m / m the engine used until F1)."""
    if isinstance(weights, (int, np.integer)) and not isinstance(weights, bool):
        return np.full(int(weights), 1.0 / int(weights)) if weights > 0 else np.zeros(0)
    return np.asarray(weights, dtype=float).ravel()


def _normal_upper(w: float) -> float:
    """P(Z > w) for a standard normal Z, accurate in both tails (as the port computes it)."""
    half = 0.5 * float(gammaincc(0.5, 0.5 * w * w))
    return half if w >= 0 else 1.0 - half


def _saddlepoint_tail(x: float, b: np.ndarray, a: float) -> float:
    """Lugannani-Rice upper tail of X with cumulant generating function
    K(t) = -a sum_k log(1 - b_k t) (b_k > 0): a = 1/2, b = 2 lambda for
    sum lambda_k chi2_1; a = 1, b = lambda for sum lambda_k chi2_2 / 2.

    The saddlepoint t solves K'(t) = x by bisection (SADDLEPOINT_STEPS steps,
    as the port); w = sign(t) sqrt(2 (t x - K(t))), u = t sqrt(K''(t)) and
    P(X > x) = 1 - Phi(w) + phi(w) (1/u - 1/w). Within SADDLEPOINT_W_MIN of
    the mean the formula's limit, 1/2 - kappa_3 / (6 sqrt(2 pi) kappa_2^3/2).
    """
    mean = a * float(np.sum(b))
    if x > mean:
        lo, hi = 0.0, 1.0 / float(np.max(b))
    else:
        # K'(t) < a n / |t| for t < 0: at t = -a n / x it is below x.
        lo, hi = -a * b.size / x, 0.0
    for _ in range(SADDLEPOINT_STEPS):
        mid = 0.5 * (lo + hi)
        if not lo < mid < hi:
            break  # the bracket is two adjacent doubles
        if a * float(np.sum(b / (1.0 - b * mid))) < x:
            lo = mid
        else:
            hi = mid
    # The root lies in [lo, hi]; lo is never past the pole 1 / max(b).
    t = lo if x > mean else 0.5 * (lo + hi)
    cgf = -a * float(np.sum(np.log1p(-b * t)))
    w = math.copysign(math.sqrt(max(2.0 * (t * x - cgf), 0.0)), t)
    if abs(w) < SADDLEPOINT_W_MIN:
        kappa2 = a * float(np.sum(b * b))
        kappa3 = 2.0 * a * float(np.sum(b * b * b))
        return 0.5 - kappa3 / (6.0 * math.sqrt(2.0 * math.pi) * kappa2 ** 1.5)
    u = t * math.sqrt(a * float(np.sum((b / (1.0 - b * t)) ** 2)))
    tail = _normal_upper(w) + math.exp(-0.5 * w * w) / math.sqrt(2.0 * math.pi) * (1.0 / u - 1.0 / w)
    return min(max(tail, 0.0), 1.0)


def null_tail(power, weights, real: bool) -> float | None:
    """P(P_null >= power) under the independent-cell null.

    ``weights`` are the null's weights (independent_cell_null: the
    eigenvalues of the element's normalised cell covariance), or an int m for
    m equal weights 1/m. ``real``: a zone-boundary point, where the Fourier
    sum is real and P = sum lambda_k chi2_1; else a general q, P = sum
    lambda_k chi2_2 / 2. One weight value: the exact gamma. Several: the
    Lugannani-Rice saddlepoint. None where undefined (power <= 0, no weight).
    """
    lam = _null_weights(weights)
    lam = lam[lam > 0]
    if power is None or not (power > 0) or not lam.size:
        return None
    a, scale = (0.5, 2.0) if real else (1.0, 1.0)
    if float(lam.max() - lam.min()) <= NULL_EQUAL_RTOL * float(lam.max()):
        # sum of n equal terms: lambda chi2_n (real) or lambda Gamma(n, 1) (complex).
        return float(gammaincc(a * lam.size, float(power) / (scale * float(lam.mean()))))
    return _saddlepoint_tail(float(power), scale * lam, a)


def zone_boundary_sigma(power, weights, family: int) -> float | None:
    """The family-wise significance (one-sided sigma) of an element's
    zone-boundary power ``power`` (normalised: its mean over q != 0 is
    N/(N - 1)) against the independent-cell null with ``weights``
    (:func:`null_tail`, real), among ``family`` points tested. None where
    undefined (power <= 0) or -inf; +inf is SIGMA_MAX."""
    p = null_tail(power, weights, True)
    if p is None:
        return None
    return _reported(_sidak_sigma(p, family))


def general_mode_sigma(power, weights, family: int) -> float | None:
    """The same at a general wavevector, where the Fourier sum is complex;
    ``family`` = the compared modes. Order there is correlation between cells
    that does not double the cell: incommensurate, or short-range (the
    2026-10-14 review's L4)."""
    if not (family >= 1):
        return None
    p = null_tail(power, weights, False)
    if p is None:
        return None
    return _reported(_sidak_sigma(p, family))


def independent_cell_null(sites: SiteDisplacements, name: str) -> np.ndarray:
    """The weights of the independent-cell null of one element's power (F1).

    With y(c) the displacements of the element's m moving sites in cell c (3m
    components; a missing copy is a zero, as on the transform's grid) and the
    power P(q) = U(q)^H W U(q), W = diag(1 / (m n_s tr C_s)), placing the
    cells at random gives U(q) the covariance N/(N - 1) sum_c y(c) y(c)^T at
    every q != 0 (U(0) = 0: each site's mean is subtracted). The weights are
    the eigenvalues of W^1/2 (that covariance) W^1/2, largest first; their sum
    is the mean power over q != 0, N/(N - 1). Weights below NULL_WEIGHT_RTOL
    of the largest are left out. Empty when no site of ``name`` moves.
    """
    supercell = np.asarray(sites.supercell, dtype=int)
    n_q = int(np.prod(supercell))
    cells = np.mod(np.asarray(sites.cell_indices, dtype=int), supercell)
    columns, weights = [], []
    for s, label in enumerate(sites.elements):
        if label != name:
            continue
        rows = np.flatnonzero(sites.site_index == s)
        u = sites.displacements[rows]
        variance = _site_variance(u)
        if not variance > 0:
            continue
        column = np.zeros((n_q, 3))
        column[np.ravel_multi_index(tuple(cells[rows].T), tuple(supercell))] = u
        columns.append(column)
        weights.append(1.0 / (rows.size * variance))
    if not columns:
        return np.zeros(0)
    scale = np.sqrt(np.repeat(np.asarray(weights) / len(columns), 3))
    z = np.concatenate(columns, axis=1) * scale
    values = np.linalg.eigvalsh((n_q / (n_q - 1.0)) * (z.T @ z))[::-1]
    return values[values > values[0] * NULL_WEIGHT_RTOL]


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
    # The diffuse test's modes: the compared modes that are not zone-boundary
    # points (tested on their own, and real where the rest are complex).
    on_boundary = np.all((grid_index == 0) | (2 * grid_index == supercell), axis=1)
    diffuse_modes = np.flatnonzero(compared & ~on_boundary)
    spectra = {}
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
        strongest_keys = []
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
            strongest_keys.append(int(np.ravel_multi_index(key, tuple(supercell))))
            if len(strongest) == top:
                break
        weights = independent_cell_null(sites, name)
        result["elements"].append({
            "element": name,
            "sites": used,
            "missingCopies": missing,
            "medianPower": median,
            "comparedModes": int(reference.size),
            # The independent-cell null's effective number of independent
            # components, (sum lambda)^2 / sum lambda^2 (3m for m independent
            # isotropic sites; fewer where sites share a cell-level variable).
            "nullDof": float(weights.sum() ** 2 / np.sum(weights * weights)),
            "zoneBoundary": zone,
            "strongest": strongest,
        })
        spectra[name] = {"weights": weights, "power": flat, "strongest": strongest_keys}
    _decide(result, family_elements, len(boundary), int(compared.sum()), spectra,
            {"modes": diffuse_modes, "q": q_frac[diffuse_modes]})
    return result


def _average_ranks(values: np.ndarray) -> np.ndarray:
    """Ranks 1..K of ``values`` in ascending order; a run of values within
    PERCENTILE_RTOL of its smallest shares the run's mean rank (as the port)."""
    order = np.argsort(values, kind="stable")
    ranks = np.empty(values.size)
    start = 0
    while start < values.size:
        first = float(values[order[start]])
        end = start + 1
        while end < values.size and float(values[order[end]]) <= first + PERCENTILE_RTOL * abs(first):
            end += 1
        ranks[order[start:end]] = 0.5 * (start + 1 + end)
        start = end
    return ranks


def _orthonormal_columns(columns: list[np.ndarray]) -> tuple[list[np.ndarray], list[int]]:
    """Modified Gram-Schmidt over ``columns`` in order: an orthonormal basis of
    their span and the indices of the columns kept. A column whose part
    orthogonal to the earlier ones is within DIFFUSE_RANK_RTOL of its own
    norm adds nothing and is left out (as the port decides it)."""
    basis, kept = [], []
    for index, column in enumerate(columns):
        vector = np.array(column, dtype=float)
        norm = math.sqrt(float(vector @ vector))
        for unit in basis:
            vector = vector - float(unit @ vector) * unit
        rest = math.sqrt(float(vector @ vector))
        if norm > 0 and rest > DIFFUSE_RANK_RTOL * norm:
            basis.append(vector / rest)
            kept.append(index)
    return basis, kept


def _solve_dense(matrix: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """A small dense system by Gaussian elimination with partial pivoting (the port's)."""
    a = np.array(matrix, dtype=float)
    b = np.array(rhs, dtype=float)
    n = b.size
    for column in range(n):
        pivot = column + int(np.argmax(np.abs(a[column:, column])))
        if pivot != column:
            a[[column, pivot]] = a[[pivot, column]]
            b[[column, pivot]] = b[[pivot, column]]
        for row in range(column + 1, n):
            factor = a[row, column] / a[column, column]
            a[row, column:] -= factor * a[column, column:]
            b[row] -= factor * b[column]
    x = np.zeros(n)
    for row in range(n - 1, -1, -1):
        x[row] = (b[row] - float(a[row, row + 1:] @ x[row + 1:])) / a[row, row]
    return x


def diffuse_order(power, q, family: int) -> dict:
    """The diffuse test (plan item 4.4): is the power spread over many q in
    the pattern correlation between neighbouring cells makes?

    ``power`` holds the element's power at the K modes compared (one per
    +q/-q pair, no long wave, no zone-boundary point, none found ordered),
    ``q`` their fractional wavevectors. Correlation Gamma(R) between a site's
    displacements in cells R apart adds 2 Gamma(R) cos(2 pi q.R) to the power
    over its mean (the power is the Fourier transform of that correlation).
    Under independent cells the power at every such q is identically
    distributed (its cross-spectral matrix is white), so the ranks of the
    power over q are exchangeable whatever the cells hold. The statistic
    regresses the normal scores a of those ranks on the lags' templates
    cos(2 pi q.R), R in DIFFUSE_LAGS (with the constant; a lag the grid does
    not resolve is left out, and its degree of freedom):
    chi^2 = (K - 1) |P a|^2 / |a|^2, P the projection onto the templates
    orthogonal to the constant, which is chi^2 with r degrees of freedom under
    the null (r the lags kept). Ranks bound what any one q can add, so no
    single mode makes it. The family is every element (Sidak, as F3).

    Parseval fixes the mean power over q != 0 (N/(N - 1)): this is about how
    the power is spread, never its total. Long waves stay out: RMC's
    box-scale waves are not short-range order (they are the longWave
    finding).

    Returns ``modes`` (K), ``df`` (r, None when the grid gives fewer than
    DIFFUSE_MODES_PER_TERM modes per term or no lag), ``chiSquare``,
    ``sigma`` (family-wise), ``correlations`` (Gamma per DIFFUSE_LAGS from
    the least-squares fit of the power itself, beta_R / (2 mean power): the
    lag's relative modulation of the power over these modes, halved; None for
    a lag left out), and the strongest: ``lag``, ``correlation``, ``phase``
    ('in_phase' for Gamma > 0, the power toward the zone centre; 'antiphase'
    for Gamma < 0, toward the zone boundary).
    """
    power = np.asarray(power, dtype=float).ravel()
    q = np.asarray(q, dtype=float).reshape(-1, 3)
    modes = int(power.size)
    out = {"modes": modes, "df": None, "chiSquare": None, "sigma": None,
           "correlations": [None] * len(DIFFUSE_LAGS), "lag": None, "correlation": None, "phase": None}
    if modes == 0:
        return out
    columns = [np.ones(modes)] + [np.cos(2.0 * np.pi * (q @ np.asarray(lag, dtype=float))) for lag in DIFFUSE_LAGS]
    basis, kept = _orthonormal_columns(columns)
    lags = [index - 1 for index in kept if index > 0]
    if not lags or kept[0] != 0 or modes < DIFFUSE_MODES_PER_TERM * (len(lags) + 1):
        return out
    out["df"] = len(lags)
    design = np.column_stack([columns[index] for index in kept])
    beta = _solve_dense(design.T @ design, design.T @ power)
    mean = float(power.mean())
    if mean > 0:
        for position, lag in enumerate(lags, start=1):
            out["correlations"][lag] = float(beta[position] / (2.0 * mean))
    scores = ndtri((_average_ranks(power) - 0.5) / modes)
    scores = scores - scores.mean()
    total = float(scores @ scores)
    if not total > 0 or not mean > 0:
        return out
    projection = sum(float(unit @ scores) ** 2 for unit in basis[1:])
    chi_square = (modes - 1) * projection / total
    out["chiSquare"] = float(chi_square)
    out["sigma"] = _reported(_sidak_sigma(float(gammaincc(len(lags) / 2.0, chi_square / 2.0)), family))
    # The strongest correlation: the largest |Gamma|, the first lag on a tie.
    magnitude = np.array([abs(out["correlations"][lag]) for lag in lags])
    strongest = lags[int(np.flatnonzero(magnitude >= magnitude.max() * (1.0 - PERCENTILE_RTOL))[0])]
    out["lag"] = list(DIFFUSE_LAGS[strongest])
    out["correlation"] = out["correlations"][strongest]
    out["phase"] = "in_phase" if out["correlation"] > 0 else "antiphase"
    return out


def _diffuse_finding(diffuse: dict) -> dict:
    """diffuseOrder: answer, because and one Evidence record (clusters.py's
    form): observed = the strongest |Gamma(R)|, null 0, z = the family-wise
    sigma, change = observed - null, bar {z: DIFFUSE_SIGMA, change:
    DIFFUSE_CORRELATION} cleared from above; clears recomputable
    (evidence_clears). The answer is the record's clears."""
    if diffuse["df"] is None:
        return {"answer": None, "because": "too_few_modes", "evidence": []}
    magnitude = abs(diffuse["correlation"]) if diffuse["correlation"] is not None else None
    record = {
        "statistic": "neighbourCellCorrelation", "against": "independentCells",
        "lag": diffuse["lag"], "correlation": diffuse["correlation"], "phase": diffuse["phase"],
        "chiSquare": diffuse["chiSquare"], "df": diffuse["df"], "modes": diffuse["modes"],
        "observed": magnitude, "null": 0.0, "z": diffuse["sigma"], "change": magnitude, "changeKind": "difference",
        "bar": {"z": DIFFUSE_SIGMA, "change": DIFFUSE_CORRELATION}, "side": "above",
    }
    record["clears"] = evidence_clears(record)
    if diffuse["chiSquare"] is None or magnitude is None:
        return {"answer": None, "because": "undefined_statistic", "evidence": [record]}
    if record["clears"]:
        because = "correlated_in_phase" if diffuse["phase"] == "in_phase" else "correlated_antiphase"
    elif diffuse["sigma"] is None or diffuse["sigma"] < DIFFUSE_SIGMA:
        because = "below_sigma"
    else:
        because = "below_correlation"
    return {"answer": record["clears"], "because": because, "evidence": [record]}


def _decide(result: dict, family_elements: list[str], points: int, modes: int, spectra: dict, diffuse: dict) -> None:
    """Annotate every tested point with its sigma, ordered, because and bar,
    and give each element its verdict, over the family: the ``points``
    zone-boundary points and ``modes`` compared modes of each of
    ``family_elements`` -- every element the configuration has with
    displacements, whichever the call reports (F3) -- and one diffuse test
    per element. The sizes are the same for every element, so this is what a
    call without ``element`` counted. ``spectra`` holds each tested element's
    null weights, power and strongest modes' grid indices; ``diffuse`` the
    diffuse test's modes and their wavevectors."""
    tested = [entry for entry in result["elements"] if entry["sites"]]
    family = {
        "zoneBoundaryPoints": points * len(family_elements),
        "comparedModes": modes * len(family_elements),
        "diffuse": len(family_elements),
        "elements": list(family_elements),
    }
    for entry in tested:
        spectrum = spectra[entry["element"]]
        for point in entry["zoneBoundary"]:
            point["kind"] = "zone_boundary"
            _test(point, "zone", entry, family["zoneBoundaryPoints"], spectrum["weights"])
        for mode in entry["strongest"]:
            if mode["kind"] == "general":
                _test(mode, "general", entry, family["comparedModes"], spectrum["weights"])
            else:
                mode.update({"sigma": None, "ordered": None, "because": mode["kind"], "bar": None})
        # The diffuse test asks what the modes found ordered leave.
        found = {key for key, mode in zip(spectrum["strongest"], entry["strongest"])
                 if mode["kind"] == "general" and mode["ordered"] is True}
        keep = np.array([int(key) not in found for key in diffuse["modes"]], dtype=bool)
        entry["diffuse"] = diffuse_order(spectrum["power"][diffuse["modes"][keep]], diffuse["q"][keep], family["diffuse"])
        entry["verdict"] = _element_verdict(entry)
    result["family"] = family
    result["rules"] = dict(MODULATION_RULES)


def _test(point: dict, test: str, entry: dict, family: int, weights: np.ndarray) -> None:
    """One point against its bar and the element's independent-cell null
    (``weights``). The power is ratio x median, as the tool formed it: sigma
    far below the background is ill-conditioned (1 - (1 - p)^F rounds), and
    this keeps it bit for bit the value the tool reported."""
    ratio = point["ratio"]
    power = ratio * entry["medianPower"] if ratio is not None else None
    defined = power is not None and power > 0 and len(weights) > 0 and (test == "zone" or family >= 1)
    sigma = None
    if not defined:
        ordered, because = None, "undefined_statistic"
    else:
        sigma = zone_boundary_sigma(power, weights, family) if test == "zone" else general_mode_sigma(power, weights, family)
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
    return {"cellDoubling": cell_doubling, "generalOrder": general_order, "longWave": long_wave,
            "diffuseOrder": _diffuse_finding(entry["diffuse"])}


def _untested_verdict() -> dict:
    """An element with no displacement to analyse (zero spread everywhere)."""
    verdict = {name: _modulation_finding(None, "no_displacements") for name in ("cellDoubling", "generalOrder", "longWave")}
    verdict["diffuseOrder"] = {"answer": None, "because": "no_displacements", "evidence": []}
    return verdict
