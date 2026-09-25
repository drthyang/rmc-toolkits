# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Keen-convention total-scattering functions and Fourier transforms.

Golden reference: D. A. Keen, J. Appl. Cryst. 34, 172-177 (2001),
doi:10.1107/S0021889800019993 — equation numbers below refer to that paper.

Conventions used throughout:

- ``S(Q)``      normalized total-scattering structure factor (Eq. 19); -> 1 at high Q.
- ``F(Q)``      ``Q [S(Q) - 1]`` (the pystog/PDF community "F(Q)"; Keen calls his
                barns-scale function F(Q) — here that one is ``FK(Q)``).
- ``FK(Q)``     Keen's F(Q) (Eq. 9), barns: ``<b>^2 [S(Q) - 1]`` with
                ``<b>^2 = (sum_i c_i b_i)^2``.
- ``g(r)``      pair distribution function; 0 below the closest approach, -> 1.
- ``G_PDF(r)``  PDFFIT-style ``4 pi rho0 r [g(r) - 1]`` (Keen Eq. 43/44 lineage).
- ``GK(r)``     Keen's G(r) (Eq. 10), barns: ``<b>^2 [g(r) - 1]``; flat ``-<b>^2``
                below the closest approach (Eq. 15).
- ``D(r)``      ``4 pi rho0 r GK(r)`` (Eq. 29); slope ``-4 pi rho0 <b>^2 r`` at low r.

The Fourier pair (validated against pystog and a classic Fortran stog run):

- ``G_PDF(r) = (2/pi) * integral F(Q) sin(Q r) dQ``
- ``F(Q)     =          integral G_PDF(r) sin(Q r) dr``

Both integrals are evaluated with the trapezoid rule on the supplied grids, which
matches the pystog ``Transformer`` discretization (and the Fortran stog to ~6e-4 rms).
"""

from __future__ import annotations

import numpy as np

try:  # numpy >= 2.0
    _trapezoid = np.trapezoid
except AttributeError:  # pragma: no cover - numpy < 2.0
    _trapezoid = np.trapz

_SINE_CHUNK = 512


def _increasing_grid(x: np.ndarray, name: str) -> np.ndarray:
    """Return ``x`` as a 1-D float array, or raise unless strictly increasing.

    The trapezoid rule integrates with signed panel widths ``x[i] - x[i-1]``: a
    descending grid silently negates every integral (and swaps Qmin/Qmax in the
    low-Q correction and the Lorch window), a repeated point gives a zero-width
    panel, and a non-monotone grid folds panels back over each other. Grids with
    fewer than two points (an empty filter section) integrate to 0 and pass.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        raise ValueError(f"the {name} grid must be one-dimensional, got shape {x.shape}")
    if x.size > 1 and not (np.all(np.isfinite(x)) and np.all(np.diff(x) > 0)):
        raise ValueError(
            f"the {name} grid must be strictly increasing (and finite); sort it "
            "and merge or remove duplicate points first"
        )
    return x


# ---------------------------------------------------------------------------
# Algebraic conversions (pure, grid-preserving)
# ---------------------------------------------------------------------------

def sq_to_fq(q: np.ndarray, sq: np.ndarray) -> np.ndarray:
    """``F(Q) = Q [S(Q) - 1]``."""
    return np.asarray(q, dtype=float) * (np.asarray(sq, dtype=float) - 1.0)


def fq_to_sq(q: np.ndarray, fq: np.ndarray) -> np.ndarray:
    """``S(Q) = F(Q)/Q + 1`` (requires Q > 0)."""
    return np.asarray(fq, dtype=float) / np.asarray(q, dtype=float) + 1.0


def sq_to_fk(sq: np.ndarray, b_avg_sq: float) -> np.ndarray:
    """Keen Eq. 19 inverted: ``FK(Q) = <b>^2 [S(Q) - 1]`` (barns)."""
    return float(b_avg_sq) * (np.asarray(sq, dtype=float) - 1.0)


def fk_to_sq(fk: np.ndarray, b_avg_sq: float) -> np.ndarray:
    """``S(Q) = FK(Q)/<b>^2 + 1``."""
    return np.asarray(fk, dtype=float) / float(b_avg_sq) + 1.0


def g_to_gpdf(r: np.ndarray, g: np.ndarray, rho0: float) -> np.ndarray:
    """``G_PDF(r) = 4 pi rho0 r [g(r) - 1]``."""
    return 4.0 * np.pi * float(rho0) * np.asarray(r, dtype=float) * (
        np.asarray(g, dtype=float) - 1.0
    )


def gpdf_to_g(r: np.ndarray, gpdf: np.ndarray, rho0: float) -> np.ndarray:
    """``g(r) = G_PDF(r) / (4 pi rho0 r) + 1`` (requires r > 0)."""
    return np.asarray(gpdf, dtype=float) / (
        4.0 * np.pi * float(rho0) * np.asarray(r, dtype=float)
    ) + 1.0


def g_to_gk(g: np.ndarray, b_avg_sq: float) -> np.ndarray:
    """Keen Eq. 16 inverted: ``GK(r) = <b>^2 [g(r) - 1]`` (barns)."""
    return float(b_avg_sq) * (np.asarray(g, dtype=float) - 1.0)


def gk_to_g(gk: np.ndarray, b_avg_sq: float) -> np.ndarray:
    """``g(r) = GK(r)/<b>^2 + 1``."""
    return np.asarray(gk, dtype=float) / float(b_avg_sq) + 1.0


def gk_to_dr(r: np.ndarray, gk: np.ndarray, rho0: float) -> np.ndarray:
    """Keen Eq. 29: ``D(r) = 4 pi rho0 r GK(r)``."""
    return 4.0 * np.pi * float(rho0) * np.asarray(r, dtype=float) * np.asarray(
        gk, dtype=float
    )


def density_line(r: np.ndarray, rho0: float, b_avg_sq: float) -> np.ndarray:
    """Theoretical low-r ``D(r) = -4 pi rho0 <b>^2 r`` (Keen Eqs. 29 + 15)."""
    return -4.0 * np.pi * float(rho0) * float(b_avg_sq) * np.asarray(r, dtype=float)


def _sinc(v: np.ndarray) -> np.ndarray:
    """``sin(v)/v`` with the limit 1 at v = 0 (unnormalized sinc)."""
    v = np.asarray(v, dtype=float)
    return np.divide(np.sin(v), v, out=np.ones_like(v), where=v != 0)


def _sinc_head(v: np.ndarray) -> np.ndarray:
    """``(v sin v + cos v - 1)/v^2 = sinc(v) - sinc(v/2)^2 / 2``, cancellation-free (1/2 at 0)."""
    return _sinc(v) - 0.5 * _sinc(0.5 * np.asarray(v, dtype=float)) ** 2


def lorch_window(q: np.ndarray, qmax: float) -> np.ndarray:
    """Lorch modification ``sin(pi Q/Qmax) / (pi Q/Qmax)`` (pystog-compatible)."""
    x = np.pi * np.asarray(q, dtype=float) / float(qmax)
    return np.divide(np.sin(x), x, out=np.ones_like(x), where=x != 0)


# ---------------------------------------------------------------------------
# Fourier transforms
# ---------------------------------------------------------------------------

def sine_transform(x: np.ndarray, y: np.ndarray, xout: np.ndarray) -> np.ndarray:
    """Trapezoid-rule ``integral y(x) sin(x * xout) dx`` for each point of ``xout``.

    ``x`` must be strictly increasing (``ValueError`` otherwise); ``xout`` may
    have any shape (a scalar included) and the result has the same shape.
    Chunked over the output grid to bound the kernel-matrix memory.
    """
    x = _increasing_grid(x, "integration (x)")
    y = np.asarray(y, dtype=float)
    xout = np.asarray(xout, dtype=float)
    flat = xout.ravel()
    out = np.empty_like(flat)
    for start in range(0, flat.size, _SINE_CHUNK):
        chunk = flat[start : start + _SINE_CHUNK]
        kernel = y[np.newaxis, :] * np.sin(np.outer(chunk, x))
        out[start : start + chunk.size] = _trapezoid(kernel, x=x, axis=1)
    return out.reshape(xout.shape)


def fq_to_gpdf(
    q: np.ndarray,
    fq: np.ndarray,
    r: np.ndarray,
    *,
    lorch: bool = False,
    low_q_correction: bool = False,
    s0_target: float = 0.0,
) -> np.ndarray:
    """Forward transform ``G_PDF(r) = (2/pi) integral F(Q) sin(Q r) dQ``.

    ``lorch`` applies the Lorch window to ``F(Q)`` before transforming.
    ``low_q_correction`` adds the analytic correction for the omitted
    ``[0, Qmin]`` range, extrapolating to ``S(0) = s0_target`` (see
    :func:`omitted_low_q_correction`).
    """
    q = _increasing_grid(q, "Q")
    fq = np.asarray(fq, dtype=float)
    weighted = fq * lorch_window(q, q[-1]) if lorch else fq
    gpdf = (2.0 / np.pi) * sine_transform(q, weighted, r)
    if low_q_correction:
        gpdf = gpdf + omitted_low_q_correction(
            q, fq, r, lorch=lorch, s0_target=s0_target
        )
    return gpdf


def gpdf_to_fq(r: np.ndarray, gpdf: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Backward transform ``F(Q) = integral G_PDF(r) sin(Q r) dr``."""
    return sine_transform(r, gpdf, q)


def low_q_correction_basis(
    q: np.ndarray,
    r: np.ndarray,
    *,
    lorch: bool = False,
    s0_target: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Affine basis of the omitted-low-Q correction on the ``r`` grid.

    The analytic correction models ``S(Q)`` on the unmeasured ``[0, Qmin]``
    as a straight line from ``S(0) = s0_target`` to the first measured point,
    and evaluates to

        ``correction(r) = coef_s0(r) * S(Qmin) - const(r)``

    (generalization of pystog ``Transformer._low_x_correction``, original
    author Jack Carpenter, in final ``G_PDF`` units; pystog's form is the
    ``s0_target = 0`` "solid-state limit" special case). The target enters the
    constant term only — ``const' = (1 - s0)*const + s0*coef`` — so the
    correction stays affine in a scale/offset applied to S(Q), which the
    auto-scaler's closed-form solve relies on. For compositions with strong
    negative-b cancellation (e.g. Mn3Sn, ``<b^2>/<b>^2`` ~ 13) the Keen
    Eq. 21 target ``S(0) = 1 - <b^2>/<b>^2`` is O(-10), and using it here
    instead of 0 removes an O(1) bias in the low-r transform.
    """
    q = _increasing_grid(q, "Q")
    r = np.asarray(r, dtype=float)
    if q[0] == 0:
        # Data starting at Q = 0 omit nothing — the [0, q[1]] panel is already
        # inside the trapezoid integral (pystog parity: F1 = F2 = 0 there).
        zero = np.zeros_like(r)
        return zero, zero
    q0 = q[0]

    v = q0 * r
    if lorch:
        # With the Lorch window M(Q) = sin(aQ)/(aQ), a = pi/qmax, the integrals
        # split into cos((r -/+ a) Q) terms. Written with v = Q0 (r -/+ a) as
        #   sin(v)/(r-a) = Q0 sinc(v),
        #   (v sin v + cos v - 1)/(r-a)^2 = Q0^2 [sinc(v) - sinc(v/2)^2 / 2]
        # (cos v - 1 = -2 sin^2(v/2)) they carry no subtraction of O(1) terms
        # and no removable singularity at r = a, so no patch is needed: the
        # previous (cos v - 1)/(r - a)^2 form lost every digit within ~1e-6 A of
        # pi/qmax (a sign-flipped coefficient at r = 0.11 for qmax = 28.56).
        a = np.pi / q[-1]
        vm = q0 * (r - a)
        vp = q0 * (r + a)
        f1 = q0**2 * (_sinc_head(vm) - _sinc_head(vp)) / (2.0 * a)
        f2 = q0 * (_sinc(vm) - _sinc(vp)) / (2.0 * a)
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            f1 = (2.0 * v * np.sin(v) - (v * v - 2.0) * np.cos(v) - 2.0) / r**3
            f2 = (np.sin(v) - v * np.cos(v)) / r**2
        f1 = np.where(r == 0.0, 0.0, f1)
        f2 = np.where(r == 0.0, 0.0, f2)
    coef = (2.0 / np.pi) * f1 / q0
    const = (2.0 / np.pi) * f2
    if s0_target != 0.0:
        const = (1.0 - s0_target) * const + s0_target * coef
    return coef, const


def omitted_low_q_correction(
    q: np.ndarray,
    fq: np.ndarray,
    r: np.ndarray,
    *,
    lorch: bool = False,
    s0_target: float = 0.0,
) -> np.ndarray:
    """Analytic ``G_PDF(r)`` contribution of the omitted ``[0, Qmin]`` range.

    See :func:`low_q_correction_basis`; here the basis is contracted with the
    measured ``S(Qmin)`` (pystog uses the first finite point).
    """
    q = np.asarray(q, dtype=float)
    fq = np.asarray(fq, dtype=float)
    if q[0] == 0:
        return np.zeros_like(np.asarray(r, dtype=float))
    s0 = fq[0] / q[0] + 1.0
    coef_s0, const = low_q_correction_basis(q, r, lorch=lorch, s0_target=s0_target)
    return coef_s0 * s0 - const


# ---------------------------------------------------------------------------
# Fourier filter (classic stog / pystog recipe)
# ---------------------------------------------------------------------------

def fourier_filter(
    q: np.ndarray,
    sq: np.ndarray,
    r: np.ndarray,
    *,
    rho0: float,
    cutoff: float,
    lorch: bool = False,
    low_q_correction: bool = False,
    s0_target: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Remove unphysical ``r < cutoff`` content from ``S(Q)``.

    Recipe (pystog ``FourierFilter``, validated to 6e-4 rms against a classic
    Fortran stog run): transform to ``g(r)``, take the raw ``g`` on
    ``[0, cutoff]``, back-transform ``4 pi rho0 r g_section`` to ``F(Q)``,
    subtract, and re-transform.

    Returns ``(sq_filtered, sq_ft, g_filtered)`` where ``sq_ft`` is the
    S(Q)-convention correction section (the classic stog ``ft.dat``), so that
    ``sq_filtered = sq - (sq_ft - 1)``.

    The ``r`` grid must be strictly increasing and non-negative; it may start
    at ``r = 0``: the section integrand ``4 pi rho0 r g`` is evaluated as
    ``G_PDF + 4 pi rho0 r`` (no division by r), and ``g_filtered`` at ``r = 0``
    is the continuous extension ``1 + G_PDF'(0) / (4 pi rho0)`` of the
    computed g(r) (:func:`gpdf_slope_at_zero`), not ``0/0``.
    """
    q = _increasing_grid(q, "Q")
    sq = np.asarray(sq, dtype=float)
    r = _increasing_grid(r, "r")
    if q[0] <= 0:
        raise ValueError(
            "fourier_filter requires a strictly positive Q grid (the S(Q) "
            "conversions divide by Q); crop Q <= 0 first"
        )
    if r.size and r[0] < 0:
        raise ValueError("fourier_filter requires a non-negative r grid")

    options = dict(lorch=lorch, low_q_correction=low_q_correction, s0_target=s0_target)
    fq = sq_to_fq(q, sq)
    gpdf = fq_to_gpdf(q, fq, r, **options)

    section = r <= float(cutoff)
    # pystog shifts the section by +1 and re-derives G_PDF, which reduces to
    # transforming 4 pi rho0 r * g(r) = G_PDF + 4 pi rho0 r over the section
    # (the identity needs no division by r, so r = 0 is harmless).
    fq_ft = gpdf_to_fq(
        r[section],
        gpdf[section] + 4.0 * np.pi * float(rho0) * r[section],
        q,
    )
    sq_ft = fq_to_sq(q, fq_ft)

    fq_filtered = fq - fq_ft
    sq_filtered = fq_to_sq(q, fq_filtered)
    gpdf_filtered = fq_to_gpdf(q, fq_filtered, r, **options)
    at_zero = r == 0
    with np.errstate(divide="ignore", invalid="ignore"):
        g_filtered = gpdf_to_g(r, gpdf_filtered, rho0)
    if np.any(at_zero):
        slope = gpdf_slope_at_zero(q, fq_filtered, **options)
        g_filtered[at_zero] = 1.0 + slope / (4.0 * np.pi * float(rho0))
    return sq_filtered, sq_ft, g_filtered


def gpdf_slope_at_zero(
    q: np.ndarray,
    fq: np.ndarray,
    *,
    lorch: bool = False,
    low_q_correction: bool = False,
    s0_target: float = 0.0,
) -> float:
    """``d G_PDF / dr`` at ``r = 0`` of :func:`fq_to_gpdf` (same options).

    ``G_PDF(r)`` is odd in r, so ``g = G_PDF/(4 pi rho0 r) + 1`` extends
    continuously to ``g(0) = 1 + G_PDF'(0)/(4 pi rho0)``. The derivative of the
    trapezoid sine transform is the trapezoid moment ``(2/pi) sum Q F(Q) M(Q)``
    (``sin(Qr)/r -> Q``); the omitted-low-Q correction contributes
    ``coef'(0) S(Q0) - const'(0)`` in closed form (``(2/pi) Q0^3/4`` and
    ``(2/pi) Q0^3/3`` unwindowed; the Lorch moments of ``sin(aQ)/(aQ)``
    otherwise).
    """
    q = _increasing_grid(q, "Q")
    fq = np.asarray(fq, dtype=float)
    weighted = fq * lorch_window(q, q[-1]) if lorch else fq
    slope = (2.0 / np.pi) * float(_trapezoid(weighted * q, x=q))
    if low_q_correction and q[0] != 0:
        q0 = float(q[0])
        if lorch:
            a = np.pi / float(q[-1])
            x = a * q0
            # (1/a) int_0^Q0 Q^k sin(aQ) dQ for k = 2 (coef) and k = 1 (const).
            f1 = (2.0 * x * np.sin(x) - (x * x - 2.0) * np.cos(x) - 2.0) / a**4
            f2 = (np.sin(x) - x * np.cos(x)) / a**3
        else:
            f1 = q0**4 / 4.0
            f2 = q0**3 / 3.0
        coef = (2.0 / np.pi) * f1 / q0
        const = (2.0 / np.pi) * f2
        if s0_target != 0.0:
            const = (1.0 - s0_target) * const + s0_target * coef
        slope += coef * (fq[0] / q0 + 1.0) - const
    return slope


def enforce_low_r(
    r: np.ndarray,
    gk: np.ndarray,
    *,
    cutoff: float,
    b_avg_sq: float,
) -> np.ndarray:
    """Classic stog parity: hard-replace ``GK(r)`` below ``cutoff`` with ``-<b>^2``.

    The Fortran stog applies this to its final RMC outputs (verified: the
    example run's ``scale_ft_rmc.gr`` equals ``-<b>^2`` to 16 digits for every
    ``r <= cutoff``). The published files therefore satisfy Keen Eq. 15 by
    construction — judge fit quality on the *pre*-enforcement residual.

    This is the special case of :func:`first_peak_zero` with the peak window
    outside ``[0, cutoff]`` (as in the validation example, where
    ``rmcmin = 2.65 > rmccut = 2.48``); use :func:`first_peak_zero` for the
    general Fortran semantics.
    """
    out = np.asarray(gk, dtype=float).copy()
    out[np.asarray(r, dtype=float) <= float(cutoff)] = -float(b_avg_sq)
    return out


def first_peak_zero(
    r: np.ndarray,
    g: np.ndarray,
    *,
    cutoff: float,
    peak_rmin: float,
    peak_rmax: float,
) -> np.ndarray:
    """Exact Fortran stog final ripple removal, in g(r) space.

    ``stog_new3.f90``: for the RMC outputs, ``g`` is zeroed wherever
    ``r <= cutoff`` **and** ``r`` lies outside the first-peak window
    ``[peak_rmin, peak_rmax]`` — removing sub-peak and inter-peak ripples
    while keeping a first peak that starts below ``cutoff``. Zeroed ``g``
    becomes ``-<b>^2`` in the written Keen G(r).
    """
    r = np.asarray(r, dtype=float)
    out = np.asarray(g, dtype=float).copy()
    zero = (r <= float(cutoff)) & (
        (r >= float(peak_rmax)) | (r <= float(peak_rmin))
    )
    out[zero] = 0.0
    return out
