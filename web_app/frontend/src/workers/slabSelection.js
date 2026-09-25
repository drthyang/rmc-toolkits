// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Slice-definition and kernel-readout helpers shared by the Structure page, the
// KDE worker and tests (pure functions, no side effects).
//
// Slab membership, shared by the KDE worker (localKdeWorker.js makeSlab) and the
// Structure page's Slab-In-Cell highlight (StructurePage.jsx inActiveSlab), and
// identical to rmc_toolkits/kde.py (SLAB_FACE_TOLERANCE): an atom is in the slab
// when its depth, normalised to [0, 1] across the unit cube's projection range,
// is within thickness / 2 + SLAB_FACE_TOLERANCE of the slice centre. Atoms of an
// ideal or unrelaxed configuration sit exactly on slider-reachable faces
// (z = 0.125 against zCenter = 0.165, thickness = 0.08), where two roundings of
// the same inequality disagree and a whole site flips in or out; the tolerance
// includes them in every runtime. This module has no side effects, so the page
// can import it without pulling in the worker.
export const SLAB_FACE_TOLERANCE = 1e-9;

export const isInSlab = (normalizedDepth, zCenter, thickness) => (
    Math.abs(normalizedDepth - zCenter) <= thickness / 2 + SLAB_FACE_TOLERANCE
);

// The custom slice is the lattice-plane family x . h = const with Miller indices
// (h k l): an atom's depth is its fractional coordinate dotted with h, so the
// slab normal is h a* + k b* + l c* -- not the real-space direction [h k l],
// which differs from it in any non-orthogonal cell (30 deg for [1 0 0] in a
// hexagonal cell). Labels therefore use parentheses, never square brackets.
export const millerPlaneLabel = (indices) => `(${indices
    .map((value) => Number(value).toLocaleString(undefined, { maximumFractionDigits: 2 }))
    .join(' ')})`;

// The same plane for a file name: locale-independent, only digits, '-', '.',
// '_' and the parentheses, e.g. "(1_1_0)" or "(-1_0.5_2)".
export const millerPlaneFileLabel = (indices) => `(${indices
    .map((value) => {
        const number = Number(value);
        return Number.isFinite(number) ? String(Math.round(number * 100) / 100) : '0';
    })
    .join('_')})`;

// The KDE kernel in Angstrom. `covariance` is H in the slice's (u, v)
// fractional coordinates (the payload's kernel.covariance); uCartesian and
// vCartesian are the Cartesian (Angstrom) images of the in-plane axes u and v.
// In real space the kernel is M H M^T with M = [uCartesian vCartesian], whose two
// nonzero eigenvalues are those of H G, G = M^T M the in-plane metric. Returns
// the principal sigmas { minor, major } in Angstrom.
export const kernelSigmaAngstrom = (covariance, uCartesian, vCartesian) => {
    const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    const g00 = dot3(uCartesian, uCartesian);
    const g01 = dot3(uCartesian, vCartesian);
    const g11 = dot3(vCartesian, vCartesian);
    const [[h00, h01], [, h11]] = covariance;
    const trace = h00 * g00 + 2 * h01 * g01 + h11 * g11;
    const determinant = (h00 * h11 - h01 * h01) * (g00 * g11 - g01 * g01);
    const major = 0.5 * trace + Math.sqrt(Math.max(0, 0.25 * trace * trace - determinant));
    const minor = major > 0 ? determinant / major : 0;
    return { minor: Math.sqrt(Math.max(minor, 0)), major: Math.sqrt(Math.max(major, 0)) };
};

// Above this aspect ratio (in Angstrom) the page notes that the kernel's shape
// is an artefact of the slab's site layout (docs/algorithms/structure.md,
// "The kernel's shape follows the slab's site layout").
export const KERNEL_ANISOTROPY_NOTE = 3;

// The slab's real thickness in Angstrom. The sliders' zCenter/thickness are
// fractions of the unit cube's projection range along the (unit) normal n
// (range = [d_min, d_max]); the depth d = n . x of a fractional position x has
// the Cartesian gradient A^-1 n = n1 a* + n2 b* + n3 c* (rows of A = unitVectors,
// in Angstrom; reciprocal vectors without the 2 pi), so a depth interval
// thickness * (d_max - d_min) is thickness * (d_max - d_min) / |A^-1 n| Angstrom
// = thickness * (|h| + |k| + |l|) * d_hkl.
export const slabThicknessAngstrom = (thickness, normal, range, unitVectors) => {
    const [a, b, c] = unitVectors;
    const cross = (u, v) => [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    const bc = cross(b, c);
    const ca = cross(c, a);
    const ab = cross(a, b);
    const volume = a[0] * bc[0] + a[1] * bc[1] + a[2] * bc[2];
    if (!(Math.abs(volume) > 0)) return Number.NaN;
    const gradient = [0, 1, 2].map((i) => (normal[0] * bc[i] + normal[1] * ca[i] + normal[2] * ab[i]) / volume);
    const length = Math.hypot(...gradient);
    const span = range[1] - range[0] || 1;
    return length > 0 ? (thickness * span) / length : Number.NaN;
};
