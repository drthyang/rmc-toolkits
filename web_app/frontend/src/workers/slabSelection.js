// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Slice-definition helpers shared by the Structure page, the KDE worker and tests.
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
