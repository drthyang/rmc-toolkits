// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Relating a site's PCA (principal-axis) frame to the crystallographic (a, b, c)
// frame. Everything in the PCA-KDE scene lives in one orthonormal Cartesian frame
// (Ångström): the displacement clouds are mapped to Cartesian through the supercell
// lattice, so the principal axes (eigenvectors of the cloud covariance) and the
// unit-cell vectors a, b, c are expressed in the SAME Cartesian basis. That shared
// basis is what makes the transformation below exact — no re-referencing is needed.
//
// Conventions (all matrices are arrays of rows, m[row][col]):
//   • Unit-cell matrix  M — rows are the unit-cell vectors a, b, c in Cartesian Å.
//     A crystallographic (fractional) vector f maps to Cartesian as  r = Mᵀ·f
//     (r_k = Σ_i f_i · M[i][k]).  Its inverse is  f = M⁻ᵀ·r.
//   • Principal-axis matrix  P — rows are the PCA axes (orthonormal, Cartesian).
//     A PCA-frame vector q maps to Cartesian as  r = Pᵀ·q, and  q = P·r  (P⁻¹ = Pᵀ).
//
// Composing the two gives the transformation matrix the UI reports:
//   • fractional → PCA :  q = (P·Mᵀ)·f          →  fracToPca = P·Mᵀ
//   • PCA → fractional :  f = (M⁻ᵀ·Pᵀ)·q = (P·Mᵀ)⁻¹·q  →  pcaToFrac = fracToPca⁻¹
//
// This is a pure module (no three.js / DOM) so the geometry can be unit-tested
// against hand-derived cases.

// --- 3×3 linear algebra (rows-of-arrays) -------------------------------------

const dot3 = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
const norm3 = (u) => Math.hypot(u[0], u[1], u[2]);

export const transpose3 = (m) => [
    [m[0][0], m[1][0], m[2][0]],
    [m[0][1], m[1][1], m[2][1]],
    [m[0][2], m[1][2], m[2][2]]
];

export const multiply3 = (a, b) => {
    const out = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < 3; i += 1) {
        for (let j = 0; j < 3; j += 1) {
            out[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
    return out;
};

export const matVec3 = (m, v) => [
    m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
    m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
    m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]
];

export const determinant3 = (m) =>
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

// Inverse by the adjugate / determinant. Returns null when the matrix is singular
// (a collapsed cell), so callers can degrade gracefully rather than emit NaNs.
export const invert3 = (m) => {
    const det = determinant3(m);
    if (!Number.isFinite(det) || Math.abs(det) < 1e-12) return null;
    const inv = 1 / det;
    return [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv
        ]
    ];
};

// --- Crystallographic frame ---------------------------------------------------

/**
 * Unit-cell edge vectors (rows a, b, c) in Cartesian Å, from the supercell lattice
 * vectors and the supercell repeat counts. This is exactly the frame the rest of
 * the app already uses (the unit-cell structure view divides the lattice by the
 * supercell), so the vectors returned here sit in the same Cartesian basis as the
 * PCA axes of every site's displacement cloud.
 */
export const unitCellVectors = (latticeVectors, supercell) => {
    if (!Array.isArray(latticeVectors) || latticeVectors.length !== 3) return null;
    if (!Array.isArray(supercell) || supercell.length !== 3) return null;
    return latticeVectors.map((row, i) => {
        const n = Math.abs(supercell[i]) > 0 ? supercell[i] : 1;
        return [row[0] / n, row[1] / n, row[2] / n];
    });
};

const RAD_TO_DEG = 180 / Math.PI;
const clamp = (v, lo, hi) => (v < lo ? lo : v > hi ? hi : v);

/**
 * Transformation matrices between the crystallographic (fractional a/b/c) frame
 * and a site's PCA frame, both derived in the shared Cartesian basis.
 *
 * @param {number[][]} pcaAxes  - principal axes as rows (unit Cartesian vectors).
 * @param {number[][]} unitCell - unit-cell vectors as rows a, b, c (Cartesian Å).
 * @returns {{fracToPca:number[][], pcaToFrac:number[][], unitCell:number[][]}|null}
 *   fracToPca maps a fractional displacement to PCA coordinates (Å); pcaToFrac is
 *   its inverse. null when the cell is singular.
 */
export const crystalPcaTransforms = (pcaAxes, unitCell) => {
    if (!Array.isArray(pcaAxes) || pcaAxes.length !== 3 || !unitCell) return null;
    const fracToPca = multiply3(pcaAxes, transpose3(unitCell)); // P·Mᵀ
    const pcaToFrac = invert3(fracToPca);                       // (P·Mᵀ)⁻¹ = M⁻ᵀ·Pᵀ
    if (!pcaToFrac) return null;
    return { fracToPca, pcaToFrac, unitCell };
};

/**
 * Orientation of each principal axis relative to the crystallographic axes.
 * For every PCA axis this returns, in the shared Cartesian frame:
 *   • cosines[i] / anglesDeg[i] — direction cosine and angle (deg) to unit-cell
 *     axis i (a=0, b=1, c=2). Angles are always well defined, even for an oblique
 *     cell, because they measure the angle between two Cartesian directions.
 *   • crystalDirection [u v w] — the crystallographic direction the PCA axis points
 *     along (the fractional components M⁻ᵀ·axis, i.e. a column of pcaToFrac),
 *     normalised so the largest component magnitude is 1.
 *   • dominant — the crystal axis it lies closest to, with that angle.
 *
 * @param {number[][]} pcaAxes  - principal axes as rows (unit Cartesian vectors).
 * @param {number[][]} unitCell - unit-cell vectors as rows a, b, c (Cartesian Å).
 */
export const principalAxisOrientation = (pcaAxes, unitCell) => {
    const transforms = crystalPcaTransforms(pcaAxes, unitCell);
    if (!transforms) return null;

    const cellLengths = unitCell.map(norm3);
    const cellUnit = unitCell.map((row, i) => {
        const len = cellLengths[i] || 1;
        return [row[0] / len, row[1] / len, row[2] / len];
    });
    const invT = transpose3(invert3(unitCell)); // M⁻ᵀ : Cartesian → fractional
    const labels = ['a', 'b', 'c'];

    const perAxis = pcaAxes.map((axis) => {
        const cosines = cellUnit.map((u) => clamp(dot3(axis, u), -1, 1));
        const anglesDeg = cosines.map((c) => Math.acos(c) * RAD_TO_DEG);
        const frac = invT ? matVec3(invT, axis) : [0, 0, 0];
        const maxComp = Math.max(Math.abs(frac[0]), Math.abs(frac[1]), Math.abs(frac[2])) || 1;
        const crystalDirection = frac.map((v) => v / maxComp);
        // Closest crystal axis = largest |cos| (smallest angle).
        let dom = 0;
        for (let i = 1; i < 3; i += 1) {
            if (Math.abs(cosines[i]) > Math.abs(cosines[dom])) dom = i;
        }
        return {
            cosines,
            anglesDeg,
            crystalDirection,
            dominant: { index: dom, label: labels[dom], angleDeg: anglesDeg[dom] }
        };
    });

    return {
        perAxis,
        cellLengths,
        fracToPca: transforms.fracToPca,
        pcaToFrac: transforms.pcaToFrac
    };
};
