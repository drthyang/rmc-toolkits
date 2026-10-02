// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/__tests__/fixtures/spaceGroups.js
//
// Ground truth for the symmetry tests: 230 space groups given by their ITA
// generators, plus the helpers that turn them into a crystal structure.
//
// Each fixture lists a few coordinate triplets which, closed under composition and
// combined with the centering translations, reproduce the group's full operation set.
// `multiplicity` is the general-position multiplicity in the conventional cell and acts
// as a checksum on the closure — a test asserts it, so a wrong generator set fails loudly
// rather than silently weakening the space-group assertions that depend on it.
//
// A structure is then built by expanding ONE generic point through the group. The orbit
// of a general position has exactly the symmetry of the group it came from — no more, no
// less — which is what makes these usable as space-group assertions.

import { itaGenerators } from '../../itaOperations.js';

/** Parse an ITA coordinate triplet ("-x+1/2,-y,z+1/2") into { R, t }. */
export function parseTriplet(text) {
    const R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    const t = [0, 0, 0];
    text.replace(/\s+/g, '').split(',').forEach((part, i) => {
        for (const token of part.match(/[+-]?[^+-]+/g) || []) {
            const sign = token.startsWith('-') ? -1 : 1;
            const body = token.replace(/^[+-]/, '');
            const axis = 'xyz'.indexOf(body[body.length - 1]);
            if (axis >= 0) {
                const coefficient = body.slice(0, -1);
                R[i][axis] = sign * (coefficient === '' ? 1 : fraction(coefficient));
            } else {
                t[i] += sign * fraction(body);
            }
        }
    });
    return { R, t: t.map(wrap) };
}

const fraction = (s) => (s.includes('/') ? Number(s.split('/')[0]) / Number(s.split('/')[1]) : Number(s));
const wrap = (v) => { const y = ((v % 1) + 1) % 1; return y > 1 - 1e-9 ? 0 : y; };

const opKey = (o) => o.R.flat().join(',') + '|' + o.t.map((v) => Math.round(wrap(v) * 12) % 12).join(',');

const compose = (a, b) => ({
    R: a.R.map((row) => [0, 1, 2].map((j) => row[0] * b.R[0][j] + row[1] * b.R[1][j] + row[2] * b.R[2][j])),
    t: [0, 1, 2].map((i) => wrap(a.R[i][0] * b.t[0] + a.R[i][1] * b.t[1] + a.R[i][2] * b.t[2] + a.t[i])),
});

/** Full operation set of a fixture: generators + centering, closed under composition. */
export function closeGroup({ generators, centering }) {
    const ops = new Map();
    const add = (o) => { const k = opKey(o); if (ops.has(k)) return false; ops.set(k, o); return true; };
    add({ R: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0, 0, 0] });
    for (const g of generators) add(parseTriplet(g));
    for (const c of centering) add({ R: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: c.split(',').map((v) => wrap(fraction(v))) });
    for (let pass = 0; pass < 12; pass += 1) {
        const snapshot = [...ops.values()];
        let grew = false;
        for (const a of snapshot) for (const b of snapshot) if (add(compose(a, b))) grew = true;
        if (!grew) break;
    }
    return [...ops.values()];
}

/** Lattice vectors (rows, Å) from the six cell parameters. */
export function cellVectors(a, b, c, alpha, beta, gamma) {
    const d = Math.PI / 180;
    const ca = Math.cos(alpha * d);
    const cb = Math.cos(beta * d);
    const cg = Math.cos(gamma * d);
    const sg = Math.sin(gamma * d);
    const cx = c * cb;
    const cy = c * (ca - cb * cg) / sg;
    return [[a, 0, 0], [b * cg, b * sg, 0], [cx, cy, Math.sqrt(Math.max(c * c - cx * cx - cy * cy, 1e-12))]];
}

// Deliberately unequal edges within each system, so the LATTICE has no accidental
// symmetry beyond the crystal system and the detected group comes from the atoms.
const CELL_FOR_SYSTEM = {
    triclinic: () => cellVectors(6.1, 7.3, 8.7, 81, 86, 94),
    monoclinic: () => cellVectors(7.1, 5.3, 9.7, 90, 104, 90),
    orthorhombic: () => cellVectors(8.1, 5.3, 6.7, 90, 90, 90),
    tetragonal: () => cellVectors(5.1, 5.1, 8.3, 90, 90, 90),
    trigonal: () => cellVectors(5.1, 5.1, 13.3, 90, 90, 120),
    hexagonal: () => cellVectors(5.1, 5.1, 13.3, 90, 90, 120),
    cubic: () => cellVectors(6.1, 6.1, 6.1, 90, 90, 90),
};

export const cellForSystem = (system) => CELL_FOR_SYSTEM[system]();

// A point on no symmetry element of any of the fixture groups.
const GENERIC_POINT = [0.137, 0.213, 0.061];

/** Expand a generic point through `ops` into a one-orbit basis for the finder. */
export function orbitBasis(ops, element = 'A', point = GENERIC_POINT) {
    const basis = [];
    for (const { R, t } of ops) {
        const q = [0, 1, 2].map((i) => wrap(R[i][0] * point[0] + R[i][1] * point[1] + R[i][2] * point[2] + t[i]));
        const dup = basis.some((s) => s.frac.every((v, i) => Math.abs((((v - q[i]) % 1) + 1.5) % 1 - 0.5) < 1e-4));
        if (!dup) basis.push({ el: element, frac: q });
    }
    return basis;
}

/**
 * Cell + basis for a fixture, ready to hand to the space-group finder.
 *
 * Two orbits: the generic one, which fixes the symmetry, plus a second element at the
 * origin. The anchor is there for speed, not for symmetry — the finder seeds candidate
 * translations from the RAREST element, and a cubic general orbit on its own is 192
 * atoms, which makes that search quadratically expensive. The origin lies on a special
 * position in every one of these groups, so the anchor orbit is small, and being a full
 * orbit it changes the structure's symmetry in neither direction.
 */
export function structureFor(fixture) {
    const ops = closeGroup(fixture);
    const basis = [...orbitBasis(ops, 'A'), ...orbitBasis(ops, 'B', [0, 0, 0])];
    return { ops, A: cellForSystem(fixture.system), basis };
}

// The 230 groups: number, symbol, crystal system, general multiplicity. The generators and
// centring are the app's own ITA table (itaOperations.js, pinned to spglib), so the tests
// and the CIF export's origin move use one copy of the same operations.
const CENTRING_TEXT = {
    P: ['0,0,0'], A: ['0,0,0', '0,1/2,1/2'], B: ['0,0,0', '1/2,0,1/2'], C: ['0,0,0', '1/2,1/2,0'],
    I: ['0,0,0', '1/2,1/2,1/2'], F: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'], R: ['0,0,0', '2/3,1/3,1/3', '1/3,2/3,2/3'],
};

const GROUPS = [
    [1, 'P1', 'triclinic', 1],
    [2, 'P-1', 'triclinic', 2],
    [3, 'P2', 'monoclinic', 2],
    [4, 'P2_1', 'monoclinic', 2],
    [5, 'C2', 'monoclinic', 4],
    [6, 'Pm', 'monoclinic', 2],
    [7, 'Pc', 'monoclinic', 2],
    [8, 'Cm', 'monoclinic', 4],
    [9, 'Cc', 'monoclinic', 4],
    [10, 'P2/m', 'monoclinic', 4],
    [11, 'P2_1/m', 'monoclinic', 4],
    [12, 'C2/m', 'monoclinic', 8],
    [13, 'P2/c', 'monoclinic', 4],
    [14, 'P2_1/c', 'monoclinic', 4],
    [15, 'C2/c', 'monoclinic', 8],
    [16, 'P222', 'orthorhombic', 4],
    [17, 'P222_1', 'orthorhombic', 4],
    [18, 'P2_12_12', 'orthorhombic', 4],
    [19, 'P2_12_12_1', 'orthorhombic', 4],
    [20, 'C222_1', 'orthorhombic', 8],
    [21, 'C222', 'orthorhombic', 8],
    [22, 'F222', 'orthorhombic', 16],
    [23, 'I222', 'orthorhombic', 8],
    [24, 'I2_12_12_1', 'orthorhombic', 8],
    [25, 'Pmm2', 'orthorhombic', 4],
    [26, 'Pmc2_1', 'orthorhombic', 4],
    [27, 'Pcc2', 'orthorhombic', 4],
    [28, 'Pma2', 'orthorhombic', 4],
    [29, 'Pca2_1', 'orthorhombic', 4],
    [30, 'Pnc2', 'orthorhombic', 4],
    [31, 'Pmn2_1', 'orthorhombic', 4],
    [32, 'Pba2', 'orthorhombic', 4],
    [33, 'Pna2_1', 'orthorhombic', 4],
    [34, 'Pnn2', 'orthorhombic', 4],
    [35, 'Cmm2', 'orthorhombic', 8],
    [36, 'Cmc2_1', 'orthorhombic', 8],
    [37, 'Ccc2', 'orthorhombic', 8],
    [38, 'Amm2', 'orthorhombic', 8],
    [39, 'Aem2', 'orthorhombic', 8],
    [40, 'Ama2', 'orthorhombic', 8],
    [41, 'Aea2', 'orthorhombic', 8],
    [42, 'Fmm2', 'orthorhombic', 16],
    [43, 'Fdd2', 'orthorhombic', 16],
    [44, 'Imm2', 'orthorhombic', 8],
    [45, 'Iba2', 'orthorhombic', 8],
    [46, 'Ima2', 'orthorhombic', 8],
    [47, 'Pmmm', 'orthorhombic', 8],
    [48, 'Pnnn', 'orthorhombic', 8],
    [49, 'Pccm', 'orthorhombic', 8],
    [50, 'Pban', 'orthorhombic', 8],
    [51, 'Pmma', 'orthorhombic', 8],
    [52, 'Pnna', 'orthorhombic', 8],
    [53, 'Pmna', 'orthorhombic', 8],
    [54, 'Pcca', 'orthorhombic', 8],
    [55, 'Pbam', 'orthorhombic', 8],
    [56, 'Pccn', 'orthorhombic', 8],
    [57, 'Pbcm', 'orthorhombic', 8],
    [58, 'Pnnm', 'orthorhombic', 8],
    [59, 'Pmmn', 'orthorhombic', 8],
    [60, 'Pbcn', 'orthorhombic', 8],
    [61, 'Pbca', 'orthorhombic', 8],
    [62, 'Pnma', 'orthorhombic', 8],
    [63, 'Cmcm', 'orthorhombic', 16],
    [64, 'Cmce', 'orthorhombic', 16],
    [65, 'Cmmm', 'orthorhombic', 16],
    [66, 'Cccm', 'orthorhombic', 16],
    [67, 'Cmme', 'orthorhombic', 16],
    [68, 'Ccce', 'orthorhombic', 16],
    [69, 'Fmmm', 'orthorhombic', 32],
    [70, 'Fddd', 'orthorhombic', 32],
    [71, 'Immm', 'orthorhombic', 16],
    [72, 'Ibam', 'orthorhombic', 16],
    [73, 'Ibca', 'orthorhombic', 16],
    [74, 'Imma', 'orthorhombic', 16],
    [75, 'P4', 'tetragonal', 4],
    [76, 'P4_1', 'tetragonal', 4],
    [77, 'P4_2', 'tetragonal', 4],
    [78, 'P4_3', 'tetragonal', 4],
    [79, 'I4', 'tetragonal', 8],
    [80, 'I4_1', 'tetragonal', 8],
    [81, 'P-4', 'tetragonal', 4],
    [82, 'I-4', 'tetragonal', 8],
    [83, 'P4/m', 'tetragonal', 8],
    [84, 'P4_2/m', 'tetragonal', 8],
    [85, 'P4/n', 'tetragonal', 8],
    [86, 'P4_2/n', 'tetragonal', 8],
    [87, 'I4/m', 'tetragonal', 16],
    [88, 'I4_1/a', 'tetragonal', 16],
    [89, 'P422', 'tetragonal', 8],
    [90, 'P42_12', 'tetragonal', 8],
    [91, 'P4_122', 'tetragonal', 8],
    [92, 'P4_12_12', 'tetragonal', 8],
    [93, 'P4_222', 'tetragonal', 8],
    [94, 'P4_22_12', 'tetragonal', 8],
    [95, 'P4_322', 'tetragonal', 8],
    [96, 'P4_32_12', 'tetragonal', 8],
    [97, 'I422', 'tetragonal', 16],
    [98, 'I4_122', 'tetragonal', 16],
    [99, 'P4mm', 'tetragonal', 8],
    [100, 'P4bm', 'tetragonal', 8],
    [101, 'P4_2cm', 'tetragonal', 8],
    [102, 'P4_2nm', 'tetragonal', 8],
    [103, 'P4cc', 'tetragonal', 8],
    [104, 'P4nc', 'tetragonal', 8],
    [105, 'P4_2mc', 'tetragonal', 8],
    [106, 'P4_2bc', 'tetragonal', 8],
    [107, 'I4mm', 'tetragonal', 16],
    [108, 'I4cm', 'tetragonal', 16],
    [109, 'I4_1md', 'tetragonal', 16],
    [110, 'I4_1cd', 'tetragonal', 16],
    [111, 'P-42m', 'tetragonal', 8],
    [112, 'P-42c', 'tetragonal', 8],
    [113, 'P-42_1m', 'tetragonal', 8],
    [114, 'P-42_1c', 'tetragonal', 8],
    [115, 'P-4m2', 'tetragonal', 8],
    [116, 'P-4c2', 'tetragonal', 8],
    [117, 'P-4b2', 'tetragonal', 8],
    [118, 'P-4n2', 'tetragonal', 8],
    [119, 'I-4m2', 'tetragonal', 16],
    [120, 'I-4c2', 'tetragonal', 16],
    [121, 'I-42m', 'tetragonal', 16],
    [122, 'I-42d', 'tetragonal', 16],
    [123, 'P4/mmm', 'tetragonal', 16],
    [124, 'P4/mcc', 'tetragonal', 16],
    [125, 'P4/nbm', 'tetragonal', 16],
    [126, 'P4/nnc', 'tetragonal', 16],
    [127, 'P4/mbm', 'tetragonal', 16],
    [128, 'P4/mnc', 'tetragonal', 16],
    [129, 'P4/nmm', 'tetragonal', 16],
    [130, 'P4/ncc', 'tetragonal', 16],
    [131, 'P4_2/mmc', 'tetragonal', 16],
    [132, 'P4_2/mcm', 'tetragonal', 16],
    [133, 'P4_2/nbc', 'tetragonal', 16],
    [134, 'P4_2/nnm', 'tetragonal', 16],
    [135, 'P4_2/mbc', 'tetragonal', 16],
    [136, 'P4_2/mnm', 'tetragonal', 16],
    [137, 'P4_2/nmc', 'tetragonal', 16],
    [138, 'P4_2/ncm', 'tetragonal', 16],
    [139, 'I4/mmm', 'tetragonal', 32],
    [140, 'I4/mcm', 'tetragonal', 32],
    [141, 'I4_1/amd', 'tetragonal', 32],
    [142, 'I4_1/acd', 'tetragonal', 32],
    [143, 'P3', 'trigonal', 3],
    [144, 'P3_1', 'trigonal', 3],
    [145, 'P3_2', 'trigonal', 3],
    [146, 'R3', 'trigonal', 9],
    [147, 'P-3', 'trigonal', 6],
    [148, 'R-3', 'trigonal', 18],
    [149, 'P312', 'trigonal', 6],
    [150, 'P321', 'trigonal', 6],
    [151, 'P3_112', 'trigonal', 6],
    [152, 'P3_121', 'trigonal', 6],
    [153, 'P3_212', 'trigonal', 6],
    [154, 'P3_221', 'trigonal', 6],
    [155, 'R32', 'trigonal', 18],
    [156, 'P3m1', 'trigonal', 6],
    [157, 'P31m', 'trigonal', 6],
    [158, 'P3c1', 'trigonal', 6],
    [159, 'P31c', 'trigonal', 6],
    [160, 'R3m', 'trigonal', 18],
    [161, 'R3c', 'trigonal', 18],
    [162, 'P-31m', 'trigonal', 12],
    [163, 'P-31c', 'trigonal', 12],
    [164, 'P-3m1', 'trigonal', 12],
    [165, 'P-3c1', 'trigonal', 12],
    [166, 'R-3m', 'trigonal', 36],
    [167, 'R-3c', 'trigonal', 36],
    [168, 'P6', 'hexagonal', 6],
    [169, 'P6_1', 'hexagonal', 6],
    [170, 'P6_5', 'hexagonal', 6],
    [171, 'P6_2', 'hexagonal', 6],
    [172, 'P6_4', 'hexagonal', 6],
    [173, 'P6_3', 'hexagonal', 6],
    [174, 'P-6', 'hexagonal', 6],
    [175, 'P6/m', 'hexagonal', 12],
    [176, 'P6_3/m', 'hexagonal', 12],
    [177, 'P622', 'hexagonal', 12],
    [178, 'P6_122', 'hexagonal', 12],
    [179, 'P6_522', 'hexagonal', 12],
    [180, 'P6_222', 'hexagonal', 12],
    [181, 'P6_422', 'hexagonal', 12],
    [182, 'P6_322', 'hexagonal', 12],
    [183, 'P6mm', 'hexagonal', 12],
    [184, 'P6cc', 'hexagonal', 12],
    [185, 'P6_3cm', 'hexagonal', 12],
    [186, 'P6_3mc', 'hexagonal', 12],
    [187, 'P-6m2', 'hexagonal', 12],
    [188, 'P-6c2', 'hexagonal', 12],
    [189, 'P-62m', 'hexagonal', 12],
    [190, 'P-62c', 'hexagonal', 12],
    [191, 'P6/mmm', 'hexagonal', 24],
    [192, 'P6/mcc', 'hexagonal', 24],
    [193, 'P6_3/mcm', 'hexagonal', 24],
    [194, 'P6_3/mmc', 'hexagonal', 24],
    [195, 'P23', 'cubic', 12],
    [196, 'F23', 'cubic', 48],
    [197, 'I23', 'cubic', 24],
    [198, 'P2_13', 'cubic', 12],
    [199, 'I2_13', 'cubic', 24],
    [200, 'Pm-3', 'cubic', 24],
    [201, 'Pn-3', 'cubic', 24],
    [202, 'Fm-3', 'cubic', 96],
    [203, 'Fd-3', 'cubic', 96],
    [204, 'Im-3', 'cubic', 48],
    [205, 'Pa-3', 'cubic', 24],
    [206, 'Ia-3', 'cubic', 48],
    [207, 'P432', 'cubic', 24],
    [208, 'P4_232', 'cubic', 24],
    [209, 'F432', 'cubic', 96],
    [210, 'F4_132', 'cubic', 96],
    [211, 'I432', 'cubic', 48],
    [212, 'P4_332', 'cubic', 24],
    [213, 'P4_132', 'cubic', 24],
    [214, 'I4_132', 'cubic', 48],
    [215, 'P-43m', 'cubic', 24],
    [216, 'F-43m', 'cubic', 96],
    [217, 'I-43m', 'cubic', 48],
    [218, 'P-43n', 'cubic', 24],
    [219, 'F-43c', 'cubic', 96],
    [220, 'I-43d', 'cubic', 48],
    [221, 'Pm-3m', 'cubic', 48],
    [222, 'Pn-3n', 'cubic', 48],
    [223, 'Pm-3n', 'cubic', 48],
    [224, 'Pn-3m', 'cubic', 48],
    [225, 'Fm-3m', 'cubic', 192],
    [226, 'Fm-3c', 'cubic', 192],
    [227, 'Fd-3m', 'cubic', 192],
    [228, 'Fd-3c', 'cubic', 192],
    [229, 'Im-3m', 'cubic', 96],
    [230, 'Ia-3d', 'cubic', 96],
];

export const SPACE_GROUP_FIXTURES = GROUPS.map(([number, symbol, system, multiplicity]) => {
    const { letter, generators } = itaGenerators(number);
    return { number, symbol, system, multiplicity, generators, centering: CENTRING_TEXT[letter] };
});
