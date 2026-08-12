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

export const SPACE_GROUP_FIXTURES = [
    {
        number: 1, symbol: 'P1', system: 'triclinic', multiplicity: 1,
        generators: [],
        centering: ['0,0,0'],
    },
    {
        number: 2, symbol: 'P-1', system: 'triclinic', multiplicity: 2,
        generators: ['-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 3, symbol: 'P2', system: 'monoclinic', multiplicity: 2,
        generators: ['-x,y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 4, symbol: 'P2_1', system: 'monoclinic', multiplicity: 2,
        generators: ['-x,y+1/2,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 5, symbol: 'C2', system: 'monoclinic', multiplicity: 4,
        generators: ['-x,y,-z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 6, symbol: 'Pm', system: 'monoclinic', multiplicity: 2,
        generators: ['x,-y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 7, symbol: 'Pc', system: 'monoclinic', multiplicity: 2,
        generators: ['x,-y,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 8, symbol: 'Cm', system: 'monoclinic', multiplicity: 4,
        generators: ['x,-y,z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 9, symbol: 'Cc', system: 'monoclinic', multiplicity: 4,
        generators: ['x,-y,z+1/2'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 10, symbol: 'P2/m', system: 'monoclinic', multiplicity: 4,
        generators: ['-x,y,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 11, symbol: 'P2_1/m', system: 'monoclinic', multiplicity: 4,
        generators: ['-x,y+1/2,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 12, symbol: 'C2/m', system: 'monoclinic', multiplicity: 8,
        generators: ['-x,y,-z', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 13, symbol: 'P2/c', system: 'monoclinic', multiplicity: 4,
        generators: ['-x,y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 14, symbol: 'P2_1/c', system: 'monoclinic', multiplicity: 4,
        generators: ['-x,y+1/2,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 15, symbol: 'C2/c', system: 'monoclinic', multiplicity: 8,
        generators: ['-x,y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 16, symbol: 'P222', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z', 'x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 17, symbol: 'P222_1', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z+1/2', 'x,-y,-z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 18, symbol: 'P2_12_12', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z', '-x+1/2,y+1/2,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 19, symbol: 'P2_12_12_1', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x+1/2,-y,z+1/2', '-x,y+1/2,-z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 20, symbol: 'C222_1', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z+1/2', 'x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 21, symbol: 'C222', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', 'x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 22, symbol: 'F222', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x,-y,z', 'x,-y,-z'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 23, symbol: 'I222', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', 'x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 24, symbol: 'I2_12_12_1', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y,z', '-x,y,-z+1/2'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 25, symbol: 'Pmm2', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z', 'x,-y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 26, symbol: 'Pmc2_1', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z+1/2', 'x,-y,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 27, symbol: 'Pcc2', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z', 'x,-y,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 28, symbol: 'Pma2', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z', 'x+1/2,-y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 29, symbol: 'Pca2_1', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z+1/2', 'x+1/2,-y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 30, symbol: 'Pnc2', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z', 'x,-y+1/2,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 31, symbol: 'Pmn2_1', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x+1/2,-y,z+1/2', '-x,y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 32, symbol: 'Pba2', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z', 'x+1/2,-y+1/2,z'],
        centering: ['0,0,0'],
    },
    {
        number: 33, symbol: 'Pna2_1', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z+1/2', 'x+1/2,-y+1/2,z'],
        centering: ['0,0,0'],
    },
    {
        number: 34, symbol: 'Pnn2', system: 'orthorhombic', multiplicity: 4,
        generators: ['-x,-y,z', 'x+1/2,-y+1/2,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 35, symbol: 'Cmm2', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', 'x,-y,z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 36, symbol: 'Cmc2_1', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z+1/2', 'x,-y,z+1/2'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 37, symbol: 'Ccc2', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', 'x,-y,z+1/2'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 38, symbol: 'Amm2', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', 'x,-y,z'],
        centering: ['0,0,0', '0,1/2,1/2'],
    },
    {
        number: 39, symbol: 'Aem2', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', 'x,-y+1/2,z'],
        centering: ['0,0,0', '0,1/2,1/2'],
    },
    {
        number: 40, symbol: 'Ama2', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', 'x+1/2,-y,z'],
        centering: ['0,0,0', '0,1/2,1/2'],
    },
    {
        number: 41, symbol: 'Aea2', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', '-x+1/2,y+1/2,z'],
        centering: ['0,0,0', '0,1/2,1/2'],
    },
    {
        number: 42, symbol: 'Fmm2', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x,-y,z', '-x,y,z'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 43, symbol: 'Fdd2', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x,-y,z', '-x+1/4,y+1/4,z+1/4'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 44, symbol: 'Imm2', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', '-x,y,z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 45, symbol: 'Iba2', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', '-x,y,z+1/2'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 46, symbol: 'Ima2', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', '-x+1/2,y,z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 47, symbol: 'Pmmm', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 48, symbol: 'Pnnn', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y+1/2,z', 'x,-y+1/2,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 49, symbol: 'Pccm', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', 'x,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 50, symbol: 'Pban', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y+1/2,z', 'x,-y+1/2,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 51, symbol: 'Pmma', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y,z', 'x+1/2,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 52, symbol: 'Pnna', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y,z', 'x,-y+1/2,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 53, symbol: 'Pmna', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y,z+1/2', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 54, symbol: 'Pcca', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y,z', 'x+1/2,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 55, symbol: 'Pbam', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', 'x+1/2,-y+1/2,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 56, symbol: 'Pccn', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y+1/2,z', 'x+1/2,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 57, symbol: 'Pbcm', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z+1/2', 'x,-y+1/2,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 58, symbol: 'Pnnm', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x,-y,z', 'x+1/2,-y+1/2,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 59, symbol: 'Pmmn', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y+1/2,z', 'x+1/2,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 60, symbol: 'Pbcn', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y+1/2,z+1/2', 'x+1/2,-y+1/2,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 61, symbol: 'Pbca', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y,z+1/2', 'x+1/2,-y+1/2,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 62, symbol: 'Pnma', system: 'orthorhombic', multiplicity: 8,
        generators: ['-x+1/2,-y,z+1/2', 'x+1/2,-y+1/2,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 63, symbol: 'Cmcm', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x,-y,z+1/2', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 64, symbol: 'Cmce', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x+1/2,-y,z+1/2', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 65, symbol: 'Cmmm', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x,-y,z', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 66, symbol: 'Cccm', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x,-y,z', 'x,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 67, symbol: 'Cmme', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x+1/2,-y,z', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 68, symbol: 'Ccce', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x+1/2,-y,z', 'x+1/2,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,0'],
    },
    {
        number: 69, symbol: 'Fmmm', system: 'orthorhombic', multiplicity: 32,
        generators: ['-x,-y,z', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 70, symbol: 'Fddd', system: 'orthorhombic', multiplicity: 32,
        generators: ['-x+1/4,-y+1/4,z', 'x,-y+1/4,-z+1/4', '-x,-y,-z'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 71, symbol: 'Immm', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x,-y,z', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 72, symbol: 'Ibam', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x,-y,z', 'x,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 73, symbol: 'Ibca', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x,-y+1/2,z', 'x,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 74, symbol: 'Imma', system: 'orthorhombic', multiplicity: 16,
        generators: ['-x,-y+1/2,z', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 75, symbol: 'P4', system: 'tetragonal', multiplicity: 4,
        generators: ['-y,x,z'],
        centering: ['0,0,0'],
    },
    {
        number: 76, symbol: 'P4_1', system: 'tetragonal', multiplicity: 4,
        generators: ['-y,x,z+1/4'],
        centering: ['0,0,0'],
    },
    {
        number: 77, symbol: 'P4_2', system: 'tetragonal', multiplicity: 4,
        generators: ['-y,x,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 78, symbol: 'P4_3', system: 'tetragonal', multiplicity: 4,
        generators: ['-y,x,z+3/4'],
        centering: ['0,0,0'],
    },
    {
        number: 79, symbol: 'I4', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 80, symbol: 'I4_1', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x+1/2,z+1/4'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 81, symbol: 'P-4', system: 'tetragonal', multiplicity: 4,
        generators: ['y,-x,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 82, symbol: 'I-4', system: 'tetragonal', multiplicity: 8,
        generators: ['y,-x,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 83, symbol: 'P4/m', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 84, symbol: 'P4_2/m', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 85, symbol: 'P4/n', system: 'tetragonal', multiplicity: 8,
        generators: ['-y+1/2,x,z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 86, symbol: 'P4_2/n', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x+1/2,z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 87, symbol: 'I4/m', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x,z', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 88, symbol: 'I4_1/a', system: 'tetragonal', multiplicity: 16,
        generators: ['-y+3/4,x+1/4,z+1/4', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 89, symbol: 'P422', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z', 'x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 90, symbol: 'P42_12', system: 'tetragonal', multiplicity: 8,
        generators: ['-y+1/2,x+1/2,z', 'y,x,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 91, symbol: 'P4_122', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z+1/4', '-x,y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 92, symbol: 'P4_12_12', system: 'tetragonal', multiplicity: 8,
        generators: ['-y+1/2,x+1/2,z+1/4', 'y,x,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 93, symbol: 'P4_222', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z+1/2', 'x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 94, symbol: 'P4_22_12', system: 'tetragonal', multiplicity: 8,
        generators: ['-y+1/2,x+1/2,z+1/2', 'y,x,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 95, symbol: 'P4_322', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z+3/4', '-x,y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 96, symbol: 'P4_32_12', system: 'tetragonal', multiplicity: 8,
        generators: ['-y+1/2,x+1/2,z+3/4', 'y,x,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 97, symbol: 'I422', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x,z', 'x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 98, symbol: 'I4_122', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x+1/2,z+1/4', '-y,-x,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 99, symbol: 'P4mm', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z', 'x,-y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 100, symbol: 'P4bm', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z', 'x+1/2,-y+1/2,z'],
        centering: ['0,0,0'],
    },
    {
        number: 101, symbol: 'P4_2cm', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z+1/2', 'x,-y,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 102, symbol: 'P4_2nm', system: 'tetragonal', multiplicity: 8,
        generators: ['-y+1/2,x+1/2,z+1/2', 'y,x,z'],
        centering: ['0,0,0'],
    },
    {
        number: 103, symbol: 'P4cc', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z', 'x,-y,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 104, symbol: 'P4nc', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z', 'x+1/2,-y+1/2,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 105, symbol: 'P4_2mc', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z+1/2', 'x,-y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 106, symbol: 'P4_2bc', system: 'tetragonal', multiplicity: 8,
        generators: ['-y,x,z+1/2', 'x+1/2,-y+1/2,z'],
        centering: ['0,0,0'],
    },
    {
        number: 107, symbol: 'I4mm', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x,z', 'x,-y,z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 108, symbol: 'I4cm', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x,z', 'x,-y,z+1/2'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 109, symbol: 'I4_1md', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x+1/2,z+1/4', '-x,y,z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 110, symbol: 'I4_1cd', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x+1/2,z+1/4', '-x,y,z+1/2'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 111, symbol: 'P-42m', system: 'tetragonal', multiplicity: 8,
        generators: ['y,-x,-z', 'x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 112, symbol: 'P-42c', system: 'tetragonal', multiplicity: 8,
        generators: ['y,-x,-z', 'x,-y,-z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 113, symbol: 'P-42_1m', system: 'tetragonal', multiplicity: 8,
        generators: ['y,-x,-z', 'x+1/2,-y+1/2,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 114, symbol: 'P-42_1c', system: 'tetragonal', multiplicity: 8,
        generators: ['y,-x,-z', 'x+1/2,-y+1/2,-z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 115, symbol: 'P-4m2', system: 'tetragonal', multiplicity: 8,
        generators: ['y,-x,-z', 'y,x,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 116, symbol: 'P-4c2', system: 'tetragonal', multiplicity: 8,
        generators: ['y,-x,-z', 'y,x,-z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 117, symbol: 'P-4b2', system: 'tetragonal', multiplicity: 8,
        generators: ['y,-x,-z', 'y+1/2,x+1/2,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 118, symbol: 'P-4n2', system: 'tetragonal', multiplicity: 8,
        generators: ['y,-x,-z', 'y+1/2,x+1/2,-z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 119, symbol: 'I-4m2', system: 'tetragonal', multiplicity: 16,
        generators: ['y,-x,-z', 'y,x,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 120, symbol: 'I-4c2', system: 'tetragonal', multiplicity: 16,
        generators: ['y,-x,-z', 'y,x,-z+1/2'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 121, symbol: 'I-42m', system: 'tetragonal', multiplicity: 16,
        generators: ['y,-x,-z', 'x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 122, symbol: 'I-42d', system: 'tetragonal', multiplicity: 16,
        generators: ['y,-x,-z', 'x,-y+1/2,-z+1/4'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 123, symbol: 'P4/mmm', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x,z', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 124, symbol: 'P4/mcc', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x,z', 'x,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 125, symbol: 'P4/nbm', system: 'tetragonal', multiplicity: 16,
        generators: ['-y+1/2,x,z', 'x,-y+1/2,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 126, symbol: 'P4/nnc', system: 'tetragonal', multiplicity: 16,
        generators: ['-y+1/2,x,z', 'x,-y+1/2,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 127, symbol: 'P4/mbm', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x,z', 'x+1/2,-y+1/2,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 128, symbol: 'P4/mnc', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x,z', 'x+1/2,-y+1/2,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 129, symbol: 'P4/nmm', system: 'tetragonal', multiplicity: 16,
        generators: ['-y+1/2,x,z', 'x+1/2,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 130, symbol: 'P4/ncc', system: 'tetragonal', multiplicity: 16,
        generators: ['-y+1/2,x,z', 'x+1/2,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 131, symbol: 'P4_2/mmc', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x,z+1/2', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 132, symbol: 'P4_2/mcm', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x,z+1/2', 'x,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 133, symbol: 'P4_2/nbc', system: 'tetragonal', multiplicity: 16,
        generators: ['-y+1/2,x,z+1/2', 'x,-y+1/2,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 134, symbol: 'P4_2/nnm', system: 'tetragonal', multiplicity: 16,
        generators: ['-y+1/2,x,z+1/2', 'x,-y+1/2,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 135, symbol: 'P4_2/mbc', system: 'tetragonal', multiplicity: 16,
        generators: ['-y,x,z+1/2', 'x+1/2,-y+1/2,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 136, symbol: 'P4_2/mnm', system: 'tetragonal', multiplicity: 16,
        generators: ['-y+1/2,x+1/2,z+1/2', 'x+1/2,-y+1/2,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 137, symbol: 'P4_2/nmc', system: 'tetragonal', multiplicity: 16,
        generators: ['-y+1/2,x,z+1/2', 'x+1/2,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 138, symbol: 'P4_2/ncm', system: 'tetragonal', multiplicity: 16,
        generators: ['-y+1/2,x,z+1/2', 'x+1/2,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 139, symbol: 'I4/mmm', system: 'tetragonal', multiplicity: 32,
        generators: ['-y,x,z', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 140, symbol: 'I4/mcm', system: 'tetragonal', multiplicity: 32,
        generators: ['-y,x,z', 'x,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 141, symbol: 'I4_1/amd', system: 'tetragonal', multiplicity: 32,
        generators: ['-y+1/4,x+3/4,z+1/4', 'x,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 142, symbol: 'I4_1/acd', system: 'tetragonal', multiplicity: 32,
        generators: ['-y+1/4,x+3/4,z+1/4', 'x,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 143, symbol: 'P3', system: 'trigonal', multiplicity: 3,
        generators: ['-y,x-y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 144, symbol: 'P3_1', system: 'trigonal', multiplicity: 3,
        generators: ['-y,x-y,z+1/3'],
        centering: ['0,0,0'],
    },
    {
        number: 145, symbol: 'P3_2', system: 'trigonal', multiplicity: 3,
        generators: ['-y,x-y,z+2/3'],
        centering: ['0,0,0'],
    },
    {
        number: 146, symbol: 'R3', system: 'trigonal', multiplicity: 9,
        generators: ['-y,x-y,z'],
        centering: ['0,0,0', '2/3,1/3,1/3', '1/3,2/3,2/3'],
    },
    {
        number: 147, symbol: 'P-3', system: 'trigonal', multiplicity: 6,
        generators: ['-y,x-y,z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 148, symbol: 'R-3', system: 'trigonal', multiplicity: 18,
        generators: ['-y,x-y,z', '-x,-y,-z'],
        centering: ['0,0,0', '2/3,1/3,1/3', '1/3,2/3,2/3'],
    },
    {
        number: 149, symbol: 'P312', system: 'trigonal', multiplicity: 6,
        generators: ['-y,x-y,z', '-y,-x,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 150, symbol: 'P321', system: 'trigonal', multiplicity: 6,
        generators: ['-y,x-y,z', 'x-y,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 151, symbol: 'P3_112', system: 'trigonal', multiplicity: 6,
        generators: ['-y,x-y,z+1/3', 'x,x-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 152, symbol: 'P3_121', system: 'trigonal', multiplicity: 6,
        generators: ['-y,x-y,z+1/3', 'y,x,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 153, symbol: 'P3_212', system: 'trigonal', multiplicity: 6,
        generators: ['-y,x-y,z+2/3', 'x,x-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 154, symbol: 'P3_221', system: 'trigonal', multiplicity: 6,
        generators: ['-y,x-y,z+2/3', 'y,x,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 155, symbol: 'R32', system: 'trigonal', multiplicity: 18,
        generators: ['-y,x-y,z', 'y,x,-z'],
        centering: ['0,0,0', '2/3,1/3,1/3', '1/3,2/3,2/3'],
    },
    {
        number: 156, symbol: 'P3m1', system: 'trigonal', multiplicity: 6,
        generators: ['-y,x-y,z', '-y,-x,z'],
        centering: ['0,0,0'],
    },
    {
        number: 157, symbol: 'P31m', system: 'trigonal', multiplicity: 6,
        generators: ['-y,x-y,z', 'y,x,z'],
        centering: ['0,0,0'],
    },
    {
        number: 158, symbol: 'P3c1', system: 'trigonal', multiplicity: 6,
        generators: ['-y,x-y,z', '-y,-x,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 159, symbol: 'P31c', system: 'trigonal', multiplicity: 6,
        generators: ['-y,x-y,z', 'y,x,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 160, symbol: 'R3m', system: 'trigonal', multiplicity: 18,
        generators: ['-y,x-y,z', '-y,-x,z'],
        centering: ['0,0,0', '2/3,1/3,1/3', '1/3,2/3,2/3'],
    },
    {
        number: 161, symbol: 'R3c', system: 'trigonal', multiplicity: 18,
        generators: ['-y,x-y,z', '-y,-x,z+1/2'],
        centering: ['0,0,0', '2/3,1/3,1/3', '1/3,2/3,2/3'],
    },
    {
        number: 162, symbol: 'P-31m', system: 'trigonal', multiplicity: 12,
        generators: ['-y,x-y,z', '-y,-x,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 163, symbol: 'P-31c', system: 'trigonal', multiplicity: 12,
        generators: ['-y,x-y,z', '-y,-x,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 164, symbol: 'P-3m1', system: 'trigonal', multiplicity: 12,
        generators: ['-y,x-y,z', 'y,x,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 165, symbol: 'P-3c1', system: 'trigonal', multiplicity: 12,
        generators: ['-y,x-y,z', 'y,x,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 166, symbol: 'R-3m', system: 'trigonal', multiplicity: 36,
        generators: ['-y,x-y,z', 'y,x,-z', '-x,-y,-z'],
        centering: ['0,0,0', '2/3,1/3,1/3', '1/3,2/3,2/3'],
    },
    {
        number: 167, symbol: 'R-3c', system: 'trigonal', multiplicity: 36,
        generators: ['-y,x-y,z', 'y,x,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0', '2/3,1/3,1/3', '1/3,2/3,2/3'],
    },
    {
        number: 168, symbol: 'P6', system: 'hexagonal', multiplicity: 6,
        generators: ['x-y,x,z'],
        centering: ['0,0,0'],
    },
    {
        number: 169, symbol: 'P6_1', system: 'hexagonal', multiplicity: 6,
        generators: ['x-y,x,z+1/6'],
        centering: ['0,0,0'],
    },
    {
        number: 170, symbol: 'P6_5', system: 'hexagonal', multiplicity: 6,
        generators: ['x-y,x,z+5/6'],
        centering: ['0,0,0'],
    },
    {
        number: 171, symbol: 'P6_2', system: 'hexagonal', multiplicity: 6,
        generators: ['x-y,x,z+1/3'],
        centering: ['0,0,0'],
    },
    {
        number: 172, symbol: 'P6_4', system: 'hexagonal', multiplicity: 6,
        generators: ['x-y,x,z+2/3'],
        centering: ['0,0,0'],
    },
    {
        number: 173, symbol: 'P6_3', system: 'hexagonal', multiplicity: 6,
        generators: ['x-y,x,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 174, symbol: 'P-6', system: 'hexagonal', multiplicity: 6,
        generators: ['-x+y,-x,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 175, symbol: 'P6/m', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 176, symbol: 'P6_3/m', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 177, symbol: 'P622', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z', 'x-y,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 178, symbol: 'P6_122', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z+1/6', 'x-y,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 179, symbol: 'P6_522', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z+5/6', 'x-y,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 180, symbol: 'P6_222', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z+1/3', 'x-y,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 181, symbol: 'P6_422', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z+2/3', 'x-y,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 182, symbol: 'P6_322', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z+1/2', 'x-y,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 183, symbol: 'P6mm', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z', 'x-y,-y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 184, symbol: 'P6cc', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z', 'x-y,-y,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 185, symbol: 'P6_3cm', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z+1/2', 'x-y,-y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 186, symbol: 'P6_3mc', system: 'hexagonal', multiplicity: 12,
        generators: ['x-y,x,z+1/2', 'x-y,-y,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 187, symbol: 'P-6m2', system: 'hexagonal', multiplicity: 12,
        generators: ['-x+y,-x,-z', '-x+y,y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 188, symbol: 'P-6c2', system: 'hexagonal', multiplicity: 12,
        generators: ['-x+y,-x,-z+1/2', '-x+y,y,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 189, symbol: 'P-62m', system: 'hexagonal', multiplicity: 12,
        generators: ['-x+y,-x,-z', 'x-y,-y,z'],
        centering: ['0,0,0'],
    },
    {
        number: 190, symbol: 'P-62c', system: 'hexagonal', multiplicity: 12,
        generators: ['-x+y,-x,-z+1/2', 'x-y,-y,z+1/2'],
        centering: ['0,0,0'],
    },
    {
        number: 191, symbol: 'P6/mmm', system: 'hexagonal', multiplicity: 24,
        generators: ['x-y,x,z', 'x-y,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 192, symbol: 'P6/mcc', system: 'hexagonal', multiplicity: 24,
        generators: ['x-y,x,z', 'x-y,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 193, symbol: 'P6_3/mcm', system: 'hexagonal', multiplicity: 24,
        generators: ['x-y,x,z+1/2', 'x-y,-y,-z+1/2', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 194, symbol: 'P6_3/mmc', system: 'hexagonal', multiplicity: 24,
        generators: ['x-y,x,z+1/2', 'x-y,-y,-z', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 195, symbol: 'P23', system: 'cubic', multiplicity: 12,
        generators: ['-x,-y,z', 'z,x,y'],
        centering: ['0,0,0'],
    },
    {
        number: 196, symbol: 'F23', system: 'cubic', multiplicity: 48,
        generators: ['-x,-y,z', 'z,x,y'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 197, symbol: 'I23', system: 'cubic', multiplicity: 24,
        generators: ['-x,-y,z', 'z,x,y'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 198, symbol: 'P2_13', system: 'cubic', multiplicity: 12,
        generators: ['-x+1/2,-y,z+1/2', 'z,x,y'],
        centering: ['0,0,0'],
    },
    {
        number: 199, symbol: 'I2_13', system: 'cubic', multiplicity: 24,
        generators: ['-x+1/2,-y,z+1/2', 'z,x,y'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 200, symbol: 'Pm-3', system: 'cubic', multiplicity: 24,
        generators: ['-x,-y,z', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 201, symbol: 'Pn-3', system: 'cubic', multiplicity: 24,
        generators: ['-x+1/2,-y+1/2,z', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 202, symbol: 'Fm-3', system: 'cubic', multiplicity: 96,
        generators: ['-x,-y,z', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 203, symbol: 'Fd-3', system: 'cubic', multiplicity: 96,
        generators: ['-x+1/4,-y+1/4,z', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 204, symbol: 'Im-3', system: 'cubic', multiplicity: 48,
        generators: ['-x,-y,z', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 205, symbol: 'Pa-3', system: 'cubic', multiplicity: 24,
        generators: ['-x+1/2,-y,z+1/2', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 206, symbol: 'Ia-3', system: 'cubic', multiplicity: 48,
        generators: ['-x+1/2,-y,z+1/2', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 207, symbol: 'P432', system: 'cubic', multiplicity: 24,
        generators: ['-y,x,z', 'z,x,y'],
        centering: ['0,0,0'],
    },
    {
        number: 208, symbol: 'P4_232', system: 'cubic', multiplicity: 24,
        generators: ['-y+1/2,x+1/2,z+1/2', 'z,x,y'],
        centering: ['0,0,0'],
    },
    {
        number: 209, symbol: 'F432', system: 'cubic', multiplicity: 96,
        generators: ['-y,x,z', 'z,x,y'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 210, symbol: 'F4_132', system: 'cubic', multiplicity: 96,
        generators: ['-y+1/4,x+1/4,z+1/4', 'z,x,y'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 211, symbol: 'I432', system: 'cubic', multiplicity: 48,
        generators: ['-y,x,z', 'z,x,y'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 212, symbol: 'P4_332', system: 'cubic', multiplicity: 24,
        generators: ['-y+3/4,x+1/4,z+3/4', 'z,x,y'],
        centering: ['0,0,0'],
    },
    {
        number: 213, symbol: 'P4_132', system: 'cubic', multiplicity: 24,
        generators: ['-y+1/4,x+3/4,z+1/4', 'z,x,y'],
        centering: ['0,0,0'],
    },
    {
        number: 214, symbol: 'I4_132', system: 'cubic', multiplicity: 48,
        generators: ['-y+1/4,x+3/4,z+1/4', 'z,x,y'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 215, symbol: 'P-43m', system: 'cubic', multiplicity: 24,
        generators: ['y,-x,-z', 'z,x,y'],
        centering: ['0,0,0'],
    },
    {
        number: 216, symbol: 'F-43m', system: 'cubic', multiplicity: 96,
        generators: ['y,-x,-z', 'z,x,y'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 217, symbol: 'I-43m', system: 'cubic', multiplicity: 48,
        generators: ['y,-x,-z', 'z,x,y'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 218, symbol: 'P-43n', system: 'cubic', multiplicity: 24,
        generators: ['y+1/2,-x+1/2,-z+1/2', 'z,x,y'],
        centering: ['0,0,0'],
    },
    {
        number: 219, symbol: 'F-43c', system: 'cubic', multiplicity: 96,
        generators: ['y+1/2,-x,-z', 'z,x,y'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 220, symbol: 'I-43d', system: 'cubic', multiplicity: 48,
        generators: ['y+1/4,-x+3/4,-z+1/4', 'z,x,y'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 221, symbol: 'Pm-3m', system: 'cubic', multiplicity: 48,
        generators: ['-y,x,z', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 222, symbol: 'Pn-3n', system: 'cubic', multiplicity: 48,
        generators: ['-y+1/2,x,z', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 223, symbol: 'Pm-3n', system: 'cubic', multiplicity: 48,
        generators: ['-y+1/2,x+1/2,z+1/2', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 224, symbol: 'Pn-3m', system: 'cubic', multiplicity: 48,
        generators: ['-y,x+1/2,z+1/2', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0'],
    },
    {
        number: 225, symbol: 'Fm-3m', system: 'cubic', multiplicity: 192,
        generators: ['-y,x,z', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 226, symbol: 'Fm-3c', system: 'cubic', multiplicity: 192,
        generators: ['-y+1/2,x,z', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 227, symbol: 'Fd-3m', system: 'cubic', multiplicity: 192,
        generators: ['-y,x+1/4,z+1/4', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 228, symbol: 'Fd-3c', system: 'cubic', multiplicity: 192,
        generators: ['-y+1/2,x+1/4,z+1/4', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0', '0,1/2,1/2', '1/2,0,1/2', '1/2,1/2,0'],
    },
    {
        number: 229, symbol: 'Im-3m', system: 'cubic', multiplicity: 96,
        generators: ['-y,x,z', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
    {
        number: 230, symbol: 'Ia-3d', system: 'cubic', multiplicity: 96,
        generators: ['-y+1/4,x+3/4,z+1/4', 'z,x,y', '-x,-y,-z'],
        centering: ['0,0,0', '1/2,1/2,1/2'],
    },
];
