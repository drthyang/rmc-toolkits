// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, it, expect } from 'vitest';

import {
    spaceGroupAtTolerance,
    symmetryLadder,
    findSpaceGroupOps,
    siteOrbits,
    classifyRotation,
    pointGroupOf,
} from '../symmetry.js';
import {
    intrinsicTranslation,
    characteristicDirection,
    classifyElement,
    glideLetter,
    rotationOrder,
    hmSymbolCandidates,
    hmSymbolInStandardSetting,
    coversAllElements,
} from '../spaceGroupSymbol.js';
import { SPACE_GROUPS, spaceGroupNumber, isStandardSymbol, canonicalSymbol } from '../spaceGroupTable.js';
import {
    SPACE_GROUP_FIXTURES,
    closeGroup,
    structureFor,
    parseTriplet,
    cellVectors,
    orbitBasis,
} from './fixtures/spaceGroups.js';

const cubic = (a) => [[a, 0, 0], [0, a, 0], [0, 0, a]];
const hexagonal = (a, c) => [[a, 0, 0], [-a / 2, (a * Math.sqrt(3)) / 2, 0], [0, 0, c]];
const at = (el, x, y, z) => ({ el, frac: [x, y, z] });
const wrap = (v) => ((v % 1) + 1) % 1;

describe('space-group table', () => {
    it('holds all 230 groups with unique numbers', () => {
        expect(SPACE_GROUPS).toHaveLength(230);
        expect(new Set(SPACE_GROUPS.map((g) => g[0])).size).toBe(230);
        expect(Math.min(...SPACE_GROUPS.map((g) => g[0]))).toBe(1);
        expect(Math.max(...SPACE_GROUPS.map((g) => g[0]))).toBe(230);
    });

    it('resolves symbols to numbers and rejects non-symbols', () => {
        expect(spaceGroupNumber('Pnma')).toBe(62);
        expect(spaceGroupNumber('I4/mcm')).toBe(140);
        expect(spaceGroupNumber('Fd-3m')).toBe(227);
        expect(spaceGroupNumber('P6_3/mmc')).toBe(194);
        expect(spaceGroupNumber('Pmmm')).toBe(47);
        expect(spaceGroupNumber('I4/mbm')).toBeNull();
        expect(isStandardSymbol('Pnma')).toBe(true);
        expect(isStandardSymbol('Pbnm')).toBe(false);
    });

    it('accepts pre-2002 spellings and canonicalises them', () => {
        expect(spaceGroupNumber('Cmca')).toBe(64);
        expect(canonicalSymbol('Cmca')).toBe('Cmce');
        expect(canonicalSymbol('Pnma')).toBe('Pnma');
    });
});

describe('symmetry elements', () => {
    it('reads the order of a rotation matrix', () => {
        expect(rotationOrder([[1, 0, 0], [0, 1, 0], [0, 0, 1]])).toBe(1);
        expect(rotationOrder([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])).toBe(2);
        expect(rotationOrder([[0, -1, 0], [1, 0, 0], [0, 0, 1]])).toBe(4);
        expect(rotationOrder([[0, -1, 0], [1, -1, 0], [0, 0, 1]])).toBe(3);
    });

    it('finds the axis of a rotation and the normal of a reflection', () => {
        expect(characteristicDirection([[0, -1, 0], [1, 0, 0], [0, 0, 1]])).toEqual([0, 0, 1]);
        expect(characteristicDirection([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])).toEqual([1, 0, 0]);
        // the identity and inversion sit on no direction
        expect(characteristicDirection([[1, 0, 0], [0, 1, 0], [0, 0, 1]])).toBeNull();
        expect(characteristicDirection([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])).toBeNull();
    });

    it('splits an operation into its intrinsic translation', () => {
        // 2-fold along c with t = (1/2, 0, 1/2): the in-axis half is intrinsic (a 2_1
        // screw), the in-plane half only moves where the axis sits.
        const twofold = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]];
        expect(intrinsicTranslation(twofold, [0.5, 0, 0.5])).toEqual([0, 0, 0.5]);
        // a pure rotation has none
        expect(intrinsicTranslation(twofold, [0.5, 0.5, 0])).toEqual([0, 0, 0]);
    });

    it('names glide planes from the intrinsic translation', () => {
        expect(glideLetter([0, 0, 0])).toBe('m');
        expect(glideLetter([0.5, 0, 0])).toBe('a');
        expect(glideLetter([0, 0, 0.5])).toBe('c');
        expect(glideLetter([0.5, 0.5, 0])).toBe('n');
        expect(glideLetter([0.25, 0.25, 0])).toBe('d');
    });

    it('tells a screw axis from a rotation', () => {
        const fourfold = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        expect(classifyElement(fourfold, [0, 0, 0])).toMatchObject({ kind: 'rotation', label: '4' });
        expect(classifyElement(fourfold, [0, 0, 0.5])).toMatchObject({ kind: 'screw', label: '4_2' });
        const mirror = [[-1, 0, 0], [0, 1, 0], [0, 0, 1]];
        expect(classifyElement(mirror, [0, 0, 0])).toMatchObject({ kind: 'mirror', label: 'm' });
        expect(classifyElement(mirror, [0, 0, 0.5])).toMatchObject({ kind: 'glide', label: 'c' });
    });

    it('distinguishes the enantiomorphic screws 4_1 and 4_3', () => {
        const fourfold = [[0, -1, 0], [1, 0, 0], [0, 0, 1]];
        const inverse = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];
        expect(classifyElement(fourfold, [0, 0, 0.25]).label).toBe('4_1');
        expect(classifyElement(inverse, [0, 0, 0.25]).label).toBe('4_3');
    });
});

describe('settings the symbol positions cannot describe', () => {
    // A 3-fold about [111] with a mirror through it: crystal class 3m, but written on
    // CUBIC axes, which is how a trigonal subgroup turns up part-way along the tolerance
    // ladder of a cubic structure. The hexagonal-axes positions expect the 3-fold on
    // [001], so nothing can be placed.
    const threefold111 = [[0, 0, 1], [1, 0, 0], [0, 1, 0]];
    const mirror110 = [[0, 1, 0], [1, 0, 0], [0, 0, 1]];
    const cubicAxes3m = [
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]], threefold111,
        [[0, 1, 0], [0, 0, 1], [1, 0, 0]], mirror110,
        [[0, 0, 1], [0, 1, 0], [1, 0, 0]], [[1, 0, 0], [0, 0, 1], [0, 1, 0]],
    ].map((R) => ({ R, t: [0, 0, 0] }));

    it('spots operations it cannot place', () => {
        expect(coversAllElements(cubicAxes3m, '3m')).toBe(false);
    });

    it('accepts the same class on hexagonal axes', () => {
        const p3m1 = closeGroup(SPACE_GROUP_FIXTURES.find((f) => f.number === 156));
        expect(coversAllElements(p3m1, '3m')).toBe(true);
    });

    it('falls back to the crystal class rather than inventing a symbol', () => {
        const found = hmSymbolInStandardSetting(cubicAxes3m, 'P', '3m', isStandardSymbol);
        expect(found).toMatchObject({ symbol: 'P3m', standard: false, placed: false });
    });

    it('keeps a correct symbol that is merely in a non-standard setting', () => {
        // Pn is a real n-glide, just not the standard spelling of #7 (Pc) — worth
        // showing as-is, unlike an unplaceable set.
        const nGlide = [
            { R: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0, 0, 0] },
            { R: [[1, 0, 0], [0, -1, 0], [0, 0, 1]], t: [0.5, 0, 0.5] },
        ];
        const found = hmSymbolInStandardSetting(nGlide, 'P', 'm', isStandardSymbol);
        expect(found).toMatchObject({ symbol: 'Pn', standard: false, placed: true });
    });
});

describe('fixture integrity', () => {
    it.each(SPACE_GROUP_FIXTURES)('$symbol (#$number) closes to $multiplicity operations', (fixture) => {
        expect(closeGroup(fixture)).toHaveLength(fixture.multiplicity);
    });
});

// Two pairs of groups this finder cannot separate, and the partner it settles on.
//
// I222/I2_12_12_1 and I23/I2_13 contain the SAME element types along the same directions:
// the I centring turns every pure 2-fold into a 2_1 screw half a cell away and vice versa,
// so both members of each pair have both. What distinguishes them is where those axes sit
// RELATIVE to one another, which naming from element types alone cannot see — telling them
// apart needs origin-aware matching against the tabulated groups, which this finder does
// not attempt. Pinned here rather than skipped so the limitation stays visible and any
// change in it fails loudly.
const LOCATION_DEGENERATE = new Map([[24, 'I222'], [199, 'I23']]);

describe('Hermann–Mauguin symbols from exact operations', () => {
    it.each(SPACE_GROUP_FIXTURES)('names $symbol (#$number)', (fixture) => {
        const ops = closeGroup(fixture);
        const pointGroup = SPACE_GROUPS.find((g) => g[0] === fixture.number)[2];
        const candidates = hmSymbolCandidates(ops, fixture.symbol[0], pointGroup);
        // The table picks the first candidate it recognises — assert that is the truth,
        // not merely that the truth appears somewhere in the list.
        expect(candidates.find(isStandardSymbol))
            .toBe(LOCATION_DEGENERATE.get(fixture.number) ?? fixture.symbol);
    });
});

describe('space-group detection end to end', () => {
    it.each(SPACE_GROUP_FIXTURES)('recovers $symbol (#$number) from its atoms', (fixture) => {
        const { A, basis } = structureFor(fixture);
        const found = spaceGroupAtTolerance(A, basis, 0.05);
        expect(found.nSpace).toBe(fixture.multiplicity);
        const degenerate = LOCATION_DEGENERATE.get(fixture.number);
        if (degenerate) {
            expect(found.spaceGroup).toBe(degenerate);
            return;
        }
        expect(found.spaceGroup).toBe(fixture.symbol);
        expect(found.spaceGroupNumber).toBe(fixture.number);
    });
});

describe('known structures', () => {
    // The cases that motivated this work: before screw and glide components were read,
    // every non-symmorphic group here was reported as its symmorphic parent — diamond as
    // Fm-3m, the tilted perovskite as I4/mmm, hcp as P6/mmm.
    const fcc = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]];

    it('names rocksalt Fm-3m', () => {
        const basis = [];
        for (const [x, y, z] of fcc) {
            basis.push(at('Na', x, y, z));
            basis.push(at('Cl', wrap(x + 0.5), y, z));
        }
        expect(spaceGroupAtTolerance(cubic(5.64), basis, 0.05)).toMatchObject({
            spaceGroup: 'Fm-3m', spaceGroupNumber: 225, nSpace: 192,
        });
    });

    it('names the cubic perovskite Pm-3m', () => {
        const basis = [
            at('Sr', 0, 0, 0), at('Ti', 0.5, 0.5, 0.5),
            at('O', 0.5, 0.5, 0), at('O', 0.5, 0, 0.5), at('O', 0, 0.5, 0.5),
        ];
        expect(spaceGroupAtTolerance(cubic(3.905), basis, 0.05)).toMatchObject({
            spaceGroup: 'Pm-3m', spaceGroupNumber: 221, nSpace: 48,
        });
    });

    it('names diamond Fd-3m, not its symmorphic parent Fm-3m', () => {
        const basis = [];
        for (const [x, y, z] of fcc) {
            basis.push(at('C', x, y, z));
            basis.push(at('C', wrap(x + 0.25), wrap(y + 0.25), wrap(z + 0.25)));
        }
        expect(spaceGroupAtTolerance(cubic(3.567), basis, 0.05)).toMatchObject({
            spaceGroup: 'Fd-3m', spaceGroupNumber: 227, nSpace: 192,
        });
    });

    it('names hcp P6_3/mmc, not P6/mmm', () => {
        const basis = [at('Mg', 1 / 3, 2 / 3, 0.25), at('Mg', 2 / 3, 1 / 3, 0.75)];
        expect(spaceGroupAtTolerance(hexagonal(3.21, 5.21), basis, 0.05)).toMatchObject({
            spaceGroup: 'P6_3/mmc', spaceGroupNumber: 194, nSpace: 24,
        });
    });

    it('names an a⁰a⁰c⁻ tilted perovskite I4/mcm, not I4/mmm', () => {
        const tilt = 0.03;
        const basis = [];
        for (const [cx, cy, cz] of [[0, 0, 0], [0.5, 0.5, 0.5]]) {
            basis.push(at('Ti', wrap(cx), wrap(cy), wrap(cz)));
            basis.push(at('Ti', wrap(cx), wrap(cy), wrap(cz + 0.5)));
            basis.push(at('Sr', wrap(cx), wrap(cy + 0.5), wrap(cz + 0.25)));
            basis.push(at('Sr', wrap(cx + 0.5), wrap(cy), wrap(cz + 0.25)));
            basis.push(at('O', wrap(cx), wrap(cy), wrap(cz + 0.25)));
            basis.push(at('O', wrap(cx), wrap(cy), wrap(cz + 0.75)));
            basis.push(at('O', wrap(cx + 0.25 + tilt), wrap(cy + 0.25 - tilt), wrap(cz)));
            basis.push(at('O', wrap(cx + 0.75 - tilt), wrap(cy + 0.75 + tilt), wrap(cz)));
            basis.push(at('O', wrap(cx + 0.75 + tilt), wrap(cy + 0.25 + tilt), wrap(cz)));
            basis.push(at('O', wrap(cx + 0.25 - tilt), wrap(cy + 0.75 - tilt), wrap(cz)));
        }
        expect(spaceGroupAtTolerance([[5.52, 0, 0], [0, 5.52, 0], [0, 0, 7.81]], basis, 0.05)).toMatchObject({
            spaceGroup: 'I4/mcm', spaceGroupNumber: 140, nSpace: 32,
        });
    });
});

describe('centering', () => {
    it('detects rhombohedral R centering on hexagonal axes', () => {
        const basis = [at('X', 0, 0, 0), at('X', 1 / 3, 2 / 3, 2 / 3), at('X', 2 / 3, 1 / 3, 1 / 3)];
        const found = spaceGroupAtTolerance(hexagonal(5, 15), basis, 0.05);
        expect(found.centering).toBe('R');
        expect(found).toMatchObject({ spaceGroup: 'R-3m', spaceGroupNumber: 166, nSpace: 36 });
    });

    it('detects the reverse R setting too', () => {
        const basis = [at('X', 0, 0, 0), at('X', 1 / 3, 2 / 3, 1 / 3), at('X', 2 / 3, 1 / 3, 2 / 3)];
        expect(spaceGroupAtTolerance(hexagonal(5, 15), basis, 0.05).centering).toBe('R');
    });

    it('reads F and I centering', () => {
        const fcc = [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]];
        expect(spaceGroupAtTolerance(cubic(4), fcc.map(([x, y, z]) => at('X', x, y, z)), 0.05).centering).toBe('F');
        const bcc = [at('X', 0, 0, 0), at('X', 0.5, 0.5, 0.5)];
        expect(spaceGroupAtTolerance(cubic(3), bcc, 0.05).centering).toBe('I');
    });
});

describe('non-standard settings', () => {
    it('reports Pnma for a structure given in the Pbnm setting', () => {
        const pnma = SPACE_GROUP_FIXTURES.find((f) => f.number === 62);
        const ops = closeGroup(pnma);
        const basis = orbitBasis(ops).map((s) => ({ el: s.el, frac: [s.frac[1], s.frac[2], s.frac[0]] }));
        // Pnma → Pbnm is the cyclic axis permutation (a, b, c) → (b, c, a).
        const A = cellVectors(5.3, 6.7, 8.1, 90, 90, 90);
        expect(spaceGroupAtTolerance(A, basis, 0.05)).toMatchObject({
            spaceGroup: 'Pnma', spaceGroupNumber: 62, nSpace: 8,
        });
    });
});

describe('tolerance ladder', () => {
    it('climbs monotonically to the full group as the tolerance loosens', () => {
        // A cubic perovskite with one oxygen nudged off its site: tight tolerance sees
        // the broken symmetry, loose tolerance recovers Pm-3m.
        const basis = [
            at('Sr', 0, 0, 0), at('Ti', 0.5, 0.5, 0.5),
            at('O', 0.5, 0.5, 0.02), at('O', 0.5, 0, 0.5), at('O', 0, 0.5, 0.5),
        ];
        const ladder = symmetryLadder(cubic(4), basis, 0.6);
        expect(ladder.length).toBeGreaterThan(1);
        expect(ladder[0].from).toBe(0);
        expect(ladder[ladder.length - 1].spaceGroup).toBe('Pm-3m');
        // rungs are contiguous and operation counts never fall as tolerance grows
        for (let i = 1; i < ladder.length; i += 1) {
            expect(ladder[i].from).toBeCloseTo(ladder[i - 1].to, 9);
            expect(ladder[i].nSpace).toBeGreaterThanOrEqual(ladder[i - 1].nSpace);
        }
    });

    it('every rung is a real space group', () => {
        const basis = [
            at('Sr', 0, 0, 0), at('Ti', 0.5, 0.5, 0.5),
            at('O', 0.5, 0.5, 0.02), at('O', 0.5, 0, 0.5), at('O', 0, 0.5, 0.5),
        ];
        for (const rung of symmetryLadder(cubic(4), basis, 0.6)) {
            expect(isStandardSymbol(rung.spaceGroup)).toBe(true);
            expect(rung.spaceGroupNumber).not.toBeNull();
        }
    });
});

describe('site orbits', () => {
    it('splits the perovskite basis into Sr, Ti and O orbits', () => {
        const basis = [
            at('Sr', 0, 0, 0), at('Ti', 0.5, 0.5, 0.5),
            at('O', 0.5, 0.5, 0), at('O', 0.5, 0, 0.5), at('O', 0, 0.5, 0.5),
        ];
        const A = cubic(3.905);
        const { ops } = spaceGroupAtTolerance(A, basis, 0.05);
        const orbits = siteOrbits(A, basis, ops, 0.05);
        expect(orbits.map((o) => [o.element, o.size, o.site])).toEqual([
            ['O', 3, '4/mmm'], ['Sr', 1, 'm-3m'], ['Ti', 1, 'm-3m'],
        ]);
    });
});

describe('rotation and point-group classification', () => {
    it('reads rotation type from determinant and trace', () => {
        expect(classifyRotation([[1, 0, 0], [0, 1, 0], [0, 0, 1]])).toBe('1');
        expect(classifyRotation([[-1, 0, 0], [0, -1, 0], [0, 0, 1]])).toBe('2');
        expect(classifyRotation([[0, -1, 0], [1, 0, 0], [0, 0, 1]])).toBe('4');
        expect(classifyRotation([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])).toBe('-1');
        expect(classifyRotation([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])).toBe('m');
    });

    it('derives the crystal class of each fixture group', { timeout: 60000 }, () => {
        for (const fixture of SPACE_GROUP_FIXTURES) {
            const rotations = [...new Map(closeGroup(fixture).map((o) => [o.R.flat().join(','), o.R])).values()];
            const expected = SPACE_GROUPS.find((g) => g[0] === fixture.number)[2];
            expect(pointGroupOf(rotations), fixture.symbol).toBe(expected);
        }
    });
});

describe('findSpaceGroupOps', () => {
    it('returns an empty group for an empty basis', () => {
        expect(findSpaceGroupOps(cubic(4), [], 0.1)).toMatchObject({ nSpace: 0, spaceGroup: 'P1' });
    });

    it('gives P-1 for a lone atom — a one-atom basis is always centrosymmetric', () => {
        const A = cellVectors(6.1, 7.3, 8.7, 81, 86, 94);
        const found = spaceGroupAtTolerance(A, [at('X', 0.137, 0.213, 0.061)], 0.05);
        expect(found).toMatchObject({ spaceGroup: 'P-1', spaceGroupNumber: 2, nSpace: 2 });
    });

    it('gives P1 once two elements sit on unrelated general positions', () => {
        const A = cellVectors(6.1, 7.3, 8.7, 81, 86, 94);
        // Inversion would have to be centred on the X atom and on the Y atom at once.
        const basis = [at('X', 0.137, 0.213, 0.061), at('Y', 0.41, 0.72, 0.23)];
        expect(spaceGroupAtTolerance(A, basis, 0.05)).toMatchObject({
            spaceGroup: 'P1', spaceGroupNumber: 1, nSpace: 1,
        });
    });

    it('parses ITA coordinate triplets', () => {
        expect(parseTriplet('-x+1/2,-y,z+1/2')).toEqual({
            R: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]], t: [0.5, 0, 0.5],
        });
        expect(parseTriplet('x-y,x,z')).toEqual({
            R: [[1, -1, 0], [1, 0, 0], [0, 0, 1]], t: [0, 0, 0],
        });
    });
});
