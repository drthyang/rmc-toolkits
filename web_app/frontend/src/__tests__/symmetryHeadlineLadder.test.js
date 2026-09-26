// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Detected SG card shows a headline group at the selected tolerance τ and a ladder of
// the groups over 0–1 Å. The ladder comes from ONE detection pass at 1 Å; the headline ran
// its own pass at τ. A pass pairs each seed within twice its tolerance while refining it,
// and merges translations of one rotation closer than its tolerance, so the two passes
// kept different seeds, landed in different local least-squares optima, and disagreed at
// some brick midpoints (P-6m2: ladder P3m1, headline Cm; a 4-site Pm: headline 'not a
// group'). The seeds also came from the FIRST site of the rarest element, so at a loose
// pairing radius the residuals — and the ladder's brick boundaries — moved with the order
// of the basis.

import { describe, it, expect } from 'vitest';

import { spaceGroupAtTolerance, symmetryLadder } from '../symmetry.js';
import { describeSymmetry } from '../symmetryModel.js';
import { SPACE_GROUP_FIXTURES, structureFor } from './fixtures/spaceGroups.js';
import { withNoise, shuffled } from './fixtures/symmetryStructures.js';
import { P6M2_NOISY, PM_NOISY } from './fixtures/headlineLadderCases.js';

const byNumber = new Map(SPACE_GROUP_FIXTURES.map((f) => [f.number, f]));
const card = ({ A, basis }, tol) => describeSymmetry({ latticeVectors: A, supercell: [1, 1, 1], basis }, tol);

// Every brick midpoint: the card's headline (describeSymmetry) names the brick's group, and
// is the finder's group at that tolerance in the ladder's pass (same operations, same fit).
// (A brick's operation count is its loosest rung's — merged rungs share a symbol — so the
// count is compared with the finder, not with the brick.)
function disagreements(st) {
    const out = [];
    for (const brick of symmetryLadder(st.A, st.basis, 1.0)) {
        const mid = (brick.from + brick.to) / 2;
        const headline = card(st, mid);
        const direct = spaceGroupAtTolerance(st.A, st.basis, mid, 1.0);
        if (headline.spaceGroup !== brick.spaceGroup) out.push(`at ${mid.toFixed(3)} Å: card ${headline.spaceGroup}, brick ${brick.spaceGroup}`);
        if (direct.spaceGroup !== headline.spaceGroup || direct.nSpace !== headline.nSpace || direct.maxResidual !== headline.maxResidual) {
            out.push(`at ${mid.toFixed(3)} Å: finder ${direct.spaceGroup} (${direct.nSpace}), card ${headline.spaceGroup} (${headline.nSpace})`);
        }
    }
    return out;
}

describe('the headline is the ladder rung at the selected tolerance', () => {
    it('on a noisy P-6m2 structure (ladder P3m1 at 0.15 Å)', () => {
        expect(symmetryLadder(P6M2_NOISY.A, P6M2_NOISY.basis, 1.0).map((b) => b.spaceGroup)).toContain('P3m1');
        expect(disagreements(P6M2_NOISY)).toEqual([]);
    });

    it('on a 4-site noisy Pm structure, never "not a group"', () => {
        expect(disagreements(PM_NOISY)).toEqual([]);
        for (const tol of [0.9, 0.95, 0.98]) expect(card(PM_NOISY, tol).spaceGroup).not.toBe('not a group');
    });

    it('on noisy trigonal and hexagonal fixtures', { timeout: 60000 }, () => {
        const out = [];
        for (const f of SPACE_GROUP_FIXTURES.filter((g) => g.number >= 143 && g.number <= 194)) {
            for (const d of disagreements(withNoise(structureFor(f), 0.03, f.number))) out.push(`${f.symbol}: ${d}`);
        }
        expect(out).toEqual([]);
    });
});

describe('the ladder does not depend on the order of the basis', () => {
    const bricks = (st) => symmetryLadder(st.A, st.basis, 1.0).map((b) => `${b.spaceGroup} ${b.from.toFixed(9)}–${b.to.toFixed(9)} (${b.nSpace})`);

    it('keeps the brick boundaries of a noisy P6_3/mmc for every shuffle', () => {
        const st = withNoise(structureFor(byNumber.get(194)), 0.03, 194);
        const reference = bricks(st);
        for (const seed of [201, 202, 203, 204, 205, 206]) {
            expect(bricks({ A: st.A, basis: shuffled(st.basis, seed) }), `seed ${seed}`).toEqual(reference);
        }
    });

    it('keeps the brick boundaries of noisy fixtures from every crystal system', { timeout: 60000 }, () => {
        for (const n of [14, 62, 63, 136, 140, 160, 166, 186, 205, 216, 221]) {
            const st = withNoise(structureFor(byNumber.get(n)), 0.03, n);
            const reference = bricks(st);
            for (let s = 1; s <= 3; s++) {
                expect(bricks({ A: st.A, basis: shuffled(st.basis, s * 7 + n) }), `No. ${n}, shuffle ${s}`).toEqual(reference);
            }
        }
    });
});
