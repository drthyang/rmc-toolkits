// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Two refined translations of one rotation closer than the pass tolerance are one
// operation at that resolution, and only one is kept. The search kept the FIRST one found,
// in seed order. At the ladder's loose 1 Å a poor near-duplicate then shadowed the true
// operation: on a noisy P4_322 the ladder held a 4-fold at 0.8 Å residual 0.7 Å away from
// the 4_3 at 0.2 Å, showed P222_1 from 0.21 to 0.96 Å — where the headline, a separate
// pass at τ, reads P4_322 — and ended in a brick labelled "not a group".

import { describe, it, expect } from 'vitest';

import { spaceGroupAtTolerance, symmetryLadder } from '../symmetry.js';
import { SPACE_GROUP_FIXTURES, structureFor } from './fixtures/spaceGroups.js';
import { withNoise, shuffled } from './fixtures/symmetryStructures.js';

const p4322 = () => withNoise(structureFor(SPACE_GROUP_FIXTURES.find((f) => f.number === 95)), 0.03, 2);
const brickAt = (ladder, tol) => ladder.find((b) => tol >= b.from && tol < b.to);

describe('near-duplicate operations keep the best fit', () => {
    it('lets the ladder reach the group the headline finds', () => {
        const s = p4322();
        const ladder = symmetryLadder(s.A, s.basis, 1.0);
        for (const tol of [0.3, 0.5, 0.8]) {
            expect(spaceGroupAtTolerance(s.A, s.basis, tol).spaceGroup).toBe('P4_322');
            expect(brickAt(ladder, tol).spaceGroup).toBe('P4_322');
        }
        expect(ladder.map((b) => b.spaceGroup)).not.toContain('not a group');
    });

    it('does not depend on the order of the sites', () => {
        const s = p4322();
        const reference = JSON.stringify(symmetryLadder(s.A, s.basis, 1.0));
        for (const seed of [1, 2, 3]) {
            expect(JSON.stringify(symmetryLadder(s.A, shuffled(s.basis, seed), 1.0))).toBe(reference);
        }
    });
});
