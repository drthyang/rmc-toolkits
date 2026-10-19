// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The two readings of a window that cuts a coordination shell, side by side
// (shellCuts.js; architecture candidate #6): the Workbench's CN margin asks,
// compute_bond_angles' counts refuse. Pinned on one dev case, where they
// agree and where they do not; `validation/run.js shell-cuts` prints every dev
// case and the Demo (docs/algorithms/copilot-validation.md).

import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES } from '../cases';
import { loadRun } from '../runFolder';
import { renderShellCuts, shellCutRows, shellCutSummary } from '../shellCuts';

describe('the two cut readings side by side (candidate #6)', () => {
    it('dev-jahn-teller: the CuF6 short bonds are a sub-shell the counts keep and the CN margin asks about', async () => {
        const item = DEV_CASES.find((entry) => entry.id === 'dev-jahn-teller');
        const run = await loadRun(buildCase(item), { ladder: false });
        const rows = await shellCutRows({ name: item.id, props: run.props, engine: run.engine });
        const row = (triplet, label) => rows.find((entry) => entry.triplet === triplet && entry.label === label);
        // The whole octahedron: neither reading cuts it.
        expect(row('F–Cu–F', 'whole shell')).toMatchObject({ holds: 'whole', asks: false, refuses: false, modal: 6 });
        // Its four short bonds: 92 % of the Cu have exactly four, so the
        // counts keep them; they hold 4 of 6 F per Cu, so the margin asks.
        expect(row('F–Cu–F', 'first peak')).toMatchObject({ holds: 'part', asks: true, refuses: false, modal: 4 });
        expect(row('F–Cu–F', 'first peak').modalShare).toBeGreaterThan(0.9);
        // Both readings cut the F–F first peak.
        expect(row('F–F–F', 'first peak')).toMatchObject({ asks: true, refuses: true, by: 'count-spread' });
        // A top 0.2 Å early that keeps the CN within the margin, its counts spread: the tool alone refuses.
        expect(row('Cu–K–Cu', 'cut 0.2 Å earlier')).toMatchObject({ asks: false, refuses: true, by: 'count-spread' });
        const summary = shellCutSummary(rows);
        expect(summary.windows + summary.empty.length).toBe(rows.length);
        expect(summary.marginOnly.map((entry) => `${entry.triplet} ${entry.label}`)).toEqual(['F–Cu–F first peak']);
        expect(renderShellCuts(rows)).toMatch(/CN margin only\)\*\* \(1\)/);
    }, 60000);
});
