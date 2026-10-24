// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The cut rule on a run's natural windows (shellCuts.js; architecture
// candidate #6): one rule since the user's approval (2026-10-19), so the
// Workbench asks exactly where compute_bond_angles refuses, beside the
// coordination counts' reading the tool refused on before. Pinned on one dev
// case, where the switch changed the verdict and where it did not;
// `validation/run.js shell-cuts` prints every dev case and the Demo
// (docs/algorithms/copilot-validation.md).

import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES } from '../cases';
import { loadRun } from '../runFolder';
import { renderShellCuts, shellCutRows, shellCutSummary } from '../shellCuts';

describe('the cut rule on the natural windows (candidate #6; one rule since 2026-10-19)', () => {
    it('dev-jahn-teller: the Workbench asks where the tool refuses; the CuF6 short bonds are refused now, two slightly cut windows no longer', async () => {
        const item = DEV_CASES.find((entry) => entry.id === 'dev-jahn-teller');
        const run = await loadRun(buildCase(item), { ladder: false });
        const rows = await shellCutRows({ name: item.id, props: run.props, engine: run.engine, tools: run.tools });
        const row = (triplet, label) => rows.find((entry) => entry.triplet === triplet && entry.label === label);
        // The whole octahedron: no rule cuts it.
        expect(row('F–Cu–F', 'whole shell')).toMatchObject({ holds: 'whole', asks: false, refuses: false, before: false, modal: 6 });
        // Its four short bonds: 92 % of the Cu have exactly four, so the
        // counts kept them; they hold 4 of 6 F per Cu, so the CN margin cuts.
        expect(row('F–Cu–F', 'first peak')).toMatchObject({ holds: 'part', asks: true, refuses: true, before: false, modal: 4 });
        expect(row('F–Cu–F', 'first peak').modalShare).toBeGreaterThan(0.9);
        // Both rules cut the F–F first peak.
        expect(row('F–F–F', 'first peak')).toMatchObject({ asks: true, refuses: true, before: true, by: 'count-spread' });
        // A top 0.2 Å early that keeps the CN within the margin, its counts spread: refused before, not now.
        expect(row('Cu–K–Cu', 'cut 0.2 Å earlier')).toMatchObject({ asks: false, refuses: false, before: true, by: 'count-spread' });
        const summary = shellCutSummary(rows);
        expect(summary.windows + summary.empty.length).toBe(rows.length);
        expect(summary.disagree).toEqual([]);
        expect(summary.nowRefused.map((entry) => `${entry.triplet} ${entry.label}`)).toEqual(['F–Cu–F first peak']);
        expect(summary.noLongerRefused.map((entry) => `${entry.triplet} ${entry.label}`)).toEqual(['Cu–K–Cu cut 0.2 Å earlier', 'F–K–F cut 0.2 Å earlier']);
        const text = renderShellCuts(rows);
        expect(text).toMatch(/they disagree on 0\. Against the counts rule the tool refused on before 2026-10-19, 3 change verdict: 2 no longer refused, 1 refused now/);
        expect(text).toMatch(/\*\*The Workbench and compute_bond_angles disagree\*\*: none/);
    }, 60000);
});
