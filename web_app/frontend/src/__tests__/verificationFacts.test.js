// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The verification's facts catalogue (llm/copilot/verification/facts.js
// FACTS) against the real tools on the Demo run: every `facts` key a tool
// returns has its row, the row names that tool among `from`, and the value
// fits the row's shape (the dataShape notation); and every row is returned
// by every tool it names. A tool that changes its facts changes its row in
// the same diff, or this fails. The engines run in-thread
// (pcaKdeWorker's inProcessEngine, as src/openRunFolder.js runs them).

import { beforeAll, describe, expect, it } from 'vitest';
import { checkShape } from '../llm/copilot/dataShape';
import { buildRunTools } from '../llm/copilot/runTools';
import { executeToolCall } from '../llm/copilot/tools';
import { FACTS } from '../llm/copilot/verification';
import { buildCopilotTools } from '../copilotTools';
import { inProcessEngine } from '../workers/pcaKdeWorker';
import { loadDemoAssistantRun } from './modelVisible/helpers';

// The calls, with the Demo's own windows and sites; the cut Ta–Se window is
// refused, and still returns its count.
const CALLS = [
    ['assess_convergence', {}],
    ['find_gr_peaks', { pair: 'Ta-Se' }],
    ['find_gr_peaks', { pair: 'Ta-Ta' }],
    ['compute_bond_angles', { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.21, r12_max: 3.23 }],
    ['compute_bond_angles', { end1: 'Se', apex: 'Ta', end2: 'Se', r12_min: 2.3, r12_max: 2.7 }],
    ['analyze_displacement_directions', { site_ref: 1 }],
    ['analyze_displacement_directions', { element: 'Ga' }],
    ['analyze_clusters', { element: 'Ta', r_min: 2.65, r_max: 3.33 }],
    ['analyze_modulations', {}],
    ['detect_symmetry', { tolerance_A: 0.2 }],
    ['measure_symmetry_lowering', {}]
];

describe('the facts catalogue against the real tools (the Demo run)', () => {
    const results = [];

    beforeAll(async () => {
        const { run, rmc6fText } = await loadDemoAssistantRun();
        const request = inProcessEngine(rmc6fText);
        const tools = [
            ...buildRunTools(run, { cards: [] }),
            ...buildCopilotTools({ request, structure: run.structure, elements: Object.keys(run.structure.elementCounts || {}).sort() })
        ];
        for (const [name, args] of CALLS) {
            const tool = tools.find((candidate) => candidate.name === name);
            results.push({ name, args, ...(await executeToolCall(tool, JSON.stringify(args), { timeoutMs: 120000 })) });
        }
    }, 120000);

    it('every key a tool returns has its row, naming that tool, and fits its shape', () => {
        for (const { name, args, error, facts } of results) {
            expect(error ?? null, `${name} ${JSON.stringify(args)}`).toBeNull();
            expect(facts, `${name} ${JSON.stringify(args)} returns facts`).toBeTruthy();
            for (const [key, value] of Object.entries(facts)) {
                const row = FACTS.find((fact) => fact.key === key);
                expect(row, `${name}: facts.${key} has no row`).toBeDefined();
                expect(row.from, `${name} is not among facts.${key}'s producers`).toContain(name);
                expect(checkShape(value, row.shape), `${name} ${JSON.stringify(args)}: facts.${key}`).toEqual([]);
            }
        }
    });

    it('every row is returned, by every tool it names', () => {
        for (const row of FACTS) {
            for (const tool of row.from) {
                expect(results.some(({ name, facts }) => name === tool && facts && row.key in facts), `${tool} → facts.${row.key}`).toBe(true);
            }
        }
    });
});
