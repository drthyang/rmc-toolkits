// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A number the model typed before its own reference is shown once in the dock
// (copilot/echoedReferences.js): gemma4:26b wrote "2.64 Å {{r1.…}}", which
// read "2.64 Å 2.64" once filled. Display only.

import { describe, expect, it } from 'vitest';
import { createReferenceStore, fillReferences } from '../copilot/references';
import { collapseEchoes, collapseFilledEchoes } from '../copilot/echoedReferences';
import { applyTurnEvent, newAssistantTurn } from '../copilot/turnState';
import { displayReferences } from '../components/referenceDisplay';

const store = () => {
    const references = createReferenceStore();
    references.put('r1', { first_shell: { mean_d_A: 2.637, cn: { Se_per_Ta: 6 } }, peaks_A: [{ r: 2.52 }] });
    return references;
};

describe('a number typed before its own reference', () => {
    it('drops the typed number and keeps the reference with its unit', () => {
        const raw = 'Mean distance: 2.64 Å {{r1.first_shell.mean_d_A}}; CN 6 {{r1.first_shell.cn.Se_per_Ta}} per Ta.';
        expect(collapseEchoes(raw, store())).toBe('Mean distance: {{r1.first_shell.mean_d_A}} Å; CN {{r1.first_shell.cn.Se_per_Ta}} per Ta.');
        expect(fillReferences(collapseEchoes(raw, store()), store()).text).toBe('Mean distance: 2.637 Å; CN 6 per Ta.');
    });

    it('leaves a different number, an unknown reference and plain text alone', () => {
        const raw = 'Shell at 3.10 Å {{r1.first_shell.mean_d_A}}, peak {{r1.peaks_A.0.r}}, missing 4 {{r1.nope}}.';
        expect(collapseEchoes(raw, store())).toBe(raw);
        expect(collapseEchoes('No references here: 2.64 Å.', store())).toBe('No references here: 2.64 Å.');
        expect(collapseEchoes(raw, null)).toBe(raw);
    });

    it('shows the filled answer\'s echoed numbers once', () => {
        const raw = 'Mean 2.64 Å {{r1.first_shell.mean_d_A}} and CN 6 {{r1.first_shell.cn.Se_per_Ta}}.';
        const filled = fillReferences(raw, store()).text;
        expect(filled).toBe('Mean 2.64 Å 2.637 and CN 6 6.');
        expect(collapseFilledEchoes(filled, raw, store())).toBe('Mean 2.637 Å and CN 6.');
    });

    it('the streaming text shows it once too', () => {
        expect(displayReferences('CN 6 {{r1.first_shell.cn.Se_per_Ta}} per Ta', { store: store(), streaming: true })).toBe('CN 6 per Ta');
    });

    it('the dock\'s ended turn shows the answer with each number once, and keeps the text the model reads', () => {
        const raw = 'Mean 2.64 Å {{r1.first_shell.mean_d_A}}.\nOutcome: done.';
        const events = [
            { type: 'round', round: 0 },
            { type: 'tool_start', id: 't1', name: 'find_gr_peaks', label: 'g(r) peaks', args: {} },
            { type: 'tool_end', id: 't1', name: 'find_gr_peaks', status: 'ok', args: {}, data: { first_shell: { mean_d_A: 2.637, cn: { Se_per_Ta: 6 } } }, checks: [], error: null, ms: 1, ref: 'r1', analysis: true },
            { type: 'round', round: 1 },
            { type: 'text', delta: raw },
            {
                type: 'done', text: 'Mean 2.64 Å 2.637.\nOutcome: done.', rawText: raw,
                outcome: { status: 'achieved', note: '1 analysis ran.', body: 'Mean 2.64 Å 2.637.', claimed: null, computed: true },
                stoppedBy: 'answer', audit: [], references: { written: 1, filled: 1, unknown: [] }, gates: []
            }
        ];
        const turn = events.reduce((current, event) => applyTurnEvent(current, event, 0), newAssistantTurn(0));
        expect(turn.items.filter((item) => item.kind === 'text').at(-1).text).toBe('Mean 2.637 Å.');
        expect(turn.text).toBe('Mean 2.64 Å 2.637.\nOutcome: done.');
    });
});
