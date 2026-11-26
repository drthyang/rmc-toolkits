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

    describe('with a short unit or label between the number and its reference', () => {
        const labelled = () => {
            const references = store();
            references.put('r2', { cn: { Se_per_Ta: 3.25 } });
            references.put('r5', { angle: { mean_deg2: 127 } });
            return references;
        };

        it('shows "3.25 Se per Ta" and "127 deg²" once, the label after the value', () => {
            const raw = 'CN 3.25 Se per Ta {{r2.cn.Se_per_Ta}}; spread 127 deg² {{r5.angle.mean_deg2}}.';
            expect(collapseEchoes(raw, labelled())).toBe('CN {{r2.cn.Se_per_Ta}} Se per Ta; spread {{r5.angle.mean_deg2}} deg².');
            const filled = fillReferences(raw, labelled()).text;
            expect(filled).toBe('CN 3.25 Se per Ta 3.25; spread 127 deg² 127.');
            expect(collapseFilledEchoes(filled, raw, labelled())).toBe('CN 3.25 Se per Ta; spread 127 deg².');
        });

        it('still collapses a plain unit', () => {
            expect(fillReferences(collapseEchoes('d = 2.64 Å {{r1.first_shell.mean_d_A}}', labelled()), labelled()).text).toBe('d = 2.637 Å');
        });

        it('keeps a mismatch visible', () => {
            const raw = 'CN 3.20 Se per Ta {{r2.cn.Se_per_Ta}}.';
            expect(collapseEchoes(raw, labelled())).toBe(raw);
            const filled = fillReferences(raw, labelled()).text;
            expect(collapseFilledEchoes(filled, raw, labelled())).toBe('CN 3.20 Se per Ta 3.25.');
        });

        it('leaves a clause, a long label, a line break and digits between them alone', () => {
            for (const raw of [
                '6 neighbours, of which {{r1.first_shell.cn.Se_per_Ta}} are Se.',
                '6 neighbours of which most {{r1.first_shell.cn.Se_per_Ta}}.',
                '6 extraordinarilylongwordlabel {{r1.first_shell.cn.Se_per_Ta}}.',
                '6 Se\n{{r1.first_shell.cn.Se_per_Ta}}.',
                '6 (Se) {{r1.first_shell.cn.Se_per_Ta}}.',
                '6 Se per Ta2 {{r1.first_shell.cn.Se_per_Ta}}.'
            ]) {
                expect(collapseEchoes(raw, labelled())).toBe(raw);
            }
        });

        it('leaves a number without its reference and a reference to another value alone', () => {
            expect(collapseEchoes('CN 3.25 Se per Ta.', labelled())).toBe('CN 3.25 Se per Ta.');
            const other = 'CN 3.25 Se per Ta {{r1.first_shell.mean_d_A}}.';
            expect(collapseEchoes(other, labelled())).toBe(other);
        });

        it('the dock\'s streaming text shows it once', () => {
            expect(displayReferences('CN 3.25 Se per Ta {{r2.cn.Se_per_Ta}}', { store: labelled(), streaming: true })).toBe('CN 3.25 Se per Ta');
        });
    });
    describe('with a thousands separator in the typed number', () => {
        const grouped = () => {
            const references = store();
            references.put('r3', { n_atoms: 215252, volume_A3: 1234.5, small: 2, tail: 234, pair: 34 });
            return references;
        };

        it('shows "215,252" and its reference 215252 once', () => {
            const raw = 'The box holds 215,252 {{r3.n_atoms}} atoms.';
            expect(collapseEchoes(raw, grouped())).toBe('The box holds {{r3.n_atoms}} atoms.');
            const filled = fillReferences(raw, grouped()).text;
            expect(filled).toBe('The box holds 215,252 215252 atoms.');
            expect(collapseFilledEchoes(filled, raw, grouped())).toBe('The box holds 215252 atoms.');
            expect(displayReferences('215,252 atoms {{r3.n_atoms}}', { store: grouped(), streaming: true })).toBe('215252 atoms');
        });

        it('reads a decimal with grouping, with its unit', () => {
            const raw = 'Volume 1,234.5 Å³ {{r3.volume_A3}}.';
            expect(collapseEchoes(raw, grouped())).toBe('Volume {{r3.volume_A3}} Å³.');
            expect(collapseFilledEchoes(fillReferences(raw, grouped()).text, raw, grouped())).toBe('Volume 1234.5 Å³.');
        });

        it('never merges a comma list into one number', () => {
            const references = grouped();
            references.put('r4', { a: 12, b: 1234, c: 12345678 });
            for (const raw of [
                'Shells 1,2 {{r4.a}}.',
                'Shells 12,34 {{r4.b}}.',
                'Counts 12345,678 {{r4.c}}.',
                'Counts 1,2345 {{r4.a}}.'
            ]) {
                expect(collapseEchoes(raw, references)).toBe(raw);
            }
            // The last item of a list is still its own echo, as before.
            expect(collapseEchoes('Shells 1,2 {{r3.small}}.', references)).toBe('Shells 1,{{r3.small}}.');
            expect(collapseEchoes('Shells 12,34 {{r3.pair}}.', references)).toBe('Shells 12,{{r3.pair}}.');
            expect(collapseEchoes('Shells 1,234 {{r3.tail}}.', references)).toBe('Shells 1,{{r3.tail}}.');
        });

        it('keeps a grouped number that differs from the engine\'s value visible', () => {
            const raw = 'The box holds 215,253 {{r3.n_atoms}} atoms.';
            expect(collapseEchoes(raw, grouped())).toBe(raw);
            expect(collapseFilledEchoes(fillReferences(raw, grouped()).text, raw, grouped())).toBe('The box holds 215,253 215252 atoms.');
        });
    });
});
