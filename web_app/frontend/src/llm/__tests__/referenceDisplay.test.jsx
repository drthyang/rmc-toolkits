// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The dock never shows a raw {{…}} (components/referenceDisplay.js): filled
// from the turn's results while the text streams, held back while one is
// still being written, else marked "number missing". Through MarkdownMessage
// and the turn view (CopilotView's AssistantTurn) as the dock renders them.

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import MarkdownMessage from '../components/MarkdownMessage';
import { AssistantTurn } from '../components/CopilotView';
import { displayReferences, MISSING_LABEL, plainReferences } from '../components/referenceDisplay';
import { createReferenceStore } from '../copilot/references';
import { applyTurnEvent, newAssistantTurn } from '../copilot/turnState';

// polyhedral_distortion's average structure on the Demo (r5), as the model read it.
const R5 = { average_structure: { angle_variance_deg2: { mean: 111, sd: 3.1 }, trans_angle_deg: { mean: 159.6, sd: 0.3 } } };
const storeWith = () => {
    const store = createReferenceStore({ context: { symmetry: { space_group: 'F-43m' } } });
    store.put('r5', R5);
    return store;
};
const USER = 'has σ²θ = {{r5.data.average_structure.angle_variance_deg2.mean}} deg² and a trans angle of '
    + '{{r5.data.average_structure.trans_angle_deg.mean}}°,';

describe('displayReferences', () => {
    it('fills what the turn\'s results fill, the user\'s sentence included', () => {
        expect(displayReferences(USER, { store: storeWith() })).toBe('has σ²θ = 111 deg² and a trans angle of 159.6°,');
        expect(displayReferences('In {{ctx.symmetry.space_group}}.', { store: storeWith() })).toBe('In F-43m.');
        expect(displayReferences('No braces.', { store: storeWith() })).toBe('No braces.');
    });

    it('marks any other reference, in prose, math and code, and leaves other braces alone', () => {
        const out = displayReferences('σ²θ {{r5.data.average_structure.angle_variance_deg2}} deg², $x = {{r9.x}}$, `{{r9.y}}`, $a^{{2}}$.', { store: storeWith() });
        expect(out).not.toMatch(/\{\{\s*r\d/);
        expect(out).toBe(`σ²θ [${MISSING_LABEL}](ref:r5.data.average_structure.angle_variance_deg2) deg², `
            + `$x = \\text{[${MISSING_LABEL}]}$, \`[${MISSING_LABEL}]\`, $a^{{2}}$.`);
        expect(displayReferences('{{r1.x}} with no store.')).toBe(`[${MISSING_LABEL}](ref:r1.x) with no store.`);
    });

    it('while streaming, holds back a reference still being written, and fills the complete ones', () => {
        const store = storeWith();
        for (const tail of ['{{', '{{r', '{{r5.data.aver', '{{ r5.data.average_structure.trans_angle_deg.mean }']) {
            expect(displayReferences(`σ²θ = {{r5.average_structure.angle_variance_deg2.mean}} deg², trans ${tail}`, { store, streaming: true }))
                .toBe('σ²θ = 111 deg², trans ');
        }
        // Done streaming: an unclosed one is no reference, shown as written.
        expect(displayReferences('trans {{r5.data.aver', { store })).toBe('trans {{r5.data.aver');
    });

    it('the screen reader\'s line reads a reference left in it as the marker\'s words', () => {
        expect(plainReferences('σ²θ {{r5.x}} deg², $a^{{2}}$')).toBe(`σ²θ ${MISSING_LABEL} deg², $a^{{2}}$`);
    });
});

describe('the dock renders no raw reference', () => {
    const render = (text, options = {}) => renderToStaticMarkup(<MarkdownMessage {...options}>{text}</MarkdownMessage>);

    it('MarkdownMessage: filled numbers, and a marked span (never a link) naming what is missing', () => {
        const html = render(`${USER} DI {{r5.distortion_index.mean}}.`, { references: storeWith() });
        expect(html).toContain('has σ²θ = 111 deg² and a trans angle of 159.6°,');
        expect(html).not.toContain('{{');
        expect(html).toMatch(/<span class="llm-ref-missing" title="Not filled: the reference r5\.distortion_index\.mean names no value this question&#x27;s analyses returned">number missing<\/span>/);
        expect(html).not.toMatch(/<a [^>]*ref:/);
        // A ref: link the model wrote itself is a mark too, never a link.
        expect(render('[x](ref:javascript:alert(1))')).not.toMatch(/<a /);
    });

    it('the turn view: a streaming answer shows the numbers the steps hold; a finished one marks what stayed missing', () => {
        let turn = newAssistantTurn(0);
        turn = applyTurnEvent(turn, { type: 'round', round: 0 });
        turn = applyTurnEvent(turn, { type: 'tool_start', id: 'c1', name: 'polyhedral_distortion', label: 'Polyhedra' });
        turn = applyTurnEvent(turn, { type: 'tool_end', id: 'c1', name: 'polyhedral_distortion', status: 'ok', resultStatus: 'ok', data: R5, checks: [], ref: 'r5', ms: 1 });
        turn = applyTurnEvent(turn, { type: 'round', round: 1 });
        turn = applyTurnEvent(turn, { type: 'text', delta: `The average octahedron ${USER} and {{r5.data.average_structure.trans_angle_d` });
        const streaming = renderToStaticMarkup(<AssistantTurn turn={turn} pages={{}} cardLinks={null} />);
        expect(streaming).toContain('has σ²θ = 111 deg² and a trans angle of 159.6°,');
        expect(streaming).not.toContain('{{');
        const unfilled = 'The average octahedron has σ²θ = {{r5.data.average_structure.angle_variance}} deg².\n\nOutcome: computed.';
        const done = applyTurnEvent(turn, {
            type: 'done', text: unfilled, rawText: unfilled, outcome: { status: 'not verified', note: 'computed.', body: 'The average octahedron has σ²θ = {{r5.data.average_structure.angle_variance}} deg².' },
            stoppedBy: 'answer', cut: false, audit: [{ rule: 'reference', text: 'The answer quotes r5.data.average_structure.angle_variance by reference, which names no value this question\'s analyses returned: that number is missing.' }],
            references: { written: 1, filled: 0, unknown: ['r5.data.average_structure.angle_variance'] }, gates: []
        });
        const html = renderToStaticMarkup(<AssistantTurn turn={done} pages={{}} cardLinks={null} />);
        expect(html).toContain(`σ²θ = <span class="llm-ref-missing"`);
        expect(html).not.toContain('{{');
        expect(html).toContain('The answer quotes r5.data.average_structure.angle_variance by reference');
    });
});
