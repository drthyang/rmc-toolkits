// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// A long answer in the dock: its lead (the conclusion) shows, the rest folds
// behind "Show details" (components/answerParts.js); a short one shows whole,
// and nothing is dropped. Through the turn view (CopilotView's AssistantTurn).

import { act } from 'react';
import { createRoot } from 'react-dom/client';
import { renderToStaticMarkup } from 'react-dom/server';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { AssistantTurn } from '../components/CopilotView';
import { answerParts, FOLD_MIN_CHARS, markdownBlocks } from '../components/answerParts';
import { applyTurnEvent, newAssistantTurn } from '../copilot/turnState';

// The kind of answer the user reported (2026-11-01), filled: a conclusion,
// then every index of polyhedral_distortion restated, a table and caveats.
const LONG = [
    'The TaSe₆ octahedra (Ta–Se 2.21–3.23 Å) are distorted in the average structure itself: σ²θ is 111 deg² there '
        + 'against 127 deg² in the configuration, so most of the distortion is static.',
    '## Configuration (one configuration)',
    '- Baur DI: mean 0.0469, sd 0.017, p10 0.0261, median 0.0459, p90 0.0683\n'
        + '- ⟨λ⟩: mean 1.0405, sd 0.013, p10 1.0247, median 1.0403, p90 1.0563\n'
        + '- σ²θ: mean 127 deg², sd 37, p10 79.4, median 128, p90 174\n'
        + '- trans angle: mean 158.8°, sd 3.7°',
    '| Index | Configuration | Average |\n| --- | --- | --- |\n| DI | 0.0469 | 0.0415 |\n| ⟨λ⟩ | 1.0405 | 1.0341 |\n| σ²θ | 127 | 111 |',
    '```\nbond ladder: 2.523 2.525 2.527 2.74 2.744 2.747\n\nsix bonds\n```',
    'These numbers come from one configuration; the coordination counts include 68 Ta atoms with 5 Se, left out.'
].join('\n\n');

describe('answerParts', () => {
    it('splits at blank lines outside fenced code', () => {
        const blocks = markdownBlocks(LONG);
        expect(blocks).toHaveLength(6);
        expect(blocks[4]).toBe('```\nbond ladder: 2.523 2.525 2.527 2.74 2.744 2.747\n\nsix bonds\n```');
    });

    it('a long answer: the conclusion leads, the rest folds, nothing is dropped', () => {
        expect(LONG.length).toBeGreaterThan(FOLD_MIN_CHARS);
        const { lead, rest } = answerParts(LONG);
        expect(lead).toBe(markdownBlocks(LONG)[0]);
        expect(rest.startsWith('## Configuration (one configuration)')).toBe(true);
        expect(`${lead}\n\n${rest}`).toBe(LONG);
    });

    it('a heading or a short first line leads with the next block', () => {
        const titled = `## Short answer\n\n${LONG}`;
        expect(answerParts(titled).lead).toBe(`## Short answer\n\n${markdownBlocks(LONG)[0]}`);
        expect(answerParts(`**Yes.**\n\n${LONG}`).lead).toBe(`**Yes.**\n\n${markdownBlocks(LONG)[0]}`);
    });

    // The shape the dock's prompt asks for (prompts/copilot.js ANSWER_RULES)
    // shows whole: the few numbers that carry the conclusion stay in view.
    const SHAPED = 'The TaSe₆ octahedra are distorted in the average structure itself, so most of the distortion is static.\n\n'
        + '- σ²θ: 111 deg² in the average structure, 127 deg² in the configuration (2.21–3.23 Å)\n'
        + '- trans angle: 159.6° against 158.8°\n\n'
        + 'From one configuration; the 68 Ta atoms with 5 Se are left out.';

    it('an answer in the prompt\'s shape shows whole; a short one whose rest is a table, headings or many numbers folds', () => {
        expect(answerParts(SHAPED)).toEqual({ lead: SHAPED, rest: '' });
        const table = `${markdownBlocks(LONG)[0]}\n\n${markdownBlocks(LONG)[3]}`;
        expect(table.length).toBeLessThan(FOLD_MIN_CHARS);
        expect(answerParts(table).rest).toBe(markdownBlocks(LONG)[3]);
        const numbers = `${markdownBlocks(LONG)[0]}\n\n${markdownBlocks(LONG)[2]}`;
        expect(numbers.length).toBeLessThan(FOLD_MIN_CHARS);
        expect(answerParts(numbers).rest).toBe(markdownBlocks(LONG)[2]);
    });

    it('a short answer, one block, or a short rest shows whole', () => {
        const short = 'The octahedra are distorted: σ²θ 111 deg² in the average structure.\n\nFrom one configuration.';
        expect(answerParts(short)).toEqual({ lead: short, rest: '' });
        const oneBlock = 'x '.repeat(500).trim();
        expect(answerParts(oneBlock)).toEqual({ lead: oneBlock, rest: '' });
        const shortRest = `${'A long conclusion sentence. '.repeat(30).trim()}\n\nOne caveat.`;
        expect(answerParts(shortRest).rest).toBe('');
        expect(answerParts(null)).toEqual({ lead: '', rest: '' });
    });
});

describe('the turn view folds a long finished answer', () => {
    let container = null;
    let root = null;
    beforeEach(() => { globalThis.IS_REACT_ACT_ENVIRONMENT = true; });
    afterEach(() => {
        act(() => root?.unmount());
        container?.remove();
        container = null;
        root = null;
    });

    const turnWith = (text, { streaming = false } = {}) => {
        let turn = newAssistantTurn(0);
        turn = applyTurnEvent(turn, { type: 'round', round: 0 });
        turn = applyTurnEvent(turn, { type: 'text', delta: text });
        if (streaming) return turn;
        return applyTurnEvent(turn, {
            type: 'done', text, rawText: text, outcome: { status: 'achieved', note: 'the indices of both.', body: text },
            stoppedBy: 'answer', cut: false, audit: [], references: null, gates: []
        });
    };

    it('shows the lead and a Show details button; the rest opens on click', () => {
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
        act(() => root.render(<AssistantTurn turn={turnWith(LONG)} pages={{}} cardLinks={null} />));
        const button = container.querySelector('.llm-answer-toggle');
        expect(button.textContent).toBe('Show details');
        expect(button.getAttribute('aria-expanded')).toBe('false');
        const details = container.querySelector(`#${CSS.escape(button.getAttribute('aria-controls'))}`);
        expect(details.hidden).toBe(true);
        expect(container.querySelector('.llm-bubble').textContent).toContain('most of the distortion is static');
        expect(container.querySelector('.llm-bubble').textContent).not.toContain('p90');
        expect(container.querySelector('table')).toBeNull();
        act(() => button.click());
        expect(button.getAttribute('aria-expanded')).toBe('true');
        expect(button.textContent).toBe('Hide details');
        expect(details.hidden).toBe(false);
        expect(details.textContent).toContain('p90 0.0683');
        expect(details.querySelector('table')).not.toBeNull();
        // The Outcome chip stays under the answer.
        expect(container.querySelector('.llm-outcome').textContent).toContain('the indices of both.');
    });

    it('a streaming answer and a short one show whole', () => {
        const streaming = renderToStaticMarkup(<AssistantTurn turn={turnWith(LONG, { streaming: true })} pages={{}} cardLinks={null} />);
        expect(streaming).not.toContain('llm-answer-toggle');
        expect(streaming).toContain('p90 0.0683');
        const short = renderToStaticMarkup(<AssistantTurn turn={turnWith('Distorted: σ²θ 111 deg².')} pages={{}} cardLinks={null} />);
        expect(short).not.toContain('llm-answer-toggle');
    });
});
