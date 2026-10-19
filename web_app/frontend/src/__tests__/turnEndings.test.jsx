// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// How a turn ended, read alike by the turn's two adapters. Every ending of
// the turn policy (llm/copilot/turnPolicy.js STOPPED_BY) is reached from the
// same scripted model through the dock (CopilotView → useCopilot → the turn)
// and through the validation harness (runLiveCase → the turn): the dock's turn
// ends with that `stoppedBy` and shows its row's note; the harness records it
// in `record.stopped_by`, and an answer's `record.final` reads as it did
// before the policy (its empty paths pinned field by field). A stall records
// 'error' and Stop records 'stopped', identically in both.
//
// The dock sets no limits (the loop's defaults: a stall after 3 min), so this
// file adds each scenario's limits to the turn's config in both adapters
// alike; nothing else of the turn is changed. No model runs and no server is
// contacted: the dock's port is providerModel on a non-local URL (never
// probed), streaming through a replaced streamCopilotChat.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';
import CopilotView from '../llm/components/CopilotView';
import { STOPPED_BY } from '../llm/copilot/turnPolicy';
import { buildCase, DEV_CASES } from '../../validation/cases';
import { runLiveCase } from '../../validation/harness';

const model = vi.hoisted(() => ({ script: null }));
vi.mock('../llm/provider/contextGuard', async (importOriginal) => ({
    ...(await importOriginal()),
    streamCopilotChat: (request) => model.script(request, { signal: request.signal })
}));
const extra = vi.hoisted(() => ({ limits: {} }));
vi.mock('../llm/copilot/turn', async (importOriginal) => {
    const turn = await importOriginal();
    return {
        ...turn,
        runTurn: (spec) => turn.runTurn({ ...spec, config: { ...spec.config, limits: { ...(spec.config?.limits || {}), ...extra.limits } } })
    };
});

const ANSWER = 'The Nb site holds its average position.\n\nOutcome: answered.';
const EMPTY_FINAL = 'The analysis budget for this question was used before the model wrote an answer. '
    + 'The analyses it ran are listed above.';
// A call to a tool neither adapter offers: answered as an error, so it adds nothing.
const BOGUS = { calls: [['no_such_tool', {}]] };

// The model of a scenario: each request takes the next reply (the last
// repeats) — { text, calls } — or 'wait': a reasoning fragment, then nothing
// until the request is aborted (Stop, or the stall clock). `stop` (optional)
// is called once that fragment is out.
const scripted = (replies, stop = null) => {
    let index = 0;
    return (request, { signal } = {}) => {
        const reply = replies[Math.min(index, replies.length - 1)];
        index += 1;
        return (async function* stream() {
            if (reply === 'wait') {
                yield { reasoning: 'Thinking about it' };
                await new Promise((resolve, reject) => {
                    signal.addEventListener('abort', () => reject(new DOMException('Aborted', 'AbortError')), { once: true });
                    stop?.();
                });
                return;
            }
            if (reply.text) yield { content: reply.text };
            for (const [position, [name, args]] of (reply.calls || []).entries()) {
                yield { toolCall: { index: position, id: `call_${index}_${position}`, name, arguments: JSON.stringify(args) } };
            }
            yield { finishReason: reply.calls ? 'tool_calls' : 'stop' };
        }());
    };
};

// One row per ending: the replies, the limits, and the note the dock shows.
const SCENARIOS = [
    { stoppedBy: 'answer', replies: [{ text: ANSWER }], limits: {} },
    { stoppedBy: 'max_rounds', replies: [BOGUS, { text: ANSWER }], limits: { maxRounds: 1 } },
    { stoppedBy: 'no_progress', replies: [BOGUS, BOGUS, { text: ANSWER }], limits: {} },
    { stoppedBy: 'empty', replies: [{ text: '' }], limits: {}, reason: 'the model wrote no answer, even after a reminder' },
    { stoppedBy: 'stopped', replies: ['wait'], limits: {}, stop: true },
    { stoppedBy: 'error', replies: ['wait'], limits: { stallMs: 20 } }
];

const RUN = Object.freeze({
    runName: 'endings-run', plotFiles: [], rValueFile: null, structure: null, symmetry: null,
    runSettings: null, pcaSites: null, liveData: false
});
const SETTINGS = Object.freeze({ baseUrl: 'https://api.example.test/v1', model: 'scripted', temperature: 0.2, apiKey: '' });
const MODEL = { name: 'scripted', label: 'scripted test model', model: 'scripted', baseUrl: 'none', temperature: 0 };
const built = () => buildCase(DEV_CASES.find((item) => item.id === 'dev-offcentre'));

const wait = (ms) => act(() => new Promise((resolve) => setTimeout(resolve, ms)));

describe('the endings, through the dock and the harness alike', () => {
    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        // No network: a non-local URL is never probed, and nothing else may fetch.
        vi.stubGlobal('fetch', vi.fn(async () => { throw new TypeError('no network in this test'); }));
    });

    afterEach(() => {
        vi.unstubAllGlobals();
        extra.limits = {};
    });

    // One question in a fresh dock: the turn it ended with, its footer and
    // its error line.
    const inTheDock = async (scenario) => {
        const ended = [];
        model.script = scripted(scenario.replies);
        const container = document.createElement('div');
        document.body.appendChild(container);
        const root = createRoot(container);
        act(() => root.render(
            <CopilotView
                context={null} run={RUN} settings={SETTINGS} connected tools={[]} toolsEnabled pages={{}}
                onOpenPage={() => {}} runKey="endings" onTurnEnd={(turn) => ended.push(turn)}
            />
        ));
        const textarea = container.querySelector('.llm-composer-input');
        act(() => {
            Object.getOwnPropertyDescriptor(HTMLTextAreaElement.prototype, 'value').set.call(textarea, 'Where is the Nb site?');
            textarea.dispatchEvent(new Event('input', { bubbles: true }));
        });
        act(() => container.querySelector('.llm-composer').dispatchEvent(new Event('submit', { bubbles: true, cancelable: true })));
        for (let i = 0; i < 100 && !ended.length; i += 1) {
            await wait(5);
            // Stop, once the reasoning is on screen.
            if (scenario.stop && container.querySelector('.llm-reasoning')) act(() => container.querySelector('.llm-send.is-stop')?.click());
        }
        const seen = {
            turn: ended[0],
            footer: container.querySelector('.llm-turn-footer')?.textContent ?? '',
            alert: container.querySelector('[role="alert"]')?.textContent ?? null
        };
        act(() => root.unmount());
        container.remove();
        expect(ended).toHaveLength(1);
        return seen;
    };

    // The same question in the harness: its record.
    const inTheHarness = async (scenario) => {
        const controller = new AbortController();
        return runLiveCase({
            built: built(), arm: 'C', model: MODEL, signal: controller.signal,
            port: scripted(scenario.replies, scenario.stop ? () => controller.abort() : null)
        });
    };

    it.each(SCENARIOS)('$stoppedBy: the dock\'s turn and the harness\'s record carry it, the dock shows its note', async (scenario) => {
        extra.limits = scenario.limits;
        const { turn, footer } = await inTheDock(scenario);
        const record = await inTheHarness(scenario);

        expect(turn.stoppedBy).toBe(scenario.stoppedBy);
        expect(record.stopped_by).toBe(scenario.stoppedBy);
        // The note under the turn is the row's (none for an answer), with an empty reply's reason.
        const { note } = STOPPED_BY[scenario.stoppedBy];
        if (note === null) {
            Object.values(STOPPED_BY).forEach((row) => { if (row.note) expect(footer).not.toContain(row.note.replace(/\.$/, '')); });
        } else {
            expect(footer).toContain(scenario.reason ? `${note.replace(/\.$/, '')}: ${scenario.reason}.` : note);
        }
        // A done ending is recorded in record.final too; a thrown one as the run's error.
        if (['stopped', 'error'].includes(scenario.stoppedBy)) {
            expect(record.final).toBeNull();
            expect(record.error.name).toBe(scenario.stoppedBy === 'stopped' ? 'AbortError' : 'ModelTimeoutError');
            expect(turn.outcome).toBeNull();
        } else {
            expect(record.error).toBeNull();
            expect(record.final.stopped_by).toBe(scenario.stoppedBy);
        }
    });

    it('a stall records error and Stop records stopped, identically in both, with nothing else recorded as an answer', async () => {
        const read = async (scenario) => {
            extra.limits = scenario.limits;
            const { turn, footer, alert } = await inTheDock(scenario);
            const record = await inTheHarness(scenario);
            return {
                dock: { stoppedBy: turn.stoppedBy, footer, alert: Boolean(alert) },
                harness: { stoppedBy: record.stopped_by, error: record.error.name, final: record.final, requests: record.requests.length }
            };
        };
        const stall = await read(SCENARIOS.find((row) => row.stoppedBy === 'error'));
        const stop = await read(SCENARIOS.find((row) => row.stoppedBy === 'stopped'));
        // The stall explained in the dock's error line; Stop shows no error.
        expect(stall).toEqual({
            dock: { stoppedBy: 'error', footer: STOPPED_BY.error.note, alert: true },
            harness: { stoppedBy: 'error', error: 'ModelTimeoutError', final: null, requests: 1 }
        });
        expect(stop).toEqual({
            dock: { stoppedBy: 'stopped', footer: STOPPED_BY.stopped.note, alert: false },
            harness: { stoppedBy: 'stopped', error: 'AbortError', final: null, requests: 1 }
        });
    });

    it('the empty paths\' record.final reads as before the policy: rawText, references and reason as the harness kept them', async () => {
        const empty = await inTheHarness({ replies: [{ text: '' }] });
        expect(empty.final).toEqual({
            text: '', raw_text: '', outcome: null, stopped_by: 'empty', cut: false, audit: [], references: null, reasoning: '',
            reason: 'the model wrote no answer, even after a reminder'
        });
        extra.limits = { maxRounds: 1 };
        const emptyFinal = await inTheHarness({ replies: [BOGUS, { text: '' }] });
        expect(emptyFinal.final).toEqual({
            text: EMPTY_FINAL, raw_text: EMPTY_FINAL,
            outcome: { status: 'not achieved', note: 'no answer was written within the analysis budget', body: EMPTY_FINAL },
            stopped_by: 'max_rounds', cut: false, audit: [], references: null, reasoning: '', reason: null
        });
        expect(emptyFinal.stopped_by).toBe('max_rounds');
        // An answer: raw_text as written, reason null, references the verification's.
        extra.limits = {};
        const answered = await inTheHarness({ replies: [{ text: ANSWER }] });
        expect(answered.final).toMatchObject({ text: ANSWER, raw_text: ANSWER, stopped_by: 'answer', reason: null, references: { written: 0, filled: 0, unknown: [] } });
    });
});
