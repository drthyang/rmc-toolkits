// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// Which models take a temperature (provider/sampling.js). Current Claude
// models refuse one with HTTP 400 ("`temperature` is deprecated for this
// model"), so the Workbench's own requests send them none: the watchdog's
// request here (the dock's: modelPort.test.js and copilotChat.test.jsx), and
// the settings drawer greys the field out. Every other model keeps its
// temperature. Nothing leaves the test: fetch is a stub.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';
import { renderToStaticMarkup } from 'react-dom/server';
import ConnectionSettings from '../components/ConnectionSettings';
import { takesNoTemperature, temperatureFor } from '../provider/sampling';
import { useWatchdog } from '../watchdog/useWatchdog';

const ANTHROPIC = 'https://api.anthropic.com/v1';
const OLLAMA = 'http://localhost:11434/v1';
const LMSTUDIO = 'http://localhost:1234/v1';

describe('takesNoTemperature and temperatureFor', () => {
    it('true for a claude- model id on any server, and for any model on Anthropic\'s API', () => {
        expect(takesNoTemperature({ baseUrl: ANTHROPIC, model: 'claude-opus-5-5' })).toBe(true);
        expect(takesNoTemperature({ baseUrl: ANTHROPIC, model: 'claude-sonnet-5-5' })).toBe(true);
        expect(takesNoTemperature({ baseUrl: ANTHROPIC, model: '' })).toBe(true);
        expect(takesNoTemperature({ baseUrl: 'https://gateway.example/v1', model: 'claude-haiku-4-5' })).toBe(true);
        expect(takesNoTemperature({ baseUrl: LMSTUDIO, model: ' Claude-Opus-5-5 ' })).toBe(true);
    });

    it('false for local and other cloud models', () => {
        expect(takesNoTemperature({ baseUrl: OLLAMA, model: 'qwen3:32b' })).toBe(false);
        expect(takesNoTemperature({ baseUrl: LMSTUDIO, model: 'muse-glimmer' })).toBe(false);
        expect(takesNoTemperature({ baseUrl: 'https://api.openai.com/v1', model: 'gpt-5' })).toBe(false);
        // Only a model id that starts with claude-; the host must be Anthropic's own.
        expect(takesNoTemperature({ baseUrl: OLLAMA, model: 'my-claude-distill' })).toBe(false);
        expect(takesNoTemperature({ baseUrl: 'https://api.anthropic.com.example/v1', model: 'm' })).toBe(false);
        expect(takesNoTemperature({ baseUrl: 'not a url', model: null })).toBe(false);
        expect(takesNoTemperature()).toBe(false);
    });

    it('temperatureFor: null for those models, else the settings\' own or the one given', () => {
        expect(temperatureFor({ baseUrl: ANTHROPIC, model: 'claude-opus-5-5', temperature: 0.2 })).toBeNull();
        expect(temperatureFor({ baseUrl: ANTHROPIC, model: 'claude-opus-5-5' }, 0)).toBeNull();
        expect(temperatureFor({ baseUrl: OLLAMA, model: 'qwen3:32b', temperature: 0.2 })).toBe(0.2);
        expect(temperatureFor({ baseUrl: OLLAMA, model: 'qwen3:32b', temperature: 0.2 }, 0)).toBe(0);
    });
});

describe('the watchdog\'s request (useWatchdog → completeChat)', () => {
    let container;
    let root;
    let bodies;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        bodies = [];
        vi.stubGlobal('fetch', vi.fn(async (url, options = {}) => {
            bodies.push({ url: String(url), body: JSON.parse(options.body) });
            return { ok: true, status: 200, json: async () => ({ choices: [{ message: { content: 'STATUS: improving\nNOTE: Still falling.' } }] }) };
        }));
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.unstubAllGlobals();
    });

    const rValueFile = { path: 'run-01.log', plotData: { chiColumn: null, series: [{ label: 'chi', x: [0, 1, 2, 3], y: [5, 4, 3, 2] }] } };
    const Probe = ({ settings }) => {
        useWatchdog({ rValueFile, settings });
        return null;
    };
    const watch = async (baseUrl, model) => {
        await act(async () => {
            root.render(<Probe settings={{ baseUrl, model, apiKey: '', watchdogEnabled: true, watchdogIntervalMin: 5 }} />);
        });
        await act(async () => { await new Promise((resolve) => setTimeout(resolve, 0)); });
        expect(bodies).toHaveLength(1);
        return bodies[0];
    };

    it('Claude: no temperature key', async () => {
        const { url, body } = await watch(ANTHROPIC, 'claude-opus-5-5');
        expect(url).toBe(`${ANTHROPIC}/chat/completions`);
        expect(body.model).toBe('claude-opus-5-5');
        expect('temperature' in body).toBe(false);
    });

    it('a local model (Ollama, LM Studio): temperature 0, as before', async () => {
        for (const [baseUrl, model] of [[OLLAMA, 'qwen3:32b'], [LMSTUDIO, 'muse-glimmer']]) {
            bodies.length = 0;
            act(() => root.unmount());
            root = createRoot(container);
            const { body } = await watch(baseUrl, model);
            expect(body.temperature).toBe(0);
        }
    });
});

describe('the settings drawer\'s temperature field', () => {
    const drawer = (settings) => renderToStaticMarkup(
        <ConnectionSettings
            settings={{ apiKey: '', temperature: 0.2, copilotTools: true, watchdogEnabled: false, watchdogIntervalMin: 5, ...settings }}
            connection={{ status: 'idle', models: [] }}
            onSave={() => {}}
            onTest={async () => null}
        />
    );
    const field = (html) => html.match(/<input[^>]*id="llm-temperature"[^>]*>/)[0];

    it('Claude: off, with one plain line why', () => {
        const html = drawer({ baseUrl: ANTHROPIC, model: 'claude-opus-5-5' });
        expect(field(html)).toContain('disabled=""');
        expect(field(html)).toContain('value=""');
        expect(html).toContain('Claude uses its own sampling');
    });

    it('a local model: the setting, editable, no note', () => {
        const html = drawer({ baseUrl: OLLAMA, model: 'qwen3:32b' });
        expect(field(html)).not.toContain('disabled');
        expect(field(html)).toContain('value="0.2"');
        expect(html).not.toContain('Claude uses its own sampling');
    });
});
