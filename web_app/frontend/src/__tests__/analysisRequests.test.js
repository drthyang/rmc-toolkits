// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The analysis request contract (fixtures/analysis_requests.json), checked at
// the seam, requestAnalysis (useSiteCloud.js), on every JS adapter:
// - every case through the worker adapter (requestAnalysis with a loopback
//   Worker that runs the PCA worker's handler) and the in-process adapter
//   (pcaKdeWorker's inProcessEngine), with the worker's expectation;
// - every case through the Flask adapter with axios mocked: each kind reaches
//   its route with { ...params, dir }, and an unknown kind sends nothing (the
//   Flask outcomes are tests/test_analysis_requests.py's);
// - the cancel table on the worker, Flask and in-process adapters and on
//   useCopilotRun's engine (its worker terminated on abort, then a fresh one);
// - the AI Copilot tool schemas (copilotTools.js, unchanged) against the
//   limits, and the fixture's own coverage: every bound has an at-bound and a
//   past-bound case, and the number of kept differences is fixed.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';
import fixture from './fixtures/analysis_requests.json';

const axiosGet = vi.hoisted(() => vi.fn());
vi.mock('axios', () => ({ default: { get: axiosGet } }));

const { requestAnalysis } = await import('../useSiteCloud');
const { default: API_BASE_URL } = await import('../api');
const { ANALYSIS_ROUTES } = await import('../workers/requestGuards');
const { handlePcaMessage, inProcessEngine, MAX_CLUSTER_WINDOW, MAX_DISTORTION_WINDOW } = await import('../workers/pcaKdeWorker');
const { buildCopilotTools } = await import('../copilotTools');
const { default: useCopilotRun } = await import('../useCopilotRun');

const ANALYSIS_TOOLS = ['list_sites', 'analyze_displacement_directions', 'compute_bond_angles', 'analyze_clusters', 'analyze_modulations'];
// The kept differences between the adapters, by kind: accepted on one and
// refused on the other (or accepted with other values), refused on both under
// another name, and read by one adapter only. A new one shows here in review.
const KEPT = { outcome: 8, naming: 10, inert: 1 };

// --- reading the fixture ---------------------------------------------------------

/** A case's expected outcome on one adapter: { ok, fields } | { error } | { aborted }. */
const outcomeOf = (testCase, adapter) => {
    const split = testCase.expect && typeof testCase.expect === 'object' && 'flask' in testCase.expect;
    const expected = split ? testCase.expect[adapter] : testCase.expect;
    if (expected === 'ok') return { ok: true, fields: split ? {} : testCase.fields ?? {} };
    if (expected === 'aborted') return { aborted: true };
    if (expected.fields) return { ok: true, fields: expected.fields };
    return { error: expected.error };
};

/** A payload value by dotted path (`.length` of a list); a missing value is null. */
const field = (payload, path) => path.split('.').reduce((value, part) => {
    if (Array.isArray(value)) return part === 'length' ? value.length : value[Number(part)];
    return value !== null && typeof value === 'object' && part in value ? value[part] : null;
}, payload);

/** Settle a request into { result } | { error }. */
const settle = (promise) => promise.then((result) => ({ result }), (error) => ({ error }));

const expectOutcome = (got, expected, label) => {
    if (expected.aborted) {
        expect(got.error?.name, label).toBe('AbortError');
        return;
    }
    if (expected.error) {
        expect(got.error, `${label}: expected a refusal naming ${expected.error}`).toBeInstanceOf(Error);
        expect(got.error.name, label).not.toBe('AbortError');
        expect(got.error.message, label).toContain(expected.error);
        return;
    }
    expect(got.error, `${label}: ${got.error?.message}`).toBeUndefined();
    for (const [path, value] of Object.entries(expected.fields)) {
        expect(field(got.result, path), `${label}: ${path}`).toEqual(value);
    }
};

// --- adapters ----------------------------------------------------------------------

/**
 * A Worker that answers like pcaKdeWorker's onmessage, in this thread:
 * handlePcaMessage over the posted message and its text.
 */
class LoopbackWorker {
    constructor({ onPost = () => {} } = {}) {
        this.listeners = new Set();
        this.posted = [];
        this.onPost = onPost;
        this.terminated = 0;
    }

    addEventListener(type, listener) { if (type === 'message') this.listeners.add(listener); }

    removeEventListener(type, listener) { if (type === 'message') this.listeners.delete(listener); }

    terminate() { this.terminated += 1; }

    postMessage(data) {
        this.posted.push(data);
        this.onPost(data);
        handlePcaMessage(data, async () => data.text).then(
            (result) => ({ id: data.id, result }),
            (error) => ({ id: data.id, error: error.message })
        ).then((reply) => { for (const listener of [...this.listeners]) listener({ data: reply }); });
    }
}

/** One case through the worker adapter (requestAnalysis' worker branch). */
const viaWorker = async (testCase) => {
    const controller = new AbortController();
    if (testCase.probe === 'aborted-before') controller.abort();
    const worker = new LoopbackWorker({ onPost: () => { if (testCase.probe === 'aborted-during') controller.abort(); } });
    const got = await settle(requestAnalysis(testCase.kind ?? undefined, testCase.params, {
        rmc6fText: fixture.rmc6f, worker, signal: controller.signal
    }));
    return { got, worker };
};

/** One case through the in-process adapter (a fresh engine, so the text is read by this case). */
const viaInProcess = async (testCase) => {
    const controller = new AbortController();
    if (testCase.probe === 'aborted-before') controller.abort();
    let reads = 0;
    const engine = inProcessEngine(() => {
        reads += 1;
        if (testCase.probe === 'aborted-during') controller.abort();
        return fixture.rmc6f;
    });
    const got = await settle(engine(testCase.kind ?? undefined, testCase.params, { signal: controller.signal }));
    return { got, reads };
};

beforeEach(() => {
    axiosGet.mockReset();
});

describe('the fixture', () => {
    it('names the kinds requestAnalysis routes by (requestGuards ANALYSIS_ROUTES)', () => {
        expect(fixture.routes).toEqual({ ...ANALYSIS_ROUTES });
    });

    it('has unique case names, and every limit names a kind and a parameter some case sends', () => {
        expect(new Set(fixture.cases.map((c) => c.name)).size).toBe(fixture.cases.length);
        for (const limit of fixture.limits) {
            expect(Object.keys(fixture.routes)).toContain(limit.kind);
            expect(fixture.cases.some((c) => c.kind === limit.kind && limit.param in c.params), `${limit.kind}.${limit.param}`).toBe(true);
        }
    });

    it('every bound has an at-bound case and a past-bound case, each on its side of it', () => {
        const ok = (c) => outcomeOf(c, 'flask').ok && outcomeOf(c, 'worker').ok;
        // Refused by Flask; the worker refuses too unless the difference is kept.
        const refused = (c) => Boolean(outcomeOf(c, 'flask').error) && (Boolean(outcomeOf(c, 'worker').error) || Boolean(c.kept));
        for (const limit of fixture.limits) {
            const label = `${limit.kind}.${limit.param}`;
            const sent = fixture.cases.filter((c) => c.kind === limit.kind && !c.probe && typeof c.params[limit.param] === 'number');
            const has = (test) => sent.some((c) => test(c.params[limit.param], c));
            if (limit.ge !== undefined) {
                expect(has((v, c) => v === limit.ge && ok(c)), `${label} at ${limit.ge}`).toBe(true);
                expect(has((v, c) => v < limit.ge && refused(c)), `${label} below ${limit.ge}`).toBe(true);
            }
            if (limit.le !== undefined) {
                expect(has((v, c) => v === limit.le && ok(c)), `${label} at ${limit.le}`).toBe(true);
                expect(has((v, c) => v > limit.le && refused(c)), `${label} above ${limit.le}`).toBe(true);
            }
            if (limit.gt !== undefined) {
                expect(has((v, c) => v === limit.gt && refused(c)), `${label} at ${limit.gt}`).toBe(true);
                expect(has((v, c) => v > limit.gt && ok(c)), `${label} above ${limit.gt}`).toBe(true);
            }
            if (limit.lt !== undefined) {
                expect(has((v, c) => v === limit.lt && refused(c)), `${label} at ${limit.lt}`).toBe(true);
                expect(has((v, c) => v < limit.lt && ok(c)), `${label} below ${limit.lt}`).toBe(true);
            }
            if (limit.clamp) {
                const [lo, hi] = limit.clamp;
                const shows = (c, value) => ok(c) && c.fields?.[limit.param] === value;
                expect(has((v, c) => v === lo && shows(c, lo)), `${label} at ${lo}`).toBe(true);
                expect(has((v, c) => v < lo && shows(c, lo)), `${label} below ${lo}, clamped`).toBe(true);
                expect(has((v, c) => v === hi && shows(c, hi)), `${label} at ${hi}`).toBe(true);
                expect(has((v, c) => v > hi && shows(c, hi)), `${label} above ${hi}, clamped`).toBe(true);
            }
            if (limit.integer) {
                expect(has((v, c) => !Number.isInteger(v) && refused(c)), `${label} fractional`).toBe(true);
            }
        }
    });

    it('every difference between the adapters is kept with a reason, and their number is fixed', () => {
        const counts = { outcome: 0, naming: 0, inert: 0 };
        for (const c of fixture.cases) {
            const split = typeof c.expect === 'object' && 'flask' in c.expect;
            if (split) expect(c.kept, `${c.name}: a split expectation needs its kept reason`).toBeTruthy();
            if (!c.kept) continue;
            expect(c.kept.length, c.name).toBeGreaterThan(40);
            if (!split) counts.inert += 1;
            else if (outcomeOf(c, 'flask').error && outcomeOf(c, 'worker').error) counts.naming += 1;
            else counts.outcome += 1;
        }
        expect(counts).toEqual(KEPT);
    });

    it('the worker\'s copy of a cap equals the contract\'s', () => {
        const rMax = fixture.limits.find((limit) => limit.kind === 'clusters' && limit.param === 'rMax');
        expect(MAX_CLUSTER_WINDOW).toBe(rMax.le);
        const distortion = fixture.limits.find((limit) => limit.kind === 'distortion' && limit.param === 'rMax');
        expect(MAX_DISTORTION_WINDOW).toBe(distortion.le);
    });
});

describe('every case on the worker and in-process adapters', () => {
    for (const testCase of fixture.cases) {
        it(testCase.name, async () => {
            const expected = outcomeOf(testCase, 'worker');
            const worker = await viaWorker(testCase);
            expectOutcome(worker.got, expected, `worker: ${testCase.name}`);
            const inProcess = await viaInProcess(testCase);
            expectOutcome(inProcess.got, expected, `in-process: ${testCase.name}`);
            if (!(testCase.kind in fixture.routes) || testCase.probe === 'aborted-before') {
                // Refused before any work: nothing posted, no text read.
                expect(worker.worker.posted).toEqual([]);
                expect(inProcess.reads).toBe(0);
            } else {
                // The posted message: the params, and kind, text and id that
                // no param can override.
                expect(worker.worker.posted).toHaveLength(1);
                expect(worker.worker.posted[0]).toEqual({
                    ...testCase.params, id: expect.any(Number), kind: testCase.kind, text: fixture.rmc6f
                });
                expect(inProcess.reads).toBe(1);
            }
        });
    }
});

describe('every case on the Flask adapter (axios mocked)', () => {
    for (const testCase of fixture.cases) {
        it(testCase.name, async () => {
            const controller = new AbortController();
            if (testCase.probe === 'aborted-before') controller.abort();
            axiosGet.mockImplementation((url, config) => {
                if (testCase.probe === 'aborted-during') {
                    // axios rejects a cancelled request with its CanceledError.
                    return new Promise((resolve, reject) => {
                        config.signal.addEventListener('abort', () => reject(Object.assign(new Error('canceled'), { code: 'ERR_CANCELED' })));
                        controller.abort();
                    });
                }
                return Promise.resolve({ data: { route: url } });
            });
            const got = await settle(requestAnalysis(testCase.kind ?? undefined, testCase.params, {
                directory: 'runs/contract', signal: controller.signal
            }));
            if (!(testCase.kind in fixture.routes)) {
                expectOutcome(got, { error: 'kind' }, testCase.name);
                expect(axiosGet).not.toHaveBeenCalled();
                return;
            }
            if (testCase.probe === 'aborted-before') {
                expectOutcome(got, { aborted: true }, testCase.name);
                expect(axiosGet).not.toHaveBeenCalled();
                return;
            }
            expect(axiosGet).toHaveBeenCalledTimes(1);
            expect(axiosGet).toHaveBeenCalledWith(`${API_BASE_URL}${fixture.routes[testCase.kind]}`, {
                params: { ...testCase.params, dir: 'runs/contract' }, signal: controller.signal
            });
            if (testCase.probe === 'aborted-during') expectOutcome(got, { aborted: true }, testCase.name);
            else expect(got.result).toEqual({ route: `${API_BASE_URL}${fixture.routes[testCase.kind]}` });
        });
    }

    it('dir is the run directory whatever params say, and a server refusal is its own message', async () => {
        axiosGet.mockResolvedValue({ data: {} });
        await requestAnalysis('sites', { dir: '/etc', probability: 0.9 }, { directory: '' });
        expect(axiosGet).toHaveBeenLastCalledWith(`${API_BASE_URL}/api/pca/sites`, { params: { dir: '.', probability: 0.9 } });
        axiosGet.mockRejectedValue(Object.assign(new Error('Request failed with status code 400'), {
            response: { data: { error: 'probability must be > 0, got 0' } }
        }));
        await expect(requestAnalysis('sites', { probability: 0 }, { directory: 'run' })).rejects.toThrow('probability must be > 0, got 0');
    });
});

// --- the cancel table ------------------------------------------------------------------

/** A Worker that only records: posts, listeners, termination; `reply` answers the last post. */
class RecordingWorker {
    static instances = [];

    constructor() {
        this.listeners = new Set();
        this.posted = [];
        this.terminated = 0;
        RecordingWorker.instances.push(this);
    }

    addEventListener(type, listener) { if (type === 'message') this.listeners.add(listener); }

    removeEventListener(type, listener) { if (type === 'message') this.listeners.delete(listener); }

    postMessage(data) { this.posted.push(data); }

    terminate() { this.terminated += 1; }

    reply(result) {
        const { id } = this.posted.at(-1);
        for (const listener of [...this.listeners]) listener({ data: { id, result } });
    }
}

describe('cancel: the worker adapter (requestAnalysis given a worker)', () => {
    it('an already aborted signal rejects with AbortError; nothing is posted', async () => {
        const worker = new RecordingWorker();
        const controller = new AbortController();
        controller.abort();
        await expect(requestAnalysis('sites', {}, { rmc6fText: 'x', worker, signal: controller.signal })).rejects.toMatchObject({ name: 'AbortError' });
        expect(worker.posted).toEqual([]);
        expect(worker.listeners.size).toBe(0);
    });

    it('an abort during the request rejects at once; the worker is not terminated (its owner decides)', async () => {
        const worker = new RecordingWorker();
        const controller = new AbortController();
        const pending = requestAnalysis('sites', {}, { rmc6fText: 'x', worker, signal: controller.signal });
        expect(worker.posted).toHaveLength(1);
        controller.abort();
        await expect(pending).rejects.toMatchObject({ name: 'AbortError' });
        expect(worker.listeners.size).toBe(0);
        expect(worker.terminated).toBe(0);
        // A late answer reaches no one.
        expect(() => worker.reply({ late: true })).not.toThrow();
    });
});

describe('cancel: the Flask adapter', () => {
    it('an already aborted signal rejects with AbortError; nothing is fetched', async () => {
        const controller = new AbortController();
        controller.abort();
        await expect(requestAnalysis('kde', {}, { directory: 'run', signal: controller.signal })).rejects.toMatchObject({ name: 'AbortError' });
        expect(axiosGet).not.toHaveBeenCalled();
    });

    it('an abort during the request cancels the HTTP request (axios holds the signal) and rejects with AbortError', async () => {
        const controller = new AbortController();
        let seen = null;
        axiosGet.mockImplementation((url, config) => new Promise((resolve, reject) => {
            seen = config.signal;
            config.signal.addEventListener('abort', () => reject(Object.assign(new Error('canceled'), { code: 'ERR_CANCELED' })));
        }));
        const pending = requestAnalysis('kde', {}, { directory: 'run', signal: controller.signal });
        controller.abort();
        await expect(pending).rejects.toMatchObject({ name: 'AbortError' });
        expect(seen).toBe(controller.signal);
    });
});

describe('cancel: the in-process adapter (inProcessEngine)', () => {
    it('an already aborted signal rejects with AbortError before the text is read', async () => {
        const source = vi.fn(() => fixture.rmc6f);
        const controller = new AbortController();
        controller.abort();
        await expect(inProcessEngine(source)('sites', {}, { signal: controller.signal })).rejects.toMatchObject({ name: 'AbortError' });
        expect(source).not.toHaveBeenCalled();
    });

    it('an abort during the request runs it to its end, then rejects with AbortError (it used to resolve)', async () => {
        const controller = new AbortController();
        const engine = inProcessEngine(async () => {
            controller.abort();
            return fixture.rmc6f;
        });
        await expect(engine('sites', {}, { signal: controller.signal })).rejects.toMatchObject({ name: 'AbortError' });
        // The text was read once and the same engine answers the next request.
        const sites = await engine('sites', {});
        expect(sites.totalAtoms).toBe(128);
    });

    it('a refusal during an aborted request is the abort', async () => {
        const controller = new AbortController();
        const engine = inProcessEngine(async () => {
            controller.abort();
            return fixture.rmc6f;
        });
        await expect(engine('sites', { probability: 2 }, { signal: controller.signal })).rejects.toMatchObject({ name: 'AbortError' });
    });
});

describe('cancel: the Copilot\'s engine (useCopilotRun) terminates its worker on abort', () => {
    let container;
    let root;
    let hook;
    const Probe = (props) => {
        hook = useCopilotRun(props);
        return null;
    };

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        RecordingWorker.instances = [];
        vi.stubGlobal('Worker', RecordingWorker);
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.unstubAllGlobals();
    });

    const localRun = { name: 'contract', runId: 'contract', structureFile: { path: 'contract/run.rmc6f', sourceFile: { text: async () => fixture.rmc6f } } };
    const posted = async (count) => {
        for (let i = 0; i < 50 && RecordingWorker.instances.flatMap((w) => w.posted).length < count; i += 1) {
            await new Promise((resolve) => setTimeout(resolve, 0));
        }
    };

    it('aborted before: AbortError, nothing posted; aborted during: terminated, then a fresh worker answers', async () => {
        await act(async () => { root.render(React.createElement(Probe, { directory: null, localRun })); });
        const engine = hook.engine;

        const before = new AbortController();
        before.abort();
        await expect(engine('sites', {}, { signal: before.signal })).rejects.toMatchObject({ name: 'AbortError' });
        expect(RecordingWorker.instances.flatMap((w) => w.posted)).toEqual([]);

        const during = new AbortController();
        const pending = engine('sites', {}, { signal: during.signal });
        await posted(1);
        const first = RecordingWorker.instances.at(-1);
        expect(first.posted).toHaveLength(1);
        during.abort();
        await expect(pending).rejects.toMatchObject({ name: 'AbortError' });
        expect(first.terminated).toBe(1);

        const next = engine('sites', {});
        await posted(2);
        const fresh = RecordingWorker.instances.at(-1);
        expect(fresh).not.toBe(first);
        expect(fresh.posted[0]).toMatchObject({ kind: 'sites', text: fixture.rmc6f });
        fresh.reply({ totalAtoms: 128 });
        await expect(next).resolves.toEqual({ totalAtoms: 128 });
        expect(fresh.terminated).toBe(0);
    });
});

// --- the AI Copilot's tool schemas -------------------------------------------------------

describe('the AI Copilot tool schemas lie inside the contract (copilotTools.js unchanged)', () => {
    const tools = buildCopilotTools({ request: async () => ({}), structure: null, elements: ['Nb', 'Se'] });
    const property = (name) => {
        const [tool, key] = name.split('.');
        return tools.find((entry) => entry.name === tool)?.parameters.properties[key];
    };

    it('every numeric property of the analysis tools is named by a limit', () => {
        expect(tools.map((tool) => tool.name)).toEqual(ANALYSIS_TOOLS);
        const named = new Set(fixture.limits.map((limit) => limit.tool).filter(Boolean));
        for (const tool of tools) {
            for (const [key, schema] of Object.entries(tool.parameters.properties)) {
                if (schema.type === 'number' || schema.type === 'integer') {
                    expect(named.has(`${tool.name}.${key}`), `${tool.name}.${key}`).toBe(true);
                }
            }
        }
    });

    it('each such property\'s [minimum, maximum] lies inside its range, and integer ↔ type integer', () => {
        for (const limit of fixture.limits.filter((entry) => entry.tool)) {
            const schema = property(limit.tool);
            expect(schema, limit.tool).toBeDefined();
            expect(schema.type === 'integer', `${limit.tool} integer`).toBe(Boolean(limit.integer));
            const lo = schema.minimum ?? -Infinity;
            const hi = schema.maximum ?? Infinity;
            if (limit.ge !== undefined) expect(lo >= limit.ge, `${limit.tool} minimum ≥ ${limit.ge}`).toBe(true);
            if (limit.gt !== undefined) expect(lo > limit.gt, `${limit.tool} minimum > ${limit.gt}`).toBe(true);
            if (limit.le !== undefined) expect(hi <= limit.le, `${limit.tool} maximum ≤ ${limit.le}`).toBe(true);
            if (limit.lt !== undefined) expect(hi < limit.lt, `${limit.tool} maximum < ${limit.lt}`).toBe(true);
        }
    });
});
