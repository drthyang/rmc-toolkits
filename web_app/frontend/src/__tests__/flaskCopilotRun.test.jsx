// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// Flask mode hands the AI Copilot the run a picked folder hands it. Until
// 2026-10-19 a server folder's plot files reached the Copilot without their
// data (each chart fetched /api/plot/data for itself), so the dock offered 8
// analyses instead of 11 (no assess_convergence, find_gr_peaks, compare_fit)
// and its context had no dataset titles or Rwp and no χ² history; and the
// partial weights were labelled Se-Ta, not Ta-Se, because /api/structure sent
// elementCounts with sorted keys.
//
// The run: the ten RMCProfile output layouts of plot_parity_fixture.json and
// the Demo's .rmc6f and three χ² logs, in one folder. The browser path reads
// them as a picked folder (openRunFolder.js's readRunFolder). The Flask path
// is the real Dashboard over a mocked axios whose payloads are what the
// server sends for the same files: each fixture file's labels, title and
// metrics as the fixture pins them, its series values as the browser parses
// them (checked against every value the fixture pins), sent through JSON
// (NaN → null); the logs' one concatenated history (the server's
// related_r_value_logs); the structure with the server's basis
// (structure_basis_fixture.json) and the species in the .rmc6f's order
// (tests/test_parsers_demo_run.py pins that order on the server's wire).
// The real InteractivePlot draws the charts, so the requests it shares with
// the Copilot's run are counted too. No server, no model.

import { afterEach, beforeAll, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import parity from './fixtures/plot_parity_fixture.json';
import basisFixture from './fixtures/structure_basis_fixture.json';

const state = vi.hoisted(() => ({
    listings: {},
    metadata: new Map(),
    payloads: new Map(),
    structures: {},
    // GET /api/plot/data, in order: the path asked for.
    dataCalls: [],
    // Set: each /api/plot/data call waits for its resolver here.
    hold: false,
    held: [],
}));

vi.mock('axios', () => ({
    default: {
        get: vi.fn(async (url, { params = {} } = {}) => {
            if (url.endsWith('/api/files')) return { data: { files: state.listings[params.dir] || [] } };
            if (url.endsWith('/api/plot/metadata')) return { data: state.metadata.get(params.path) };
            if (url.endsWith('/api/structure')) {
                const structure = state.structures[params.dir];
                if (!structure) throw Object.assign(new Error('404'), { response: { data: { error: 'No .rmc6f file found' } } });
                return { data: structure };
            }
            if (url.endsWith('/api/plot/data')) {
                state.dataCalls.push(params.path);
                const answer = () => {
                    const payload = state.payloads.get(params.path);
                    if (payload instanceof Error) throw payload;
                    return { data: payload };
                };
                if (state.hold) return new Promise((done, fail) => state.held.push(() => { try { done(answer()); } catch (error) { fail(error); } }));
                return answer();
            }
            throw new Error(`unexpected ${url}`);
        }),
    },
}));

// Flask mode (the test environment would otherwise count as static).
vi.mock('../browserData', async (importOriginal) => ({
    ...(await importOriginal()),
    isStaticMode: () => false,
}));
vi.mock('../llm', async (importOriginal) => ({ ...(await importOriginal()), WatchdogBadge: () => null }));

const { default: Dashboard } = await import('../components/Dashboard');
const {
    CHI_HISTORY_Y_LABEL, chiHistoryLabels, combineChiTerms, detectPlotKind, plotDataFromText, structureFromRmc6f, WATCH_INTERVAL_MS
} = await import('../browserData');
const { readRunFolder } = await import('../openRunFolder');
const { copilotRun, DEFAULT_VIEW } = await import('../copilotRun');
const { buildCopilotCards } = await import('../copilotCards');
const { buildRunContext, contextToJson } = await import('../llm/context/runContext');
const { buildRunTools } = await import('../llm/copilot/runTools');
const { clearPlotDataCache } = await import('../plotDataCache');
const { IssueScope, IssueStoreProvider } = await import('../ui');

const HERE = dirname(fileURLToPath(import.meta.url));
const DEMO = join(HERE, '../../public/demo');
const TIME = Date.UTC(2026, 0, 1);
const FOLDER = 'Run';
const SERVER_DIR = `/srv/runs/${FOLDER}`;
const LOGS = ['GTS_250K-00.log', 'GTS_250K-01.log', 'GTS_250K-02.log'];
const RMC6F = 'GTS_250K.rmc6f';

// The run's files: the fixture's ten layouts and the Demo's model and logs.
const SOURCES = [
    ...parity.cases.map(({ name, text }) => ({ name, text })),
    ...[RMC6F, ...LOGS].map((name) => ({ name, text: readFileSync(join(DEMO, name), 'utf8') })),
];
const entries = () => SOURCES.map(({ name, text }) => ({
    path: `${FOLDER}/${name}`, file: new File([text], name, { lastModified: TIME }),
}));
const textOf = (name) => SOURCES.find((source) => source.name === name).text;
const serverPath = (name, dir = SERVER_DIR) => `${dir}/${name}`;

// /api/files: what list_files() sends (path, kind, mtime in seconds, size).
const listing = (dir = SERVER_DIR, modified = {}) => SOURCES.map(({ name, text }) => ({
    name,
    path: serverPath(name, dir),
    type: 'file',
    plotKind: detectPlotKind(name),
    modified: modified[name] ?? TIME / 1000,
    size: new TextEncoder().encode(text).length,
}));

const throughJson = (value) => JSON.parse(JSON.stringify(value));

// A fixture file's /api/plot/data payload: the fixture's (Flask's) kind,
// title, labels and metrics; the values as the browser parses them, which
// agree with every value the fixture pins (checked below).
const fixturePayload = ({ name, text, expected }) => {
    const parsed = plotDataFromText({ name, plotKind: detectPlotKind(name), text });
    return throughJson({
        kind: expected.kind,
        title: expected.title,
        metrics: expected.metrics,
        xLabel: expected.xLabel,
        yLabel: expected.yLabel,
        series: expected.series.map((series, index) => ({ label: series.label, x: parsed.series[index].x, y: parsed.series[index].y })),
    });
};

// A log's /api/plot/data payload: the history of every log of its run, in
// sequence (related_r_value_logs + read_chi_log), as app.py's plot_data sends it,
// with every χ² term and their total (chi_terms_payload; its values pinned by
// test_parsers_demo_run.py's test_log_data_carries_every_chi_term_and_their_total).
// test_parsers_demo_run.py's log_last_column: the last token of each row below the two header lines.
const lastColumn = (text) => text.split(/\r\n|\r|\n/).slice(2).filter((line) => line.trim()).map((line) => Number(line.trim().split(/\s+/).pop()));
const logPayload = () => {
    const parsed = LOGS.map((name) => plotDataFromText({ name, plotKind: 'r_value', text: textOf(name) }));
    const column = parsed[0].chiColumn;
    const { title, label } = chiHistoryLabels(column);
    const y = parsed.flatMap((plot) => plot.series[0].y);
    return throughJson({
        kind: 'r_value',
        title,
        metrics: { final_chi_r: parsed[parsed.length - 1].metrics.final_chi_r },
        xLabel: 'Time steps',
        yLabel: CHI_HISTORY_Y_LABEL,
        chiColumn: column,
        chiTerms: combineChiTerms(parsed.map((plot) => plot.chiTerms)),
        series: [{ label, x: y.map((_, index) => index), y }],
    });
};

// /api/structure: the model as the parsers read it, the server's basis rows,
// and the species in the .rmc6f's order ('Atom types present', as the server
// builds elementCounts and atomIndices, and keeps them on the wire).
const header = (key) => textOf(RMC6F).split('\n').find((line) => line.startsWith(key)).split(':')[1].trim().split(/\s+/);
const ATOM_TYPES = header('Atom types present');
const serverStructure = (dir = SERVER_DIR) => {
    const { basis: _browserRows, elementCounts, atomIndices, ...model } = structureFromRmc6f({ path: serverPath(RMC6F, dir), text: textOf(RMC6F) });
    return throughJson({
        ...model,
        elementCounts: Object.fromEntries(ATOM_TYPES.map((element) => [element, elementCounts[element]])),
        atomIndices: Object.fromEntries(ATOM_TYPES.map((element) => [element, atomIndices[element]])),
        basis: basisFixture.cases[0].basis,
        basisOmitted: null,
    });
};

const serveFolder = (dir = SERVER_DIR) => {
    state.listings[dir] = listing(dir);
    state.structures[dir] = serverStructure(dir);
    parity.cases.forEach((entry) => {
        const payload = fixturePayload(entry);
        state.payloads.set(serverPath(entry.name, dir), payload);
        state.metadata.set(serverPath(entry.name, dir), { kind: payload.kind, title: entry.expected.metadataTitle, metrics: payload.metrics });
    });
    const history = logPayload();
    LOGS.forEach((name) => {
        state.payloads.set(serverPath(name, dir), history);
        state.metadata.set(serverPath(name, dir), { kind: 'r_value', title: history.title, metrics: history.metrics });
    });
};

const ENGINE = async () => ({});
const definitions = (tools) => tools.map(({ name, description, parameters }) => ({ name, description, parameters }));
// What the dock is handed for a run: the host's tools and the module's run
// tools over its props (the 11 analyses), the cards, and the run context.
const offered = (run) => {
    const cards = buildCopilotCards({ plotFiles: run.props.plotFiles, rValueFile: run.props.rValueFile, symmetry: run.props.symmetry });
    return {
        cards,
        tools: [...run.tools, ...buildRunTools(run.props, { cards })],
        context: buildRunContext(run.props),
    };
};

const flush = () => act(async () => {
    for (let i = 0; i < 10; i += 1) await new Promise((resolve) => setTimeout(resolve, 0));
});
const until = async (condition, what) => {
    for (let i = 0; i < 50 && !condition(); i += 1) await flush();
    if (!condition()) throw new Error(`timed out waiting for ${what}`);
};

const plotFileNames = parity.cases.map(({ name }) => name);
const allFetched = (published) => Boolean(published?.files && published.structure) && plotFileNames.every(
    (name) => published.files.find((file) => file.name === name)?.plotData
) && published.files.find((file) => file.name === LOGS[0])?.plotData;

describe('the fixture\'s Flask payloads are what the server sends', () => {
    it.each(parity.cases.map((entry) => [entry.name, entry]))('%s: every pinned value', (name, entry) => {
        const payload = fixturePayload(entry);
        payload.series.forEach((series, index) => {
            const finite = series.y.filter(Number.isFinite);
            const pinned = entry.expected.series[index];
            expect([series.label, series.y.length, finite.length]).toEqual([pinned.label, pinned.n, pinned.finite]);
            if (pinned.first === null) expect(finite).toEqual([]);
            else {
                expect(finite[0]).toBeCloseTo(pinned.first, 12);
                expect(finite[finite.length - 1]).toBeCloseTo(pinned.last, 12);
            }
        });
    });

    it('the logs: one history of the three restarts, as test_parsers_demo_run.py pins it', () => {
        const history = logPayload();
        const rows = LOGS.flatMap((name) => lastColumn(textOf(name)));
        expect(history.title).toBe('χ² history: X_ray_(R)1');
        expect(history.series[0].y).toHaveLength(rows.length);
        expect(history.metrics.final_chi_r).toBe(rows[rows.length - 1]);
        expect(history.series[0].y[rows.length - 1]).toBeCloseTo(Math.log(rows[rows.length - 1]), 12);
    });

    it('the structure: the species in the .rmc6f\'s order, as the browser builds them', () => {
        const browser = structureFromRmc6f({ path: `${FOLDER}/${RMC6F}`, text: textOf(RMC6F) });
        expect(ATOM_TYPES).toEqual(['Ga', 'Ta', 'Se']);
        expect(Object.keys(browser.elementCounts)).toEqual(ATOM_TYPES);
        expect(Object.keys(serverStructure().elementCounts)).toEqual(ATOM_TYPES);
    });
});

describe('Flask mode: the Copilot gets the run a picked folder gives it', () => {
    let container;
    let root;
    let published;
    let browserLoaded;

    beforeAll(async () => {
        browserLoaded = await readRunFolder(entries());
    });

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        clearPlotDataCache();
        state.listings = {};
        state.metadata = new Map();
        state.payloads = new Map();
        state.structures = {};
        state.dataCalls = [];
        state.hold = false;
        state.held = [];
        serveFolder();
        published = null;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.useRealTimers();
    });

    const render = async (props = {}) => {
        act(() => root.render(
            <IssueStoreProvider>
                <IssueScope page="dashboard">
                    <Dashboard
                        directory={SERVER_DIR}
                        localRun={null}
                        wantAssistantData
                        onRunChange={(run) => { published = run; }}
                        {...props}
                    />
                </IssueScope>
            </IssueStoreProvider>
        ));
        await flush();
    };
    const flaskRun = () => copilotRun({
        name: FOLDER, files: published.files, structure: published.structure,
        structurePath: published.structurePath, settings: null, engine: ENGINE,
    }, DEFAULT_VIEW);
    const browserRun = () => copilotRun({ ...browserLoaded, name: FOLDER, engine: ENGINE }, DEFAULT_VIEW);
    const callsFor = (name) => state.dataCalls.filter((path) => path.endsWith(`/${name}`)).length;

    it('the same tools, cards and run context as the browser path', async () => {
        await render();
        await until(() => allFetched(published), 'the published plot data');
        const flask = offered(flaskRun());
        const browser = offered(browserRun());

        const names = flask.tools.map((tool) => tool.name);
        expect(names).toEqual(browser.tools.map((tool) => tool.name));
        expect(names).toEqual(expect.arrayContaining(['assess_convergence', 'find_gr_peaks', 'compare_fit', 'read_run_data']));
        expect(names).toHaveLength(11);
        expect(definitions(flask.tools)).toEqual(definitions(browser.tools));
        expect(flask.cards).toEqual(browser.cards);

        // The context keys, top level and nested, then the whole text.
        const keys = (context) => Object.fromEntries(Object.entries(context).map(([key, value]) => [
            key, value && typeof value === 'object' && !Array.isArray(value) ? Object.keys(value) : typeof value,
        ]));
        expect(keys(flask.context)).toEqual(keys(browser.context));
        expect(Object.keys(flask.context)).toEqual(expect.arrayContaining(['datasets', 'convergence', 'partial_weights', 'pair_correlations', 'symmetry']));
        expect(flask.context.datasets.map(Object.keys)).toEqual(browser.context.datasets.map(Object.keys));
        expect(flask.context.partial_weights).toEqual(browser.context.partial_weights);
        expect(Object.keys(flask.context.partial_weights.xray)).toContain('Ta-Se');
        expect(Object.keys(flask.context.partial_weights.neutron)).toContain('Ta-Se');
        expect(contextToJson(flask.context)).toBe(contextToJson(browser.context));

        // The values the fixture pins: each dataset's title and Rwp.
        const roundSig = (value) => Number(value.toPrecision(3));
        const byKind = new Map();
        parity.cases.filter(({ expected }) => expected.kind !== 'pdf_partials').forEach(({ expected }) => {
            if (!byKind.has(expected.kind)) byKind.set(expected.kind, []);
            byKind.get(expected.kind).push(expected);
        });
        byKind.forEach((pinned, kind) => {
            const entries = flask.context.datasets.filter((entry) => entry.kind === kind);
            expect(entries.map((entry) => entry.title ?? kind).sort()).toEqual(pinned.map((expected) => expected.title).sort());
            expect(entries.map((entry) => entry.rwp ?? null).sort()).toEqual(
                pinned.map((expected) => (Number.isFinite(expected.metrics.rwp) ? roundSig(expected.metrics.rwp) : null)).sort()
            );
        });
        expect(flask.context.convergence).toMatchObject({ n_steps: logPayload().series[0].y.length, column: 'X_ray_(R)1' });
        // Every χ² term and their total, as the browser reads the three logs (null → NaN).
        const terms = flaskRun().props.rValueFile.plotData.chiTerms;
        expect(terms).toEqual(browserRun().props.rValueFile.plotData.chiTerms);
        expect(terms.segments).toEqual([1320, 510, 432]);
    });

    it('one request per file and version, shared with the charts; none before the Copilot opens', async () => {
        await render({ wantAssistantData: false });
        await until(() => plotFileNames.every((name) => callsFor(name) === 1), 'the charts\' requests');
        // Only the charts on screen asked: the χ² panel starts collapsed.
        expect(callsFor(LOGS[0])).toBe(0);
        expect(published).toBeNull();

        // The Copilot opens: the run reuses every chart's answer and asks only
        // for the χ² history (its run's first log, not each log).
        await render();
        await until(() => allFetched(published), 'the published plot data');
        plotFileNames.forEach((name) => expect(callsFor(name)).toBe(1));
        expect(LOGS.map(callsFor)).toEqual([1, 0, 0]);

        // Expanding the χ² panel draws the run's answer: no second request.
        const show = [...container.querySelectorAll('[data-copilot-card="chi2-history"] button')].find((button) => button.textContent === 'Show');
        act(() => show.dispatchEvent(new MouseEvent('click', { bubbles: true })));
        await flush();
        expect(container.querySelector('[data-copilot-card="chi2-history"] .interactive-plot svg')).not.toBeNull();
        expect(LOGS.map(callsFor)).toEqual([1, 0, 0]);
    });

    it('Live Data: only the files whose version moved are asked again', async () => {
        vi.useFakeTimers({ toFake: ['setInterval', 'clearInterval'] });
        await render({ watchFiles: true });
        await until(() => allFetched(published), 'the published plot data');
        const before = published.files;
        state.dataCalls = [];

        // RMCProfile rewrites one fit and appends to its current (last) log.
        const changed = 'NPDF_PDF1.csv';
        state.listings[SERVER_DIR] = listing(SERVER_DIR, { [changed]: TIME / 1000 + 60, [LOGS[2]]: TIME / 1000 + 60 });
        const fresh = { ...state.payloads.get(serverPath(changed)), title: 'PDF1 (rewritten)' };
        state.payloads.set(serverPath(changed), fresh);
        await act(async () => { vi.advanceTimersByTime(WATCH_INTERVAL_MS); });
        await until(() => published.files.find((file) => file.name === changed)?.plotData?.title === 'PDF1 (rewritten)', 'the new answer');

        // The fit, and the χ² history through its run's first log (the run's
        // version moved with its last log); nothing else.
        expect([...state.dataCalls].sort()).toEqual([serverPath(changed), serverPath(LOGS[0])].sort());
        const plotDataOf = (files, name) => files.find((file) => file.name === name).plotData;
        plotFileNames.filter((name) => name !== changed).forEach((name) => {
            expect(plotDataOf(published.files, name)).toBe(plotDataOf(before, name));
        });
    });

    it('a folder switch drops the folder before\'s answers', async () => {
        const OTHER = '/srv/runs/Other';
        serveFolder(OTHER);
        state.hold = true;
        await render();
        await until(() => state.held.length > 0, 'the first folder\'s requests');
        const firstFolder = state.held.splice(0);

        // Another folder opens before the first folder's answers land.
        await render({ directory: OTHER });
        await until(() => published?.files?.[0]?.path.startsWith(OTHER) && state.held.length > 0, 'the other folder\'s listing');
        await act(async () => { firstFolder.forEach((answer) => answer()); });
        await flush();
        expect(published.files.every((file) => file.path.startsWith(OTHER))).toBe(true);
        expect(published.files.some((file) => file.plotData)).toBe(false);

        // Its own answers do land.
        state.hold = false;
        await act(async () => { state.held.splice(0).forEach((answer) => answer()); });
        await until(() => allFetched(published), 'the other folder\'s plot data');
        expect(published.files.every((file) => file.path.startsWith(OTHER))).toBe(true);
    });

    it('a request that fails is in the Problems section and in the run as a parseError', async () => {
        // The χ² history: its panel is collapsed, so no chart lists the failure.
        state.payloads.set(serverPath(LOGS[0]), Object.assign(new Error('500'), { response: { data: { error: 'GTS_250K-00.log does not contain chi values' } } }));
        await render();
        await until(() => published?.files?.find((file) => file.name === LOGS[0])?.parseError, 'the failure');
        expect(container.querySelector('[role="alert"]').textContent).toBe('GTS_250K-00.log: GTS_250K-00.log does not contain chi values');
        const run = flaskRun();
        expect(run.props.rValueFile.parseError).toBe('GTS_250K-00.log does not contain chi values');
        expect(offered(run).tools.map((tool) => tool.name)).not.toContain('assess_convergence');
    });
});
