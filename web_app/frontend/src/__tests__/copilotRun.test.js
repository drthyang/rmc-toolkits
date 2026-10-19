// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The run as the Copilot sees it, at its interface: the pure core
// (copilotRun.js) and its in-process adapter (openRunFolder.js) on the Demo,
// a Flask-shaped run, the view, the errors, the drift between the adapter and
// the Dashboard's own parses, the import guard that keeps the adapter out of
// the app, and every copy of the Detected SG card's default tolerance. The
// browser adapter (useCopilotRun) is checked in useCopilotRun.test.jsx.

import { readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { ESLint } from 'eslint';
import { beforeAll, describe, expect, it } from 'vitest';
import { ARM_A_RULES } from '../../validation/armA';
import { COPILOT_RUN_PROPS, DEFAULT_VIEW, copilotRun, normalizeView } from '../copilotRun';
import { openRunFolder, readRunFolder } from '../openRunFolder';
import basisFixture from './fixtures/structure_basis_fixture.json';
import { demoEntries, loadDemoAssistantRun } from './modelVisible/helpers';

const HERE = dirname(fileURLToPath(import.meta.url));
const FRONTEND = join(HERE, '../..');
const read = (path) => readFileSync(join(FRONTEND, path), 'utf8');

// A tool's definition, as the model is sent it (its functions dropped).
const definitions = (tools) => tools.map(({ name, description, parameters }) => ({ name, description, parameters }));

const DEMO_PLOTS = [
    'GTS_250K-00.log', 'GTS_250K-01.log', 'GTS_250K-02.log', 'GTS_250K_FQ1.csv', 'GTS_250K_FT_XFQ1.csv', 'GTS_250K_PDFpartials.csv'
];

describe('the Demo through openRunFolder (the in-process adapter)', () => {
    let run;
    let loaded;
    beforeAll(async () => {
        loaded = await readRunFolder(demoEntries());
        run = await openRunFolder(demoEntries());
    });

    it('names the run after its folder and lists the Dashboard\'s plot files in its order', () => {
        expect(run.props.runName).toBe('Demo');
        expect(run.props.plotFiles.map((file) => file.name)).toEqual(DEMO_PLOTS);
        expect(run.props.plotFiles.every((file) => file.plotData && !file.parseError)).toBe(true);
        expect(run.props.plotFiles.find((file) => file.name === 'GTS_250K_PDFpartials.csv').plotKind).toBe('pdf_partials');
    });

    it('charts the three logs as one χ² history, and reads the .dat', () => {
        const history = run.props.rValueFile;
        expect(history.sourceNames).toEqual(['GTS_250K-00.log', 'GTS_250K-01.log', 'GTS_250K-02.log']);
        expect(history.plotData.series[0].y.length).toBeGreaterThan(10);
        expect(run.props.runSettings).not.toBeNull();
        expect(run.props.runSettings).toEqual(loaded.settings);
    });

    it('gives the Detected SG card at 0.2 Å with its ladder, the elements and every engine tool', () => {
        expect(run.props.symmetry).toMatchObject({ spaceGroup: 'F-43m', toleranceA: 0.2 });
        expect(run.props.symmetry.ladder.length).toBeGreaterThan(1);
        expect(run.props.structure.source).toBe('Demo/GTS_250K.rmc6f');
        expect(run.elements).toEqual(['Ga', 'Se', 'Ta']);
        expect(run.tools.map((tool) => tool.name)).toEqual([
            'list_sites', 'analyze_displacement_directions', 'compute_bond_angles', 'analyze_clusters',
            'analyze_modulations', 'detect_symmetry', 'measure_symmetry_lowering'
        ]);
        expect(run.view).toEqual(DEFAULT_VIEW);
    });

    it('runs the engine in-thread with the worker\'s message, and refuses an aborted request', async () => {
        const sites = await run.engine('sites', {});
        expect(sites.sites.length).toBeGreaterThan(0);
        // The worker's precedence: params never override the kind.
        expect(await run.engine('sites', { kind: 'kde' })).toEqual(sites);
        const controller = new AbortController();
        controller.abort();
        await expect(run.engine('sites', {}, { signal: controller.signal })).rejects.toMatchObject({ name: 'AbortError' });
    });

    it('drift: the adapter hands the Copilot what the Dashboard\'s own parses give the core', async () => {
        // The Demo button's files as the Dashboard parses them (helpers.js:
        // parseDashboardPlotFile, the structure worker's handler, the .dat).
        const { copilot: dashboard } = await loadDemoAssistantRun();
        expect(run.props).toEqual(dashboard.props);
        expect(run.elements).toEqual(dashboard.elements);
        expect(run.view).toEqual(dashboard.view);
        expect(definitions(run.tools)).toEqual(definitions(dashboard.tools));
    });
});

describe('copilotRun (the pure core)', () => {
    let loaded;
    beforeAll(async () => {
        loaded = await readRunFolder(demoEntries());
    });

    it('hands CopilotChat exactly its eight run props, frozen', () => {
        const run = copilotRun(loaded);
        // The list AGENTS.md states for the src/llm/ import boundary.
        const agents = readFileSync(join(FRONTEND, '../../AGENTS.md'), 'utf8').replace(/\s+/g, ' ');
        const listed = /receives run data \*\*only as props\*\* \(([^)]*)\)/.exec(agents)?.[1] ?? '';
        const names = [...listed.matchAll(/`(\w+)`/g)].map((match) => match[1]);
        expect(names).toHaveLength(8);
        expect(Object.keys(run.props).sort()).toEqual([...names].sort());
        expect([...COPILOT_RUN_PROPS].sort()).toEqual([...names].sort());
        expect(Object.isFrozen(run.props)).toBe(true);
        expect(Object.keys(run).sort()).toEqual(['elements', 'engine', 'props', 'tools', 'view']);
    });

    it('a tolerance changes the symmetry and not the ladder', () => {
        const at02 = copilotRun(loaded).props.symmetry;
        // A tighter rung of the ladder, where another group holds.
        const rung = at02.ladder.find((brick) => brick.spaceGroup && brick.spaceGroup !== at02.spaceGroup);
        expect(rung).toBeDefined();
        const tolerance = (rung.from + rung.to) / 2;
        const tighter = copilotRun(loaded, { toleranceA: tolerance }).props.symmetry;
        expect(tighter.toleranceA).toBe(tolerance);
        expect(tighter.spaceGroup).toBe(rung.spaceGroup);
        expect(tighter.ladder).toEqual(at02.ladder);
    });

    it('a hidden log drops out of the χ² history, as on screen', () => {
        const run = copilotRun(loaded, { hiddenPlots: ['Demo/GTS_250K-01.log'] });
        expect(run.props.rValueFile.sourceNames).toEqual(['GTS_250K-00.log', 'GTS_250K-02.log']);
        // The plot list still names every file: the Copilot may read a hidden one.
        expect(run.props.plotFiles.map((file) => file.name)).toEqual(DEMO_PLOTS);
        expect(run.view.hiddenPlots).toEqual(['Demo/GTS_250K-01.log']);
    });

    it('the ladder off gives an empty ladder and the same group', () => {
        const run = copilotRun(loaded, { ladder: false });
        expect(run.props.symmetry.ladder).toEqual([]);
        expect(run.props.symmetry.spaceGroup).toBe(copilotRun(loaded).props.symmetry.spaceGroup);
    });

    // A server folder: /api/files' listing (plot kinds from Flask, nothing
    // parsed) and /api/structure's model.
    const serverFiles = [
        { name: 'run.rmc6f', path: '/runs/a/run.rmc6f', type: 'file', plotKind: null },
        { name: 'run_PDFpartials.csv', path: '/runs/a/run_PDFpartials.csv', type: 'file', plotKind: 'pdf_partials' },
        { name: 'run-01.log', path: '/runs/a/run-01.log', type: 'file', plotKind: 'r_value' }
    ];
    const serverRun = (basis) => {
        const { basis: _basis, ...model } = loaded.structure;
        const structure = basis === undefined ? model : { ...model, basis };
        const engine = async () => ({});
        return { engine, run: copilotRun({ name: '/runs/a', files: serverFiles, structure, structurePath: '/runs/a/run.rmc6f', settings: null, engine }) };
    };
    const tool = (run, name) => run.tools.find((entry) => entry.name === name);

    it('a server folder: the basis of /api/structure gives the symmetry and both symmetry tools', async () => {
        // The basis Flask computes (rmc_toolkits/structure_basis.py: the golden's demo rows).
        const { run, engine } = serverRun(basisFixture.cases[0].basis);
        const browser = copilotRun(loaded);
        expect(run.props).toMatchObject({ runName: '/runs/a', runSettings: null });
        expect(run.props.plotFiles.map((file) => file.name)).toEqual(['run-01.log', 'run_PDFpartials.csv']);
        expect(run.props.rValueFile.name).toBe('run-01.log');
        expect(run.engine).toBe(engine);
        expect(run.props.symmetry).toMatchObject({ spaceGroup: 'F-43m', toleranceA: 0.2 });
        expect(run.props.symmetry.ladder.map((brick) => brick.spaceGroup)).toEqual(browser.props.symmetry.ladder.map((brick) => brick.spaceGroup));
        // Every tool the same files picked in the browser get, as the model is sent them...
        expect(definitions(run.tools)).toEqual(definitions(browser.tools));
        expect(run.tools.map((entry) => entry.name)).toEqual(expect.arrayContaining(['detect_symmetry', 'measure_symmetry_lowering']));
        // ...and the symmetry tools' answers.
        for (const args of [{ tolerance_A: 0.2 }, { tolerance_A: 0.02 }]) {
            expect(await tool(run, 'detect_symmetry').run(args)).toEqual(await tool(browser, 'detect_symmetry').run(args));
        }
        expect(await tool(run, 'measure_symmetry_lowering').run({})).toEqual(await tool(browser, 'measure_symmetry_lowering').run({}));
    });

    it('a server folder without a basis (a coords-only file, an older server): no symmetry, no symmetry tool', () => {
        // Flask answers `basis: []` for a file with no reference numbers; a
        // server before 0.8 sent no basis at all.
        for (const basis of [[], undefined]) {
            const { run, engine } = serverRun(basis);
            expect(run.props).toMatchObject({ runName: '/runs/a', symmetry: null, runSettings: null });
            expect(run.engine).toBe(engine);
            const names = run.tools.map((entry) => entry.name);
            expect(names).toContain('list_sites');
            expect(names).not.toContain('detect_symmetry');
            expect(names).not.toContain('measure_symmetry_lowering');
        }
    });

    it('never throws on a missing piece; a non-positive tolerance is a RangeError', () => {
        for (const missing of [null, {}, { name: 'empty', files: [] }]) {
            const run = copilotRun(missing);
            expect(run.props).toMatchObject({ plotFiles: [], rValueFile: null, structure: null, symmetry: null, runSettings: null });
            expect(run.tools).toEqual([]);
            expect(run.engine).toBeNull();
        }
        for (const toleranceA of [0, -0.1, Number.NaN, Infinity, '0.2']) {
            expect(() => copilotRun(loaded, { toleranceA })).toThrow(RangeError);
            expect(() => normalizeView({ toleranceA })).toThrow(RangeError);
        }
    });

    it('the view comes back whole, normalised and JSON-safe', () => {
        const view = normalizeView({ toleranceA: 0.1, hiddenPlots: ['b', 'a', 'b'], liveData: 1, pcaSites: { sites: [] } });
        expect(view).toEqual({ toleranceA: 0.1, hiddenPlots: ['a', 'b'], pcaSites: { sites: [] }, liveData: true, ladder: true });
        expect(JSON.parse(JSON.stringify(view))).toEqual(view);
        expect(JSON.parse(JSON.stringify(copilotRun(loaded).view))).toEqual(DEFAULT_VIEW);
        expect(Object.isFrozen(view)).toBe(true);
        // Normalising a normalised view keeps its parts.
        expect(normalizeView(view).hiddenPlots).toBe(view.hiddenPlots);
        expect(Object.isFrozen(DEFAULT_VIEW)).toBe(true);
    });

    it('a PCA table and Live Data pass through to the props, nothing else', () => {
        const table = { sites: [{ referenceNumber: 1 }] };
        const plain = copilotRun(loaded);
        const run = copilotRun(loaded, { pcaSites: table, liveData: true });
        expect(run.props.pcaSites).toBe(table);
        expect(run.props.liveData).toBe(true);
        expect({ ...run.props, pcaSites: null, liveData: false }).toEqual(plain.props);
    });
});

describe('openRunFolder errors', () => {
    const file = (path, text) => ({ path, file: new File([text], path.split('/').pop(), { lastModified: 0 }) });

    it('rejects a folder with no RMCProfile file, and a bad view before reading', async () => {
        await expect(openRunFolder([file('run/notes.txt', 'hello')])).rejects.toThrow(/^No supported RMCProfile files/);
        await expect(openRunFolder(demoEntries(), { toleranceA: 0 })).rejects.toThrow(RangeError);
    });

    it('rejects with the .rmc6f\'s parse error', async () => {
        const broken = 'Number of atoms: 1\nCell (Ang/deg): 5 5 5 90 90 90\nAtoms:\n   1   Ga  not a number\n';
        await expect(openRunFolder([file('run/run.rmc6f', broken)])).rejects.toThrow('run/run.rmc6f is missing lattice or supercell metadata');
    });

    it('keeps a plot file\'s parse error, as the Dashboard does', async () => {
        const [rmc6f] = demoEntries().filter((entry) => entry.path.endsWith('.rmc6f'));
        const run = await openRunFolder([
            { path: 'run/GTS_250K.rmc6f', file: rmc6f.file },
            file('run/GTS_250K_PDFpartials.csv', 'not,a\nplot,file\n')
        ]);
        const [partials] = run.props.plotFiles;
        expect(partials.plotData).toBeUndefined();
        expect(partials.parseError).toEqual(expect.any(String));
        expect(run.props.symmetry).not.toBeNull();
    });
});

describe('the import guard: no app module imports the in-process adapter', () => {
    const eslint = new ESLint({ cwd: FRONTEND });
    const lint = async (filePath, code) => {
        const [result] = await eslint.lintText(code, { filePath: join(FRONTEND, filePath) });
        return result.messages.filter((message) => ['no-restricted-imports', 'no-restricted-syntax'].includes(message.ruleId));
    };

    it('rejects an import of openRunFolder from src/', async () => {
        expect(await lint('src/someModule.js', 'import { openRunFolder } from \'./openRunFolder\';\nexport default openRunFolder;\n'))
            .toEqual([expect.objectContaining({ ruleId: 'no-restricted-imports', severity: 2 })]);
        expect(await lint('src/components/Page.jsx', 'export { openRunFolder } from \'../openRunFolder.js\';\n'))
            .toEqual([expect.objectContaining({ ruleId: 'no-restricted-imports', severity: 2 })]);
        expect(await lint('src/someModule.js', 'export const load = () => import(\'./openRunFolder\');\n'))
            .toEqual([expect.objectContaining({ ruleId: 'no-restricted-syntax', severity: 2 })]);
    });

    it('allows it in tests and in validation/', async () => {
        const code = 'import { openRunFolder } from \'../openRunFolder\';\nexport default openRunFolder;\n';
        expect(await lint('src/__tests__/someTest.test.js', code)).toEqual([]);
        expect(await lint('validation/someHarness.js', 'import { openRunFolder } from \'../src/openRunFolder\';\nexport default openRunFolder;\n')).toEqual([]);
    });
});

describe('the Detected SG card\'s default tolerance is one number: DEFAULT_VIEW.toleranceA', () => {
    const tolerance = DEFAULT_VIEW.toleranceA;

    it('the app starts the shared tolerance there, and the copies the module cannot import agree', async () => {
        // ...and back there when another run opens (appStateLifecycle.test.jsx).
        expect(read('src/App.jsx')).toContain('const symTolState = useRunTolerance(runKey, DEFAULT_VIEW.toleranceA);');
        // Arm A reads the card at its default (a reading rule fixed before any run).
        expect(ARM_A_RULES.CARD_TOLERANCE_A).toBe(tolerance);
        // ModelSummary's fallback without a provider, and the Workbench's
        // tolerance question (llm module: it may not import the host), name
        // the card's default as a literal.
        expect(read('src/components/ModelSummary.jsx')).toContain(`const localSymTol = useState(${tolerance});`);
        expect(read('src/llm/copilot/clarify.js')).toContain(`brick.from <= ${tolerance} && ${tolerance} < brick.to`);
    });

    it('the symmetry tools state it to the model and default to it', async () => {
        const run = await openRunFolder(demoEntries());
        const tools = Object.fromEntries(run.tools.map((tool) => [tool.name, tool]));
        expect(tools.detect_symmetry.parameters.properties.tolerance_A.description).toContain(`the card default is ${tolerance}.`);
        const lowering = tools.measure_symmetry_lowering;
        expect(lowering.description).toContain(`default ${tolerance} Å, the Detected SG card's`);
        expect(lowering.parameters.properties.parent_tolerance_A.description).toContain(`Default ${tolerance}.`);
        expect(await lowering.run({})).toEqual(await lowering.run({ parent_tolerance_A: tolerance }));
    });
});
