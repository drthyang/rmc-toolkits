// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// Above the space-group finder's site cap (MAX_SYMMETRY_SITES), Flask's
// /api/structure sends no basis, only `basisOmitted: { sites, cap }`
// (rmc_toolkits.structure_basis.MAX_BASIS_SITES, the same value). A server
// folder must then look exactly like the same files picked in the browser,
// whose structureFromRmc6f keeps every row: the Detected SG card's 'not
// analysed' state (no ladder, no CIF button), the Copilot's symmetry and run
// context, and the two symmetry tools offered with the same definitions and
// the same answers. The configuration is the case the cap is for: a 1×1×1 box
// with one reference site, of one copy, per atom.

import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';
import { structureFromRmc6f } from '../browserData';
import ModelSummary from '../components/ModelSummary';
import { copilotRun } from '../copilotRun';
import { buildRunContext, contextToJson } from '../llm/context/runContext';
import { describeSymmetry, MAX_SYMMETRY_SITES, toleranceLadder } from '../symmetryModel';

const SITES = MAX_SYMMETRY_SITES + 1;

// The box tests/test_structure_basis.py's route test sends (_one_copy_box).
const oneCopyBox = (sites) => {
    const lines = [];
    for (let number = 1; number <= sites; number += 1) {
        const [x, y, z] = [0.6180339887, 0.4142135623, 0.7320508075].map((step) => ((number * step) % 1).toFixed(10));
        lines.push(`${number} ${number % 3 ? 'Si' : 'O'} [${number}] ${x} ${y} ${z} ${number} 0 0 0`);
    }
    return [
        '(Version 6f format configuration file)',
        `Number of atoms: ${sites}`,
        'Supercell dimensions: 1 1 1',
        'Lattice vectors (Ang):',
        '40.0 0.0 0.0',
        '0.0 40.0 0.0',
        '0.0 0.0 40.0',
        'Atoms:',
        ...lines
    ].join('\n');
};

const browser = structureFromRmc6f({ path: '/runs/glass/glass.rmc6f', text: oneCopyBox(SITES) });
const { basis: _rows, ...model } = browser;
// What the Dashboard holds for the same file from /api/structure.
const server = { ...model, basis: null, basisOmitted: { sites: SITES, cap: MAX_SYMMETRY_SITES } };

const definitions = (tools) => tools.map(({ name, description, parameters }) => ({ name, description, parameters }));
const loaded = (structure) => ({
    name: '/runs/glass', files: [], structure, structurePath: structure.source, settings: null, engine: async () => ({})
});

describe('above the finder\'s site cap: a server folder looks like the browser', () => {
    it('the browser\'s basis is above the cap, and the server sent none', () => {
        expect(browser.basis).toHaveLength(SITES);
        expect(server.basis).toBeNull();
    });

    it('the same \'not analysed\' description at every tolerance, and no ladder', () => {
        for (const tolerance of [0.02, 0.2, 1]) {
            const fromBrowser = describeSymmetry(browser, tolerance);
            expect(fromBrowser).toMatchObject({ skipped: true, spaceGroup: 'not analysed', pointGroup: `${SITES} sites > ${MAX_SYMMETRY_SITES} limit` });
            expect(describeSymmetry(server, tolerance)).toEqual(fromBrowser);
        }
        expect(toleranceLadder(browser, 1.0)).toEqual([]);
        expect(toleranceLadder(server, 1.0)).toEqual([]);
    });

    it('the Copilot gets the same symmetry, run context and symmetry tools, with the same answers', async () => {
        const fromBrowser = copilotRun(loaded(browser));
        const fromServer = copilotRun(loaded(server));
        expect(fromBrowser.props.symmetry).toMatchObject({ skipped: true, ladder: [] });
        expect(fromServer.props.symmetry).toEqual(fromBrowser.props.symmetry);
        expect(contextToJson(buildRunContext(fromServer.props))).toBe(contextToJson(buildRunContext(fromBrowser.props)));

        expect(definitions(fromServer.tools)).toEqual(definitions(fromBrowser.tools));
        const names = fromServer.tools.map((tool) => tool.name);
        expect(names).toEqual(expect.arrayContaining(['detect_symmetry', 'measure_symmetry_lowering']));
        const run = (copilot, name, args) => copilot.tools.find((tool) => tool.name === name).run(args);
        for (const [name, args] of [['detect_symmetry', { tolerance_A: 0.2 }], ['measure_symmetry_lowering', {}]]) {
            const answer = await run(fromBrowser, name, args);
            expect(answer.checks).toEqual([expect.objectContaining({ level: 'warn', text: expect.stringMatching(/^Not analysed: /) })]);
            expect(await run(fromServer, name, args)).toEqual(answer);
        }
    });
});

describe('above the cap: the Detected SG card', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const card = (structure) => {
        act(() => root.render(<ModelSummary structure={structure} />));
        return container.querySelector('[aria-label="Detected space group"]');
    };

    it('renders the same \'not analysed\' card for both, with no ladder and no CIF button', () => {
        const fromBrowser = card(browser).outerHTML;
        expect(fromBrowser).toContain('not analysed');
        expect(fromBrowser).toContain(`${SITES} sites &gt; ${MAX_SYMMETRY_SITES} limit`);
        expect(container.querySelector('.sym-brick')).toBeNull();
        expect(container.querySelector('.sym-cif-action')).toBeNull();

        expect(card(server).outerHTML).toBe(fromBrowser);
        expect(container.querySelector('.sym-brick')).toBeNull();
        expect(container.querySelector('.sym-cif-action')).toBeNull();
    });
});
