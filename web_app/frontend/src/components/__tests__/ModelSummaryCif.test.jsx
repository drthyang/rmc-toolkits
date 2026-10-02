// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */
/* global process */

// The Detected SG ladder's CIF download: a click on a brick only selects it (the card shows
// that group); a double-click downloads the group's symmetry-averaged CIF. The label row
// says so ("double-click for CIF"); a failed export is listed in the page's Problems
// section until an export succeeds.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { act } from 'react';
import { createRoot } from 'react-dom/client';

const downloads = vi.hoisted(() => []);
vi.mock('../../figureExport', async (importOriginal) => ({
    ...(await importOriginal()),
    downloadBlob: (blob, filename) => downloads.push({ blob, filename }),
}));

const { default: ModelSummary } = await import('../ModelSummary');
const { IssueScope, IssueStoreProvider, PageIssues } = await import('../../ui');
const { structureFromRmc6f } = await import('../../browserData');

// The jsdom environment gives import.meta.url an http scheme: read from the project root.
const RMC6F = readFileSync(resolve(process.cwd(), 'public/demo/GTS_250K.rmc6f'), 'utf8');
const demo = structureFromRmc6f({ path: 'Demo/GTS_250K.rmc6f', text: RMC6F });

describe('Detected SG CIF download', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        downloads.length = 0;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const bricks = () => [...container.querySelectorAll('.sym-brick')];
    const headline = () => container.querySelector('[aria-label="Detected space group"] dd').textContent;
    const fire = (el, type) => act(() => el.dispatchEvent(new MouseEvent(type, { bubbles: true })));
    // A double-click: two clicks, then the dblclick event.
    const doubleClick = (el) => { fire(el, 'click'); fire(el, 'click'); fire(el, 'dblclick'); };

    it('selects on a click and downloads only on a double-click', async () => {
        act(() => root.render(<ModelSummary structure={demo} />));
        expect(container.querySelector('.sym-cif-hint').textContent).toBe('double-click for CIF');

        const tetragonal = bricks().find((b) => b.textContent === 'P-42_1m');
        fire(tetragonal, 'click');
        expect(tetragonal.classList.contains('is-active')).toBe(true);
        expect(headline()).toMatch(/^P-42_1m/);
        expect(downloads).toHaveLength(0);

        const cubic = bricks().at(-1);
        doubleClick(cubic);
        expect(cubic.classList.contains('is-active')).toBe(true);
        expect(downloads).toHaveLength(1);
        expect(downloads[0].filename).toBe('GTS_250K_F-43m.cif');
        expect(downloads[0].blob.type).toBe('chemical/x-cif');
        const text = await downloads[0].blob.text();
        expect(text).toMatch(/^data_GTS_250K_F-43m$/m);
        expect(text).toMatch(/_space_group_IT_number\s+216/);
    });

    it('lists a failed export in the Problems section until one succeeds', () => {
        act(() => root.render(
            <IssueStoreProvider>
                <IssueScope page="dashboard">
                    <PageIssues />
                    <ModelSummary structure={demo} />
                </IssueScope>
            </IssueStoreProvider>,
        ));
        const problems = () => container.querySelector('.ui-issues')?.textContent ?? '';
        const cubic = bricks().at(-1);

        // Empty the basis behind the card's back so the export throws.
        const saved = demo.basis;
        demo.basis = [];
        try {
            fire(cubic, 'dblclick');
        } finally {
            demo.basis = saved;
        }
        expect(downloads).toHaveLength(0);
        expect(problems()).toMatch(/CIF · F-43m.*no average-structure basis/);

        fire(cubic, 'dblclick');
        expect(downloads).toHaveLength(1);
        expect(problems()).not.toMatch(/CIF/);
    });
});
