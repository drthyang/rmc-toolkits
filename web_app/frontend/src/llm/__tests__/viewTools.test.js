// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// highlight_cards: built from the host's cards, it names every card for the
// model, resolves the spellings a model uses, rejects unknown ids with the
// list, keeps marks for plot cards only, and reports what was shown.

import { describe, expect, it, vi } from 'vitest';
import { buildFocusTools, describeCards } from '../copilot/viewTools';
import { validateToolArguments } from '../copilot/toolCalls';

const CARDS = [
    { id: 'plot:GTS_PDFpartials.csv', page: 'dashboard', title: 'g(r) partials', marks: true, about: 'partial g(r)' },
    { id: 'pca-sites', page: 'ellipsoids', title: 'Site ellipsoids (3D)', about: 'sites' }
];
const PAGES = { dashboard: 'Dashboard', ellipsoids: 'PCA Ellipsoid' };

const fakeFocus = (result = {}) => ({
    show: vi.fn(async (ids) => ({ shown: ids, missing: [], page: null, switched: false, scrolled: true, ...result }))
});

describe('highlight_cards', () => {
    const [tool] = buildFocusTools({ cards: CARDS, pages: PAGES });

    it('is offered only with cards, as a view tool that lists them', () => {
        expect(buildFocusTools({ cards: [] })).toEqual([]);
        expect(tool.name).toBe('highlight_cards');
        expect(tool.view).toBe(true);
        expect(tool.description).toContain('plot:GTS_PDFpartials.csv — g(r) partials (Dashboard): partial g(r)');
        expect(describeCards(CARDS, PAGES)).toContain('pca-sites — Site ellipsoids (3D) (PCA Ellipsoid)');
        expect(tool.label({ cards: ['GTS_PDFpartials.csv'] })).toBe('Pointing at g(r) partials');
        expect(tool.label({ cards: [] })).toBe('Clear the highlight');
    });

    it('validates a list of ids and marks, also when a model sends them as text', () => {
        const checked = validateToolArguments(tool.parameters, { cards: '["pca-sites"]', marks: [{ x: '2.52', label: 'Ta–Se' }] });
        expect(checked.ok).toBe(true);
        expect(checked.value).toEqual({ cards: ['pca-sites'], marks: [{ x: 2.52, label: 'Ta–Se' }] });
        expect(validateToolArguments(tool.parameters, { cards: 'pca-sites' }).value.cards).toEqual(['pca-sites']);
        expect(validateToolArguments(tool.parameters, { cards: ['a', 'b', 'c', 'd', 'e'] }).errors[0]).toMatch(/at most 4/);
        expect(validateToolArguments(tool.parameters, { cards: [], marks: [{ label: 'x' }] }).errors[0]).toBe('marks[0].x is required');
    });

    it('shows the cards, with marks on a plot card, and reports what it pointed at', async () => {
        const focus = fakeFocus({ switched: true });
        const result = await tool.run({ cards: ['card:plot:GTS_PDFpartials.csv', 'pca-sites'], marks: [{ x: 2.52, label: 'Ta–Se' }] }, { focus });
        expect(focus.show).toHaveBeenCalledWith(['plot:GTS_PDFpartials.csv', 'pca-sites'], { marks: [{ x: 2.52, label: 'Ta–Se' }] });
        expect(result.data).toMatchObject({
            highlighted: [
                { id: 'plot:GTS_PDFpartials.csv', title: 'g(r) partials', page: 'Dashboard' },
                { id: 'pca-sites', title: 'Site ellipsoids (3D)', page: 'PCA Ellipsoid' }
            ],
            page: 'Dashboard',
            switched_page: true,
            marks: [{ x: 2.52, label: 'Ta–Se' }]
        });
        expect(result.checks).toEqual([{ level: 'ok', text: 'Pointing at g(r) partials, Site ellipsoids (3D) on the Dashboard page, 1 mark.' }]);
    });

    it('drops marks on a card that is not a plot, with a note', async () => {
        const focus = fakeFocus();
        const result = await tool.run({ cards: ['pca-sites'], marks: [{ x: 1 }] }, { focus });
        expect(focus.show).toHaveBeenCalledWith(['pca-sites'], { marks: null });
        expect(result.checks[0]).toEqual({ level: 'note', text: 'Marks are drawn on plot cards only; Site ellipsoids (3D) has none.' });
    });

    it('rejects an unknown card with the list, warns when the main card is not on screen, and clears on []', async () => {
        await expect(tool.run({ cards: ['rietveld'] }, { focus: fakeFocus() }))
            .rejects.toThrow('unknown card "rietveld"; the cards are: plot:GTS_PDFpartials.csv, pca-sites');
        await expect(tool.run({ cards: ['pca-sites'] }, {})).rejects.toThrow('the page view is not available here');

        const missing = await tool.run({ cards: ['pca-sites'] }, { focus: fakeFocus({ shown: [], missing: ['pca-sites'] }) });
        expect(missing.checks).toEqual([{ level: 'warn', text: expect.stringMatching(/^Not on screen: Site ellipsoids \(3D\)/) }]);
        expect(missing.data.not_on_screen).toEqual(['pca-sites']);

        const focus = fakeFocus();
        const cleared = await tool.run({ cards: [] }, { focus });
        expect(focus.show).toHaveBeenCalledWith([], { marks: null });
        expect(cleared).toEqual({ data: { highlighted: [] }, checks: [{ level: 'ok', text: 'Highlight cleared.' }] });
    });
});
