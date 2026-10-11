// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// orient_3d_view and current_view against fake 3D views in the registry:
// the card is brought on screen first, the direction is read in the view's
// cell (with c or b up), sites are selected (waiting for the view to settle),
// centred and ringed — the rings ending with the card's focus — and every
// refusal names what was wrong. current_view reports the page on screen.

import { afterEach, describe, expect, it, vi } from 'vitest';
import { buildViewTools } from '../copilotViewTools';
import { buildCopilotCards } from '../copilotCards';
import { normalize, registerViewer } from '../viewer3d';

const CUBIC = [[10.4, 0, 0], [0, 10.4, 0], [0, 0, 10.4]];
const PAGES = { orientation: 'Displacement Directions', ellipsoids: 'PCA Ellipsoid', structure: 'Atomic Density' };

const fakeViewer = (overrides = {}) => {
    const sites = [{ ref: 1, element: 'Ga' }, { ref: 5, element: 'Ta' }];
    let selected = 1;
    const viewer = {
        ready: () => true,
        has: (name) => typeof viewer[name] === 'function',
        cell: () => CUBIC,
        pcaAxes: () => [[0, 1, 0], [1, 0, 0], [0, 0, 1]],
        describe: () => ({ view_along: '[1 1 1]', site: selected }),
        lookAlong: vi.fn(() => true),
        reset: vi.fn(() => true),
        sites: () => sites,
        centerOn: vi.fn((ref) => sites.some((site) => site.ref === ref)),
        highlightSites: vi.fn((refs) => ({ found: refs.filter((ref) => sites.some((site) => site.ref === ref)), missing: refs.filter((ref) => !sites.some((site) => site.ref === ref)) })),
        clearHighlight: vi.fn(),
        selectSite: vi.fn((ref) => { if (!sites.some((site) => site.ref === ref)) return false; selected = ref; return true; }),
        ...overrides
    };
    return viewer;
};

const fakeFocus = () => {
    const undo = [];
    return {
        show: vi.fn(async () => ({})),
        onClear: vi.fn((fn) => undo.push(fn)),
        state: () => ({ cards: ['direction-sphere'] }),
        clear: () => undo.splice(0).forEach((fn) => fn())
    };
};

const tools = (activePage = 'orientation') => Object.fromEntries(
    buildViewTools({ cards: buildCopilotCards(), pages: PAGES, activePage: () => activePage, waitMs: 60, settleMs: 200 })
        .map((tool) => [tool.name, tool])
);

describe('orient_3d_view', () => {
    let unregister = [];
    afterEach(() => {
        unregister.forEach((fn) => fn());
        unregister = [];
    });
    const register = (id, viewer) => {
        unregister.push(registerViewer(id, viewer));
        return viewer;
    };

    it('is a view tool over the 3D cards, and needs something to do', async () => {
        const { orient_3d_view: tool, current_view: current } = tools();
        expect(tool.view).toBe(true);
        expect(current.view).toBe(true);
        expect(tool.parameters.properties.card.enum).toEqual(['density-cell', 'bond-cell', 'pca-density', 'pca-sites', 'direction-sphere', 'direction-sites']);
        expect(tool.cards({ card: 'pca-sites' })).toEqual(['pca-sites']);
        expect(tool.label({ card: 'direction-sphere', along: '[111]', select_site: 5 }))
            .toBe('Displacement-direction sphere (3D) · site 5 · along [111]');
        await expect(tool.run({ card: 'direction-sphere' }, {})).rejects.toThrow('nothing to do');
        await expect(tool.run({ card: 'kde-slice', along: 'c' }, {})).rejects.toThrow('"kde-slice" is not a 3D view');
        await expect(tool.run({ card: 'pca-sites', along: '[11]' }, {})).rejects.toThrow('could not read along="[11]"');
    });

    it('brings the card on screen, then looks down [111] with c up', async () => {
        const viewer = register('direction-sphere', fakeViewer());
        const focus = fakeFocus();
        const { data, checks } = await tools().orient_3d_view.run({ card: 'direction-sphere', along: '[111]' }, { focus });
        expect(focus.show).toHaveBeenCalledWith(['direction-sphere']);
        const [direction, up] = viewer.lookAlong.mock.calls[0];
        normalize([1, 1, 1]).forEach((value, i) => expect(direction[i]).toBeCloseTo(value, 12));
        normalize([-1, -1, 2]).forEach((value, i) => expect(up[i]).toBeCloseTo(value, 12));
        expect(data).toEqual({
            card: 'direction-sphere',
            page: 'Displacement Directions',
            requested: '[1 1 1]',
            view_along: '[1 1 1]',
            screen_up: '[-1 -1 2]'
        });
        expect(checks).toEqual([{ level: 'ok', text: 'Displacement-direction sphere (3D): looking down [1 1 1] with [-1 -1 2] up.' }]);
    });

    it('reads axes and principal axes', async () => {
        const viewer = register('pca-density', fakeViewer());
        await tools().orient_3d_view.run({ card: 'pca-density', along: 'b' }, { focus: fakeFocus() });
        expect(viewer.lookAlong.mock.calls[0][0]).toEqual([0, 1, 0]);
        expect(viewer.lookAlong.mock.calls[0][1]).toEqual([0, 0, 1]);
        // PC1 of this site is along y.
        const { data } = await tools().orient_3d_view.run({ card: 'pca-density', along: 'PC1' }, { focus: fakeFocus() });
        expect(viewer.lookAlong.mock.calls[1][0]).toEqual([0, 1, 0]);
        expect(data.requested).toBe('PC1');
        viewer.pcaAxes = () => null;
        await expect(tools().orient_3d_view.run({ card: 'pca-density', along: 'PC2' }, { focus: fakeFocus() }))
            .rejects.toThrow('no principal axes');
    });

    it('selects a site and waits for it to settle, then resets, centres and rings sites until the focus ends', async () => {
        let settledRef = null;
        const viewer = register('pca-sites', fakeViewer({ settled: (ref) => ref === settledRef }));
        setTimeout(() => { settledRef = 5; }, 50);
        const focus = fakeFocus();
        const { data, checks } = await tools().orient_3d_view.run({
            card: 'pca-sites', select_site: 5, reset: true, along: 'c', center_site: 5, highlight_sites: [5, 1, 9]
        }, { focus });
        expect(viewer.selectSite).toHaveBeenCalledWith(5);
        expect(viewer.reset).toHaveBeenCalled();
        expect(viewer.centerOn).toHaveBeenCalledWith(5);
        expect(viewer.highlightSites).toHaveBeenCalledWith([5, 1, 9]);
        expect(data).toMatchObject({ selected_site: 5, centred_on: 5, highlighted_sites: [5, 1], view_along: '[0 0 1]', screen_up: '[0 1 0]' });
        expect(checks).toEqual([{ level: 'warn', text: 'No site 9 on Site ellipsoids (3D).' }]);
        // The rings end with the card's focus.
        expect(focus.onClear).toHaveBeenCalledTimes(1);
        focus.clear();
        expect(viewer.clearHighlight).toHaveBeenCalledTimes(1);
    });

    it('warns when a selected site is still loading, and refuses an unknown site with the list', async () => {
        register('direction-sphere', fakeViewer({ settled: () => false }));
        const { checks } = await tools().orient_3d_view.run({ card: 'direction-sphere', select_site: 5 }, { focus: fakeFocus() });
        expect(checks).toEqual([{ level: 'warn', text: 'Site 5 was still loading; the view may change once it has drawn.' }]);
        await expect(tools().orient_3d_view.run({ card: 'direction-sphere', select_site: 42 }, { focus: fakeFocus() }))
            .rejects.toThrow('no site 42 here; sites: 1, 5');
    });

    it('notes what a view cannot do, and says when a view has not drawn', async () => {
        const viewer = register('density-cell', fakeViewer({ centerOn: undefined, highlightSites: undefined, selectSite: undefined }));
        const { checks } = await tools('structure').orient_3d_view.run({ card: 'density-cell', center_site: 1, highlight_sites: [1], select_site: 1 }, { focus: fakeFocus() });
        expect(checks.map((check) => check.level)).toEqual(['note', 'note', 'note', 'ok']);
        expect(viewer.lookAlong).not.toHaveBeenCalled();
        await expect(tools().orient_3d_view.run({ card: 'bond-cell', along: 'a' }, { focus: fakeFocus() }))
            .rejects.toThrow('the Folded cell with bonds (3D) view has not drawn yet');
    });
});

describe('current_view', () => {
    it('reports the page on screen, its cards, the highlight and its 3D views', async () => {
        const unregister = [
            registerViewer('direction-sphere', fakeViewer()),
            registerViewer('direction-sites', fakeViewer({ ready: () => false }))
        ];
        try {
            const { data } = await tools('orientation').current_view.run({}, { focus: fakeFocus() });
            expect(data).toEqual({
                page: 'Displacement Directions',
                cards: ['direction-sphere', 'direction-views', 'direction-sites'],
                highlighted: ['direction-sphere'],
                views_3d: [
                    { card: 'direction-sphere', view_along: '[1 1 1]', site: 1 },
                    { card: 'direction-sites', drawn: false }
                ]
            });
        } finally {
            unregister.forEach((fn) => fn());
        }
    });
});
