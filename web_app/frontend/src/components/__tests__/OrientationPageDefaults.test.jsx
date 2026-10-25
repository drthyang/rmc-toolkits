// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The Displacement Directions page opens at ν = 10 with the map smoothed 4×
// (the user's choice, 2026-10-24); Auto, the engine's resolution (the tiling
// analyze_displacement_directions reads), shows the ν the engine chose; a
// picked resolution stays when the site changes; and every readout comes from
// the raw histogram (smoothing 0) while the sphere draws the map at the
// page's smoothing (plan item 2.5), so the smoothing never moves a readout.
// The engine is a recording fake; the Three.js panels are stubbed.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { act } from 'react';
import { createRoot } from 'react-dom/client';

const state = vi.hoisted(() => ({ requests: [], view: null }));

// The engine's choice for every site here: ν = 2 (42 cells).
const RECOMMENDED = 2;
// A payload that says which request it answers.
const payload = (params) => {
    const frequency = params.frequency ?? RECOMMENDED;
    return {
        referenceNumber: params.referenceNumber,
        frequency,
        cellCount: 10 * frequency * frequency + 2,
        recommendedFrequency: RECOMMENDED,
        smoothing: params.smoothing,
        weight: params.weight,
        // Raw and smoothed maps disagree, as the engines' peak cells can.
        peakSignificance: params.smoothing > 0 ? 1.5 : 4.5,
        peakDirection: params.smoothing > 0 ? [0, 0, -1] : [0, 0, 1]
    };
};
const requestPca = async (kind, params) => {
    state.requests.push({ kind, params });
    return payload(params);
};

const SITES = {
    elements: ['K', 'Nb', 'O'],
    sites: [
        { referenceNumber: 1, element: 'K', count: 512, uIso: 0.004 },
        { referenceNumber: 2, element: 'Nb', count: 512, uIso: 0.006 }
    ]
};

vi.mock('../../useSiteCloud', async () => {
    const { useState } = await import('react');
    // The hook's contract, with the picked site as its state.
    const useSiteCloudFake = () => {
        const [selectedRef, setSelectedRef] = useState(1);
        return {
            sites: SITES,
            sitesError: null,
            loadingSites: false,
            selectedRef,
            setSelectedRef,
            selectedEllipsoid: SITES.sites.find((site) => site.referenceNumber === selectedRef) ?? null,
            requestPca,
            localFile: null,
            rmc6fText: null,
            ready: true,
            unitCell: null,
            datasetKey: 'run'
        };
    };
    return { default: useSiteCloudFake };
});

// The sphere is Three.js: record what the page hands it instead.
vi.mock('../OrientationView', () => ({
    default: (props) => {
        state.view = props;
        return null;
    }
}));
vi.mock('../SiteStructurePanel', () => ({ default: () => null }));

const { default: OrientationPage } = await import('../OrientationPage');

describe('the Displacement Directions page opens at ν = 10, smoothed 4× (2026-10-24)', () => {
    let container;
    let root;

    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        state.requests = [];
        state.view = null;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });

    afterEach(() => {
        act(() => root.unmount());
        container.remove();
    });

    const render = async (props = {}) => {
        await act(async () => {
            root.render(<OrientationPage directory="." localRun={null} {...props} />);
        });
    };
    const select = (label) => container.querySelector(`select[aria-label="${label}"]`);
    const change = async (element, value) => {
        await act(async () => {
            element.value = String(value);
            element.dispatchEvent(new Event('change', { bubbles: true }));
        });
    };
    // React tracks a range input's value: set it through the native setter.
    const slide = async (label, value) => {
        const input = container.querySelector(`input[aria-label="${label}"]`);
        await act(async () => {
            Object.getOwnPropertyDescriptor(window.HTMLInputElement.prototype, 'value').set.call(input, String(value));
            input.dispatchEvent(new Event('input', { bubbles: true }));
        });
    };
    const resolution = () => select('Sphere resolution');
    const selectedText = (element) => element.options[element.selectedIndex].textContent;
    const orientationRequests = () => state.requests.filter((request) => request.kind === 'orientation').map((request) => request.params);
    // The latest readout request (every readout request is at smoothing 0).
    const lastReadout = () => orientationRequests().filter((params) => params.smoothing === 0).at(-1);

    it('opens at ν = 10 with the map smoothed 4×, and reads its readouts from the raw histogram', async () => {
        await render();
        expect(resolution().value).toBe('10');
        expect(selectedText(resolution())).toBe('ν=10 (1002)');
        // Auto still names the engine's choice for the site on screen.
        expect([...resolution().options].find((option) => option.value === 'auto').textContent).toBe('Auto · ν=2 (42)');
        // The map opens smoothed 4× (four neighbour passes).
        expect(container.querySelector('input[aria-label="Neighbour smoothing passes"]').value).toBe('4');
        expect(container.textContent).toContain('4×');

        // Two requests landing together: the readout at smoothing 0 and the
        // drawn map, the same request smoothed 4×.
        const readoutRequest = { referenceNumber: 1, frequency: 10, weight: 'count', frame: 'cartesian', minAmplitudeQuantile: 0, smoothing: 0 };
        expect(orientationRequests()).toEqual([
            expect.objectContaining(readoutRequest),
            expect.objectContaining({ ...readoutRequest, smoothing: 4 })
        ]);
        expect(state.view.readout).toMatchObject({ smoothing: 0, peakSignificance: 4.5, frequency: 10 });
        expect(state.view.result).toMatchObject({ smoothing: 4, peakSignificance: 1.5, frequency: 10 });
        expect(state.view.loadedRef).toBe(1);
    });

    it('the smoothing redraws the map only: the readout is neither asked again nor changed', async () => {
        await render();
        const readout = state.view.readout;
        const [readoutRequest] = orientationRequests();
        await slide('Neighbour smoothing passes', 5);
        // One request: the map, smoothed; the readout is the raw one it was.
        expect(orientationRequests().slice(2)).toEqual([{ ...readoutRequest, smoothing: 5 }]);
        expect(state.view.readout).toBe(readout);
        expect(state.view.result).toMatchObject({ smoothing: 5, peakSignificance: 1.5 });
        // At 0 the readout is the map: nothing to ask.
        await slide('Neighbour smoothing passes', 0);
        expect(orientationRequests()).toHaveLength(3);
        expect(state.view.result).toBe(readout);
        expect(state.view.readout).toBe(readout);
    });

    it('a picked resolution overrides the default and stays when the site changes; Auto names the engine\'s ν meanwhile', async () => {
        await render();
        await change(resolution(), 6);
        expect(resolution().value).toBe('6');
        expect(lastReadout()).toMatchObject({ referenceNumber: 1, frequency: 6, smoothing: 0 });
        // The Auto choice still says what Auto would pick for this site.
        expect([...resolution().options].find((option) => option.value === 'auto').textContent).toBe('Auto · ν=2 (42)');

        // Another site: the pick is page state and stays.
        await change(select('Site'), 2);
        expect(resolution().value).toBe('6');
        expect(lastReadout()).toMatchObject({ referenceNumber: 2, frequency: 6, smoothing: 0 });
        expect(state.view.loadedRef).toBe(2);

        // Back to Auto: the engine chooses again.
        await change(resolution(), 'auto');
        expect(lastReadout()).toMatchObject({ referenceNumber: 2, frequency: null, smoothing: 0 });
        expect(state.view.readout.frequency).toBe(RECOMMENDED);
    });

    it('a Copilot step\'s Show in sets its site and resolution', async () => {
        await render();
        await change(resolution(), 6);
        await render({ request: { id: 1, page: 'orientation', params: { referenceNumber: 2, frequency: 'auto' } } });
        expect(select('Site').value).toBe('2');
        expect(resolution().value).toBe('auto');
        expect(lastReadout()).toMatchObject({ referenceNumber: 2, frequency: null, smoothing: 0 });
    });
});
