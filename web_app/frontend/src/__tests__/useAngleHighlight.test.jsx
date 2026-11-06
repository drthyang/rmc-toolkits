// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */
/* global process */

// The page's hover highlight hook on the demo run: after a Compute the
// worker's triplets arrive, a hovered angle picks its bin's triplets, leaving
// clears it, and an answer whose counts are not the histogram's (another
// configuration or request) is never shown. A stale request is stopped.

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { afterEach, beforeEach, describe, expect, it } from 'vitest';
import { act, useEffect } from 'react';
import { createRoot } from 'react-dom/client';
import useAngleHighlight from '../useAngleHighlight';
import { foldedCellNetwork } from '../bondNetwork';
import { handleAngleTriplets } from '../workers/angleTripletsWorker';
import { siteDisplacementsFromRmc6f } from '../workers/pcaKde';
import { bondAngleSummary, tripletOptionsFromRequest } from '../workers/triplets';
import { BOND_COLORS } from '../plotPalette';

const text = readFileSync(resolve(process.cwd(), 'public/demo/GTS_250K.rmc6f'), 'utf8');
const parsed = siteDisplacementsFromRmc6f(text);
const sites = {
    latticeVectors: parsed.latticeVectors,
    supercell: parsed.supercell,
    sites: parsed.sites.map((site) => ({ referenceNumber: site.referenceNumber, element: site.element, siteFractional: site.siteFractional }))
};
const request = { end1: 'Se', apex: 'Ta', end2: 'Se', r12Min: 2.3, r12Max: 2.9, binWidth: 1 };
const result = bondAngleSummary(parsed.atomList.fractional, parsed.atomList.elements, parsed.latticeVectors, tripletOptionsFromRequest(request));
const bondSets = [{ elements: ['Se', 'Ta'], window: result.bond12, color: BOND_COLORS.ab }];
const network = foldedCellNetwork(null, sites, bondSets);

// A worker that answers in-process, on demand.
const fakeWorkers = [];
const createWorker = () => {
    const listeners = new Set();
    const worker = {
        posted: [],
        terminated: false,
        addEventListener: (type, fn) => listeners.add(fn),
        removeEventListener: (type, fn) => listeners.delete(fn),
        postMessage: (data) => worker.posted.push(data),
        terminate: () => { worker.terminated = true; },
        answer: (tamper = (value) => value) => {
            const data = worker.posted.shift();
            const message = { data: { id: data.id, result: tamper(handleAngleTriplets(data)) } };
            [...listeners].forEach((fn) => fn(message));
        }
    };
    fakeWorkers.push(worker);
    return worker;
};

const probe = { current: null };
const Probe = (props) => {
    const value = useAngleHighlight({ createWorker, sites, network, ...props });
    useEffect(() => { probe.current = value; });
    return null;
};

describe('useAngleHighlight', () => {
    let container;
    let root;
    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        fakeWorkers.length = 0;
        container = document.createElement('div');
        root = createRoot(container);
    });
    afterEach(() => act(() => root.unmount()));

    const render = (props) => act(() => root.render(<Probe rmc6fText={text} result={result} request={request} {...props} />));

    it('highlights the hovered bin with the histogram\'s count, and clears on leave', () => {
        render();
        expect(fakeWorkers).toHaveLength(1);
        act(() => fakeWorkers[0].answer());
        expect(probe.current.ready).toBe(true);
        const busiest = result.counts.indexOf(Math.max(...result.counts));
        act(() => probe.current.onHoverX(busiest + 0.5));
        expect(probe.current.highlight.total).toBe(result.counts[busiest]);
        expect(probe.current.highlight.shown).toBe(result.counts[busiest]);
        expect(probe.current.highlight.ab.length).toBeGreaterThan(0);
        act(() => probe.current.onHoverX(null));
        expect(probe.current.highlight).toBeNull();
    });

    it('highlights a picked range over the hovered bin, with the range\'s count', () => {
        render({ range: [85, 95] });
        act(() => fakeWorkers[0].answer());
        act(() => probe.current.onHoverX(30.5));
        const inRange = result.counts.slice(85, 95).reduce((a, b) => a + b, 0);
        expect(probe.current.highlight.total).toBe(inRange);
        expect(probe.current.highlight.range).toBe(true);
        render({ range: null });
        expect(probe.current.highlight.total).toBe(result.counts[30]);
    });

    it('never shows an answer whose counts are not the histogram\'s', () => {
        render();
        act(() => fakeWorkers[0].answer((data) => ({ ...data, bins: data.bins.map((bin, i) => (i === 90 ? { ...bin, total: bin.total + 1 } : bin)) })));
        act(() => probe.current.onHoverX(90.5));
        expect(probe.current.ready).toBe(false);
        expect(probe.current.highlight).toBeNull();
    });

    it('stops a request a newer one made stale', () => {
        render();
        const first = fakeWorkers[0];
        const other = { ...result };
        render({ result: other });
        expect(first.terminated).toBe(true);
        expect(fakeWorkers).toHaveLength(2);
        act(() => fakeWorkers[1].answer());
        expect(probe.current.ready).toBe(true);
    });

    it('asks nothing without the structure text (a server folder)', () => {
        render({ rmc6fText: null });
        expect(fakeWorkers).toHaveLength(0);
        act(() => probe.current.onHoverX(90.5));
        expect(probe.current.highlight).toBeNull();
        expect(probe.current.available).toBe(false);
    });
});
