// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// A panel's registration in the Copilot's view registry: the controller
// appears while the panel is mounted, its orbit functions turn the panel's
// camera (at once under prefers-reduced-motion), a panel's own function
// overrides or extends them with the orbit as its last argument, it always
// reads the panel's latest props, and it leaves on unmount.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import React, { act } from 'react';
import { createRoot } from 'react-dom/client';
import * as THREE from 'three';
import useCopilotViewer from '../useCopilotViewer';
import { getViewer } from '../viewer3d';

const CUBIC = [[5, 0, 0], [0, 5, 0], [0, 0, 5]];

const Panel = ({ id, handle, site, onSelect }) => {
    useCopilotViewer(id, {
        orbit: { getHandle: () => handle, getCell: () => CUBIC, onReset: handle.reset },
        ready: () => site != null,
        describe: (orbit) => ({ ...orbit.describe(), site }),
        selectSite: (ref) => { onSelect(ref); return true; }
    });
    return null;
};

describe('useCopilotViewer', () => {
    let container;
    let root;
    beforeEach(() => {
        globalThis.IS_REACT_ACT_ENVIRONMENT = true;
        vi.stubGlobal('matchMedia', (query) => ({ matches: query.includes('reduce'), addEventListener() {}, removeEventListener() {} }));
        window.matchMedia = globalThis.matchMedia;
        container = document.createElement('div');
        document.body.appendChild(container);
        root = createRoot(container);
    });
    afterEach(() => {
        act(() => root.unmount());
        container.remove();
        vi.unstubAllGlobals();
        delete window.matchMedia;
    });

    it('registers while mounted and turns the camera through the orbit', () => {
        const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 100);
        camera.position.set(0, 0, 9);
        const handle = { camera, controls: { target: new THREE.Vector3(), enableDamping: true, update: () => {} }, reset: vi.fn() };
        const onSelect = vi.fn();
        act(() => root.render(<Panel id="test-panel" handle={handle} site={null} onSelect={onSelect} />));
        const viewer = getViewer('test-panel');
        expect(viewer.ready()).toBe(false);
        expect(viewer.has('lookAlong')).toBe(true);
        expect(viewer.has('centerOn')).toBe(false);
        expect(viewer.cell()).toEqual(CUBIC);

        act(() => root.render(<Panel id="test-panel" handle={handle} site={3} onSelect={onSelect} />));
        expect(viewer.ready()).toBe(true);
        expect(viewer.describe()).toEqual({ view_along: '[0 0 1]', screen_up: '[0 1 0]', site: 3 });

        expect(viewer.lookAlong([1, 0, 0], [0, 0, 1])).toBe(true);
        [9, 0, 0].forEach((value, i) => expect(camera.position.toArray()[i]).toBeCloseTo(value, 9));
        expect(viewer.describe()).toEqual({ view_along: '[1 0 0]', screen_up: '[0 0 1]', site: 3 });
        expect(viewer.centerAt([1, 1, 1])).toBe(true);
        expect(handle.controls.target.toArray()).toEqual([1, 1, 1]);
        viewer.reset();
        expect(handle.reset).toHaveBeenCalled();
        expect(viewer.selectSite(7)).toBe(true);
        expect(onSelect).toHaveBeenCalledWith(7);

        act(() => root.render(<div />));
        expect(getViewer('test-panel')).toBeNull();
    });
});
