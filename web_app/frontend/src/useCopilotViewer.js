// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A Three.js panel joins the Copilot's view registry (viewer3d.js) while it
// is mounted, under its card id. `spec` holds the panel's own functions and,
// as `spec.orbit`, how to reach its OrbitControls camera:
//
//   orbit: {
//     getHandle: () => ({ camera, controls }),   // the live pair
//     getCell: () => [a, b, c] | null,          // unit-cell rows (Cartesian)
//     toScene / fromScene: (v) => v,            // Cartesian ↔ the scene's frame
//     onReset: () => …                          // the panel's own Reset view
//   }
//
// From the orbit the hook builds the common functions — describe the view,
// look along a direction (an animated turn), centre a point, reset — which a
// panel may override or extend: every spec function receives the orbit as
// its last argument. The spec is read through a ref (updated after each
// render), so its functions always see the panel's latest state, and the
// controller registers once per id.

import { useEffect, useLayoutEffect, useRef } from 'react';
import { readPose, turnTo } from './cameraPose';
import { describeDirection, posePlan, registerViewer } from './viewer3d';

const METHODS = [
    'describe', 'lookAlong', 'reset', 'centerOn', 'centerAt', 'highlightSites', 'clearHighlight',
    'selectSite', 'sites', 'cell', 'pcaAxes', 'settled'
];

/**
 * The common controller functions over an OrbitControls scene handle (see
 * the `orbit` config above). Turns are animated (cameraPose.turnTo); a new
 * one cancels the last.
 */
export const orbitViewer = ({ getHandle, getCell, toScene = (v) => v, fromScene = (v) => v, onReset }) => {
    let cancel = null;
    const stop = () => {
        cancel?.();
        cancel = null;
    };
    const live = () => {
        const handle = getHandle?.();
        return handle?.camera && handle?.controls ? handle : null;
    };
    return {
        describe: () => {
            const handle = live();
            if (!handle) return null;
            const pose = readPose(handle.camera, handle.controls);
            // The line of sight runs from the target toward the camera.
            return describeDirection(fromScene(pose.direction), fromScene(pose.up), getCell?.() ?? null);
        },
        lookAlong: (direction, up) => {
            const handle = live();
            if (!handle) return false;
            stop();
            const current = readPose(handle.camera, handle.controls);
            const plan = posePlan({ target: current.target, direction: toScene(direction), distance: current.distance, up: toScene(up) });
            if (!plan) return false;
            cancel = turnTo(handle.camera, handle.controls, { ...plan, distance: current.distance }, { onDone: () => { cancel = null; } });
            return true;
        },
        // Move the orbit centre to `point` (scene frame), keeping the view.
        centerAt: (point) => {
            const handle = live();
            if (!handle) return false;
            stop();
            const current = readPose(handle.camera, handle.controls);
            cancel = turnTo(handle.camera, handle.controls, { ...current, target: point }, { onDone: () => { cancel = null; } });
            return true;
        },
        reset: () => {
            stop();
            onReset?.();
            return true;
        },
        cell: () => getCell?.() ?? null
    };
};

export default function useCopilotViewer(id, spec) {
    const specRef = useRef(spec);
    useLayoutEffect(() => {
        specRef.current = spec;
    });
    useEffect(() => {
        if (!id) return undefined;
        const config = () => specRef.current?.orbit || null;
        const orbit = orbitViewer({
            getHandle: () => config()?.getHandle?.() ?? null,
            getCell: () => config()?.getCell?.() ?? null,
            toScene: (v) => (config()?.toScene ? config().toScene(v) : v),
            fromScene: (v) => (config()?.fromScene ? config().fromScene(v) : v),
            onReset: () => config()?.onReset?.()
        });
        const own = (name) => (typeof specRef.current?.[name] === 'function' ? specRef.current[name] : null);
        const fromOrbit = (name) => (config() && typeof orbit[name] === 'function' ? orbit[name] : null);
        const controller = {
            id,
            ready: () => Boolean(own('ready')?.(orbit)),
            has: (name) => Boolean(own(name) || fromOrbit(name))
        };
        METHODS.forEach((name) => {
            controller[name] = (...args) => {
                const fn = own(name);
                if (fn) return fn(...args, orbit);
                const shared = fromOrbit(name);
                return shared ? shared(...args) : undefined;
            };
        });
        return registerViewer(id, controller);
    }, [id]);
}
