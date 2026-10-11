// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Camera turns for the Copilot's 3D view tools: a pose is set exactly, a turn
// passes through in-between orientations and lands exactly on the pose, a
// half-turn (looking from the opposite side) is well defined, and leftover
// OrbitControls momentum is flushed first.

import { describe, expect, it, vi } from 'vitest';
import * as THREE from 'three';
import { applyPose, readPose, turnTo } from '../cameraPose';

// OrbitControls needs a DOM element; its part in a pose is the target, the
// damping flag and update().
const rig = () => {
    const camera = new THREE.PerspectiveCamera(45, 1, 0.01, 100);
    camera.position.set(3, 2, 4);
    const controls = { target: new THREE.Vector3(0, 0, 0), enableDamping: true, update: vi.fn() };
    camera.lookAt(controls.target);
    return { camera, controls };
};

const close = (actual, expected, digits = 9) => actual.forEach((value, i) => expect(value).toBeCloseTo(expected[i], digits));

// A scripted clock and frame queue, so a turn can be stepped through.
const clock = () => {
    let time = 0;
    const queue = [];
    return {
        now: () => time,
        frame: (fn) => { queue.push(fn); return queue.length; },
        cancelFrame: vi.fn(),
        advance(ms) {
            time += ms;
            const due = queue.splice(0);
            due.forEach((fn) => fn());
        }
    };
};

describe('camera poses', () => {
    it('reads and sets a pose exactly, flushing the damped velocity first', () => {
        const { camera, controls } = rig();
        const pose = readPose(camera, controls);
        expect(pose.distance).toBeCloseTo(Math.hypot(3, 2, 4), 12);
        close(pose.direction, [3, 2, 4].map((v) => v / Math.hypot(3, 2, 4)));

        applyPose(camera, controls, { target: [1, 1, 1], direction: [0, 0, 2], distance: 5, up: [0, 1, 0] });
        close(camera.position.toArray(), [1, 1, 6]);
        close(controls.target.toArray(), [1, 1, 1]);
        close(camera.up.toArray(), [0, 1, 0]);
        // The camera looks at the target: its −z axis points there.
        close(new THREE.Vector3(0, 0, -1).applyQuaternion(camera.quaternion).toArray(), [0, 0, -1]);
        expect(controls.update).toHaveBeenCalled();
        expect(controls.enableDamping).toBe(true);
    });

    it('turns through in-between poses and lands exactly, keeping the distance', () => {
        const { camera, controls } = rig();
        const goal = { target: [0, 0, 0], direction: [1, 0, 0], distance: 6, up: [0, 0, 1] };
        const time = clock();
        const done = vi.fn();
        turnTo(camera, controls, goal, { ms: 600, now: time.now, frame: time.frame, cancelFrame: time.cancelFrame, onDone: done });
        time.advance(300);
        const midway = readPose(camera, controls).direction;
        // Halfway, neither where it started nor where it goes.
        expect(new THREE.Vector3(...midway).angleTo(new THREE.Vector3(1, 0, 0))).toBeGreaterThan(0.05);
        expect(done).not.toHaveBeenCalled();
        time.advance(400);
        expect(done).toHaveBeenCalledTimes(1);
        close(camera.position.toArray(), [6, 0, 0]);
        close(camera.up.toArray(), [0, 0, 1]);
    });

    it('can be cancelled, and turns at once when ms is 0', () => {
        const { camera, controls } = rig();
        const time = clock();
        const cancel = turnTo(camera, controls, { target: [0, 0, 0], direction: [0, 1, 0], distance: 2, up: [0, 0, 1] },
            { ms: 600, now: time.now, frame: time.frame, cancelFrame: time.cancelFrame });
        cancel();
        time.advance(700);
        // Cancelled before its first frame: the camera never moved.
        close(camera.position.toArray(), [3, 2, 4]);

        // Straight down the opposite side (a half-turn) is still well defined.
        turnTo(camera, controls, { target: [0, 0, 0], direction: [-3, -2, -4], distance: 2, up: [0, 1, 0] }, { ms: 0 });
        close(camera.position.toArray(), [-3, -2, -4].map((v) => (2 * v) / Math.hypot(3, 2, 4)));
        // An up along the line of sight falls back to a perpendicular.
        turnTo(camera, controls, { target: [0, 0, 0], direction: [0, 1, 0], distance: 1, up: [0, 1, 0] }, { ms: 0 });
        close(camera.position.toArray(), [0, 1, 0]);
    });
});
