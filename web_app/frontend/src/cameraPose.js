// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Camera moves for the Copilot's 3D view tools (viewer3d.js): read the pose
// of an OrbitControls camera, and turn it to a new pose — a short slerp of the
// camera's orientation (the model visibly rotates into place), or at once
// under prefers-reduced-motion. The final pose is set exactly, whatever the
// frame timing. OrbitControls' leftover damped velocity is flushed first, as
// PcaKdePage's placeMainCamera does, so no drag momentum bends the result.

import * as THREE from 'three';

export const TURN_MS = 650;

const reducedMotion = () => typeof window !== 'undefined'
    && typeof window.matchMedia === 'function'
    && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

/** `{ target, direction, distance, up }` of a camera orbiting controls.target. */
export const readPose = (camera, controls) => {
    const target = controls.target.clone();
    const offset = camera.position.clone().sub(target);
    const distance = offset.length() || 1;
    return {
        target: target.toArray(),
        direction: offset.normalize().toArray(),
        distance,
        up: camera.up.clone().normalize().toArray()
    };
};

// The camera orientation whose +z (toward the viewer) is `direction` and
// whose +y is `up`, as a quaternion.
const orientation = (direction, up) => {
    const z = new THREE.Vector3(...direction).normalize();
    const y = new THREE.Vector3(...up);
    y.addScaledVector(z, -y.dot(z));
    if (y.lengthSq() < 1e-12) {
        // `up` along the line of sight: any perpendicular will do.
        y.set(...(Math.abs(z.y) < 0.9 ? [0, 1, 0] : [1, 0, 0]));
        y.addScaledVector(z, -y.dot(z));
    }
    y.normalize();
    const x = new THREE.Vector3().crossVectors(y, z).normalize();
    return new THREE.Quaternion().setFromRotationMatrix(new THREE.Matrix4().makeBasis(x, y, z));
};

const flushDamping = (controls) => {
    const damping = controls.enableDamping;
    controls.enableDamping = false;
    controls.update();
    controls.enableDamping = damping;
};

/** Set a pose exactly: target, position = target + direction·distance, up. */
export const applyPose = (camera, controls, { target, direction, distance, up }) => {
    flushDamping(controls);
    const t = new THREE.Vector3(...target);
    const d = new THREE.Vector3(...direction).normalize();
    controls.target.copy(t);
    camera.up.set(...up).normalize();
    camera.position.copy(t).addScaledVector(d, distance);
    camera.lookAt(t);
    camera.updateProjectionMatrix();
    const damping = controls.enableDamping;
    controls.enableDamping = false;
    controls.update();
    controls.enableDamping = damping;
};

/**
 * Turn the camera to `pose` over `ms` (a slerp of orientation, a lerp of
 * target and distance); returns a cancel function. Immediate when `ms` is 0
 * or the user prefers reduced motion. `onDone` runs once the pose is set.
 */
const raf = typeof requestAnimationFrame === 'function' ? (fn) => requestAnimationFrame(fn) : null;
const caf = typeof cancelAnimationFrame === 'function' ? (id) => cancelAnimationFrame(id) : () => {};

export const turnTo = (camera, controls, pose, { ms = TURN_MS, onDone = null, now = () => performance.now(), frame = raf, cancelFrame = caf } = {}) => {
    const goal = {
        target: pose.target,
        direction: pose.direction,
        distance: pose.distance,
        up: pose.up
    };
    if (!ms || reducedMotion() || typeof frame !== 'function') {
        applyPose(camera, controls, goal);
        onDone?.();
        return () => {};
    }
    flushDamping(controls);
    const start = readPose(camera, controls);
    const q0 = orientation(start.direction, start.up);
    const q1 = orientation(goal.direction, goal.up);
    const t0 = new THREE.Vector3(...start.target);
    const t1 = new THREE.Vector3(...goal.target);
    const began = now();
    let id = null;
    let cancelled = false;
    const step = () => {
        if (cancelled) return;
        const s = Math.min(1, (now() - began) / ms);
        if (s >= 1) {
            applyPose(camera, controls, goal);
            onDone?.();
            return;
        }
        // Ease in-out, so the turn starts and lands softly.
        const e = s < 0.5 ? 2 * s * s : 1 - ((-2 * s + 2) ** 2) / 2;
        const q = q0.clone().slerp(q1, e);
        const direction = new THREE.Vector3(0, 0, 1).applyQuaternion(q);
        const up = new THREE.Vector3(0, 1, 0).applyQuaternion(q);
        applyPose(camera, controls, {
            target: t0.clone().lerp(t1, e).toArray(),
            direction: direction.toArray(),
            distance: start.distance + (goal.distance - start.distance) * e,
            up: up.toArray()
        });
        id = frame(step);
    };
    id = frame(step);
    return () => {
        cancelled = true;
        if (id != null) cancelFrame(id);
    };
};
