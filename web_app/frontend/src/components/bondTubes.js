// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The folded cell's bond tubes (FoldedCellPanel): one instanced mesh per
// bond set, as before, but drawn like a ball-and-stick figure's sticks.
//
// - Round: a lit Phong material with a soft highlight, smooth normals and
//   rounded ends (a capsule whose caps keep their shape at any bond length:
//   the cap vertices are pushed out along the bond in the shader, by the
//   tube's radius, after the instance matrix stretches the body).
// - Edges darken where the tube turns away from the eye, so crossing bonds
//   stay apart without a theme-dependent fog.
// - Each end stops BOND_END_GAP short of its site, at the edge of the site's
//   thermal cloud, so the cloud stays visible.
// - A per-instance opacity (and radius) lets the angle highlight draw each
//   bond as strongly as its share.
//
// On a software renderer (SwiftShader, llvmpipe) the tubes take fewer
// segments and the canvas skips multisampling and draws at 1× pixel ratio.

import * as THREE from 'three';

// Å each tube stops short of its site centre (at most a quarter of the
// bond), about a thermal cloud's radius.
export const BOND_END_GAP = 0.22;
const MAX_GAP_FRACTION = 0.25;
// How much the edges darken (1: not at all).
const RIM_SHADE = 0.55;

let qualityCache = null;

/**
 * Whether this browser renders WebGL in software, read once from a probe
 * context; `{ software, radialSegments, capRings, antialias, maxPixelRatio }`.
 */
export const renderQuality = () => {
    if (qualityCache) return qualityCache;
    let software = false;
    try {
        const gl = document.createElement('canvas').getContext('webgl2') || document.createElement('canvas').getContext('webgl');
        const info = gl?.getExtension('WEBGL_debug_renderer_info');
        const name = String(gl && info ? gl.getParameter(info.UNMASKED_RENDERER_WEBGL) : gl?.getParameter(gl.RENDERER) ?? '');
        software = /swiftshader|llvmpipe|softpipe|software/i.test(name);
        gl?.getExtension('WEBGL_lose_context')?.loseContext();
    } catch {
        // No probe: assume a GPU.
    }
    qualityCache = software
        ? { software, radialSegments: 8, capRings: 2, antialias: false, maxPixelRatio: 1 }
        : { software, radialSegments: 16, capRings: 3, antialias: true, maxPixelRatio: 2 };
    return qualityCache;
};

/**
 * A unit capsule: radius 1, body y ∈ [−0.5, 0.5], each cap's vertices on the
 * body's end ring with `capOffset` the distance (in radii) the shader pushes
 * them out along the bond. Rings run from the bottom pole to the top one.
 */
export const capsuleGeometry = (radialSegments = 16, capRings = 3) => {
    const rings = [];
    for (let k = capRings; k >= 0; k -= 1) rings.push({ y: -0.5, phi: -(k / capRings) * Math.PI / 2 });
    for (let k = 0; k <= capRings; k += 1) rings.push({ y: 0.5, phi: (k / capRings) * Math.PI / 2 });
    const n = radialSegments;
    const positions = new Float32Array(rings.length * n * 3);
    const normals = new Float32Array(rings.length * n * 3);
    const offsets = new Float32Array(rings.length * n);
    rings.forEach(({ y, phi }, r) => {
        const ring = Math.cos(phi);
        for (let i = 0; i < n; i += 1) {
            const theta = (i / n) * Math.PI * 2;
            const v = r * n + i;
            positions.set([Math.cos(theta) * ring, y, Math.sin(theta) * ring], v * 3);
            normals.set([Math.cos(theta) * ring, Math.sin(phi), Math.sin(theta) * ring], v * 3);
            offsets[v] = Math.sin(phi);
        }
    });
    const index = [];
    for (let r = 0; r + 1 < rings.length; r += 1) {
        for (let i = 0; i < n; i += 1) {
            const a = r * n + i;
            const b = r * n + ((i + 1) % n);
            const c = a + n;
            const d = b + n;
            index.push(a, c, b, b, c, d);
        }
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('normal', new THREE.BufferAttribute(normals, 3));
    geometry.setAttribute('capOffset', new THREE.BufferAttribute(offsets, 1));
    geometry.setIndex(index);
    return geometry;
};

/**
 * The tubes' material: Phong with a soft specular, rounded caps, darkened
 * edges and a per-instance opacity (attribute `instanceOpacity`).
 */
export const bondMaterial = ({ color, emissiveIntensity = 0, opacity = 1, transparent = false } = {}) => {
    const material = new THREE.MeshPhongMaterial({
        color,
        emissive: color,
        emissiveIntensity,
        specular: 0x3a3a3a,
        shininess: 48,
        opacity,
        transparent
    });
    material.onBeforeCompile = (shader) => {
        shader.vertexShader = shader.vertexShader
            .replace('#include <common>', `#include <common>
attribute float capOffset;
attribute float instanceOpacity;
varying float vBondOpacity;`)
            // The instance matrix is rotation × scale (radius, length,
            // radius): its rotation alone turns the normals, caps included.
            .replace('#include <defaultnormal_vertex>', `vec3 transformedNormal = objectNormal;
#ifdef USE_INSTANCING
    mat3 bondFrame = mat3(instanceMatrix);
    transformedNormal = mat3(normalize(bondFrame[0]), normalize(bondFrame[1]), normalize(bondFrame[2])) * transformedNormal;
#endif
transformedNormal = normalMatrix * transformedNormal;`)
            .replace('#include <project_vertex>', `vec4 mvPosition = vec4(transformed, 1.0);
#ifdef USE_INSTANCING
    mvPosition = instanceMatrix * mvPosition;
    mvPosition.xyz += normalize(instanceMatrix[1].xyz) * capOffset * length(instanceMatrix[0].xyz);
    vBondOpacity = instanceOpacity;
#else
    vBondOpacity = 1.0;
#endif
mvPosition = modelViewMatrix * mvPosition;
gl_Position = projectionMatrix * mvPosition;`);
        shader.fragmentShader = shader.fragmentShader
            .replace('#include <common>', `#include <common>
varying float vBondOpacity;`)
            .replace('#include <alphamap_fragment>', `#include <alphamap_fragment>
diffuseColor.a *= vBondOpacity;`)
            .replace('#include <opaque_fragment>', `outgoingLight *= mix(${RIM_SHADE.toFixed(2)}, 1.0,
    smoothstep(0.0, 0.75, abs(dot(normalize(normal), normalize(vViewPosition)))));
#include <opaque_fragment>`);
    };
    material.customProgramCacheKey = () => 'bond-tube-1';
    return material;
};

/**
 * One instanced capsule per segment (six numbers each: from xyz, to xyz),
 * each end `gap` short of its point (at most a quarter of the bond, caps
 * included). `scales` (per segment, × radius) and `opacities` (per
 * segment) default to 1.
 */
export const tubeMesh = (segments, radius, material, {
    gap = BOND_END_GAP, scales = null, opacities = null, quality = renderQuality()
} = {}) => {
    const count = segments.length / 6;
    const geometry = capsuleGeometry(quality.radialSegments, quality.capRings);
    const opacity = new Float32Array(count).fill(1);
    if (opacities) opacity.set(opacities.subarray ? opacities.subarray(0, count) : opacities.slice(0, count));
    geometry.setAttribute('instanceOpacity', new THREE.InstancedBufferAttribute(opacity, 1));
    const mesh = new THREE.InstancedMesh(geometry, material, count);
    const up = new THREE.Vector3(0, 1, 0);
    const start = new THREE.Vector3();
    const end = new THREE.Vector3();
    const direction = new THREE.Vector3();
    const rotation = new THREE.Quaternion();
    const scale = new THREE.Vector3();
    const matrix = new THREE.Matrix4();
    for (let i = 0; i < count; i += 1) {
        start.fromArray(segments, i * 6);
        end.fromArray(segments, i * 6 + 3);
        direction.subVectors(end, start);
        const length = direction.length();
        const r = radius * (scales ? scales[i] : 1);
        const trim = Math.min(gap, MAX_GAP_FRACTION * length);
        // The body between the caps; never zero, so its axis stays defined.
        const body = Math.max(length - 2 * trim - 2 * r, 1e-6 * (length || 1));
        rotation.setFromUnitVectors(up, length > 0 ? direction.divideScalar(length) : up);
        scale.set(r, body, r);
        matrix.compose(start.add(end).multiplyScalar(0.5), rotation, scale);
        mesh.setMatrixAt(i, matrix);
    }
    mesh.instanceMatrix.needsUpdate = true;
    mesh.computeBoundingSphere();
    return mesh;
};
