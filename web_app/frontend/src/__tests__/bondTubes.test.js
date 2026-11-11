// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The folded cell's bond tubes (components/bondTubes.js): a capsule whose
// rounded ends stop BOND_END_GAP short of each site, one instance per bond,
// with a per-bond radius, colour and opacity for the angle highlight, and a
// fade with distance for the bonds behind it.

import * as THREE from 'three';
import { describe, expect, it } from 'vitest';
import { BOND_END_GAP, bondMaterial, capsuleGeometry, setFadeRange, tubeMesh } from '../components/bondTubes';

const GPU = { radialSegments: 16, capRings: 3 };

// Where a tube's instance puts its body and its caps: centre, axis, radius,
// body length, and the full drawn length (body + both caps).
const instanceOf = (mesh, index) => {
    const matrix = new THREE.Matrix4();
    mesh.getMatrixAt(index, matrix);
    const position = new THREE.Vector3();
    const rotation = new THREE.Quaternion();
    const scale = new THREE.Vector3();
    matrix.decompose(position, rotation, scale);
    const axis = new THREE.Vector3(0, 1, 0).applyQuaternion(rotation);
    return { position, axis, radius: scale.x, body: scale.y, drawn: scale.y + 2 * scale.x };
};

describe('the capsule', () => {
    it('has a ring per cap step, unit normals, and cap offsets only on the caps', () => {
        const geometry = capsuleGeometry(16, 3);
        expect(geometry.attributes.position.count).toBe(2 * 4 * 16);
        const normals = geometry.attributes.normal;
        for (let i = 0; i < normals.count; i += 1) {
            expect(Math.hypot(normals.getX(i), normals.getY(i), normals.getZ(i))).toBeCloseTo(1, 6);
        }
        const offsets = [...geometry.attributes.capOffset.array];
        expect(Math.min(...offsets)).toBeCloseTo(-1, 6);
        expect(Math.max(...offsets)).toBeCloseTo(1, 6);
        // The body's two rings sit at y = ±0.5 with no offset.
        const ys = geometry.attributes.position;
        for (let i = 0; i < ys.count; i += 1) {
            expect(Math.abs(ys.getY(i))).toBe(0.5);
            if (offsets[i] === 0) expect(Math.hypot(ys.getX(i), ys.getZ(i))).toBeCloseTo(1, 6);
        }
    });
});

describe('a tube per bond', () => {
    it('stops BOND_END_GAP short of each site, caps included', () => {
        const segments = new Float64Array([0, 0, 0, 3, 0, 0, 1, 1, 1, 1, 1, 0.6]);
        const mesh = tubeMesh(segments, 0.1, bondMaterial({ color: '#000' }), { quality: GPU });
        expect(mesh.count).toBe(2);
        const long = instanceOf(mesh, 0);
        expect(long.position.toArray()).toEqual([1.5, 0, 0]);
        expect(Math.abs(long.axis.x)).toBeCloseTo(1, 6);
        expect(long.radius).toBeCloseTo(0.1, 6);
        expect(long.drawn).toBeCloseTo(3 - 2 * BOND_END_GAP, 6);
        // A 0.4 Å bond keeps half its length: the gap is at most a quarter.
        const short = instanceOf(mesh, 1);
        expect(short.drawn).toBeCloseTo(0.4 - 2 * 0.1, 5);
        expect(short.body).toBeGreaterThan(0);
    });

    it('takes a radius scale and an opacity per bond', () => {
        const segments = new Float64Array([0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0]);
        const mesh = tubeMesh(segments, 0.1, bondMaterial({ color: '#000', transparent: true }), {
            quality: GPU, scales: [1.8, 1.2], opacities: new Float64Array([1, 0.4])
        });
        expect(instanceOf(mesh, 0).radius).toBeCloseTo(0.18, 6);
        expect(instanceOf(mesh, 1).radius).toBeCloseTo(0.12, 6);
        // The drawn length does not depend on the radius.
        expect(instanceOf(mesh, 0).drawn).toBeCloseTo(2 - 2 * BOND_END_GAP, 6);
        expect([...mesh.geometry.attributes.instanceOpacity.array].map((v) => Number(v.toFixed(6)))).toEqual([1, 0.4]);
        // Without opacities every bond is opaque.
        const plain = tubeMesh(segments, 0.1, bondMaterial({ color: '#000' }), { quality: GPU });
        expect([...plain.geometry.attributes.instanceOpacity.array]).toEqual([1, 1]);
    });

    it('takes a colour per bond, and none unless asked', () => {
        const segments = new Float64Array([0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 2, 0]);
        const colors = [new THREE.Color('#1f6fd6'), new THREE.Color('#dce4ec')];
        const mesh = tubeMesh(segments, 0.1, bondMaterial({ color: '#ffffff' }), { quality: GPU, colors });
        const read = new THREE.Color();
        colors.forEach((color, index) => {
            mesh.getColorAt(index, read);
            expect(read.getHexString()).toBe(color.getHexString());
        });
        expect(tubeMesh(segments, 0.1, bondMaterial({ color: '#000' }), { quality: GPU }).instanceColor).toBe(null);
    });

    it('fades towards a colour between two view depths, only when asked', () => {
        const compiled = (material) => {
            const shader = {
                uniforms: {},
                vertexShader: THREE.ShaderLib.phong.vertexShader,
                fragmentShader: THREE.ShaderLib.phong.fragmentShader
            };
            material.onBeforeCompile(shader);
            return shader;
        };
        const material = bondMaterial({ color: '#aab3bf', specular: 0x000000, fade: '#e3e7ec' });
        const shader = compiled(material);
        // After the edge shading, before the fragment is written.
        expect(shader.fragmentShader).toMatch(/outgoingLight \*= mix\([\s\S]*outgoingLight = mix\(outgoingLight, bondFadeColor,[\s\S]*#include <opaque_fragment>/);
        expect(shader.uniforms.bondFadeColor.value.getHexString()).toBe('e3e7ec');
        setFadeRange(material, 12, 20);
        expect(shader.uniforms.bondFadeNear.value).toBe(12);
        expect(shader.uniforms.bondFadeFar.value).toBe(20);
        expect(material.specular.getHex()).toBe(0x000000);
        // Its own shader program, apart from the plain tubes'.
        expect(material.customProgramCacheKey()).toBe('bond-tube-1-fade');
        const plain = bondMaterial({ color: '#336699' });
        expect(compiled(plain).fragmentShader).not.toContain('bondFadeColor');
        expect(() => setFadeRange(plain, 1, 2)).not.toThrow();
    });

    it('patches the material for rounded caps, darker edges and the per-bond opacity', () => {
        const material = bondMaterial({ color: '#336699' });
        const shader = {
            vertexShader: THREE.ShaderLib.phong.vertexShader,
            fragmentShader: THREE.ShaderLib.phong.fragmentShader
        };
        material.onBeforeCompile(shader);
        // Every chunk it replaces was found (the three.js version still has it).
        expect(shader.vertexShader).toContain('capOffset * length(instanceMatrix[0].xyz)');
        expect(shader.vertexShader).toContain('vBondOpacity = instanceOpacity');
        expect(shader.vertexShader).not.toContain('#include <project_vertex>');
        expect(shader.vertexShader).not.toContain('#include <defaultnormal_vertex>');
        expect(shader.fragmentShader).toContain('diffuseColor.a *= vBondOpacity');
        expect(shader.fragmentShader).toContain('outgoingLight *= mix(');
        expect(material.customProgramCacheKey()).toBe('bond-tube-1');
    });
});
