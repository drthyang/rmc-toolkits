// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// WebGPU is not available under node, so the WGSL density shader is checked by
// emulation: the test packs the kernel with the same helpers computeDensityGpu
// uses (packKdeParams / packKdeSamples), replays KDE_WGSL's arithmetic in float32
// (Math.fround after every operation, as the shader rounds), and compares the
// result with the float64 CPU loop on the same kernel. A real GPU may fuse
// multiply-adds and its exp() is not correctly rounded, so this bounds the
// float32 error of the formula, not of a particular device.

import { describe, expect, it } from 'vitest';
import { KDE_WGSL, packKdeParams, packKdeSamples } from '../gpuKde';
import { computeDensityCpu, covariance, makeKernel } from '../localKdeWorker';

const f32 = Math.fround;

const makeRandom = (seed) => {
    let value = seed >>> 0;
    return () => {
        value += 0x6D2B79F5;
        let mixed = value;
        mixed = Math.imul(mixed ^ (mixed >>> 15), mixed | 1);
        mixed ^= mixed + Math.imul(mixed ^ (mixed >>> 7), mixed | 61);
        return ((mixed ^ (mixed >>> 14)) >>> 0) / 4294967296;
    };
};
const gaussian = (random) => Math.sqrt(-2 * Math.log(1 - random())) * Math.cos(2 * Math.PI * random());

// KDE_WGSL's main(), one invocation per node, on the packed float32 buffers.
const emulateShader = (paramData, sampleData) => {
    const floats = new Float32Array(paramData);
    const uints = new Uint32Array(paramData);
    const [w00, w10, w11, normalizer, xMin, yMin, xStep, yStep] = floats;
    const grid = uints[8];
    const count = uints[9];
    const density = [];
    for (let gidY = 0; gidY < grid; gidY += 1) {
        const row = [];
        const gy = f32(yMin + f32(f32(gidY) * yStep));
        for (let gidX = 0; gidX < grid; gidX += 1) {
            const gx = f32(xMin + f32(f32(gidX) * xStep));
            let sum = f32(0);
            for (let i = 0; i < count; i += 1) {
                const dx = f32(gx - sampleData[2 * i]);
                const dy = f32(gy - sampleData[2 * i + 1]);
                const w0 = f32(w00 * dx);
                const w1 = f32(f32(w10 * dx) + f32(w11 * dy));
                const e = f32(-0.5 * f32(f32(w0 * w0) + f32(w1 * w1)));
                if (e > -60) sum = f32(sum + f32(Math.exp(e)));
            }
            row.push(f32(sum * normalizer));
        }
        density.push(row);
    }
    return density;
};

const compare = (samples, factor, grid = 48) => {
    const kernel = makeKernel(covariance(samples), factor, 1 / samples.length);
    const args = { samples, kernel, grid, xMin: 0, yMin: 0, xStep: 1 / (grid - 1), yStep: 1 / (grid - 1) };
    const cpu = computeDensityCpu(args);
    const paramData = packKdeParams({ ...args, sampleCount: samples.length });
    const gpu = emulateShader(paramData, packKdeSamples(samples));
    const peak = Math.max(...cpu.flat());
    let worst = 0;
    cpu.forEach((row, y) => row.forEach((value, x) => {
        worst = Math.max(worst, Math.abs(gpu[y][x] - value) / peak);
    }));
    return { kernel, paramData, worst, peak };
};

describe('WGSL KDE shader (float32 emulation) vs the CPU loop', () => {
    it('consumes exactly the kernel fields the CPU loop uses', () => {
        const random = makeRandom(1);
        const samples = Array.from({ length: 200 }, () => [random(), random()]);
        const { kernel, paramData } = compare(samples, 0.03, 16);
        const floats = new Float32Array(paramData);
        expect([...floats.slice(0, 4)]).toEqual([kernel.w00, kernel.w10, kernel.w11, kernel.normalizer].map(f32));
        // The shader text evaluates the same whitened form as computeDensityCpu.
        expect(KDE_WGSL).toContain('let w0 = P.kern.x * dx;');
        expect(KDE_WGSL).toContain('let w1 = P.kern.y * dx + P.kern.z * dy;');
        expect(KDE_WGSL).toContain('let e = -0.5 * (w0 * w0 + w1 * w1);');
        expect(KDE_WGSL).toContain('if (e > -60.0)');
    });

    it('matches the CPU loop to float32 precision for a cell-filling slab', () => {
        const random = makeRandom(2);
        const samples = Array.from({ length: 1500 }, () => [random(), random()]);
        // Measured 2.8e-6 of the peak.
        expect(compare(samples, 0.03).worst).toBeLessThan(1e-5);
    });

    it('stays within float32 position rounding for a needle kernel', () => {
        // Two sites on the anti-diagonal at bw = 0.01: minor kernel sigma ~6e-5,
        // so the ~6e-8 float32 rounding of node and atom positions is ~1e-3 sigma.
        const random = makeRandom(3);
        const site = (u, v) => Array.from({ length: 300 }, () => [u + 0.006 * gaussian(random), v + 0.006 * gaussian(random)]);
        const samples = [...site(0.25, 0.75), ...site(0.75, 0.25)];
        const { kernel, worst } = compare(samples, 0.01, 64);
        expect(kernel.sigmaMinor).toBeLessThan(1e-4);
        // Measured 2.0e-4 of the peak.
        expect(worst).toBeLessThan(2e-3);
    });
});
