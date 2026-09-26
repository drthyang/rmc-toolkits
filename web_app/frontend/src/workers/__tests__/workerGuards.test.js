// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Static-mode workers refuse what the Flask routes refuse, and never leave the
// caller waiting:
// - pcaKdeWorker: orientation `frequency` must be an integer and `smoothing` an
//   integer in [0, 64] (/api/pca/orientation); kde and orientation results that
//   hold NaN/Infinity are an error with _strict_result_response's message; an
//   unknown `kind` is an error, not a silent KDE.
// - autoScaleWorker: manual mode needs a finite, non-zero `a` and a finite `b`
//   (/api/scaling/*), and a non-finite scaling result is an error
//   (_require_finite_scaling), never `ok: true` with NaN curves.
// - localStructureWorker: `maxPoints` is an integer, clamped to [100, 1e6]
//   (/api/structure).
// - every worker answers a null or non-object message with an error.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';

import { handlePcaMessage } from '../pcaKdeWorker.js';
import { NON_FINITE_RESULT_MESSAGE, requestNumber } from '../requestGuards.js';
import { makeConfig } from '../autoScale.js';

const buildRmc6f = (elements, { supercell = [3, 3, 3], edge = 8, seed = 1 } = {}) => {
    let state = seed >>> 0;
    const rand = () => {
        state = (state * 1103515245 + 12345) & 0x7fffffff;
        return state / 0x7fffffff - 0.5;
    };
    const [sx, sy, sz] = supercell;
    const lines = [
        `Supercell dimensions ${sx} ${sy} ${sz}`,
        'Lattice vectors (Ang):',
        `${edge * sx} 0 0`,
        `0 ${edge * sy} 0`,
        `0 0 ${edge * sz}`,
        'Atoms:'
    ];
    let atom = 0;
    elements.forEach((element, ref) => {
        for (let ix = 0; ix < sx; ix += 1) {
            for (let iy = 0; iy < sy; iy += 1) {
                for (let iz = 0; iz < sz; iz += 1) {
                    atom += 1;
                    const c = [ix / sx, iy / sy, iz / sz].map((v) => v + 0.02 * rand());
                    lines.push(`${atom} ${element} [1] ${c[0].toFixed(8)} ${c[1].toFixed(8)} `
                        + `${c[2].toFixed(8)} ${ref + 1} ${ix} ${iy} ${iz}`);
                }
            }
        }
    });
    return lines.join('\n');
};

// Load a worker module with a stubbed `self`, returning its onmessage and the posts.
const loadWorker = async (path) => {
    const posted = [];
    vi.stubGlobal('self', { postMessage: (message) => posted.push(message) });
    vi.resetModules();
    await import(path);
    const handler = globalThis.self.onmessage;
    const send = async (data) => {
        posted.length = 0;
        await handler({ data });
        // Let any awaited work inside the handler finish.
        for (let i = 0; i < 20 && posted.length === 0; i += 1) await new Promise((r) => setTimeout(r, 5));
        return posted[0];
    };
    return { send, posted };
};

afterEach(() => {
    vi.unstubAllGlobals();
});

describe('requestNumber mirrors app._number', () => {
    it('reads numbers and numeric strings, blank means the fallback', () => {
        expect(requestNumber(3, 'x')).toBe(3);
        expect(requestNumber(' 2.5 ', 'x')).toBe(2.5);
        expect(requestNumber('', 'x', { fallback: 7 })).toBe(7);
        expect(requestNumber(null, 'x', { fallback: null })).toBeNull();
        expect(requestNumber(50, 'x', { clamp: [100, 1000] })).toBe(100);
    });

    it('refuses what _number refuses, naming the parameter', () => {
        expect(() => requestNumber('abc', 'x')).toThrow("x must be a number, got 'abc'");
        expect(() => requestNumber(true, 'x')).toThrow('x must be a number, got true');
        expect(() => requestNumber([1], 'x')).toThrow('x must be a number');
        expect(() => requestNumber(NaN, 'x')).toThrow('x must be a finite number, got NaN');
        expect(() => requestNumber('inf', 'x')).toThrow("x must be a finite number, got 'inf'");
        expect(() => requestNumber(2.5, 'x', { integer: true })).toThrow('x must be an integer, got 2.5');
        expect(() => requestNumber(65, 'x', { le: 64 })).toThrow('x must be <= 64, got 65');
        expect(() => requestNumber(-1, 'x', { ge: 0 })).toThrow('x must be >= 0, got -1');
        expect(() => requestNumber(0, 'x', { gt: 0 })).toThrow('x must be > 0, got 0');
    });
});

describe('pcaKdeWorker request rules (as /api/pca/*)', () => {
    const dataset = buildRmc6f(['Ga', 'Se'], { seed: 21 });
    const get = async () => dataset;

    it('orientation smoothing must be an integer in [0, 64]', async () => {
        await expect(handlePcaMessage({ kind: 'orientation', referenceNumber: 1, smoothing: 1e9 }, get))
            .rejects.toThrow('smoothing must be <= 64, got 1000000000');
        await expect(handlePcaMessage({ kind: 'orientation', referenceNumber: 1, smoothing: 1.5 }, get))
            .rejects.toThrow('smoothing must be an integer, got 1.5');
        await expect(handlePcaMessage({ kind: 'orientation', referenceNumber: 1, smoothing: -1 }, get))
            .rejects.toThrow('smoothing must be >= 0, got -1');
    });

    it('orientation frequency must be an integer', async () => {
        await expect(handlePcaMessage({ kind: 'orientation', referenceNumber: 1, frequency: 2.5 }, get))
            .rejects.toThrow('frequency must be an integer, got 2.5');
        await expect(handlePcaMessage({ kind: 'orientation', referenceNumber: 1, frequency: 'many' }, get))
            .rejects.toThrow("frequency must be a number, got 'many'");
        const ok = await handlePcaMessage(
            { kind: 'orientation', referenceNumber: 1, frequency: '3', geometry: false }, get);
        expect(ok.frequency).toBe(3);
    });

    it('a kde result with NaN/Infinity is refused with the Flask message', async () => {
        for (const extreme of [{ bw: 1e-200 }, { extent: 1e300 }, { bw: 1e200 }]) {
            await expect(handlePcaMessage(
                { kind: 'kde', referenceNumber: 1, grid: 12, projections: false, ...extreme }, get
            )).rejects.toThrow(NON_FINITE_RESULT_MESSAGE);
        }
        const ok = await handlePcaMessage({ kind: 'kde', referenceNumber: 1, grid: 12, projections: false }, get);
        expect(ok.density.length).toBe(12 ** 3);
    });

    it('an unknown kind is an error, not a silent KDE', async () => {
        await expect(handlePcaMessage({ kind: 'bogus' }, get))
            .rejects.toThrow("unknown request kind 'bogus'");
    });

    it('a null message is an error the caller receives', async () => {
        await expect(handlePcaMessage(null, get)).rejects.toThrow('worker request must be an object, got null');
    });
});

describe('worker entry points always post an answer', () => {
    beforeEach(() => {
        vi.resetModules();
    });

    it('pcaKdeWorker posts an error for a null message and for a bad request', async () => {
        const { send } = await loadWorker('../pcaKdeWorker.js');
        const empty = await send(null);
        expect(empty.error).toMatch(/worker request must be an object/);
        const bad = await send({ id: 7, kind: 'bogus', text: 'x' });
        expect(bad).toEqual({ id: 7, error: expect.stringContaining("unknown request kind 'bogus'") });
    });

    it('localKdeWorker posts an error for a null message', async () => {
        const { send } = await loadWorker('../localKdeWorker.js');
        const empty = await send(null);
        expect(empty.error).toMatch(/worker request must be an object/);
    });

    it('localStructureWorker validates maxPoints and posts an error for a null message', async () => {
        const { send } = await loadWorker('../localStructureWorker.js');
        const empty = await send(null);
        expect(empty.error).toMatch(/worker request must be an object/);
        const file = { name: 'run.rmc6f', sourceFile: { text: async () => buildRmc6f(['Ga', 'Se'], { seed: 3 }) } };
        for (const [maxPoints, message] of [
            [NaN, 'maxPoints must be a finite number, got NaN'],
            [Infinity, 'maxPoints must be a finite number, got Infinity'],
            [2.5, 'maxPoints must be an integer, got 2.5'],
            ['lots', "maxPoints must be a number, got 'lots'"],
        ]) {
            const answer = await send({ id: 1, file, maxPoints });
            expect(answer).toEqual({ id: 1, error: message });
        }
        // Clamped to [100, 1e6] like /api/structure: -5 and 0 read as 100.
        for (const maxPoints of [-5, 0]) {
            const answer = await send({ id: 2, file, maxPoints });
            expect(answer.result.sampledAtoms ?? answer.result.points.length).toBe(54);
            expect(Number.isFinite(answer.result.sampleStride)).toBe(true);
        }
    });

    it('autoScaleWorker: null message, manual a/b rules and non-finite results', async () => {
        const { send } = await loadWorker('../autoScaleWorker.js');
        const empty = await send(null);
        expect(empty).toMatchObject({ ok: false, error: expect.stringMatching(/worker request must be an object/) });

        // A smooth synthetic S(Q) -> 1 on Q 0.5..25.
        const q = Array.from({ length: 1000 }, (_, i) => 0.5 + i * 0.0245);
        const sq = q.map((x) => 1 + 0.3 * Math.sin(2.7 * x) * Math.exp(-0.05 * x * x));
        const config = makeConfig({ qmin: 0.5, qmax: 25, rho0: 0.07, bAvgSq: 0.1, rmax: 20, nr: 400 });
        const job = (extra) => ({ id: 5, config, q, sq, mode: 'manual', enforcement: null, ...extra });

        for (const [extra, message] of [
            [{ a: undefined }, "manual mode requires a scale 'a'"],
            [{ a: NaN }, 'a must be a finite number, got NaN'],
            [{ a: 'abc' }, "a must be a number, got 'abc'"],
            [{ a: 0 }, "manual mode requires a finite, non-zero scale 'a', got 0"],
            [{ a: 1, b: NaN }, 'b must be a finite number, got NaN'],
        ]) {
            const answer = await send(job(extra));
            expect(answer).toEqual({ id: 5, ok: false, error: message });
        }

        const overflow = await send(job({ a: 1e308, b: 1e308 }));
        expect(overflow.ok).toBe(false);
        expect(overflow.error).toMatch(/^the scaling result contains NaN or Infinity \(a = 1e\+308, b = 1e\+308\)/);

        const fine = await send(job({ a: 1, b: 0 }));
        expect(fine.ok).toBe(true);
        expect(new Float64Array(fine.result.gk).every(Number.isFinite)).toBe(true);

        const unknownMode = await send(job({ mode: 'sideways', a: 1 }));
        expect(unknownMode).toEqual({ id: 5, ok: false, error: "mode must be 'auto' or 'manual', got 'sideways'" });
    });
});
