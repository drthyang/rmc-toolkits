// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, it, expect } from 'vitest';

import { handlePcaMessage } from '../pcaKdeWorker.js';

// Minimal .rmc6f text with one reference site per element, n^3 box copies each.
const buildRmc6f = (elements, { supercell = [4, 4, 4], edge = 8, seed = 1 } = {}) => {
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
                    const base = [ix / sx, iy / sy, iz / sz];
                    const c = base.map((v) => v + 0.02 * rand());
                    lines.push(`${atom} ${element} [1] ${c[0].toFixed(8)} ${c[1].toFixed(8)} `
                        + `${c[2].toFixed(8)} ${ref + 1} ${ix} ${iy} ${iz}`);
                }
            }
        }
    });
    return lines.join('\n');
};

describe('pcaKdeWorker cache is content-addressed', () => {
    it('re-parses when the .rmc6f text changes (does not return the previous dataset)', async () => {
        const datasetA = buildRmc6f(['Ga', 'Se'], { seed: 1 });
        const datasetB = buildRmc6f(['Nb', 'Ta', 'O'], { seed: 2 });

        const a1 = await handlePcaMessage({ kind: 'sites' }, async () => datasetA);
        expect(a1.elements).toEqual(['Ga', 'Se']);

        // Loading a different dataset must reflect the new model, not the cached one.
        const b = await handlePcaMessage({ kind: 'sites' }, async () => datasetB);
        expect(b.elements).toEqual(['Nb', 'O', 'Ta']);
        expect(b.sites).toHaveLength(3);

        // And switching back returns the first dataset again (not a stale mix).
        const a2 = await handlePcaMessage({ kind: 'sites' }, async () => datasetA);
        expect(a2.elements).toEqual(['Ga', 'Se']);
        expect(a2.sites).toHaveLength(2);
    });

    it('serves a repeated identical request from cache (same result)', async () => {
        const dataset = buildRmc6f(['Ga', 'Se', 'Se'], { seed: 5 });
        const first = await handlePcaMessage({ kind: 'sites' }, async () => dataset);
        const second = await handlePcaMessage({ kind: 'sites' }, async () => dataset);
        expect(second.totalAtoms).toBe(first.totalAtoms);
        expect(second.referenceNumbers).toEqual(first.referenceNumbers);
    });

    it('kde requests also follow the current dataset', async () => {
        const datasetA = buildRmc6f(['Ga'], { seed: 7 });
        const datasetB = buildRmc6f(['Nb'], { seed: 8 });
        const kdeA = await handlePcaMessage({ kind: 'kde', referenceNumber: 1, grid: 12, projections: false }, async () => datasetA);
        expect(kdeA.element).toBe('Ga');
        const kdeB = await handlePcaMessage({ kind: 'kde', referenceNumber: 1, grid: 12, projections: false }, async () => datasetB);
        expect(kdeB.element).toBe('Nb');
    });
});

// Integrated 1.0 rule, the same in both runtimes: an atom line with a
// non-finite coordinate is skipped and counted by the shared grammar
// (parseRmc6fAtoms / iter_rmc6f_atoms), and the sites and orientation
// responses carry the report's warning -- the text /api/pca/sites and
// /api/pca/orientation return as parseWarning (tests/test_pca_api.py).
describe('pcaKdeWorker surfaces the .rmc6f parse report', () => {
    const withLine = (text, index, edit) => {
        const lines = text.split('\n');
        const atomsAt = lines.findIndex((line) => line.startsWith('Atoms'));
        lines[atomsAt + 1 + index] = edit(lines[atomsAt + 1 + index]);
        return { text: lines.join('\n'), line: lines[atomsAt + 1 + index] };
    };

    it('reports null for a clean file', async () => {
        const sites = await handlePcaMessage({ kind: 'sites' }, async () => buildRmc6f(['Ga'], { seed: 11 }));
        expect(sites.parseWarning).toBeNull();
    });

    it('skips a NaN line, keeps the other atoms, and names the line in sites and orientation', async () => {
        const { text, line } = withLine(buildRmc6f(['Se'], { seed: 12 }), 9, (row) => {
            const parts = row.split(' ');
            parts[3] = 'NaN';
            return parts.join(' ');
        });
        const sites = await handlePcaMessage({ kind: 'sites' }, async () => text);
        expect(sites.totalAtoms).toBe(63);
        expect(sites.parseWarning).toBe(`1 atom lines skipped for non-finite coordinates (first: '${line}')`);
        expect(sites.sites[0].uIso).toBeGreaterThan(0);
        expect(Number.isFinite(sites.sites[0].uIso)).toBe(true);

        const orientation = await handlePcaMessage(
            { kind: 'orientation', referenceNumber: 1, frequency: 3, geometry: false }, async () => text);
        expect(orientation.parseWarning).toBe(sites.parseWarning);
    });

    it('reads any Atoms-marker spelling and bare-CR files like the Dashboard parser', async () => {
        const clean = buildRmc6f(['Ga', 'Se'], { seed: 13 });
        const variants = [clean.replace('Atoms:', 'Atoms :'), clean.replace('Atoms:', 'atoms:'), clean.replace(/\n/g, '\r')];
        for (const text of variants) {
            const sites = await handlePcaMessage({ kind: 'sites' }, async () => text);
            expect(sites.totalAtoms).toBe(128);
            expect(sites.parseWarning).toBeNull();
        }
    });

    it('rejects a cell index outside the declared supercell, as iter_rmc6f_atoms does', async () => {
        const { text } = withLine(buildRmc6f(['Se'], { seed: 14 }), 0, (row) => {
            const parts = row.split(' ');
            parts[parts.length - 1] = '4';
            return parts.join(' ');
        });
        const sites = await handlePcaMessage({ kind: 'sites' }, async () => text);
        expect(sites.totalAtoms).toBe(63);
        expect(sites.parseWarning).toMatch(/^1 of 64 atom lines unparsed \(first: '1 Se \[1\] /);
    });

    it('says why when no atom can be parsed', async () => {
        const text = buildRmc6f(['Se'], { supercell: [2, 2, 2], seed: 15 })
            .split('\n')
            .map((row) => (/^\d+ Se /.test(row) ? row.replace(/^(\d+ Se \[1\] )\S+/, '$1inf') : row))
            .join('\n');
        await expect(handlePcaMessage({ kind: 'sites', file: { name: 'blown.rmc6f' } }, async () => text))
            .rejects.toThrow("blown.rmc6f: no atoms could be parsed — 8 atom lines skipped for non-finite coordinates (first: '1 Se [1] inf ");
    });
});
