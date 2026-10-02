// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The symmetry-averaged structure behind the Detected SG card's CIF download
// (averageStructure.js). What it must guarantee, whatever cell or origin the .rmc6f uses:
//   • exact — the operations close exactly, every translation is a 1/48 fraction, and each
//     site's position and U are invariant under its site symmetry to round-off;
//   • faithful — expanding the written sites through the operations puts an atom within the
//     reported largest shift of every site mean of the model (no site lost or invented);
//   • measured — free coordinates keep the averaged value, U is the pooled second moment.

import { describe, expect, it } from 'vitest';

import { structureFromRmc6f } from '../browserData.js';
import { describeSymmetry, toleranceLadder } from '../symmetryModel.js';
import { cellRows, exactGroup, niceDenominator, symmetryAveragedStructure, wrapTidy } from '../averageStructure.js';
import { inv3 } from '../symmetry.js';
import { wyckoffPositions, fitsForm } from '../wyckoff.js';
import {
    STRUCTURES, lacunarSpinel, redescribe, closureDefects, demoStructure, orbits,
} from './fixtures/symmetryStructures.js';

const SUPER = 4;
const cyc = (x) => x - Math.round(x);
const mulV = (M, v) => M.map((row) => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
const mul = (X, Y) => X.map((row) => [0, 1, 2].map((j) => row[0] * Y[0][j] + row[1] * Y[1][j] + row[2] * Y[2][j]));
const transpose = (M) => [0, 1, 2].map((i) => [0, 1, 2].map((j) => M[j][i]));
const cartLength = (A, d) => Math.hypot(...[0, 1, 2].map((k) => d[0] * A[0][k] + d[1] * A[1][k] + d[2] * A[2][k]));

const seeded = (seed) => {
    let s = seed >>> 0;
    const rand = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
    const gauss = () => { let u = 0; while (u === 0) u = rand(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rand()); };
    return { rand, gauss };
};

/**
 * A fixture ({A, basis}) as structureFromRmc6f returns it: a SUPER³ box whose sites carry
 * the per-site moments — the mean moved by `offset` (cell fractions) and Gaussian noise of
 * `sigma` Å, and an anisotropic Cartesian covariance (u² on the diagonal plus a random
 * symmetric part of `aniso` Å²), different for every site.
 */
function asStructure({ A, basis }, { sigma = 0, u = 0.09, aniso = 0.002, offset = [0, 0, 0], seed = 7 } = {}) {
    const { gauss } = seeded(seed);
    const Ainv = inv3(A);
    const n = SUPER ** 3;
    const toFrac = (d) => [0, 1, 2].map((j) => d[0] * Ainv[0][j] + d[1] * Ainv[1][j] + d[2] * Ainv[2][j]);
    return {
        source: 'run/test.rmc6f',
        latticeVectors: A.map((row) => row.map((v) => v * SUPER)),
        supercell: [SUPER, SUPER, SUPER],
        totalAtoms: basis.length * n,
        basis: basis.map((site, index) => {
            const noise = toFrac([gauss() * sigma, gauss() * sigma, gauss() * sigma]);
            const mean = site.frac.map((v, i) => v + offset[i] + noise[i]);
            const B = [0, 1, 2].map(() => [gauss(), gauss(), gauss()]);
            const C = [0, 1, 2].map((i) => [0, 1, 2].map((j) => (i === j ? u * u : 0)
                + aniso * (B[i][0] * B[j][0] + B[i][1] * B[j][1] + B[i][2] * B[j][2]) / 3));
            // Cartesian → fractional: V = A⁻ᵀ·C·A⁻¹.
            const covFrac = mul(mul(transpose(Ainv), C), Ainv);
            return {
                el: site.el,
                referenceNumber: index + 1,
                frac: mean.map(wrapTidy),
                mean,
                covFrac,
                count: n,
                elementCounts: site.elementCounts ?? { [site.el]: n },
                dispA: 0,
            };
        }),
    };
}

// Distinct images (mod 1) of a point under the operations.
const images = (ops, x) => {
    const out = [];
    for (const { R, t } of ops) {
        const p = mulV(R, x).map((v, i) => wrapTidy(v + t[i]));
        if (!out.some((q) => q.every((v, i) => Math.abs(cyc(v - p[i])) < 1e-7))) out.push(p);
    }
    return out;
};

// Metric of the written cell and the fractional covariance of a site (U_ij·a*_i·a*_j).
const metricOf = ({ a, b, c, alpha, beta, gamma }) => {
    const cos = (deg) => Math.cos((deg * Math.PI) / 180);
    return [[a * a, a * b * cos(gamma), a * c * cos(beta)], [a * b * cos(gamma), b * b, b * c * cos(alpha)], [a * c * cos(beta), b * c * cos(alpha), c * c]];
};
const fractionalCovariance = (site, G) => {
    const Gs = inv3(G);
    const aStar = [0, 1, 2].map((i) => Math.sqrt(Gs[i][i]));
    return site.U.map((row, i) => row.map((v, j) => v * aStar[i] * aStar[j]));
};

/** The exactness and faithfulness checks every export must pass. */
function expectSoundModel(model, structure) {
    const ops = model.operations;
    // Exact group, translations on the 1/48 grid.
    expect(closureDefects(ops, 1e-9)).toEqual([]);
    for (const { t } of ops) for (const v of t) expect(niceDenominator(v), `translation ${v}`).not.toBeNull();
    // Site symmetry: position and U invariant under every operation fixing the site.
    const G = metricOf(model.cell);
    let written = 0;
    for (const site of model.sites) {
        const V = fractionalCovariance(site, G);
        const stabilizer = ops.filter(({ R, t }) => mulV(R, site.x).every((v, i) => Math.abs(cyc(v + t[i] - site.x[i])) < 1e-9));
        expect(stabilizer.length * site.multiplicity).toBe(ops.length);
        expect(images(ops, site.x)).toHaveLength(site.multiplicity);
        const scale = Math.max(...V.flat().map(Math.abs));
        for (const { R } of stabilizer) {
            const W = mul(mul(R, V), transpose(R));
            W.forEach((row, i) => row.forEach((v, j) => expect(Math.abs(v - V[i][j])).toBeLessThan(1e-12 * scale + 1e-18)));
        }
        written += site.multiplicity;
    }
    // Faithful: every site mean of the model sits within the reported shift of a written
    // atom of its element, and the cell holds as many atoms as the model's sites map to.
    const p = model.provenance;
    expect(written).toBeCloseTo(structure.basis.length * p.ratio, 9);
    const Qinv = inv3(p.Q);
    const A = structure.latticeVectors.map((row, i) => row.map((v) => v / structure.supercell[i]));
    const Aout = cellRows(A, p.Q);
    const expanded = model.sites.flatMap((site) => images(ops, site.x).map((x) => ({ element: site.element, x })));
    for (const b of structure.basis) {
        const x = mulV(Qinv, (b.mean ?? b.frac).map((v, i) => v + p.originShift[i]));
        const nearest = Math.min(...expanded.filter((e) => e.element === b.el)
            .map((e) => cartLength(Aout, e.x.map((v, i) => cyc(v - x[i])))));
        expect(nearest).toBeLessThanOrEqual(p.maxShiftA + 1e-9);
    }
}

describe('exact operations', () => {
    it('makes a noisy group exact without moving it', () => {
        // A hand-perturbed P2_1/c: exact closure, translations moved by the noise only.
        const A = [[7.1, 0, 0], [0, 5.3, 0], [0, 0, 9.4]];
        const p21c = [
            { R: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0.0002, -0.0001, 0.0003] },
            { R: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]], t: [0.0004, 0.5003, 0.4996] },
            { R: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]], t: [0.9997, 0.0002, 0.0001] },
            { R: [[1, 0, 0], [0, -1, 0], [0, 0, 1]], t: [0.0001, 0.4998, 0.5002] },
        ];
        const exact = exactGroup(p21c, A);
        expect(exact.defect).toBeLessThan(2e-3);
        expect(closureDefects(exact.ops, 1e-12)).toEqual([]);
        exact.ops.forEach((o, g) => o.t.forEach((v, i) => expect(Math.abs(cyc(v - p21c[g].t[i]))).toBeLessThan(1e-3)));
    });

    it('refuses a set that does not close', () => {
        const A = [[5, 0, 0], [0, 5, 0], [0, 0, 5]];
        const notGroup = [
            { R: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0, 0, 0] },
            { R: [[0, -1, 0], [1, 0, 0], [0, 0, 1]], t: [0, 0, 0] },
        ];
        expect(exactGroup(notGroup, A)).toBeNull();
    });
});

describe('symmetry-averaged structure', () => {
    it('rocksalt on an arbitrary origin: exact special positions, isotropic U, nice operations', () => {
        const structure = asStructure(STRUCTURES.rocksalt(), { sigma: 0.01, offset: [0.0123, 0.0456, 0.0789] });
        const model = symmetryAveragedStructure(structure, 0.06);
        expect(model.spaceGroup.symbol).toBe('Fm-3m');
        expect(model.operations).toHaveLength(192);
        expectSoundModel(model, structure);
        expect(model.provenance.niceOrigin).toBe(true);
        expect(model.sites.map((s) => [s.element, s.multiplicity])).toEqual([['Cl', 4], ['Na', 4]]);
        for (const site of model.sites) {
            // On a special position of Fm-3m: every coordinate a multiple of 1/4.
            site.x.forEach((v) => expect(Math.abs(cyc(v * 4))).toBeLessThan(1e-9));
            // m-3m site symmetry: U isotropic, at least the within-site u² (0.09² Å²).
            const U = site.U;
            expect(U[0][1]).toBeCloseTo(0, 12);
            expect(U[1][1]).toBeCloseTo(U[0][0], 12);
            expect(U[2][2]).toBeCloseTo(U[0][0], 12);
            // u² plus the random anisotropic part (0.002 Å² on average) and the noise.
            expect(site.Ueq).toBeGreaterThan(0.0081);
            expect(site.Ueq).toBeLessThan(0.0081 + 0.0045);
        }
    });

    it('finds an origin on the screw axes of a group with no special position (P2_12_12_1)', () => {
        const ops = [
            (x, y, z) => [x, y, z], (x, y, z) => [-x + 0.5, -y, z + 0.5],
            (x, y, z) => [-x, y + 0.5, -z + 0.5], (x, y, z) => [x + 0.5, -y + 0.5, -z],
        ];
        const fixture = {
            A: [[5.1, 0, 0], [0, 6.3, 0], [0, 0, 7.7]],
            basis: orbits(ops, [['K', 0.137, 0.213, 0.061], ['Br', 0.412, 0.073, 0.311], ['Br', 0.29, 0.38, 0.83]]),
        };
        const structure = asStructure(fixture, { sigma: 0.006, offset: [0.0731, 0.2113, 0.0457] });
        const model = symmetryAveragedStructure(structure, 0.04);
        expect(model.spaceGroup.symbol).toBe('P2_12_12_1');
        expect(model.provenance.niceOrigin).toBe(true);
        expectSoundModel(model, structure);
    });

    it('puts the origin of a hexagonal box on a 3-fold axis wherever the box was cut', () => {
        const structure = asStructure(STRUCTURES.wurtzite(), { sigma: 0.004, offset: [0.1712, 0.0533, 0.2901] });
        const model = symmetryAveragedStructure(structure, 0.04);
        expect(model.spaceGroup.symbol).toBe('P6_3mc');
        expect(model.provenance.niceOrigin).toBe(true);
        expectSoundModel(model, structure);
        for (const site of model.sites) for (const v of site.x.slice(0, 2)) expect(Math.abs(cyc(v * 3))).toBeLessThan(1e-9);
    });

    it('keeps a box on a standard origin where it is (shift of the order of the noise)', () => {
        const structure = asStructure(STRUCTURES.rocksalt(), { sigma: 0.005 });
        const model = symmetryAveragedStructure(structure, 0.03);
        expect(model.provenance.originShiftA).toBeLessThan(0.01);
        const na = model.sites.find((s) => s.element === 'Na');
        na.x.forEach((v) => expect(v).toBe(0));
    });

    it('wurtzite: the polar z is measured, x and y exact, U of 3m form', () => {
        const structure = asStructure(STRUCTURES.wurtzite(), { sigma: 0.008, offset: [0, 0, 0.0371] });
        const model = symmetryAveragedStructure(structure, 0.05);
        expect(model.spaceGroup.symbol).toBe('P6_3mc');
        expectSoundModel(model, structure);
        const zn = model.sites.find((s) => s.element === 'Zn');
        const o = model.sites.find((s) => s.element === 'O');
        for (const v of [zn.x[0], zn.x[1], o.x[0], o.x[1]]) expect(Math.abs(cyc(v * 3))).toBeLessThan(1e-9);
        // The O–Zn separation along c is the model's u = 0.382 (± the noise), not snapped.
        expect(Math.abs(cyc(o.x[2] - zn.x[2]) - 0.382)).toBeLessThan(0.004);
        expect(Math.abs(cyc(o.x[2] - zn.x[2]) - 0.382)).toBeGreaterThan(0);
        for (const { U } of [zn, o]) {
            expect(U[1][1]).toBeCloseTo(U[0][0], 12);
            expect(U[0][1]).toBeCloseTo(U[0][0] / 2, 12);
            expect(U[0][2]).toBeCloseTo(0, 12);
            expect(U[1][2]).toBeCloseTo(0, 12);
        }
        expect(model.cell.gamma).toBeCloseTo(120, 9);
    });

    it('rutile: O on 4f (x,x,0) with its measured x and the m.mm form of U', () => {
        const structure = asStructure(STRUCTURES.rutile(), { sigma: 0.01 });
        const model = symmetryAveragedStructure(structure, 0.06);
        expect(model.spaceGroup.symbol).toBe('P4_2/mnm');
        expectSoundModel(model, structure);
        const o = model.sites.find((s) => s.element === 'O');
        expect(o.wyckoff).toBe('f');
        const f = wyckoffPositions(136).find((r) => r.letter === 'f');
        expect(fitsForm(f.form, o.x, 1e-9)).toBe(true);
        expect(Math.abs(o.x[0] - 0.305)).toBeLessThan(0.006);
        expect(o.U[1][1]).toBeCloseTo(o.U[0][0], 12);
        expect(o.U[0][2]).toBeCloseTo(0, 12);
        expect(o.U[1][2]).toBeCloseTo(0, 12);
    });

    it('carries an F-cubic R3m distortion into its hexagonal cell', () => {
        const structure = asStructure(lacunarSpinel(0.01), { sigma: 0.003 });
        const model = symmetryAveragedStructure(structure, 0.02);
        expect(model.spaceGroup.symbol).toBe('R3m');
        expect(model.spaceGroup.centring).toBe('R');
        expect(model.operations).toHaveLength(18);
        expect(model.cell.gamma).toBeCloseTo(120, 9);
        expect(model.cell.a).toBeCloseTo(model.cell.b, 9);
        expectSoundModel(model, structure);
        expect(model.sites.find((s) => s.element === 'Ga')).toMatchObject({ multiplicity: 3, wyckoff: 'a' });
    });

    it('reduces a supercell to the true cell', () => {
        // Rutile in a 1x1x2 cell: the standard cell is half of it.
        const structure = asStructure(redescribe(STRUCTURES.rutile(), [[1, 0, 0], [0, 1, 0], [0, 0, 2]]), { sigma: 0.005 });
        const model = symmetryAveragedStructure(structure, 0.03);
        expect(model.spaceGroup.symbol).toBe('P4_2/mnm');
        expect(model.provenance.ratio).toBeCloseTo(0.5, 9);
        expect(model.cell.c).toBeCloseTo(2.959, 9);
        expectSoundModel(model, structure);
        // Each orbit of the box (twice the sites) maps to one site of the true cell, occupancy 1.
        expect(model.sites.map((s) => [s.element, s.multiplicity, s.elements[0].occupancy]))
            .toEqual([['O', 4, 1], ['Ti', 2, 1]]);
        expect(model.formula).toEqual({ counts: { O: 2, Ti: 1 }, Z: 2 });
    });

    it('writes a group with no standard cell (a lower bound) in the .rmc6f cell', () => {
        // Rocksalt in a 2x1x1 cell: the cubic 3-folds are not lattice rotations of that cell,
        // so the card reports a lower bound with no number; the export keeps its operations.
        const structure = asStructure(redescribe(STRUCTURES.rocksalt(), [[2, 0, 0], [0, 1, 0], [0, 0, 1]]), { sigma: 0.005 });
        const label = describeSymmetry(structure, 0.03).spaceGroup;
        expect(label.startsWith('≥')).toBe(true);
        const model = symmetryAveragedStructure(structure, 0.03);
        expect(model.spaceGroup).toMatchObject({ label, symbol: null, number: null, standard: false });
        expect(model.provenance.ratio).toBe(1);
        expect(model.cell.a).toBeCloseTo(11.28, 9);
        expectSoundModel(model, structure);
    });

    it('builds the conventional cell from a primitive one', () => {
        const primitive = [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]];
        const structure = asStructure(redescribe(STRUCTURES.rocksalt(), primitive), { sigma: 0.005 });
        const model = symmetryAveragedStructure(structure, 0.03);
        expect(model.spaceGroup.symbol).toBe('Fm-3m');
        expect(model.provenance.ratio).toBeCloseTo(4, 9);
        expect(model.cell.a).toBeCloseTo(5.64, 9);
        expect(model.cell.alpha).toBeCloseTo(90, 9);
        expectSoundModel(model, structure);
        expect(model.sites.map((s) => [s.element, s.multiplicity, s.elements[0].occupancy]))
            .toEqual([['Cl', 4, 1], ['Na', 4, 1]]);
        expect(model.formula).toEqual({ counts: { Cl: 1, Na: 1 }, Z: 4 });
    });

    it('writes the occupancies of a mixed site and keeps its majority element', () => {
        const base = STRUCTURES.rocksalt();
        const n = SUPER ** 3;
        base.basis = base.basis.map((site) => (site.el === 'Na'
            ? { ...site, elementCounts: { Na: n - 4, K: 4 } } : site));
        const structure = asStructure(base, { sigma: 0.004 });
        const model = symmetryAveragedStructure(structure, 0.03);
        const mixed = model.sites.find((s) => s.elements.length === 2);
        expect(mixed.element).toBe('Na');
        expect(mixed.elements).toEqual([{ element: 'Na', occupancy: (n - 4) / n }, { element: 'K', occupancy: 4 / n }]);
        expect(model.formula.Z).toBe(4);
        expect(model.formula.counts.K).toBeCloseTo(4 / n, 12);
        expect(model.formula.counts.Na).toBeCloseTo((n - 4) / n, 12);
    });

    it('in P1 writes every site mean and its own covariance unchanged', () => {
        const structure = asStructure(STRUCTURES.pnmaPerovskite(), { sigma: 0.02 });
        const model = symmetryAveragedStructure(structure, 1e-6);
        expect(model.spaceGroup.symbol).toBe(describeSymmetry(structure, 1e-6).spaceGroup);
        expect(model.spaceGroup.symbol).toBe('P1');
        expect(model.sites).toHaveLength(structure.basis.length);
        expect(model.provenance.maxShiftA).toBe(0);
        expectSoundModel(model, structure);
        const G = metricOf(model.cell);
        for (const b of structure.basis) {
            const site = model.sites.find((s) => s.x.every((v, i) => Math.abs(cyc(v - b.mean[i])) < 1e-12));
            expect(site).toBeDefined();
            fractionalCovariance(site, G).forEach((row, i) => row.forEach((v, j) => expect(v).toBeCloseTo(b.covFrac[i][j], 14)));
        }
    });

    it('pools a distortion the group averages away into U', () => {
        // Tight: R3m keeps the Ga shift; loose: F-43m averages the four Ga copies, and their
        // spread about the cubic position appears as extra mean-square displacement. The
        // finder places the symmetry elements by least squares over all 52 sites, so the
        // cubic origin moves 4/52 of the way along the Ga shift s (0.01 of the 10.4 Å edge
        // on each axis): Ga is left (48/52)·s off its 4a site, isotropic over the four
        // tetrahedral images.
        const structure = asStructure(lacunarSpinel(0.01), { sigma: 0, aniso: 0 });
        const ladder = toleranceLadder(structure, 1.0);
        const cubic = ladder.find((b) => b.spaceGroup === 'F-43m');
        const model = symmetryAveragedStructure(structure, cubic.from + 1e-6);
        expectSoundModel(model, structure);
        const ga = model.sites.find((s) => s.element === 'Ga');
        expect(ga.Ueq).toBeCloseTo(0.09 ** 2 + ((48 / 52) * 0.01 * 10.4) ** 2, 10);
        expect(ga.U[0][1]).toBeCloseTo(0, 12);
    });
});

describe('the demo run (GaTa4Se8, 250 K)', () => {
    const structure = demoStructure();
    const ladder = toleranceLadder(structure, 1.0);

    it('exports every rung of the ladder as a sound model of that group', () => {
        expect(ladder.length).toBeGreaterThan(1);
        for (const brick of ladder) {
            const model = symmetryAveragedStructure(structure, brick.from + 1e-6);
            expect(model.spaceGroup.label).toBe(brick.spaceGroup);
            expectSoundModel(model, structure);
        }
    });

    it('F-43m: Ga on 4c, Ta and Se on 16e with measured x, cubic site U', () => {
        const brick = ladder[ladder.length - 1];
        expect(brick.spaceGroup).toBe('F-43m');
        const model = symmetryAveragedStructure(structure, brick.from + 1e-6);
        expect(model.formula).toEqual({ counts: { Ga: 1, Se: 8, Ta: 4 }, Z: 4 });
        expect(model.operations).toHaveLength(96);
        expect(model.sites.map((s) => `${s.element} ${s.multiplicity}${s.wyckoff}`)).toEqual(['Ga 4c', 'Se 16e', 'Se 16e', 'Ta 16e']);
        for (const site of model.sites.filter((s) => s.multiplicity === 16)) {
            expect(site.x[1]).toBeCloseTo(site.x[0], 12);
            expect(site.x[2]).toBeCloseTo(site.x[0], 12);
            expect(site.U[1][1]).toBeCloseTo(site.U[0][0], 12);
            expect(site.U[0][2]).toBeCloseTo(site.U[0][1], 12);
        }
        // Box built on the standard origin: the shift is noise.
        expect(model.provenance.originShiftA).toBeLessThan(0.005);
    });
});

describe('structureFromRmc6f per-site moments', () => {
    it('gives the arithmetic mean, the population covariance and the element counts', () => {
        // One site, three copies in a 3x1x1 box of a 10 Å cell: x = 0.40, 0.40, 0.55 (skewed),
        // two Ga and one Ta (a mixed site).
        const text = [
            'Supercell dimensions: 3 1 1',
            'Lattice vectors (Ang):',
            ' 30.0 0.0 0.0',
            ' 0.0 10.0 0.0',
            ' 0.0 0.0 10.0',
            'Atoms:',
            ` 1 Ga [1] ${0.40 / 3} 0.25 0.75 1 0 0 0`,
            ` 2 Ga [1] ${1.40 / 3} 0.25 0.75 1 1 0 0`,
            ` 3 Ta [1] ${2.55 / 3} 0.25 0.75 1 2 0 0`,
        ].join('\n');
        const site = structureFromRmc6f({ path: 'test.rmc6f', text }).basis[0];
        expect(site.count).toBe(3);
        expect(site.elementCounts).toEqual({ Ga: 2, Ta: 1 });
        expect(site.mean[0]).toBeCloseTo(0.45, 12);
        expect(Math.abs(site.frac[0] - 0.45)).toBeGreaterThan(1e-4);   // the circular mean leans to the mode
        expect(site.covFrac[0][0]).toBeCloseTo((0.05 ** 2 * 2 + 0.1 ** 2) / 3, 12);
        expect(site.covFrac[1][1]).toBeCloseTo(0, 12);
        expect(site.covFrac[0][1]).toBeCloseTo(0, 12);
    });
});
