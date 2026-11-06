// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* global process */

// The angle-histogram hover highlight (angleHighlight.js): the triplets it
// walks are exactly the ones the histogram counts — checked against an
// independent all-pairs search over every periodic image, bin by bin — they
// fold onto the right average sites and images, and a large cell keeps to the
// drawn region and says how many it shows.

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import {
    binIndexOf, highlightForBin, highlightLabel, MAX_HIGHLIGHT_TRIPLETS, siteAngleTriplets, unitCellOf
} from '../angleHighlight';
import { buildBondNetwork } from '../bondNetwork';
import { siteDisplacementsFromRmc6f } from '../workers/pcaKde';
import { bondAngleSummary, tripletOptionsFromRequest, visitAngleTriplets } from '../workers/triplets';

// Seeded uniform numbers (mulberry32), so every run builds the same model.
const rng = (seed) => () => {
    seed |= 0;
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
};

// An .rmc6f of `basis` (unit-cell sites) repeated over `supercell` with a
// random displacement of up to `jitter` (unit-cell fractions) per atom.
const rmc6fText = ({ basis, unit, supercell, jitter, seed = 7 }) => {
    const random = rng(seed);
    const lattice = unit.map((row, i) => row.map((value) => value * supercell[i]));
    const lines = [];
    let index = 0;
    for (let i = 0; i < supercell[0]; i += 1) {
        for (let j = 0; j < supercell[1]; j += 1) {
            for (let k = 0; k < supercell[2]; k += 1) {
                basis.forEach(([element, f], site) => {
                    index += 1;
                    const coords = [i, j, k].map((cell, axis) => (
                        (cell + f[axis] + (random() - 0.5) * 2 * jitter) / supercell[axis]
                    ));
                    lines.push(`${index} ${element} [1] ${coords.map((v) => v.toFixed(12)).join(' ')} ${site + 1} ${i} ${j} ${k}`);
                });
            }
        }
    }
    return [
        '(Version 6f format configuration file)',
        `Number of atoms: ${index}`,
        `Supercell dimensions: ${supercell.join(' ')}`,
        'Lattice vectors (Ang):',
        ...lattice.map((row) => row.map((v) => v.toFixed(6)).join(' ')),
        'Atoms:',
        ...lines
    ].join('\n');
};

const parse = (text) => {
    const parsed = siteDisplacementsFromRmc6f(text);
    return {
        parsed,
        atoms: parsed.atomList,
        siteFractional: new Map(parsed.sites.map((site) => [site.referenceNumber, site.siteFractional]))
    };
};

// A small perovskite-like ABX3 cell with a triclinic tilt.
const PEROVSKITE = {
    basis: [
        ['Sr', [0, 0, 0]],
        ['Ti', [0.5, 0.5, 0.5]],
        ['O', [0.5, 0.5, 0]],
        ['O', [0.5, 0, 0.5]],
        ['O', [0, 0.5, 0.5]]
    ],
    unit: [[3.9, 0, 0], [0.25, 3.92, 0], [0.1, -0.15, 3.95]]
};

// Independent reference: every apex against every end atom over all images
// within two supercells, each angle binned from its plain arccos (random
// displacements put no angle on a bin edge).
const bruteForceTriplets = (atoms, lattice, options) => {
    const [a, b, c] = options.triplet;
    const [lo12, hi12] = options.bond12.map(Number);
    const [lo23, hi23] = (options.bond23 ?? options.bond12).map(Number);
    const nbins = Math.max(1, Math.round(180 / options.binWidth));
    const width = 180 / nbins;
    const wrapped = atoms.fractional.map((row) => row.map((v) => v - Math.floor(v)));
    const bondsTo = (apex, element, lo, hi) => {
        const out = [];
        atoms.elements.forEach((el, row) => {
            if (el !== element) return;
            for (let ix = -2; ix <= 2; ix += 1) {
                for (let iy = -2; iy <= 2; iy += 1) {
                    for (let iz = -2; iz <= 2; iz += 1) {
                        if (row === apex && !ix && !iy && !iz) continue;
                        const d = [ix, iy, iz].map((shift, axis) => wrapped[row][axis] - wrapped[apex][axis] + shift);
                        const v = [0, 1, 2].map((axis) => d[0] * lattice[0][axis] + d[1] * lattice[1][axis] + d[2] * lattice[2][axis]);
                        const length = Math.hypot(...v);
                        if (length > 0 && length >= lo && length <= hi) out.push({ key: `${row}@${ix},${iy},${iz}`, v, length });
                    }
                }
            }
        });
        return out;
    };
    const angle = (p, q) => Math.acos(Math.min(1, Math.max(-1,
        (p.v[0] * q.v[0] + p.v[1] * q.v[1] + p.v[2] * q.v[2]) / (p.length * q.length)))) * 180 / Math.PI;
    const keys = [];
    atoms.elements.forEach((el, apex) => {
        if (el !== b) return;
        const first = bondsTo(apex, a, lo12, hi12);
        const second = bondsTo(apex, c, lo23, hi23);
        if (a !== c) {
            first.forEach((p) => second.forEach((q) => {
                keys.push(`${Math.floor(angle(p, q) / width)}|${apex}|${p.key}|${q.key}`);
            }));
            return;
        }
        // Same end element: each unordered pair of distinct bond images
        // once, when one bond can take each window.
        const all = new Map();
        [...first, ...second].forEach((bond) => all.set(bond.key, bond));
        const in12 = new Set(first.map((bond) => bond.key));
        const in23 = new Set(second.map((bond) => bond.key));
        const list = [...all.values()];
        for (let i = 0; i < list.length; i += 1) {
            for (let j = i + 1; j < list.length; j += 1) {
                const [p, q] = [list[i], list[j]];
                if (!((in12.has(p.key) && in23.has(q.key)) || (in12.has(q.key) && in23.has(p.key)))) continue;
                const pair = [p.key, q.key].sort().join('|');
                keys.push(`${Math.floor(angle(p, q) / width)}|${apex}|${pair}`);
            }
        }
    });
    return keys.sort();
};

const engineTriplets = (atoms, lattice, options) => {
    const keys = [];
    const sameEnd = options.triplet[0] === options.triplet[2];
    visitAngleTriplets(atoms.fractional, atoms.elements, lattice, options, (bin, apex, one, two) => {
        const p = `${one.row}@${one.ix},${one.iy},${one.iz}`;
        const q = `${two.row}@${two.ix},${two.iy},${two.iz}`;
        keys.push(`${bin}|${apex}|${sameEnd ? [p, q].sort().join('|') : `${p}|${q}`}`);
    });
    return keys.sort();
};

const histogramOf = (atoms, lattice, options) => bondAngleSummary(atoms.fractional, atoms.elements, lattice, options).counts;

const binCounts = (keys, nbins) => {
    const counts = new Array(nbins).fill(0);
    keys.forEach((key) => { counts[Number(key.split('|')[0])] += 1; });
    return counts;
};

describe('visitAngleTriplets walks the histogram\'s triplets', () => {
    const { parsed, atoms } = parse(rmc6fText({ ...PEROVSKITE, supercell: [2, 2, 2], jitter: 0.03 }));
    const lattice = parsed.latticeVectors;

    const cases = [
        ['different ends (Sr–Ti–O)', { end1: 'Sr', apex: 'Ti', end2: 'O', r12Min: 3, r12Max: 3.8, r23Min: 1.7, r23Max: 2.3, binWidth: 1 }],
        ['same end, one window (O–Ti–O)', { end1: 'O', apex: 'Ti', end2: 'O', r12Min: 1.6, r12Max: 2.4, binWidth: 2 }],
        ['same end, two windows (O–Ti–O)', { end1: 'O', apex: 'Ti', end2: 'O', r12Min: 1.6, r12Max: 2.0, r23Min: 1.9, r23Max: 2.4, binWidth: 1 }],
        // rmax past half the 2-cell box: several images of one atom bond.
        ['homonuclear, long window (Sr–Sr–Sr)', { end1: 'Sr', apex: 'Sr', end2: 'Sr', r12Min: 3.5, r12Max: 5.8, binWidth: 1 }]
    ];
    cases.forEach(([name, request]) => {
        it(`${name}: the same triplets, in the same bins, as an all-pairs search and the histogram`, () => {
            const options = tripletOptionsFromRequest(request);
            const engine = engineTriplets(atoms, lattice, options);
            const reference = bruteForceTriplets(atoms, lattice, options);
            expect(engine.length).toBeGreaterThan(0);
            expect(engine).toEqual(reference);
            const counts = histogramOf(atoms, lattice, options);
            expect(binCounts(engine, counts.length)).toEqual(counts);
        });
    });

    it('names the A–B bond first when one end element has two windows', () => {
        const options = tripletOptionsFromRequest(cases[2][1]);
        let seen = 0;
        visitAngleTriplets(atoms.fractional, atoms.elements, lattice, options, (bin, apex, one, two) => {
            seen += 1;
            expect(one.member & 1).toBe(1);
            expect(two.member & 2).toBe(2);
        });
        expect(seen).toBeGreaterThan(0);
    });
});

// The angle at the apex site between two folded bonds.
const siteAngle = (siteFractional, unit, triplet) => {
    const cart = (f) => [0, 1, 2].map((axis) => f[0] * unit[0][axis] + f[1] * unit[1][axis] + f[2] * unit[2][axis]);
    const apex = cart(siteFractional.get(triplet.apex));
    const end = ([site, ...image]) => {
        const f = siteFractional.get(site);
        return cart([f[0] + image[0], f[1] + image[1], f[2] + image[2]]).map((value, axis) => value - apex[axis]);
    };
    const p = end(triplet.ab);
    const q = end(triplet.bc);
    return Math.acos((p[0] * q[0] + p[1] * q[1] + p[2] * q[2]) / (Math.hypot(...p) * Math.hypot(...q))) * 180 / Math.PI;
};

describe('siteAngleTriplets folds them onto the average sites', () => {
    it('counts every angle of every bin and, uncapped, shows them all', () => {
        const { parsed, atoms, siteFractional } = parse(rmc6fText({ ...PEROVSKITE, supercell: [3, 3, 3], jitter: 0.03 }));
        const options = tripletOptionsFromRequest({ end1: 'O', apex: 'Ti', end2: 'O', r12Min: 1.6, r12Max: 2.4, binWidth: 1 });
        const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options);
        const counts = histogramOf(atoms, parsed.latticeVectors, options);
        expect(data.bins.map((bin) => bin.total)).toEqual(counts);
        data.bins.forEach((bin) => {
            expect(bin.triplets.reduce((sum, triplet) => sum + triplet.count, 0)).toBe(bin.total);
        });
        // One Ti site and its six O neighbours: 12 right angles and 3 straight
        // ones per Ti, folded onto 15 site triplets (each O pair once).
        const folded = data.bins.flatMap((bin) => bin.triplets);
        const key = (triplet) => [triplet.ab, triplet.bc].map((end) => end.join()).sort().join('|');
        expect(new Set(folded.map(key)).size).toBe(15);
        expect(new Set(folded.map((triplet) => triplet.apex))).toEqual(new Set([2]));
        // Within a bin each site triplet appears once.
        data.bins.forEach((bin) => expect(new Set(bin.triplets.map(key)).size).toBe(bin.triplets.length));
    });

    it('puts each site triplet at the average angle of the bin it came from', () => {
        // Displacements far below a degree: a wrong site or image would move
        // an angle by tens of degrees.
        const { parsed, atoms, siteFractional } = parse(rmc6fText({ ...PEROVSKITE, supercell: [3, 3, 3], jitter: 1e-4 }));
        const unit = unitCellOf(parsed.latticeVectors, parsed.supercell);
        const options = tripletOptionsFromRequest({ end1: 'Sr', apex: 'Ti', end2: 'O', r12Min: 3, r12Max: 3.8, r23Min: 1.7, r23Max: 2.3, binWidth: 1 });
        const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options);
        let checked = 0;
        data.bins.forEach((bin, index) => bin.triplets.forEach((triplet) => {
            expect(Math.abs(siteAngle(siteFractional, unit, triplet) - (index + 0.5))).toBeLessThan(0.6);
            checked += 1;
        }));
        expect(checked).toBeGreaterThan(0);
    });

    it('agrees with the histogram on the demo run, every angle shown', () => {
        const text = readFileSync(resolve(process.cwd(), 'public/demo/GTS_250K.rmc6f'), 'utf8');
        const { parsed, atoms, siteFractional } = parse(text);
        const options = tripletOptionsFromRequest({ end1: 'Se', apex: 'Ta', end2: 'Se', r12Min: 2.3, r12Max: 2.9, binWidth: 1 });
        const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options);
        expect(data.bins.map((bin) => bin.total)).toEqual(histogramOf(atoms, parsed.latticeVectors, options));
        data.bins.forEach((bin) => {
            expect(bin.triplets.reduce((sum, triplet) => sum + triplet.count, 0)).toBe(bin.total);
        });
    });
});

describe('a large cell keeps to the drawn region', () => {
    // A P1 box: every atom is its own site (one unit cell), 6×6×6 perovskite
    // cells as one, so the bond network is capped.
    const big = (() => {
        const n = 7;
        const basis = [];
        for (let i = 0; i < n; i += 1) {
            for (let j = 0; j < n; j += 1) {
                for (let k = 0; k < n; k += 1) {
                    PEROVSKITE.basis.forEach(([element, f]) => basis.push([element, [(i + f[0]) / n, (j + f[1]) / n, (k + f[2]) / n]]));
                }
            }
        }
        const unit = PEROVSKITE.unit.map((row) => row.map((value) => value * n));
        return parse(rmc6fText({ basis, unit, supercell: [1, 1, 1], jitter: 0.004 }));
    })();

    it('shows only triplets whose tubes sit inside the drawn bonds, and counts the rest', () => {
        const { parsed, atoms, siteFractional } = big;
        const unit = unitCellOf(parsed.latticeVectors, parsed.supercell);
        const sites = parsed.sites.map((site) => ({ referenceNumber: site.referenceNumber, element: site.element, siteFractional: site.siteFractional }));
        const center = [0, 1, 2].map((axis) => 0.5 * (unit[0][axis] + unit[1][axis] + unit[2][axis]));
        const network = buildBondNetwork(sites, unit, [{ elements: ['O', 'Ti'], window: [1.6, 2.4], color: '#000' }], { center });
        expect(network.limited).toBe(true);
        const options = tripletOptionsFromRequest({ end1: 'O', apex: 'Ti', end2: 'O', r12Min: 1.6, r12Max: 2.4, binWidth: 1 });
        const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options,
            { center: network.center, radiusSq: network.radiusSq });
        const counts = histogramOf(atoms, parsed.latticeVectors, options);
        expect(data.bins.map((bin) => bin.total)).toEqual(counts);

        const cart = (f) => [0, 1, 2].map((axis) => f[0] * unit[0][axis] + f[1] * unit[1][axis] + f[2] * unit[2][axis]);
        const busiest = counts.indexOf(Math.max(...counts));
        const highlight = highlightForBin(data, busiest, siteFractional, unit);
        expect(highlight.total).toBe(counts[busiest]);
        expect(highlight.shown).toBeGreaterThan(0);
        expect(highlight.shown).toBeLessThan(highlight.total);
        for (let s = 0; s < highlight.ab.length; s += 6) {
            const mid = [0, 1, 2].map((axis) => (highlight.ab[s + axis] + highlight.ab[s + 3 + axis]) / 2);
            const d = mid.reduce((sum, value, axis) => sum + (value - network.center[axis]) ** 2, 0);
            expect(d).toBeLessThanOrEqual(network.radiusSq + 1e-9);
        }
        // Each drawn tube runs from its apex site to an O site in reach.
        const apexes = new Set(parsed.sites.filter((site) => site.element === 'Ti').map((site) => cart(site.siteFractional).map((v) => v.toFixed(6)).join()));
        expect(apexes.has([...highlight.ab.subarray(0, 3)].map((v) => v.toFixed(6)).join())).toBe(true);
        expect(highlightLabel(highlight)).toMatch(/^[\d,]+ of [\d,]+ angles at \d+–\d+° shown$/);
    });

    it('draws at most MAX_HIGHLIGHT_TRIPLETS, nearest the centre first', () => {
        const { parsed, atoms, siteFractional } = big;
        const unit = unitCellOf(parsed.latticeVectors, parsed.supercell);
        const options = tripletOptionsFromRequest({ end1: 'O', apex: 'Ti', end2: 'O', r12Min: 1.6, r12Max: 2.4, binWidth: 180 });
        const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options);
        expect(data.bins[0].triplets.length).toBeGreaterThan(MAX_HIGHLIGHT_TRIPLETS);
        const highlight = highlightForBin(data, 0, siteFractional, unit);
        // One shared window: both tubes of a triplet are A–B.
        expect(highlight.ab.length).toBe(MAX_HIGHLIGHT_TRIPLETS * 12);
        expect(highlight.bc.length).toBe(0);
        expect(highlight.shown).toBe(data.bins[0].triplets.slice(0, MAX_HIGHLIGHT_TRIPLETS).reduce((s, t) => s + t.count, 0));
        const distances = data.bins[0].triplets.map((triplet) => {
            const p = siteFractional.get(triplet.apex);
            return [0, 1, 2].reduce((sum, axis) => sum + (p[axis] - 0.5) ** 2, 0);
        });
        // Sorted by Cartesian distance; the cell is near-cubic, so the
        // fractional distances are near-sorted too: the drawn ones are central.
        expect(Math.max(...distances.slice(0, 50))).toBeLessThan(Math.min(...distances.slice(-50)));
    });

    it('keeps at most maxPerBin site triplets per bin, still counting every angle', () => {
        const { parsed, atoms, siteFractional } = big;
        const options = tripletOptionsFromRequest({ end1: 'O', apex: 'Ti', end2: 'O', r12Min: 1.6, r12Max: 2.4, binWidth: 180 });
        const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options, null, { maxPerBin: 100 });
        expect(data.bins[0].triplets).toHaveLength(100);
        expect(data.bins[0].total).toBe(histogramOf(atoms, parsed.latticeVectors, options)[0]);
    });
});

describe('the hovered bin and its legend line', () => {
    it('maps a bin centre, or any x inside the bin, to the bin', () => {
        expect(binIndexOf(88.5, 1, 180)).toBe(88);
        expect(binIndexOf(88.01, 1, 180)).toBe(88);
        expect(binIndexOf(180, 1, 180)).toBe(179);
        expect(binIndexOf(1.25, 0.5, 360)).toBe(2);
        expect(binIndexOf(null, 1, 180)).toBe(null);
    });

    it('says how many angles it shows, plainly', () => {
        expect(highlightLabel({ shown: 830, total: 830, lo: 88, hi: 89 })).toBe('830 angles at 88–89°');
        expect(highlightLabel({ shown: 12, total: 830, lo: 88, hi: 89 })).toBe('12 of 830 angles at 88–89° shown');
        expect(highlightLabel({ shown: 1, total: 1, lo: 88.5, hi: 89 })).toBe('1 angle at 88.5–89°');
        expect(highlightLabel({ shown: 1200, total: 15000, lo: 0.25, hi: 0.5 })).toBe('1,200 of 15,000 angles at 0.25–0.5° shown');
    });
});
