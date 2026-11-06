// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* global process */

// The angle-histogram hover highlight (angleHighlight.js): the triplets it
// walks are exactly the ones the histogram counts — checked against an
// independent all-pairs search over every periodic image, bin by bin — they
// fold onto the right average sites and images, and a large cell keeps to the
// drawn region and says how many it shows. Each site triplet is weighted by
// its share (its angles in the range over all of its angles): a cell whose
// distorted octahedra are known shows exactly those in the distorted range.

import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import {
    binIndexOf, binsOfRange, highlightForBin, highlightForBins, highlightLabel, highlightStrength, MAX_HIGHLIGHT_TRIPLETS,
    MIN_HIGHLIGHT_SHARE, rangeText, siteAngleTriplets, unitCellOf
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

    it('keeps at most MAX_HIGHLIGHT_TRIPLETS, equal shares nearest the centre first', () => {
        const { parsed, atoms, siteFractional } = big;
        const unit = unitCellOf(parsed.latticeVectors, parsed.supercell);
        const options = tripletOptionsFromRequest({ end1: 'O', apex: 'Ti', end2: 'O', r12Min: 1.6, r12Max: 2.4, binWidth: 180 });
        const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options);
        expect(data.bins[0].triplets.length).toBeGreaterThan(MAX_HIGHLIGHT_TRIPLETS);
        const highlight = highlightForBin(data, 0, siteFractional, unit);
        // One bin holds every angle: every share is 1. One shared window:
        // both tubes of a triplet are A–B, each Ti–O bond drawn once (750
        // triplets: 50 whole octahedra, 300 bonds).
        expect(highlight.ab.length).toBe(300 * 6);
        expect([...highlight.abShare].every((share) => share === 1)).toBe(true);
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
        // A triplet's total counts the angles of a full bin too.
        expect(data.tripletTotals.reduce((a, b) => a + b, 0)).toBe(data.bins[0].total);
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
        expect(highlightLabel({ shown: 0, total: 0, lo: 59, hi: 60 })).toBe('No angles at 59–60°');
        expect(highlightLabel({ shown: 1200, total: 15000, lo: 0.25, hi: 0.5 })).toBe('1,200 of 15,000 angles at 0.25–0.5° shown');
    });
});

describe('an angle range', () => {
    const { parsed, atoms, siteFractional } = parse(rmc6fText({ ...PEROVSKITE, supercell: [3, 3, 3], jitter: 0.03 }));
    const unit = unitCellOf(parsed.latticeVectors, parsed.supercell);
    const request = { end1: 'O', apex: 'Ti', end2: 'O', r12Min: 1.6, r12Max: 2.4, binWidth: 1 };
    const options = tripletOptionsFromRequest(request);
    const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options);
    const counts = histogramOf(atoms, parsed.latticeVectors, options);

    it('holds exactly the triplets whose angle bins fall inside it, edges at bin edges', () => {
        const [first, last] = binsOfRange([85, 95], data.binWidth, data.nbins);
        expect([first, last]).toEqual([85, 94]);
        const inside = bruteForceTriplets(atoms, parsed.latticeVectors, options)
            .filter((key) => { const bin = Number(key.split('|')[0]); return bin >= first && bin <= last; });
        const highlight = highlightForBins(data, first, last, siteFractional, unit, { range: true });
        expect(inside.length).toBeGreaterThan(0);
        expect(highlight.total).toBe(inside.length);
        expect(highlight.total).toBe(counts.slice(85, 95).reduce((a, b) => a + b, 0));
        expect(highlight.shown).toBe(highlight.total);
        expect([highlight.lo, highlight.hi]).toEqual([85, 95]);
        expect(highlightLabel(highlight)).toBe(`${highlight.total.toLocaleString('en-US')} angles in 85–95°`);
    });

    it('draws a site triplet found in several of its bins once, its counts summed', () => {
        const highlight = highlightForBins(data, 0, data.nbins - 1, siteFractional, unit, { range: true });
        // One Ti site, six O neighbours: 15 site triplets over 6 bonds, each
        // bond drawn once, every share 1 (the range is the whole histogram).
        expect(highlight.ab.length).toBe(6 * 6);
        expect([...highlight.abShare]).toEqual([1, 1, 1, 1, 1, 1]);
        expect(highlight.total).toBe(counts.reduce((a, b) => a + b, 0));
        expect(highlight.shown).toBe(highlight.total);
        // A single bin is the hovered one's highlight.
        expect(highlightForBins(data, 90, 90, siteFractional, unit)).toEqual(highlightForBin(data, 90, siteFractional, unit));
    });

    it('maps a range to whole bins, off-grid edges to the nearest edge', () => {
        expect(binsOfRange([85, 95], 1, 180)).toEqual([85, 94]);
        expect(binsOfRange([85.2, 94.8], 1, 180)).toEqual([85, 94]);
        expect(binsOfRange([0, 180], 0.5, 360)).toEqual([0, 359]);
        expect(binsOfRange([90, 90.2], 1, 180)).toBe(null);
        expect(binsOfRange(null, 1, 180)).toBe(null);
        expect(rangeText(85, 95)).toBe('85–95°');
        expect(highlightLabel({ shown: 40, total: 2310, lo: 85, hi: 95, range: true })).toBe('40 of 2,310 angles in 85–95° shown');
    });
});


describe('the share weighting', () => {
    // Hand-made: sites on a line, unit cell 1 Å, two bins. Site 1 is the
    // apex; ends 2, 3, 4 sit at x = 2, 3, 4; distinct windows (A–B, B–C).
    const siteFractional = new Map([[1, [0, 0, 0]], [2, [2, 0, 0]], [3, [3, 0, 0]], [4, [4, 0, 0]], [5, [0, 5, 0]]]);
    const unit = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    const t = (id, ab, bc, count, d) => ({ id, apex: 1, ab: [ab, 0, 0, 0], bc: [bc, 0, 0, 0], count, d });
    // id 0: 8 of its 10 angles in bin 0 (share 0.8), far from the centre.
    // id 1: 2 of 40 (0.05), nearest the centre, shares the A–B bond 1→2.
    // id 2: 1 of 100 (0.01), below MIN_HIGHLIGHT_SHARE.
    const data = {
        binWidth: 1,
        sharedEnds: false,
        bins: [
            { total: 11, triplets: [t(1, 2, 4, 2, 0), t(0, 2, 3, 8, 9), t(2, 5, 4, 1, 1)] },
            { total: 140, triplets: [t(0, 2, 3, 2, 9), t(1, 2, 4, 38, 0), t(2, 5, 4, 99, 1)] }
        ],
        tripletTotals: [10, 40, 100]
    };
    const segmentEnds = (segments) => {
        const ends = [];
        for (let s = 0; s < segments.length; s += 6) ends.push(segments[s + 3] || segments[s + 4]);
        return ends;
    };

    it('weights each site triplet by its angles in the range over all its angles', () => {
        const highlight = highlightForBin(data, 0, siteFractional, unit);
        // The A–B bond 1→2 is used by id 0 (0.8) and id 1 (0.05): the larger.
        expect(segmentEnds(highlight.ab)).toEqual([2]);
        expect([...highlight.abShare]).toEqual([0.8]);
        expect(segmentEnds(highlight.bc)).toEqual([3, 4]);
        expect([...highlight.bcShare]).toEqual([0.8, 0.05]);
        // The legend counts are the histogram's, faint and undrawn included.
        expect(highlight.total).toBe(11);
        expect(highlight.shown).toBe(11);
        expect(highlightLabel(highlight)).toBe('11 angles at 0–1°');
    });

    it('does not draw a triplet below the threshold, but counts it', () => {
        const highlight = highlightForBin(data, 0, siteFractional, unit);
        // id 2's bonds (1→5 A–B, 1→4 B–C at 0.01) add nothing: 1→4 keeps id 1's 0.05.
        expect(segmentEnds(highlight.ab)).not.toContain(5);
        expect(MIN_HIGHLIGHT_SHARE).toBeGreaterThan(0.01);
        expect(MIN_HIGHLIGHT_SHARE).toBeLessThanOrEqual(0.05);
    });

    it('sums a range\'s bins before dividing, so the whole histogram is share 1', () => {
        const highlight = highlightForBins(data, 0, 1, siteFractional, unit, { range: true });
        expect([...highlight.abShare]).toEqual([1, 1]);
        expect(highlight.total).toBe(151);
        const second = highlightForBin(data, 1, siteFractional, unit);
        // Bin 1: id 2 at 0.99 now leads; id 1 at 0.95, id 0 at 0.2.
        expect(segmentEnds(second.ab)).toEqual([5, 2]);
        expect([...second.abShare]).toEqual([0.99, 0.95]);
    });

    it('keeps the largest shares when capped, not the nearest the centre', () => {
        const highlight = highlightForBin(data, 0, siteFractional, unit, { max: 1 });
        expect(highlight.shown).toBe(8);
        expect(highlightLabel(highlight)).toBe('8 of 11 angles at 0–1° shown');
        expect(segmentEnds(highlight.bc)).toEqual([3]);
        // Equal shares: nearest the centre first, as before.
        const even = { ...data, tripletTotals: [8, 2, 1] };
        expect(segmentEnds(highlightForBin(even, 0, siteFractional, unit, { max: 1 }).bc)).toEqual([4]);
    });

    it('draws a share by its square root, nothing below the threshold', () => {
        expect(highlightStrength(1)).toBe(1);
        expect(highlightStrength(0.25)).toBe(0.5);
        expect(highlightStrength(0.04)).toBeCloseTo(0.2, 12);
        expect(highlightStrength(MIN_HIGHLIGHT_SHARE / 2)).toBe(0);
        expect(highlightStrength(1.2)).toBe(1);
    });

    it('draws a homonuclear bond seen from both its atoms once', () => {
        const line = new Map([[1, [0, 0, 0]], [2, [0.5, 0, 0]], [3, [0, 0.5, 0]]]);
        const both = {
            binWidth: 1,
            sharedEnds: true,
            bins: [{ total: 2, triplets: [
                { id: 0, apex: 1, ab: [2, 0, 0, 0], bc: [3, 0, 0, 0], count: 1, d: 0 },
                { id: 1, apex: 2, ab: [1, 0, 0, 0], bc: [3, 0, 0, 0], count: 1, d: 0 }
            ] }],
            tripletTotals: [1, 4]
        };
        // id 1's first bond, 2 → 1 at image 0, is 1 → 2 at image 0 seen
        // from site 2: one tube, at the larger share (id 0's 1, not 0.25).
        const shared = highlightForBin(both, 0, line, unit);
        expect(shared.ab.length).toBe(3 * 6);
        expect(shared.abShare[0]).toBe(1);
        // 2 → 1 at image (1,0,0) is another bond (0.5 → 1 along a).
        both.bins[0].triplets[1].ab = [1, 1, 0, 0];
        expect(highlightForBin(both, 0, line, unit).ab.length).toBe(4 * 6);
    });

    it('agrees with each triplet\'s own bins on the demo run', () => {
        const text = readFileSync(resolve(process.cwd(), 'public/demo/GTS_250K.rmc6f'), 'utf8');
        const { parsed, atoms, siteFractional } = parse(text);
        const options = tripletOptionsFromRequest({ end1: 'Se', apex: 'Ta', end2: 'Se', r12Min: 2.3, r12Max: 2.9, binWidth: 1 });
        const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options);
        const perTriplet = new Map();
        data.bins.forEach((bin) => bin.triplets.forEach(({ id, count }) => perTriplet.set(id, (perTriplet.get(id) ?? 0) + count)));
        perTriplet.forEach((count, id) => expect(data.tripletTotals[id]).toBe(count));
        // 16 Ta sites × 15 Se pairs.
        expect(data.tripletTotals).toHaveLength(240);
    });
});

describe('a cell with a known subset of distorted octahedra', () => {
    // Four Ti sites in a 2×2×1 perovskite cell; Ti site 1 sits 0.3 Å off
    // its octahedron's centre along c, so four of its cis angles open to
    // ~99° and four close to ~81°. The other three are regular (90°).
    const shift = 0.3 / 3.95;
    const basis = [];
    [[0, 0], [1, 0], [0, 1], [1, 1]].forEach(([i, j], n) => {
        const at = (f) => [(i + f[0]) / 2, (j + f[1]) / 2, f[2]];
        basis.push(['Ti', at([0.5, 0.5, 0.5 + (n === 0 ? shift : 0)])]);
        basis.push(['Sr', at([0, 0, 0])], ['O', at([0.5, 0.5, 0])], ['O', at([0.5, 0, 0.5])], ['O', at([0, 0.5, 0.5])]);
    });
    const unit = [[7.8, 0, 0], [0, 7.8, 0], [0, 0, 3.95]];
    const { parsed, atoms, siteFractional } = parse(rmc6fText({ basis, unit, supercell: [5, 5, 5], jitter: 0.012 }));
    const cell = unitCellOf(parsed.latticeVectors, parsed.supercell);
    const options = tripletOptionsFromRequest({ end1: 'O', apex: 'Ti', end2: 'O', r12Min: 1.6, r12Max: 2.4, binWidth: 1 });
    const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options);
    const distortedApex = parsed.sites.find((site) => site.element === 'Ti' && site.siteFractional[2] > 0.55);
    const cart = (f) => [0, 1, 2].map((axis) => f[0] * cell[0][axis] + f[1] * cell[1][axis] + f[2] * cell[2][axis]);
    const apexAt = cart(distortedApex.siteFractional);
    // Per drawn bond, its share: on the distorted octahedron (its apex end) or not.
    const split = (highlight) => {
        const out = { distorted: [], regular: [] };
        for (let s = 0; s < highlight.ab.length; s += 6) {
            const fromDistorted = [0, 1, 2].every((axis) => Math.abs(highlight.ab[s + axis] - apexAt[axis]) < 1e-9);
            out[fromDistorted ? 'distorted' : 'regular'].push(highlight.abShare[s / 6]);
        }
        return out;
    };
    const mean = (values) => values.reduce((a, b) => a + b, 0) / values.length;
    const [first, last] = binsOfRange([78, 85], data.binWidth, data.nbins);

    it('draws the distorted octahedron strongly, the regular ones faintly at most', () => {
        expect(distortedApex).toBeTruthy();
        const highlight = highlightForBins(data, first, last, siteFractional, cell, { range: true });
        const { distorted, regular } = split(highlight);
        // Its four equatorial bonds and the short apical one: the four ~81°
        // angles lie (almost) wholly in the range; the long apical bond
        // takes part only in ~99° and 180° angles, so it is not drawn.
        expect(distorted).toHaveLength(5);
        expect(Math.min(...distorted)).toBeGreaterThan(0.8);
        // The regular octahedra reach 78–85° only with their tails.
        expect(Math.max(...regular)).toBeLessThan(0.1);
        expect(mean(distorted) / mean(regular)).toBeGreaterThan(10);
        // As drawn (square root): at least three times brighter.
        expect(Math.min(...distorted.map(highlightStrength)) / Math.max(...regular.map(highlightStrength))).toBeGreaterThan(3);
        expect(highlight.shown).toBe(highlight.total);
    });

    it('marks only the distorted octahedron\'s angles in the range strongly', () => {
        const atDistorted = (arc) => [0, 1, 2].every((axis) => Math.abs(arc.apex[axis] - apexAt[axis]) < 1e-9);
        [[[78, 85], 81.2], [[95, 102], 98.8]].forEach(([range, at]) => {
            const [lo, hi] = binsOfRange(range, data.binWidth, data.nbins);
            const arcs = arcList(highlightForBins(data, lo, hi, siteFractional, cell, { range: true }));
            const strong = arcs.filter((arc) => highlightStrength(arc.share) >= 0.85);
            // Its four ~81° (or ~99°) angles, between the shifted Ti and its
            // four equatorial O.
            expect(strong).toHaveLength(4);
            strong.forEach((arc) => {
                expect(atDistorted(arc)).toBe(true);
                expect(Math.abs(arcDegrees(arc) - at)).toBeLessThan(1);
            });
            // The rest are the tails of ~89–90° angles, faint.
            arcs.filter((arc) => !strong.includes(arc)).forEach((arc) => {
                expect(Math.abs(arcDegrees(arc) - 89.5)).toBeLessThan(1.5);
                expect(highlightStrength(arc.share)).toBeLessThan(0.4);
            });
        });
    });

    it('with little thermal spread, gives arcs to its in-range angles only', () => {
        const still = parse(rmc6fText({ basis, unit, supercell: [5, 5, 5], jitter: 0.002 }));
        const stillCell = unitCellOf(still.parsed.latticeVectors, still.parsed.supercell);
        const stillData = siteAngleTriplets(still.atoms, still.parsed.latticeVectors, still.parsed.supercell, still.siteFractional, options);
        const [lo, hi] = binsOfRange([78, 85], stillData.binWidth, stillData.nbins);
        const arcs = arcList(highlightForBins(stillData, lo, hi, still.siteFractional, stillCell, { range: true }));
        const shifted = still.parsed.sites.find((site) => site.element === 'Ti' && site.siteFractional[2] > 0.55);
        const shiftedAt = pointKey([0, 1, 2].map((axis) => [0, 1, 2].reduce(
            (sum, k) => sum + shifted.siteFractional[k] * stillCell[k][axis], 0
        )));
        expect(arcs).toHaveLength(4);
        arcs.forEach((arc) => {
            expect(pointKey(arc.apex)).toBe(shiftedAt);
            expect(Math.abs(arcDegrees(arc) - 81.2)).toBeLessThan(0.5);
            expect(arc.share).toBe(1);
        });
    });

    it('where an unweighted highlight lit the regular octahedra as brightly', () => {
        // Every site triplet with any count in the range (the old rule):
        // the regular octahedra have angles there too.
        const touched = new Set();
        for (let bin = first; bin <= last; bin += 1) {
            data.bins[bin].triplets.forEach((triplet) => {
                if (triplet.apex !== distortedApex.referenceNumber) touched.add(triplet.apex);
            });
        }
        expect(touched.size).toBe(3);
    });
});

// The arcs (`arcs`, nine numbers each: apex, A end, C end) mark each drawn
// triplet's angle at its apex; they are the highlight's own triplets.
const arcList = (highlight) => Array.from(highlight.arcShare, (share, index) => ({
    apex: [...highlight.arcs.subarray(index * 9, index * 9 + 3)],
    a: [...highlight.arcs.subarray(index * 9 + 3, index * 9 + 6)],
    c: [...highlight.arcs.subarray(index * 9 + 6, index * 9 + 9)],
    share
}));
const arcDegrees = ({ apex, a, c }) => {
    const p = a.map((value, axis) => value - apex[axis]);
    const q = c.map((value, axis) => value - apex[axis]);
    return Math.acos((p[0] * q[0] + p[1] * q[1] + p[2] * q[2]) / (Math.hypot(...p) * Math.hypot(...q))) * 180 / Math.PI;
};
const pointKey = (p) => p.map((value) => value.toFixed(6)).join();

// Every drawn bond (apex → end) of a highlight, with its share.
const drawnBonds = (highlight) => {
    const bonds = new Map();
    [[highlight.ab, highlight.abShare], [highlight.bc, highlight.bcShare]].forEach(([segments, shares]) => {
        for (let s = 0; s < segments.length; s += 6) {
            const from = pointKey([...segments.subarray(s, s + 3)]);
            const to = pointKey([...segments.subarray(s + 3, s + 6)]);
            bonds.set(`${from}|${to}`, shares[s / 6]);
        }
    });
    return bonds;
};

describe('the angle arcs', () => {
    it('mark each drawn triplet once, at its share, the largest first', () => {
        const siteFractional = new Map([[1, [0, 0, 0]], [2, [2, 0, 0]], [3, [0, 3, 0]], [4, [0, 0, 4]], [5, [0, 5, 0]]]);
        const unit = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        const t = (id, ab, bc, count) => ({ id, apex: 1, ab: [ab, 0, 0, 0], bc: [bc, 0, 0, 0], count, d: 0 });
        const data = {
            binWidth: 1,
            sharedEnds: false,
            bins: [{ total: 11, triplets: [t(1, 2, 4, 2), t(0, 2, 3, 8), t(2, 5, 4, 1)] }],
            tripletTotals: [10, 40, 100]
        };
        const arcs = arcList(highlightForBin(data, 0, siteFractional, unit));
        // id 0 (0.8) then id 1 (0.05); id 2 (0.01) is under the threshold.
        expect(arcs.map((arc) => arc.share)).toEqual([0.8, 0.05]);
        expect(arcs[0]).toEqual({ apex: [0, 0, 0], a: [2, 0, 0], c: [0, 3, 0], share: 0.8 });
        expect(arcs[1]).toEqual({ apex: [0, 0, 0], a: [2, 0, 0], c: [0, 0, 4], share: 0.05 });
        // Capped: only the kept triplets have arcs.
        expect(arcList(highlightForBin(data, 0, siteFractional, unit, { max: 1 })).map((arc) => arc.share)).toEqual([0.8]);
        // id 3 is id 0's two bonds in the other windows' order, the same
        // angle: one arc, at the larger share (0.8, then its own 1).
        data.bins[0].triplets.push(t(3, 3, 2, 1));
        data.tripletTotals.push(2);
        expect(arcList(highlightForBin(data, 0, siteFractional, unit)).map((arc) => arc.share)).toEqual([0.8, 0.05]);
        data.tripletTotals[3] = 1;
        expect(arcList(highlightForBin(data, 0, siteFractional, unit)).map((arc) => arc.share)).toEqual([1, 0.05]);
    });

    it('are exactly the triplets whose bonds are drawn, on a capped cell', () => {
        // A P1 box of 6×6×6 perovskite cells (one site per atom).
        const n = 6;
        const basis = [];
        for (let i = 0; i < n; i += 1) {
            for (let j = 0; j < n; j += 1) {
                for (let k = 0; k < n; k += 1) {
                    PEROVSKITE.basis.forEach(([element, f]) => basis.push([element, [(i + f[0]) / n, (j + f[1]) / n, (k + f[2]) / n]]));
                }
            }
        }
        const unitOfBox = PEROVSKITE.unit.map((row) => row.map((value) => value * n));
        const { parsed, atoms, siteFractional } = parse(rmc6fText({ basis, unit: unitOfBox, supercell: [1, 1, 1], jitter: 0.004 }));
        const unit = unitCellOf(parsed.latticeVectors, parsed.supercell);
        const options = tripletOptionsFromRequest({ end1: 'O', apex: 'Ti', end2: 'O', r12Min: 1.6, r12Max: 2.4, binWidth: 1 });
        const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options);
        const [first, last] = binsOfRange([86, 94], data.binWidth, data.nbins);
        const highlight = highlightForBins(data, first, last, siteFractional, unit, { range: true });

        // The reference: the bins' site triplets by share, capped, thresholded.
        const counts = new Map();
        for (let bin = first; bin <= last; bin += 1) {
            data.bins[bin].triplets.forEach(({ id, count }) => counts.set(id, (counts.get(id) ?? 0) + count));
        }
        const shares = [...counts].map(([id, count]) => count / data.tripletTotals[id]).sort((p, q) => q - p);
        expect(shares.length).toBeGreaterThan(MAX_HIGHLIGHT_TRIPLETS);
        const drawn = shares.slice(0, MAX_HIGHLIGHT_TRIPLETS).filter((share) => share >= MIN_HIGHLIGHT_SHARE);
        const arcs = arcList(highlight);
        expect(arcs.map((arc) => arc.share)).toEqual(drawn);

        // Each arc's two legs are drawn bonds, and each bond's share is the
        // largest of the arcs that use it.
        const bonds = drawnBonds(highlight);
        const largest = new Map();
        arcs.forEach((arc) => [arc.a, arc.c].forEach((end) => {
            const key = `${pointKey(arc.apex)}|${pointKey(end)}`;
            expect(bonds.has(key)).toBe(true);
            largest.set(key, Math.max(largest.get(key) ?? 0, arc.share));
        }));
        expect(largest).toEqual(bonds);
        // Cis angles only: every arc is near 90°.
        arcs.forEach((arc) => expect(Math.abs(arcDegrees(arc) - 90)).toBeLessThan(4));
    });

    it('tell on the demo which angles at each Ta are in the range', () => {
        const text = readFileSync(resolve(process.cwd(), 'public/demo/GTS_250K.rmc6f'), 'utf8');
        const { parsed, atoms, siteFractional } = parse(text);
        const unit = unitCellOf(parsed.latticeVectors, parsed.supercell);
        const options = tripletOptionsFromRequest({ end1: 'Se', apex: 'Ta', end2: 'Se', r12Min: 2.0, r12Max: 3.0, binWidth: 1 });
        const data = siteAngleTriplets(atoms, parsed.latticeVectors, parsed.supercell, siteFractional, options);
        // Per range: per Ta site (16), its arcs by angle, strong (strength
        // ≥ 0.9) or faint.
        const arcsIn = (range) => {
            const [first, last] = binsOfRange(range, data.binWidth, data.nbins);
            const arcs = arcList(highlightForBins(data, first, last, siteFractional, unit, { range: true }));
            const apexes = new Map();
            arcs.forEach((arc) => {
                const key = pointKey(arc.apex);
                if (!apexes.has(key)) apexes.set(key, []);
                apexes.get(key).push(arc);
            });
            expect(apexes.size).toBe(16);
            return { arcs, apexes: [...apexes.values()] };
        };
        const strong = (arc) => highlightStrength(arc.share) >= 0.9;

        // 80–100°: each Ta's six ~88.5° cis angles are strong; its three ~75°
        // and three ~104° ones only faint (their spreads' tails).
        const cis = arcsIn([80, 100]);
        cis.apexes.forEach((arcs) => {
            expect(arcs).toHaveLength(12);
            const bright = arcs.filter(strong);
            expect(bright).toHaveLength(6);
            bright.forEach((arc) => expect(arcDegrees(arc)).toBeGreaterThan(87));
            bright.forEach((arc) => expect(arcDegrees(arc)).toBeLessThan(90));
            arcs.filter((arc) => !strong(arc)).forEach((arc) => {
                expect(Math.abs(arcDegrees(arc) - 88.5)).toBeGreaterThan(10);
                expect(highlightStrength(arc.share)).toBeLessThan(0.45);
            });
        });
        // 70–80° and 100–110°: the three ~75° and the three ~104° angles.
        [[[70, 80], 75.5], [[100, 110], 104]].forEach(([range, at]) => {
            arcsIn(range).apexes.forEach((arcs) => {
                expect(arcs).toHaveLength(3);
                arcs.forEach((arc) => {
                    expect(Math.abs(arcDegrees(arc) - at)).toBeLessThan(1.5);
                    expect(highlightStrength(arc.share)).toBeGreaterThan(0.85);
                });
            });
        });
        // Where the bonds could not tell: 80–100° lights all 96 Ta–Se bonds.
        const [first, last] = binsOfRange([80, 100], data.binWidth, data.nbins);
        expect(highlightForBins(data, first, last, siteFractional, unit, { range: true }).ab.length).toBe(96 * 6);
    });
});
