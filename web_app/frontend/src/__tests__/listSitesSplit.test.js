// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// list_sites shows the split (plan 4.1): known-answer tests. Synthetic clouds
// go through the engine's site rows (siteEllipsoids, the rows of the 'sites'
// request in both runtimes) and the tool's summary; the planted splits of the
// validation's dev cases go through the tool as the Copilot runs it. The
// held-out test cases are never built here.
//
// What is pinned: a split site's static displacement δ = √(rms1² − rms⊥²)
// with its sampling error reads back the planted δ, (δ̂ − δ)/SE ≈ N(0, 1)
// over many seeds; PC1 as the crystal direction it was planted along, in
// oblique cells too; the per-axis kurtosis on resolved axes only; and none of
// it where it does not apply (a Gaussian ellipsoid, however anisotropic, an
// isotropic or oblate site, a site with zero spread).

import { describe, expect, it } from 'vitest';
import { siteEllipsoids } from '../workers/pcaKde';
import { isSplitSite, splitDisplacement, splitMemberSites, summarizeSites } from '../copilotTools';
import { buildCase, DEV_CASES, mulberry32 } from '../../validation/cases';
import { loadRun } from '../../validation/runFolder';

const normal = (rng) => {
    let u = 0;
    while (u === 0) u = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * rng());
};
const norm = (v) => Math.hypot(v[0], v[1], v[2]);
const unit = (v) => v.map((x) => x / norm(v));
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
// Fractional [u v w] → Cartesian, cell rows a, b, c (r = Mᵀ·f).
const cartesian = (uvw, cell) => [0, 1, 2].map((k) => uvw[0] * cell[0][k] + uvw[1] * cell[1][k] + uvw[2] * cell[2][k]);
const angleDeg = (a, b) => (Math.acos(Math.min(1, Math.abs(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (norm(a) * norm(b)))) * 180) / Math.PI;

const CUBIC = [[4, 0, 0], [0, 4, 0], [0, 0, 4]];
const HEXAGONAL = [[4, 0, 0], [-2, 2 * Math.sqrt(3), 0], [0, 0, 6]];
const MONOCLINIC = [[5, 0, 0], [0, 6, 0], [-1.5, 0, 7]];

/**
 * One site's displacement cloud (Cartesian Å): a static displacement along
 * `axis` (a crystal [u v w] of `cell`) taken from `positions` ({ at, p },
 * mean zero), plus Gaussian thermal motion of widths `sigma` along that axis
 * and two axes across it.
 */
const plantedCloud = ({ n, cell, axis, positions, sigma }, seed) => {
    const rng = mulberry32(seed);
    const e1 = unit(cartesian(axis, cell));
    const e2 = unit(cross(e1, Math.abs(e1[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]));
    const e3 = cross(e1, e2);
    return Array.from({ length: n }, () => {
        let r = rng();
        let at = positions.at(-1).at;
        for (const position of positions) {
            if (r < position.p) { at = position.at; break; }
            r -= position.p;
        }
        const c = [at + sigma[0] * normal(rng), sigma[1] * normal(rng), sigma[2] * normal(rng)];
        return [0, 1, 2].map((k) => c[0] * e1[k] + c[1] * e2[k] + c[2] * e3[k]);
    });
};

// The engine row of one site, and the 'sites' payload around it (a 1×1×1
// supercell, so the unit cell is the lattice).
const siteRow = (displacements, element = 'Nb') => siteEllipsoids([{
    referenceNumber: 1, element, elementCounts: { [element]: displacements.length }, mixed: false,
    count: displacements.length, displacements, siteFractional: [0, 0, 0], copiesPerCell: null
}])[0];
const payloadOf = (rows, cell) => ({
    sites: rows, totalAtoms: rows.reduce((sum, row) => sum + row.count, 0), latticeVectors: cell, supercell: [1, 1, 1], parseWarning: null
});

const twoWells = (delta) => [{ at: delta, p: 0.5 }, { at: -delta, p: 0.5 }];
// [name, config, planted δ]: equal wells far apart and close (δ = 1.5σ), a
// small site, unequal occupancy (1/3 at +2b, 2/3 at −b: δ = b√2 is their
// RMS), and anisotropic thermal motion whose width along the split is the
// mean across it; along cubic, hexagonal and monoclinic directions.
const PLANTED = [
    ['equal wells, δ = 3.3σ, cubic [0 0 1]', { n: 512, cell: CUBIC, axis: [0, 0, 1], positions: twoWells(0.2), sigma: [0.06, 0.06, 0.06] }, 0.2],
    ['equal wells, δ = 1.5σ, cubic [1 1 1]', { n: 512, cell: CUBIC, axis: [1, 1, 1], positions: twoWells(0.09), sigma: [0.06, 0.06, 0.06] }, 0.09],
    ['216 copies, hexagonal [2 1 0]', { n: 216, cell: HEXAGONAL, axis: [2, 1, 0], positions: twoWells(0.15), sigma: [0.05, 0.05, 0.05] }, 0.15],
    ['unequal occupancy, monoclinic [1 0 2]', {
        n: 512, cell: MONOCLINIC, axis: [1, 0, 2], positions: [{ at: 0.21, p: 1 / 3 }, { at: -0.105, p: 2 / 3 }], sigma: [0.05, 0.05, 0.05]
    }, Math.sqrt(0.21 * 0.105)],
    ['anisotropic thermal motion, cubic [1 2 3]', {
        n: 512, cell: CUBIC, axis: [1, 2, 3], positions: twoWells(0.15), sigma: [Math.sqrt((0.07 ** 2 + 0.05 ** 2) / 2), 0.07, 0.05]
    }, 0.15]
];
const SEEDS = 1000;

describe('list_sites: a planted split site (synthetic clouds)', () => {
    it.each(PLANTED)('%s: δ within its sampling error over many seeds, PC1 kurtosis ≤ −0.5, along the planted [u v w]', (_, config, delta) => {
        const z = [];
        const offsets = [];
        let covered = 0;
        const stream = PLANTED.findIndex((entry) => entry[1] === config);
        for (let seed = 1; seed <= SEEDS; seed += 1) {
            const row = siteRow(plantedCloud(config, 7919 * seed + 104729 * stream));
            const { data } = summarizeSites(payloadOf([row], config.cell));
            const [site] = data.sites;
            // A kurtosis on each resolved axis, and only there.
            site.kurtosis_axes.forEach((value, axis) => expect(value === null).toBe(!row.axisResolved[axis]));
            // Split by the one rule the cluster summary uses too; δ only then.
            if (!isSplitSite(row)) {
                expect(site.delta_A).toBeUndefined();
                continue;
            }
            expect(site.kurtosis_axes[0]).toBeLessThanOrEqual(-0.5);
            const exact = splitDisplacement(row);
            expect(site.delta_A).toBe(Number(exact.delta.toPrecision(3)));
            expect(site.delta_se_A).toBe(Number(exact.se.toPrecision(2)));
            const score = (exact.delta - delta) / exact.se;
            z.push(score);
            if (Math.abs(score) <= 1.96) covered += 1;
            offsets.push(angleDeg(cartesian(site.pc1_uvw, config.cell), cartesian(config.axis, config.cell)));
        }
        // The rule finds the planted split in (almost) every seed.
        expect(z.length / SEEDS).toBeGreaterThanOrEqual(0.99);
        const mean = z.reduce((sum, value) => sum + value, 0) / z.length;
        const sd = Math.sqrt(z.reduce((sum, value) => sum + (value - mean) ** 2, 0) / (z.length - 1));
        // ≈ N(0, 1): bands of about three standard errors over 1000 seeds
        // (mean ±0.095; sd ±0.07; 95 % coverage ±0.021), the mean's widened
        // by the small upward bias of the largest sample eigenvalue (+0.02 to
        // +0.08 SE over 2000 seeds).
        expect(Math.abs(mean)).toBeLessThan(0.15);
        expect(sd).toBeGreaterThan(0.92);
        expect(sd).toBeLessThan(1.08);
        expect(covered / z.length).toBeGreaterThan(0.929);
        expect(covered / z.length).toBeLessThan(0.971);
        // PC1 is the planted direction: every seed within the scoring's 15°, most within 3°.
        offsets.sort((a, b) => a - b);
        expect(offsets.at(-1)).toBeLessThan(15);
        expect(offsets[Math.floor(offsets.length / 2)]).toBeLessThan(3);
    });

    it('[u v w] is a direct-lattice direction (not a plane normal), in the sense whose closest crystal axis is acute', () => {
        // Large clouds, so PC1 lies within a fraction of a degree of the planted axis.
        const large = (index, axis) => siteRow(plantedCloud({ ...PLANTED[index][1], n: 5000, ...(axis ? { axis } : {}) }, 11));
        const tenths = (uvw) => uvw.map((value) => Math.round(value * 10) / 10 + 0);
        const hexagonal = summarizeSites(payloadOf([large(2)], HEXAGONAL)).data.sites[0];
        expect(tenths(hexagonal.pc1_uvw)).toEqual([1, 0.5, 0]);
        // Not the normal of the (2 1 0) plane, which in this cell lies along [5 4 0] (G⁻¹·h).
        const monoclinic = summarizeSites(payloadOf([large(3)], MONOCLINIC)).data.sites[0];
        expect(tenths(monoclinic.pc1_uvw)).toEqual([0.5, 0, 1]);
        // Planted along −[1 0 2], it is still quoted as [0.5 0 1]: c is the closest axis.
        const reversed = summarizeSites(payloadOf([large(3, [-1, 0, -2])], MONOCLINIC)).data.sites[0];
        expect(tenths(reversed.pc1_uvw)).toEqual([0.5, 0, 1]);
        // The unit cell given to the tool is the one used: the same cloud read in a cubic cell is [1 0.58 0].
        const cubic = summarizeSites(payloadOf([large(2)], HEXAGONAL), { unitCell: CUBIC }).data.sites[0];
        expect(tenths(cubic.pc1_uvw)).toEqual([1, 0.6, 0]);
    });

    it('the sampling error is the two-well closed form, √((4δ²σ² + 3σ⁴)/n)/(2δ), on a large sample', () => {
        const n = 20000;
        const row = siteRow(plantedCloud({ n, cell: CUBIC, axis: [0, 0, 1], positions: twoWells(0.2), sigma: [0.06, 0.06, 0.06] }, 5));
        const { delta, se } = splitDisplacement(row);
        expect(Math.abs(delta - 0.2)).toBeLessThan(3 * se);
        const closed = Math.sqrt((4 * 0.2 ** 2 * 0.06 ** 2 + 3 * 0.06 ** 4) / n) / (2 * 0.2);
        expect(Math.abs(se / closed - 1)).toBeLessThan(0.05);
    });
});

describe('list_sites: no split where there is none', () => {
    it('an anisotropic Gaussian ellipsoid is not split: no delta_A, PC1 kurtosis near 0, PC1 along its long axis', () => {
        const prolate = { n: 512, cell: CUBIC, axis: [1, 1, 1], positions: [{ at: 0, p: 1 }], sigma: [0.1, 0.05, 0.05] };
        let split = 0;
        for (let seed = 1; seed <= SEEDS; seed += 1) {
            const row = siteRow(plantedCloud(prolate, 104729 * seed));
            if (isSplitSite(row)) split += 1;
        }
        // The rule's false-positive rate for a Gaussian is far below 0.1 % (κ̂ ≤ −3 SE).
        expect(split).toBe(0);
        const { data, checks } = summarizeSites(payloadOf([siteRow(plantedCloud(prolate, 3))], CUBIC));
        const [site] = data.sites;
        expect(site.delta_A).toBeUndefined();
        expect(site.delta_se_A).toBeUndefined();
        expect(Math.abs(site.kurtosis_axes[0])).toBeLessThan(0.5);
        expect(site.kurtosis_axes.slice(1)).toEqual([null, null]);
        expect(angleDeg(cartesian(site.pc1_uvw, CUBIC), [1, 1, 1])).toBeLessThan(5);
        expect(data.note).not.toMatch(/delta_A/);
        expect(data.note).toMatch(/pc1_uvw = PC1 as a direct-lattice direction \[u v w\] \(not a plane normal\)/);
        expect(checks.some((check) => /^Split/.test(check.text))).toBe(false);
    });

    it('an isotropic site resolves no axis: no per-axis kurtosis, no PC1 direction, no delta_A', () => {
        const row = siteRow(plantedCloud({ n: 512, cell: CUBIC, axis: [0, 0, 1], positions: [{ at: 0, p: 1 }], sigma: [0.07, 0.07, 0.07] }, 21));
        expect(row.axisResolved).toEqual([false, false, false]);
        const { data } = summarizeSites(payloadOf([row], CUBIC));
        expect(Object.keys(data.sites[0]).sort()).toEqual(['U_iso_A2', 'anisotropy', 'count', 'kurtosis', 'ref', 'rms_axes_A', 'site']);
        expect(data.note).toBe('kurtosis = Mardia multivariate excess kurtosis (0 Gaussian; >0 peaked/heavy-tailed; <0 flat or split); '
            + 'anisotropy = rms1/rms3; sites sorted by |kurtosis|, largest first');
    });

    it('an oblate site (a flattened Gaussian): only PC3, the resolved axis, has its kurtosis; PC1 has no direction', () => {
        const row = siteRow(plantedCloud({ n: 512, cell: CUBIC, axis: [0, 0, 1], positions: [{ at: 0, p: 1 }], sigma: [0.05, 0.09, 0.09] }, 22));
        expect(row.axisResolved).toEqual([false, false, true]);
        const [site] = summarizeSites(payloadOf([row], CUBIC)).data.sites;
        expect(site.kurtosis_axes[0]).toBeNull();
        expect(site.kurtosis_axes[1]).toBeNull();
        expect(Number.isFinite(site.kurtosis_axes[2])).toBe(true);
        expect(site.pc1_uvw).toBeUndefined();
        expect(site.delta_A).toBeUndefined();
    });

    it('zero spread (an averaged configuration): no axes, kurtosis or delta_A, and the warning', () => {
        const row = siteRow(Array.from({ length: 64 }, () => [0, 0, 0]));
        expect(row.zeroSpread).toBe(true);
        expect(isSplitSite(row)).toBe(false);
        expect(splitDisplacement(row)).toBeNull();
        const { data, checks } = summarizeSites(payloadOf([row], CUBIC));
        expect(data.sites[0]).toEqual({
            ref: 1, site: 'Nb', count: 64, U_iso_A2: 0, rms_axes_A: [0, 0, 0], anisotropy: null, kurtosis: null, zero_spread: true
        });
        expect(checks).toEqual([{ level: 'warn', text: '1 site(s) have zero spread (an averaged configuration?): no axes or kurtosis.' }]);
    });

    it('split sites come first, then |kurtosis|: trimming to the budget never drops a split site before another site', () => {
        const heavy = siteRow(plantedCloud({ n: 512, cell: CUBIC, axis: [0, 0, 1], positions: [{ at: 0, p: 0.9 }, { at: 0.3, p: 0.05 }, { at: -0.3, p: 0.05 }], sigma: [0.05, 0.05, 0.05] }, 31));
        const wide = siteRow(plantedCloud({ n: 512, cell: CUBIC, axis: [0, 0, 1], positions: twoWells(0.2), sigma: [0.06, 0.06, 0.06] }, 32));
        const close = siteRow(plantedCloud({ n: 512, cell: CUBIC, axis: [1, 1, 1], positions: twoWells(0.1), sigma: [0.06, 0.06, 0.06] }, 33));
        const rows = [heavy, close, wide].map((row, index) => ({ ...row, referenceNumber: index + 1 }));
        // The heavy-tailed site has the largest |Mardia kurtosis|, yet ranks last.
        expect(Math.abs(heavy.nonGaussianity)).toBeGreaterThan(Math.abs(wide.nonGaussianity));
        const { data, checks } = summarizeSites(payloadOf(rows, CUBIC));
        expect(data.sites.map((site) => site.ref)).toEqual([3, 2, 1]);
        expect(data.largest_abs_kurtosis_refs[0]).toBe(1);
        expect(data.note).toMatch(/split sites first \(largest delta_A first\), then by \|kurtosis\|, largest first$/);
        expect(checks).toEqual([
            {
                level: 'note',
                text: 'Split: sites 3 (Nb), 2 (Nb) have two positions along PC1 (PC1 excess kurtosis ≤ −0.5, at least 3 sampling errors '
                    + 'below 0); delta_A is the static displacement, along pc1_uvw.'
            },
            { level: 'ok', text: '3 sites, 1536 atoms.' }
        ]);
    });
});

describe('list_sites on the dev cases', () => {
    const listSites = async (id) => {
        const run = await loadRun(buildCase(DEV_CASES.find((item) => item.id === id), { partials: false }));
        const tool = run.tools.find((entry) => entry.name === 'list_sites');
        return { result: await tool.run({}), rows: (await run.engine('sites', {})).sites };
    };
    const splitRows = (data) => data.sites.filter((site) => site.delta_A != null);

    it('dev-offcentre: the planted two-site disorder reads back (±0.2 Å along [0 0 1]), and nothing else is split', async () => {
        const { result: { data, checks } } = await listSites('dev-offcentre');
        const [site, ...others] = splitRows(data);
        expect(others).toEqual([]);
        expect(site.site).toBe('Nb');
        expect(angleDeg(cartesian(site.pc1_uvw, CUBIC), [0, 0, 1])).toBeLessThan(1);
        expect(Math.abs(site.delta_A - 0.2)).toBeLessThan(3 * site.delta_se_A);
        expect(site.kurtosis_axes[0]).toBeLessThan(-0.5);
        expect(data.sites[0]).toBe(site);
        expect(checks.some((check) => check.text.startsWith(`Split: site ${site.ref} (Nb) has two positions along PC1`))).toBe(true);
        // The O sites are flattened Gaussians: PC3 (their bond) is the resolved axis.
        data.sites.filter((row) => row.site === 'O').forEach((row) => {
            expect(row.kurtosis_axes.slice(0, 2)).toEqual([null, null]);
            expect(row.pc1_uvw).toBeUndefined();
        });
    });

    it('dev-zone-boundary: the antiphase displacements give the same single-site split (±0.1 Å along [0 0 1])', async () => {
        const { result: { data } } = await listSites('dev-zone-boundary');
        const site = splitRows(data).find((row) => row.site === 'Ti');
        expect(Math.abs(site.delta_A - 0.1)).toBeLessThan(3 * site.delta_se_A);
        expect(angleDeg(cartesian(site.pc1_uvw, CUBIC), [0, 0, 1])).toBeLessThan(15);
    });

    it('dev-jahn-teller: list_sites and the cluster summary call the same sites split (one rule)', async () => {
        const { result: { data }, rows } = await listSites('dev-jahn-teller');
        const listed = splitRows(data).map((site) => site.ref).sort((a, b) => a - b);
        const member = splitMemberSites(rows, rows.map((row) => row.referenceNumber)).map((row) => row.referenceNumber);
        expect(listed).toEqual(member);
        expect(listed.length).toBeGreaterThanOrEqual(3);
        expect(splitRows(data).every((site) => site.site === 'F')).toBe(true);
    });

    it('the dev cases without a split site have none: anisotropic Gaussians (breathing Mo, perovskite O) and thermal motion only', async () => {
        for (const id of ['dev-breathing', 'dev-cell-wrap', 'dev-thermal-null', 'dev-symmetry-lowering', 'dev-chirality-trap', 'dev-chiral']) {
            const { result: { data } } = await listSites(id);
            expect(splitRows(data), id).toEqual([]);
        }
        // The breathing Mo: a resolved PC1 along its radial direction, Gaussian (κ near 0).
        const { result: { data } } = await listSites('dev-breathing');
        const mo = data.sites.filter((site) => site.site === 'Mo');
        expect(mo).toHaveLength(3);
        mo.forEach((site) => {
            expect(Math.abs(site.kurtosis_axes[0])).toBeLessThan(0.5);
            expect(site.pc1_uvw).toHaveLength(3);
        });
    });
});
