// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { existsSync, readFileSync } from 'node:fs';
import { describe, it, expect } from 'vitest';
import { parseAtomLine } from '../rmc6f.js';
import { siteDisplacementsFromRmc6f, siteEllipsoids } from '../workers/pcaKde.js';

describe('parseAtomLine', () => {
    it('parses the current 10-field format (with a bracketed type label)', () => {
        const parts = '1 Ga [1] 0.024375 0.073409 0.074344 4 0 0 0'.trim().split(/\s+/);
        expect(parseAtomLine(parts)).toEqual({
            element: 'Ga',
            coords: [0.024375, 0.073409, 0.074344],
            referenceNumber: 4,
            cellIndices: [0, 0, 0],
        });
    });

    it('parses the older 9-field format (no type label)', () => {
        const parts = '3 Pb 0.998456155893590 0.001266296130703 0.049993349135976 5 0 0 2'
            .trim().split(/\s+/);
        const atom = parseAtomLine(parts);
        expect(atom.element).toBe('Pb');
        expect(atom.referenceNumber).toBe(5);
        expect(atom.cellIndices).toEqual([0, 0, 2]);
        expect(atom.coords[0]).toBeCloseTo(0.998456155893590, 12);
        expect(atom.coords[2]).toBeCloseTo(0.049993349135976, 12);
    });

    it('parses the coords-only format with a bracket label ( id el [n] x y z )', () => {
        expect(parseAtomLine('1 S [1] 0.979763 0.994536 0.987370'.split(/\s+/))).toEqual({
            element: 'S',
            coords: [0.979763, 0.994536, 0.987370],
            referenceNumber: null,
            cellIndices: null,
        });
    });

    it('parses the coords-only format without a label ( id el x y z )', () => {
        expect(parseAtomLine('7 F 0.10 0.20 0.30'.split(/\s+/))).toEqual({
            element: 'F',
            coords: [0.10, 0.20, 0.30],
            referenceNumber: null,
            cellIndices: null,
        });
    });

    it('rejects too-few fields and non-numeric coords / ref / cell', () => {
        expect(parseAtomLine('1 Pb 0.1 0.2 0.3 5 0 0'.split(/\s+/))).toBeNull(); // 8 fields (truncated full)
        expect(parseAtomLine('1 Pb x y z 5 0 0 0'.split(/\s+/))).toBeNull();      // coords NaN
        expect(parseAtomLine('1 Pb 0.1 0.2 0.3 x 0 0 0'.split(/\s+/))).toBeNull(); // ref NaN
        expect(parseAtomLine('1 Pb 0.1 0.2'.split(/\s+/))).toBeNull();            // < 5 fields
    });
});

// Reconstruction path: an even older file carries neither the reference-site nor the
// per-atom cell columns, so sites are recovered by folding into one unit cell and
// clustering. A synthetic 2×2×2 supercell with two Pb sites must come back as exactly
// two sites, each with one box copy per cell.
describe('siteDisplacementsFromRmc6f — coords-only (fold + cluster reconstruction)', () => {
    const SITE_A = [0.2, 0.2, 0.2];
    const SITE_B = [0.7, 0.7, 0.7];
    const buildCoordsOnly = () => {
        const N = 2;
        const lines = [];
        let id = 1;
        for (let cx = 0; cx < N; cx += 1) {
            for (let cy = 0; cy < N; cy += 1) {
                for (let cz = 0; cz < N; cz += 1) {
                    [SITE_A, SITE_B].forEach((uf) => {
                        // Small deterministic per-copy displacement (fraction of a cell)
                        // so each site's cloud has spread on every axis.
                        const jit = [(id % 3) - 1, ((id >> 1) % 3) - 1, ((id >> 2) % 3) - 1].map((d) => d * 0.02);
                        const cell = [cx, cy, cz];
                        const sf = uf.map((v, k) => (v + jit[k] + cell[k]) / N);
                        lines.push(`  ${id}  Pb  [1]  ${sf[0].toFixed(6)}  ${sf[1].toFixed(6)}  ${sf[2].toFixed(6)}`);
                        id += 1;
                    });
                }
            }
        }
        return [
            '(Version 6f format configuration file)',
            'Number of atoms: 16',
            'Supercell dimensions: 2 2 2',
            'Lattice vectors (Ang):',
            ' 10.0 0.0 0.0',
            ' 0.0 10.0 0.0',
            ' 0.0 0.0 10.0',
            'Atoms:',
            ...lines,
        ].join('\n');
    };

    it('recovers two sites of 8 copies each, flagged reconstructed', () => {
        const parsed = siteDisplacementsFromRmc6f(buildCoordsOnly(), { clusterThreshold: 1.5 });
        expect(parsed.reconstructed).toBe(true);
        expect(parsed.supercell).toEqual([2, 2, 2]);
        expect(parsed.sites).toHaveLength(2);
        parsed.sites.forEach((site) => {
            expect(site.element).toBe('Pb');
            expect(site.count).toBe(8);
            expect(site.copiesPerCell).toBe(8);   // 2×2×2 images per site
            expect(site.displacements).toHaveLength(8);
        });
        // Ordered by folded centroid, so site 1 is A (0.2) and site 2 is B (0.7).
        expect(parsed.sites[0].referenceNumber).toBe(1);
        expect(parsed.sites[1].referenceNumber).toBe(2);
        SITE_A.forEach((v, k) => expect(parsed.sites[0].siteFractional[k]).toBeCloseTo(v, 1));
        SITE_B.forEach((v, k) => expect(parsed.sites[1].siteFractional[k]).toBeCloseTo(v, 1));
        // Displacements are small local offsets (well under half a cell edge).
        parsed.sites.forEach((site) => site.displacements.forEach((d) => {
            expect(Math.hypot(...d)).toBeLessThan(1.0);
        }));
    });

    it('a coarser threshold cannot merge sites 4.3 Å apart', () => {
        const parsed = siteDisplacementsFromRmc6f(buildCoordsOnly(), { clusterThreshold: 2.5 });
        expect(parsed.sites).toHaveLength(2);
    });
});

// Clustering distances must use the full cell metric, not per-axis edge lengths. In an
// oblique cell the closest image of a site can be diagonal: here A=(0,0,0) and
// B=(½,0,½) are 2.24 Å apart via the (a,c) diagonal but 4.47 Å by per-axis reduction.
// A per-axis metric would mis-cluster at a threshold between the two.
describe('siteDisplacementsFromRmc6f — oblique (monoclinic) cell metric', () => {
    // 2×1×1 supercell; unit vectors a=(5,0,0) b=(0,5,0) c=(3,0,4) (a^c not 90°).
    const OBLIQUE = [
        '(Version 6f format configuration file)',
        'Number of atoms: 4',
        'Supercell dimensions: 2 1 1',
        'Lattice vectors (Ang):',
        ' 10.0 0.0 0.0',
        ' 0.0 5.0 0.0',
        ' 3.0 0.0 4.0',
        'Atoms:',
        '  1  Pb  [1]  0.000  0.000  0.000',
        '  2  Pb  [1]  0.500  0.000  0.000',
        '  3  Pb  [1]  0.250  0.000  0.500',
        '  4  Pb  [1]  0.750  0.000  0.500',
    ].join('\n');

    it('keeps the two sites apart at 1.5 Å (true separation 2.24 Å)', () => {
        const parsed = siteDisplacementsFromRmc6f(OBLIQUE, { clusterThreshold: 1.5 });
        expect(parsed.reconstructed).toBe(true);
        expect(parsed.sites).toHaveLength(2);
        parsed.sites.forEach((site) => expect(site.count).toBe(2));
    });

    it('merges them at 3.0 Å via the diagonal image (per-axis metric would not)', () => {
        const parsed = siteDisplacementsFromRmc6f(OBLIQUE, { clusterThreshold: 3.0 });
        expect(parsed.sites).toHaveLength(1);
        expect(parsed.sites[0].count).toBe(4);
    });
});

// Regression on the real shared file: SF6 at 190 K (plastic/rotator phase). The S
// molecular centres fold into two clean BCC sites (27 copies each); the orientationally
// disordered F fold into two shells of 162 (= 6 F × 27) at a bond-length threshold.
// The file lives under the gitignored data/ dir (user research data), so this block
// runs locally where the file is present and is skipped in CI / a fresh clone.
const SF6_PATH = new URL('../../../../data/old_rmc6f/rmcsf6_190k.rmc6f', import.meta.url);
const describeSf6 = existsSync(SF6_PATH) ? describe : describe.skip;
describeSf6('siteDisplacementsFromRmc6f — real coords-only SF6 file', () => {
    const text = readFileSync(SF6_PATH, 'utf8');

    it('reconstructs 2 S sites (27/27) and 2 F shells (162/27)', () => {
        const parsed = siteDisplacementsFromRmc6f(text, { clusterThreshold: 1.5 });
        expect(parsed.reconstructed).toBe(true);
        expect(parsed.supercell).toEqual([3, 3, 3]);
        const byElement = (el) => parsed.sites.filter((s) => s.element === el);
        const s = byElement('S');
        const f = byElement('F');
        expect(s).toHaveLength(2);
        expect(f).toHaveLength(2);
        s.forEach((site) => { expect(site.count).toBe(27); expect(site.copiesPerCell).toBe(27); });
        f.forEach((site) => { expect(site.count).toBe(162); });
        // Ellipsoids compute for every reconstructed site without error.
        const ellipsoids = siteEllipsoids(parsed.sites);
        expect(ellipsoids).toHaveLength(4);
        ellipsoids.forEach((e) => expect(Number.isFinite(e.uIso)).toBe(true));
        // The disordered F shell reads much larger than the S centre motion.
        const sRms = siteEllipsoids(s)[0].rmsIso;
        const fRms = siteEllipsoids(f)[0].rmsIso;
        expect(fRms).toBeGreaterThan(sRms * 2);
        // Whole-cluster (circular-mean) unwrap: the two physically-identical F shells
        // must reconstruct to nearly the same U_iso, and no displacement may exceed the
        // ~3 Å shell diameter. A member-relative unwrap wraps the shell's far cap a full
        // cell out — it made the two shells disagree (0.95 vs 1.41 Å²) and scattered 22
        // displacements past 3 Å.
        const fEll = siteEllipsoids(f);
        expect(Math.abs(fEll[0].uIso - fEll[1].uIso) / fEll[0].uIso).toBeLessThan(0.1);
        f.forEach((site) => site.displacements.forEach((d) => {
            expect(Math.hypot(...d)).toBeLessThan(3);
        }));
    });
});

// End-to-end: the site extractor must load an older, label-less `.rmc6f` file
// (this is the 2018-era format shared by users) the same as the current one.
describe('siteDisplacementsFromRmc6f — older 9-field .rmc6f', () => {
    const OLD_FORMAT = `(Version 6f format configuration file)
Metadata owner:     Maxim
Metadata date:      17-07-2018
Number of atoms:                     4
Supercell dimensions:                2  2  2
Cell (Ang/deg):   10.000000  10.000000  10.000000   90.000000   90.000000   90.000000
Lattice vectors (Ang):
  10.000000   0.000000   0.000000
   0.000000  10.000000   0.000000
   0.000000   0.000000  10.000000
Atoms:
       1   Pb  0.010  0.010  0.010     1   0   0   0
       2   Pb  0.020  0.480  0.020     1   0   1   0
       3   Pb  0.490  0.010  0.480     1   1   0   1
       4   Pb  0.510  0.510  0.510     1   1   1   1
`;

    it('reads supercell, lattice, and the reference site from an old-format file', () => {
        const parsed = siteDisplacementsFromRmc6f(OLD_FORMAT);
        expect(parsed.supercell).toEqual([2, 2, 2]);
        expect(parsed.latticeVectors[0]).toEqual([10, 0, 0]);
        expect(parsed.sites).toHaveLength(1);
        expect(parsed.sites[0].referenceNumber).toBe(1);
        expect(parsed.sites[0].element).toBe('Pb');
        expect(parsed.sites[0].count).toBe(4);
        expect(parsed.sites[0].displacements).toHaveLength(4);
    });
});
