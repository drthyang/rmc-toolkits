import { describe, expect, it } from 'vitest';
import { structureFromRmc6f } from '../browserData';

// Minimal .rmc6f fixture: 2×1×1 supercell of a 10 Å cubic cell, one Ga
// reference site whose two supercell copies sit at within-cell fraction
// 0.5 ± delta on x (exact 0.25 / 0.75 on y / z). Box coords are
// (withinFrac + cellIndex) / supercell.
const rmc6fText = (delta) => [
    '(Version 6f format configuration file)',
    'Number of atoms: 2',
    'Supercell dimensions: 2 1 1',
    'Lattice vectors (Ang):',
    ' 20.0 0.0 0.0',
    ' 0.0 10.0 0.0',
    ' 0.0 0.0 10.0',
    'Atoms:',
    ` 1 Ga [1] ${(0.5 + delta) / 2} 0.25 0.75 1 0 0 0`,
    ` 2 Ga [1] ${(0.5 - delta + 1) / 2} 0.25 0.75 1 1 0 0`
].join('\n');

describe('structureFromRmc6f site displacement (dispA)', () => {
    it('recovers the mean site and the rms displacement about it', () => {
        const structure = structureFromRmc6f({ path: 'test.rmc6f', text: rmc6fText(0.02) });
        expect(structure.basis).toHaveLength(1);
        const site = structure.basis[0];
        expect(site.el).toBe('Ga');
        expect(site.frac[0]).toBeCloseTo(0.5, 6);
        expect(site.frac[1]).toBeCloseTo(0.25, 6);
        expect(site.frac[2]).toBeCloseTo(0.75, 6);
        // Two points at ±0.02 in cell fraction on a 10 Å edge: circular std
        // √(−2 ln cos(2πδ))/2π · 10 ≈ 0.2003 Å; y and z contribute zero.
        expect(site.dispA).toBeCloseTo(0.2003, 3);
    });

    it('reports zero displacement for coincident copies', () => {
        const structure = structureFromRmc6f({ path: 'test.rmc6f', text: rmc6fText(0) });
        expect(structure.basis[0].dispA).toBe(0);
    });

    it('handles a boundary-wrapping site (mean at 0 ≡ 1)', () => {
        const text = [
            'Supercell dimensions: 2 1 1',
            'Lattice vectors (Ang):',
            ' 20.0 0.0 0.0',
            ' 0.0 10.0 0.0',
            ' 0.0 0.0 10.0',
            'Atoms:',
            ' 1 Ga [1] 0.495 0.25 0.75 1 0 0 0',   // within-cell 0.99
            ' 2 Ga [1] 0.505 0.25 0.75 1 1 0 0'    // within-cell 0.01
        ].join('\n');
        const site = structureFromRmc6f({ path: 'test.rmc6f', text }).basis[0];
        const wrapped = Math.min(site.frac[0], 1 - site.frac[0]);
        expect(wrapped).toBeLessThan(1e-6);
        expect(site.dispA).toBeCloseTo(0.1, 3);   // ±0.01 frac on a 10 Å edge
    });
});
