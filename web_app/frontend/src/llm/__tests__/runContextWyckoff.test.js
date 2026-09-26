// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A Wyckoff label pairs a letter with the multiplicity of the cell the letter is read in.
// describeSymmetry reads letters in the standard cell the group is named in, which need not
// be the given cell: the four Ga of a lacunar spinel in its F-cubic cell are R3m's 3a (on
// hexagonal axes), rocksalt on its primitive cell is Fm-3m's 4a/4b. Each orbit carries both
// multiplicities: `size` in the given cell, `wyckoffMultiplicity` in the naming cell.

import { describe, expect, it } from 'vitest';
import { buildRunContext } from '../context/runContext';

const symmetryOf = (orbits) => ({ spaceGroup: 'R3m', spaceGroupNumber: 160, pointGroup: '3m', nSpace: 24, maxResidual: 0.01, orbits });
const sitesOf = (symmetry) => Object.fromEntries(buildRunContext({ runName: 'x', symmetry }).symmetry.sites.map((s) => [s.element, s]));

describe('Wyckoff labels in the symmetry context', () => {
    it('pairs the letter with the naming cell multiplicity, and keeps the given-cell size apart', () => {
        const sites = sitesOf(symmetryOf([
            { element: 'Nb', size: 12, site: '3m', rep: [0.6, 0.6, 0.6], wyckoff: 'b', wyckoffMultiplicity: 9, members: [1] },
            { element: 'Ga', size: 4, site: '3m', rep: [0.01, 0.01, 0.01], wyckoff: 'a', wyckoffMultiplicity: 3, members: [0] },
            { element: 'Se', size: 4, site: '3m', rep: [0.37, 0.37, 0.37], wyckoff: null, wyckoffMultiplicity: null, members: [2] },
        ]));
        expect(sites.Ga).toMatchObject({ multiplicity: 4, wyckoff: '3a' });
        expect(sites.Nb).toMatchObject({ multiplicity: 12, wyckoff: '9b' });
        expect(sites.Se.wyckoff).toBeUndefined();              // no letter: no label
        expect(sites.Se.multiplicity).toBe(4);
    });

    it('falls back to the given-cell size for an orbit without wyckoffMultiplicity', () => {
        const sites = sitesOf(symmetryOf([{ element: 'Ga', size: 4, site: '-43m', rep: [0, 0, 0], wyckoff: 'a', members: [0] }]));
        expect(sites.Ga.wyckoff).toBe('4a');
    });

    it('says why a structure was not analysed', () => {
        const reason = 'The average structure has 900 reference sites; symmetry detection runs in the browser.';
        const block = buildRunContext({
            runName: 'x',
            symmetry: { skipped: true, reason, spaceGroup: 'not analysed', pointGroup: '900 sites', maxResidual: Number.NaN, orbits: [] },
        }).symmetry;
        expect(block.note).toBe(reason);
        expect(block.max_residual_A).toBeUndefined();
        expect(buildRunContext({ runName: 'x', symmetry: symmetryOf([]) }).symmetry.note).toBeUndefined();
    });
});
