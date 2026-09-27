// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Detected SG card: a structure the finder does not analyse has no fit
// (maxResidual NaN) and a reason, and the card must show the reason -- never
// 'fits to NaN Å'. Its help text must describe the 0.6.0 naming (standard
// setting from the group's own elements), not the old 'given cell, up to an
// axis permutation'.

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import ModelSummary from '../ModelSummary';
import { MAX_SYMMETRY_SITES } from '../../symmetryModel';

const decode = (html) => html
    .replace(/&#x27;/g, "'")
    .replace(/&quot;/g, '"')
    .replace(/&amp;/g, '&');

// A cubic box with more reference sites than the finder analyses.
const tooManySites = () => {
    const n = MAX_SYMMETRY_SITES + 1;
    const basis = Array.from({ length: n }, (_, i) => ({
        el: 'Se', referenceNumber: i + 1, frac: [(i % 50) / 50, Math.floor(i / 50) / 50, 0.5], dispA: 0.1
    }));
    return {
        source: 'run/big.rmc6f',
        totalAtoms: n,
        supercell: [1, 1, 1],
        latticeVectors: [[10, 0, 0], [0, 10, 0], [0, 0, 10]],
        elementCounts: { Se: n },
        atomIndices: { Se: basis.map((site) => site.referenceNumber) },
        basis,
    };
};

describe('Detected SG card', () => {
    it('shows why a structure was not analysed instead of a NaN fit', () => {
        const html = decode(renderToStaticMarkup(<ModelSummary structure={tooManySites()} />));
        expect(html).not.toMatch(/NaN/);
        expect(html).toMatch(new RegExp(`title="The average structure has ${MAX_SYMMETRY_SITES + 1} reference sites`));
    });

    it('describes the standard-setting naming in its help text', () => {
        const text = decode(renderToStaticMarkup(<ModelSummary structure={tooManySites()} />))
            .replace(/<[^>]+>/g, '').replace(/\s+/g, ' ');
        expect(text).not.toMatch(/up to an axis permutation/);
        expect(text).toMatch(/the symbol is reported in its standard setting/);
        expect(text).toMatch(/a symbol marked ≥ is a lower bound/);
    });
});
