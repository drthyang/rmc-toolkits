// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Card references in answers: the spellings a model uses for a card id, the
// short [[id]] form rewritten as a card: link, and code left alone.

import { describe, expect, it } from 'vitest';
import { cardIdFromHref, citedCards, linkCardReferences, resolveCardId } from '../focus/cardRefs';

const CARDS = [
    { id: 'chi2-history', title: 'χ² history' },
    { id: 'plot:GTS_250K_PDFpartials.csv', title: 'g(r) partials' },
    { id: 'pca-sites', title: 'Site ellipsoids (3D)' }
];

describe('resolveCardId', () => {
    it('accepts the id, the scheme, brackets, case, the title and a bare plot file name', () => {
        expect(resolveCardId('chi2-history', CARDS)).toBe('chi2-history');
        expect(resolveCardId('card:chi2-history', CARDS)).toBe('chi2-history');
        expect(resolveCardId('`CHI2-History`', CARDS)).toBe('chi2-history');
        expect(resolveCardId('[χ² history]', CARDS)).toBe('chi2-history');
        expect(resolveCardId('GTS_250K_PDFpartials.csv', CARDS)).toBe('plot:GTS_250K_PDFpartials.csv');
        expect(resolveCardId('card:plot:GTS_250K_PDFpartials.csv', CARDS)).toBe('plot:GTS_250K_PDFpartials.csv');
    });

    it('returns null for anything else', () => {
        expect(resolveCardId('rietveld', CARDS)).toBeNull();
        expect(resolveCardId('', CARDS)).toBeNull();
        expect(resolveCardId('chi2-history', [])).toBeNull();
        expect(resolveCardId(null, CARDS)).toBeNull();
    });
});

describe('cardIdFromHref', () => {
    it('reads card: targets only', () => {
        expect(cardIdFromHref('card:pca-sites', CARDS)).toBe('pca-sites');
        expect(cardIdFromHref('card:Site%20ellipsoids%20(3D)', CARDS)).toBe('pca-sites');
        expect(cardIdFromHref('https://example.org', CARDS)).toBeNull();
        expect(cardIdFromHref('card:unknown', CARDS)).toBeNull();
    });
});

describe('linkCardReferences', () => {
    it('rewrites [[id]] and [[card:id]] for known cards, and leaves others and code alone', () => {
        const text = 'See [[chi2-history]] and [[card:pca-sites]], not [[nope]]; `[[chi2-history]]` stays.';
        expect(linkCardReferences(text, CARDS)).toBe(
            'See [χ² history](<card:chi2-history>) and [Site ellipsoids (3D)](<card:pca-sites>), not [[nope]]; `[[chi2-history]]` stays.'
        );
        expect(linkCardReferences('```\n[[chi2-history]]\n```', CARDS)).toBe('```\n[[chi2-history]]\n```');
    });

    it('wraps a card: target with spaces in angle brackets', () => {
        expect(linkCardReferences('[the history](card:χ² history)', CARDS))
            .toBe('[the history](<card:χ² history>)');
        expect(linkCardReferences('[g](card:chi2-history)', CARDS)).toBe('[g](card:chi2-history)');
    });

    it('lists the cards an answer cites, in order, once', () => {
        expect(citedCards('[a](card:pca-sites) then [[chi2-history]] and [b](card:pca-sites)', CARDS))
            .toEqual(['pca-sites', 'chi2-history']);
        expect(citedCards('no cards', CARDS)).toEqual([]);
    });
});
