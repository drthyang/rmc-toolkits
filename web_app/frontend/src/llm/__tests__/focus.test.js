// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The card-focus controller: which cards breathe (data-copilot-focus), how a
// card on another page is reached (openPage, then wait for it to mount), the
// marks handed to the host, and every way a highlight ends — the next show,
// a click outside the Copilot, clear() — with the tools' undo callbacks.

import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { cardSelector, createFocusController, findCard, isCardShown } from '../focus/focus';

const CARDS = [
    { id: 'chi2-history', page: 'dashboard', title: 'χ² history' },
    { id: 'plot:GTS_250K_PDFpartials.csv', page: 'dashboard', title: 'g(r) partials' },
    { id: 'bond-angles', page: 'geometry', title: 'Bond-angle distribution' },
    { id: 'pca-sites', page: 'ellipsoids', title: 'Site ellipsoids' }
];

const focusOf = (root, id) => root.querySelector(cardSelector(id))?.getAttribute('data-copilot-focus') ?? null;

describe('the card-focus controller', () => {
    let root;
    let dock;

    beforeEach(() => {
        root = document.createElement('div');
        root.innerHTML = `
            <main>
              <section id="dashboard">
                <article data-copilot-card="chi2-history"></article>
                <article data-copilot-card="plot:GTS_250K_PDFpartials.csv"></article>
              </section>
              <section id="geometry" aria-hidden="true">
                <article data-copilot-card="bond-angles"></article>
              </section>
            </main>
            <aside class="llm-dock"><button id="chip">chip</button></aside>`;
        document.body.appendChild(root);
        dock = root.querySelector('.llm-dock');
    });

    afterEach(() => root.remove());

    const controller = (options = {}) => createFocusController({ root, getCards: () => CARDS, waitMs: 300, ...options });

    it('finds only cards on screen, quoting ids with colons and dots', () => {
        expect(findCard(root, 'plot:GTS_250K_PDFpartials.csv')).not.toBeNull();
        expect(findCard(root, 'bond-angles')).toBeNull();
        expect(isCardShown(root.querySelector(cardSelector('bond-angles')))).toBe(false);
        expect(findCard(root, 'nothing')).toBeNull();
    });

    it('lets the main card breathe as primary and the others as secondary, and hands marks to the host', async () => {
        const onMarks = vi.fn();
        const scroll = vi.fn(() => true);
        const focus = controller({ onMarks, scroll });
        const result = await focus.show(['plot:GTS_250K_PDFpartials.csv', 'chi2-history'], { marks: [{ x: 2.52, label: 'Ta–Se' }] });
        expect(result).toMatchObject({ shown: ['plot:GTS_250K_PDFpartials.csv', 'chi2-history'], missing: [], page: 'dashboard', switched: false, scrolled: true });
        expect(focusOf(root, 'plot:GTS_250K_PDFpartials.csv')).toBe('primary');
        expect(focusOf(root, 'chi2-history')).toBe('secondary');
        expect(scroll).toHaveBeenCalledWith(root.querySelector(cardSelector('plot:GTS_250K_PDFpartials.csv')));
        expect(onMarks).toHaveBeenLastCalledWith({ 'plot:GTS_250K_PDFpartials.csv': [{ x: 2.52, label: 'Ta–Se' }] });
        expect(focus.state().cards).toEqual(['plot:GTS_250K_PDFpartials.csv', 'chi2-history']);

        // The next show replaces it.
        await focus.show(['chi2-history']);
        expect(focusOf(root, 'plot:GTS_250K_PDFpartials.csv')).toBeNull();
        expect(focusOf(root, 'chi2-history')).toBe('primary');
        expect(onMarks).toHaveBeenLastCalledWith({});
    });

    it('opens the page of a card that is not on screen and waits for it to mount', async () => {
        const openPage = vi.fn((page) => {
            if (page === 'geometry') root.querySelector('#geometry').removeAttribute('aria-hidden');
        });
        const focus = controller({ openPage, scroll: () => false });
        const result = await focus.show(['bond-angles']);
        expect(openPage).toHaveBeenCalledWith('geometry');
        expect(result).toMatchObject({ shown: ['bond-angles'], switched: true });
        expect(focusOf(root, 'bond-angles')).toBe('primary');

        // A page visited for the first time mounts its cards later.
        const late = vi.fn((page) => {
            if (page !== 'ellipsoids') return;
            setTimeout(() => {
                const card = document.createElement('article');
                card.setAttribute('data-copilot-card', 'pca-sites');
                root.querySelector('main').appendChild(card);
            }, 80);
        });
        const later = controller({ openPage: late, scroll: () => false });
        const shown = await later.show(['pca-sites']);
        expect(shown).toMatchObject({ shown: ['pca-sites'], missing: [], switched: true });
        expect(focusOf(root, 'pca-sites')).toBe('primary');
    });

    it('reports a card that never appears, and a show replaced before it landed', async () => {
        const focus = controller({ openPage: () => {}, scroll: () => false, waitMs: 120 });
        const missing = await focus.show(['bond-angles']);
        expect(missing).toMatchObject({ shown: [], missing: ['bond-angles'], switched: true });

        const first = focus.show(['bond-angles']);
        await focus.show(['chi2-history']);
        expect((await first).superseded).toBe(true);
        expect(focusOf(root, 'chi2-history')).toBe('primary');
    });

    it('clears on a click outside the Copilot, never on its own clicks, and runs the tools\' undo', async () => {
        const undo = vi.fn();
        const focus = controller({ scroll: () => false });
        const detach = focus.attach();
        await focus.show(['chi2-history']);
        focus.onClear(undo);

        dock.querySelector('#chip').dispatchEvent(new Event('click', { bubbles: true }));
        expect(focusOf(root, 'chi2-history')).toBe('primary');
        expect(undo).not.toHaveBeenCalled();

        root.querySelector('#dashboard').dispatchEvent(new Event('pointerdown', { bubbles: true }));
        // A press alone (a touch-scroll toward the card) keeps it.
        expect(focusOf(root, 'chi2-history')).toBe('primary');
        root.querySelector('#dashboard').dispatchEvent(new Event('click', { bubbles: true }));
        expect(focusOf(root, 'chi2-history')).toBeNull();
        expect(undo).toHaveBeenCalledTimes(1);
        expect(focus.state()).toMatchObject({ cards: [], reason: 'user' });
        detach();

        // A new show runs the last show's undo too ("the agent moves on").
        const undo2 = vi.fn();
        await focus.show(['chi2-history']);
        focus.onClear(undo2);
        await focus.show(['plot:GTS_250K_PDFpartials.csv']);
        expect(undo2).toHaveBeenCalledTimes(1);
    });

    it('lets running analyses breathe softly without scrolling, below an explicit highlight', async () => {
        const scroll = vi.fn(() => true);
        const focus = controller({ scroll });
        focus.running('step-1', ['chi2-history', 'bond-angles']);
        expect(focusOf(root, 'chi2-history')).toBe('running');
        // Hidden pages are marked too, so the card breathes when the user opens it.
        expect(focusOf(root, 'bond-angles')).toBe('running');
        expect(scroll).not.toHaveBeenCalled();

        await focus.show(['chi2-history']);
        expect(focusOf(root, 'chi2-history')).toBe('primary');
        focus.done('step-1');
        expect(focusOf(root, 'bond-angles')).toBeNull();
        expect(focusOf(root, 'chi2-history')).toBe('primary');

        focus.preview('plot:GTS_250K_PDFpartials.csv');
        expect(focusOf(root, 'plot:GTS_250K_PDFpartials.csv')).toBe('preview');
        focus.preview(null);
        expect(focusOf(root, 'plot:GTS_250K_PDFpartials.csv')).toBeNull();
    });

    it('re-applies the focus to a card the host re-renders', async () => {
        const focus = controller({ scroll: () => false });
        await focus.show(['chi2-history']);
        const old = root.querySelector(cardSelector('chi2-history'));
        const fresh = document.createElement('article');
        fresh.setAttribute('data-copilot-card', 'chi2-history');
        old.replaceWith(fresh);
        await new Promise((resolve) => setTimeout(resolve, 40));
        expect(fresh.getAttribute('data-copilot-focus')).toBe('primary');
        focus.dispose();
        expect(fresh.hasAttribute('data-copilot-focus')).toBe(false);
    });

    it('lowers a tall sheet (beforeScroll) before it scrolls', async () => {
        const order = [];
        const focus = controller({
            beforeScroll: async () => { order.push('room'); },
            scroll: () => { order.push('scroll'); return true; }
        });
        await focus.show(['chi2-history']);
        expect(order).toEqual(['room', 'scroll']);
    });
});
