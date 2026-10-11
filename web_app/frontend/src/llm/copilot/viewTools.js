// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// highlight_cards: the Copilot points at the host's cards while it explains
// them. Built from the host's `cards` prop (id, page, title, about, and
// `marks: true` for a plot whose x axis can carry marks); run() hands the ids
// to the card-focus controller (focus/focus.js), which switches page, scrolls
// the main card into view and lets it breathe until the agent moves on or the
// user clicks. A view tool: it changes what is on screen, never data, a page's
// inputs or a file.

import { resolveCardId } from '../focus/cardRefs';

export const MAX_FOCUS_CARDS = 4;
export const MAX_MARKS = 8;

const pageName = (pages, key) => (key && pages?.[key]) || key || 'this page';

// One line per card for the model: id — title (page): what it shows.
export const describeCards = (cards, pages) => cards
    .map((card) => `${card.id} — ${card.title}${card.page ? ` (${pageName(pages, card.page)})` : ''}`
        + `${card.about ? `: ${card.about}` : ''}`)
    .join('; ');

export const buildFocusTools = ({ cards, pages } = {}) => {
    if (!Array.isArray(cards) || !cards.length) return [];
    const byId = new Map(cards.map((card) => [card.id, card]));
    const titled = (ids) => ids.map((id) => byId.get(id)?.title || id);
    return [{
        name: 'highlight_cards',
        view: true,
        description: 'Point at the Workbench cards you are talking about: the first card is brought on screen '
            + '(its page opens, it scrolls into view) and the cards get a breathing blue edge until you highlight '
            + 'something else or the user clicks. Call it before explaining what a card shows; [] clears. On a '
            + 'plot card (marks: yes) you can also mark x positions, e.g. the g(r) peaks you name (r in Å). '
            + `Cards: ${describeCards(cards, pages)}.`,
        parameters: {
            type: 'object',
            properties: {
                cards: {
                    type: 'array',
                    items: { type: 'string' },
                    maxItems: MAX_FOCUS_CARDS,
                    description: 'Card ids, the main one first.'
                },
                marks: {
                    type: 'array',
                    maxItems: MAX_MARKS,
                    items: {
                        type: 'object',
                        properties: {
                            x: { type: 'number', description: 'Position on the plot\'s x axis (r in Å, Q in Å⁻¹, …).' },
                            label: { type: 'string', description: 'A short label, e.g. "Ta–Se 2.52 Å".' }
                        },
                        required: ['x']
                    },
                    description: 'Optional marks on the first card when it is a plot.'
                }
            },
            required: ['cards']
        },
        label: ({ cards: ids }) => {
            const list = Array.isArray(ids) ? ids.map((id) => resolveCardId(String(id), cards) || id) : [];
            return list.length ? `Pointing at ${titled(list).join(', ')}` : 'Clear the highlight';
        },
        cards: () => [],
        run: async ({ cards: raw = [], marks = [] }, { focus } = {}) => {
            if (!focus) throw new Error('the page view is not available here');
            const unknown = [];
            const ids = [];
            raw.forEach((value) => {
                const id = resolveCardId(String(value), cards);
                if (!id) unknown.push(value);
                else if (!ids.includes(id)) ids.push(id);
            });
            if (unknown.length) {
                throw new Error(`unknown card ${unknown.map((value) => JSON.stringify(value)).join(', ')}; `
                    + `the cards are: ${cards.map((card) => card.id).join(', ')}`);
            }
            const checks = [];
            const main = byId.get(ids[0]);
            let kept = [];
            if (marks.length) {
                if (!main?.marks) {
                    checks.push({ level: 'note', text: `Marks are drawn on plot cards only; ${main ? main.title : 'no card'} has none.` });
                } else {
                    kept = marks.filter((mark) => Number.isFinite(mark.x))
                        .map((mark) => ({ x: mark.x, ...(mark.label ? { label: String(mark.label).slice(0, 40) } : {}) }));
                }
            }
            const result = await focus.show(ids, { marks: kept.length ? kept : null });
            if (!ids.length) {
                return { data: { highlighted: [] }, checks: [{ level: 'ok', text: 'Highlight cleared.' }] };
            }
            const page = main?.page ? pageName(pages, main.page) : null;
            const data = {
                highlighted: ids.map((id) => ({ id, title: byId.get(id)?.title || id, page: pageName(pages, byId.get(id)?.page) })),
                ...(page ? { page } : {}),
                ...(result?.switched ? { switched_page: true } : {}),
                ...(kept.length ? { marks: kept } : {})
            };
            if (result?.superseded) {
                checks.push({ level: 'note', text: 'Another highlight replaced this one before it was shown.' });
            } else if (result?.missing?.length) {
                const missing = titled(result.missing);
                data.not_on_screen = result.missing;
                checks.push({
                    level: result.missing.includes(ids[0]) ? 'warn' : 'note',
                    text: `Not on screen: ${missing.join(', ')} — its page has not drawn it (no data yet?) or it is on another page.`
                });
            }
            if (!checks.some((check) => check.level === 'warn')) {
                checks.push({
                    level: 'ok',
                    text: `Pointing at ${titled(ids).join(', ')}${page ? ` on the ${page} page` : ''}`
                        + `${kept.length ? `, ${kept.length} mark${kept.length === 1 ? '' : 's'}` : ''}.`
                });
            }
            return { data, checks };
        }
    }];
};
