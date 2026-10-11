// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The AI Copilot's view tools on the host side: orient_3d_view turns one of
// the Workbench's 3D views (the panels registered in viewer3d.js) to look
// down a crystal direction, centres, rings or selects sites, or resets it;
// current_view says what the user has on screen. Both are view tools
// (`view: true` in the llm tool contract): they change only what is on
// screen — the page, the camera, a highlight, the selected site as a click
// would — never data, an analysis input or a file. The card itself is
// brought on screen and breathes through the dock's card focus (`focus` in
// the run context), and a site highlight ends with that focus.

import {
    describeDirection, directionLabel, directionVector, getViewer, parseDirection, viewUp, waitForViewer
} from './viewer3d';

export const VIEWER_WAIT_MS = 6000;
// How long a selected site may take to load (its KDE volume, its sphere map).
export const SETTLE_WAIT_MS = 15000;
export const MAX_HIGHLIGHT_SITES = 24;

const AXIS_NAMES = ['a', 'b', 'c'];

const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

const waitUntil = async (test, { timeoutMs, signal, pollMs = 100 }) => {
    const started = Date.now();
    for (;;) {
        if (signal?.aborted) throw new DOMException('Aborted', 'AbortError');
        if (test()) return true;
        if (Date.now() - started >= timeoutMs) return false;
        await sleep(pollMs);
    }
};

const parsedLabel = (parsed) => {
    if (parsed.uvw) return directionLabel(parsed.uvw);
    if (parsed.axis != null) return `${AXIS_NAMES[parsed.axis]} axis`;
    if (parsed.pc != null) return `PC${parsed.pc + 1}`;
    return `Cartesian ${directionLabel(parsed.cartesian)}`;
};

/**
 * The view tools for the host's 3D cards (`cards` entries with `view3d`).
 * `pages` maps page keys to labels; `activePage()` names the page on screen.
 */
export const buildViewTools = ({ cards = [], pages = {}, activePage = () => null, waitMs = VIEWER_WAIT_MS, settleMs = SETTLE_WAIT_MS } = {}) => {
    const views = cards.filter((card) => card.view3d);
    if (!views.length) return [];
    const byId = new Map(cards.map((card) => [card.id, card]));
    const pageLabel = (key) => pages[key] || key;

    const orient = {
        name: 'orient_3d_view',
        view: true,
        description: 'Turn one of the Workbench\'s 3D views so the user sees what you describe: look down a crystal '
            + 'direction [u v w] (or an axis a, b, c; PC1–PC3 for a site\'s principal axes), centre a site, ring sites '
            + 'in blue, select a site as a click would, or reset. The view\'s page opens and its card breathes. '
            + '"Looking down [u v w]" puts that direction toward the viewer, with c (or b) up. Views: '
            + views.map((card) => `${card.id} (${pageLabel(card.page)}: ${card.about || card.title}`
                + `${card.sites ? '; sites' : ''})`).join('; ') + '.',
        parameters: {
            type: 'object',
            properties: {
                card: { type: 'string', enum: views.map((card) => card.id), description: 'Which 3D view.' },
                along: {
                    type: 'string',
                    description: 'Direction to look down: "[111]", "[1 -1 0]", "[0.71 0.71 1]", or "a", "b", "c", "PC1"…'
                },
                center_site: { type: 'integer', minimum: 1, description: 'Site ref to centre (views with sites).' },
                highlight_sites: {
                    type: 'array',
                    items: { type: 'integer', minimum: 1 },
                    maxItems: MAX_HIGHLIGHT_SITES,
                    description: 'Site refs to ring in blue (views with sites).'
                },
                select_site: {
                    type: 'integer',
                    minimum: 1,
                    description: 'Select this site on the page (the direction sphere and the density then show it).'
                },
                reset: { type: 'boolean', description: 'Back to the default view first.' }
            },
            required: ['card']
        },
        label: ({ card, along, select_site: select }) => {
            const title = byId.get(card)?.title || card || '3D view';
            const parts = [title];
            if (select != null) parts.push(`site ${select}`);
            if (along) parts.push(`along ${String(along)}`);
            return parts.join(' · ');
        },
        cards: ({ card }) => (card ? [card] : []),
        run: async (args, { signal, focus } = {}) => {
            const card = byId.get(args.card);
            if (!card?.view3d) throw new Error(`"${args.card}" is not a 3D view; views: ${views.map((view) => view.id).join(', ')}`);
            const parsed = args.along != null ? parseDirection(args.along) : null;
            if (args.along != null && !parsed) {
                throw new Error(`could not read along=${JSON.stringify(args.along)}: give a direction like "[111]", "[1 -1 0]" or an axis "a", "b", "c"`);
            }
            const actions = ['along', 'center_site', 'highlight_sites', 'select_site', 'reset'].filter((key) => args[key] != null);
            if (!actions.length) throw new Error('nothing to do: give along, center_site, highlight_sites, select_site or reset');

            if (focus) await focus.show([card.id]);
            const viewer = await waitForViewer(card.id, { signal, timeoutMs: waitMs });
            if (!viewer) {
                throw new Error(`the ${card.title} view has not drawn yet — its page needs data `
                    + '(a run with a structure; on the site pages, sites to load)');
            }
            const data = { card: card.id, page: pageLabel(card.page) };
            const checks = [];

            if (args.select_site != null) {
                if (!viewer.has('selectSite')) {
                    checks.push({ level: 'note', text: `${card.title} selects no site; select_site was ignored.` });
                } else if (!viewer.selectSite(args.select_site)) {
                    const refs = (viewer.sites?.() || []).map((site) => site.ref);
                    throw new Error(`no site ${args.select_site} here${refs.length ? `; sites: ${refs.join(', ')}` : ''}`);
                } else {
                    data.selected_site = args.select_site;
                    if (viewer.has('settled')) {
                        const settled = await waitUntil(() => viewer.settled(args.select_site), { timeoutMs: settleMs, signal });
                        if (!settled) checks.push({ level: 'warn', text: `Site ${args.select_site} was still loading; the view may change once it has drawn.` });
                    }
                }
            }
            if (args.reset) viewer.reset();

            const cell = viewer.cell?.() || null;
            if (parsed) {
                const direction = directionVector(parsed, cell, viewer.pcaAxes?.() || null);
                if (!direction) {
                    throw new Error(parsed.pc != null
                        ? 'this view has no principal axes (no site selected?)'
                        : 'this view has no unit cell to read a crystal direction in');
                }
                const up = viewUp(direction, cell);
                viewer.lookAlong(direction, up);
                Object.assign(data, { requested: parsedLabel(parsed), ...describeDirection(direction, up, cell) });
            }
            if (args.center_site != null) {
                if (!viewer.has('centerOn')) {
                    checks.push({ level: 'note', text: `${card.title} has no sites to centre; center_site was ignored.` });
                } else if (!viewer.centerOn(args.center_site)) {
                    throw new Error(`no site ${args.center_site} on ${card.title}`);
                } else {
                    data.centred_on = args.center_site;
                }
            }
            if (args.highlight_sites != null) {
                if (!viewer.has('highlightSites')) {
                    checks.push({ level: 'note', text: `${card.title} cannot ring sites; highlight_sites was ignored.` });
                } else {
                    const { found = [], missing = [] } = viewer.highlightSites(args.highlight_sites) || {};
                    // The rings end with the card's focus: the next highlight,
                    // a new question, or the user's click.
                    if (found.length) focus?.onClear(() => viewer.clearHighlight());
                    data.highlighted_sites = found;
                    if (missing.length) checks.push({ level: 'warn', text: `No site ${missing.join(', ')} on ${card.title}.` });
                }
            }

            const did = [];
            if (data.selected_site != null) did.push(`site ${data.selected_site} selected`);
            if (args.reset) did.push('reset');
            if (data.view_along) did.push(`looking down ${data.view_along}${data.screen_up ? ` with ${data.screen_up} up` : ''}`);
            else if (data.view_along_cartesian) did.push(`looking down Cartesian ${directionLabel(data.view_along_cartesian)}`);
            if (data.centred_on != null) did.push(`centred on site ${data.centred_on}`);
            if (data.highlighted_sites?.length) did.push(`ringed site${data.highlighted_sites.length === 1 ? '' : 's'} ${data.highlighted_sites.join(', ')}`);
            if (!checks.some((check) => check.level === 'warn')) {
                checks.push({ level: 'ok', text: `${card.title}: ${did.join(', ') || 'shown'}.` });
            }
            return { data, checks };
        }
    };

    const current = {
        name: 'current_view',
        view: true,
        description: 'What the user has on screen: the page, its cards, the cards highlighted, and for each 3D view '
            + 'on the page the direction it looks down and the selected site. Use it when the user says "this", '
            + '"here" or "what I\'m looking at".',
        parameters: { type: 'object', properties: {} },
        label: () => 'What is on screen',
        run: async (args, { focus } = {}) => {
            const key = activePage();
            const onPage = cards.filter((card) => card.page === key);
            const data = {
                page: pageLabel(key),
                cards: onPage.map((card) => card.id),
                highlighted: focus?.state?.().cards || []
            };
            const viewsOnPage = onPage.filter((card) => card.view3d).map((card) => {
                const viewer = getViewer(card.id);
                if (!viewer?.ready()) return { card: card.id, drawn: false };
                return { card: card.id, ...(viewer.describe() || {}) };
            });
            if (viewsOnPage.length) data.views_3d = viewsOnPage;
            return { data, checks: [{ level: 'ok', text: `On the ${data.page} page.` }] };
        }
    };

    return [orient, current];
};
