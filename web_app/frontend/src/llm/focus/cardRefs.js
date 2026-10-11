// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Card references in the Copilot's answers. The model cites a card as a
// Markdown link with the card: scheme, `[the g(r) partials](card:plot:GTS_PDFpartials.csv)`,
// or in the short form `[[card:chi2-history]]` / `[[chi2-history]]`; the
// answer renders each known card as a chip that focuses it (MarkdownMessage).
// Unknown ids render as their plain label: a reference never becomes a link
// to anything but a card the host listed. Pure functions.

export const CARD_SCHEME = 'card:';

// A model's spelling of an id → the catalog's id, or null. Accepts the
// scheme, surrounding brackets or backticks, any letter case, and the card's
// title ("χ² history").
export const resolveCardId = (raw, cards) => {
    if (typeof raw !== 'string' || !Array.isArray(cards) || !cards.length) return null;
    const unscheme = (text) => (text.toLowerCase().startsWith(CARD_SCHEME) ? text.slice(CARD_SCHEME.length) : text).trim();
    const plain = unscheme(raw.trim());
    // As written first (a title may end in a parenthesis), then without the
    // brackets, backticks or quotes a model wraps it in.
    const variants = [plain, unscheme(raw.trim().replace(/^[[(`'"<]+|[\])`'">]+$/g, '').trim())]
        .filter((text, index, all) => text && all.indexOf(text) === index);
    for (const text of variants) {
        const lower = text.toLowerCase();
        const card = cards.find((entry) => entry.id === text)
            || cards.find((entry) => entry.id.toLowerCase() === lower)
            || cards.find((entry) => typeof entry.title === 'string' && entry.title.toLowerCase() === lower)
            // A plot card named by its file alone ("GTS_PDFpartials.csv").
            || cards.find((entry) => entry.id.toLowerCase() === `plot:${lower}`);
        if (card) return card.id;
    }
    return null;
};

/** The card id of a `card:` link target, or null for any other href. */
export const cardIdFromHref = (href, cards) => {
    if (typeof href !== 'string') return null;
    const text = decodeURIComponentSafe(href.trim());
    if (!text.toLowerCase().startsWith(CARD_SCHEME)) return null;
    return resolveCardId(text, cards);
};

function decodeURIComponentSafe(text) {
    try {
        return decodeURIComponent(text);
    } catch {
        return text;
    }
}

// Markdown link text must not close early or open a nested link.
const linkText = (text) => String(text).replace(/[[\]]/g, '');

/**
 * Rewrite the short forms `[[card:id]]` and `[[id]]` (known ids only) as
 * `[title](card:id)` links, and give a `card:` link target spaces-safe
 * angle brackets, so Markdown parses every reference as one link. Text inside
 * code spans and fences is left alone.
 */
export const linkCardReferences = (text, cards) => {
    if (typeof text !== 'string' || !text || !Array.isArray(cards) || !cards.length) return text;
    const parts = text.split(/(```[\s\S]*?```|`[^`\n]*`)/g);
    return parts.map((part, index) => {
        if (index % 2 === 1) return part;
        return part
            .replace(/\[\[([^[\]\n]{1,160})\]\]/g, (whole, inner) => {
                const id = resolveCardId(inner, cards);
                if (!id) return whole;
                const card = cards.find((entry) => entry.id === id);
                return `[${linkText(card.title || id)}](<${CARD_SCHEME}${id}>)`;
            })
            // [label](card:some id with spaces) → [label](<card:...>)
            .replace(/\]\((card:[^)<>\n]+)\)/gi, (whole, target) => (/\s/.test(target.trim()) ? `](<${target.trim()}>)` : whole));
    }).join('');
};

/** The ids an answer cites, in order of first appearance (known cards only). */
export const citedCards = (text, cards) => {
    if (typeof text !== 'string' || !text) return [];
    const linked = linkCardReferences(text, cards);
    const ids = [];
    const pattern = /\]\(<?(card:[^)>\n]+)>?\)/gi;
    let match = pattern.exec(linked);
    while (match) {
        const id = cardIdFromHref(match[1], cards);
        if (id && !ids.includes(id)) ids.push(id);
        match = pattern.exec(linked);
    }
    return ids;
};
