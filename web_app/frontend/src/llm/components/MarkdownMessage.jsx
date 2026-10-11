// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useMemo } from 'react';
import Markdown, { defaultUrlTransform } from 'react-markdown';
import remarkGfm from 'remark-gfm';
import remarkMath from 'remark-math';
import rehypeRaw from 'rehype-raw';
import rehypeSanitize, { defaultSchema } from 'rehype-sanitize';
import rehypeKatex from 'rehype-katex';
import 'katex/dist/katex.min.css';
import { CARD_SCHEME, cardIdFromHref, linkCardReferences } from '../focus/cardRefs';

// Model replies render as Markdown (GitHub-flavored: tables, lists, code) via
// react-markdown, which builds a React tree — there is no innerHTML.
//
// Models often use inline HTML (e.g. <br> for line breaks inside table cells),
// so raw HTML is parsed (rehypeRaw) and then sanitized (rehypeSanitize): safe
// tags like <br>/<sub> render, while scripts, event handlers, and images (no
// external loads from model output) are stripped — the no-injection guarantee
// holds.
//
// Models also emit LaTeX (e.g. $R_{wp}$, $\text{\AA}$) when summarising results,
// often inside table cells. remark-math turns $...$ into placeholder
// <span class="math math-inline"> / <div class="math math-display"> nodes and
// rehype-katex renders those to HTML + MathML. rehype-katex runs *after*
// rehypeSanitize on purpose: KaTeX markup relies on class names, inline `style`
// (glyph heights/offsets) and MathML tags that our locked-down schema strips —
// sanitising KaTeX's output would break the math, and loosening the schema
// enough to keep it would also re-permit `style`/tags on the model's own raw
// HTML (re-opening the external-load vector we close by dropping <img>). Running
// KaTeX last keeps model HTML fully sanitised while KaTeX emits its own trusted
// markup; KaTeX's default trust:false means LaTeX can't smuggle HTML in either.
// So the schema only has to whitelist the math* placeholder classes so the spans
// survive sanitize and reach KaTeX. Order: raw → sanitize → katex.
//
// Card references: the Copilot cites a host card as `[label](card:<id>)` (or
// `[[id]]`, rewritten first by linkCardReferences). The `card:` scheme passes
// react-markdown's URL filter and the sanitizer's protocol list, and the link
// renders as a chip button that focuses the card — only for a card the host
// listed; any other card: target renders as its plain label.
const SANITIZE_SCHEMA = {
    ...defaultSchema,
    protocols: {
        ...defaultSchema.protocols,
        href: [...(defaultSchema.protocols?.href || []), CARD_SCHEME.slice(0, -1)]
    },
    tagNames: (defaultSchema.tagNames || []).filter((tag) => tag !== 'img'),
    attributes: {
        ...defaultSchema.attributes,
        span: [...(defaultSchema.attributes?.span || []), ['className', 'math', 'math-inline', 'math-display']],
        div: [...(defaultSchema.attributes?.div || []), ['className', 'math', 'math-inline', 'math-display']]
    }
};
const REMARK_PLUGINS = [remarkGfm, remarkMath];
const REHYPE_PLUGINS = [
    rehypeRaw,
    [rehypeSanitize, SANITIZE_SCHEMA],
    rehypeKatex
];

// The plain text of a rendered label (strings and nested elements).
const plainText = (node) => {
    if (node == null || typeof node === 'boolean') return '';
    if (typeof node === 'string' || typeof node === 'number') return String(node);
    if (Array.isArray(node)) return node.map(plainText).join('');
    return plainText(node.props?.children);
};

const urlTransform = (url) => (
    typeof url === 'string' && url.trim().toLowerCase().startsWith(CARD_SCHEME) ? url : defaultUrlTransform(url)
);

const TABLE = ({ node, ...props }) => <div className="llm-md-table"><table {...props} /></div>;

// Open links in a new tab and wrap tables so wide ones scroll instead of
// breaking the bubble. `node` is dropped so it is not spread onto the DOM.
// A card: link is a chip that focuses the card (onCardClick), and previews it
// while hovered (onCardHover).
const buildComponents = (cards, onCardClick, onCardHover) => ({
    a: ({ node, href, children, ...props }) => {
        if (typeof href === 'string' && href.toLowerCase().startsWith(CARD_SCHEME)) {
            const id = cardIdFromHref(href, cards);
            if (!id || !onCardClick) return <span className="llm-card-ref">{children}</span>;
            const card = cards.find((entry) => entry.id === id);
            // A label that is only the reference itself (`<card:id>`, an
            // autolink, or `[card:id](card:id)`) reads as the card's title.
            const label = plainText(children).trim();
            const bare = !label || label.toLowerCase().startsWith(CARD_SCHEME) || label === id;
            return (
                <button
                    type="button"
                    className="llm-card-chip"
                    data-card={id}
                    onClick={() => onCardClick(id)}
                    onMouseEnter={() => onCardHover?.(id)}
                    onMouseLeave={() => onCardHover?.(null)}
                    onFocus={() => onCardHover?.(id)}
                    onBlur={() => onCardHover?.(null)}
                    title={`Show ${card?.title || id}`}
                >
                    <span className="llm-card-chip-dot" aria-hidden="true" />
                    {bare ? (card?.title || id) : children}
                </button>
            );
        }
        return <a {...props} href={href} target="_blank" rel="noreferrer">{children}</a>;
    },
    table: TABLE
});

// `inline`: one line of running text (the Outcome note) — paragraphs render
// as spans, so the math and emphasis models put there render in place.
const INLINE_PARAGRAPH = ({ node, ...props }) => <span {...props} />;

const MarkdownMessage = ({ children, cards = null, onCardClick = null, onCardHover = null, inline = false }) => {
    const components = useMemo(() => ({
        ...buildComponents(cards || [], onCardClick, onCardHover),
        ...(inline ? { p: INLINE_PARAGRAPH } : {})
    }), [cards, onCardClick, onCardHover, inline]);
    const text = useMemo(() => linkCardReferences(children, cards || []), [children, cards]);
    return (
        <Markdown
            remarkPlugins={REMARK_PLUGINS}
            rehypePlugins={REHYPE_PLUGINS}
            components={components}
            urlTransform={urlTransform}
        >
            {text}
        </Markdown>
    );
};

export default MarkdownMessage;
