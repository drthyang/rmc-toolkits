// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';
import Markdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import remarkMath from 'remark-math';
import rehypeRaw from 'rehype-raw';
import rehypeSanitize, { defaultSchema } from 'rehype-sanitize';
import rehypeKatex from 'rehype-katex';
import 'katex/dist/katex.min.css';

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
const SANITIZE_SCHEMA = {
    ...defaultSchema,
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

// Open links in a new tab and wrap tables so wide ones scroll instead of
// breaking the bubble. `node` is dropped so it is not spread onto the DOM.
const MARKDOWN_COMPONENTS = {
    a: ({ node, ...props }) => <a {...props} target="_blank" rel="noreferrer" />,
    table: ({ node, ...props }) => <div className="llm-md-table"><table {...props} /></div>
};

const MarkdownMessage = ({ children }) => (
    <Markdown remarkPlugins={REMARK_PLUGINS} rehypePlugins={REHYPE_PLUGINS} components={MARKDOWN_COMPONENTS}>
        {children}
    </Markdown>
);

export default MarkdownMessage;
