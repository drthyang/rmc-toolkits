// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import Markdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import remarkMath from 'remark-math';
import rehypeRaw from 'rehype-raw';
import rehypeSanitize, { defaultSchema } from 'rehype-sanitize';
import rehypeKatex from 'rehype-katex';

// Same plugin pipeline ChatView uses: GFM + math on the remark side; on the
// rehype side raw HTML is parsed, then sanitized (images dropped so model output
// can't load external resources), then KaTeX renders $...$ into trusted markup.
// Order matters: raw → sanitize → katex (see ChatView.jsx for the rationale).
const sanitizeSchema = {
    ...defaultSchema,
    tagNames: (defaultSchema.tagNames || []).filter((tag) => tag !== 'img'),
    attributes: {
        ...defaultSchema.attributes,
        span: [...(defaultSchema.attributes?.span || []), ['className', 'math', 'math-inline', 'math-display']],
        div: [...(defaultSchema.attributes?.div || []), ['className', 'math', 'math-inline', 'math-display']]
    }
};
const remarkPlugins = [remarkGfm, remarkMath];
const rehypePlugins = [
    rehypeRaw,
    [rehypeSanitize, sanitizeSchema],
    rehypeKatex
];

const render = (md) => renderToStaticMarkup(
    <Markdown remarkPlugins={remarkPlugins} rehypePlugins={rehypePlugins}>{md}</Markdown>
);

describe('assistant markdown rendering', () => {
    it('renders a GFM pipe table as a real <table>', () => {
        const html = render('| Dataset | Rwp |\n| --- | --- |\n| x-ray S(Q) | 0.046 |\n| xPDF | 0.012 |');
        expect(html).toContain('<table>');
        expect(html).toContain('<th>Dataset</th>');
        expect(html).toContain('<td>0.046</td>');
    });

    it('renders inline <br> line breaks (as models emit in table cells)', () => {
        const html = render('Space group: F-43m<br>Point group: -43m');
        expect(html).toContain('<br');
        expect(html).toContain('Point group');
    });

    it('strips dangerous HTML — no scripts, event handlers, or images', () => {
        const html = render('ok <img src=x onerror=alert(1)> <script>alert(2)</script> <a href="javascript:alert(3)">x</a>');
        expect(html).not.toContain('<img');
        expect(html).not.toContain('<script');
        expect(html).not.toContain('onerror');
        expect(html).not.toContain('javascript:');
    });

    it('renders inline LaTeX ($...$) as KaTeX, not literal dollar text', () => {
        const html = render('The weighted profile factor is $R_{wp}$ for the fit.');
        expect(html).toContain('class="katex"');
        // KaTeX keeps the source in a MathML <annotation> — the raw $…$ is gone.
        expect(html).toContain('R_{wp}');
        expect(html).not.toContain('$R_{wp}$');
    });

    it('renders inline math inside GFM table cells', () => {
        const html = render('| Metric | Value |\n| --- | --- |\n| $R_{wp}$ | 0.046 |\n| $r$ (Å) | 1.98 |');
        expect(html).toContain('<table>');
        expect(html).toContain('class="katex"');
        expect(html).toContain('<td>0.046</td>');
    });

    it('renders display math (block $$...$$) in a scrollable block', () => {
        const html = render('$$\nR_{wp} = \\sqrt{\\frac{\\sum w_i (y_i - y_{c,i})^2}{\\sum w_i y_i^2}}\n$$');
        expect(html).toContain('katex-display');
    });

    it('keeps KaTeX inline styles that survive because katex runs after sanitize', () => {
        // If sanitize ran last it would strip KaTeX's style attrs and break layout.
        const html = render('$\\frac{a}{b}$');
        expect(html).toContain('class="katex"');
        expect(html).toContain('style="');
    });
});
