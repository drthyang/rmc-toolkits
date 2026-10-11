// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import MarkdownMessage from '../components/MarkdownMessage';

// The Copilot's reply renderer itself: GFM + math on the remark side; on the
// rehype side raw HTML is parsed, then sanitized (images dropped so model
// output can't load external resources), then KaTeX renders $...$ into trusted
// markup. Order matters: raw → sanitize → katex (see MarkdownMessage.jsx).
const render = (md) => renderToStaticMarkup(<MarkdownMessage>{md}</MarkdownMessage>);

describe('Copilot markdown rendering', () => {
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

describe('card references', () => {
    const cards = [{ id: 'chi2-history', title: 'χ² history' }, { id: 'plot:a_PDFpartials.csv', title: 'g(r) partials' }];
    const withCards = (md) => renderToStaticMarkup(
        <MarkdownMessage cards={cards} onCardClick={() => {}}>{md}</MarkdownMessage>
    );

    it('renders a cited card as a chip button, never as a link', () => {
        const html = withCards('The [χ² history](card:chi2-history) and [[plot:a_PDFpartials.csv]].');
        expect(html).toContain('<button type="button" class="llm-card-chip" data-card="chi2-history" title="Show χ² history">');
        expect(html).toContain('data-card="plot:a_PDFpartials.csv"');
        expect(html).toContain('g(r) partials</button>');
        expect(html).not.toContain('href="card:');
    });

    it('labels a chip written as the bare reference with the card\'s title', () => {
        // gemma4:26b wrote "<card:direction-sphere>"-style references in the
        // real-model check; the chip read "card:direction-sphere".
        const html = withCards('See <card:chi2-history>, [card:chi2-history](card:chi2-history) and [chi2-history](card:chi2-history).');
        expect(html.match(/class="llm-card-chip"/g)).toHaveLength(3);
        expect(html).not.toContain('>card:chi2-history<');
        expect(html.match(/χ² history<\/button>/g)).toHaveLength(3);
        // A real label is kept.
        expect(withCards('[the last run](card:chi2-history)')).toContain('the last run</button>');
    });

    it('shows an unknown card as its label, and leaves ordinary links alone', () => {
        const html = withCards('[gone](card:nothing) and [docs](https://example.org)');
        expect(html).toContain('<span class="llm-card-ref">gone</span>');
        expect(html).toContain('<a href="https://example.org" target="_blank" rel="noreferrer">docs</a>');
    });

    it('without cards, a card: link is plain text and javascript: stays blocked', () => {
        const html = render('[the history](card:chi2-history) <a href="javascript:alert(1)">x</a>');
        expect(html).toContain('<span class="llm-card-ref">the history</span>');
        expect(html).not.toContain('javascript:');
    });
});

describe('inline notes', () => {
    it('renders math and emphasis in place, without a paragraph', () => {
        const html = renderToStaticMarkup(<MarkdownMessage inline>{'The site with the largest $U_{iso}$ is **23**.'}</MarkdownMessage>);
        expect(html.startsWith('<span>')).toBe(true);
        expect(html).not.toContain('<p>');
        expect(html).toContain('class="katex"');
        expect(html).not.toContain('$U_{iso}$');
        expect(html).toContain('<strong>23</strong>');
    });
});
