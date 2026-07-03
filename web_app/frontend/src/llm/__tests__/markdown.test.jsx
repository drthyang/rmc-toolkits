import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import Markdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

describe('markdown table rendering (remark-gfm)', () => {
    it('renders a GFM pipe table as a real <table>', () => {
        const md = '| Dataset | Rwp |\n| --- | --- |\n| x-ray S(Q) | 0.046 |\n| xPDF | 0.012 |';
        const html = renderToStaticMarkup(<Markdown remarkPlugins={[remarkGfm]}>{md}</Markdown>);
        expect(html).toContain('<table>');
        expect(html).toContain('<th>Dataset</th>');
        expect(html).toContain('<td>0.046</td>');
    });

    it('does not emit raw HTML from the markdown source', () => {
        const html = renderToStaticMarkup(
            <Markdown remarkPlugins={[remarkGfm]}>{'hi <img src=x onerror=alert(1)>'}</Markdown>
        );
        expect(html).not.toContain('<img');
    });
});
