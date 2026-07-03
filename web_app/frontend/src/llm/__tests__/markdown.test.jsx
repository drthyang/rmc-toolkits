import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import Markdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import rehypeRaw from 'rehype-raw';
import rehypeSanitize, { defaultSchema } from 'rehype-sanitize';

// Same plugin pipeline ChatView uses: GFM + raw HTML parsed then sanitized
// (images dropped so model output can't load external resources).
const remarkPlugins = [remarkGfm];
const rehypePlugins = [
    rehypeRaw,
    [rehypeSanitize, { ...defaultSchema, tagNames: (defaultSchema.tagNames || []).filter((tag) => tag !== 'img') }]
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
});
