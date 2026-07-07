import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import InfoBadge from '../InfoBadge';

const render = (props) => renderToStaticMarkup(<InfoBadge {...props} />);

describe('InfoBadge', () => {
    it('renders an accessible trigger described by the popover', () => {
        const html = render({ label: 'How it works', children: <p>Detail text.</p> });
        expect(html).toContain('aria-label="How it works"');
        expect(html).toContain('role="tooltip"');
        expect(html).toContain('Detail text.');
        // Trigger's aria-describedby points at the popover's id.
        const describedBy = html.match(/aria-describedby="([^"]+)"/)?.[1];
        expect(describedBy).toBeTruthy();
        expect(html).toContain(`id="${describedBy}"`);
    });

    it('aligns the popover to the requested edge', () => {
        expect(render({ align: 'end', children: 'x' })).toContain('info-badge-popover--end');
        expect(render({ children: 'x' })).toContain('info-badge-popover--start');
    });
});
