// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

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
        expect(render({ align: 'end', children: 'x' })).toContain('ui-info__popover--end');
        expect(render({ children: 'x' })).toContain('ui-info__popover--start');
    });

    it('opens above on request, below by default', () => {
        expect(render({ side: 'above', align: 'end', children: 'x' }))
            .toContain('class="ui-info__popover ui-info__popover--end ui-info__popover--above"');
        expect(render({ children: 'x' })).toContain('class="ui-info__popover ui-info__popover--start"');
    });
});
