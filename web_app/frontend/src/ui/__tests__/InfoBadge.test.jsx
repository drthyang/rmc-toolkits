// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { renderToStaticMarkup } from 'react-dom/server';
import InfoBadge from '../InfoBadge';
import { popoverEdge } from '../popoverEdge';

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

describe('popover edge', () => {
    const at = (left, align = 'start') => popoverEdge({
        align, left, right: left + 16, viewportWidth: 820, popoverWidth: 300, gap: 11
    });

    it('keeps the requested edge when the popover fits beside it', () => {
        expect(at(100)).toBe('start');
        expect(at(600, 'end')).toBe('end');
    });

    it('turns to the other edge when the requested one runs off the window', () => {
        // Displacement Directions at 820 px: "About smoothing" sat at x = 532,
        // and its start-aligned 300 px popover reached 832.
        expect(at(532)).toBe('end');
        expect(at(100, 'end')).toBe('start');
    });

    it('keeps the requested edge when neither fits', () => {
        const narrow = { left: 150, right: 166, viewportWidth: 320, popoverWidth: 300, gap: 11 };
        expect(popoverEdge({ ...narrow, align: 'start' })).toBe('start');
        expect(popoverEdge({ ...narrow, align: 'end' })).toBe('end');
    });
});
