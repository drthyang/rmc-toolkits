// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useId, useState } from 'react';
import { measurePopoverEdge } from './popoverEdge';

/**
 * Small "?" badge that reveals a short explanation on hover or keyboard focus.
 * Accessible: the trigger is a real button labelled by `label`, and the popover
 * is associated via aria-describedby so screen readers announce it too.
 * Styled by the kit's `ui-info` classes (ui.css). Takes only the props below —
 * no `className` and no `...rest` pass-through, unlike the other kit components.
 *
 * A closed popover takes no space (display: none), so a badge near the right
 * edge never makes its page scroll sideways; an open one aligns to the other
 * edge when the requested one would run off the window (popoverEdge.js).
 *
 * @param {string}  label   - accessible name for the trigger (e.g. "About …").
 * @param {React.ReactNode} children - the description shown in the popover.
 * @param {'start'|'end'} [align='start'] - horizontal edge the popover aligns to
 *        when it fits the window.
 * @param {'below'|'above'} [side='below'] - 'above' for a badge near the
 *        bottom of the page, so the open popover stays on the page.
 */
const InfoBadge = ({ label = 'More information', children, align = 'start', side = 'below' }) => {
    const id = useId();
    const [edge, setEdge] = useState(null);
    const place = (event) => setEdge(measurePopoverEdge(event.currentTarget, align));
    return (
        <span className="ui-info" onPointerEnter={place} onFocus={place}>
            <button
                type="button"
                className="ui-info__trigger"
                aria-label={label}
                aria-describedby={id}
            >
                ?
            </button>
            <span
                id={id}
                role="tooltip"
                className={`ui-info__popover ui-info__popover--${edge || align}${side === 'above' ? ' ui-info__popover--above' : ''}`}
            >
                {children}
            </span>
        </span>
    );
};

export default InfoBadge;
