// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useId } from 'react';

/**
 * Small "?" badge that reveals a short explanation on hover or keyboard focus.
 * Accessible: the trigger is a real button labelled by `label`, and the popover
 * is associated via aria-describedby so screen readers announce it too.
 * Styled by the kit's `ui-info` classes (ui.css). Takes only the props below —
 * no `className` and no `...rest` pass-through, unlike the other kit components.
 *
 * @param {string}  label   - accessible name for the trigger (e.g. "About …").
 * @param {React.ReactNode} children - the description shown in the popover.
 * @param {'start'|'end'} [align='start'] - horizontal edge the popover aligns to.
 */
const InfoBadge = ({ label = 'More information', children, align = 'start' }) => {
    const id = useId();
    return (
        <span className="ui-info">
            <button
                type="button"
                className="ui-info__trigger"
                aria-label={label}
                aria-describedby={id}
            >
                ?
            </button>
            <span id={id} role="tooltip" className={`ui-info__popover ui-info__popover--${align}`}>
                {children}
            </span>
        </span>
    );
};

export default InfoBadge;
