import React, { useId } from 'react';
import './InfoBadge.css';

/**
 * Small "?" badge that reveals a short explanation on hover or keyboard focus.
 * Accessible: the trigger is a real button labelled by `label`, and the popover
 * is associated via aria-describedby so screen readers announce it too.
 *
 * @param {string}  label   - accessible name for the trigger (e.g. "About …").
 * @param {React.ReactNode} children - the description shown in the popover.
 * @param {'start'|'end'} [align='start'] - horizontal edge the popover aligns to.
 */
const InfoBadge = ({ label = 'More information', children, align = 'start' }) => {
    const id = useId();
    return (
        <span className="info-badge">
            <button
                type="button"
                className="info-badge-trigger"
                aria-label={label}
                aria-describedby={id}
            >
                ?
            </button>
            <span id={id} role="tooltip" className={`info-badge-popover info-badge-popover--${align}`}>
                {children}
            </span>
        </span>
    );
};

export default InfoBadge;
