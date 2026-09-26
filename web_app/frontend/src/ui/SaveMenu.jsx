// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useEffect, useRef, useState } from 'react';
import cx from './cx';

const DEFAULT_OPTIONS = [{ id: 'png', label: 'PNG image', hint: '.png' }];

// A compact "save" badge. With several formats it opens a small menu; with a
// single format it saves directly on click. Shared by the chart toolbars and
// the KDE panels so every figure offers the same control. Styled by the kit's
// `ui-save` / `ui-menu` classes (ui.css); pass className="ui-save--accent" for
// the accent look.
const SaveMenu = ({ onSave, options = DEFAULT_OPTIONS, label = 'Save', align = 'right', disabled = false, busy = false, className = '' }) => {
    const [open, setOpen] = useState(false);
    const rootRef = useRef(null);

    useEffect(() => {
        if (!open) return undefined;
        const handlePointer = (event) => {
            if (rootRef.current && !rootRef.current.contains(event.target)) setOpen(false);
        };
        const handleKey = (event) => {
            if (event.key === 'Escape') setOpen(false);
        };
        document.addEventListener('pointerdown', handlePointer);
        document.addEventListener('keydown', handleKey);
        return () => {
            document.removeEventListener('pointerdown', handlePointer);
            document.removeEventListener('keydown', handleKey);
        };
    }, [open]);

    const choose = (id) => {
        setOpen(false);
        onSave(id);
    };

    const handleTrigger = () => {
        if (options.length === 1) choose(options[0].id);
        else setOpen((value) => !value);
    };

    const multiple = options.length > 1;

    return (
        <div className={cx('ui-save', className)} ref={rootRef}>
            <button
                type="button"
                className="ui-save__trigger"
                onClick={handleTrigger}
                disabled={disabled || busy}
                aria-haspopup={multiple ? 'menu' : undefined}
                aria-expanded={multiple ? open : undefined}
                title="Save figure"
            >
                <span className="ui-save__icon" aria-hidden="true">⤓</span>
                {busy ? 'Saving…' : label}
            </button>
            {open && multiple && (
                <div className={`ui-menu ui-menu--${align}`} role="menu">
                    {options.map((option) => (
                        <button
                            key={option.id}
                            type="button"
                            role="menuitem"
                            className="ui-menu__item"
                            onClick={() => choose(option.id)}
                        >
                            {option.label}
                            {option.hint && <span>{option.hint}</span>}
                        </button>
                    ))}
                </div>
            )}
        </div>
    );
};

export default SaveMenu;
