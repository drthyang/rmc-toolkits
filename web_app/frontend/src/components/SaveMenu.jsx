// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useEffect, useRef, useState } from 'react';
import './SaveMenu.css';

const DEFAULT_OPTIONS = [{ id: 'png', label: 'PNG image', hint: '.png' }];

// A compact "save" badge. With several formats it opens a small menu; with a
// single format it saves directly on click. Shared by the chart toolbars and
// the KDE panels so every figure offers the same control.
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
        <div className={`save-menu ${className}`} ref={rootRef}>
            <button
                type="button"
                className="save-menu-trigger"
                onClick={handleTrigger}
                disabled={disabled || busy}
                aria-haspopup={multiple ? 'menu' : undefined}
                aria-expanded={multiple ? open : undefined}
                title="Save figure"
            >
                <span className="save-menu-icon" aria-hidden="true">⤓</span>
                {busy ? 'Saving…' : label}
            </button>
            {open && multiple && (
                <div className={`save-menu-list save-menu-list--${align}`} role="menu">
                    {options.map((option) => (
                        <button
                            key={option.id}
                            type="button"
                            role="menuitem"
                            className="save-menu-item"
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
