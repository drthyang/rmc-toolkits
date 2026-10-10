// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';

// The Copilot's inline SVG icons (currentColor strokes, aria-hidden).

export const ArrowUpIcon = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.2" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <path d="M12 19V5M6 11l6-6 6 6" />
    </svg>
);

export const StopIcon = () => (
    <svg viewBox="0 0 24 24" fill="currentColor" aria-hidden="true">
        <rect x="7" y="7" width="10" height="10" rx="2" />
    </svg>
);

export const UserGlyph = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.9" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <circle cx="12" cy="8" r="3.5" />
        <path d="M5 20c0-3.6 3.1-5.5 7-5.5s7 1.9 7 5.5" />
    </svg>
);

// The app's pair-distribution-wave mark, shared with the header and favicon.
export const WaveMark = () => (
    <svg viewBox="0 0 100 100" fill="none" aria-hidden="true">
        <path
            d="M14 62 C22 62 27 24 37 24 C47 24 45 66 54 66 C62 66 61 44 69 44 C76 44 78 57 86 57"
            transform="translate(0,5)"
            stroke="currentColor"
            strokeWidth="9"
            strokeLinecap="round"
            strokeLinejoin="round"
        />
    </svg>
);

export const ChevronIcon = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <path d="M6 9l6 6 6-6" />
    </svg>
);

export const CloseIcon = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <path d="M6 6l12 12M18 6L6 18" />
    </svg>
);

export const ExpandIcon = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.9" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <path d="M14 4h6v6M10 20H4v-6M20 4l-7 7M4 20l7-7" />
    </svg>
);

export const ShrinkIcon = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.9" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <path d="M20 10h-6V4M4 14h6v6M14 10l7-7M10 14l-7 7" />
    </svg>
);
