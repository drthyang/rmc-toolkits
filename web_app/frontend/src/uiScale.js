// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The root type size in index.css is 15px up to a 1080p-class window and
// grows on 2K / 4K monitors (17px at 1440p, 21px at 2160p). Everything sized
// in rem follows it on its own; code that draws in CSS pixels — canvas
// overlay text, the 'fit' plot whose viewBox is one user unit per pixel —
// multiplies its sizes by uiScale() so it grows in step with the chrome.

export const BASE_ROOT_FONT_PX = 15;

// Current root font size over the 15px design size (1 on laptops, tablets and
// phones; 1 when there is no document to measure).
export const uiScale = () => {
    if (typeof window === 'undefined' || typeof document === 'undefined') return 1;
    const size = parseFloat(window.getComputedStyle(document.documentElement).fontSize);
    return Number.isFinite(size) && size > 0 ? size / BASE_ROOT_FONT_PX : 1;
};
