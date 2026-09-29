// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// uiScale() is the root type size over the 15px design size: code drawing in
// CSS pixels (canvas text, the 'fit' plot viewBox) multiplies by it so it grows
// with the rem-sized UI on 2K / 4K monitors.

// @vitest-environment jsdom

import { afterEach, describe, expect, it, vi } from 'vitest';
import { BASE_ROOT_FONT_PX, uiScale } from '../uiScale';

describe('uiScale', () => {
    afterEach(() => {
        document.documentElement.style.fontSize = '';
    });

    it('is 1 at the 15px design size', () => {
        document.documentElement.style.fontSize = `${BASE_ROOT_FONT_PX}px`;
        expect(uiScale()).toBe(1);
    });

    it('follows the root size on 2K and 4K monitors', () => {
        document.documentElement.style.fontSize = '17px';
        expect(uiScale()).toBeCloseTo(17 / 15, 12);
        document.documentElement.style.fontSize = '21px';
        expect(uiScale()).toBeCloseTo(1.4, 12);
    });

    it('falls back to 1 when the root size cannot be read', () => {
        const spy = vi.spyOn(window, 'getComputedStyle').mockReturnValue({ fontSize: '' });
        try {
            expect(uiScale()).toBe(1);
        } finally {
            spy.mockRestore();
        }
    });
});
