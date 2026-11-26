// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang
/* @vitest-environment jsdom */

// The Three.js panels draw nothing while their page is hidden
// (pageVisibility.js): App keeps a visited page mounted with `is-hidden`.

import { describe, expect, it } from 'vitest';
import { isDisplayed } from '../pageVisibility';

describe('isDisplayed', () => {
    it('is false inside a hidden workspace page and true on the page in view', () => {
        document.body.innerHTML = `
            <div class="workspace-page is-active"><div><canvas id="shown"></canvas></div></div>
            <div class="workspace-page is-hidden"><div><canvas id="hidden"></canvas></div></div>
            <canvas id="loose"></canvas>`;
        expect(isDisplayed(document.getElementById('shown'))).toBe(true);
        expect(isDisplayed(document.getElementById('hidden'))).toBe(false);
        expect(isDisplayed(document.getElementById('loose'))).toBe(true);
        expect(isDisplayed(null)).toBe(true);
    });
});
