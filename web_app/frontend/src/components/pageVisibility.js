// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The app keeps a visited page mounted and hides it (App.jsx: the
// workspace page's `is-hidden`, display: none), so its Three.js panels keep
// their scenes. Their render loops ask this each frame and draw nothing
// while their page is hidden: five hidden panels otherwise kept ~10 % of
// the main thread busy on the page in view. A DOM walk, no layout read.

const HIDDEN_PAGE = '.workspace-page.is-hidden';

/** Whether `element` sits on the page in view (not in a hidden workspace page). */
export const isDisplayed = (element) => !element?.closest?.(HIDDEN_PAGE);
