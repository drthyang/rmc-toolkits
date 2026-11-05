// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Bond Geometry hero card's tabs: the angle distribution (sin-corrected,
// the default, or density) and the partial g(r) window helper, which used to
// be a card of its own. The Copilot's card ids stay as they were
// (copilotCards.js): 'bond-angles' is the card, 'bond-windows' the g(r) tab
// panel inside it. A highlight on either brings its tab forward, so a
// highlight always lands on what it names.

const pointed = (level) => level === 'primary' || level === 'secondary';

/**
 * Which tab a Copilot highlight asks for, from the focus levels
 * (data-copilot-focus) on the card ('bond-angles') and on its g(r) panel
 * ('bond-windows'). Returns true for the g(r) tab, false for the angle view,
 * null for no change. Only the agent's or the user's pointing (primary,
 * secondary) switches; a running analysis or a hovered chip does not. When
 * both are named, the primary one wins.
 */
export const grTabForFocus = ({ card = null, panel = null } = {}) => {
    if (panel === 'primary') return true;
    if (card === 'primary') return false;
    if (pointed(panel) && !pointed(card)) return true;
    if (pointed(card) && !pointed(panel)) return false;
    return null;
};
