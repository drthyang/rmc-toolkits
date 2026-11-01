// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Whether a card that asks the user something (a question, a change to
// approve) may take the focus when it appears: from nowhere in particular or
// from inside the conversation (the empty composer, where it rests after a
// question is sent), never out of a field that holds text the user is typing.
export const mayTakeFocus = (card) => {
    if (typeof document === 'undefined' || !card) return false;
    const active = document.activeElement;
    if (!active || active === document.body) return true;
    if (active.matches?.('input, textarea, select, [contenteditable="true"]') && (active.value || active.isContentEditable)) {
        return false;
    }
    return Boolean(card.closest('.llm-chat-panel')?.contains(active));
};
