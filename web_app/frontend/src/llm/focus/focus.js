// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The cards the Copilot points at. The host marks each card it offers with
// `data-copilot-card="<id>"` and lists them in the `cards` prop (id, page,
// title); this controller gives the cards in focus a breathing blue edge
// (`data-copilot-focus` = primary | secondary | running | preview, styled in
// Copilot.css), switches to the card's page when it is not on screen
// (through the host's openPage), waits for it to mount, and scrolls it into
// view. Its lifecycle follows the conversation:
//
//   - show(ids)        the agent (highlight_cards, a 3D view tool) or the user
//                      (a card chip in an answer) points at cards; replaces
//                      the previous focus.
//   - running(step)    the card of an analysis that is running breathes
//                      softly, on the page on screen only (no scroll, no page
//                      switch), until the step ends.
//   - clear(reason)    a new question, the user clicking anywhere outside the
//                      Copilot (a click, not a touch-scroll toward the card),
//                      the Copilot collapsing, or another run.
//
// `onClear(fn)` lets a tool undo what it drew outside the cards (a 3D view's
// site highlight) when its focus ends. Marks (x positions inside a plot card)
// travel to the host through `onMarks`. No React: the dock owns one instance.

export const CARD_ATTR = 'data-copilot-card';
export const FOCUS_ATTR = 'data-copilot-focus';
// How long show() waits for a card to mount after a page switch (a page
// visited for the first time mounts, then loads its data).
export const FOCUS_WAIT_MS = 3000;
const POLL_MS = 60;

const quoteAttr = (value) => String(value).replace(/["\\]/g, '\\$&');
export const cardSelector = (id) => `[${CARD_ATTR}="${quoteAttr(id)}"]`;

// A card is on screen unless it sits in a hidden subtree: the host hides its
// inactive pages with aria-hidden (App.jsx), and jsdom has no layout, so this
// is the test that works in both.
export const isCardShown = (element) => Boolean(
    element?.isConnected && !element.closest('[hidden], [aria-hidden="true"], [inert]')
);

export const findCard = (root, id) => {
    if (!root || id == null) return null;
    const all = [...root.querySelectorAll(cardSelector(id))];
    return all.find(isCardShown) || null;
};

const prefersReducedMotion = () => typeof window !== 'undefined'
    && typeof window.matchMedia === 'function'
    && window.matchMedia('(prefers-reduced-motion: reduce)').matches;

// The nearest ancestor that scrolls vertically (the page, or the phone-size
// shell), else the document's scrolling element.
const scrollParent = (element) => {
    for (let node = element.parentElement; node; node = node.parentElement) {
        const style = window.getComputedStyle(node);
        if (/(auto|scroll|overlay)/.test(style.overflowY) && node.scrollHeight > node.clientHeight) return node;
    }
    return document.scrollingElement || document.documentElement;
};

/**
 * Bring a card into view in its own scroll container: nothing when it is
 * fully visible, its top when it is taller than the container, else centred.
 * Smooth unless the user prefers reduced motion.
 */
export const scrollCardIntoView = (element) => {
    if (!element || typeof element.scrollIntoView !== 'function') return false;
    const box = element.getBoundingClientRect();
    const container = scrollParent(element);
    const view = container === document.scrollingElement || container === document.documentElement
        ? { top: 0, bottom: window.innerHeight }
        : container.getBoundingClientRect();
    if (box.height > 0 && box.top >= view.top && box.bottom <= view.bottom) return false;
    const tall = box.height > (view.bottom - view.top);
    element.scrollIntoView({
        block: tall ? 'start' : 'center',
        inline: 'nearest',
        behavior: prefersReducedMotion() ? 'auto' : 'smooth'
    });
    return true;
};

const unique = (ids) => [...new Set((ids || []).filter((id) => typeof id === 'string' && id))];

export function createFocusController({
    root = typeof document !== 'undefined' ? document : null,
    getCards = () => [],
    openPage = null,
    onMarks = null,
    onShow = null,
    // Runs before a shown card is scrolled to (the dock lowers a tall sheet
    // so the card has room above it); may return a promise.
    beforeScroll = null,
    waitMs = FOCUS_WAIT_MS,
    // Whether a click target belongs to the Copilot itself (its own clicks —
    // a chip, the composer — never clear the focus).
    isOwnTarget = (target) => Boolean(target?.closest?.('.llm-dock')),
    scroll = scrollCardIntoView
} = {}) {
    let state = { cards: [], primary: null, marks: {}, source: null };
    let preview = null;
    const running = new Map();
    const listeners = new Set();
    let cleanups = [];
    let generation = 0;
    let observer = null;
    let applyQueued = false;

    const desired = () => {
        const wanted = new Map();
        running.forEach((ids) => ids.forEach((id) => wanted.set(id, 'running')));
        state.cards.forEach((id) => wanted.set(id, id === state.primary ? 'primary' : 'secondary'));
        if (preview && !wanted.has(preview)) wanted.set(preview, 'preview');
        return wanted;
    };

    // Set the attribute on every element of every wanted card (hidden pages
    // included, so a page the user opens later shows it too) and remove it
    // everywhere else.
    const apply = () => {
        if (!root) return;
        const wanted = desired();
        root.querySelectorAll(`[${FOCUS_ATTR}]`).forEach((element) => {
            const id = element.getAttribute(CARD_ATTR);
            if (!wanted.has(id)) element.removeAttribute(FOCUS_ATTR);
        });
        wanted.forEach((level, id) => {
            root.querySelectorAll(cardSelector(id)).forEach((element) => {
                if (element.getAttribute(FOCUS_ATTR) !== level) element.setAttribute(FOCUS_ATTR, level);
            });
        });
        watch(wanted.size > 0);
    };

    // While anything is in focus, re-apply when the host re-renders a card
    // (a page mounts, a plot card appears once its file is parsed).
    const watch = (on) => {
        if (typeof MutationObserver === 'undefined' || !root) return;
        if (on && !observer) {
            observer = new MutationObserver(() => {
                if (applyQueued) return;
                applyQueued = true;
                const run = () => { applyQueued = false; apply(); };
                if (typeof requestAnimationFrame === 'function') requestAnimationFrame(run);
                else setTimeout(run, 0);
            });
            observer.observe(root === document ? document.body : root, { childList: true, subtree: true });
        } else if (!on && observer) {
            observer.disconnect();
            observer = null;
        }
    };

    const notify = () => listeners.forEach((listener) => listener(state));

    const runCleanups = () => {
        const pending = cleanups;
        cleanups = [];
        pending.forEach((fn) => {
            try {
                fn();
            } catch {
                // A cleanup that fails must not keep the focus from clearing.
            }
        });
    };

    const clear = (reason = 'clear') => {
        generation += 1;
        const had = state.cards.length > 0 || Object.keys(state.marks).length > 0 || cleanups.length > 0;
        runCleanups();
        state = { cards: [], primary: null, marks: {}, source: null, reason };
        if (had) onMarks?.({});
        apply();
        if (had) notify();
    };

    const waitForCard = (id, token) => new Promise((resolve) => {
        const started = Date.now();
        const poll = () => {
            if (token !== generation) return resolve(null);
            const element = findCard(root, id);
            if (element) return resolve(element);
            if (Date.now() - started >= waitMs) return resolve(null);
            setTimeout(poll, POLL_MS);
        };
        poll();
    });

    /**
     * Point at `ids` (the first is the main one): switch to its page when it
     * is not on screen, wait for it, scroll it into view. Resolves to what was
     * shown: `{ shown, missing, page, switched, scrolled, superseded }`.
     */
    const show = async (ids, { marks = null, scroll: doScroll = true, source = 'agent' } = {}) => {
        const cards = unique(ids);
        runCleanups();
        generation += 1;
        const token = generation;
        const markMap = marks && cards.length && Array.isArray(marks) && marks.length
            ? { [cards[0]]: marks }
            : {};
        state = { cards, primary: cards[0] || null, marks: markMap, source };
        onMarks?.(markMap);
        apply();
        notify();
        if (!cards.length) return { shown: [], missing: [], page: null, switched: false, scrolled: false };

        const catalog = getCards() || [];
        const entry = catalog.find((card) => card.id === state.primary) || null;
        let element = findCard(root, state.primary);
        let switched = false;
        if (!element && entry?.page && openPage) {
            openPage(entry.page);
            switched = true;
        }
        if (!element) element = await waitForCard(state.primary, token);
        if (token !== generation) return { shown: [], missing: cards, page: entry?.page ?? null, switched, scrolled: false, superseded: true };
        apply();
        if (element && doScroll && beforeScroll) {
            await beforeScroll();
            if (token !== generation) return { shown: [], missing: cards, page: entry?.page ?? null, switched, scrolled: false, superseded: true };
        }
        const scrolled = Boolean(element && doScroll && scroll(element));
        onShow?.({ ids: cards, element, switched });
        const shown = cards.filter((id) => findCard(root, id));
        return {
            shown,
            missing: cards.filter((id) => !shown.includes(id)),
            page: entry?.page ?? null,
            switched,
            scrolled
        };
    };

    const onClick = (event) => {
        if (!state.cards.length && !cleanups.length) return;
        if (isOwnTarget(event.target)) return;
        clear('user');
    };

    return {
        show,
        clear,
        /** Register an undo for something drawn outside the cards; runs on the next clear or show. */
        onClear: (fn) => { if (typeof fn === 'function') cleanups.push(fn); },
        running: (stepId, ids) => {
            const list = unique(ids);
            if (!list.length) return;
            running.set(stepId, list);
            apply();
        },
        done: (stepId) => {
            if (running.delete(stepId)) apply();
        },
        doneAll: () => {
            if (!running.size) return;
            running.clear();
            apply();
        },
        preview: (id) => {
            preview = id || null;
            apply();
        },
        state: () => state,
        subscribe: (listener) => {
            listeners.add(listener);
            return () => listeners.delete(listener);
        },
        attach: () => {
            if (!root) return () => {};
            // A click, not pointerdown: scrolling a touch screen toward the
            // card must not end its highlight.
            root.addEventListener('click', onClick, true);
            return () => root.removeEventListener('click', onClick, true);
        },
        dispose: () => {
            generation += 1;
            running.clear();
            preview = null;
            runCleanups();
            state = { cards: [], primary: null, marks: {}, source: null };
            apply();
            watch(false);
        }
    };
}
