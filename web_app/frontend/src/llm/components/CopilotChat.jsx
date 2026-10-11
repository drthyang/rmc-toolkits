// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useCallback, useEffect, useId, useMemo, useRef, useState } from 'react';
import { buildRunTools } from '../copilot/runTools';
import { buildFocusTools } from '../copilot/viewTools';
import { createFocusController } from '../focus/focus';
import { useAssistant } from '../useAssistant';
import AssistantConnectionBar from './AssistantConnectionBar';
import ConnectionSettings from './ConnectionSettings';
import CopilotView from './CopilotView';
import {
    DOCK_QUERY, DOCK_STORAGE_KEY, SHEET_MAX, SHEET_MIN,
    clampDockWidth, clampSheet, defaultDockWidth, dragDockWidth, dragSheet, keyResize, maxDockWidth,
    readDockState, sheetForCard, writeDockState, DOCK_MIN
} from './dockLayout';
import { ChevronIcon, CloseIcon, WaveMark } from './icons';
import './Copilot.css';
import './CopilotChat.css';

// The AI Copilot docked beside the workbench pages instead of over them. The
// host renders it once, as the last child of its app grid, and places it by
// its `data-mode`: "docked" (a column to the right of the pages, resizable
// from its left edge) or "sheet" (a row under them, on portrait tablets, phones and
// short windows, resizable from its top edge). Collapsed, it is a slim rail
// (or bar) with an unread mark; the panel mounts on first open — no network
// request before that — and stays mounted while collapsed, so a running
// answer continues and the conversation is kept. The open state, width and
// sheet height persist (dockLayout.js).
//
// The dock owns the card-focus controller (focus/focus.js): the agent's
// highlight_cards tool, the cards of running analyses, and the card chips in
// its answers all go through it; `cards` lists what the host offers
// (id, page, title, about, kind/file, marks), `onOpenPage(key, params)`
// switches the page, `onMarks(marksByCard)` hands marks inside a plot card to
// the host. Run data and the host's tools arrive as props; the tools over the
// run data itself are built here.

const hasMatchMedia = () => typeof window !== 'undefined' && typeof window.matchMedia === 'function';
const viewportWidth = () => (typeof window !== 'undefined' ? window.innerWidth : 1440);
const viewportHeight = () => (typeof window !== 'undefined' ? window.innerHeight : 900);
const rootFontSize = () => {
    if (typeof window === 'undefined') return 15;
    const size = parseFloat(window.getComputedStyle(document.documentElement).fontSize);
    return Number.isFinite(size) && size > 0 ? size : 15;
};
const storage = () => {
    try {
        return typeof window !== 'undefined' ? window.localStorage : null;
    } catch {
        return null;
    }
};
// The card's nearest vertically scrolling ancestor (the page, or the shell).
const scrollAncestor = (element) => {
    for (let node = element.parentElement; node; node = node.parentElement) {
        if (/(auto|scroll|overlay)/.test(window.getComputedStyle(node).overflowY)) return node;
    }
    return null;
};
const nextFrame = () => new Promise((resolve) => {
    if (typeof requestAnimationFrame === 'function') requestAnimationFrame(() => requestAnimationFrame(() => resolve()));
    else setTimeout(resolve, 0);
});

// "docked" when the window is wide enough for the pages and the panel side by
// side, else "sheet". Without matchMedia (tests), docked.
const useDockMode = () => {
    const [docked, setDocked] = useState(() => !hasMatchMedia() || window.matchMedia(DOCK_QUERY).matches);
    useEffect(() => {
        if (!hasMatchMedia()) return undefined;
        const query = window.matchMedia(DOCK_QUERY);
        const update = () => setDocked(query.matches);
        update();
        query.addEventListener?.('change', update);
        return () => query.removeEventListener?.('change', update);
    }, []);
    return docked ? 'docked' : 'sheet';
};

// The edge the panel is resized from: a vertical separator on the docked
// panel's left edge, a grab bar on the sheet's top edge. Pointer drag (with
// capture), arrow keys, Home/End; double-click or Enter restores the default.
// A drag moves the size variable on the dock element directly and commits
// the state once, on release, so the conversation and the host's pages are
// not re-rendered on every pointer move.
const ResizeHandle = ({ mode, width, sheet, dockRef, onResize, onCommit, onReset }) => {
    const dragRef = useRef(null);
    const docked = mode === 'docked';
    const vw = viewportWidth();
    const value = docked ? width : Math.round(sheet * 100);
    const onPointerDown = (event) => {
        if (event.button !== 0) return;
        event.preventDefault();
        event.currentTarget.setPointerCapture?.(event.pointerId);
        dragRef.current = { x: event.clientX, width, last: null };
    };
    const onPointerMove = (event) => {
        const drag = dragRef.current;
        if (!drag) return;
        drag.last = docked
            ? { width: dragDockWidth(drag.width, drag.x, event.clientX, viewportWidth()) }
            : { sheet: dragSheet(event.clientY, viewportHeight()) };
        const dock = dockRef.current;
        if (!dock) return;
        if (docked) dock.style.setProperty('--llm-dock-w', `${drag.last.width}px`);
        else dock.style.setProperty('--llm-sheet-f', String(drag.last.sheet));
    };
    const onPointerUp = () => {
        const drag = dragRef.current;
        dragRef.current = null;
        if (!drag?.last) return;
        onResize(drag.last);
        onCommit(drag.last);
    };
    return (
        <div
            className={`llm-dock-resize is-${mode}`}
            role="separator"
            aria-orientation={docked ? 'vertical' : 'horizontal'}
            aria-label={docked ? 'Resize the Copilot panel' : 'Resize the Copilot sheet'}
            aria-valuemin={docked ? DOCK_MIN : Math.round(SHEET_MIN * 100)}
            aria-valuemax={docked ? maxDockWidth(vw) : Math.round(SHEET_MAX * 100)}
            aria-valuenow={value}
            tabIndex={0}
            title="Drag to resize · double-click to reset"
            onPointerDown={onPointerDown}
            onPointerMove={onPointerMove}
            onPointerUp={onPointerUp}
            onPointerCancel={onPointerUp}
            onDoubleClick={onReset}
            onKeyDown={(event) => {
                if (event.key === 'Enter') {
                    event.preventDefault();
                    onReset();
                    return;
                }
                const next = keyResize(mode, event.key, docked ? width : sheet, viewportWidth());
                if (next == null) return;
                event.preventDefault();
                onResize(docked ? { width: next } : { sheet: next });
                onCommit(docked ? { width: next } : { sheet: next });
            }}
        >
            <span className="llm-dock-resize-grip" aria-hidden="true" />
        </div>
    );
};

const CopilotPanel = ({
    id, hidden, mode, onCollapse, run, hostTools, pages, onOpenPage, runKey, cards, focus, onTurnEnd, onBusyChange
}) => {
    const assistant = useAssistant({ ...run, enabled: true });
    const hasRun = Boolean((run.plotFiles && run.plotFiles.length) || run.structure);
    const { runName, plotFiles, rValueFile, structure, symmetry, runSettings, pcaSites, liveData } = run;

    const tools = useMemo(() => [
        ...buildRunTools({ runName, plotFiles, rValueFile, structure, symmetry, runSettings, pcaSites, liveData }, { cards }),
        ...buildFocusTools({ cards, pages }),
        ...(hostTools || [])
    ], [runName, plotFiles, rValueFile, structure, symmetry, runSettings, pcaSites, liveData, hostTools, cards, pages]);

    return (
        <section
            id={id}
            className={`llm-chat-panel${hidden ? ' is-closed' : ''}`}
            aria-label="AI Copilot conversation"
            aria-hidden={hidden}
            tabIndex={-1}
        >
            <header className="llm-chat-panel-head">
                <div className="llm-chat-panel-title">
                    <span className="llm-chat-panel-mark" aria-hidden="true"><WaveMark /></span>
                    <h2>AI Copilot</h2>
                    <span className="llm-beta-tag">Beta</span>
                </div>
                <div className="llm-chat-panel-actions">
                    <button
                        type="button"
                        className={`llm-icon-btn llm-dock-collapse is-${mode}`}
                        onClick={onCollapse}
                        aria-label="Collapse the Copilot"
                        title="Collapse (Esc) — the conversation is kept"
                    >
                        {mode === 'docked' ? <ChevronIcon /> : <CloseIcon />}
                    </button>
                </div>
            </header>

            <div className="llm-chat-panel-conn">
                <AssistantConnectionBar
                    connection={assistant.connection}
                    settings={assistant.settings}
                    onSave={assistant.saveSettings}
                    showSettings={assistant.showSettings}
                    onToggleSettings={() => assistant.setShowSettings((value) => !value)}
                />
            </div>

            <div className="llm-chat-panel-body">
                {assistant.showSettings && (
                    <div className="llm-chat-panel-settings">
                        <ConnectionSettings
                            settings={assistant.settings}
                            connection={assistant.connection}
                            onSave={assistant.saveSettings}
                            onTest={assistant.runTest}
                        />
                    </div>
                )}
                {/* Kept mounted under the settings, so the conversation survives. */}
                <div className={`llm-chat-panel-conversation${assistant.showSettings ? ' is-covered' : ''}`}>
                    {hasRun ? (
                        <CopilotView
                            context={assistant.context}
                            settings={assistant.settings}
                            connected={assistant.connected}
                            tools={tools}
                            toolsEnabled={assistant.settings.copilotTools !== false}
                            pages={pages}
                            onOpenPage={onOpenPage}
                            runKey={runKey}
                            cards={cards}
                            focus={focus}
                            onTurnEnd={onTurnEnd}
                            onBusyChange={onBusyChange}
                        />
                    ) : (
                        <p className="llm-chat-norun">Open a run folder, or load the Demo, to start.</p>
                    )}
                </div>
            </div>
        </section>
    );
};

const CopilotChat = ({
    runName, plotFiles, rValueFile, structure, symmetry, runSettings, pcaSites, liveData,
    tools = [], pages = {}, onOpenPage, runKey = null, onOpenChange, defaultOpen = false,
    cards = null, onMarks = null, onLayoutChange = null, className = '', storageKey = DOCK_STORAGE_KEY
}) => {
    const mode = useDockMode();
    const [saved] = useState(() => readDockState(storage(), storageKey));
    const initialOpen = defaultOpen || saved.open;
    const [open, setOpen] = useState(initialOpen);
    const [mounted, setMounted] = useState(initialOpen);
    const [width, setWidth] = useState(() => clampDockWidth(saved.width ?? defaultDockWidth(viewportWidth(), rootFontSize()), viewportWidth()));
    const [sheet, setSheet] = useState(() => clampSheet(saved.sheet));
    const [unread, setUnread] = useState(false);
    const [busy, setBusy] = useState(false);
    const panelId = useId();
    const railRef = useRef(null);
    const shellRef = useRef(null);
    // Where focus goes after the next open/close: the composer, or the rail.
    const focusTargetRef = useRef(null);
    const openRef = useRef(open);
    openRef.current = open;
    const modeRef = useRef(mode);
    modeRef.current = mode;
    const sheetRef = useRef(sheet);
    sheetRef.current = sheet;

    // The card-focus controller, one per dock; it reads the latest props
    // through refs.
    const cardsRef = useRef(cards);
    cardsRef.current = cards;
    const openPageRef = useRef(onOpenPage);
    openPageRef.current = onOpenPage;
    const onMarksRef = useRef(onMarks);
    onMarksRef.current = onMarks;
    const [focus] = useState(() => createFocusController({
        getCards: () => cardsRef.current || [],
        openPage: (page) => openPageRef.current?.(page, null),
        onMarks: (marks) => onMarksRef.current?.(marks),
        beforeScroll: async (card) => {
            // On a sheet, lower it until the card fits whole above it (not
            // remembered — the user's height is). The space above the card's
            // scroll area is its own top, or the pinned header on a phone
            // (the card's scroll-margin-top).
            if (modeRef.current !== 'sheet' || !card) return;
            const scroller = scrollAncestor(card);
            const top = scroller ? Math.max(0, scroller.getBoundingClientRect().top) : 0;
            const margin = parseFloat(window.getComputedStyle(card).scrollMarginTop) || 0;
            const next = sheetForCard({
                viewportHeight: viewportHeight(),
                topOffset: top + margin,
                cardHeight: card.getBoundingClientRect().height,
                current: sheetRef.current
            });
            if (next >= sheetRef.current) return;
            setSheet(next);
            await nextFrame();
        }
    }));
    useEffect(() => {
        const detach = focus.attach();
        return () => {
            detach();
            focus.dispose();
        };
    }, [focus]);

    // Another run: its highlight means nothing on the new one.
    const lastRunKeyRef = useRef(runKey);
    useEffect(() => {
        if (lastRunKeyRef.current === runKey) return;
        lastRunKeyRef.current = runKey;
        focus.clear('run');
    }, [runKey, focus]);

    // An initially open dock (remembered) tells the host once.
    useEffect(() => {
        if (initialOpen) onOpenChange?.(true);
    // Mount only.
    // eslint-disable-next-line react-hooks/exhaustive-deps
    }, []);

    const setOpenState = useCallback((next) => {
        if (next) {
            setMounted(true);
            setUnread(false);
        } else {
            // Collapsed, the agent is not on screen: its highlight ends.
            focus.clear('collapse');
        }
        focusTargetRef.current = next ? 'composer' : 'rail';
        setOpen(next);
        writeDockState(storage(), { open: next }, storageKey);
        onOpenChange?.(next);
    }, [focus, onOpenChange, storageKey]);

    useEffect(() => {
        const target = focusTargetRef.current;
        if (!target) return;
        focusTargetRef.current = null;
        if (target === 'rail') {
            railRef.current?.focus();
        } else {
            const composer = shellRef.current?.querySelector('.llm-composer-input:not(:disabled)');
            (composer || shellRef.current?.querySelector('.llm-chat-panel'))?.focus?.();
        }
    }, [open]);

    // Keep the docked width inside the window as it resizes.
    useEffect(() => {
        if (typeof window === 'undefined') return undefined;
        const onResize = () => setWidth((current) => clampDockWidth(current, viewportWidth()));
        window.addEventListener('resize', onResize);
        return () => window.removeEventListener('resize', onResize);
    }, []);

    useEffect(() => {
        onLayoutChange?.({ mode, open, width, sheet });
    }, [onLayoutChange, mode, open, width, sheet]);

    const resize = useCallback(({ width: nextWidth, sheet: nextSheet }) => {
        if (nextWidth != null) setWidth(nextWidth);
        if (nextSheet != null) setSheet(nextSheet);
    }, []);
    const commit = useCallback((patch) => {
        writeDockState(storage(), patch, storageKey);
    }, [storageKey]);
    const reset = useCallback(() => {
        if (modeRef.current === 'docked') {
            const next = clampDockWidth(defaultDockWidth(viewportWidth(), rootFontSize()), viewportWidth());
            setWidth(next);
            writeDockState(storage(), { width: null }, storageKey);
        } else {
            setSheet(clampSheet(null));
            writeDockState(storage(), { sheet: clampSheet(null) }, storageKey);
        }
    }, [storageKey]);

    const onTurnEnd = useCallback(() => {
        if (!openRef.current) setUnread(true);
    }, []);

    const run = { runName, plotFiles, rValueFile, structure, symmetry, runSettings, pcaSites, liveData };
    const railLabel = unread ? 'AI Copilot — new answer' : (busy ? 'AI Copilot — working' : 'AI Copilot');

    return (
        <aside
            ref={shellRef}
            className={`llm-dock is-${mode} ${open ? 'is-open' : 'is-collapsed'}${className ? ` ${className}` : ''}`}
            data-mode={mode}
            data-open={open ? 'true' : 'false'}
            aria-label="AI Copilot"
            style={{ '--llm-dock-w': `${width}px`, '--llm-sheet-f': String(sheet) }}
            onKeyDown={(event) => {
                if (event.key === 'Escape' && open && !event.defaultPrevented) {
                    event.stopPropagation();
                    setOpenState(false);
                }
            }}
        >
            {open && (
                <ResizeHandle
                    mode={mode}
                    width={width}
                    sheet={sheet}
                    dockRef={shellRef}
                    onResize={resize}
                    onCommit={commit}
                    onReset={reset}
                />
            )}
            {!open && (
                <button
                    ref={railRef}
                    type="button"
                    className={`llm-dock-rail llm-chat-launcher is-${mode}${unread ? ' has-unread' : ''}${busy ? ' is-busy' : ''}`}
                    onClick={() => setOpenState(true)}
                    aria-expanded="false"
                    aria-controls={mounted ? panelId : undefined}
                    aria-label={railLabel}
                    title="Open the Copilot: it runs and checks analyses on this run, and points at the cards it explains"
                >
                    <span className="llm-dock-rail-icon" aria-hidden="true"><WaveMark /></span>
                    <span className="llm-dock-rail-label">AI Copilot</span>
                    {unread && <span className="llm-dock-rail-dot" aria-hidden="true" />}
                    {mode === 'sheet' && <span className="llm-dock-rail-chevron" aria-hidden="true"><ChevronIcon /></span>}
                </button>
            )}
            {mounted && (
                <CopilotPanel
                    id={panelId}
                    hidden={!open}
                    mode={mode}
                    onCollapse={() => setOpenState(false)}
                    run={run}
                    hostTools={tools}
                    pages={pages}
                    onOpenPage={onOpenPage}
                    runKey={runKey}
                    cards={cards}
                    focus={focus}
                    onTurnEnd={onTurnEnd}
                    onBusyChange={setBusy}
                />
            )}
        </aside>
    );
};

export default CopilotChat;
