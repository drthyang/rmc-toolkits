// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useEffect, useId, useMemo, useRef, useState } from 'react';
import { buildRunTools } from '../copilot/runTools';
import { useAssistant } from '../useAssistant';
import AssistantConnectionBar from './AssistantConnectionBar';
import ConnectionSettings from './ConnectionSettings';
import CopilotView from './CopilotView';
import { CloseIcon, ExpandIcon, ShrinkIcon, WaveMark } from './icons';
import './Copilot.css';
import './CopilotChat.css';

// The AI Copilot as a chatbox over every workbench page: a labelled launcher
// in the bottom-right corner opens a floating panel with the conversation, so
// the analysis the Copilot ran can be shown on its page ("Show in <page>")
// while the chat stays in view. The panel mounts on first open — no network
// request before that — and then stays mounted while closed, so a running
// answer continues and the conversation is kept. Run data and the host's
// analysis tools arrive as props; the tools over the run data itself are
// built here. `pages` maps a tool's page key to its tab label;
// `onOpenPage(key, params)` switches the workbench there and loads `params`
// (on a phone-size screen, where the panel is a full sheet, it then closes);
// `runKey` names the open run (a change restarts the conversation);
// `onOpenChange(open)` tells the host when the chatbox opens or closes.

// Where the panel is a full-screen sheet (CopilotChat.css): a page shown
// from a step would sit hidden behind it, so the sheet closes instead.
const SHEET_QUERY = '(max-width: 760px), (max-height: 540px)';
const isSheet = () => typeof window !== 'undefined' && typeof window.matchMedia === 'function'
    && window.matchMedia(SHEET_QUERY).matches;

const CopilotPanel = ({ id, hidden, expanded, onToggleExpand, onClose, run, hostTools, pages, onOpenPage, runKey }) => {
    const assistant = useAssistant({ ...run, enabled: true });
    const hasRun = Boolean((run.plotFiles && run.plotFiles.length) || run.structure);
    const { runName, plotFiles, rValueFile, structure, symmetry, runSettings, pcaSites, liveData } = run;

    const tools = useMemo(() => [
        ...buildRunTools({ runName, plotFiles, rValueFile, structure, symmetry, runSettings, pcaSites, liveData }),
        ...(hostTools || [])
    ], [runName, plotFiles, rValueFile, structure, symmetry, runSettings, pcaSites, liveData, hostTools]);

    return (
        <section
            id={id}
            className={`llm-chat-panel${expanded ? ' is-expanded' : ''}${hidden ? ' is-closed' : ''}`}
            role="dialog"
            aria-label="AI Copilot"
            aria-hidden={hidden}
            tabIndex={-1}
            onKeyDown={(event) => {
                if (event.key === 'Escape' && !event.defaultPrevented) {
                    event.stopPropagation();
                    onClose();
                }
            }}
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
                        className="llm-icon-btn"
                        onClick={onToggleExpand}
                        aria-pressed={expanded}
                        aria-label={expanded ? 'Make the chat smaller' : 'Make the chat larger'}
                        title={expanded ? 'Smaller' : 'Larger'}
                    >
                        {expanded ? <ShrinkIcon /> : <ExpandIcon />}
                    </button>
                    <button
                        type="button"
                        className="llm-icon-btn"
                        onClick={onClose}
                        aria-label="Close the Copilot"
                        title="Close (Esc) — the conversation is kept"
                    >
                        <CloseIcon />
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
    tools = [], pages = {}, onOpenPage, runKey = null, onOpenChange, defaultOpen = false
}) => {
    const [open, setOpen] = useState(defaultOpen);
    const [mounted, setMounted] = useState(defaultOpen);
    const [expanded, setExpanded] = useState(false);
    const panelId = useId();
    const launcherRef = useRef(null);
    const shellRef = useRef(null);
    // Where focus goes after the next open/close: the composer, or the launcher.
    const focusRef = useRef(null);

    const setOpenState = (next) => {
        if (next) setMounted(true);
        focusRef.current = next ? 'composer' : 'launcher';
        setOpen(next);
        onOpenChange?.(next);
    };

    useEffect(() => {
        const target = focusRef.current;
        if (!target) return;
        focusRef.current = null;
        if (target === 'launcher') {
            launcherRef.current?.focus();
        } else {
            const composer = shellRef.current?.querySelector('.llm-composer-input:not(:disabled)');
            (composer || shellRef.current?.querySelector('.llm-chat-panel'))?.focus?.();
        }
    }, [open]);

    const run = { runName, plotFiles, rValueFile, structure, symmetry, runSettings, pcaSites, liveData };

    const openPage = onOpenPage
        ? (page, params) => {
            onOpenPage(page, params);
            if (isSheet()) setOpenState(false);
        }
        : undefined;

    return (
        <div className="llm-chat-shell" ref={shellRef}>
            {!open && (
                <button
                    ref={launcherRef}
                    type="button"
                    className="llm-chat-launcher"
                    onClick={() => setOpenState(true)}
                    aria-expanded="false"
                    aria-controls={mounted ? panelId : undefined}
                    title="Ask the Copilot to run and check analyses on this run"
                >
                    <WaveMark />
                    <span>AI Copilot</span>
                </button>
            )}
            {mounted && (
                <CopilotPanel
                    id={panelId}
                    hidden={!open}
                    expanded={expanded}
                    onToggleExpand={() => setExpanded((value) => !value)}
                    onClose={() => setOpenState(false)}
                    run={run}
                    hostTools={tools}
                    pages={pages}
                    onOpenPage={openPage}
                    runKey={runKey}
                />
            )}
        </div>
    );
};

export default CopilotChat;
