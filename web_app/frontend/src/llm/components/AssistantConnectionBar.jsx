// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';

// The always-visible connection controls: a status indicator, a model switcher,
// and a gear that toggles the settings drawer. Rendered in the page header and
// in the dashboard card, driven by the shared useAssistant hook.

// The status pill's label + state class. Before the user has actively tested —
// on first load, or when the automatic probe fails — show a gentle prompt in
// the neutral (idle) style rather than a red "Connection failed"; that error is
// reserved for a failed manual test.
const statusDisplay = (connection) => {
    if (connection.status === 'ok') return { state: 'ok', text: 'Connected' };
    if (connection.status === 'testing') return { state: 'testing', text: 'Connecting…' };
    if (connection.status === 'error' && connection.manual) {
        return { state: 'error', text: 'Connection failed' };
    }
    return { state: 'idle', text: 'Connect a model to start' };
};

const GearIcon = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <circle cx="12" cy="12" r="3" />
        <path d="M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 1 1-2.83 2.83l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 0 1-4 0v-.09A1.65 1.65 0 0 0 9 19.4a1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 1 1-2.83-2.83l.06-.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 0 1 0-4h.09A1.65 1.65 0 0 0 4.6 9a1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 1 1 2.83-2.83l.06.06a1.65 1.65 0 0 0 1.82.33H9a1.65 1.65 0 0 0 1-1.51V3a2 2 0 0 1 4 0v.09a1.65 1.65 0 0 0 1 1.51 1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 1 1 2.83 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82V9a1.65 1.65 0 0 0 1.51 1H21a2 2 0 0 1 0 4h-.09a1.65 1.65 0 0 0-1.51 1z" />
    </svg>
);

const AssistantConnectionBar = ({ connection, settings, onSave, showSettings, onToggleSettings }) => {
    const hasModels = connection.models.length > 0;
    const status = statusDisplay(connection);

    return (
        <div className="llm-connbar">
            <span className={`llm-status is-${status.state}`} role="status">
                <span className="llm-status-dot" aria-hidden="true" />
                {status.text}
            </span>
            <div className="llm-connbar-controls">
                {hasModels ? (
                    <select
                        className="llm-model-select"
                        value={settings.model}
                        onChange={(event) => onSave({ model: event.target.value })}
                        aria-label="Model"
                    >
                        <option value="">Select a model…</option>
                        {connection.models.map((model) => (
                            <option key={model} value={model}>{model}</option>
                        ))}
                    </select>
                ) : (
                    <input
                        className="llm-model-input"
                        type="text"
                        placeholder="model name"
                        value={settings.model}
                        spellCheck="false"
                        onChange={(event) => onSave({ model: event.target.value })}
                        aria-label="Model"
                    />
                )}
                <button
                    type="button"
                    className={`llm-icon-btn${showSettings ? ' is-active' : ''}`}
                    onClick={onToggleSettings}
                    aria-expanded={showSettings}
                    aria-label="Connection settings"
                    title="Connection settings"
                >
                    <GearIcon />
                </button>
            </div>
        </div>
    );
};

export default AssistantConnectionBar;
