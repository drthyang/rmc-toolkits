// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useState } from 'react';
import { ChevronIcon } from './icons';

// One analysis the Copilot ran: status, title, duration and the verdict that
// matters most, expandable to every check, the arguments and the exact data
// the model read. Data is shown as JSON text, never as HTML. "Show in <page>"
// switches the workbench to that analysis's page and loads the same analysis
// there ("Open <page>" when the step has nothing to load) — a user action;
// the Copilot itself never navigates.

const STATUS_TEXT = {
    running: 'Running…',
    ok: 'Done',
    warn: 'Warning',
    failed: 'Failed',
    stopped: 'Stopped'
};

const LEVEL_RANK = { fail: 0, warn: 1, ok: 2, note: 3 };

const formatMs = (ms) => {
    if (!Number.isFinite(ms)) return null;
    if (ms < 100) return '<0.1 s';
    return ms < 10000 ? `${(ms / 1000).toFixed(1)} s` : `${Math.round(ms / 1000)} s`;
};

const pretty = (value) => {
    if (value == null) return '—';
    return typeof value === 'string' ? value : JSON.stringify(value, null, 1);
};

// The single line shown under the title: the error, else the most severe check.
const headline = (step) => {
    if (step.error) return { level: 'fail', text: step.error };
    const ranked = [...(step.checks || [])].sort((a, b) => LEVEL_RANK[a.level] - LEVEL_RANK[b.level]);
    return ranked[0] || null;
};

const ToolStep = ({ step, pageLabel, onOpenPage }) => {
    const [open, setOpen] = useState(false);
    const running = step.status === 'running';
    const top = running ? null : headline(step);
    const meta = [STATUS_TEXT[step.status] || step.status, step.cached ? 'reused' : formatMs(step.ms)]
        .filter(Boolean)
        .join(' · ');

    return (
        <li className={`llm-step is-${step.status}`}>
            <button
                type="button"
                className="llm-step-head"
                onClick={() => setOpen((value) => !value)}
                aria-expanded={open}
                disabled={running}
            >
                <span className="llm-step-dot" aria-hidden="true" />
                <span className="llm-step-label">{step.label || step.name}</span>
                <span className="llm-step-meta">{meta}</span>
                {!running && <span className="llm-step-chevron" aria-hidden="true"><ChevronIcon /></span>}
            </button>
            {top && !open && <p className={`llm-step-headline is-${top.level}`}>{top.text}</p>}
            {open && (
                <div className="llm-step-body">
                    {step.error && <p className="llm-step-headline is-fail">{step.error}</p>}
                    {(step.checks || []).length > 0 && (
                        <ul className="llm-step-checks" aria-label="Checks">
                            {step.checks.map((check, index) => (
                                <li key={index} className={`llm-step-check is-${check.level}`}>
                                    <span className="llm-step-check-level">{check.level}</span>
                                    {check.text}
                                </li>
                            ))}
                        </ul>
                    )}
                    <div className="llm-step-section">
                        <span className="llm-step-section-title">{step.name} arguments</span>
                        <pre>{pretty(step.args)}</pre>
                    </div>
                    <div className="llm-step-section">
                        <span className="llm-step-section-title">{step.error ? 'Error sent to the model' : 'Result sent to the model'}</span>
                        <pre>{step.error ? step.error : pretty(step.data)}</pre>
                    </div>
                </div>
            )}
            {pageLabel && onOpenPage && !running && (
                <button
                    type="button"
                    className="ui-pill llm-step-open"
                    onClick={() => onOpenPage(step.page, step.pageParams ?? null)}
                    title={step.pageParams ? `Show this analysis on the ${pageLabel} page` : `Go to the ${pageLabel} page`}
                >
                    {step.pageParams ? `Show in ${pageLabel}` : `Open ${pageLabel}`}
                </button>
            )}
        </li>
    );
};

export default ToolStep;
