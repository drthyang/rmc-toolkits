// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React from 'react';
import InfoBadge from '../../ui/InfoBadge';

// The context meter by the chat box (THY-76): "Context 62%" of the active
// model's window, the reply reserve counted, in the guard's own numbers
// (useCopilot.js meter, provider/contextMeter.js). Neutral below 70 %, amber
// to 90 %, red above; at amber or red it offers Start fresh and, when it
// would shrink the next request, Summarize and continue. The details sit
// behind "?". On a phone it stays folded away until amber (Copilot.css).

const tokens = (count) => Math.round(count).toLocaleString('en-US');

const shortTokens = (count) => (count >= 10000 ? `${Math.round(count / 1000)}k` : tokens(count));

const percent = (share) => `${Math.round(share * 100)}%`;

const usd = (amount) => (amount < 0.01 ? `$${amount.toFixed(4)}` : `$${amount.toFixed(2)}`);

const SERVER_WINDOW = {
    ollama_max: (meter) => `the model's maximum on Ollama${meter.limits.loadedContext ? ` (loaded now with ${tokens(meter.limits.loadedContext)})` : ''}`,
    lmstudio_loaded: () => 'the context LM Studio loaded the model with',
    published: (meter) => `${meter.model}'s published window`
};

// Where the window is not known, why.
const unknownWindow = (meter) => {
    if (meter.limits.server === 'lmstudio') return 'LM Studio has not loaded the model yet, so its window is not known.';
    if (meter.limits.server === 'ollama') return 'Ollama does not report this model\'s maximum context.';
    return `The window of ${meter.model || 'this model'} is not known here.`;
};

const ContextMeter = ({ meter, busy = false, onStartFresh = null, onSummarize = null }) => {
    if (!meter) return null;
    const { level, share, window } = meter;
    const label = Number.isFinite(share) ? `Context ${percent(share)}` : `Context ${shortTokens(meter.used)} tokens`;
    const act = (level === 'amber' || level === 'red') && !busy;
    return (
        <span className={`llm-meter is-${level}`} data-level={level}>
            <span className="llm-meter-main">
                {Number.isFinite(share) && (
                    <span className="llm-meter-bar" aria-hidden="true">
                        <span className="llm-meter-fill" style={{ width: `${Math.min(100, share * 100)}%` }} />
                    </span>
                )}
                <span className="llm-meter-label">{label}</span>
                <InfoBadge label="About the context meter" align="end" side="above">
                    <span className="llm-meter-details">
                        {meter.refused ? (
                            <span>Refused: this request needs about {tokens(meter.used)} tokens; the model holds {window ? tokens(window) : 'less'}. Nothing was sent.</span>
                        ) : (
                            <span>
                                {tokens(meter.used)} {window ? `of ${tokens(window)} ` : ''}tokens: the prompt {tokens(meter.prompt)} + the reply
                                reserve {tokens(meter.reserve)}.
                            </span>
                        )}
                        {window && meter.limits.source
                            ? <span>Window: {SERVER_WINDOW[meter.limits.source](meter)}.</span>
                            : (!window && <span>{unknownWindow(meter)}</span>)}
                        {meter.running ? (
                            <span>Shown: the latest request of this answer.</span>
                        ) : (meter.next && meter.last && !meter.refused && (
                            <span>
                                Next question starts at {tokens(meter.next.used)}; the last answer peaked at {tokens(meter.last.used)}.
                            </span>
                        ))}
                        <span>
                            Estimated high (2.5 characters a token, 1.4 for analysis results; about 30% above a real count),
                            the same count the Copilot uses to refuse a request above 100%.
                        </span>
                        {meter.cost != null && <span>Cost so far: about {usd(meter.cost)} (same estimate; hidden thinking not included).</span>}
                        {meter.summarized && <span>Earlier answers are sent summarized.</span>}
                    </span>
                </InfoBadge>
            </span>
            {act && (
                <span className="llm-meter-actions">
                    {onStartFresh && <button type="button" className="llm-link-btn" onClick={onStartFresh}>Start fresh</button>}
                    {meter.canSummarize && onSummarize && (
                        <button type="button" className="llm-link-btn" onClick={onSummarize}>Summarize and continue</button>
                    )}
                </span>
            )}
        </span>
    );
};

export default ContextMeter;
