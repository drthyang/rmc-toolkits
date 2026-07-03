import React from 'react';
import { contextToJson } from '../context/runContext';
import { buildSummaryMessages } from '../prompts/templates';
import { useStreamedReply } from '../useStreamedReply';

// One-click natural-language assessment of the loaded run, streamed from the
// local model. The collapsible inspector shows the exact context JSON sent to
// the model — the point of this module is to make the pipeline legible.

const RunSummaryView = ({ context, settings, connected }) => {
    const { text, streaming, error, start, stop } = useStreamedReply();

    const handleAssess = () => start({
        baseUrl: settings.baseUrl,
        model: settings.model,
        temperature: settings.temperature,
        messages: buildSummaryMessages(context)
    });

    return (
        <div className="llm-view">
            <div className="llm-view-actions">
                <button type="button" onClick={handleAssess} disabled={!connected || streaming}>
                    {streaming ? 'Assessing…' : 'Assess this run'}
                </button>
                {streaming && (
                    <button type="button" className="llm-stop" onClick={stop}>Stop</button>
                )}
            </div>
            {!connected && (
                <p className="llm-view-hint">Connect to a local model in Settings to generate an assessment.</p>
            )}
            {error && <div className="llm-connection-status is-error" role="alert">{error}</div>}
            {text && <div className="llm-reply-text">{text}</div>}
            <details className="llm-context-inspector">
                <summary>Context sent to the model</summary>
                <pre>{contextToJson(context)}</pre>
            </details>
        </div>
    );
};

export default RunSummaryView;
