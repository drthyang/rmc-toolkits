import React from 'react';
import { buildReportMessages } from '../prompts/templates';
import { buildReportMarkdown, downloadReport } from '../report/buildReport';
import { useStreamedReply } from '../useStreamedReply';

// One-click Markdown report: deterministic tables (model summary, fit metrics,
// convergence) always work; the AI assessment section streams in when a model
// is connected, and the report can be downloaded without it.

const ReportView = ({ context, settings, connected }) => {
    const { text: narrative, streaming, error, start, stop } = useStreamedReply();

    const handleGenerate = () => start({
        baseUrl: settings.baseUrl,
        model: settings.model,
        temperature: settings.temperature,
        apiKey: settings.apiKey,
        messages: buildReportMessages(context)
    });

    const handleDownload = () => {
        const markdown = buildReportMarkdown(context, narrative || null, {
            model: narrative ? settings.model : null,
            temperature: narrative ? settings.temperature : null
        });
        downloadReport(markdown, context.run);
    };

    return (
        <div className="llm-view">
            <div className="llm-view-actions">
                <button type="button" className="llm-primary" onClick={handleGenerate} disabled={!connected || streaming}>
                    {streaming ? 'Writing assessment…' : narrative ? 'Regenerate AI assessment' : 'Generate AI assessment'}
                </button>
                {streaming && (
                    <button type="button" className="llm-stop" onClick={stop}>Stop</button>
                )}
                <button type="button" onClick={handleDownload} disabled={streaming}>
                    {narrative ? 'Download report (.md)' : 'Download report without AI assessment (.md)'}
                </button>
            </div>
            {!connected && (
                <p className="llm-view-hint">
                    The metrics tables work offline; connect a local model in Settings to add a narrative assessment.
                </p>
            )}
            {error && <div className="llm-connection-status is-error" role="alert">{error}</div>}
            {narrative && <div className="llm-reply-text">{narrative}</div>}
        </div>
    );
};

export default ReportView;
