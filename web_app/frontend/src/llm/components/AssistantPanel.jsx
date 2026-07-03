import React, { useMemo, useState } from 'react';
import { buildRunContext } from '../context/runContext';
import { checkConnection } from '../provider/client';
import { saveSettings, useLlmSettings } from '../settings';
import ChatView from './ChatView';
import ConnectionSettings from './ConnectionSettings';
import ReportView from './ReportView';
import RunSummaryView from './RunSummaryView';
import './AssistantPanel.css';

// The assistant's shell: a collapsible dashboard card with Summary / Chat /
// Report tabs and a settings drawer. Collapsed by default and makes zero
// network requests until the user opens it and acts — users without a local
// LLM see nothing but a closed card.

const TABS = [
    { id: 'summary', label: 'Summary' },
    { id: 'chat', label: 'Chat' },
    { id: 'report', label: 'Report' }
];

const AssistantPanel = ({ runName, plotFiles, rValueFile, structure, symmetry, liveData }) => {
    const settings = useLlmSettings();
    const [expanded, setExpanded] = useState(false);
    const [activeTab, setActiveTab] = useState('summary');
    // Open by default so first-time users land on the connection setup; a
    // successful test keeps it visible (confirmation included) until hidden.
    const [showSettings, setShowSettings] = useState(true);
    const [connection, setConnection] = useState({ status: 'idle', models: [], error: null, hint: null });

    const context = useMemo(() => (
        expanded
            ? buildRunContext({ runName, plotFiles, rValueFile, structure, symmetry, liveData })
            : null
    ), [expanded, runName, plotFiles, rValueFile, structure, symmetry, liveData]);

    const handleTest = async () => {
        setConnection({ status: 'testing', models: [], error: null, hint: null });
        const result = await checkConnection(settings.baseUrl);
        setConnection({
            status: result.ok ? 'ok' : 'error',
            models: result.models,
            error: result.error,
            hint: result.hint
        });
        // Auto-pick a model when none is chosen (or the chosen one is gone).
        if (result.ok && result.models.length && !result.models.includes(settings.model)) {
            saveSettings({ model: result.models[0] });
        }
    };

    const connected = connection.status === 'ok' && Boolean(settings.model);

    return (
        <article className={`plot-card llm-assistant-card${expanded ? '' : ' is-collapsed'}`}>
            <div className="plot-card-header">
                <div className="llm-assistant-heading">
                    <h3>AI Assistant</h3>
                    <span className="llm-experimental-tag">Experimental</span>
                    {connection.status === 'ok' && (
                        <span className="llm-connection-dot" title={`Connected to ${settings.baseUrl}`} />
                    )}
                </div>
                <div className="plot-card-header-actions">
                    {expanded && (
                        <button
                            type="button"
                            className="panel-toggle"
                            onClick={() => setShowSettings((value) => !value)}
                            aria-expanded={showSettings}
                        >
                            Settings
                        </button>
                    )}
                    <button
                        type="button"
                        className="panel-toggle"
                        onClick={() => setExpanded((value) => !value)}
                        aria-expanded={expanded}
                    >
                        {expanded ? 'Hide' : 'Show'}
                    </button>
                </div>
            </div>
            {expanded && (
                <div className="llm-assistant-body">
                    {(showSettings || connection.status !== 'ok') && (
                        <ConnectionSettings
                            settings={settings}
                            connection={connection}
                            onSave={saveSettings}
                            onTest={handleTest}
                        />
                    )}
                    <nav className="llm-tabs" aria-label="Assistant features">
                        {TABS.map((tab) => (
                            <button
                                key={tab.id}
                                type="button"
                                className={activeTab === tab.id ? 'active' : ''}
                                onClick={() => setActiveTab(tab.id)}
                            >
                                {tab.label}
                            </button>
                        ))}
                    </nav>
                    {activeTab === 'summary' && (
                        <RunSummaryView context={context} settings={settings} connected={connected} />
                    )}
                    {activeTab === 'chat' && (
                        <ChatView context={context} settings={settings} connected={connected} />
                    )}
                    {activeTab === 'report' && (
                        <ReportView context={context} settings={settings} connected={connected} />
                    )}
                </div>
            )}
        </article>
    );
};

export default AssistantPanel;
