import React from 'react';
import { useAssistant } from '../useAssistant';
import AssistantConnectionBar from './AssistantConnectionBar';
import ChatView from './ChatView';
import ConnectionSettings from './ConnectionSettings';
import './AssistantPanel.css';
import './AssistantPage.css';

// Full-page presentation of the assistant, matching the app's other workspace
// pages (Dashboard, KDE/3D). The connection controls (status, model switcher,
// settings) live in the page header beside the title; the settings drawer opens
// below it; the tabs and active view sit in the content box.

const AssistantPage = ({ runName, plotFiles, rValueFile, structure, symmetry, runSettings, liveData }) => {
    const assistant = useAssistant({ runName, plotFiles, rValueFile, structure, symmetry, runSettings, liveData, enabled: true });
    const hasRun = Boolean((plotFiles && plotFiles.length) || structure);

    return (
        <section className="assistant-page">
            <div className="assistant-page-inner">
                <header className="assistant-page-head">
                    <div className="assistant-page-title">
                        <h2>AI Assistant</h2>
                        <span className="llm-experimental-tag">Experimental</span>
                    </div>
                    <AssistantConnectionBar
                        connection={assistant.connection}
                        settings={assistant.settings}
                        onSave={assistant.saveSettings}
                        showSettings={assistant.showSettings}
                        onToggleSettings={() => assistant.setShowSettings((value) => !value)}
                    />
                </header>

                {assistant.showSettings && (
                    <ConnectionSettings
                        settings={assistant.settings}
                        connection={assistant.connection}
                        onSave={assistant.saveSettings}
                        onTest={assistant.runTest}
                    />
                )}

                {hasRun ? (
                    <div className="plot-card assistant-surface">
                        <ChatView
                            context={assistant.context}
                            settings={assistant.settings}
                            connected={assistant.connected}
                        />
                    </div>
                ) : (
                    <div className="assistant-empty">Open a run folder to use the AI Assistant.</div>
                )}
            </div>
        </section>
    );
};

export default AssistantPage;
