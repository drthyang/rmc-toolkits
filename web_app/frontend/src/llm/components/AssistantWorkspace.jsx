import React from 'react';
import { useAssistant } from '../useAssistant';
import AssistantConnectionBar from './AssistantConnectionBar';
import AssistantViews from './AssistantViews';
import ConnectionSettings from './ConnectionSettings';
import './AssistantPanel.css';

// The assistant's inner UI for the dashboard card (AssistantPanel): the
// connection bar, the optional settings drawer, and the feature tabs/views. The
// dedicated page (AssistantPage) composes the same pieces with the connection
// bar lifted into the page header instead.
const AssistantWorkspace = (props) => {
    const assistant = useAssistant(props);

    return (
        <div className="llm-assistant-body">
            <AssistantConnectionBar
                connection={assistant.connection}
                settings={assistant.settings}
                onSave={assistant.saveSettings}
                showSettings={assistant.showSettings}
                onToggleSettings={() => assistant.setShowSettings((value) => !value)}
            />
            {assistant.showSettings && (
                <ConnectionSettings
                    settings={assistant.settings}
                    connection={assistant.connection}
                    onSave={assistant.saveSettings}
                    onTest={assistant.runTest}
                />
            )}
            <AssistantViews
                activeTab={assistant.activeTab}
                onTabChange={assistant.setActiveTab}
                context={assistant.context}
                settings={assistant.settings}
                connected={assistant.connected}
            />
        </div>
    );
};

export default AssistantWorkspace;
