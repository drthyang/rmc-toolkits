import React from 'react';
import ChatView from './ChatView';
import ReportView from './ReportView';
import RunSummaryView from './RunSummaryView';

// The feature tabs and the active view. Shared by the page and the card.

const TABS = [
    { id: 'chat', label: 'Chat' },
    { id: 'summary', label: 'Summary' },
    { id: 'report', label: 'Report' }
];

const AssistantViews = ({ activeTab, onTabChange, context, settings, connected }) => (
    <>
        <nav className="llm-tabs" aria-label="Assistant features">
            {TABS.map((tab) => (
                <button
                    key={tab.id}
                    type="button"
                    className={activeTab === tab.id ? 'active' : ''}
                    onClick={() => onTabChange(tab.id)}
                >
                    {tab.label}
                </button>
            ))}
        </nav>
        {activeTab === 'chat' && (
            <ChatView context={context} settings={settings} connected={connected} />
        )}
        {activeTab === 'summary' && (
            <RunSummaryView context={context} settings={settings} connected={connected} />
        )}
        {activeTab === 'report' && (
            <ReportView context={context} settings={settings} connected={connected} />
        )}
    </>
);

export default AssistantViews;
