// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useState } from 'react';
import AssistantWorkspace from './AssistantWorkspace';
import './AssistantPanel.css';

// Collapsible dashboard-card presentation of the assistant. The host app now
// hosts the assistant on its own page (AssistantPage); this card is kept as a
// self-contained, extractable surface (see README) that renders the same
// AssistantWorkspace body. Collapsed by default and makes zero network requests
// until the user opens it and acts.

const AssistantPanel = (props) => {
    const [expanded, setExpanded] = useState(false);

    return (
        <article className={`plot-card llm-assistant-card${expanded ? '' : ' is-collapsed'}`}>
            <div className="plot-card-header">
                <div className="llm-assistant-heading">
                    <h3>AI Assistant</h3>
                    <span className="llm-beta-tag">Beta</span>
                </div>
                <div className="plot-card-header-actions">
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
            {expanded && <AssistantWorkspace {...props} enabled={expanded} />}
        </article>
    );
};

export default AssistantPanel;
