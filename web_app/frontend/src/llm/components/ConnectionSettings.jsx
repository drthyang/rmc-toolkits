import React, { useState } from 'react';
import { PROVIDER_PRESETS } from '../provider/presets';

// Connection + model settings drawer: provider presets, base URL, connection
// test with actionable failure hints, model picker (from GET /models),
// temperature, and the watchdog toggle. All persisted via saveSettings.

const ConnectionSettings = ({ settings, connection, onSave, onTest }) => {
    const [showGuide, setShowGuide] = useState(false);
    const testing = connection.status === 'testing';

    return (
        <div className="llm-settings">
            <div className="llm-settings-row">
                <label htmlFor="llm-base-url">Server</label>
                <div className="llm-settings-controls">
                    {PROVIDER_PRESETS.map((preset) => (
                        <button
                            key={preset.id}
                            type="button"
                            className={`llm-preset${settings.baseUrl === preset.baseUrl ? ' is-active' : ''}`}
                            onClick={() => onSave({ baseUrl: preset.baseUrl })}
                        >
                            {preset.label}
                        </button>
                    ))}
                    <input
                        id="llm-base-url"
                        type="text"
                        value={settings.baseUrl}
                        spellCheck="false"
                        onChange={(event) => onSave({ baseUrl: event.target.value })}
                    />
                    <button type="button" onClick={onTest} disabled={testing}>
                        {testing ? 'Testing…' : 'Test connection'}
                    </button>
                </div>
            </div>

            {connection.status === 'ok' && (
                <div className="llm-connection-status is-ok" role="status">
                    Connected — {connection.models.length} model{connection.models.length === 1 ? '' : 's'} available.
                </div>
            )}
            {connection.status === 'error' && (
                <div className="llm-connection-status is-error" role="alert">
                    <strong>{connection.error}</strong>
                    {connection.hint && <p>{connection.hint}</p>}
                </div>
            )}

            <div className="llm-settings-row">
                <label htmlFor="llm-model">Model</label>
                <div className="llm-settings-controls">
                    {connection.models.length ? (
                        <select
                            id="llm-model"
                            value={settings.model}
                            onChange={(event) => onSave({ model: event.target.value })}
                        >
                            <option value="">Select a model…</option>
                            {connection.models.map((model) => (
                                <option key={model} value={model}>{model}</option>
                            ))}
                        </select>
                    ) : (
                        <input
                            id="llm-model"
                            type="text"
                            placeholder="e.g. llama3.2"
                            value={settings.model}
                            spellCheck="false"
                            onChange={(event) => onSave({ model: event.target.value })}
                        />
                    )}
                    <label className="llm-inline-label" htmlFor="llm-temperature">Temperature</label>
                    <input
                        id="llm-temperature"
                        type="number"
                        min="0"
                        max="2"
                        step="0.1"
                        value={settings.temperature}
                        onChange={(event) => onSave({ temperature: Number(event.target.value) })}
                    />
                </div>
            </div>

            <div className="llm-settings-row">
                <label htmlFor="llm-watchdog">Watchdog</label>
                <div className="llm-settings-controls">
                    <label className="llm-checkbox">
                        <input
                            id="llm-watchdog"
                            type="checkbox"
                            checked={settings.watchdogEnabled}
                            onChange={(event) => onSave({ watchdogEnabled: event.target.checked })}
                        />
                        Flag stalled or diverging refinements on the dashboard
                    </label>
                    <label className="llm-inline-label" htmlFor="llm-watchdog-interval">LLM check every</label>
                    <input
                        id="llm-watchdog-interval"
                        type="number"
                        min="1"
                        max="120"
                        step="1"
                        value={settings.watchdogIntervalMin}
                        onChange={(event) => onSave({ watchdogIntervalMin: Math.max(1, Number(event.target.value) || 1) })}
                    />
                    <span className="llm-inline-label">min</span>
                </div>
            </div>

            <button
                type="button"
                className="llm-guide-toggle"
                onClick={() => setShowGuide((value) => !value)}
                aria-expanded={showGuide}
            >
                {showGuide ? 'Hide setup guide' : 'How to set up a local LLM'}
            </button>
            {showGuide && (
                <div className="llm-setup-guide">
                    <p>
                        The assistant talks straight from your browser to an OpenAI-compatible server running
                        on <strong>your own machine</strong> — run-derived summaries never leave your device
                        except to that local server.
                    </p>
                    <ol>
                        <li>
                            Install <a href="https://ollama.com" target="_blank" rel="noreferrer">Ollama</a> or{' '}
                            <a href="https://lmstudio.ai" target="_blank" rel="noreferrer">LM Studio</a> and pull a
                            model (e.g. <code>ollama pull llama3.2</code>).
                        </li>
                        {PROVIDER_PRESETS.map((preset) => (
                            <li key={preset.id}><strong>{preset.label}:</strong> {preset.corsHint}</li>
                        ))}
                        <li>Pick the matching preset above and press <em>Test connection</em>.</li>
                    </ol>
                    <p>
                        The hosted app is HTTPS but browsers allow requests to <code>http://localhost</code> as a
                        trustworthy origin in Chrome, Edge, and Firefox; Safari may block them — if the test keeps
                        failing there, run the app locally instead.
                    </p>
                </div>
            )}
        </div>
    );
};

export default ConnectionSettings;
