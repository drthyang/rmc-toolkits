import React, { useState } from 'react';
import { PROVIDER_PRESETS, isLocalUrl, providerForUrl } from '../provider/presets';

// Connection + model settings drawer: provider presets (local + cloud), base
// URL, an API-key field and data-leaves-your-device warning for cloud
// providers, connection test with actionable hints, temperature, and the
// watchdog toggle. The model itself is picked from the toolbar.

const CheckIcon = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.6" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <path d="M5 12.5l4.5 4.5L19 6.5" />
    </svg>
);

const HelpIcon = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.9" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <circle cx="12" cy="12" r="9.2" />
        <path d="M9.6 9.4a2.4 2.4 0 0 1 4.5 1.1c0 1.6-2.1 2-2.1 2.4" />
        <path d="M12 16.4h.01" />
    </svg>
);

const WarnIcon = () => (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.9" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true">
        <path d="M10.3 3.9 1.8 18a2 2 0 0 0 1.7 3h17a2 2 0 0 0 1.7-3L13.7 3.9a2 2 0 0 0-3.4 0z" />
        <path d="M12 9v4" />
        <path d="M12 17h.01" />
    </svg>
);

const hostLabel = (url) => {
    try {
        return new URL(url).host;
    } catch {
        return url;
    }
};

const ConnectionSettings = ({ settings, connection, onSave, onTest }) => {
    const [showGuide, setShowGuide] = useState(false);
    const testing = connection.status === 'testing';

    const activeProvider = providerForUrl(settings.baseUrl);
    const isCloud = activeProvider ? activeProvider.cloud : !isLocalUrl(settings.baseUrl);
    const providerLabel = activeProvider?.label || 'this remote server';
    // This page's origin — the value a local server must allow via CORS.
    const pageOrigin = typeof window !== 'undefined' ? window.location.origin : 'https://drthyang.github.io';

    return (
        <div className="llm-settings">
            <div className="llm-settings-head">
                <span className="llm-settings-title">Connection</span>
                <span className="llm-settings-note">
                    Local models (Ollama, LM Studio) keep your run data on your device. Cloud providers send it to
                    their API using your key.
                </span>
            </div>

            <div className="llm-field">
                <span className="llm-field-label">Provider</span>
                <div className="llm-preset-grid">
                    {PROVIDER_PRESETS.map((preset) => {
                        const active = settings.baseUrl === preset.baseUrl;
                        return (
                            <button
                                key={preset.id}
                                type="button"
                                className={`llm-preset${active ? ' is-active' : ''}`}
                                onClick={() => onSave({ baseUrl: preset.baseUrl })}
                                aria-pressed={active}
                            >
                                <span className="llm-preset-text">
                                    <span className="llm-preset-name">
                                        {preset.label}
                                        {preset.cloud && <span className="llm-preset-cloud">Cloud</span>}
                                    </span>
                                    <span className="llm-preset-host">{hostLabel(preset.baseUrl)}</span>
                                </span>
                                {active && <span className="llm-preset-check"><CheckIcon /></span>}
                            </button>
                        );
                    })}
                </div>
            </div>

            {isCloud && (
                <div className="llm-cloud-warning" role="note">
                    <WarnIcon />
                    <span>
                        Sends your run data to <strong>{providerLabel}</strong> — it leaves your device. Use a local
                        provider to keep everything on your machine.
                    </span>
                </div>
            )}

            <div className="llm-field">
                <label className="llm-field-label" htmlFor="llm-base-url">Server URL</label>
                <div className="llm-url-row">
                    <input
                        id="llm-base-url"
                        className="llm-input llm-url-input"
                        type="text"
                        value={settings.baseUrl}
                        spellCheck="false"
                        onChange={(event) => onSave({ baseUrl: event.target.value })}
                    />
                    <button type="button" className="llm-primary llm-test-btn" onClick={onTest} disabled={testing}>
                        {testing ? 'Testing…' : 'Test'}
                    </button>
                </div>
            </div>

            {isCloud && (
                <div className="llm-field">
                    <label className="llm-field-label" htmlFor="llm-api-key">API key</label>
                    <input
                        id="llm-api-key"
                        className="llm-input llm-key-input"
                        type="password"
                        placeholder="Paste your API key"
                        value={settings.apiKey}
                        spellCheck="false"
                        autoComplete="off"
                        onChange={(event) => onSave({ apiKey: event.target.value })}
                    />
                    <span className="llm-field-hint">
                        {activeProvider?.keyUrl && (
                            <>Get a key at <a href={activeProvider.keyUrl} target="_blank" rel="noreferrer">{hostLabel(activeProvider.keyUrl)}</a>. </>
                        )}
                        Stored only in this browser.
                    </span>
                </div>
            )}

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

            <div className="llm-settings-divider" />

            <div className="llm-field">
                <label className="llm-field-label" htmlFor="llm-temperature">Temperature</label>
                <div className="llm-inline-row">
                    <input
                        id="llm-temperature"
                        className="llm-input llm-num"
                        type="number"
                        min="0"
                        max="2"
                        step="0.1"
                        value={settings.temperature}
                        onChange={(event) => onSave({ temperature: Number(event.target.value) })}
                    />
                    <span className="llm-field-hint">Lower is more factual — 0.2 is a good default.</span>
                </div>
            </div>

            <div className="llm-field">
                <span className="llm-field-label">Watchdog</span>
                <label className="llm-check" htmlFor="llm-watchdog">
                    <input
                        id="llm-watchdog"
                        type="checkbox"
                        checked={settings.watchdogEnabled}
                        onChange={(event) => onSave({ watchdogEnabled: event.target.checked })}
                    />
                    <span>Flag stalled or diverging refinements on the dashboard.</span>
                </label>
                <div className="llm-watch-row">
                    <span className="llm-field-hint">Re-check with the model every</span>
                    <input
                        className="llm-input llm-num"
                        type="number"
                        min="1"
                        max="120"
                        step="1"
                        value={settings.watchdogIntervalMin}
                        onChange={(event) => onSave({ watchdogIntervalMin: Math.max(1, Number(event.target.value) || 1) })}
                        aria-label="Watchdog interval in minutes"
                    />
                    <span className="llm-field-hint">min</span>
                </div>
            </div>

            <button
                type="button"
                className="llm-guide-toggle"
                onClick={() => setShowGuide((value) => !value)}
                aria-expanded={showGuide}
            >
                <HelpIcon />
                {showGuide ? 'Hide setup guide' : 'How to connect a model'}
            </button>
            {showGuide && (
                <div className="llm-setup-guide">
                    <p>
                        <strong>Local (private):</strong> the assistant talks straight from your browser to an
                        OpenAI-compatible server on <strong>your own machine</strong> — run-derived summaries never
                        leave your device.
                    </p>
                    <ol>
                        <li>
                            Install <a href="https://ollama.com" target="_blank" rel="noreferrer">Ollama</a> or{' '}
                            <a href="https://lmstudio.ai" target="_blank" rel="noreferrer">LM Studio</a> and pull a
                            model (e.g. <code>ollama pull llama3.2</code>).
                        </li>
                        {PROVIDER_PRESETS.filter((preset) => !preset.cloud).map((preset) => (
                            <li key={preset.id}><strong>{preset.label}:</strong> {preset.corsHint}</li>
                        ))}
                        <li>Pick the matching provider above and press <em>Test</em>.</li>
                    </ol>
                    <p>
                        <strong>Cloud:</strong> pick OpenAI or Gemini and paste an API key — no server to run, but your
                        run context is sent to that provider. Some providers or networks may block direct browser
                        calls (CORS); if the test fails, run the app locally.
                    </p>
                    <p>
                        <strong>Hosted app + Ollama:</strong> Ollama must allow this page's origin
                        (<code>{pageOrigin}</code>). Quit Ollama, then start it with the origin allowed:{' '}
                        <code>OLLAMA_ORIGINS="{pageOrigin}" ollama serve</code>. On the macOS menu-bar app run{' '}
                        <code>launchctl setenv OLLAMA_ORIGINS "{pageOrigin}"</code> and reopen Ollama (this resets on
                        logout).
                    </p>
                    <p>
                        <strong>Safari</strong> blocks HTTPS pages from calling <code>http://localhost</code> — Chrome,
                        Edge, and Firefox allow it. If the test fails only in Safari, use one of those, or run the app
                        locally over http (<code>VITE_STATIC_MODE=true npm run dev</code>) so the page and the model
                        share the localhost origin.
                    </p>
                </div>
            )}
        </div>
    );
};

export default ConnectionSettings;
