import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { buildRunContext } from './context/runContext';
import { checkConnection } from './provider/client';
import { saveSettings, useLlmSettings } from './settings';

// Shared assistant state for both the page (AssistantPage) and the dashboard
// card (AssistantWorkspace): settings, the connection probe (+ auto-connect),
// the settings-drawer toggle, the active tab, and the compact run context.
// `enabled` gates the (local-only) context build and the auto-connect so an
// inactive/collapsed host does no work; nothing leaves the machine but a GET to
// the user's own localhost model server.
export const useAssistant = ({
    runName,
    plotFiles,
    rValueFile,
    structure,
    symmetry,
    runSettings,
    liveData,
    enabled = true
} = {}) => {
    const settings = useLlmSettings();
    const [showSettings, setShowSettings] = useState(false);
    const [connection, setConnection] = useState({ status: 'idle', models: [], error: null, hint: null });
    const autoTestedRef = useRef(null);

    const context = useMemo(() => (
        enabled
            ? buildRunContext({ runName, plotFiles, rValueFile, structure, symmetry, runSettings, liveData })
            : null
    ), [enabled, runName, plotFiles, rValueFile, structure, symmetry, runSettings, liveData]);

    // `manual` records whether the user actively pressed Test (vs. the automatic
    // on-load probe): a failed auto-probe is shown as a gentle "connect a model"
    // prompt, while a failed manual test surfaces the full error + hint.
    const probe = useCallback(async (manual) => {
        setConnection({ status: 'testing', models: [], error: null, hint: null, manual });
        const result = await checkConnection(settings.baseUrl, { apiKey: settings.apiKey });
        setConnection({
            status: result.ok ? 'ok' : 'error',
            models: result.models,
            error: result.error,
            hint: result.hint,
            manual
        });
        // Auto-pick a model when none is chosen (or the chosen one is gone).
        if (result.ok && result.models.length && !result.models.includes(settings.model)) {
            saveSettings({ model: result.models[0] });
        }
    }, [settings.baseUrl, settings.apiKey, settings.model]);

    const runTest = useCallback(() => probe(true), [probe]);

    // Probe the server once per base URL when active, so the status indicator
    // and model list populate without a manual click. The ref guard keeps a
    // model switch (which changes probe) from re-probing.
    useEffect(() => {
        if (!enabled) return;
        if (autoTestedRef.current === settings.baseUrl) return;
        autoTestedRef.current = settings.baseUrl;
        // One-shot probe of the user's local server — an intentional sync with
        // an external system, not React-derived state.
        // eslint-disable-next-line react-hooks/set-state-in-effect
        probe(false);
    }, [enabled, settings.baseUrl, probe]);

    const connected = connection.status === 'ok' && Boolean(settings.model);

    return {
        settings,
        saveSettings,
        connection,
        connected,
        showSettings,
        setShowSettings,
        runTest,
        context
    };
};
