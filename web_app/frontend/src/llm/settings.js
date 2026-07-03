import { useSyncExternalStore } from 'react';

// Assistant settings live in localStorage so they survive reloads without any
// backend. A tiny external store (subscribe + snapshot) keeps the panel and the
// dashboard watchdog badge in sync when either edits the settings.

const STORAGE_KEY = 'rmc-llm-settings-v1';

export const DEFAULT_SETTINGS = {
    baseUrl: 'http://localhost:11434/v1',
    model: '',
    // Sent as a Bearer token to cloud providers (OpenAI, Gemini, …). Empty for
    // local servers, which need no key. Persisted in this browser only.
    apiKey: '',
    temperature: 0.2,
    watchdogEnabled: false,
    watchdogIntervalMin: 5
};

const readStorage = () => {
    try {
        const raw = window.localStorage.getItem(STORAGE_KEY);
        const parsed = raw ? JSON.parse(raw) : null;
        return {
            ...DEFAULT_SETTINGS,
            ...(parsed && typeof parsed === 'object' ? parsed : {})
        };
    } catch {
        return { ...DEFAULT_SETTINGS };
    }
};

let snapshot = null;
const listeners = new Set();

export const loadSettings = () => {
    if (!snapshot) snapshot = readStorage();
    return snapshot;
};

export const saveSettings = (patch) => {
    snapshot = { ...loadSettings(), ...patch };
    try {
        window.localStorage.setItem(STORAGE_KEY, JSON.stringify(snapshot));
    } catch {
        // Private mode or quota errors just mean settings do not persist.
    }
    listeners.forEach((listener) => listener());
    return snapshot;
};

const subscribe = (listener) => {
    listeners.add(listener);
    return () => listeners.delete(listener);
};

export const useLlmSettings = () => useSyncExternalStore(subscribe, loadSettings);
