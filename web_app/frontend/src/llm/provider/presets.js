// OpenAI-compatible servers this module can talk to. Local servers (Ollama, LM
// Studio) expose GET /models and POST /chat/completions under a /v1 prefix with
// no API key, so the browser calls them directly and run data never leaves the
// machine. Cloud providers speak the same dialect but need a Bearer API key —
// and, unlike the local ones, they receive the run-derived context, so they are
// flagged `cloud` and surfaced with a data-leaves-your-device warning.
export const PROVIDER_PRESETS = [
    {
        id: 'ollama',
        label: 'Ollama',
        baseUrl: 'http://localhost:11434/v1',
        cloud: false,
        corsHint: 'Allow this site before starting the server: OLLAMA_ORIGINS="https://drthyang.github.io" ollama serve (or OLLAMA_ORIGINS="*").'
    },
    {
        id: 'lmstudio',
        label: 'LM Studio',
        baseUrl: 'http://localhost:1234/v1',
        cloud: false,
        corsHint: 'Enable CORS in LM Studio: Developer tab → server settings → Enable CORS, then start the local server.'
    },
    {
        id: 'openai',
        label: 'OpenAI',
        baseUrl: 'https://api.openai.com/v1',
        cloud: true,
        keyUrl: 'https://platform.openai.com/api-keys',
        corsHint: 'Needs an OpenAI API key. Your run data is sent to OpenAI — it does not stay on your device.'
    },
    {
        id: 'gemini',
        label: 'Gemini',
        baseUrl: 'https://generativelanguage.googleapis.com/v1beta/openai',
        cloud: true,
        keyUrl: 'https://aistudio.google.com/apikey',
        corsHint: 'Needs a Google AI Studio (Gemini) API key. Your run data is sent to Google — it does not stay on your device.'
    }
];

export const providerForUrl = (baseUrl) => (
    PROVIDER_PRESETS.find((preset) => preset.baseUrl === baseUrl) || null
);

const LOCAL_HOSTNAMES = new Set(['localhost', '127.0.0.1', '[::1]', '::1', '0.0.0.0']);

// A base URL is treated as local (private, no key) when its host is loopback.
// Anything else is a remote/cloud endpoint that should collect an API key.
export const isLocalUrl = (baseUrl) => {
    try {
        return LOCAL_HOSTNAMES.has(new URL(baseUrl).hostname);
    } catch {
        return false;
    }
};
