// Local OpenAI-compatible servers this module is designed to talk to. Both
// expose GET /models and POST /chat/completions under a /v1 prefix and need no
// API key, so the browser can call them directly and run data never leaves the
// machine the model runs on.
export const PROVIDER_PRESETS = [
    {
        id: 'ollama',
        label: 'Ollama',
        baseUrl: 'http://localhost:11434/v1',
        corsHint: 'Allow this site before starting the server: OLLAMA_ORIGINS="https://drthyang.github.io" ollama serve (or OLLAMA_ORIGINS="*").'
    },
    {
        id: 'lmstudio',
        label: 'LM Studio',
        baseUrl: 'http://localhost:1234/v1',
        corsHint: 'Enable CORS in LM Studio: Developer tab → server settings → Enable CORS, then start the local server.'
    }
];
