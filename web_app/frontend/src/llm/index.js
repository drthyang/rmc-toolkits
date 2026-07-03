// Public surface of the experimental LLM assistant module. The host app should
// import only from here; everything else in src/llm/ is internal. See
// src/llm/README.md for the architecture and the extraction recipe.
export { default as AssistantPanel } from './components/AssistantPanel';
export { default as WatchdogBadge } from './components/WatchdogBadge';
