// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Public surface of the experimental LLM module (the AI Copilot and the
// convergence watchdog badge). The host app should import only from here;
// everything else in src/llm/ is internal. See src/llm/README.md for the
// architecture, the tool contract and the extraction recipe.
export { default as CopilotChat } from './components/CopilotChat';
export { default as WatchdogBadge } from './components/WatchdogBadge';
// What a point group implies (inversion, chirality, polarity), for host tools
// that report a space group (copilotTools.js detect_symmetry).
export { describePointGroup, pointGroupFacts } from './context/pointGroups';
