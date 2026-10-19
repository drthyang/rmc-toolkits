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
// The run's one shell reading (context/shells.js): when two coordination
// numbers agree, and compute_bond_angles' counts rule for a window that cuts
// a shell, with its thresholds (copilotTools.js).
export {
    cnMatches, countsCut, MODAL_COORDINATION_SHARE, shellCatalogue, UNBONDED_APEX_FRACTION
} from './context/shells';
// A check whose text holds method advice (copilot/method.js): the turn sends
// it without the advice under method 'off', the validation's arm C
// (copilotTools.js).
export { methodCheck } from './copilot/method';
