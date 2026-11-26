// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What each Copilot tool does to the Workbench: its effect. One list,
// COPILOT_TOOL_EFFECTS, names every tool the Copilot can be offered — the
// turn's own, the tools over the run props, the module's view tool and the
// host's engine, view and change tools — with its effect, one of
// TOOL_EFFECTS:
//
//   'compute'  computes or reads and changes nothing the user sees: the
//              analyses, the run data, a skill, a question to the user;
//   'view'     changes only what is on screen (a card's edge, a page, a 3D
//              camera, a site selected as a click would): never data, an
//              analysis input or a file;
//   'change'   changes an analysis input the user can set on a page (the Bond
//              Geometry inputs). It asks the user first on an approval card
//              unless they chose Auto, and can be undone.
//
// The dock and the validation harness read their tools through this list
// (toolsWithEffects): the dock offers every effect, the live arms 'compute'
// only. A change tool also needs someone to approve it: the turn offers one
// only with an approve port (turn.js), which the harness has none of. The
// list holds the host's tool names as the verification's catalogues do
// (claims.js, facts.js). A tool it does not name (a test's) falls back on its
// own declaration: `effect`, else `view: true` → 'view', else 'compute'.
// src/__tests__/toolEffects.test.js checks that every tool the builders make
// is named here, that `view: true` is set exactly on the 'view' tools, and
// what the dock and the harness offer.

export const TOOL_EFFECTS = Object.freeze(['compute', 'view', 'change']);

export const COPILOT_TOOL_EFFECTS = Object.freeze({
    // The turn's own (askUser.js, skills.js): a question to the user, a skill read.
    ask_user: 'compute',
    read_skill: 'compute',
    // The tools over the run props (runTools.js).
    read_run_data: 'compute',
    assess_convergence: 'compute',
    find_gr_peaks: 'compute',
    compare_fit: 'compute',
    // The module's view tool over the host's cards (viewTools.js).
    highlight_cards: 'view',
    // The host's engine tools (src/copilotTools.js).
    list_sites: 'compute',
    analyze_displacement_directions: 'compute',
    compute_bond_angles: 'compute',
    polyhedral_distortion: 'compute',
    analyze_clusters: 'compute',
    analyze_modulations: 'compute',
    detect_symmetry: 'compute',
    measure_symmetry_lowering: 'compute',
    // Termination ripples (src/copilotRipples.js), offered in the dock only for
    // now: validation/harness.js LIVE_ARM_EXCLUDED leaves it out of the arms.
    check_ripples: 'compute',
    // The host's view tools (src/copilotViewTools.js).
    orient_3d_view: 'view',
    current_view: 'view',
    // The host's change tool (src/copilotChangeTools.js).
    set_bond_geometry: 'change'
});

/** What the dock offers: every effect (components/CopilotChat.jsx). */
export const DOCK_EFFECTS = TOOL_EFFECTS;

/** A tool's effect: the list's, else its own declaration (see the header). */
export const effectOf = (tool) => {
    const listed = Object.hasOwn(COPILOT_TOOL_EFFECTS, tool?.name ?? '') ? COPILOT_TOOL_EFFECTS[tool.name] : null;
    if (listed) return listed;
    if (TOOL_EFFECTS.includes(tool?.effect)) return tool.effect;
    return tool?.view ? 'view' : 'compute';
};

/** Whether a tool changes an analysis input (an approval card first, undoable). */
export const isChangeTool = (tool) => effectOf(tool) === 'change';

/**
 * The tools whose effect is one of `effects`, in their order (falsy entries
 * dropped): what an adapter offers from the tools it built. A TypeError for
 * an effect that is not one of TOOL_EFFECTS.
 */
export const toolsWithEffects = (tools, effects) => {
    if (!Array.isArray(effects) || !effects.length || effects.some((effect) => !TOOL_EFFECTS.includes(effect))) {
        throw new TypeError(`toolsWithEffects: effects must be a non-empty list of ${TOOL_EFFECTS.join(', ')}`);
    }
    return (Array.isArray(tools) ? tools : []).filter((tool) => tool && effects.includes(effectOf(tool)));
};
