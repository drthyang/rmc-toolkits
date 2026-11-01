// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The host's tools for the AI Copilot dock, as App.jsx passes them (the
// `tools` prop of CopilotChat): the run's engine tools (the run module's,
// copilotRun.js), the view tools over the host's 3D cards
// (copilotViewTools.js) and the change tool for the Bond Geometry inputs
// (copilotChangeTools.js). One function, so the app and the tests that
// capture what the dock sends build the same list. The dock adds the tools
// over the run props and highlight_cards before them and reads them all
// through the one list of tool effects (llm/copilot/effects.js); the
// validation harness takes the engine tools only (validation/harness.js
// liveArmTools).

import { buildChangeTools } from './copilotChangeTools';
import { buildViewTools } from './copilotViewTools';

/**
 * `engineTools` and `elements` are the run's (CopilotRun's `tools` and
 * `elements`), `cards` / `pages` / `activePage` what the view tools point at,
 * `history` the app's change history (bondGeometryControl.js) and
 * `openPage(key)` a page switch for the change tool.
 */
export const copilotHostTools = ({
    engineTools = [], elements = [], cards = [], pages = {}, activePage = () => null, history, openPage = null
} = {}) => [
    ...engineTools,
    ...buildViewTools({ cards, pages, activePage }),
    ...buildChangeTools({ elements, history, openPage })
];
