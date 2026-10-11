// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { SYSTEM_PROMPT } from './system';
import { contextToJson } from '../context/runContext';

// Messages for one Copilot turn. The domain rules (SYSTEM_PROMPT) are shared
// with the watchdog; the Copilot adds how to work with tools and the verdict
// line its page turns into the Outcome chip (copilot/outcome.js). The run
// context still travels inside a fenced JSON block, so the model can tell
// data from instructions, and tool results arrive as `tool` messages.

const OUTCOME_RULE = [
    'Finish every answer with one final line, exactly:',
    '`Outcome: achieved|partial|not achieved — <one sentence>`.',
    'achieved: the question is answered by numbers you quoted; partial: some',
    'evidence is missing or a check warned; not achieved: the run or the tools',
    'could not provide the evidence. Say why in the sentence.'
].join(' ');

const TOOL_RULES = [
    'You are the Copilot: you can run the Workbench analyses with the tools',
    'provided, and you should, whenever the question needs numbers the context',
    'does not hold. Do not call a tool when the context already answers.',
    'Work in steps:',
    '(1) decide what evidence the question needs;',
    '(2) call the tools, saying in one short sentence why before each call;',
    '(3) read each result and its `checks` — `fail` means the result does not',
    'answer as asked: change the arguments (e.g. widen a bond window, pick',
    'another site) or the approach, and never repeat an identical call;',
    '`warn` means use the result with the stated caveat;',
    '(4) answer, quoting the numbers from the tool results.',
    'peakSignificance, mapSignificance and anisotropy significances are',
    'calibrated one-sided sigma; ≥ 3 is significant.',
    'A g(r) peak position or width cannot tell a distortion from thermal',
    'spread: for a cluster or polyhedron distortion (Jahn–Teller, dimers,',
    'trimers) run analyze_clusters or compute_bond_angles, and for where the',
    'fit fails run compare_fit. A z-score with a tiny change is not a finding:',
    'quote both.'
].join(' ');

const COMPUTE_ONLY = 'Tools only compute: they never change what the user sees.';

// With ask_user (copilot/askUser.js): ask instead of guessing a parameter that
// changes the science, and only then.
export const ASK_RULES = [
    'Ask instead of guessing a parameter that changes the science. When the',
    'question leaves open a choice that changes which analysis you run or its',
    'parameters — which coordination shell or r-window, which site or element',
    'pair, which symmetry tolerance, which run or dataset — and neither the',
    'context nor an analysis can settle it, call ask_user with 2–4 short',
    'options, each saying what it means (the window in Å, the site), and your',
    'default. Settle what the data can first: e.g. run find_gr_peaks, then offer',
    'its shells as the options. Never ask when the question is clear, when the',
    'context or an analysis answers it, or for permission to run an analysis.',
    'One question at a time, at most two per request. Then use exactly the',
    'option the user chose and say in the answer which one you used.'
].join(' ');

// With the view tools (highlight_cards from the host's cards, and the host's
// own orient_3d_view / current_view): explain by pointing at the page.
export const VIEW_RULES = [
    'You can also point at the Workbench the user is looking at. Its results sit',
    'on cards (the ids are listed in highlight_cards). When you explain an',
    'observation, show it: call highlight_cards with the card first (on a plot,',
    'with marks at the x positions you name, e.g. the g(r) peaks), then describe',
    'what the card shows. In the answer, cite a card as a link',
    '[short label](card:<id>), e.g. [the χ² history](card:chi2-history); the user',
    'clicks it to see the card. For directions in 3D, call orient_3d_view to look',
    'down the crystal direction you describe (e.g. along [111] on the',
    'displacement sphere) or to centre and highlight the sites you name, before',
    'you describe them. When the user says "this", "here" or "what I see", call',
    'current_view. The analysis tools only compute; the view tools change only',
    'what is on screen (never data or a file). One or two view calls per answer',
    'is enough: point at what matters.'
].join(' ');

export const COPILOT_PROMPT = `${SYSTEM_PROMPT} ${TOOL_RULES} ${COMPUTE_ONLY} ${OUTCOME_RULE}`;
// With the view tools: the same rules, pointing at the page instead of the
// compute-only sentence.
export const COPILOT_PROMPT_VIEW = `${SYSTEM_PROMPT} ${TOOL_RULES} ${VIEW_RULES} ${OUTCOME_RULE}`;
// For a model that cannot call tools, or with the tools switched off.
export const COPILOT_PROMPT_NO_TOOLS = `${SYSTEM_PROMPT} ${OUTCOME_RULE}`;

export const contextBlock = (context) => (
    `Context for this RMC modeling run:\n\`\`\`json\n${contextToJson(context)}\n\`\`\``
);

// Keep only the newest turns so long conversations stay inside small local
// models' context windows; the run context is re-sent every call anyway.
// Earlier turns travel as their question and final answer only — never their
// tool exchanges, which are re-run when needed.
export const CHAT_HISTORY_TURNS = 8;

/**
 * The system prompt for the tools the turn offers: the tool rules, the ask
 * rules when ask_user is among them, then pointing (view tools) or the
 * compute-only sentence, then the Outcome rule.
 */
export const copilotSystemPrompt = ({ toolsEnabled = true, viewTools = false, askUser = false } = {}) => {
    if (!toolsEnabled) return COPILOT_PROMPT_NO_TOOLS;
    if (!askUser) return viewTools ? COPILOT_PROMPT_VIEW : COPILOT_PROMPT;
    return `${SYSTEM_PROMPT} ${TOOL_RULES} ${ASK_RULES} ${viewTools ? VIEW_RULES : COMPUTE_ONLY} ${OUTCOME_RULE}`;
};

export const buildCopilotMessages = ({
    context, history = [], userText, toolsEnabled = true, viewTools = false, askUser = false
}) => [
    {
        role: 'system',
        content: copilotSystemPrompt({ toolsEnabled, viewTools, askUser })
    },
    { role: 'user', content: contextBlock(context) },
    {
        role: 'assistant',
        content: 'Understood. I will answer questions about this RMC modeling run from the context above'
            + `${toolsEnabled ? ' and the analyses I run' : ''}.`
    },
    ...(history || []).slice(-CHAT_HISTORY_TURNS).map(({ role, content }) => ({ role, content })),
    { role: 'user', content: userText }
];

// Appended when the turn has used its tool rounds: the next (last) request
// sends tool_choice 'none', and this note tells the model why.
export const toolBudgetNote = (rounds) => ({
    role: 'user',
    content: `You have used the analysis budget for this question (${rounds} rounds of tool calls). `
        + 'Answer now from the results above, say what is still missing, and end with the Outcome line.'
});
