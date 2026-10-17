// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { SYSTEM_PROMPT } from './system';
import { contextToJson } from '../context/runContext';

// Messages for one Copilot turn. The domain rules (SYSTEM_PROMPT) are shared
// with the watchdog; the Copilot adds how to work with tools and the verdict
// line its page turns into the Outcome chip (copilot/outcome.js). The run
// context still travels inside a fenced JSON block, so the model can tell
// data from instructions, and tool results arrive as `tool` messages.

// The Outcome line explains; the app computes the verdict from its checks
// (copilot/verdict.js), so the model is not asked to grade itself.
export const OUTCOME_RULE = [
    'Finish every answer with one final line, exactly:',
    '`Outcome: <one sentence>` — what the evidence shows, or what is missing.',
    'Do not grade your answer: the app sets the verdict (achieved, partial, not',
    'verified, not achieved) from its own checks.'
].join(' ');

// Numbers by reference (copilot/references.js) and the evidence gate
// (copilot/evidence.js), with the tools.
export const REFERENCE_RULES = [
    'Quote numbers from tool results by reference: each result has a `ref`',
    '(r1, r2, …); write {{r2.peak.significance_sigma}} — the ref, then the path of',
    'the value in that result\'s data, array items by index ({{r3.sites.0.shift_A}})',
    '— and the app writes in the engine\'s value. The run context is',
    '{{ctx.<path>}}. Type a number yourself only when you derived it, and say how.'
].join(' ');

export const EVIDENCE_RULES = [
    'Every conclusion needs its analysis in this question: order between cells',
    'or a split site called disorder, analyze_modulations; a preferred direction,',
    'analyze_displacement_directions; clusters, analyze_clusters; a coordination',
    'number or split shell, find_gr_peaks; a lowering of the average structure,',
    'measure_symmetry_lowering; "no order", every one of them. The app sends an',
    'answer back when its analysis did not run.'
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
    'default. A request that names an element but not which neighbours,',
    'triplet or site (e.g. "the bond angles around Ta") is such a choice: ask,',
    'offering the likely ones. Settle what the data can first: e.g. run',
    'find_gr_peaks, then offer its shells as the options. Never ask when the',
    'question is clear, when the context or an analysis answers it, or for',
    'permission to run an analysis. One question at a time, at most two per',
    'request. Then use exactly the option the user chose and say in the answer',
    'which one you used; when you choose without asking, say which choice you',
    'made.'
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

export const COPILOT_PROMPT = `${SYSTEM_PROMPT} ${TOOL_RULES} ${REFERENCE_RULES} ${EVIDENCE_RULES} ${COMPUTE_ONLY} ${OUTCOME_RULE}`;
// With the view tools: the same rules, pointing at the page instead of the
// compute-only sentence.
export const COPILOT_PROMPT_VIEW = `${SYSTEM_PROMPT} ${TOOL_RULES} ${REFERENCE_RULES} ${EVIDENCE_RULES} ${VIEW_RULES} ${OUTCOME_RULE}`;
// For a model that cannot call tools, or with the tools switched off.
export const COPILOT_PROMPT_NO_TOOLS = `${SYSTEM_PROMPT} ${OUTCOME_RULE}`;

export const contextBlock = (context) => (
    `Context for this RMC modeling run:\n\`\`\`json\n${contextToJson(context)}\n\`\`\``
);

/**
 * The run context as the model reads it (contextBlock trims it to its
 * budget: a long ladder to its two ends, a long history downsampled), for
 * {{ctx.…}} references: a reference names what the model saw, never a value
 * the trim removed.
 */
export const shownContext = (context) => (context ? JSON.parse(contextToJson(context)) : null);

// Keep only the newest turns so long conversations stay inside small local
// models' context windows; the run context is re-sent every call anyway.
// Earlier turns travel as their question and final answer only — never their
// tool exchanges, which are re-run when needed.
export const CHAT_HISTORY_TURNS = 8;

// With skills (copilot/skills.js): their names and one-line descriptions
// only; the agent reads a skill's text with read_skill when it needs it.
export const skillsBlock = (index) => [
    'Skills — your methods for these questions. Before an analysis a skill covers, read it with read_skill and',
    'follow its steps and its checklist. Where a draft skill says "to confirm with the user", say in the answer',
    'which choice you made, or ask with ask_user when the choice changes the analysis.',
    `Your skills:\n${index}\n`
].join(' ');

/**
 * The system prompt for the tools the turn offers: the tool rules, the ask
 * rules when ask_user is among them, the skills' index when read_skill is,
 * then pointing (view tools) or the compute-only sentence, then the Outcome
 * rule.
 */
export const copilotSystemPrompt = ({ toolsEnabled = true, viewTools = false, askUser = false, skills = '' } = {}) => {
    if (!toolsEnabled) return COPILOT_PROMPT_NO_TOOLS;
    return [
        SYSTEM_PROMPT,
        TOOL_RULES,
        REFERENCE_RULES,
        EVIDENCE_RULES,
        askUser ? ASK_RULES : null,
        skills ? skillsBlock(skills) : null,
        viewTools ? VIEW_RULES : COMPUTE_ONLY,
        OUTCOME_RULE
    ].filter(Boolean).join(' ');
};

export const buildCopilotMessages = ({
    context, history = [], userText, toolsEnabled = true, viewTools = false, askUser = false, skills = ''
}) => [
    {
        role: 'system',
        content: copilotSystemPrompt({ toolsEnabled, viewTools, askUser, skills })
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
        + 'Answer now from the results above (quote their numbers by reference), say what is still missing, '
        + 'and end with the Outcome line.'
});
