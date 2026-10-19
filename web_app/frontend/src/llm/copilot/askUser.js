// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// ask_user: the Copilot asks the user one short clarifying question instead
// of guessing a parameter that changes the science — which coordination shell
// or r-window, which site or pair, which tolerance, which run. It renders as a
// question card in the conversation (components/QuestionCard.jsx) with 2–4
// option buttons and an "Other…" field; the turn pauses until the user
// answers, and the answer returns to the model as this tool's result. An
// interactive tool (`interactive: true`): never answered from the call cache,
// no time limit while it waits, and at most MAX_QUESTIONS per question from the
// user (agentLoop.js). It shows nothing on the page and changes no data.
//
// The conversation supplies `askUser(request, { signal })` in the tools'
// context (useCopilot.js); it resolves with the user's reply
// `{ answer, source: 'option' | 'other' | 'skipped', index }` and rejects with
// an AbortError when the turn is stopped or the conversation restarts. The
// tool returns the reply beside its result as `choice` (replyChoice), which
// the tool contract keeps on the step (tools.js): the dock's card and the
// validation's record read whether and how the user answered from it, never
// from the text the model reads.

export const ASK_USER = 'ask_user';
// Per question from the user, and one per model step.
export const MAX_QUESTIONS = 2;
export const MIN_OPTIONS = 2;
export const MAX_OPTIONS = 4;
const MAX_QUESTION_CHARS = 300;
const MAX_LABEL_CHARS = 80;
const MAX_DESCRIPTION_CHARS = 200;

const clip = (text, limit) => (text.length > limit ? `${text.slice(0, limit - 1).trimEnd()}…` : text);

// A label and an optional one-line description from one option. A model may
// write the description into the label ("Both peaks — 2.33–2.96 Å"): split it
// at a spaced dash only when no description was given, and never inside a
// range ("2.2 – 2.7 Å" is one window).
const readOption = (option) => {
    const raw = typeof option === 'string' ? { label: option } : (option || {});
    let label = String(raw.label ?? '').replace(/\s+/g, ' ').trim();
    let description = String(raw.description ?? '').replace(/\s+/g, ' ').trim();
    if (!description) {
        for (const dash of label.matchAll(/\s+[—–-]\s+/g)) {
            const head = label.slice(0, dash.index);
            const tail = label.slice(dash.index + dash[0].length);
            if (/\d$/.test(head) && /^[-−]?\d/.test(tail)) continue;
            if (head.length >= 2 && head.length <= MAX_LABEL_CHARS && tail) [label, description] = [head, tail];
            break;
        }
    }
    return { label: clip(label, MAX_LABEL_CHARS), ...(description ? { description: clip(description, MAX_DESCRIPTION_CHARS) } : {}) };
};

const sameLabel = (a, b) => a.toLocaleLowerCase() === b.toLocaleLowerCase();

/**
 * The validated question the card shows: `{ question, options, allowOther,
 * defaultIndex }`. Throws an error the model can read and fix (no question,
 * fewer than two distinct options).
 */
export const normalizeAskRequest = (args = {}) => {
    const question = clip(String(args.question ?? '').replace(/\s+/g, ' ').trim(), MAX_QUESTION_CHARS);
    if (!question) throw new Error('question is empty: write one short question');
    const options = [];
    for (const raw of Array.isArray(args.options) ? args.options : []) {
        const option = readOption(raw);
        if (option.label && !options.some((kept) => sameLabel(kept.label, option.label))) options.push(option);
    }
    if (options.length < MIN_OPTIONS) {
        throw new Error(`give ${MIN_OPTIONS}–${MAX_OPTIONS} distinct options, each a short label`);
    }
    if (options.length > MAX_OPTIONS) {
        throw new Error(`give at most ${MAX_OPTIONS} options; drop the least likely`);
    }
    const wanted = String(args.default ?? '').trim();
    let defaultIndex = wanted ? options.findIndex((option) => sameLabel(option.label, wanted)) : -1;
    // "2" or "option 2" for the second option.
    if (defaultIndex < 0 && /^(option\s*)?[1-4]$/i.test(wanted)) {
        const index = Number(wanted.replace(/\D/g, '')) - 1;
        if (index < options.length) defaultIndex = index;
    }
    return {
        question,
        options,
        allowOther: args.allow_other !== false,
        defaultIndex: defaultIndex >= 0 ? defaultIndex : null
    };
};

/**
 * The user's reply as a step's choice: where it came from (an option, a typed
 * answer, a skip) and the option's index (null for a typed answer or a skip).
 * A reply of no known source reads as an option, as answerResult reads it.
 */
export const replyChoice = (reply) => ({
    source: reply?.source === 'other' || reply?.source === 'skipped' ? reply.source : 'option',
    index: Number.isInteger(reply?.index) ? reply.index : null
});

/** The tool result the model reads after the user replied. */
export const answerResult = (request, reply) => {
    const chosen = Number.isInteger(reply?.index) ? request.options[reply.index] : null;
    const fallback = request.defaultIndex != null ? request.options[request.defaultIndex] : null;
    if (reply?.source === 'skipped') {
        return {
            data: {
                question: request.question,
                answer: fallback ? fallback.label : null,
                skipped: true,
                ...(fallback?.description ? { means: fallback.description } : {})
            },
            checks: [{
                level: 'note',
                text: fallback
                    ? `The user skipped the question: use your default, "${fallback.label}", and say in the answer that you assumed it.`
                    : 'The user skipped the question: choose the most standard option and say in the answer what you assumed.'
            }]
        };
    }
    const answer = String(reply?.answer ?? '').trim();
    if (reply?.source === 'other') {
        return {
            data: { question: request.question, answer, typed: true },
            checks: [{ level: 'ok', text: `The user answered in their own words: "${answer}". Follow it exactly.` }]
        };
    }
    return {
        data: {
            question: request.question,
            answer: chosen?.label ?? answer,
            ...(chosen?.description ? { means: chosen.description } : {})
        },
        checks: [{
            level: 'ok',
            text: `The user chose "${chosen?.label ?? answer}"${chosen?.description ? ` (${chosen.description})` : ''}. `
                + 'Use exactly this choice and say in the answer which one you used.'
        }]
    };
};

export const buildAskUserTool = () => ({
    name: ASK_USER,
    interactive: true,
    description: 'Ask the user ONE short clarifying question and wait for the answer. Only when the answer changes '
        + 'which analysis you run or its parameters (which coordination shell or r-window, which site or pair, which '
        + 'tolerance, which run or dataset) and neither the run context nor an analysis can settle it. Offer 2–4 '
        + 'short options, each with a one-line description (e.g. the window in Å), and your default.',
    parameters: {
        type: 'object',
        properties: {
            question: { type: 'string', description: 'One short question.' },
            options: {
                type: 'array',
                minItems: MIN_OPTIONS,
                maxItems: MAX_OPTIONS,
                items: {
                    type: 'object',
                    properties: {
                        label: { type: 'string', description: 'A short choice, a few words.' },
                        description: { type: 'string', description: 'Optional: one line saying what it means, e.g. "2.33–2.96 Å, both Ta–Se peaks".' }
                    },
                    required: ['label']
                },
                description: 'The choices, the most likely first.'
            },
            allow_other: { type: 'boolean', description: 'Let the user type another answer. Default true.' },
            default: { type: 'string', description: 'The label of the option you would pick yourself.' }
        },
        required: ['question', 'options']
    },
    // answerResult: skipped (answer = the default's label, or null), typed
    // (the user's own words) or an option chosen; `means` is that option's description.
    dataShape: {
        question: 'string',
        answer: 'string|null',
        'skipped?': 'boolean',
        'typed?': 'boolean',
        'means?': 'string'
    },
    label: () => 'Question for you',
    run: async (args, { signal, askUser } = {}) => {
        const request = normalizeAskRequest(args);
        if (typeof askUser !== 'function') {
            throw new Error('the user cannot be asked here: proceed with your default and say what you assumed');
        }
        const reply = await askUser(request, { signal });
        return { ...answerResult(request, reply), choice: replyChoice(reply) };
    }
});
