// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The method a turn follows (turn.js `config.method`), and what it changes in
// what the model reads.
//
//   'full'  the Copilot as it ships (the dock; the validation's arm B): the
//           skills (read_skill and their index), the prompt's method rules,
//           and the method advice in the tools' descriptions and checks;
//   'off'   none of the method (the validation's arm C, plan item 1.8,
//           decided 2026-10-19): no skills, no EVIDENCE_RULES, no method
//           sentence of TOOL_RULES or ASK_RULES (prompts/copilot.js), no
//           method advice in a tool's description or in a check. The tools,
//           their refusals, the Workbench's questions, the evidence gate and
//           the answer checks still run; only the advice goes.
//
// Method advice is text that tells the model which analysis answers which
// question, where to take an analysis's parameter from (a window from
// find_gr_peaks), which window is the right one, or what to run next: what
// the skills teach. An inference in a check — what a pattern is evidence of,
// as the domain rules' inference parts (prompts/system.js GROUND_RULES) —
// is method too (decided 2026-10-19: analyze_modulations' long-wave and
// in-phase notes). A description keeps what the tool does and returns; a
// check keeps its fact, its verdict and its limits, and a refusal still
// refuses.
//
// A tool spec carries its advice apart from its description, in an optional
// `method` field: a string, sent after the description and a space, or a
// list of parts in order, each a string (after the description and a space)
// or `{ after, text }` (`text` right after the one place `after` appears in
// the description: advice inside a sentence). With method 'full' the
// description is sent with its advice put back — byte for byte the
// description the tool had before its advice was split out — and with 'off'
// alone. A check carries its advice inside its text instead, listed in its
// own `method` (methodCheck), because a check's text is read as written by
// more than the turn (the step cards, the record, arm A+'s readings): its
// text stays what the dock reads, and method 'off' sends it without the
// advice.

export const TURN_METHODS = Object.freeze(['full', 'off']);

const invalid = (message) => {
    throw new TypeError(`method: ${message}`);
};

// A tool's method parts, checked: [{ after?, text }] in order (a string part
// has no `after`).
const methodParts = (tool) => {
    const advice = tool?.method;
    if (advice == null) return [];
    const list = Array.isArray(advice) ? advice : [advice];
    return list.map((part) => {
        if (typeof part === 'string' && part) return { text: part };
        if (part && typeof part === 'object' && typeof part.after === 'string' && part.after
            && typeof part.text === 'string' && part.text) {
            return { after: part.after, text: part.text };
        }
        return invalid(`${tool?.name}: a method part is a non-empty string or { after, text } with both non-empty`);
    });
};

/**
 * The description of `tool` as the model reads it under `method` ('full' |
 * 'off'): with 'off' its description alone; with 'full' its `method` advice
 * put back (see the header). A TypeError for a malformed `method` field, or
 * an `after` that does not appear exactly once in the description.
 */
export const toolDescription = (tool, method) => {
    const parts = methodParts(tool);
    const description = tool.description;
    // Checked whatever the method, so a malformed spec fails in every arm.
    const anchored = parts.filter((part) => part.after).map((part) => {
        const at = description.indexOf(part.after);
        if (at < 0 || description.indexOf(part.after, at + 1) >= 0) {
            invalid(`${tool.name}: "${part.after}" must appear exactly once in its description`);
        }
        return { at: at + part.after.length, text: part.text };
    });
    if (method === 'off' || !parts.length) return description;
    // Inserted from the end, so each place is read in the description as written.
    let text = description;
    [...anchored].sort((a, b) => b.at - a.at).forEach(({ at, text: inserted }) => {
        text = `${text.slice(0, at)}${inserted}${text.slice(at)}`;
    });
    return [text, ...parts.filter((part) => !part.after).map((part) => part.text)].join(' ');
};

/**
 * A check whose text holds method advice. `parts` in order: a string (the
 * fact, the verdict, the refusal) or `{ method: '…' }` (advice, with its own
 * leading space or punctuation). `text` is every part, as the dock and arm B
 * read it; `method` lists the advice, which method 'off' leaves out
 * (checkText). `level` is the check's.
 */
export const methodCheck = (level, ...parts) => {
    const pieces = parts.map((part) => {
        if (typeof part === 'string') return { text: part };
        if (part && typeof part.method === 'string' && part.method) return { text: part.method, advice: true };
        return invalid('a check part is a string or { method: \'…\' }');
    });
    return {
        level,
        text: pieces.map((piece) => piece.text).join(''),
        method: pieces.filter((piece) => piece.advice).map((piece) => piece.text)
    };
};

/**
 * A check's text as the model reads it under `method`: its text, without the
 * advice its `method` lists when the method is 'off'. A TypeError when that
 * advice is not in the text, in order (a check not built by methodCheck):
 * advice is never sent in arm C by mistake, and the error (which the model
 * reads as the call's failure) does not quote it.
 */
export const checkText = (check, method) => {
    const advice = Array.isArray(check?.method) ? check.method : [];
    if (method !== 'off' || !advice.length) return check.text;
    let text = check.text;
    let from = 0;
    for (const piece of advice) {
        const at = typeof piece === 'string' && piece ? text.indexOf(piece, from) : -1;
        if (at < 0) invalid('a check\'s method advice is not in its text, in order (build it with methodCheck)');
        text = `${text.slice(0, at)}${text.slice(at + piece.length)}`;
        from = at;
    }
    return text;
};

// A result as the model reads it under method 'off': each check's text
// without its advice; data, facts and anything else unchanged.
const resultWithoutAdvice = (result) => (result && Array.isArray(result.checks)
    ? { ...result, checks: result.checks.map((check) => (check && typeof check === 'object' ? { level: check.level, text: checkText(check, 'off') } : check)) }
    : result);

/**
 * `tool` as the turn offers it under `method`: its description as the model
 * reads it (toolDescription), no `method` field, and with method 'off' a
 * run() whose checks come without their advice. With 'full' run() is the
 * tool's own: the loop keeps a check's level and text only, as written.
 */
export const toolAsSent = (tool, method) => {
    // eslint's ignoreRestSiblings: `advice` is named only to leave it out.
    const { method: advice, ...rest } = tool;
    const sent = { ...rest, description: toolDescription(tool, method) };
    if (method !== 'off') return sent;
    return { ...sent, run: async (args, context) => resultWithoutAdvice(await tool.run(args, context)) };
};
