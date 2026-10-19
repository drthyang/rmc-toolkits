// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The skill lint (docs/AGENT_PLAN.md, item 5.1): pure functions over the text
// of the Copilot's skills (.claude/skills/<name>/SKILL.md and the reference
// files beside them). It reads the skills and never changes them: the skill
// text is the user's. Three checks, each a list of findings:
//
//   tool       every tool a skill names is one the Copilot can offer. Every
//              snake_case identifier in a skill (backticked, or plain in the
//              prose) must be a tool, or else a parameter of one or a field of
//              what the model sees (a tool's result, the run context); one
//              that is neither names a tool that does not exist.
//   gate       every gate number a skill states equals the constant of the
//              module that decides it (BINDINGS: a phrase that finds the
//              sentence, the number as written, the constants it must equal,
//              optionally the statistic the constant gates). Every number a
//              skill puts after a comparator (≥, below, fewer than …) or
//              before a decimal × must be bound, or listed as no gate with why
//              (NOT_GATES), so a new threshold in a skill cannot go unbound.
//   firm-rule  every "*Firm rule in the app:*" names a check or refusal that
//              exists: each such sentence is mapped (FIRM_RULES) to the code
//              that enforces it, and that code is exercised (a probe), not
//              grepped.
//   value      every value a skill quotes as measured on the Demo that a
//              table binds (VALUES: a phrase, the number as written, and what
//              the tool measures there now) equals the tool's, rounded as the
//              skill writes it (candidate #6: the shell readings).
//   rule       every rule a skill states in words that the app and arm A+
//              apply (RULES, architecture candidate 8: the split site, the
//              cluster elements, the lower group's tolerance) is stated as
//              they apply it: a phrase finds the sentence, a pattern says
//              whether it states the app's rule.
//
// The tables live in tables.js; skillLint.test.js runs the lint on the skills
// the app bundles, against an allow-list of today's violations.

// --- text -------------------------------------------------------------------------

/**
 * A text with every whitespace run collapsed to one space (the skills wrap
 * their lines by hand), and the original line of each character, so a phrase
 * matches across a line break and a finding still names its line.
 */
export const flatten = (text) => {
    const source = String(text || '').replace(/\r/g, '');
    let flat = '';
    const lines = [];
    let line = 1;
    let space = false;
    for (let index = 0; index < source.length; index += 1) {
        const char = source[index];
        if (/\s/.test(char)) {
            if (!space && flat.length) {
                flat += ' ';
                lines.push(line);
            }
            space = true;
            if (char === '\n') line += 1;
            continue;
        }
        space = false;
        flat += char;
        lines.push(line);
    }
    return { flat, lineAt: (at) => lines[Math.max(0, Math.min(at, lines.length - 1))] ?? 1 };
};

/** At most `max` characters of `text` (an allow-list fragment's length). */
export const clip = (text, max = 60) => (text.length > max ? `${text.slice(0, max - 1)}…` : text);

// --- check 1: the tools a skill names ----------------------------------------------

// A snake_case identifier: lower-case first, at least one underscore; a dotted
// path (first_shell.window_A) or a glob (largest_*_refs) inside backticks.
const CODE_SPAN = /`([^`]+)`/g;
const SPAN_PIECE = /[A-Za-z0-9_.*]+/g;
const PROSE_IDENTIFIER = /(?<![\w./*-])([a-z][a-z0-9]*(?:_[A-Za-z0-9]+)+)(?![\w/*-])/g;

const identifierShaped = (token) => /^[a-z]/.test(token) && token.split('.').some((segment) => /[_*]/.test(segment));

/**
 * The identifiers a skill text names: `{ token, at, coded }` (at = index in
 * the flat text; coded = inside backticks). A code span that holds a path
 * (`docs/…`, `data/5K_try1`) or starts with a glob (`*_FT_*`, a file name) is
 * not an identifier. In the prose, a single-letter first part (q_i, U_iso) is
 * mathematics, not code.
 */
export const identifierMentions = (flat) => {
    const found = [];
    const masked = flat.replace(CODE_SPAN, (span, inner, offset) => {
        if (!inner.includes('/') && !inner.trim().startsWith('*')) {
            for (const piece of inner.matchAll(SPAN_PIECE)) {
                const token = piece[0].replace(/\.+$/, '');
                if (identifierShaped(token)) found.push({ token, at: offset + 1 + piece.index, coded: true });
            }
        }
        return ' '.repeat(span.length);
    });
    for (const match of masked.matchAll(PROSE_IDENTIFIER)) {
        if (/^[a-z]_/.test(match[1])) continue;
        found.push({ token: match[1], at: match.index, coded: false });
    }
    return found.sort((a, b) => a.at - b.at);
};

const globPattern = (segment) => new RegExp(`^${segment.split('*').map((part) => part.replace(/[.+?^${}()|[\]\\]/g, '\\$&')).join('[A-Za-z0-9_]*')}$`);

/**
 * An identifier as a tool (its first part is a tool's name), a field (every
 * part is a tool parameter or a result or context field) or unknown.
 */
export const classifyIdentifier = (token, { tools, fields }) => {
    const segments = token.split('.');
    if (tools.has(segments[0])) return { kind: 'tool', tool: segments[0] };
    const known = (segment) => (segment.includes('*')
        ? [...fields].some((field) => globPattern(segment).test(field))
        : fields.has(segment) || tools.has(segment));
    return segments.every(known) ? { kind: 'field' } : { kind: 'unknown' };
};

// --- check 2: gate numbers ------------------------------------------------------------

const NUMBER = '[−–-]?\\d+(?:\\.\\d+)?';
// A number after a comparator (not a <100> direction family), or a decimal
// factor before ×.
const COMPARATOR_NUMBER = new RegExp(
    `(?:≥|≤|(?<![<\\d])>|<(?!\\d+>)|\\b(?:under|below|above|at least|at most|fewer than|more than|less than|within))`
        + `\\s?(?:~|≈|about\\s)?\\s?±?\\s?(${NUMBER})`,
    'g'
);
const FACTOR_NUMBER = /(?<![\d.])(\d+\.\d+)\s?×/g;

/** The numbers a skill states as thresholds: `{ text, at, end }` (the number itself). */
export const gateNumberCandidates = (flat) => {
    const found = [];
    for (const pattern of [COMPARATOR_NUMBER, FACTOR_NUMBER]) {
        for (const match of flat.matchAll(pattern)) {
            const at = match.index + match[0].lastIndexOf(match[1]);
            found.push({ text: match[0], at, end: at + match[1].length });
        }
    }
    // One number once (1.2× is a comparator's number and a factor).
    return found.sort((a, b) => a.at - b.at).filter((item, index, all) => index === 0 || item.at !== all[index - 1].at);
};

const WORD_VALUES = { half: 0.5, most: 0.5 };

/** A number as written in a skill: −0.5, 60 (with unit '%': 0.6), "half", a list "4, 6, 8, 12". */
export const parseStated = (text, { unit = null, list = false } = {}) => {
    const one = (piece) => {
        const word = piece.trim().toLowerCase();
        const value = Object.hasOwn(WORD_VALUES, word) ? WORD_VALUES[word] : Number(word.replace(/[−–]/g, '-'));
        return unit === '%' ? value / 100 : value;
    };
    return list ? String(text).split(/\s*,\s*/).map(one) : one(String(text));
};

const sameValue = (stated, constant) => (Array.isArray(constant)
    ? Array.isArray(stated) && stated.length === constant.length && stated.every((value, index) => Math.abs(value - constant[index]) <= 1e-12)
    : Number.isFinite(stated) && Math.abs(stated - constant) <= 1e-12 * Math.max(1, Math.abs(constant)));

const constantNote = (constants) => constants.map(({ name, value, module }) => `${name} = ${JSON.stringify(value)} (${module})`).join(', ');

/**
 * One binding against one flat text: `{ problem, at, end, stated, written,
 * context }` (`written`: the number as the skill writes it, "60 %", "most");
 * problem null when the sentence states the constant's number (and names the
 * statistic it gates, when the binding says which).
 */
export const evaluateBinding = (binding, flat) => {
    const at = flat.indexOf(binding.find);
    if (at < 0) return { problem: 'stale', at: 0, end: 0, context: binding.find };
    const start = Math.max(0, at - (binding.before ?? 0));
    const window = flat.slice(start, at + binding.find.length + (binding.after ?? 120));
    const match = binding.number.exec(window);
    const context = window.trim();
    if (!match) return { problem: 'no-number', at, end: at + binding.find.length, context };
    const numberAt = start + match.index + match[0].lastIndexOf(match.groups.n);
    const stated = parseStated(match.groups.n, { unit: binding.unit, list: Array.isArray(binding.constants[0].value) });
    const end = numberAt + match.groups.n.length;
    const sentence = flat.slice(Math.min(start + match.index, at), start + match.index + match[0].length);
    const problems = [
        binding.constants.every(({ value }) => sameValue(stated, value)) ? null : 'differs',
        binding.statistic && !binding.statistic.pattern.test(sentence) ? 'statistic' : null
    ].filter(Boolean);
    const written = `${match.groups.n}${binding.unit === '%' ? ' %' : ''}`;
    return { problem: problems[0] ?? null, problems, at: numberAt, end, stated, written, context: sentence };
};

/** The constants a binding names, as an allow-list entry writes the app's side: "CN_CHANGE_REL = 0.15". */
export const decidedBy = (constants) => constants.map(({ name, value }) => `${name} = ${JSON.stringify(value)}`).join(', ');

const bindingDetail = (binding, result) => {
    const constants = constantNote(binding.constants);
    switch (result.problem) {
        case 'stale': return `the binding "${binding.gate}" no longer finds "${clip(binding.find, 50)}": rebind it to the new text, or remove it`;
        case 'no-number': return `states no number for ${binding.gate}, which the app decides with ${constants}`;
        case 'differs': return `states ${JSON.stringify(result.stated)} for ${binding.gate}; the app decides with ${constants}`
            + (result.problems.includes('statistic') ? `, applied to ${binding.statistic.says}` : '');
        case 'statistic': return `${binding.gate}: the app applies ${constants} to ${binding.statistic.says}, which the sentence does not name`;
        default: return '';
    }
};

// --- check 4: measured values ------------------------------------------------------------

const decimalsOf = (text) => (String(text).split('.')[1] || '').length;

/** A measured value as a skill would write it: 3 significant digits (as the tools show numbers). */
export const showValue = (value) => (Number.isFinite(value) ? String(Number(value.toPrecision(3))) : String(value));
const showValues = (values, unit) => values.map((value) => (unit === '%' && Number.isFinite(value) ? `${showValue(value * 100)} %` : showValue(value))).join('–');

/**
 * One measured value against one flat text: `{ problem, at, end, stated,
 * measured, context }`. `entry.number` reads the value as written (group n,
 * and m for a pair such as a window), `entry.unit: '%'` reads 47 as 0.47;
 * `measured` is what the tool shows on the Demo (a number, or [lo, hi]). They
 * agree when the measurement, rounded to the decimals the skill writes,
 * equals the stated number: "3.9" against 3.25 differs, "3" against 2.75 does
 * not (a count written as a whole number, within half an atom). `stated` and
 * `measured` come back as text ("2.33–2.71", "2.21–2.65").
 */
export const evaluateValue = (entry, flat, measured) => {
    const at = flat.indexOf(entry.find);
    if (at < 0) return { problem: 'stale', at: 0, end: 0, context: entry.find };
    const start = Math.max(0, at - (entry.before ?? 0));
    const window = flat.slice(start, at + entry.find.length + (entry.after ?? 120));
    // With group indices ('d'), so a finding names the number's own line.
    const match = new RegExp(entry.number.source, `${entry.number.flags.replace('d', '')}d`).exec(window);
    if (!match) return { problem: 'no-number', at, end: at + entry.find.length, context: window.trim() };
    const keys = ['n', 'm'].filter((key) => match.groups?.[key] != null);
    const texts = keys.map((key) => match.groups[key]);
    const values = Array.isArray(measured) ? measured : [measured];
    const scale = entry.unit === '%' ? 100 : 1;
    const agree = values.length === texts.length && values.every((value, index) => (
        Number.isFinite(value) && Number((value * scale).toFixed(decimalsOf(texts[index]))) === Number(texts[index].replace(/[−–]/g, '-'))
    ));
    const numberAt = start + match.indices.groups[keys[0]][0];
    const end = start + match.indices.groups[keys.at(-1)][1];
    return {
        problem: agree ? null : 'differs',
        at: numberAt,
        end,
        stated: texts.join('–') + (entry.unit === '%' ? ' %' : ''),
        measured: showValues(values, entry.unit),
        context: flat.slice(Math.min(start + match.index, at), start + match.index + match[0].length)
    };
};

// --- check 5: rules the app and arm A+ apply -----------------------------------------------

/**
 * One rule entry against one flat text: `{ problem, at, end, context }`.
 * `entry.find` locates the rule's statement (the finding names its line),
 * `before` / `after` widen the sentence read around it, and `entry.states`
 * is a pattern the sentence matches when it states the rule as the app
 * applies it. problem: null (stated as applied), 'differs', or 'stale' (the
 * phrase is gone).
 */
export const evaluateRule = (entry, flat) => {
    const at = flat.indexOf(entry.find);
    if (at < 0) return { problem: 'stale', at: 0, end: 0, context: entry.find };
    const start = Math.max(0, at - (entry.before ?? 0));
    const end = Math.min(flat.length, at + entry.find.length + (entry.after ?? 0));
    const context = flat.slice(start, end).trim();
    return { problem: entry.states.test(context) ? null : 'differs', at, end, context };
};

// --- check 3: firm rules ----------------------------------------------------------------

// "*Firm rule in the app:*", also "*Firm rule in the app (THY-18):*"; the
// sentence runs to the end of its list item or paragraph.
const FIRM_RULE = /\*Firm rule in the app(?: \([^)]*\))?:\*/g;
const ITEM_END = /\s(?:- |#{1,6} )/g;

/** The firm-rule sentences of a flat text: `{ at, text }` (text after the marker). */
export const firmRuleSentences = (flat) => [...flat.matchAll(FIRM_RULE)].map((match) => {
    const from = match.index + match[0].length;
    ITEM_END.lastIndex = from;
    const end = ITEM_END.exec(flat)?.index ?? flat.length;
    return { at: match.index, marker: match[0], text: flat.slice(from, end).trim() };
});

// --- the lint -------------------------------------------------------------------------

/**
 * Every finding on the skills: `[{ check, kind, file, line, fragment, context,
 * detail }]` (a value finding also carries `stated` and `measured`, a gate
 * finding whose number disagrees `stated` and `decided`). `files`
 * maps a repository path (.claude/skills/<name>/<file>.md) to its text;
 * `tools` is the set of tool names the Copilot can offer, `fields` the
 * parameter and field names the model sees; `bindings`, `notGates`,
 * `firmRules`, `values` and `rules` are the tables (tables.js). Firm-rule
 * probes and value measurements run here, each given `probeContext`: one that
 * throws is a finding. A rule finding also carries the skill's rule (`skill`)
 * and the app's (`app`).
 */
export const lintSkills = async ({
    files, tools, fields, bindings = [], notGates = [], firmRules = [], values = [], rules = [], probeContext = {}
}) => {
    const findings = [];
    const texts = Object.entries(files)
        .sort(([a], [b]) => (a < b ? -1 : a > b ? 1 : 0))
        .map(([path, text]) => [path.replace(/^\.claude\/skills\//, ''), flatten(text)]);
    const named = new Map();

    for (const [file, { flat, lineAt }] of texts) {
        const add = (check, kind, at, context, detail, extra = {}) => findings.push({
            check, kind, file, line: lineAt(at), fragment: clip(context), context, detail, ...extra
        });

        // 1. Tools.
        for (const mention of identifierMentions(flat)) {
            const kind = classifyIdentifier(mention.token, { tools, fields });
            if (kind.kind === 'tool') {
                if (!named.has(file)) named.set(file, new Set());
                named.get(file).add(kind.tool);
            } else if (kind.kind === 'unknown') {
                add('tool', 'unknown', mention.at, mention.token,
                    `\`${mention.token}\` is no tool the Copilot offers, nor a parameter or a result or context field of one`);
            }
        }

        // 2. Gate numbers: the bindings, then every comparator number covered.
        const covered = [];
        for (const binding of bindings.filter((item) => item.file === file)) {
            const result = evaluateBinding(binding, flat);
            if (result.problem !== 'stale' && result.problem !== 'no-number') covered.push([result.at, result.end]);
            // A number that disagrees also says both sides: the skill's as written, the app's constants.
            const sides = result.problem === 'differs' || result.problem === 'statistic'
                ? { stated: result.written, decided: decidedBy(binding.constants) }
                : {};
            if (result.problem) add('gate', result.problem, result.at, result.context, bindingDetail(binding, result), sides);
        }
        for (const item of notGates.filter((entry) => entry.file === file)) {
            const at = flat.indexOf(item.find);
            if (at < 0) add('gate', 'stale', 0, item.find, `the no-gate entry no longer finds "${clip(item.find, 50)}": update or remove it`);
            else covered.push([at, at + item.find.length]);
        }
        for (const candidate of gateNumberCandidates(flat)) {
            if (covered.some(([from, to]) => candidate.at >= from && candidate.end <= to)) continue;
            add('gate', 'unbound', candidate.at, flat.slice(Math.max(0, candidate.at - 30), candidate.end + 20).trim(),
                `"${candidate.text}" is stated as a threshold but bound to no constant: bind it (BINDINGS) or say why it is no gate (NOT_GATES)`);
        }

        // 3. Firm rules.
        const sentences = firmRuleSentences(flat);
        const entries = firmRules.filter((rule) => rule.file === file);
        for (const sentence of sentences) {
            if (!entries.some((rule) => sentence.text.includes(rule.find))) {
                add('firm-rule', 'unmapped', sentence.at, sentence.text,
                    `${sentence.marker} names no check the lint maps: map it (FIRM_RULES) to the answer check, gate rule, tool refusal or question that enforces it`);
            }
        }
        for (const rule of entries) {
            const sentence = sentences.find((item) => item.text.includes(rule.find));
            if (!sentence) {
                add('firm-rule', 'stale', 0, rule.find, `the firm-rule entry no longer finds "${clip(rule.find, 50)}" after a firm-rule marker: update or remove it`);
                continue;
            }
            try {
                await rule.probe(probeContext);
            } catch (failure) {
                add('firm-rule', 'broken', sentence.at, sentence.text,
                    `${rule.enforcedBy} does not enforce it any more: ${String(failure?.message || failure).split('\n')[0]}`);
            }
        }

        // 4. Measured values.
        for (const entry of values.filter((item) => item.file === file)) {
            let measured;
            try {
                measured = await entry.measure(probeContext);
            } catch (failure) {
                add('value', 'broken', Math.max(0, flat.indexOf(entry.find)), entry.find,
                    `${entry.tool} could not be measured: ${String(failure?.message || failure).split('\n')[0]}`);
                continue;
            }
            const result = evaluateValue(entry, flat, measured);
            if (!result.problem) continue;
            const detail = result.problem === 'differs'
                ? `states ${result.stated} for ${entry.what}; ${entry.tool} gives ${result.measured} on the Demo`
                : result.problem === 'stale'
                    ? `the value entry "${entry.what}" no longer finds "${clip(entry.find, 50)}": rebind it to the new text, or remove it`
                    : `states no number for ${entry.what}, which ${entry.tool} measures`;
            add('value', result.problem, result.at, result.context, detail,
                result.problem === 'differs' ? { stated: result.stated, measured: result.measured } : {});
        }

        // 5. Rules the app and arm A+ apply.
        for (const entry of rules.filter((item) => item.file === file)) {
            const result = evaluateRule(entry, flat);
            if (!result.problem) continue;
            const detail = result.problem === 'stale'
                ? `the rule entry "${entry.rule}" no longer finds "${clip(entry.find, 50)}": rebind it to the new text, or remove it`
                : `states ${entry.skill} for ${entry.rule}; the app and arm A+ apply ${entry.app}`;
            add('rule', result.problem, result.at, result.context, detail,
                result.problem === 'differs' ? { skill: entry.skill, app: entry.app } : {});
        }
    }
    return { findings, named };
};

/**
 * The findings an allow-list does not hold (`unexpected`) and the entries no
 * finding needs (`stale`). An entry `{ check, file, fragment, why }` holds a
 * finding of that check in that file whose context contains the fragment.
 */
export const compareAllowList = (findings, allowList) => {
    const holds = (entry, finding) => entry.check === finding.check && entry.file === finding.file && finding.context.includes(entry.fragment);
    return {
        unexpected: findings.filter((finding) => !allowList.some((entry) => holds(entry, finding))),
        stale: allowList.filter((entry) => !findings.some((finding) => holds(entry, finding)))
    };
};

/** A finding as one line of a test message. */
export const describeFinding = (finding) => `${finding.file}:${finding.line} [${finding.check}/${finding.kind}] "${finding.fragment}" — ${finding.detail}`;
