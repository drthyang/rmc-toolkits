// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Checks on the answer itself. A tool's checks judge its result; these judge
// what the model then wrote, against facts the tools and the run context
// established — for the failures the reviews found and that must not recur
// unnoticed (the skills' firm rules):
//
//   - convergence: "converged" while the watchdog or assess_convergence says
//     improving or diverging (rmc-fit-assessment);
//   - chirality, centrosymmetry, a lost inversion centre: claims a group's
//     point group contradicts (local-symmetry; the spot check called P-42₁m
//     chiral and said F-43m lost its centrosymmetry);
//   - a preferred displacement direction when every directions analysis of
//     the question found none at ≥ 3σ (displacement-anisotropy).
//
// Deliberately narrow: a rule fires only on a clear positive claim, not on a
// negated one ("has not converged", "achiral", "no preferred direction"), not
// in a question, and not when the same sentence already says the right thing.
// A false alarm under a correct answer would teach the user to ignore these.
// Each finding is `{ rule, text }`; pure, no React.

const NEGATION = /\b(not|no|never|neither|nor|none|without|lacks?|lacking|absence|absent|cannot|non|yet)\b|n't\b/i;

// Sentences, then clauses: a contrast ("…, but …") starts a new claim.
const sentencesOf = (text) => String(text || '')
    .replace(/\r/g, '')
    .split(/(?<=[.!?])\s+|\n+/)
    .map((sentence) => sentence.trim())
    .filter(Boolean);
const clausesOf = (sentence) => sentence
    .split(/;\s+|:\s+|,\s+(?=(?:but|while|whereas|unlike|although|though|however)\b)/i)
    .map((clause) => clause.trim())
    .filter(Boolean);

// The positive matches of `pattern` in a clause: not in a question, not after
// a negation in the 60 characters before it.
const positiveMatches = (clause, pattern) => {
    if (/\?\s*$/.test(clause)) return [];
    const out = [];
    for (const match of clause.matchAll(pattern)) {
        const before = clause.slice(Math.max(0, match.index - 60), match.index + (match[1] || '').length);
        if (!NEGATION.test(before)) out.push(match);
    }
    return out;
};

// --- space-group symbols ----------------------------------------------------------

const SUBSCRIPTS = { '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4', '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9' };

/** A symbol as written anywhere (P-42₁m, P−4 2_1 m, $P\bar{4}2_1m$, P4̅2₁m) → p-421m. */
export const normalizeSymbol = (text) => String(text || '')
    .replace(/\\bar\{?(\d)\}?/g, '-$1')
    .replace(/(\d)̅/g, '-$1')
    .replace(/[₀-₉]/g, (digit) => SUBSCRIPTS[digit])
    .replace(/[−–—]/g, '-')
    .replace(/\(no\.?\s*\d+\)/gi, '')
    .replace(/[\s_${}\\]/g, '')
    .toLowerCase();

// Known groups (symbol → facts) from the tool facts and the context.
const knownGroups = (facts, context) => {
    const groups = new Map();
    const add = (symbol, details) => {
        const key = normalizeSymbol(symbol);
        if (key && details && !groups.has(key)) groups.set(key, { symbol, ...details });
    };
    for (const fact of facts) for (const group of fact.groups || []) add(group.symbol, group);
    for (const group of context?.groups || []) add(group.symbol, group);
    return groups;
};

// Text with symbols in one spelling but words kept apart: overbars and
// subscripts written as digits, every dash a hyphen.
const symbolText = (text) => String(text || '')
    .replace(/\\bar\{?(\d)\}?/g, '-$1')
    .replace(/(\d)̅/g, '-$1')
    .replace(/[₀-₉]/g, (digit) => SUBSCRIPTS[digit])
    .replace(/[−–—]/g, '-')
    .toLowerCase();

// A normalized symbol as a pattern that allows spacing, underscores, braces
// and $ inside it (P−4 2_1 m) but needs a word boundary around it.
const symbolPattern = (key) => {
    const filler = '[\\s_{}$\\\\]*';
    const body = [...key].map((char) => char.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')).join(filler);
    return new RegExp(`(?<![a-z0-9])${body}(?![a-z0-9])`);
};

// The known groups a piece of text names, each matched as a whole symbol.
const groupsIn = (text, groups) => {
    const prepared = symbolText(text);
    return [...groups].filter(([key]) => symbolPattern(key).test(prepared)).map(([, group]) => group);
};

const groupNote = (group) => `${group.symbol} (${group.point_group})`;

// --- the rules ----------------------------------------------------------------------

const CONVERGED = /(^|[^\w-])(converged|reached convergence|at convergence|plateau(?:ed|s)?)\b/gi;
const STILL_MOVING = /\b(improving|diverg\w*|still (?:falling|decreasing|dropping|declining|descending|going down|moving))\b/i;

const convergenceRule = (text, facts, context) => {
    const reported = [...facts].reverse().find((fact) => fact.convergence)?.convergence ?? context?.convergence ?? null;
    if (reported !== 'improving' && reported !== 'diverging') return [];
    for (const sentence of sentencesOf(text)) {
        if (STILL_MOVING.test(sentence)) continue;
        if (clausesOf(sentence).some((clause) => positiveMatches(clause, CONVERGED).length)) {
            return [{
                rule: 'convergence',
                text: `The answer says the run has converged, but the χ² history is still ${reported} `
                    + '(the watchdog and assess_convergence): it has not converged.'
            }];
        }
    }
    return [];
};

const CHIRAL = /(^|[^\w-])chiral\b/gi;
const CENTROSYMMETRIC = /(^|[^\w-])centrosymmetric\b/gi;
const LOST_INVERSION = [
    /\b(los(?:e|es|t|ing)|break(?:s|ing)?|broke|broken|remov(?:e|es|ed|ing)|destroy(?:s|ed|ing)?)\b[^.;]{0,50}\b(inversion|centrosymmetr\w*|cent(?:re|er) of symmetry)/gi,
    /\b(inversion|centrosymmetr\w*|cent(?:re|er) of symmetry)\b[^.;]{0,30}\b(is|was|are|were|gets|got)\s+(lost|broken|removed|destroyed)/gi
];

const symmetryRules = (text, facts, context) => {
    const groups = knownGroups(facts, context);
    if (!groups.size) return [];
    const found = [];
    const seen = new Set();
    const push = (rule, message) => {
        if (seen.has(rule)) return;
        seen.add(rule);
        found.push({ rule, text: message });
    };
    const answerGroups = groupsIn(text, groups);
    for (const sentence of sentencesOf(text)) {
        for (const clause of clausesOf(sentence)) {
            const named = groupsIn(clause, groups).length ? groupsIn(clause, groups) : groupsIn(sentence, groups);
            if (named.length && named.every((group) => !group.sohncke) && positiveMatches(clause, CHIRAL).length) {
                push('chirality', `The answer calls ${named.map(groupNote).join(', ')} chiral, but a group with mirrors, glides or `
                    + 'rotoinversions is achiral: only a Sohncke group (proper rotations only) allows a chiral structure.');
            }
            if (named.length && named.every((group) => !group.centrosymmetric) && positiveMatches(clause, CENTROSYMMETRIC).length) {
                push('centrosymmetry', `The answer calls ${named.map(groupNote).join(', ')} centrosymmetric, but `
                    + `${named.length > 1 ? 'none of these point groups has' : 'its point group has no'} inversion centre.`);
            }
            const lost = LOST_INVERSION.some((pattern) => positiveMatches(clause, pattern).length);
            const among = answerGroups.length ? answerGroups : [...groups.values()].slice(0, 1);
            if (lost && among.length && among.every((group) => !group.centrosymmetric)) {
                push('inversion-lost', `The answer says an inversion centre is lost, but none of ${among.map(groupNote).join(', ')} `
                    + 'has one: there is no inversion to lose.');
            }
        }
    }
    return found;
};

const PREFERRED = /(^|[^\w-])(preferred|preferential)\s+(?:displacement\s+)?(direction|orientation|axis)\b/gi;

const directionRule = (text, facts) => {
    const analyses = facts.map((fact) => fact.direction).filter(Boolean);
    if (!analyses.length || analyses.some((analysis) => analysis.significant)) return [];
    for (const sentence of sentencesOf(text)) {
        if (clausesOf(sentence).some((clause) => positiveMatches(clause, PREFERRED).length)) {
            const detail = analyses.map((analysis) => `${analysis.subject}: peak ${analysis.peak_sigma}σ`
                + `${analysis.map_sigma != null ? `, map ${analysis.map_sigma}σ` : ''}`).join('; ');
            return [{
                rule: 'direction',
                text: `The answer claims a preferred displacement direction, but no directions analysis in it found one at ≥ 3σ (${detail}).`
            }];
        }
    }
    return [];
};

/**
 * The findings on one answer: `[{ rule, text }]`, empty when nothing
 * contradicts. `facts` are the turn's tool facts; `context` the run's own
 * (`{ convergence, groups }`, contextAuditFacts).
 */
export const auditAnswer = (text, { facts = [], context = null } = {}) => {
    if (!text || !String(text).trim()) return [];
    return [
        ...convergenceRule(text, facts, context),
        ...symmetryRules(text, facts, context),
        ...directionRule(text, facts)
    ];
};

/**
 * The verdict after the checks: a finding contradicts the answer, so a
 * self-assessed 'achieved' or 'partial' (or no verdict at all) becomes
 * 'not verified', keeping what the model claimed; 'not achieved' is already
 * the lowest and stays. The model's own verdict never overrides a computed
 * contradiction. `outcome` is parseOutcome's `{ status, note, body }` or null.
 */
export const lowerOutcome = (outcome, findings, text = '') => {
    if (!Array.isArray(findings) || !findings.length) return outcome;
    if (outcome?.status === 'not achieved') return outcome;
    return {
        status: 'not verified',
        claimed: outcome?.status ?? null,
        note: `${findings.length === 1 ? 'a check below contradicts' : `${findings.length} checks below contradict`} the answer`
            + `${outcome?.status ? ` (the model said "${outcome.status}"${outcome.note ? `: ${outcome.note}` : ''})` : ''}.`,
        body: outcome?.body ?? text
    };
};

/** The run context's facts for the audit: the watchdog's status and the ladder's groups. */
export const contextAuditFacts = (context, pointGroupFacts) => {
    if (!context) return null;
    const groups = [];
    const symmetry = context.symmetry;
    if (symmetry?.space_group && symmetry.point_group) {
        const facts = pointGroupFacts(symmetry.point_group);
        if (facts) groups.push({ symbol: String(symmetry.space_group).replace(/\s*\(No\.\s*\d+\)\s*$/, ''), ...facts });
    }
    for (const rung of symmetry?.ladder || []) {
        const facts = rung.pg ? pointGroupFacts(rung.pg) : null;
        if (facts && rung.sg) groups.push({ symbol: rung.sg, ...facts });
    }
    return { convergence: context.convergence?.watchdog_status ?? null, groups };
};
