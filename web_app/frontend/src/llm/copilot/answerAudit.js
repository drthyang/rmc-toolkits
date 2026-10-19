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
//   - a preferred displacement direction when no directions analysis of the
//     question found a peak at ≥ 3σ (displacement-anisotropy; a significant
//     whole-map test alone is anisotropy, not a direction: review L3);
//   - a polar group the point group says is not polar (local-symmetry);
//   - a superstructure (cell doubling, zone-boundary order) when every
//     analyze_modulations of the question found no ordered point (a general
//     wavevector's order allows "superstructure", never "cell doubling");
//   - cluster claims (distorted, breathing, rigid) that every analyze_clusters
//     of the element contradicts, and a coordination number that every
//     whole-shell count of the element contradicts (review L5).
//
// Deliberately narrow: a rule fires only on a clear positive claim, not on a
// negated one ("has not converged", "achiral", "no preferred direction"), not
// on a hypothetical, conditional or a question ("to check whether it has
// converged", "until χ² plateaus", "I tested whether site 23 has a preferred
// direction", "a converged run"), and not when the same sentence already
// says the right thing.
// A negation or such a word counts only in the stretch of the clause just
// before the claim, after its last comma or "and"/"so": "With no remaining
// drift in χ², the run has converged" is a claim.
// A JSON key ("chiral": false in a structured block) is not a claim either.
// A false alarm under a correct answer would teach the user to ignore these.
// Each finding is `{ rule, text }`; pure, no React.

import { cnMatches } from '../context/shells';

const NEGATION = /\b(not|no|never|neither|nor|none|without|lacks?|lacking|absence|absent|cannot|non|yet|stopped|ceased)\b|n't\b/i;
// What makes a mention hypothetical, conditional or a question about it.
const NOT_ASSERTED = new RegExp('\\b(whether|if|unless|until|once|when(?:ever)?|before|would|could|might|'
    + 'to (?:check|test|see|confirm|verify|determine|assess|judge|tell|decide)|'
    + 'test(?:s|ed|ing)? (?:for|whether|if)|check(?:s|ed|ing)? (?:for|whether|if|that)|'
    + 'look(?:s|ed|ing)? for|make sure|ensure|verify that|confirm that)\\b', 'i');
// An attributive use of a state names a kind, it claims nothing: "a
// converged run" (but "has a preferred direction" is a claim).
const ATTRIBUTIVE = /\b(a|an|any|every|each)\s+$/i;
// The stretch of a clause a claim's negation or condition must sit in: after
// its last comma, "and", "so", "therefore", "thus" or "hence".
const leadOf = (text) => text.split(/,\s*|\s+(?:and|so|therefore|thus|hence)\s+/i).pop();

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

// The asserted matches of `pattern` in a clause: not in a question, and not
// after a negation or a hypothetical word in the stretch just before it (nor
// after an article, for a `state` such as "converged").
const positiveMatches = (clause, pattern, { state = false } = {}) => {
    if (/\?\s*$/.test(clause)) return [];
    const out = [];
    for (const match of clause.matchAll(pattern)) {
        // A JSON key ("chiral": false) names a field, it claims nothing.
        if (/^"\s*(?::|$)/.test(clause.slice(match.index + match[0].length))) continue;
        const lead = leadOf(clause.slice(0, match.index + (match[1] || '').length));
        if (!NEGATION.test(lead) && !NOT_ASSERTED.test(lead) && !(state && ATTRIBUTIVE.test(lead))) out.push(match);
    }
    return out;
};

/**
 * Whether a text asserts `pattern` somewhere: a positive match in one of its
 * clauses, read as the checks read claims (not negated, not hypothetical, not
 * a question, not a JSON key). Group 1 of the pattern is the character before
 * the claim. `within(sentence)`, when given, keeps only the sentences it
 * accepts. Shared with the evidence gate (evidence.js).
 */
export const asserts = (text, pattern, { within = null } = {}) => sentencesOf(text)
    .some((sentence) => (!within || within(sentence))
        && clausesOf(sentence).some((clause) => positiveMatches(clause, pattern).length > 0));

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

// Something written like a space-group symbol (a lattice letter, then axes
// with at least one digit: P2₁2₁2, I4/mmm, R-3c), known or not.
const SYMBOL_LIKE = /(?<![a-z0-9])[pifcabr](?=[a-z/-]*\d)(?:-?\d(?:_?\d)?|\/?[mnabcde]){1,5}(?![a-z0-9])/;
const namesSymbol = (text) => SYMBOL_LIKE.test(symbolText(text).replace(/[\s{}$\\]/g, (char) => (/\s/.test(char) ? ' ' : '')));

// The groups a claim is about, before a symmetry lowering: the context's
// average group and the top of its ladder (`parent`), or, for groups given
// without that flag, the context's first.
const parentGroups = (context) => {
    const groups = context?.groups || [];
    const flagged = groups.filter((group) => group.parent);
    return flagged.length ? flagged : groups.slice(0, 1);
};

// --- the rules ----------------------------------------------------------------------

const CONVERGED = /(^|[^\w-])(converged|reached convergence|at convergence|plateau(?:ed|s)?)\b/gi;
const STILL_MOVING = /(^|\W)(improving|diverg\w*|still (?:falling|decreasing|dropping|declining|descending|going down|moving))\b/gi;

const convergenceRule = (text, facts, context) => {
    const reported = [...facts].reverse().find((fact) => fact.convergence)?.convergence ?? context?.convergence ?? null;
    if (reported !== 'improving' && reported !== 'diverging') return [];
    for (const sentence of sentencesOf(text)) {
        // A sentence that itself says the run still moves already says the
        // right thing; "stopped improving" or "no longer diverging" does not.
        if (clausesOf(sentence).some((clause) => positiveMatches(clause, STILL_MOVING).length)) continue;
        if (clausesOf(sentence).some((clause) => positiveMatches(clause, CONVERGED, { state: true }).length)) {
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
// "polar" as a property ("non-polar", "nonpolar" and "polarisation" do not match).
const POLAR = /(^|[^\w-])polar\b(?!\s+(?:plot|coordinates?|angle))/gi;
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
    const parents = parentGroups(context);
    for (const sentence of sentencesOf(text)) {
        for (const clause of clausesOf(sentence)) {
            // A clause that names no group is about the sentence's ("…; it is
            // chiral"), but not one that names a group the checks do not know
            // ("its subgroup P2₁2₁2 is chiral" is not about the parent).
            const inClause = groupsIn(clause, groups);
            const named = inClause.length ? inClause : (namesSymbol(clause) ? [] : groupsIn(sentence, groups));
            if (named.length && named.every((group) => !group.sohncke) && positiveMatches(clause, CHIRAL).length) {
                push('chirality', `The answer calls ${named.map(groupNote).join(', ')} chiral, but a group with mirrors, glides or `
                    + 'rotoinversions is achiral: only a Sohncke group (proper rotations only) allows a chiral structure.');
            }
            if (named.length && named.every((group) => !group.centrosymmetric) && positiveMatches(clause, CENTROSYMMETRIC).length) {
                push('centrosymmetry', `The answer calls ${named.map(groupNote).join(', ')} centrosymmetric, but `
                    + `${named.length > 1 ? 'none of these point groups has' : 'its point group has no'} inversion centre.`);
            }
            if (named.length && named.every((group) => group.polar === false) && positiveMatches(clause, POLAR).length) {
                push('polar', `The answer calls ${named.map(groupNote).join(', ')} polar, but `
                    + `${named.length > 1 ? 'none of these point groups is' : 'its point group is not'} polar: an operation reverses every direction.`);
            }
            // Losing inversion is wrong only when no group it could be lost
            // from has it: neither a group the answer names nor the parent.
            const lost = LOST_INVERSION.some((pattern) => positiveMatches(clause, pattern).length);
            const among = [...new Map([...answerGroups, ...parents].map((group) => [normalizeSymbol(group.symbol), group])).values()];
            if (lost && among.length && among.every((group) => !group.centrosymmetric)) {
                push('inversion-lost', `The answer says an inversion centre is lost, but none of ${among.map(groupNote).join(', ')} `
                    + 'has one: there is no inversion to lose.');
            }
        }
    }
    return found;
};

export const PREFERRED = /(^|[^\w-])(preferred|preferential)\s+(?:displacement\s+)?(direction|orientation|axis)\b/gi;

const directionRule = (text, facts) => {
    const analyses = facts.map((fact) => fact.direction).filter(Boolean);
    if (!analyses.length || analyses.some((analysis) => analysis.significant)) return [];
    for (const sentence of sentencesOf(text)) {
        if (clausesOf(sentence).some((clause) => positiveMatches(clause, PREFERRED).length)) {
            const detail = analyses.map((analysis) => `${analysis.subject}: peak ${analysis.peak_sigma}σ`
                + `${analysis.map_sigma != null ? `, map ${analysis.map_sigma}σ` : ''}`).join('; ');
            const anisotropic = analyses.some((analysis) => analysis.anisotropic);
            return [{
                rule: 'direction',
                text: `The answer claims a preferred displacement direction, but no directions analysis in it found a peak at ≥ 3σ (${detail}).`
                    + (anisotropic ? ' A significant map test is anisotropy, which any thermal ellipsoid shows, not a preferred direction.' : '')
            }];
        }
    }
    return [];
};

// A claimed superstructure: cell doubling or ordered zone-boundary
// displacements (needs a zone-boundary point), or order between cells in
// general (a zone-boundary point or a general wavevector will do).
export const CELL_DOUBLING = /(^|[^\w-])(cell[- ]doubling|doubl(?:es|ed|ing) (?:of )?the (?:unit )?cell|zone[- ]boundary (?:order|modulation|mode|distortion)|commensurate modulation)\b/gi;
export const SUPERSTRUCTURE = /(^|[^\w-])(superstructure|supercell order|ordered (?:\w+ )?displacements? at q|incommensurate (?:order|modulation))\b/gi;

const superstructureRule = (text, facts) => {
    const analyses = facts.map((fact) => fact.modulation).filter(Boolean);
    if (!analyses.length) return [];
    const doubling = analyses.some((analysis) => analysis.ordered?.length);
    const any = doubling || analyses.some((analysis) => analysis.general?.length);
    for (const sentence of sentencesOf(text)) {
        const clauses = clausesOf(sentence);
        if (!doubling && clauses.some((clause) => positiveMatches(clause, CELL_DOUBLING).length)) {
            return [{
                rule: 'superstructure',
                text: 'The answer claims a cell doubling, but no analyze_modulations of this question found an ordered '
                    + 'zone-boundary point (each was within the independent-displacement background after the correction).'
            }];
        }
        if (!any && clauses.some((clause) => positiveMatches(clause, SUPERSTRUCTURE).length)) {
            return [{
                rule: 'superstructure',
                text: 'The answer claims a superstructure, but no analyze_modulations of this question found order between cells '
                    + '(no zone-boundary point or general wavevector beyond the independent-displacement background after the correction).'
            }];
        }
    }
    return [];
};

// --- cluster and coordination claims (the review's L5, approved 2026-10-14) ---------

// An element symbol as written in a formula or alone (Mo, Mo3, Ta₄), not inside a word.
const elementPattern = (element) => new RegExp(`(?<![A-Za-z])${element}(?:[0-9₀-₉]+)?(?![a-z])`);
const namesElement = (text, element) => elementPattern(element).test(text);

const CLUSTER_NOUN = '(?:clusters?|tetrahedr(?:on|a)|trimers?|dimers?|triangles?|octahedr(?:on|a)|polyhedr(?:on|a))';
// An adverb may stand between the verb and "distorted", but not a negation.
export const DISTORTED = new RegExp(`(^|[^\\w-])(distorted\\s+(?:\\w+\\s+)?${CLUSTER_NOUN}|${CLUSTER_NOUN}\\s+(?:are|is|appear|look|seem)\\s+`
    + '(?:(?!not\\b|never\\b|hardly\\b|barely\\b)\\w+\\s+)?(?:distorted|deformed)|cluster\\s+distortions?)\\b', 'gi');
export const BREATHING = /(^|[^\w-])(breath(?:e|es|ing))\b/gi;
export const RIGID = /(^|[^\w-])(rigid(?:ly)?\b|move(?:s|d)?\s+(?:together\s+)?as\s+(?:a\s+)?(?:rigid\s+)?(?:units?|bod(?:y|ies))\b)/gi;
const CLUSTER_CLAIMS = [
    { key: 'distorted', pattern: DISTORTED, holds: (fact) => fact.distorted || fact.split || fact.static, says: 'are distorted', found: 'no shape distortion beyond independent atoms, no split sites and a regular average cluster' },
    { key: 'breathing', pattern: BREATHING, holds: (fact) => fact.breathing, says: 'breathe', found: 'no size fluctuation beyond independent atoms' },
    { key: 'rigid', pattern: RIGID, holds: (fact) => fact.rigid, says: 'are rigid', found: 'no rigidity beyond independent atoms (or split sites, where the null already holds the distortion)' }
];

// A claim about clusters is checked against the analyses of the element it
// names (or, naming none but saying "cluster", all of them); a polyhedron word
// without an element or "cluster" may be about another polyhedron and is left.
const clusterRule = (text, facts) => {
    const analyses = facts.map((fact) => fact.cluster).filter(Boolean);
    if (!analyses.length) return [];
    const found = [];
    for (const claim of CLUSTER_CLAIMS) {
        let hit = null;
        for (const sentence of sentencesOf(text)) {
            for (const clause of clausesOf(sentence)) {
                if (!positiveMatches(clause, claim.pattern).length) continue;
                const named = analyses.filter((fact) => namesElement(clause, fact.element) || namesElement(sentence, fact.element));
                const about = named.length ? named : (/\bclusters?\b/i.test(clause) ? analyses : []);
                if (about.length && about.every((fact) => !claim.holds(fact))) hit = about;
            }
            if (hit) break;
        }
        if (hit) {
            const which = [...new Set(hit.map((fact) => `${fact.element}${fact.members}`))].join(', ');
            found.push({
                rule: 'cluster',
                text: `The answer says the ${which} clusters ${claim.says}, but analyze_clusters found ${claim.found}.`
            });
        }
    }
    return found;
};

const WORD_NUMBER = { two: 2, three: 3, four: 4, five: 5, six: 6, seven: 7, eight: 8, nine: 9, ten: 10, eleven: 11, twelve: 12 };
const COUNT = '(\\d+(?:\\.\\d+)?|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve)';
export const CN_CLAIMS = [
    new RegExp(`(^|[^\\w.-])${COUNT}[- ]?(?:fold[- ]?)?coordinat(?:e|ed|ion)\\b`, 'gi'),
    new RegExp(`(^|[^\\w-])(?:coordination(?:\\s+number)?|CN)\\s*(?:of|is|=|equals?|:)?\\s*${COUNT}(?!\\d|\\.\\d)`, 'gi')
];
export const POLYHEDRON_CN = [
    [/(^|[^\w-])(tetrahedral(?:ly)?)\s+(?:\w+\s+)?coordinat/gi, 4],
    [/(^|[^\w-])(octahedral(?:ly)?)\s+(?:\w+\s+)?coordinat/gi, 6],
    [/(^|[^\w-])(square[- ]planar)\s+(?:\w+\s+)?coordinat/gi, 4],
    [/(^|[^\w-])(cub(?:ic|ically))\s+(?:\w+\s+)?coordinat/gi, 8]
];
const countOf = (word) => (WORD_NUMBER[word.toLowerCase()] ?? Number(word));
// A clause about part of the shell (the first peak, the short bonds, a window
// in Å) is not a claim about the whole shell's count.
const PARTIAL_SHELL = /\bfirst[- ]peak\b|\bshort(?:est|er)?\s+bonds?\b|\bwithin\s+[\d.]+\s*(?:[–—−-]|to)\s*[\d.]+\s*(?:Å|A\b)|\d+(?:\.\d+)?\s*[–—−-]\s*\d+(?:\.\d+)?\s*(?:Å|A\b)/i;

// A coordination number claimed for an element, against every whole-shell
// count of it (find_gr_peaks' first shell, compute_bond_angles in a whole
// shell) — not against a window the user chose, whose count stands. A claim
// matches a count within half an atom (shells.js cnMatches).
const coordinationRule = (text, facts) => {
    const counts = facts.flatMap((fact) => (Array.isArray(fact.coordination) ? fact.coordination : (fact.coordination ? [fact.coordination] : [])))
        .filter((count) => Number.isFinite(count.cn));
    if (!counts.length) return [];
    const centres = [...new Set(counts.map((count) => count.centre))];
    // The centre a claim is about: the named element nearest to it.
    const nearestCentre = (clause, at) => {
        let best = null;
        centres.forEach((centre) => {
            for (const match of clause.matchAll(new RegExp(elementPattern(centre).source, 'g'))) {
                const distance = Math.abs(match.index - at);
                if (!best || distance < best.distance) best = { centre, distance };
            }
        });
        return best?.centre ?? null;
    };
    for (const sentence of sentencesOf(text)) {
        for (const clause of clausesOf(sentence)) {
            if (PARTIAL_SHELL.test(clause)) continue;
            const claimed = [
                ...CN_CLAIMS.flatMap((pattern) => positiveMatches(clause, pattern).map((match) => ({ n: countOf(match[2]), at: match.index }))),
                ...POLYHEDRON_CN.flatMap(([pattern, cn]) => positiveMatches(clause, pattern).map((match) => ({ n: cn, at: match.index })))
            ].filter((claim) => Number.isFinite(claim.n));
            for (const { n, at } of claimed) {
                const centre = nearestCentre(clause, at);
                if (!centre) continue;
                const mine = counts.filter((count) => count.centre === centre);
                const partners = [...new Set(mine.map((count) => count.partner))].filter((partner) => partner !== centre && namesElement(clause, partner));
                const relevant = partners.length === 1 ? mine.filter((count) => count.partner === partners[0]) : mine;
                const whole = relevant.filter((count) => count.whole_shell);
                if (!whole.length) continue;
                if (relevant.some((count) => (count.whole_shell || count.confirmed) && cnMatches(count.cn, n))) continue;
                const shells = whole.map((count) => `${count.cn} ${count.partner} per ${centre}${count.window ? ` in ${count.window[0]}–${count.window[1]} Å` : ''}`);
                return [{
                    rule: 'coordination',
                    text: `The answer gives ${centre} a coordination of ${n}, but its whole first shell holds ${[...new Set(shells)].join('; ')} `
                        + '(a window that cuts the shell counts fewer: say so if that window was meant).'
                }];
            }
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
        ...directionRule(text, facts),
        ...superstructureRule(text, facts),
        ...clusterRule(text, facts),
        ...coordinationRule(text, facts)
    ];
};

/**
 * The run context's facts for the audit: the watchdog's status and the
 * ladder's groups. The average group and the ladder's top (its rung with the
 * most operations) are the `parent` groups a lowering starts from.
 */
export const contextAuditFacts = (context, pointGroupFacts) => {
    if (!context) return null;
    const groups = [];
    const symmetry = context.symmetry;
    if (symmetry?.space_group && symmetry.point_group) {
        const facts = pointGroupFacts(symmetry.point_group);
        if (facts) groups.push({ symbol: String(symmetry.space_group).replace(/\s*\(No\.\s*\d+\)\s*$/, ''), ...facts, parent: true });
    }
    const ladder = symmetry?.ladder || [];
    const most = Math.max(...ladder.map((rung) => (Number.isFinite(rung.n_ops) ? rung.n_ops : -Infinity)));
    ladder.forEach((rung, index) => {
        const facts = rung.pg ? pointGroupFacts(rung.pg) : null;
        // Without operation counts, the ladder's two ends.
        const top = Number.isFinite(most) ? rung.n_ops === most : index === 0 || index === ladder.length - 1;
        if (facts && rung.sg) groups.push({ symbol: rung.sg, ...facts, ...(top ? { parent: true } : {}) });
    });
    return { convergence: context.convergence?.watchdog_status ?? null, groups };
};
