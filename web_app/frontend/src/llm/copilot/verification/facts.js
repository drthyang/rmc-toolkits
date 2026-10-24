// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The facts catalogue (FACTS, the checks on the answer): one row per key of
// the `facts` a tool returns beside its data (never sent to the model), in
// the order the checks run. A tool's checks judge its result; these judge
// what the model then wrote, against what the tools and the run context
// established — the failures the reviews found and that must not recur
// unnoticed (the skills' firm rules):
//
//   convergence   "converged" while the watchdog or assess_convergence says
//                 improving or diverging — their verdict on every χ² term of
//                 the .log and their total (rmc-fit-assessment; the user's
//                 rule, 2026-10-24);
//   groups        chirality, centrosymmetry, polarity and a lost inversion
//                 centre that a group's point group contradicts
//                 (local-symmetry; the spot check called P-42₁m chiral and said
//                 F-43m lost its centrosymmetry);
//   direction     a preferred displacement direction when no directions
//                 analysis of the question found a peak at ≥ 3σ
//                 (displacement-anisotropy; a significant whole-map test alone
//                 is anisotropy, not a direction: review L3);
//   modulation    a superstructure (cell doubling, zone-boundary order) when
//                 every analyze_modulations of the question found no ordered
//                 point (a general wavevector's order allows "superstructure",
//                 never "cell doubling"); and order claimed at a box-scale
//                 wavevector, a long wave (every |q_i| ≤ the engine's
//                 long_wave_q), whatever its power (the user's rule,
//                 2026-10-24: a long wave says nothing about order at its q);
//   cluster       cluster claims (distorted, breathing, rigid) every
//                 analyze_clusters of the element contradicts (review L5);
//   coordination  a coordination number every whole-shell count of the
//                 element contradicts (review L5);
//   distribution  a distribution's numbers (bond angles, coordination counts,
//                 displacement directions, cluster shapes) quoted without
//                 saying they come from one configuration (the user's rule,
//                 2026-10-24: the app always analyses one);
//   lowering      measure_symmetry_lowering's shifts: no check reads them yet
//                 (contradicts: null, so the missing consumer shows).
//
// Each row: `key`; `from`, the tools that return it; `shape`, its shape in
// the dataShape notation (copilot/dataShape.js; checked against the real
// tools on the Demo by src/__tests__/verificationFacts.test.js); `reads`, the
// phrasings its check reads (phrasings.js); and `contradicts(reading, facts,
// context)` → findings, where `reading` is the answer read once ({ text,
// sentences }), `facts` the turn's tools' facts in the order they ran, and
// `context` the run context's own (contextAuditFacts). Each finding is
// `{ rule, text }`. A facts key no row names is kept but read by nothing.
// Pure, no React.

import { oneOf } from '../dataShape';
import { cnMatches } from '../../context/shells';
import { pointGroupFacts } from '../../context/pointGroups';
import { clausesOf, elementPattern, groupsIn, namesElement, namesSymbol, normalizeSymbol, positiveMatches, qVectorsIn } from './reading';
import {
    BREATHING, CELL_DOUBLING, CENTROSYMMETRIC, CHIRAL, CN_CLAIMS, CONVERGED, DISTORTED, DISTRIBUTION_NOUN, DISTRIBUTION_STAT, LONG_WAVE_WORDS,
    LOST_INVERSION, ORDER_CLAIM, PARTIAL_SHELL, POLAR, POLYHEDRON_CN, PREFERRED, RIGID, SAYS_NOTHING, SINGLE_CONFIGURATION, STILL_MOVING,
    SUPERSTRUCTURE, WORD_NUMBER
} from './phrasings';

// --- the run context's own facts ------------------------------------------------------

/**
 * The run context's facts for the checks: the watchdog's status and the
 * ladder's groups. The average group and the ladder's top (its rung with the
 * most operations) are the `parent` groups a lowering starts from. Null
 * without a context.
 */
export const contextAuditFacts = (context) => {
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

// --- convergence ----------------------------------------------------------------------

const convergenceRule = (reading, facts, context) => {
    const reported = [...facts].reverse().find((fact) => fact.convergence)?.convergence ?? context?.convergence ?? null;
    if (reported !== 'improving' && reported !== 'diverging') return [];
    for (const sentence of reading.sentences) {
        // A sentence that itself says the run still moves already says the
        // right thing; "stopped improving" or "no longer diverging" does not.
        if (clausesOf(sentence).some((clause) => positiveMatches(clause, STILL_MOVING).length)) continue;
        if (clausesOf(sentence).some((clause) => positiveMatches(clause, CONVERGED, { state: true }).length)) {
            return [{
                rule: 'convergence',
                text: `The answer says the run has converged, but the χ² history is still ${reported} `
                    + '(the watchdog and assess_convergence, judging every χ² term and their total): it has not converged.'
            }];
        }
    }
    return [];
};

// --- groups ---------------------------------------------------------------------------

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

const groupNote = (group) => `${group.symbol} (${group.point_group})`;

// The groups a claim is about, before a symmetry lowering: the context's
// average group and the top of its ladder (`parent`), or, for groups given
// without that flag, the context's first.
const parentGroups = (context) => {
    const groups = context?.groups || [];
    const flagged = groups.filter((group) => group.parent);
    return flagged.length ? flagged : groups.slice(0, 1);
};

const symmetryRules = (reading, facts, context) => {
    const groups = knownGroups(facts, context);
    if (!groups.size) return [];
    const found = [];
    const seen = new Set();
    const push = (rule, message) => {
        if (seen.has(rule)) return;
        seen.add(rule);
        found.push({ rule, text: message });
    };
    const answerGroups = groupsIn(reading.text, groups);
    const parents = parentGroups(context);
    for (const sentence of reading.sentences) {
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

// --- direction ------------------------------------------------------------------------

const directionRule = (reading, facts) => {
    const analyses = facts.map((fact) => fact.direction).filter(Boolean);
    if (!analyses.length || analyses.some((analysis) => analysis.significant)) return [];
    for (const sentence of reading.sentences) {
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

// --- modulation -----------------------------------------------------------------------

const superstructureRule = (reading, facts) => {
    const analyses = facts.map((fact) => fact.modulation).filter(Boolean);
    if (!analyses.length) return [];
    const doubling = analyses.some((analysis) => analysis.ordered?.length);
    const any = doubling || analyses.some((analysis) => analysis.general?.length);
    for (const sentence of reading.sentences) {
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

// --- long waves (the user's rule, 2026-10-24) -------------------------------------------

// "A long wave is a correlation the finite box can't represent, so it says
// nothing about order at its own q; a point or mode that passes its own test
// at another q is judged on that test alone" (the user, 2026-10-24). A clause
// that claims order (ORDER_CLAIM, asserted) at a box-scale wavevector — a q
// it writes whose every |q_i| is at most the engine's bound (facts'
// long_wave_q, the modulation engine's isLongWave), or a long or box-scale
// wave in words — is flagged, whatever the mode's power; a clause that
// denies order anywhere ("a long wave, not order", "says nothing about
// order") claims none. Only once an analyze_modulations ran, as the
// superstructure rule.
const isLongWaveQ = (q, bound) => q.some((value) => value !== 0) && q.every((value) => Math.abs(value) <= bound + 1e-9);
const qText = (q) => `(${q.map((value) => (Math.abs(value - 0.5) < 1e-9 ? '½' : String(Number(value.toFixed(3))))).join(', ')})`;

const longWaveRule = (reading, facts) => {
    const analyses = facts.map((fact) => fact.modulation).filter((fact) => fact && Number.isFinite(fact.long_wave_q));
    if (!analyses.length) return [];
    const bound = Math.max(...analyses.map((fact) => fact.long_wave_q));
    for (const sentence of reading.sentences) {
        for (const clause of clausesOf(sentence)) {
            const claims = positiveMatches(clause, ORDER_CLAIM).filter((match) => !SAYS_NOTHING.test(clause.slice(0, match.index)));
            if (!claims.length) continue;
            // A clause that also denies order is about the long wave, not a claim of it.
            if ([...clause.matchAll(ORDER_CLAIM)].length > claims.length) continue;
            const wave = qVectorsIn(clause).find((q) => isLongWaveQ(q, bound));
            if (!wave && !LONG_WAVE_WORDS.test(clause)) continue;
            const at = wave ? `order at q = ${qText(wave)}, a box-scale wavevector (every |q_i| ≤ ${bound}: a long wave)` : 'order at a long (box-scale) wave';
            return [{
                rule: 'long-wave',
                text: `The answer claims ${at}, but a long wave is a correlation the finite box can't represent, so it says nothing `
                    + 'about order at its own q; a point or mode that passes its own test at another q is judged on that test alone.'
            }];
        }
    }
    return [];
};

const modulationRules = (reading, facts, context) => [...superstructureRule(reading, facts, context), ...longWaveRule(reading, facts)];

// --- distributions from one configuration (the user's rule, 2026-10-24) -----------------

// Every analysis in the app reads one configuration (the run's .rmc6f): a
// distribution it gives — bond angles and the coordination counts, the
// displacement directions, the clusters' shapes (facts.distribution) — is one
// configuration's, with no spread between runs to judge it by. An answer that
// quotes such a distribution's numbers (a sentence with a statistic of it,
// DISTRIBUTION_STAT, that names a distribution or gives degrees,
// DISTRIBUTION_NOUN) and nowhere says they come from one configuration
// (SINGLE_CONFIGURATION) is flagged. Only once such an analysis ran.
const singleConfigurationRule = (reading, facts) => {
    const distributions = facts.map((fact) => fact.distribution).filter(Boolean);
    if (!distributions.length || SINGLE_CONFIGURATION.test(reading.text)) return [];
    const quotes = reading.sentences.some((sentence) => DISTRIBUTION_STAT.test(sentence) && DISTRIBUTION_NOUN.test(sentence));
    if (!quotes) return [];
    const subjects = [...new Set(distributions.map((fact) => fact.subject))];
    return [{
        rule: 'single-configuration',
        text: `The answer quotes the numbers of ${subjects.join(', ')} without saying they come from one configuration: every `
            + 'analysis here reads the run\'s one configuration, so its distribution carries no spread between runs; say so.'
    }];
};

// --- cluster (the review's L5, approved 2026-10-14) -------------------------------------

const CLUSTER_CLAIMS = [
    { key: 'distorted', pattern: DISTORTED, holds: (fact) => fact.distorted || fact.split || fact.static, says: 'are distorted', found: 'no shape distortion beyond independent atoms, no split sites and a regular average cluster' },
    { key: 'breathing', pattern: BREATHING, holds: (fact) => fact.breathing, says: 'breathe', found: 'no size fluctuation beyond independent atoms' },
    { key: 'rigid', pattern: RIGID, holds: (fact) => fact.rigid, says: 'are rigid', found: 'no rigidity beyond independent atoms (or split sites, where the null already holds the distortion)' }
];

// A claim about clusters is checked against the analyses of the element it
// names (or, naming none but saying "cluster", all of them); a polyhedron word
// without an element or "cluster" may be about another polyhedron and is left.
const clusterRule = (reading, facts) => {
    const analyses = facts.map((fact) => fact.cluster).filter(Boolean);
    if (!analyses.length) return [];
    const found = [];
    for (const claim of CLUSTER_CLAIMS) {
        let hit = null;
        for (const sentence of reading.sentences) {
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

// --- coordination (the review's L5) -----------------------------------------------------

const countOf = (word) => (WORD_NUMBER[word.toLowerCase()] ?? Number(word));

// A coordination number claimed for an element, against every whole-shell
// count of it (find_gr_peaks' first shell, compute_bond_angles in a whole
// shell) — not against a window the user chose, whose count stands. A claim
// matches a count within half an atom (shells.js cnMatches).
const coordinationRule = (reading, facts) => {
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
    for (const sentence of reading.sentences) {
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

// --- the catalogue --------------------------------------------------------------------

const COORDINATION_COUNT = {
    centre: 'string', partner: 'string', window: oneOf(['number', 'number'], 'null'), cn: 'number|null', whole_shell: 'boolean', confirmed: 'boolean'
};

export const FACTS = Object.freeze([
    Object.freeze({
        key: 'convergence',
        from: Object.freeze(['assess_convergence']),
        // The watchdog's status (watchdog/heuristics.js classifyConvergence);
        // the run context's watchdog_status is read too.
        shape: 'string',
        reads: Object.freeze([CONVERGED, STILL_MOVING]),
        contradicts: convergenceRule
    }),
    Object.freeze({
        key: 'groups',
        from: Object.freeze(['detect_symmetry']),
        // The group and its point group's facts (context/pointGroups.js); the
        // run context's ladder adds its groups, the parents flagged.
        shape: [{ symbol: 'string', point_group: 'string', centrosymmetric: 'boolean', sohncke: 'boolean', polar: 'boolean', crystal_system: 'string' }],
        reads: Object.freeze([CHIRAL, CENTROSYMMETRIC, POLAR, ...LOST_INVERSION]),
        contradicts: symmetryRules
    }),
    Object.freeze({
        key: 'direction',
        from: Object.freeze(['analyze_displacement_directions']),
        shape: { subject: 'string', significant: 'boolean', anisotropic: 'boolean', peak_sigma: 'number', map_sigma: 'number|null' },
        reads: Object.freeze([PREFERRED]),
        contradicts: directionRule
    }),
    Object.freeze({
        key: 'modulation',
        from: Object.freeze(['analyze_modulations']),
        // long_wave_q: the engine's long-wave bound (rules.longWaveQ);
        // long_waves: the long waves among each element's strongest modes.
        shape: {
            elements: ['string'],
            ordered: [{ element: 'string', q: ['number', 'number', 'number'] }],
            general: [{ element: 'string', q: ['number', 'number', 'number'] }],
            long_wave_q: 'number',
            long_waves: [{ element: 'string', q: ['number', 'number', 'number'] }]
        },
        reads: Object.freeze([CELL_DOUBLING, SUPERSTRUCTURE, ORDER_CLAIM, LONG_WAVE_WORDS, SAYS_NOTHING]),
        // The superstructure rule, then the long-wave rule (the user's, 2026-10-24).
        contradicts: modulationRules
    }),
    Object.freeze({
        key: 'cluster',
        from: Object.freeze(['analyze_clusters']),
        shape: {
            element: 'string', members: 'number', window: ['number', 'number'], distorted: 'boolean|null', static: 'boolean|null',
            breathing: 'boolean|null', rigid: 'boolean|null', split: 'boolean|null', families: ['string']
        },
        reads: Object.freeze([DISTORTED, BREATHING, RIGID]),
        contradicts: clusterRule
    }),
    Object.freeze({
        key: 'coordination',
        from: Object.freeze(['find_gr_peaks', 'compute_bond_angles']),
        // find_gr_peaks: the whole first shell, each way; compute_bond_angles:
        // the count in its window, whole shell or the user's.
        shape: oneOf([COORDINATION_COUNT], COORDINATION_COUNT),
        reads: Object.freeze([...CN_CLAIMS, ...POLYHEDRON_CN.map(([pattern]) => pattern), PARTIAL_SHELL]),
        contradicts: coordinationRule
    }),
    Object.freeze({
        key: 'distribution',
        from: Object.freeze(['compute_bond_angles', 'analyze_displacement_directions', 'analyze_clusters']),
        // A distribution read from the run's one configuration (the user's
        // rule, 2026-10-24): its kind, what it is of, and the configurations
        // it comes from (always 1: every analysis reads the one .rmc6f).
        shape: { kind: 'string', subject: 'string', configurations: 'number' },
        reads: Object.freeze([DISTRIBUTION_STAT, DISTRIBUTION_NOUN, SINGLE_CONFIGURATION]),
        contradicts: singleConfigurationRule
    }),
    Object.freeze({
        key: 'lowering',
        from: Object.freeze(['measure_symmetry_lowering']),
        shape: {
            parent: 'string', lower: 'string|null',
            shifted: [{ element: 'string', shift_A: 'number', direction_uvw: oneOf(['number', 'number', 'number'], 'null') }]
        },
        reads: Object.freeze([]),
        // No check reads a lowering's shifts yet.
        contradicts: null
    })
]);

/** The facts keys, in the order the checks run. */
export const FACT_KEYS = Object.freeze(FACTS.map((fact) => fact.key));
