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
//                 F-43m lost its centrosymmetry); and a group–subgroup relation
//                 that does not hold — not a subgroup, the wrong index, t or k
//                 (subgroup; the user's rule, 2026-10-24);
//   direction     a preferred displacement direction when no directions
//                 analysis of the question found one — a finding of the
//                 user's tiers (2026-10-24: ≥ 3σ, ≥ 1.2× isotropic, held at the
//                 neighbouring resolution), a hint called a hint passing
//                 (displacement-anisotropy; a significant whole-map test alone
//                 is anisotropy, not a direction: review L3);
//   sites         a split or flat-topped site called disorder when no
//                 analyze_modulations of the question tested its element
//                 (disorder; the user's rule, 2026-10-24);
//   modulation    a superstructure (cell doubling, zone-boundary order) when
//                 every analyze_modulations of the question found no ordered
//                 point (a general wavevector's order allows "superstructure",
//                 never "cell doubling"); and order claimed at a box-scale
//                 wavevector, a long wave (every |q_i| ≤ the engine's
//                 long_wave_q), whatever its power (the user's rule,
//                 2026-10-24: a long wave says nothing about order at its q);
//                 and a cell doubling the box cannot hold — odd along an axis
//                 it doubles — or along axes no ordered zone-boundary point
//                 doubles (cell-doubling; the user's rule, 2026-10-24); and
//                 the diffuse verdict read as the user reads it (2026-10-24:
//                 only antiphase or general-q order is order): "no order"
//                 when order was found (no-order), order between cells when
//                 only in-phase correlation — correlated motion — was
//                 (correlated-motion), no correlation when any was
//                 (no-correlation);
//   cluster       cluster claims (distorted, breathing, rigid) every
//                 analyze_clusters of the element contradicts (review L5);
//   coordination  a coordination number every whole-shell count of the
//                 element contradicts (review L5); and a coordination number
//                 or a bond angle quoted without its window in Å (window; the
//                 user's rule, 2026-10-24: every CN or angle carries it);
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
import { isPointSubgroup, pointGroupFacts, pointGroupOfSymbol, pointGroupOrder } from '../../context/pointGroups';
import {
    clausesOf, elementPattern, groupsIn, namesElement, namesSymbol, normalizeSymbol, positiveMatches, qVectorsIn, symbolSpelling, symbolTokens, windowsIn
} from './reading';
import {
    ALL_AXES, ALONG_AXES, ARROW, BOND_ANGLE_WORDS, BREATHING, CANNOT_HOLD, CELL_DOUBLING, CELL_ORDER, CENTROSYMMETRIC, CHIRAL, CN_CLAIMS, CONVERGED,
    CORRELATION_SUBJECT, DEGREES, DISORDER, DISTORTED, DISTRIBUTION_NOUN, DISTRIBUTION_STAT, DOUBLED_AXIS, DOUBLED_CELL, HALF_TUPLE, HAS_AS_SUBGROUP,
    INDEX_CLAIM, K_TYPE, LONG_WAVE_WORDS, LOST_INVERSION, LOWERING_TO, MAXIMAL, NEIGHBOUR_COUNTS, NOT_A_BOND_ANGLE, NO_CORRELATION, NO_ORDER,
    ORDER_CLAIM, PARTIAL_SHELL, PER_COUNT, POLAR, POLYHEDRON_CN, PREFERRED, REFERENCE_VALUE, RIGID, RUN_NOUN, SAYS_NOTHING, SINGLE_CONFIGURATION,
    SITE_REFS, SPLIT_SITE, STILL_MOVING, SUBGROUP_OF, SUPERGROUP_OF, SUPERSTRUCTURE, T_TYPE, WINDOW_CUTOFF, WINDOW_RANGES, WORD_NUMBER
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
    const ladder = symmetry?.ladder || [];
    // A group's operations in the cell and its place on the ladder (tight →
    // loose), for the subgroup check.
    const counted = (rung) => ({
        ...(Number.isFinite(rung?.n_ops) ? { n_ops: rung.n_ops } : {}),
        ...(ladder.includes(rung) ? { rank: ladder.indexOf(rung) } : {})
    });
    if (symmetry?.space_group && symmetry.point_group) {
        const facts = pointGroupFacts(symmetry.point_group);
        const symbol = String(symmetry.space_group).replace(/\s*\(No\.\s*\d+\)\s*$/, '');
        const rung = ladder.find((item) => item.sg && normalizeSymbol(item.sg) === normalizeSymbol(symbol));
        if (facts) groups.push({ symbol, ...facts, ...counted(rung ?? symmetry), parent: true });
    }
    const most = Math.max(...ladder.map((rung) => (Number.isFinite(rung.n_ops) ? rung.n_ops : -Infinity)));
    ladder.forEach((rung, index) => {
        const facts = rung.pg ? pointGroupFacts(rung.pg) : null;
        // Without operation counts, the ladder's two ends.
        const top = Number.isFinite(most) ? rung.n_ops === most : index === 0 || index === ladder.length - 1;
        if (facts && rung.sg) groups.push({ symbol: rung.sg, ...facts, ...counted(rung), ...(top ? { parent: true } : {}) });
    });
    return {
        convergence: context.convergence?.watchdog_status ?? null,
        groups,
        // The box (unit cells along a, b, c) and the run's elements, for the box check.
        supercell: Array.isArray(context.structure?.supercell) ? context.structure.supercell : null,
        elements: Object.keys(context.structure?.elements || {})
    };
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

// --- group–subgroup relations (the local-symmetry interview, 2026-10-24) ----------------

// "An answer that claims a group–subgroup relation that doesn't hold gets a
// warning under it" (the user, 2026-10-24, answer 4a). A relation the answer
// asserts — "H is a (maximal) (t-/k-) subgroup of G", "G is a supergroup of
// H", "G has H as a subgroup", or a lowering from G to H ("lowers G to H",
// "from G to H", "G → H" in a clause that says lower, reduce, descend or
// break down) — is read against the run's own groups where it knows both:
// the ladder's rungs are nested (each holds every operation of the tighter
// ones, symmetry.js groupsByThreshold), so H is a subgroup of G when it is
// tighter on the ladder with fewer operations of the cell, which G's count
// divides; the index is |G|/|H| in the same cell, and the type is t
// (translationengleiche) when the translations per cell (operations over the
// point group's order) stay, k (klassengleiche) when the point group stays,
// else neither (a general subgroup, which is never maximal). A group the run
// does not know (a published phase) is read by its point group alone (its
// symbol's, pointGroupOfSymbol): a subgroup's point group is one of the
// parent's subgroups (isPointSubgroup), a k-subgroup keeps it, a t-subgroup
// lowers it by its index. A negated, hypothetical or asked relation claims
// nothing.
const translationsPerCell = (group) => (Number.isFinite(group.n_ops) && pointGroupOrder(group.point_group) ? group.n_ops / pointGroupOrder(group.point_group) : null);

// The run's groups for the relation: each with its point group, and, where
// the ladder or detect_symmetry gives them, its operations in the cell and
// its rank on the ladder (the context's entry and a tool's merged).
const relationGroups = (facts, context) => {
    const groups = new Map();
    const add = (group) => {
        const key = normalizeSymbol(group.symbol);
        if (!key || !/^[pifcabr]/.test(key) || !group.point_group) return;
        const known = groups.get(key) || {};
        groups.set(key, {
            ...known,
            symbol: known.symbol ?? group.symbol,
            point_group: known.point_group ?? group.point_group,
            ...(Number.isFinite(group.n_ops) && !Number.isFinite(known.n_ops) ? { n_ops: group.n_ops } : {}),
            ...(Number.isFinite(group.rank) && !Number.isFinite(known.rank) ? { rank: group.rank } : {})
        });
    };
    for (const group of context?.groups || []) add(group);
    for (const fact of facts) for (const group of fact.groups || []) add(group);
    return groups;
};

// The relations a clause asserts: [{ sub, parent, said }], each group a
// symbol token (reading.js symbolTokens) with its group or its point group,
// `said` the words that state its index and type: the stretch from the first
// of the two symbols to a little past the second ("H is a t-subgroup of G of
// index 3"), or, for a step of a lowering, what follows the step's subgroup
// up to the next symbol ("G → H (t3) → K (k2)").
const SAID_AFTER = 30;
const relationsIn = (clause, groups, sentence = clause) => {
    const tokens = symbolTokens(clause, groups);
    if (tokens.length < 2) return [];
    const text = symbolSpelling(clause);
    const before = (at) => [...tokens].reverse().find((token) => token.end <= at) ?? null;
    const after = (at) => tokens.find((token) => token.at >= at) ?? null;
    const stretch = (a, b) => (a && b ? text.slice(Math.min(a.at, b.at), Math.max(a.end, b.end) + SAID_AFTER) : '');
    const pairs = [];
    positiveMatches(text, SUBGROUP_OF).forEach((match) => {
        const at = match.index + match[1].length;
        const [sub, parent] = [before(at), after(at + match[2].length)];
        pairs.push({ sub, parent, said: stretch(sub, parent) });
    });
    positiveMatches(text, SUPERGROUP_OF).forEach((match) => {
        const at = match.index + match[1].length;
        const [parent, sub] = [before(at), after(at + match[2].length)];
        pairs.push({ sub, parent, said: stretch(sub, parent) });
    });
    positiveMatches(text, HAS_AS_SUBGROUP).forEach((match) => {
        const at = match.index + match[1].length;
        const parent = before(at);
        const found = after(at + match[2].length);
        const sub = found && found.at < match.index + match[0].length ? found : null;
        pairs.push({ sub, parent, said: stretch(sub, parent) });
    });
    // A lowering from one symbol to the next, in a sentence that says lower,
    // reduce, descend or break down: only an arrow or "to" (with a tolerance
    // or an index in brackets, say) between them.
    if (clausesOf(symbolSpelling(sentence)).some((part) => positiveMatches(part, LOWERING_TO).length)) {
        tokens.slice(1).forEach((token, k) => {
            const between = text.slice(tokens[k].end, token.at).replace(/\([^()]{0,30}\)/g, ' ');
            if (!ARROW.test(between)) return;
            const next = tokens[k + 2];
            pairs.push({ parent: tokens[k], sub: token, said: text.slice(token.end, next ? next.at : token.end + SAID_AFTER) });
        });
    }
    const asGroup = (token) => ({
        ...(token.group || {}), symbol: token.group?.symbol ?? token.written, point_group: token.group?.point_group ?? pointGroupOfSymbol(token.written)
    });
    return pairs.filter((pair) => pair.sub && pair.parent && pair.sub.key !== pair.parent.key)
        .map(({ sub, parent, said }) => ({ sub: asGroup(sub), parent: asGroup(parent), said }));
};

const subgroupRule = (reading, facts, context) => {
    const groups = relationGroups(facts, context);
    for (const sentence of reading.sentences) {
        for (const clause of clausesOf(sentence)) {
            for (const { sub, parent, said } of relationsIn(clause, groups, sentence)) {
                const claimedIndex = INDEX_CLAIM.exec(said);
                const index = claimedIndex ? Number(claimedIndex[1] ?? claimedIndex[2]) : null;
                const claimsT = T_TYPE.test(said);
                const claimsK = K_TYPE.test(said) && !claimsT;
                const claimsMaximal = MAXIMAL.test(said);
                const pgSub = sub.point_group;
                const pgParent = parent.point_group;
                const says = `${sub.symbol} as a subgroup of ${parent.symbol}`;
                const finding = (text2) => [{ rule: 'subgroup', text: `The answer has ${says}, but ${text2}` }];
                if (Number.isFinite(sub.n_ops) && Number.isFinite(parent.n_ops)) {
                    // Both the run's: the ladder's nesting, its operations in the cell.
                    const ladderOrder = Number.isFinite(sub.rank) && Number.isFinite(parent.rank) ? sub.rank < parent.rank : true;
                    const holds = sub.n_ops < parent.n_ops && parent.n_ops % sub.n_ops === 0 && ladderOrder && isPointSubgroup(pgSub, pgParent) !== false;
                    if (!holds) {
                        return finding(`on the ladder ${sub.symbol} holds ${sub.n_ops} operations of the cell and ${parent.symbol} ${parent.n_ops}: `
                            + `${sub.n_ops >= parent.n_ops ? `${parent.symbol} is the subgroup, not ${sub.symbol}` : `${parent.symbol} does not contain ${sub.symbol}`}.`);
                    }
                    const actual = parent.n_ops / sub.n_ops;
                    const pgRatio = pointGroupOrder(pgParent) / pointGroupOrder(pgSub);
                    const [tSub, tParent] = [translationsPerCell(sub), translationsPerCell(parent)];
                    const type = tSub === tParent ? 't' : (pgSub === pgParent ? 'k' : 'general');
                    const how = type === 't' ? `the translations stay (${tParent} per cell) and the point group lowers from ${pgParent} to ${pgSub}: a t-subgroup`
                        : type === 'k' ? `the point group stays (${pgParent}) and the translations drop from ${tParent} to ${tSub} per cell: a k-subgroup`
                            : `the translations drop from ${tParent} to ${tSub} per cell and the point group from ${pgParent} to ${pgSub}: a general `
                                + `subgroup, neither t nor k (index ${actual} = ${pgRatio} × ${tParent / tSub})`;
                    if (index != null && index !== actual) return finding(`its index is ${parent.n_ops}/${sub.n_ops} = ${actual} operations of the cell, not ${index}: ${how}.`);
                    if (claimsT && type !== 't') return finding(`it is not translationengleiche: ${how}.`);
                    if (claimsK && type !== 'k') return finding(`it is not klassengleiche: ${how}.`);
                    if (claimsMaximal && type === 'general') return finding(`${how}, and a subgroup that is neither t nor k is never maximal (it lies below an intermediate group).`);
                    continue;
                }
                // A group the run does not know: its point group.
                const contained = isPointSubgroup(pgSub, pgParent);
                if (contained === false) {
                    return finding(`${sub.symbol}'s point group, ${pgSub}, is not a subgroup of ${parent.symbol}'s, ${pgParent}: no subgroup of ${parent.symbol} has it.`);
                }
                if (contained && claimsK && pgSub !== pgParent) {
                    return finding(`a klassengleiche subgroup keeps the point group, and ${sub.symbol}'s (${pgSub}) is not ${parent.symbol}'s (${pgParent}).`);
                }
                if (contained && claimsT && pgSub === pgParent) {
                    return finding(`a translationengleiche subgroup lowers the point group, and ${sub.symbol} keeps ${parent.symbol}'s (${pgParent}).`);
                }
                const pgRatio = contained ? pointGroupOrder(pgParent) / pointGroupOrder(pgSub) : null;
                if (contained && claimsT && index != null && index !== pgRatio) {
                    return finding(`a translationengleiche subgroup's index is the ratio of the point groups' orders, ${pointGroupOrder(pgParent)}/${pointGroupOrder(pgSub)} = ${pgRatio}, not ${index}.`);
                }
            }
        }
    }
    return [];
};

// --- direction ------------------------------------------------------------------------

// A preferred direction needs a finding of the user's tiers (2026-10-24: ≥ 3σ,
// ≥ 1.2× isotropic, and it holds at the neighbouring resolution;
// facts.direction.tier, analyze_displacement_directions' directionTier); facts
// without a tier read their peak test (≥ 3σ). A clause that calls it a hint
// claims no finding.
const HINT_WORDS = /\bhints?\b|\btentative(?:ly)?\b|\bnot\s+(?:yet\s+)?a\s+finding\b/i;
const isFinding = (analysis) => (analysis.tier ? analysis.tier === 'finding' : Boolean(analysis.significant));
const TIER_TEXT = {
    below_sigma: 'between 2.5σ and 3σ', small_enhancement: 'under 1.2× isotropic', not_held: 'not held at the neighbouring resolution',
    not_auto: 'neither at the engine\'s resolution nor beside it', unchecked: 'not checked at a second resolution', below_hint_sigma: 'below 2.5σ'
};

const directionRule = (reading, facts) => {
    const analyses = facts.map((fact) => fact.direction).filter(Boolean);
    if (!analyses.length || analyses.some(isFinding)) return [];
    for (const sentence of reading.sentences) {
        if (clausesOf(sentence).some((clause) => !HINT_WORDS.test(clause) && positiveMatches(clause, PREFERRED).length)) {
            const detail = analyses.map((analysis) => `${analysis.subject}: peak ${analysis.peak_sigma}σ`
                + `${Number.isFinite(analysis.enhancement) ? ` at ${analysis.enhancement}× isotropic` : ''}`
                + `${analysis.map_sigma != null ? `, map ${analysis.map_sigma}σ` : ''}`
                + `${analysis.tier ? ` — ${analysis.tier === 'hint' ? `a hint, ${TIER_TEXT[analysis.because] ?? 'not a finding'}` : 'nothing to report'}` : ''}`).join('; ');
            const anisotropic = analyses.some((analysis) => analysis.anisotropic);
            const tiered = analyses.some((analysis) => analysis.tier);
            return [{
                rule: 'direction',
                text: (tiered
                    ? 'The answer claims a preferred displacement direction, but no directions analysis in it found one (a finding: ≥ 3σ, '
                        + `≥ 1.2× isotropic, and held at the neighbouring resolution; ${detail}).`
                    : `The answer claims a preferred displacement direction, but no directions analysis in it found a peak at ≥ 3σ (${detail}).`)
                    + (anisotropic ? ' A significant map test is anisotropy, which any thermal ellipsoid shows, not a preferred direction.' : '')
            }];
        }
    }
    return [];
};

// --- a split or flat-topped site called disorder (the displacement-anisotropy interview, 2026-10-24) ---

// "An answer that calls a split or flat-topped site 'disorder' when
// analyze_modulations has not run for its element gets a warning under it"
// (the user, 2026-10-24, answer 8): ordered antiphase displacements give a
// site the same split cloud as random two-site disorder; only the
// correlation between cells tells them apart. A clause that asserts disorder
// (DISORDER — not a ranking, "the most disordered site", nor thermal disorder
// — or the gate's split-site words, SPLIT_SITE) about sites — named by number
// ("site 3"), else by an element the clause (or its sentence) names — of
// which list_sites found one split (the app's rule) or flat-topped (a
// negative Mardia kurtosis that misses it) is flagged when no
// analyze_modulations of the question tested that site's element. Only once
// list_sites ran (its facts name the split and flat-topped sites).
const siteRefsIn = (text) => [...String(text).matchAll(SITE_REFS)].flatMap((match) => match[1].split(/\s*(?:,|and|&)\s*/).map(Number));

const disorderRule = (reading, facts) => {
    const sites = facts.flatMap((fact) => fact.sites || []);
    const shaped = sites.filter((site) => site.split || site.flat);
    if (!shaped.length) return [];
    const tested = new Set(facts.flatMap((fact) => fact.modulation?.elements || []));
    const elements = [...new Set(sites.flatMap((site) => site.elements || [site.element]))];
    for (const sentence of reading.sentences) {
        for (const clause of clausesOf(sentence)) {
            if (!positiveMatches(clause, DISORDER).length && !positiveMatches(clause, SPLIT_SITE).length) continue;
            const refs = siteRefsIn(clause);
            const named = elements.filter((element) => namesElement(clause, element));
            const about = refs.length
                ? shaped.filter((site) => refs.includes(site.ref))
                : shaped.filter((site) => (named.length ? named : elements.filter((element) => namesElement(sentence, element)))
                    .some((element) => (site.elements || [site.element]).includes(element)));
            const untested = about.filter((site) => !(site.elements || [site.element]).some((element) => tested.has(element)));
            if (!untested.length) continue;
            const which = untested.slice(0, 3).map((site) => `site ${site.ref} (${site.element}, ${site.split ? 'split' : 'flat-topped'})`).join(', ');
            const missing = [...new Set(untested.map((site) => site.element))];
            return [{
                rule: 'disorder',
                text: `The answer calls ${which} disorder, but analyze_modulations has not run for ${missing.join(', ')} in this question: ordered `
                    + 'antiphase displacements (a superstructure) give a site the same split cloud as random disorder, and only the '
                    + 'correlation between cells tells them apart.'
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

// --- a cell doubling the box cannot hold (the local-symmetry interview, 2026-10-24) -----

// "An answer that claims a cell doubling the box's average structure can't
// support gets a warning under it" (the user, 2026-10-24, answer 4b). A
// clause that asserts a cell doubling (CELL_DOUBLING, DOUBLED_CELL) along
// some axes — "along a and b", "the c axis is doubled", "2a × 2a × c", a ½
// component of a wavevector — is flagged when the box (the run's supercell,
// unit cells along a, b, c) is odd along one of them, or along every axis
// when the clause names none: a doubled cell does not fit in it (its ½
// wavevector is not on the box's grid), so the configuration can neither show
// nor rule out that doubling — in a sentence about this run (an element of
// it, or the run itself), unless the answer says the box cannot hold it. And,
// once an analyze_modulations found ordered zone-boundary points, a doubling
// along axes none of them is ½ on (at that q): the superstructure check
// flags a doubling with no ordered point at all.
const AXIS_INDEX = { a: 0, b: 1, c: 2 };
const isHalf = (value) => Math.abs(Math.abs(value) % 1 - 0.5) < 1e-6;
const HALF_VALUE = { '½': 0.5, '1/2': 0.5, '.5': 0.5, '0.5': 0.5, 0: 0 };
const doublingAxes = (clause) => {
    const axes = new Set();
    const letter = (word) => (word ? /[abc]/.exec(word.toLowerCase())?.[0] : null);
    if (ALL_AXES.test(clause)) [0, 1, 2].forEach((k) => axes.add(k));
    const along = ALONG_AXES.exec(clause);
    if (along) along.slice(1).map(letter).filter(Boolean).forEach((axis) => axes.add(AXIS_INDEX[axis]));
    for (const match of clause.matchAll(DOUBLED_AXIS)) [match[1], match[2], match[3]].map(letter).filter(Boolean).forEach((axis) => axes.add(AXIS_INDEX[axis]));
    for (const q of qVectorsIn(clause)) q.forEach((value, k) => { if (isHalf(value)) axes.add(k); });
    for (const match of clause.matchAll(HALF_TUPLE)) match.slice(1, 4).forEach((value, k) => { if (HALF_VALUE[value] === 0.5) axes.add(k); });
    return axes;
};
const axisNames = (axes) => [...axes].sort().map((k) => 'abc'[k]).join(', ');

const cellDoublingRule = (reading, facts, context) => {
    const supercell = context?.supercell;
    const ordered = facts.map((fact) => fact.modulation).filter(Boolean).flatMap((fact) => fact.ordered || []);
    if (!Array.isArray(supercell) && !ordered.length) return [];
    const elements = context?.elements || [];
    const aboutRun = (sentence) => !elements.length || RUN_NOUN.test(sentence) || elements.some((element) => namesElement(sentence, element));
    const cannotHold = CANNOT_HOLD.test(reading.text);
    for (const sentence of reading.sentences) {
        for (const clause of clausesOf(sentence)) {
            if (!positiveMatches(clause, CELL_DOUBLING).length && !positiveMatches(clause, DOUBLED_CELL).length) continue;
            const axes = doublingAxes(clause);
            if (Array.isArray(supercell) && !cannotHold && aboutRun(sentence)) {
                const odd = [0, 1, 2].filter((k) => !(supercell[k] % 2 === 0));
                const box = supercell.join(' × ');
                if (axes.size && [...axes].some((k) => odd.includes(k))) {
                    const along = [...axes].filter((k) => odd.includes(k));
                    return [{
                        rule: 'cell-doubling',
                        text: `The answer claims a cell doubling along ${axisNames(axes)}, but the box (${box} unit cells) is odd along `
                            + `${axisNames(along)}: a doubled cell there does not fit in it (its ½ wavevector is not on the box's grid), so `
                            + 'the configuration can neither show nor rule out that doubling.'
                    }];
                }
                if (!axes.size && odd.length === 3) {
                    return [{
                        rule: 'cell-doubling',
                        text: `The answer claims a cell doubling, but the box (${box} unit cells) is odd along every axis: no doubled cell `
                            + 'fits in it, so the configuration can neither show nor rule out a doubling.'
                    }];
                }
            }
            if (axes.size && ordered.length && !ordered.some(({ q }) => [...axes].every((k) => isHalf(q[k])))) {
                return [{
                    rule: 'cell-doubling',
                    text: `The answer claims a cell doubling along ${axisNames(axes)}, but analyze_modulations found ordered zone-boundary `
                        + `points only at ${[...new Set(ordered.map(({ q }) => qText(q)))].join(', ')}: none doubles the cell along ${axisNames(axes)}.`
                }];
            }
        }
    }
    return [];
};

// --- order between cells as the diffuse verdict decides it (the local-order interview, 2026-10-24) ---

// "Only antiphase or general-q diffuse order past the bar counts as order;
// in-phase correlation is neighbouring cells moving together (correlated
// motion), not order, and never supports a 'no correlation' statement" (the
// user, 2026-10-24, answers 5 and 6). facts.modulation.diffuse holds each
// element whose diffuse verdict passed its bar, with its phase. Order found is
// an ordered zone-boundary point, an ordered general mode, or antiphase
// diffuse order; a correlation found is any of them, or in-phase diffuse
// correlation. Three checks, once an analyze_modulations ran:
//   no-order           "no order" (NO_ORDER, the gate's words) when order was
//                      found — unless the sentence itself names order between
//                      cells or antiphase ("no order beyond the antiphase
//                      correlation");
//   correlated-motion  order between cells (CELL_ORDER) when none was found:
//                      in phase only, correlated motion, not order;
//   no-correlation     no correlation of displacements or cells
//                      (NO_CORRELATION; "uncorrelated" about
//                      CORRELATION_SUBJECT) when a correlation was found.
const diffuseOf = (analyses) => analyses.flatMap((analysis) => analysis.diffuse || []);
const lagText = (lag) => `[${lag.join(' ')}]`;
const signedCorrelation = (value) => (value < 0 ? `−${Math.abs(value)}` : `+${value}`);
const diffuseText = ({ element, phase, lag, correlation, sigma }) => `${element} ${phase === 'antiphase' ? 'in antiphase' : 'in phase'} between `
    + `cells ${lagText(lag)} apart (${signedCorrelation(correlation)}, ${sigma >= 8 ? '> 8' : sigma}σ)`;
const orderFound = (analyses) => [
    ...analyses.flatMap((analysis) => (analysis.ordered || []).map(({ element, q }) => `${element} ordered at q = ${qText(q)}`)),
    ...analyses.flatMap((analysis) => (analysis.general || []).map(({ element, q }) => `${element} ordered at the general q = ${qText(q)}`)),
    ...diffuseOf(analyses).filter((item) => item.phase === 'antiphase').map(diffuseText)
];
const NAMES_CELL_ORDER = /\bantiphase\b/i;

const noOrderRule = (reading, facts) => {
    const found = orderFound(facts.map((fact) => fact.modulation).filter(Boolean));
    if (!found.length) return [];
    for (const sentence of reading.sentences) {
        if (positiveMatches(sentence, CELL_ORDER).length || NAMES_CELL_ORDER.test(sentence)) continue;
        if (clausesOf(sentence).some((clause) => positiveMatches(clause, NO_ORDER).length)) {
            return [{
                rule: 'no-order',
                text: `The answer says there is no order, but analyze_modulations found order between cells in this question: ${found.slice(0, 3).join('; ')}.`
            }];
        }
    }
    return [];
};

const correlatedMotionRule = (reading, facts) => {
    const analyses = facts.map((fact) => fact.modulation).filter(Boolean);
    if (!analyses.length || orderFound(analyses).length) return [];
    const inPhase = diffuseOf(analyses).filter((item) => item.phase === 'in_phase');
    for (const sentence of reading.sentences) {
        if (clausesOf(sentence).some((clause) => positiveMatches(clause, CELL_ORDER).length)) {
            return [{
                rule: 'correlated-motion',
                text: inPhase.length
                    ? `The answer claims order between cells, but the only correlation between cells analyze_modulations found in this question is in phase (${inPhase.slice(0, 3).map(diffuseText).join('; ')}): `
                        + 'neighbouring cells moving together (correlated motion), not order; only antiphase or general-q order past its bar is order.'
                    : 'The answer claims order between cells, but no analyze_modulations of this question found any: no ordered point or mode, and no '
                        + 'correlation between neighbouring cells past its bar.'
            }];
        }
    }
    return [];
};

const noCorrelationRule = (reading, facts) => {
    const analyses = facts.map((fact) => fact.modulation).filter(Boolean);
    const found = [...orderFound(analyses), ...diffuseOf(analyses).filter((item) => item.phase === 'in_phase').map(diffuseText)];
    if (!found.length) return [];
    for (const sentence of reading.sentences) {
        for (const clause of clausesOf(sentence)) {
            const claims = positiveMatches(clause, NO_CORRELATION).filter((match) => !/^uncorrelated$/i.test(match[2] ?? '') || CORRELATION_SUBJECT.test(clause));
            if (!claims.length) continue;
            return [{
                rule: 'no-correlation',
                text: `The answer says the displacements are not correlated, but analyze_modulations found a correlation between cells past its bar in this question: ${found.slice(0, 3).join('; ')}`
                    + (found.length && !orderFound(analyses).length ? ' (in phase: correlated motion, not order, but correlation all the same).' : '.')
            }];
        }
    }
    return [];
};

const modulationRules = (reading, facts, context) => [
    ...superstructureRule(reading, facts, context), ...longWaveRule(reading, facts), ...cellDoublingRule(reading, facts, context),
    ...noOrderRule(reading, facts), ...correlatedMotionRule(reading, facts), ...noCorrelationRule(reading, facts)
];

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

// The coordination counts of the question's analyses (find_gr_peaks: each
// way; compute_bond_angles: its window), one list.
const coordinationCounts = (facts) => facts
    .flatMap((fact) => (Array.isArray(fact.coordination) ? fact.coordination : (fact.coordination ? [fact.coordination] : [])))
    .filter((count) => Number.isFinite(count.cn));

// A window the answer states is the shell's when its top agrees within this
// (Å; the tools round to 0.01 Å, an answer may round to 0.1 Å) and it starts
// no later than the shell's onset (a cutoff "within 2.52 Å" counts what the
// shell holds: g is nought below the onset).
const WINDOW_MATCH_A = 0.05;
const sameWindow = (stated, shell) => Array.isArray(stated) && Array.isArray(shell)
    && Math.abs(stated[1] - shell[1]) <= WINDOW_MATCH_A && stated[0] <= shell[0] + WINDOW_MATCH_A;

// A coordination number claimed for an element, against every whole-shell
// count of it (find_gr_peaks' first shell, compute_bond_angles in a whole
// shell) — not against a window the user chose, whose count stands: a clause
// about part of the shell (the first peak, the short bonds) or one stating a
// window in Å other than the whole shell's is left (one stating the whole
// shell's window is a claim about it, which the window check asks every CN
// to carry). A claim matches a count within half an atom (shells.js
// cnMatches).
const coordinationRule = (reading, facts) => {
    const counts = coordinationCounts(facts);
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
            const stated = windowsIn(clause);
            // "n A per B" names its centre (B) and partner (A); a CN phrase over the
            // same number ("the coordination is 3.99 Se per Ga") is read through it.
            const perCounts = positiveMatches(clause, PER_COUNT).map((match) => ({
                n: countOf(match[2]), at: match.index, centre: match[4], partner: match[3], from: match.index + match[1].length, to: match.index + match[0].length
            }));
            const overlapsPer = (match) => perCounts.some((per) => per.from >= match.index && per.from < match.index + match[0].length);
            const claimed = [
                ...CN_CLAIMS.flatMap((pattern) => positiveMatches(clause, pattern).filter((match) => !overlapsPer(match))
                    .map((match) => ({ n: countOf(match[2]), at: match.index }))),
                ...POLYHEDRON_CN.flatMap(([pattern, cn]) => positiveMatches(clause, pattern).map((match) => ({ n: cn, at: match.index }))),
                ...perCounts
            ].filter((claim) => Number.isFinite(claim.n));
            for (const { n, at, centre: named = null, partner: given = null } of claimed) {
                const centre = named ?? nearestCentre(clause, at);
                if (!centre) continue;
                const mine = counts.filter((count) => count.centre === centre);
                const partners = given ? [given] : [...new Set(mine.map((count) => count.partner))].filter((partner) => partner !== centre && namesElement(clause, partner));
                const relevant = partners.length === 1 ? mine.filter((count) => count.partner === partners[0]) : mine;
                const whole = relevant.filter((count) => count.whole_shell);
                if (!whole.length) continue;
                // A window of the answer's own that is not the whole shell's: its count stands.
                if (stated.length && !stated.some((window) => whole.some((count) => sameWindow(window, count.window)))) continue;
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

// --- every CN and angle with its window (the coordination-shells interview, 2026-10-24) -

// "Every CN or angle in an answer carries its window in Å; an answer that gives
// one without it gets a warning under it" (the user, 2026-10-24, answer 8). A
// coordination number (CN_CLAIMS' "6-coordinated" and "CN 6", NEIGHBOUR_COUNTS'
// "3.99 Se per Ga" and "4 neighbours"; a polyhedron named in words names no
// number) or a bond angle (degrees in a sentence about angles, in a clause
// that is not about a direction's cell radius or cone, nor an ideal value:
// "Se–Ga–Se angles: mean 109.4°") carries its window when its
// sentence states a window in Å (reading.js windowsIn), or when another
// sentence of the answer does that names no element or every element the
// quoting sentence names — the window stated once for the pair, as an answer
// that opens with it reads. A quoting sentence that names no element is
// carried by any window in the answer. Only once an analysis that counts a
// shell or measures angles ran (find_gr_peaks with a first shell,
// compute_bond_angles: their coordination facts), whose windows the finding
// names.
const quotedCounts = (sentence) => {
    const aboutAngles = BOND_ANGLE_WORDS.test(sentence);
    return clausesOf(sentence).flatMap((clause) => [
        ...[...CN_CLAIMS, ...NEIGHBOUR_COUNTS].flatMap((pattern) => positiveMatches(clause, pattern)
            .map((match) => `a coordination number ("${match[0].slice(match[1].length).trim()}")`)),
        ...(aboutAngles && !NOT_A_BOND_ANGLE.test(clause) && !REFERENCE_VALUE.test(clause)
            ? positiveMatches(clause, DEGREES).map((match) => `a bond angle (${match[2]}°)`)
            : [])
    ]);
};

const pairWindowText = (count) => (Array.isArray(count.window) ? `${count.centre}–${count.partner} ${count.window[0]}–${count.window[1]} Å` : null);

const windowRule = (reading, facts) => {
    const counts = coordinationCounts(facts);
    if (!counts.length) return [];
    const elements = [...new Set(counts.flatMap((count) => [count.centre, count.partner]))];
    const sentences = reading.sentences.map((sentence) => ({
        sentence, windows: windowsIn(sentence), named: elements.filter((element) => namesElement(sentence, element))
    }));
    const stating = sentences.filter((item) => item.windows.length);
    for (const { sentence, windows, named } of sentences) {
        if (windows.length) continue;
        const quoted = quotedCounts(sentence);
        if (!quoted.length) continue;
        const carried = named.length
            ? stating.some((item) => !item.named.length || named.every((element) => item.named.includes(element)))
            : stating.length > 0;
        if (carried) continue;
        const windowsUsed = [...new Set(counts.map(pairWindowText).filter(Boolean))];
        return [{
            rule: 'window',
            text: `The answer gives ${quoted[0]} without its window in Å: every coordination number and bond angle carries the window `
                + `it was counted in${windowsUsed.length ? ` (this question's: ${windowsUsed.join('; ')})` : ''}.`
        }];
    }
    return [];
};

// The coordination rule, then the window check (the user's, 2026-10-24).
const coordinationRules = (reading, facts) => [...coordinationRule(reading, facts), ...windowRule(reading, facts)];

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
        // The group and its point group's facts (context/pointGroups.js), its
        // operations in the cell and the tolerance it holds at; the run
        // context's ladder adds its groups (the parents flagged, each with its
        // operations and its rank on the ladder).
        shape: [{
            symbol: 'string', point_group: 'string', centrosymmetric: 'boolean', sohncke: 'boolean', polar: 'boolean', crystal_system: 'string',
            'n_ops?': 'number', 'tolerance_A?': 'number'
        }],
        reads: Object.freeze([
            CHIRAL, CENTROSYMMETRIC, POLAR, ...LOST_INVERSION, SUBGROUP_OF, SUPERGROUP_OF, HAS_AS_SUBGROUP, LOWERING_TO, ARROW, INDEX_CLAIM, T_TYPE, K_TYPE,
            MAXIMAL
        ]),
        // The point-group facts, then the group–subgroup relation (the user's, 2026-10-24).
        contradicts: (reading, facts, context) => [...symmetryRules(reading, facts, context), ...subgroupRule(reading, facts, context)]
    }),
    Object.freeze({
        key: 'direction',
        from: Object.freeze(['analyze_displacement_directions']),
        // The peak's tier (the user's, 2026-10-24: directionTier) and why,
        // its enhancement and the resolution it was read at.
        shape: {
            subject: 'string', significant: 'boolean', anisotropic: 'boolean', peak_sigma: 'number', map_sigma: 'number|null', tier: 'string',
            because: 'string', enhancement: 'number|null', frequency: 'number'
        },
        reads: Object.freeze([PREFERRED]),
        contradicts: directionRule
    }),
    Object.freeze({
        key: 'sites',
        from: Object.freeze(['list_sites']),
        // Each site's species, and whether it is split (the app's rule) or
        // flat-topped (a negative Mardia kurtosis that misses it); with
        // analyze_modulations' tested elements, the disorder check (the
        // user's rule, 2026-10-24).
        shape: [{ ref: 'number', element: 'string', elements: ['string'], split: 'boolean', flat: 'boolean' }],
        reads: Object.freeze([DISORDER, SPLIT_SITE, SITE_REFS]),
        contradicts: disorderRule
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
            long_waves: [{ element: 'string', q: ['number', 'number', 'number'] }],
            // The elements whose diffuse verdict passed its bar, in phase or antiphase (the local-order interview, 2026-10-24).
            'diffuse?': [{ element: 'string', phase: 'string', lag: ['number', 'number', 'number'], correlation: 'number', sigma: 'number' }]
        },
        reads: Object.freeze([
            CELL_DOUBLING, SUPERSTRUCTURE, ORDER_CLAIM, LONG_WAVE_WORDS, SAYS_NOTHING, DOUBLED_CELL, ALONG_AXES, DOUBLED_AXIS, ALL_AXES, HALF_TUPLE,
            CANNOT_HOLD, RUN_NOUN, NO_ORDER, CELL_ORDER, NO_CORRELATION, CORRELATION_SUBJECT
        ]),
        // The superstructure rule, the long-wave rule, the box rule and the diffuse verdict's three (the user's, 2026-10-24).
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
        reads: Object.freeze([
            ...CN_CLAIMS, ...POLYHEDRON_CN.map(([pattern]) => pattern), PARTIAL_SHELL, PER_COUNT, ...WINDOW_RANGES, WINDOW_CUTOFF,
            ...NEIGHBOUR_COUNTS, DEGREES, BOND_ANGLE_WORDS, NOT_A_BOND_ANGLE, REFERENCE_VALUE
        ]),
        // The coordination rule, then the window check (the user's, 2026-10-24).
        contradicts: coordinationRules
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
