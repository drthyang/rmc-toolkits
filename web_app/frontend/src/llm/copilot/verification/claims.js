// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The claim catalogue (CLAIMS, THY-18's evidence gate): one row per kind of
// conclusion that needs an analysis in this question, in the order the gate
// sends an answer back. A model cannot be forced to call a tool (tool_choice
// is not available on every model the Copilot runs), so this is read after
// the answer: an answer that concludes something whose analysis did not run
// goes back to the model (verification index.js judge), and an answer that
// still lacks it is "not verified" with the missing analysis named under it.
//
// Each row is data:
//   id, what     the claim, as the send-back and the finding say it;
//   scope        'run': read only in a sentence about this run — one that
//                names an element of the run (Ti, Mo3, Ta₄, TaSe6) or the run
//                itself (this structure, the configuration, the clusters, site
//                3); a definition ("A superstructure is order that enlarges the
//                cell") names neither and needs no analysis; without the run's
//                elements every sentence counts. 'always': about the run in
//                any sentence ("no order");
//   phrasings    { asserted, literal }: `asserted` patterns are read as claims
//                (reading.js asserts: not negated, hypothetical or asked);
//                `literal` ones match anywhere in the text. A caller's own
//                phrasings join them (createVerification's `phrasings`: the
//                validation's verdict words, validation/verdict.js);
//   backedBy     the analyses that back it, through one combinator — anyOf
//                (met when any ran; asks for the first the turn offers),
//                firstOffered (the first the turn offers, else the last, is
//                the one needed) or allOf (each requirement, every `each`
//                value its own). A
//                requirement is a tool name or { tool, with, each, say }:
//                `with` names the argument the call must carry (`each`: equal
//                to each value), `say` how a missing one is named. A row is
//                checked only when the turn offers one of its tools, and never
//                asks for a tool the turn does not offer;
//   prompt       { phrase, analysis }: how EVIDENCE_RULES (prompts/copilot.js,
//                still written by hand) states it; evidenceSentence renders
//                the sentence from the rows and promptGaps lists what it does
//                not say (the binding test).
//
// Under the turn's method 'off' (the validation's arm C) the gate checks the
// same rows, but its send-back names no analysis and no tool (OFF_NUDGE):
// which analysis a conclusion needs is the method. The notice and `missing`
// (the step card, the record) still name them; the model never reads those.
// Pure, no React.

import { asserts } from './reading';
import {
    BREATHING, CELL_DOUBLING, CN_CLAIMS, DISTORTED_CLUSTERS, DISTORTED_POLYHEDRA, LOWERING, MODULATION, NO_ORDER, POLYHEDRON_CN, PREFERRED, RIGID,
    RUN_NOUN, SHELL_CLAIMS, SPLIT_OR_FLAT_DISORDER, SPLIT_SITE, SUPERSTRUCTURE
} from './phrasings';

// --- the combinators ----------------------------------------------------------------

const requirement = (spec) => Object.freeze(typeof spec === 'string' ? { tool: spec, say: spec } : { ...spec, say: spec.say ?? spec.tool });
const combinator = (need) => (...specs) => Object.freeze({ need, requirements: Object.freeze(specs.map(requirement)) });

/** Met when any requirement ran; a missing one is named by the first. */
export const anyOf = combinator('any');
/** The first requirement the turn offers (else the last) is the one needed. */
export const firstOffered = combinator('first');
/** Every requirement is needed, and a requirement with `each` once per value. */
export const allOf = combinator('all');

// --- the catalogue ------------------------------------------------------------------

const row = (spec) => Object.freeze({
    ...spec,
    phrasings: Object.freeze({ asserted: Object.freeze(spec.phrasings.asserted || []), literal: Object.freeze(spec.phrasings.literal || []) }),
    prompt: spec.prompt ? Object.freeze(spec.prompt) : null
});

// An element whose own first shell is a finite cluster's (shells.js
// clusterElementsOf) needs analyze_clusters before "no order"; one that
// centres octahedra or tetrahedra no finite cluster holds (networkCentresOf:
// an extended network, or isolated polyhedra) polyhedral_distortion, where
// the turn offers it (the local-order interview, 2026-10-24).
const clusterSay = (element) => `analyze_clusters for ${element} (its first ${element}–${element} shell is a finite cluster's)`;
const networkSay = (element) => `polyhedral_distortion for ${element} (centre ${element}: its polyhedra form no finite cluster)`;

export const CLAIMS = Object.freeze([
    // Order between cells — a superstructure, cell doubling, modulation or
    // incommensurate order.
    row({
        id: 'order-between-cells',
        what: 'order between cells',
        scope: 'run',
        phrasings: { asserted: [SUPERSTRUCTURE, CELL_DOUBLING, MODULATION] },
        backedBy: anyOf('analyze_modulations'),
        prompt: { phrase: 'order between cells', analysis: 'analyze_modulations' }
    }),
    // Split sites — random two-site disorder or off-centring: ordered
    // antiphase displacements give a site the same cloud (the
    // displacement-anisotropy skill); a split or flat-topped site called
    // disorder too (its interview, 2026-10-24).
    row({
        id: 'split-site',
        what: 'random (two-site) disorder of a split site',
        scope: 'run',
        phrasings: { asserted: [SPLIT_SITE, SPLIT_OR_FLAT_DISORDER] },
        backedBy: anyOf('analyze_modulations'),
        prompt: { phrase: 'a split site called disorder', analysis: 'analyze_modulations' }
    }),
    row({
        id: 'direction',
        what: 'a preferred displacement direction',
        scope: 'run',
        phrasings: { asserted: [PREFERRED] },
        backedBy: anyOf('analyze_displacement_directions'),
        prompt: { phrase: 'a preferred direction', analysis: 'analyze_displacement_directions' }
    }),
    row({
        id: 'clusters',
        what: 'the clusters\' distortion, breathing or rigidity',
        scope: 'run',
        phrasings: { asserted: [DISTORTED_CLUSTERS, BREATHING, RIGID] },
        backedBy: anyOf('analyze_clusters'),
        prompt: { phrase: 'clusters', analysis: 'analyze_clusters' }
    }),
    // Distorted coordination polyhedra (octahedra, tetrahedra): a distortion
    // index (polyhedral_distortion, where the run offers it) or the cluster
    // null — a corner-sharing polyhedron or another extended network has no
    // finite clusters to test (the coordination-shells interview, 2026-10-24,
    // answer 7). Without polyhedral_distortion the gate asks for analyze_clusters.
    row({
        id: 'polyhedra',
        what: 'a distortion of the coordination polyhedra',
        scope: 'run',
        phrasings: { asserted: [DISTORTED_POLYHEDRA] },
        backedBy: anyOf('polyhedral_distortion', 'analyze_clusters'),
        prompt: { phrase: 'distorted polyhedra', analysis: 'polyhedral_distortion' }
    }),
    // A coordination number, a split shell or a Jahn–Teller distortion.
    row({
        id: 'coordination',
        what: 'a coordination number or a split shell',
        scope: 'run',
        phrasings: { asserted: [...CN_CLAIMS, ...POLYHEDRON_CN.map(([pattern]) => pattern), SHELL_CLAIMS] },
        backedBy: anyOf('find_gr_peaks', 'compute_bond_angles'),
        prompt: { phrase: 'a coordination number or split shell', analysis: 'find_gr_peaks' }
    }),
    // A lowering of the average structure (a coherent shift):
    // measure_symmetry_lowering, detect_symmetry where it is not offered.
    row({
        id: 'lowering',
        what: 'a lowering of the average structure',
        scope: 'run',
        phrasings: { asserted: [LOWERING] },
        backedBy: firstOffered('measure_symmetry_lowering', 'detect_symmetry'),
        prompt: { phrase: 'a lowering of the average structure', analysis: 'measure_symmetry_lowering' }
    }),
    // "No order beyond thermal motion": the superstructure test (its diffuse
    // verdict included: analyze_modulations), a per-site directions map,
    // analyze_clusters for every element whose own first shell is a finite
    // cluster's, and polyhedral_distortion for every centre of polyhedra no
    // finite cluster holds, counted only where the turn offers it (the
    // local-order skill; its interview, 2026-10-24).
    row({
        id: 'no-order',
        what: '"no order beyond thermal motion"',
        scope: 'always',
        phrasings: { asserted: [NO_ORDER] },
        backedBy: allOf(
            'analyze_modulations',
            { tool: 'analyze_displacement_directions', with: 'site_ref', say: 'analyze_displacement_directions on each kind of site (site_ref)' },
            { tool: 'analyze_clusters', with: 'element', each: 'clusterElement', say: clusterSay },
            { tool: 'polyhedral_distortion', with: 'centre', each: 'networkCentre', say: networkSay }
        ),
        prompt: { phrase: '"no order"', analysis: 'analyze_modulations' }
    })
]);

/** The claim ids, in send-back order. */
export const CLAIM_IDS = Object.freeze(CLAIMS.map((claim) => claim.id));

// --- the send-back ------------------------------------------------------------------

// The send-back under method 'off', after "Your answer concludes <what>, but":
// the evidence a conclusion needs, with no analysis or tool named (decided by
// the manager on the user's delegation, 2026-10-19).
export const OFF_NUDGE = 'this conclusion needs its evidence from an analysis run in this question. Answer from '
    + 'the analysis\'s result (quote its numbers by reference). If none can run, say so in the answer instead of concluding.';

/** An unmet claim ({ id, what, missing }) as the gate sends it back under `method`. */
export const sendBackOf = ({ id, what, missing }, method) => ({
    id,
    what,
    missing,
    nudge: method === 'off'
        ? `(From the app, not the user.) Your answer concludes ${what}, but ${OFF_NUDGE}`
        : `(From the app, not the user.) Your answer concludes ${what}, but this question has not run `
            + `${missing.join(', ')}. Run it, then answer from its result (quote its numbers by reference). If it cannot `
            + 'run, say so in the answer instead of concluding.',
    notice: `The answer concluded ${what} without ${missing.join(', ')}; the app asked for it first.`
});

/** An unmet claim as the finding under the answer. */
export const evidenceFinding = ({ id, what, missing }) => ({
    rule: 'evidence',
    id,
    missing,
    text: `The answer concludes ${what}, but ${missing.join(', ')} did not run in this question.`
});

// --- reading the claims of one turn -------------------------------------------------

// Whether a call (one that ran: { name, args }) meets a requirement instance.
const meets = (call, item) => call.name === item.tool
    && (!item.with || (item.each ? call.args?.[item.with] === item.value : call.args?.[item.with] != null));

/**
 * The claims one turn checks: the rows a tool it offers can back (`offered`:
 * tool names, null for every tool), each with `extra` phrasings joined
 * ({ [id]: { asserted, literal } }), read in sentences about the run
 * (`elements`, the run's), with `each` values from `values`
 * ({ clusterElement: [...], networkCentre: [...] }). Returns { rows, unmet(text, ran) }: unmet is
 * [{ id, what, missing }] in row order, `missing` naming only offered tools.
 */
export const claimsOf = ({ offered = null, extra = null, elements = [], values = {} } = {}) => {
    const has = (name) => !offered || offered.includes(name);
    // A sentence about this run (above); every sentence without the elements.
    const symbols = elements.length ? elements.map((element) => new RegExp(`(?<![A-Za-z])${element}(?![a-z])`)) : null;
    const within = symbols ? (sentence) => RUN_NOUN.test(sentence) || symbols.some((pattern) => pattern.test(sentence)) : null;
    const rows = CLAIMS
        .filter((claim) => claim.backedBy.requirements.some((item) => has(item.tool)))
        .map((claim) => ({
            claim,
            asserted: [...claim.phrasings.asserted, ...(extra?.[claim.id]?.asserted || [])],
            literal: [...claim.phrasings.literal, ...(extra?.[claim.id]?.literal || [])]
        }));
    const applies = ({ claim, asserted, literal }, text) => {
        const scope = claim.scope === 'run' ? within : null;
        return asserted.some((pattern) => asserts(text, pattern, { within: scope })) || literal.some((pattern) => text.search(pattern) >= 0);
    };
    // The requirement instances of a row: an `each` requirement once per value.
    const instances = ({ requirements }) => requirements.flatMap((item) => (item.each
        ? (values[item.each] || []).map((value) => ({ ...item, value, say: typeof item.say === 'function' ? item.say(value) : item.say }))
        : [item]));
    const missingOf = ({ backedBy }, ran) => {
        const ranAny = (item) => ran.some((call) => meets(call, item));
        let missing;
        if (backedBy.need === 'any') {
            // Asks for the first the turn offers (a row is read only when one is).
            missing = backedBy.requirements.some(ranAny) ? [] : [backedBy.requirements.find((item) => has(item.tool)) ?? backedBy.requirements[0]];
        } else if (backedBy.need === 'first') {
            const needed = backedBy.requirements.find((item) => has(item.tool)) ?? backedBy.requirements.at(-1);
            missing = ranAny(needed) ? [] : [needed];
        } else {
            missing = instances(backedBy).filter((item) => !ranAny(item));
        }
        return missing.filter((item) => has(item.tool)).map((item) => item.say);
    };
    return {
        rows: rows.map(({ claim }) => claim),
        unmet: (text, ran = []) => rows
            .filter((entry) => applies(entry, String(text || '')))
            .map(({ claim }) => ({ id: claim.id, what: claim.what, missing: missingOf(claim, ran) }))
            .filter((claim) => claim.missing.length)
    };
};

// --- the prompt sentence and the binding test ---------------------------------------

const need = (claim) => claim.need ?? claim.backedBy?.need;

const EVIDENCE_HEAD = 'Every conclusion needs its analysis in this question: ';
const EVIDENCE_TAIL = '. The app sends an answer back when its analysis did not run.';

/**
 * The prompt sentence the rows state (EVIDENCE_RULES' form): consecutive rows
 * with the same analysis are joined ("order between cells or a split site
 * called disorder, analyze_modulations"), and an allOf row reads "…, every one
 * of them". Rows without a prompt are left out.
 */
export const evidenceSentence = (claims) => {
    const groups = [];
    claims.filter((claim) => claim.prompt).forEach((claim) => {
        const analysis = need(claim) === 'all' ? 'every one of them' : claim.prompt.analysis;
        const last = groups.at(-1);
        if (last && last.analysis === analysis) last.phrases.push(claim.prompt.phrase);
        else groups.push({ analysis, phrases: [claim.prompt.phrase] });
    });
    return `${EVIDENCE_HEAD}${groups.map((group) => `${group.phrases.join(' or ')}, ${group.analysis}`).join('; ')}${EVIDENCE_TAIL}`;
};

const label = (item) => `${item.tool}${item.with ? ` with ${item.with}` : ''}${item.each ? ` for each ${item.each}` : ''}`;

/**
 * What the prompt sentence does not say of each row's backing, by claim id
 * (rows with nothing unsaid left out): an alternative or a fallback it does
 * not name (anyOf, firstOffered); for an allOf row ("every one of them": every
 * analysis the sentence names) each requirement with an argument it does not
 * state, and "not <tool>" for each analysis it names that the row does not
 * need. A row without a prompt lists every requirement.
 */
export const promptGaps = (claims) => {
    const named = [...new Set(claims.filter((claim) => claim.prompt && need(claim) !== 'all').map((claim) => claim.prompt.analysis))];
    const gaps = {};
    claims.forEach((claim) => {
        const { requirements } = claim.backedBy;
        let gap;
        if (!claim.prompt) {
            gap = requirements.map(label);
        } else if (need(claim) === 'all') {
            gap = [
                ...requirements.filter((item) => item.with || !named.includes(item.tool)).map(label),
                ...named.filter((tool) => !requirements.some((item) => item.tool === tool)).map((tool) => `not ${tool}`)
            ];
        } else {
            gap = requirements.filter((item) => item.tool !== claim.prompt.analysis || item.with).map(label);
        }
        if (gap.length) gaps[claim.id] = gap;
    });
    return gaps;
};
