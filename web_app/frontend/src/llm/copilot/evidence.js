// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The evidence gate (THY-18): an answer that concludes something must have
// run, in this turn, the analysis that gives the evidence for it. A model
// cannot be forced to call a tool (tool_choice is not available on every
// model the Copilot runs), so this is a check after the answer: the agent
// loop sends an answer whose evidence did not run back to the model once per
// rule, with what to run (agentLoop.js), and an answer that still lacks it is
// "not verified" (verdict.js), with the missing analysis named under it.
//
// The rules, each read as the answer checks read claims (answerAudit.js
// `asserts`: not negated, not hypothetical, not a question):
//
//   order between cells — a superstructure, cell doubling, modulation or
//     incommensurate order needs analyze_modulations;
//   split sites — random two-site disorder or off-centring needs
//     analyze_modulations, since ordered antiphase displacements give a site
//     the same cloud (the displacement-anisotropy skill);
//   direction — a preferred displacement direction needs
//     analyze_displacement_directions;
//   clusters — distorted, breathing or rigid clusters need analyze_clusters;
//   coordination — a coordination number, a split shell or a Jahn–Teller
//     distortion needs find_gr_peaks or compute_bond_angles;
//   lowering — a lowering of the average structure (a coherent shift)
//     needs measure_symmetry_lowering (detect_symmetry where it is not offered);
//   no order — "no order beyond thermal motion" needs analyze_modulations,
//     analyze_displacement_directions on a site, and analyze_clusters for every
//     element whose own first shell is a finite cluster's (the local-order skill).
//
// A conclusion is about this run when its sentence names one of the run's
// elements (a symbol, or in a formula: Ti, Mo3, Ta₄, TaSe6) or the run itself
// (this structure, the configuration, the clusters, site 3); a definition ("A
// superstructure is order that enlarges the cell") names neither and needs no
// analysis. "No order" is always about the run. Without the run's elements
// every sentence counts.
//
// Pure, no React. Words of the validation's verdict schema (off_centre,
// cluster_breathing, "features": [] …) are read as the claims they are.

import { BREATHING, CELL_DOUBLING, CN_CLAIMS, DISTORTED, POLYHEDRON_CN, PREFERRED, RIGID, SUPERSTRUCTURE, asserts } from './answerAudit';

const MODULATION = /(^|[^\w-])(modulat(?:ion|ions|ed)|incommensurate|order(?:ed)? between cells|correlated displacements between cells)\b/gi;
const SPLIT_SITE = /(^|[^\w-])(split[- ]sites?|two[- ]site(?:\s+disorder)?|disordered over two|random(?:ly)?\s+(?:\w+\s+)?(?:off[- ]cent\w*|displac\w*|distribut\w*)|off[- ]cent(?:re|er|ring|ering)\s+(?:dis)?order|off_centre)\b/gi;
const CLUSTER_KINDS = /(^|[^\w-])(cluster_distortion|cluster_breathing|rigid_cluster)\b/gi;
const SHELL_CLAIMS = /(^|[^\w-])(jahn[-–_ ]teller|split (?:first )?(?:coordination )?shell|short and long bonds|jahn_teller)\b/gi;
const LOWERING = /(^|[^\w-])(symmetry[-_ ]lowering|lower(?:s|ed|ing)? (?:the )?(?:average )?symmetry|coherent(?:ly)? (?:shift\w*|off[- ]cent\w*|displace\w*)|polar shift)\b/gi;
const NO_ORDER = /(^|[^\w-])(no (?:local |short-range |additional )?order(?!\s+(?:that|which)\s+doubles)|no (?:ordered|correlated) (?:displacements?|distortions?)|(?:only|purely) thermal motion|thermal motion only|nothing beyond (?:the )?(?:average structure|thermal motion)|no (?:sign|evidence) of (?:local )?order)\b/gi;
const EMPTY_VERDICT = /"features"\s*:\s*\[\s*\]/;
const RUN_NOUN = /\b(?:th(?:is|e) (?:run|configuration|structure|model|box|supercell|sample|compound|crystal|material|system)|the (?:sites?|atoms|clusters|octahedra|tetrahedra|polyhedra|displacements)|sites? \d+)\b/i;

// Every element whose own first shell holds at most this many neighbours may
// form finite clusters (a dimer 1, a trimer 2, a tetrahedron 3, an octahedron
// 4 per member): "no order" needs analyze_clusters on it.
export const CLUSTER_CN_MAX = 4.5;

// The tool calls that ran: [{ name, args }] → helpers.
const ranTool = (ran, name) => ran.some((call) => call.name === name);

/**
 * The gate for one run. `catalog` (clarify.js buildShellCatalog, optional)
 * names the elements whose own first shell is a finite cluster's; `offered`
 * the tool names this turn offers (a rule never asks for a tool that is not
 * there).
 */
export const buildEvidenceGate = ({ catalog = null, offered = null, elements = null } = {}) => {
    const has = (name) => !offered || offered.includes(name);
    // A sentence about this run (above); every sentence without the elements.
    const symbols = elements ? elements.map((element) => new RegExp(`(?<![A-Za-z])${element}(?![a-z])`)) : null;
    const within = symbols ? (sentence) => RUN_NOUN.test(sentence) || symbols.some((pattern) => pattern.test(sentence)) : null;
    const claims = (text, pattern) => asserts(text, pattern, { within });
    const clusterElements = catalog
        ? catalog.elements.filter((element) => {
            const shell = catalog.shell(element, element);
            return shell && Number.isFinite(shell.cn) && shell.cn >= 0.5 && shell.cn <= CLUSTER_CN_MAX;
        })
        : [];
    const rules = [
        {
            id: 'order-between-cells',
            what: 'order between cells',
            tools: ['analyze_modulations'],
            applies: (text) => claims(text, SUPERSTRUCTURE) || claims(text, CELL_DOUBLING) || claims(text, MODULATION),
            missing: (ran) => (ranTool(ran, 'analyze_modulations') ? [] : ['analyze_modulations'])
        },
        {
            id: 'split-site',
            what: 'random (two-site) disorder of a split site',
            tools: ['analyze_modulations'],
            applies: (text) => claims(text, SPLIT_SITE),
            missing: (ran) => (ranTool(ran, 'analyze_modulations') ? [] : ['analyze_modulations'])
        },
        {
            id: 'direction',
            what: 'a preferred displacement direction',
            tools: ['analyze_displacement_directions'],
            applies: (text) => claims(text, PREFERRED),
            missing: (ran) => (ranTool(ran, 'analyze_displacement_directions') ? [] : ['analyze_displacement_directions'])
        },
        {
            id: 'clusters',
            what: 'the clusters\' distortion, breathing or rigidity',
            tools: ['analyze_clusters'],
            applies: (text) => claims(text, DISTORTED) || claims(text, BREATHING) || claims(text, RIGID) || claims(text, CLUSTER_KINDS),
            missing: (ran) => (ranTool(ran, 'analyze_clusters') ? [] : ['analyze_clusters'])
        },
        {
            id: 'coordination',
            what: 'a coordination number or a split shell',
            tools: ['find_gr_peaks', 'compute_bond_angles'],
            applies: (text) => CN_CLAIMS.some((pattern) => claims(text, pattern))
                || POLYHEDRON_CN.some(([pattern]) => claims(text, pattern)) || claims(text, SHELL_CLAIMS),
            missing: (ran) => (ranTool(ran, 'find_gr_peaks') || ranTool(ran, 'compute_bond_angles') ? [] : ['find_gr_peaks'])
        },
        {
            id: 'lowering',
            what: 'a lowering of the average structure',
            tools: ['measure_symmetry_lowering', 'detect_symmetry'],
            applies: (text) => claims(text, LOWERING),
            missing: (ran) => {
                const needed = has('measure_symmetry_lowering') ? 'measure_symmetry_lowering' : 'detect_symmetry';
                return ranTool(ran, needed) ? [] : [needed];
            }
        },
        {
            id: 'no-order',
            what: '"no order beyond thermal motion"',
            tools: ['analyze_modulations', 'analyze_displacement_directions', 'analyze_clusters'],
            applies: (text) => asserts(text, NO_ORDER) || EMPTY_VERDICT.test(text),
            missing: (ran) => [
                ranTool(ran, 'analyze_modulations') ? null : 'analyze_modulations',
                ran.some((call) => call.name === 'analyze_displacement_directions' && call.args?.site_ref != null)
                    ? null : 'analyze_displacement_directions on each kind of site (site_ref)',
                ...clusterElements
                    .filter((element) => !ran.some((call) => call.name === 'analyze_clusters' && call.args?.element === element))
                    .map((element) => `analyze_clusters for ${element} (its first ${element}–${element} shell is a finite cluster's)`)
            ].filter(Boolean)
        }
    ].filter((rule) => rule.tools.some(has));

    // What an answer concludes without its evidence: [{ id, what, missing }].
    const unmet = (text, ran) => rules
        .filter((rule) => rule.applies(String(text || '')))
        .map((rule) => ({ id: rule.id, what: rule.what, missing: rule.missing(ran || []).filter((item) => has(item.split(' ')[0])) }))
        .filter((rule) => rule.missing.length);
    return {
        rules,
        clusterElements,
        unmet,
        /** The first unmet rule not sent back yet this turn, as the note for the model; or null. */
        check: (text, { ran = [], fired = new Set() } = {}) => {
            const rule = unmet(text, ran).find((item) => !fired.has(item.id));
            if (!rule) return null;
            return {
                ...rule,
                nudge: `(From the app, not the user.) Your answer concludes ${rule.what}, but this question has not run `
                    + `${rule.missing.join(', ')}. Run it, then answer from its result (quote its numbers by reference). If it cannot `
                    + 'run, say so in the answer instead of concluding.',
                notice: `The answer concluded ${rule.what} without ${rule.missing.join(', ')}; the app asked for it first.`
            };
        }
    };
};
