// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The skill lint's tables (lint.js; docs/AGENT_PLAN.md item 5.1):
//
//   BINDINGS    a gate number a skill states → the constant(s) that decide it;
//   NOT_GATES   a number after a comparator that is no gate, with why;
//   FIRM_RULES  a "*Firm rule in the app:*" → the code that enforces it, and a
//               probe that exercises that code;
//   VALUES      a value a skill quotes as measured on the Demo → what the tool
//               measures there now (the shell readings, candidate #6);
//   RULES       a rule a skill states in words → the rule the app and arm A+
//               apply (candidate 8: the split site, the cluster elements, the
//               lower group's tolerance), and a probe that the code applies it;
//   ALLOW_LIST  today's violations, each with why and what clears it: a FACT
//               (the skill states a number or a rule the tool decides
//               otherwise; the skill fact-sync PR clears it), a JUDGMENT (a
//               method, an order or a judgment call: it stays for the user's
//               interviews, docs/SKILL_INTERVIEWS.md) or AWAITS CODE (an
//               interviewed skill states the user's rule, which the named
//               code PR still has to build).
//
// When a skill's text changes (the user's review of a draft), its entries here
// change in the same PR: a binding or a firm rule whose sentence moved is
// "stale", a violation the change fixes leaves the allow-list (the test fails
// on an entry nothing needs), and a new threshold or firm rule gets its entry.
// Files are named by their path under .claude/skills/.

import assert from 'node:assert/strict';
import {
    directionResolved, LOW_INDEX_CONE_DEG, SIGNIFICANCE_SIGMA, summarizeBondAngles, summarizeOrientation, summarizeSymmetry
} from '../../copilotTools';
import { buildClarifier } from '../../llm/copilot/clarify';
import { buildRunTools } from '../../llm/copilot/runTools';
import { toolResult } from '../../llm/copilot/tools';
import { createVerification, FACTS } from '../../llm/copilot/verification';
import { buildRunContext, lowerGroupTolerance, NOISE_FLOOR_MARGIN } from '../../llm/context/runContext';
import {
    CLUSTER_CN_MAX, CN_CHANGE_ABS, CN_CHANGE_REL, CN_MATCH_TOLERANCE, clusterElementsOf, cnMarginCut, cnMatches, CONSTRAINT_LIMIT_FRACTION,
    IDEAL_COORDINATION, isIdealCoordination, MODAL_COORDINATION_SHARE, SEPARATION_G_FRACTION, shellCatalogue, SPLIT_SHELL_REACH_A,
    UNBONDED_APEX_FRACTION
} from '../../llm/context/shells';
import {
    classifyConvergence, MIN_TOTAL_DROP, TREND_MIN_ROWS, TREND_NOISE_FACTOR, TREND_WINDOW_FRACTION
} from '../../llm/watchdog/heuristics';
import { SWEEP_RULES } from '../../../validation/armA';
import {
    AXIS_CONE_DEG, AXIS_ENHANCEMENT, CLUSTER_EFFECT, CLUSTER_Z, SPLIT_SITE_KURTOSIS, SPLIT_SITE_SE, TRIGONALITY_EFFECT, isSplitSite
} from '../../workers/clusters';
import { isLongWave, LONG_WAVE_Q, MODULATION_RATIO, MODULATION_RULES, MODULATION_SIGMA } from '../../workers/modulations';
import { ASYMMETRY_FLAG_SIGMA, MAP_TEST_MIN_PAIRS } from '../../workers/orientation';
import { AXIS_RESOLUTION_SIGMAS } from '../../workers/pcaKde';
import { coordinationSummaryOf } from '../../workers/triplets';

const ANISOTROPY = 'displacement-anisotropy/SKILL.md';
const SIGNIFICANCE = 'displacement-anisotropy/significance-reading.md';
const SHELLS = 'coordination-shells/SKILL.md';
const SHELL_TABLE = 'coordination-shells/lacunar-spinel-shells.md';
const LOCAL_ORDER = 'local-order/SKILL.md';
const ORDER_SIGNATURES = 'local-order/order-signatures.md';
const SYMMETRY = 'local-symmetry/SKILL.md';
const FIT = 'rmc-fit-assessment/SKILL.md';
const RUN_NUMBERS = 'rmc-fit-assessment/run-numbers.md';

const constant = (name, value, module) => ({ name, value, module });
const SIGMA = constant('SIGNIFICANCE_SIGMA', SIGNIFICANCE_SIGMA, 'copilotTools.js');
const ASYMMETRY_SIGMA = constant('ASYMMETRY_FLAG_SIGMA', ASYMMETRY_FLAG_SIGMA, 'workers/orientation.js');
const CLUSTER_SIGMA = constant('CLUSTER_Z', CLUSTER_Z, 'workers/clusters.js');
// compute_bond_angles' cut rule, the one rule the user chose on 2026-10-19
// (shells.js windowCut): the CN margin (cnMarginCut) wherever a window meets
// its shell, bound by both its parts; the coordination counts (countsCut) only
// where the run reads no shell to compare, the skill's fallback clause, whose
// thresholds keep their own constants.
const CN_MARGIN = constant('CN_CHANGE_REL', CN_CHANGE_REL, 'llm/context/shells.js, cnMarginCut via windowCut');
const CN_MARGIN_ABS = constant('CN_CHANGE_ABS', CN_CHANGE_ABS, 'llm/context/shells.js, cnMarginCut via windowCut');
const MODAL_SHARE = constant('MODAL_COORDINATION_SHARE', MODAL_COORDINATION_SHARE, 'llm/context/shells.js, countsCut (the fallback with no shell to compare)');
const UNBONDED_SHARE = constant('UNBONDED_APEX_FRACTION', UNBONDED_APEX_FRACTION, 'llm/context/shells.js, countsCut (the fallback with no shell to compare)');
// The split-site rule (workers/clusters.js isSplitSite: list_sites' delta_A, the clusters' split member sites).
const SPLIT_KAPPA = constant('SPLIT_SITE_KURTOSIS', SPLIT_SITE_KURTOSIS, 'workers/clusters.js, isSplitSite (list_sites; the cluster verdict\'s split member sites)');
const SPLIT_ERRORS = constant('SPLIT_SITE_SE', SPLIT_SITE_SE, 'workers/clusters.js, isSplitSite');
const PC1_RULE = { pattern: /\bPC1\b/, says: `the PC1 excess kurtosis on a resolved axis, at least ${SPLIT_SITE_SE} sampling errors below zero` };
// The order rule analyze_modulations renders (the engine decides it: each
// tested point's `ordered`, its rules in the payload's MODULATION_RULES).
const ORDER_SIGMA = constant('MODULATION_SIGMA', MODULATION_SIGMA, 'workers/modulations.js, the order rule');
const ORDER_RATIO = constant('MODULATION_RATIO', MODULATION_RATIO, 'workers/modulations.js, the order rule');
const LONG_WAVE = constant('LONG_WAVE_Q', LONG_WAVE_Q, 'workers/modulations.js, isLongWave');
// find_gr_peaks' join of a split shell (shells.js readShell).
const SPLIT_REACH = constant('SPLIT_SHELL_REACH_A', SPLIT_SHELL_REACH_A, 'llm/context/shells.js, readShell (find_gr_peaks\' join)');
const IDEAL_COUNTS = constant('IDEAL_COORDINATION', [...IDEAL_COORDINATION], 'llm/context/shells.js, isIdealCoordination');
const CN_MATCH = constant('CN_MATCH_TOLERANCE', CN_MATCH_TOLERANCE, 'llm/context/shells.js, cnMatches');
const FLOOR_MARGIN = constant('NOISE_FLOOR_MARGIN', NOISE_FLOOR_MARGIN, 'llm/context/runContext.js, lowerGroupTolerance and the Workbench\'s tolerance question');
const LOW_INDEX_CONE = constant('LOW_INDEX_CONE_DEG', LOW_INDEX_CONE_DEG, 'copilotTools.js, measure_symmetry_lowering and analyze_modulations\' polarisation');
const AXIS_CONE = constant('AXIS_CONE_DEG', AXIS_CONE_DEG, 'workers/clusters.js, analyze_clusters\' axis families');
// Both cones are inclusive: nearestLowIndex skips a direction only when it is
// more than the cone off; the axis families count |cos| ≥ cos 15°.
const INCLUSIVE_CONE = { pattern: /\bwithin\b|≤|\bat most\b/, says: 'an inclusive cone, a family named at 15° or less (nearestLowIndex: off > toleranceDeg is left out; clusters: |cos| ≥ cos 15°)' };
const PEAK_TEST = { pattern: /\bpeak\b/, says: 'the peak test only (analyze_displacement_directions facts.direction.significant; the map test is anisotropy)' };
// The CN margin's sentences name the CN they compare (#113).
const CN_MARGIN_STATISTIC = {
    pattern: /\bCN\b|coordination number|partner atoms/,
    says: 'the window\'s CN against the whole shell\'s (fewer by at least max(0.5, 15 % of its CN))'
};

// --- probes: the code applies a constant, or enforces a firm rule, as stated ----------

// An analyze_displacement_directions engine payload (no unit cell: no [u v w]).
const orientationPayload = ({ peak, map }) => ({
    referenceNumber: 1, element: 'Ta', frequency: 6, cellCount: 362, usedPoints: 4000, rejectedPoints: 0,
    vmin: 0.8, vmax: 1.3, peakEnhancement: 1.3, peakCount: 20, peakExpected: 11, peakSignificance: peak,
    peakDirection: null, mapSignificance: map, antipodalAsymmetry: 0, antipodalAsymmetryZ: 0,
    antipodalAsymmetrySignificant: false, orientationEigenvalues: [0.4, 0.3, 0.3], orientationAnisotropy: 0.1,
    orientationAnisotropySignificance: 1, recommendedFrequency: 6
});

// A compute_bond_angles engine payload with the coordination counts given
// (counts[n] = central atoms with n A–B bonds in the window) and the engine's
// summary of them.
const bondAnglePayload = (coordination) => ({
    triplet: ['Se', 'Ta', 'Se'], bond12: [2.3, 2.7], bond23: [2.3, 2.7], sharedEnds: true,
    angleCount: 1000, apexCount: coordination.reduce((sum, count) => sum + count, 0), coordination,
    coordinationSummary: coordinationSummaryOf(coordination),
    meanAngle: 100, stdAngle: 20, binCenters: [45, 90, 135, 180], counts: [100, 600, 200, 100], lengths12: null
});
const refused = (result) => typeof result.data.refused === 'string' && !('mean_deg' in result.data);
// A bond window's reading against its first shell (a 6-atom shell at
// 2.2–3.2 Å by default), as the shell catalogue's window() gives it, with the
// CN-margin verdict from the shell module itself.
const shellReading = (cn, shellCn = 6) => {
    const shell = { window: [2.2, 3.2], cn: shellCn };
    return { cn, holds: 'part', cnMargin: { cut: cnMarginCut(cn, shell, [2.3, 2.7]) }, shell };
};

// The answer judged by the verification (llm/copilot/verification/): the
// tools' `facts` recorded as steps, the run `context` given, no claim offered
// (only the checks on the answer speak). Its `done`: the findings and the
// Outcome they compute.
const judgeAnswer = (text, { facts = [], context = null } = {}) => {
    const verification = createVerification({ context, offered: [] });
    facts.forEach((item) => verification.record({ name: 'analysis', ...toolResult({ facts: item }) }, () => ({ data: null })));
    return verification.judge(text).done;
};
// One message of the answer check `rule` on `text`.
const findsRule = (text, options, rule) => {
    const rules = judgeAnswer(text, options).audit.map((finding) => finding.rule);
    assert.ok(rules.includes(rule), `the answer checks gave [${rules.join(', ')}] on "${text}", not "${rule}"`);
};
// A finding under the answer lowers the computed Outcome (the "warning under it").
const lowersOutcome = (text, options) => assert.equal(judgeAnswer(text, options).outcome.status, 'not verified');

// A χ² history still falling at the end: the watchdog says improving.
const fallingHistory = { name: 'run-01.log', plotData: { chiColumn: 'X_ray_(R)1', series: [{ y: Array.from({ length: 60 }, (_, i) => 5 - 0.05 * i) }] } };

// The Demo's ladder (local-symmetry's reference): F-43m over P-42₁m, neither centrosymmetric.
const demoSymmetry = {
    spaceGroup: 'F-43m', spaceGroupNumber: 216, pointGroup: '-43m', nSpace: 96, toleranceA: 0.2,
    ladder: [
        { spaceGroup: 'P-42_1m', pointGroup: '-42m', from: 0.02, to: 0.031, nSpace: 8 },
        { spaceGroup: 'F-43m', pointGroup: '-43m', from: 0.031, to: 1, nSpace: 96 }
    ]
};

// A run whose Cu–F partial g(r) is a CuF6-type split shell: 4 F at 1.93 Å and
// `second` F at `r2` (σ 0.05 Å), 200 Cu and 1200 F in a 20 Å cube — the
// Jahn–Teller octahedron of the coordination-shells skill (as in
// llm/__tests__/runTools.test.js). The trough between the peaks is deep, so
// only the join makes the shell whole.
const splitShellRun = ({ r2 = 2.25, second = 2 } = {}) => {
    const L = 20;
    const elementCounts = { Cu: 200, F: 1200 };
    const rhoF = elementCounts.F / L ** 3;
    const x = Array.from({ length: 600 }, (_, i) => i * 0.005);
    const normal = (r, r0) => Math.exp(-((r - r0) ** 2) / (2 * 0.05 ** 2)) / (0.05 * Math.sqrt(2 * Math.PI));
    const y = x.map((r) => (r > 0 ? (4 * normal(r, 1.93) + second * normal(r, r2)) / (4 * Math.PI * rhoF * r * r) : 0));
    return {
        structure: { latticeVectors: [[L, 0, 0], [0, L, 0], [0, 0, L]], elementCounts },
        plotFiles: [{ name: 'cuf.csv', plotKind: 'pdf_partials', plotData: { series: [{ label: 'Cu-F', x, y }] } }]
    };
};
const cuF = (options) => shellCatalogue(splitShellRun(options)).shell('Cu', 'F');
// An analysis with one A–B window (compute_bond_angles' shellWindows, in brief).
const windowedAnalysis = {
    shellWindows: (args) => [{ centre: 'Cu', partner: 'F', window: [args.lo, args.hi], apply: ([lo, hi]) => ({ ...args, lo, hi }) }]
};

// A Ga–Se peak at 2.42 Å (σ 0.1 Å, height 8) cut at its closest-approach
// constraint `limit`: g is 0 below it, and just above it the share of the
// first peak the cut's position sets (2.13 Å: 3 %; 2.18 Å: 9 %).
const constrainedGaSe = (limit) => {
    const x = Array.from({ length: 1000 }, (_, i) => Number((i * 0.01).toFixed(6)));
    const y = x.map((r) => (r <= limit ? 0 : 8 * Math.exp(-((r - 2.42) ** 2) / (2 * 0.1 ** 2))));
    return {
        plotFiles: [{ name: 'gase.csv', plotKind: 'pdf_partials', plotData: { series: [{ label: 'Ga-Se', x, y }] } }],
        runSettings: { atoms: ['Ga', 'Se'], minimumDistancesA: [3.0, limit, 3.0] }
    };
};
// The constraint-limited test (shells.js, the user's 5 %): 3 % of the peak at the
// cut is not limited (it was under the separation's 2 %), 9 % is, and
// find_gr_peaks says so with the user's advice.
const verifyConstraintLimit = async () => {
    const reading = (limit) => shellCatalogue(constrainedGaSe(limit)).shell('Ga', 'Se').constraint;
    const low = reading(2.13);
    assert.ok(low.gFraction > SEPARATION_G_FRACTION && low.gFraction < CONSTRAINT_LIMIT_FRACTION, `3 %: g_fraction ${low.gFraction}`);
    assert.equal(low.limited, false, '3 % of the peak at the cut was called constraint-limited');
    const high = reading(2.18);
    assert.ok(high.gFraction > CONSTRAINT_LIMIT_FRACTION, `9 %: g_fraction ${high.gFraction}`);
    assert.equal(high.limited, true, '9 % of the peak at the cut was not called constraint-limited');
    const tool = buildRunTools(constrainedGaSe(2.18)).find((item) => item.name === 'find_gr_peaks');
    const { checks } = await tool.run({ pair: 'Ga-Se' });
    assert.ok(checks.some((check) => check.level === 'warn'
        && /constraint-limited: interpret nothing on its low-r side; consider rerunning with the constraint loosened/.test(check.text)),
    'find_gr_peaks does not say constraint-limited with the user\'s advice');
};

// A run whose small F(Q)_1 term still falls (0.01 → 0.006 in χ² over the run)
// while a large flat constraint term, charted last, holds the total flat
// within its noise: the convergence rule on every term and the total.
const termStillFalling = () => {
    let seed = 17;
    const uniform = () => {
        seed = (seed * 16807) % 2147483647;
        return seed / 2147483647;
    };
    const noise = () => Math.sqrt(-2 * Math.log(uniform())) * Math.cos(2 * Math.PI * uniform());
    const n = 400;
    const small = Array.from({ length: n }, (_, t) => 0.01 * Math.exp(-t / 800) * Math.exp(0.002 * noise()));
    const large = Array.from({ length: n }, (_, t) => (t < 60 ? 20 - 0.15 * t : 11) * Math.exp(0.01 * noise()));
    const ln = (values) => values.map((value) => Math.log(value));
    return {
        name: 'run-00.log',
        plotData: {
            chiColumn: 'Curvature',
            series: [{ y: ln(large) }],
            chiTerms: { names: ['F(Q)_1', 'Curvature'], series: [ln(small), ln(large)], total: ln(small.map((value, t) => value + large[t])), segments: [n] }
        }
    };
};

// --- BINDINGS -------------------------------------------------------------------------

// `find` locates the sentence (in the text with its line breaks collapsed);
// `number` reads the number as written (group n) from `before` characters
// ahead of it to `after` beyond; `unit: '%'` reads 60 as 0.6; `statistic`,
// when given, is what the constant gates and the sentence must name;
// `verify` checks the code applies the constant as bound.
export const BINDINGS = [
    {
        file: ANISOTROPY,
        gate: 'the significance threshold',
        find: 'Read the significances, all calibrated one-sided σ',
        number: /\(≥ (?<n>\d+(?:\.\d+)?) is significant\)/,
        constants: [SIGMA, ASYMMETRY_SIGMA]
    },
    {
        file: ANISOTROPY,
        gate: 'a preferred direction',
        find: '"Preferred direction" only for a finding (step 4:',
        number: /significance ≥ (?<n>\d+(?:\.\d+)?)σ/,
        constants: [SIGMA],
        statistic: { pattern: /\bpeak\b/, says: 'the peak test only (analyze_displacement_directions facts.direction.significant; the map test is anisotropy)' },
        verify: () => {
            const thermal = summarizeOrientation(orientationPayload({ peak: SIGNIFICANCE_SIGMA - 0.7, map: 7.2 }), null).facts.direction;
            assert.equal(thermal.significant, false, 'a map test alone made a direction preferred');
            assert.equal(thermal.anisotropic, true);
            assert.equal(summarizeOrientation(orientationPayload({ peak: SIGNIFICANCE_SIGMA, map: 0.5 }), null).facts.direction.significant, true);
        }
    },
    {
        file: ANISOTROPY,
        gate: 'the low-index family cone',
        find: 'Name the nearest low-index family (<100>, <110>, <111>) only when the angle to it is',
        number: /only when the angle to it is (?:under|within|below|less than|≤|<) ?(?<n>\d+(?:\.\d+)?) ?°/,
        constants: [LOW_INDEX_CONE, AXIS_CONE],
        statistic: INCLUSIVE_CONE
    },
    {
        file: ANISOTROPY,
        gate: 'a resolved peak direction',
        find: 'and the cell is no wider (`peak.radius_deg`',
        number: /`peak\.radius_deg` ≤ (?<n>\d+(?:\.\d+)?)°/,
        constants: [constant('LOW_INDEX_CONE_DEG', LOW_INDEX_CONE_DEG, 'copilotTools.js, directionResolved (the peak cell\'s radius)')],
        verify: () => {
            assert.equal(directionResolved(LOW_INDEX_CONE_DEG), true, 'a 15° cell did not resolve the direction');
            assert.equal(directionResolved(LOW_INDEX_CONE_DEG + 0.1), false, 'a cell wider than 15° resolved the direction');
        }
    },
    {
        file: ANISOTROPY,
        gate: 'the low-index family cone (checklist)',
        find: 'the low-index family named only within',
        number: /named only within (?<n>\d+(?:\.\d+)?) ?°/,
        constants: [LOW_INDEX_CONE, AXIS_CONE]
    },
    {
        file: ANISOTROPY,
        gate: 'a split site',
        find: 'A site is split (two wells) only by the app\'s rule',
        after: 160,
        number: /whose excess kurtosis is ≤ (?<n>[−–-]?\d+(?:\.\d+)?) and/,
        constants: [SPLIT_KAPPA],
        statistic: PC1_RULE,
        verify: () => {
            const site = (kappa) => ({ referenceNumber: 1, count: 4000, axisResolved: [true, true, true], excessKurtosis: [kappa, 0, 0] });
            assert.equal(isSplitSite(site(SPLIT_SITE_KURTOSIS - 0.01)), true);
            assert.equal(isSplitSite(site(SPLIT_SITE_KURTOSIS + 0.01)), false, 'a kurtosis above the bound counted as split');
        }
    },
    {
        file: ANISOTROPY,
        gate: 'a split site\'s sampling errors',
        find: 'A site is split (two wells) only by the app\'s rule',
        after: 160,
        number: /and at least (?<n>\d+) sampling errors below zero/,
        constants: [SPLIT_ERRORS],
        verify: () => {
            // κ₁ −0.93 at 216 copies is 2.8 errors below zero: not split; at 4000 copies it is.
            const site = (count) => ({ referenceNumber: 1, count, axisResolved: [true, true, true], excessKurtosis: [-0.93, 0, 0] });
            assert.equal(isSplitSite(site(216)), false, '2.8 sampling errors made a site split');
            assert.equal(isSplitSite(site(4000)), true);
        }
    },
    {
        file: ANISOTROPY,
        gate: 'a preferred direction (step 3)',
        find: 'A preferred direction needs the peak test:',
        number: /`peak\.significance_sigma` ≥ (?<n>\d+(?:\.\d+)?)σ, unrounded/,
        constants: [SIGMA],
        statistic: PEAK_TEST,
        verify: () => {
            // On the unrounded value: 2.996σ prints as 3 at three digits, yet it is no preferred direction.
            const below = summarizeOrientation(orientationPayload({ peak: SIGNIFICANCE_SIGMA - 0.004, map: 0.5 }), null).facts.direction;
            assert.equal(below.significant, false, 'a peak just below 3σ made a direction preferred');
            assert.equal(below.peak_sigma, SIGNIFICANCE_SIGMA, 'the rounded σ no longer shows the case this guards');
        }
    },
    {
        file: ANISOTROPY,
        gate: 'the direction firm rule',
        find: 'every directions analysis of the question found none at',
        number: /found none at ≥ (?<n>\d+(?:\.\d+)?)σ/,
        constants: [SIGMA]
    },
    // The user's three tiers (interview 2026-10-24, answer 3): the σ of each
    // is the peak test's; the 1.2× and 2.5σ lines await their code (ALLOW_LIST).
    {
        file: ANISOTROPY,
        gate: 'a finding (step 4)',
        find: '**a finding** (a preferred direction):',
        number: /\(a preferred direction\): ≥ (?<n>\d+(?:\.\d+)?)σ and/,
        constants: [SIGMA]
    },
    {
        file: ANISOTROPY,
        gate: 'a hint (step 4)',
        find: '**a hint**:',
        number: /\*\*a hint\*\*: ≥ (?<n>\d+(?:\.\d+)?)σ under/,
        constants: [SIGMA]
    },
    {
        file: ANISOTROPY,
        gate: 'a cluster distortion\'s z (step 7)',
        find: 'gyration-axis families against the shuffle null, with effect size and z, at',
        number: /, at z ≥ (?<n>\d+(?:\.\d+)?) with the same sign/,
        constants: [CLUSTER_SIGMA]
    },
    {
        file: ANISOTROPY,
        gate: 'a cluster distortion\'s z (checklist)',
        find: 'Cluster-level claims rest on `analyze_clusters` (effect size and',
        number: /effect size and z ≥ (?<n>\d+(?:\.\d+)?),/,
        constants: [CLUSTER_SIGMA]
    },
    {
        file: SIGNIFICANCE,
        gate: 'the significance threshold',
        find: 'counts as significant',
        before: 14,
        number: /≥ (?<n>\d+(?:\.\d+)?)σ counts as significant/,
        // The direction tool's σ; the cluster tool's z is bound in its own paragraph.
        constants: [SIGMA, ASYMMETRY_SIGMA]
    },
    {
        file: SIGNIFICANCE,
        gate: 'a split site',
        find: 'A site is split (two positions) only by the app\'s rule on its PC1',
        number: /its excess kurtosis ≤ (?<n>[−–-]?\d+(?:\.\d+)?) and/,
        constants: [SPLIT_KAPPA],
        statistic: PC1_RULE
    },
    {
        file: SIGNIFICANCE,
        gate: 'a split site\'s sampling errors',
        find: 'A site is split (two positions) only by the app\'s rule on its PC1',
        number: /and at least (?<n>\d+) sampling errors below zero/,
        constants: [SPLIT_ERRORS]
    },
    {
        file: SIGNIFICANCE,
        gate: 'the map test\'s floor',
        find: 'expected coincident pairs',
        before: 24,
        number: /below ~?(?<n>\d+(?:\.\d+)?) expected coincident pairs/,
        constants: [constant('MAP_TEST_MIN_PAIRS', MAP_TEST_MIN_PAIRS, 'workers/orientation.js')]
    },
    {
        file: SIGNIFICANCE,
        gate: 'a resolved principal axis',
        find: 'an eigenvalue gap above',
        number: /gap above (?<n>\d+(?:\.\d+)?) standard errors/,
        constants: [constant('AXIS_RESOLUTION_SIGMAS', AXIS_RESOLUTION_SIGMAS, 'workers/pcaKde.js')]
    },
    {
        file: SIGNIFICANCE,
        gate: 'a cluster statistic\'s z',
        find: 'clear their thresholds (z',
        number: /\(z ≥ (?<n>\d+(?:\.\d+)?) and/,
        constants: [constant('CLUSTER_Z', CLUSTER_Z, 'workers/clusters.js')]
    },
    {
        file: SIGNIFICANCE,
        gate: 'a cluster statistic\'s effect size',
        find: 'a relative change',
        number: /a relative change ≥ (?<n>\d+(?:\.\d+)?) %/,
        unit: '%',
        constants: [constant('CLUSTER_EFFECT', CLUSTER_EFFECT, 'workers/clusters.js')]
    },
    {
        file: SIGNIFICANCE,
        gate: 'the T2 trigonality effect',
        find: 'T2 trigonality',
        number: /T2 trigonality ≥ (?<n>\d+(?:\.\d+)?)/,
        constants: [constant('TRIGONALITY_EFFECT', TRIGONALITY_EFFECT, 'workers/clusters.js')]
    },
    {
        file: SIGNIFICANCE,
        gate: 'the axis-family enhancement',
        find: 'axis enhancement',
        number: /axis enhancement ≥ (?<n>\d+(?:\.\d+)?) ?×/,
        constants: [constant('AXIS_ENHANCEMENT', AXIS_ENHANCEMENT, 'workers/clusters.js')]
    },
    {
        file: SHELLS,
        gate: 'a cut window: the CN margin\'s absolute part (firm rule)',
        find: 'holds fewer partner atoms than the whole first shell by at least max(',
        number: /by at least max\((?<n>\d+(?:\.\d+)?), \d+ % of the shell's CN\)/,
        constants: [CN_MARGIN_ABS],
        statistic: CN_MARGIN_STATISTIC,
        verify: () => {
            // Half an atom of a 2-atom shell (15 % of it is less), in the shell module and in compute_bond_angles.
            const shell = { window: [2.2, 3.2], cn: 2 };
            assert.equal(cnMarginCut(2 - CN_CHANGE_ABS - 0.01, shell, [2.3, 2.7]), true, 'a window 0.51 atoms short of a 2-atom shell was not cut');
            assert.equal(cnMarginCut(2 - CN_CHANGE_ABS + 0.01, shell, [2.3, 2.7]), false, 'a window 0.49 atoms short of a 2-atom shell was cut');
            const uniform = bondAnglePayload([0, 0, 0, 0, 0, 100]);
            assert.ok(refused(summarizeBondAngles(uniform, { shells: [shellReading(2 - CN_CHANGE_ABS - 0.01, 2)] })), 'a window 0.51 atoms short of a 2-atom shell was not refused');
            assert.ok(!refused(summarizeBondAngles(uniform, { shells: [shellReading(2 - CN_CHANGE_ABS + 0.01, 2)] })), 'a window 0.49 atoms short of a 2-atom shell was refused');
        }
    },
    {
        file: SHELLS,
        gate: 'a cut window: the CN margin\'s relative part (firm rule)',
        find: 'holds fewer partner atoms than the whole first shell by at least max(',
        number: /by at least max\(\d+(?:\.\d+)?, (?<n>\d+) % of the shell's CN\)/,
        unit: '%',
        constants: [CN_MARGIN],
        statistic: CN_MARGIN_STATISTIC,
        verify: () => {
            // 15 % of a 6-atom shell is 0.9 atoms; a window apart from the shell is not cut on the margin.
            const shell = { window: [2.2, 3.2], cn: 6 };
            const margin = Math.max(CN_CHANGE_ABS, CN_CHANGE_REL * 6);
            assert.equal(cnMarginCut(6 - margin - 0.01, shell, [2.3, 2.7]), true, 'a window 0.91 atoms short of a 6-atom shell was not cut');
            assert.equal(cnMarginCut(6 - margin + 0.01, shell, [2.3, 2.7]), false, 'a window 0.89 atoms short of a 6-atom shell was cut');
            assert.equal(cnMarginCut(0, shell, [3.3, 3.8]), false, 'a window apart from the shell was cut on the margin');
            // compute_bond_angles refuses on the margin wherever a shell is given, whatever its counts.
            const uniform = bondAnglePayload([0, 0, 0, 0, 0, 100]);
            const spread = bondAnglePayload([0, 0, 0, 0, 41, 0, 59]);
            assert.ok(refused(summarizeBondAngles(uniform, { shells: [shellReading(6 - margin - 0.01)] })), 'a window 0.91 atoms short of a 6-atom shell was not refused');
            assert.ok(!refused(summarizeBondAngles(spread, { shells: [shellReading(6 - margin + 0.01)] })), 'a window 0.89 atoms short of a 6-atom shell was refused on its counts');
            assert.ok(refused(summarizeBondAngles(bondAnglePayload([40, 0, 60]), { shells: [shellReading(1.2)] })), 'a window holding 1.2 of 6 was not refused');
        }
    },
    {
        file: SHELLS,
        gate: 'a cut window with no shell to compare (firm rule)',
        find: 'of the central atoms share the most common bond count',
        before: 24,
        number: /fewer than (?<n>\d+(?:\.\d+)?) % of the central atoms/,
        unit: '%',
        constants: [MODAL_SHARE],
        verify: () => {
            // With no shell given, the counts decide: 59 of 100 central atoms share the modal count: refused; 61: kept.
            assert.ok(refused(summarizeBondAngles(bondAnglePayload([0, 0, 0, 0, 41, 0, 59]))), 'a 59 % modal share was not refused');
            assert.ok(!refused(summarizeBondAngles(bondAnglePayload([0, 0, 0, 0, 39, 0, 61]))), 'a 61 % modal share was refused');
        }
    },
    {
        file: SHELLS,
        gate: 'a window that misses the shell, with no shell to compare (firm rule)',
        find: 'or most have none',
        number: /or (?<n>most) have none/,
        constants: [UNBONDED_SHARE],
        verify: () => {
            // With no shell given, "most" is more than half: 51 of 100 central atoms with no bond → refused.
            assert.ok(refused(summarizeBondAngles(bondAnglePayload([51, 0, 49]))), 'a window that misses 51 % was not refused');
        }
    },
    {
        file: SHELLS,
        gate: 'a distortion\'s z (step 7)',
        find: 'A distortion needs `analyze_clusters`\' effect size at',
        number: /effect size at z ≥ (?<n>\d+(?:\.\d+)?)/,
        constants: [CLUSTER_SIGMA]
    },
    {
        file: SHELLS,
        gate: 'a distortion\'s z (checklist)',
        find: 'A distortion claim rests on the cluster effect size at',
        number: /effect size at z ≥ (?<n>\d+(?:\.\d+)?)/,
        constants: [CLUSTER_SIGMA]
    },
    {
        file: SHELLS,
        gate: 'the cluster elements: the fewest neighbours (step 7)',
        find: 'own first shell is a finite cluster\'s (',
        number: /a finite cluster's \((?<n>\d+(?:\.\d+)?)–\d+(?:\.\d+)? neighbours/,
        constants: [CN_MATCH]
    },
    {
        file: SHELLS,
        gate: 'the cluster elements: the most neighbours (step 7)',
        find: 'own first shell is a finite cluster\'s (',
        number: /a finite cluster's \(\d+(?:\.\d+)?–(?<n>\d+(?:\.\d+)?) neighbours/,
        constants: [constant('CLUSTER_CN_MAX', CLUSTER_CN_MAX, 'llm/context/shells.js, clusterElementsOf')],
        verify: () => verifyClusterRule()
    },
    {
        file: SHELLS,
        gate: 'the split-shell join\'s reach (step 2)',
        find: 'then joins a later peak within',
        number: /joins a later peak within (?<n>\d+(?:\.\d+)?) Å of the first/,
        constants: [SPLIT_REACH],
        verify: () => {
            // A second peak 0.45 Å out joins the shell (CN 6); 0.55 Å out it is the next shell's.
            const near = cuF({ r2: 1.93 + SPLIT_SHELL_REACH_A - 0.05 });
            assert.equal(near.joined, true, 'a second peak 0.45 Å from the first was not joined');
            assert.ok(cnMatches(near.cn, 6), `the joined shell holds ${near.cn} F per Cu, not 6`);
            assert.equal(cuF({ r2: 1.93 + SPLIT_SHELL_REACH_A + 0.05 }).joined, false, 'a second peak 0.55 Å from the first was joined');
        }
    },
    {
        file: SHELLS,
        gate: 'the split-shell join\'s ideal counts (step 2)',
        find: 'whose summed CN is an ideal polyhedral count',
        number: /an ideal polyhedral count \((?<n>\d+(?:, \d+)*)\)/,
        constants: [IDEAL_COUNTS],
        verify: () => {
            // 4 + 2 makes an octahedron: joined; 4 + 1 makes no ideal count either way: not joined.
            assert.equal(cuF().joined, true, '4 + 2 F per Cu was not joined');
            assert.equal(cuF({ second: 1 }).joined, false, '4 + 1 F per Cu (no ideal count) was joined');
            assert.deepEqual(IDEAL_COORDINATION.filter((count) => isIdealCoordination(count)), [...IDEAL_COORDINATION]);
        }
    },
    {
        file: LOCAL_ORDER,
        gate: 'a split site',
        find: 'Note the split sites (two positions, by the app\'s rule',
        number: /whose excess kurtosis is ≤ (?<n>[−–-]?\d+(?:\.\d+)?) and/,
        constants: [SPLIT_KAPPA],
        statistic: PC1_RULE
    },
    {
        file: LOCAL_ORDER,
        gate: 'a split site\'s sampling errors',
        find: 'Note the split sites (two positions, by the app\'s rule',
        number: /and at least (?<n>\d+) sampling errors below zero/,
        constants: [SPLIT_ERRORS]
    },
    {
        file: LOCAL_ORDER,
        gate: 'the cluster elements: the fewest neighbours',
        find: 'own first shell is closed and holds',
        number: /closed and holds (?<n>\d+(?:\.\d+)?)–\d+(?:\.\d+)? neighbours/,
        constants: [CN_MATCH]
    },
    {
        file: LOCAL_ORDER,
        gate: 'the cluster elements: the most neighbours',
        find: 'own first shell is closed and holds',
        number: /closed and holds \d+(?:\.\d+)?–(?<n>\d+(?:\.\d+)?) neighbours/,
        constants: [constant('CLUSTER_CN_MAX', CLUSTER_CN_MAX, 'llm/context/shells.js, clusterElementsOf')]
    },
    {
        file: LOCAL_ORDER,
        gate: 'the lower group\'s margin above the noise floor',
        find: 'lowest rung clearly above the noise floor (from',
        number: /\(from (?<n>\d+(?:\.\d+)?) × `noise_floor_A`\)/,
        constants: [FLOOR_MARGIN]
    },
    {
        file: SYMMETRY,
        gate: 'the lower group\'s margin above the noise floor (step 4)',
        find: 'A rung below `noise_floor_A` is within the noise of the site means. A',
        number: /rung at or above (?<n>\d+(?:\.\d+)?) × `noise_floor_A` \(the app's margin\)/,
        constants: [FLOOR_MARGIN]
    },
    {
        file: SYMMETRY,
        gate: 'the limit of the sampling (step 4)',
        find: 'A rung between the floor and',
        number: /between the floor and (?<n>\d+(?:\.\d+)?) × it/,
        constants: [FLOOR_MARGIN],
        verify: () => verifyToleranceRule()
    },
    {
        file: SYMMETRY,
        gate: 'a lower-symmetry CIF (step 7)',
        find: '(the Detected SG card\'s export) only for a rung at or above',
        number: /only for a rung at or above (?<n>\d+(?:\.\d+)?) ×/,
        constants: [FLOOR_MARGIN]
    },
    {
        file: LOCAL_ORDER,
        gate: 'an empty directions analysis',
        find: 'Empty means no peak at',
        number: /no peak at ≥ (?<n>\d+(?:\.\d+)?)σ/,
        constants: [SIGMA],
        statistic: PEAK_TEST
    },
    {
        file: LOCAL_ORDER,
        gate: 'a long wave',
        find: 'A long wave (every |q_i|',
        number: /every \|q_i\| ≤ (?<n>\d+(?:\.\d+)?)\)/,
        constants: [LONG_WAVE],
        verify: () => {
            assert.equal(isLongWave([LONG_WAVE_Q, -LONG_WAVE_Q, 0]), true, '|q_i| = 0.25 is no long wave');
            assert.equal(isLongWave([0.3, 0, 0]), false, '|q_i| = 0.3 is a long wave');
            assert.equal(isLongWave([0.5, 0, 0]), false, 'a zone-boundary point is a long wave');
        }
    },
    {
        file: LOCAL_ORDER,
        gate: 'the reported order',
        find: 'the order found, or "none at',
        number: /"none at ≥ (?<n>\d+(?:\.\d+)?)σ"/,
        constants: [SIGMA, ORDER_SIGMA, CLUSTER_SIGMA]
    },
    {
        file: ORDER_SIGNATURES,
        gate: 'an ordered zone-boundary point: the ratio',
        find: '`analyze_modulations`: the zone-boundary point,',
        number: /the zone-boundary point, ratio ≥ (?<n>\d+(?:\.\d+)?) and/,
        constants: [ORDER_RATIO]
    },
    {
        file: ORDER_SIGNATURES,
        gate: 'an ordered zone-boundary point: the significance',
        find: '`analyze_modulations`: the zone-boundary point,',
        number: /the zone-boundary point, ratio ≥ \d+(?:\.\d+)? and ≥ (?<n>\d+(?:\.\d+)?)σ/,
        constants: [ORDER_SIGMA]
    },
    {
        file: ORDER_SIGNATURES,
        gate: 'order at a general wavevector: the ratio',
        find: 'neither zone boundary nor long wave,',
        number: /nor long wave, ratio ≥ (?<n>\d+(?:\.\d+)?) and/,
        constants: [ORDER_RATIO]
    },
    {
        file: ORDER_SIGNATURES,
        gate: 'order at a general wavevector: the significance',
        find: 'neither zone boundary nor long wave,',
        number: /nor long wave, ratio ≥ \d+(?:\.\d+)? and ≥ (?<n>\d+(?:\.\d+)?)σ/,
        constants: [ORDER_SIGMA]
    },
    {
        file: SYMMETRY,
        gate: 'the lower group\'s margin above the noise floor',
        find: '"The" group is the lowest rung at or above',
        number: /lowest rung at or above (?<n>\d+(?:\.\d+)?) × `noise_floor_A`/,
        constants: [FLOOR_MARGIN],
        verify: () => verifyToleranceRule()
    },
    {
        file: FIT,
        gate: 'a long wave (step 7)',
        find: 'A long wave has every |q_i|',
        number: /every \|q_i\| ≤ (?<n>\d+(?:\.\d+)?) of the cell/,
        constants: [LONG_WAVE]
    },
    {
        file: FIT,
        gate: 'the order rule\'s ratio',
        find: 'The test, as `analyze_modulations` applies it:',
        number: /it: ratio ≥ (?<n>\d+(?:\.\d+)?) and/,
        constants: [ORDER_RATIO],
        verify: () => assert.deepEqual([MODULATION_RULES.ratio, MODULATION_RULES.sigma], [MODULATION_RATIO, MODULATION_SIGMA])
    },
    {
        file: FIT,
        gate: 'the order rule\'s significance',
        find: 'The test, as `analyze_modulations` applies it:',
        number: /and ≥ (?<n>\d+(?:\.\d+)?)σ against independent cells/,
        constants: [ORDER_SIGMA]
    },
    {
        file: FIT,
        gate: 'a long wave (checklist)',
        find: 'No order is claimed at a box-scale wavevector',
        number: /every \|q_i\| ≤ (?<n>\d+(?:\.\d+)?)\)/,
        constants: [LONG_WAVE]
    },
    {
        file: FIT,
        gate: 'a long wave (firm rule)',
        find: 'an answer that claims order at a box-scale wavevector',
        number: /every \|q_i\| ≤ (?<n>\d+(?:\.\d+)?)\)/,
        constants: [LONG_WAVE]
    },
    {
        // The user's line (interview 2026-10-24): 5 %, the shell module's own constant since the
        // fit-assessment code PR (the shell separation keeps SEPARATION_G_FRACTION, 2 %).
        file: FIT,
        gate: 'a constraint-limited peak (step 4)',
        find: 'When g(r) just above the cut is above',
        number: /is above (?<n>\d+(?:\.\d+)?) % of the peak/,
        unit: '%',
        constants: [constant('CONSTRAINT_LIMIT_FRACTION', CONSTRAINT_LIMIT_FRACTION, 'llm/context/shells.js, the constraint-limited test (find_gr_peaks\' warning and g_at_constraint_fraction, the context\'s cut_at_min_distance)')],
        verify: () => verifyConstraintLimit()
    }
];

// --- NOT_GATES ------------------------------------------------------------------------

// A number after a comparator that no module decides with: `find` holds it.
export const NOT_GATES = [
    { file: ANISOTROPY, find: '> 0 heavy-tailed or peaked', why: 'the sign convention of Mardia\'s kurtosis (0 for a Gaussian); the app decides nothing at it' },
    { file: ANISOTROPY, find: '< 0 flat-topped', why: 'the sign convention of Mardia\'s kurtosis; a split site is the app\'s PC1 rule, bound below' },
    {
        file: ANISOTROPY,
        find: 'the same sign in at least 3 independent runs (the',
        why: 'the user\'s bar for a cluster distortion (interview 2026-10-24, answer 5; the coordination-shells skill\'s bar): every tool reads one run, no module decides it'
    },
    { file: ANISOTROPY, find: 'the same sign in at least 3 independent runs).', why: 'the user\'s bar for a cluster distortion (interview 2026-10-24, answer 5), in the checklist: no module decides it' },
    {
        file: ANISOTROPY,
        find: 'finds the order at (0, 0, ½) at > 8σ',
        why: 'measured on the zone-boundary dev case, written as analyze_modulations writes a σ above 8 ("> 8σ"): a value, not a threshold'
    },
    { file: SHELLS, find: 'a broad shoulder below 180°', why: 'geometry (the trans angle), not a threshold' },
    { file: SHELLS, find: 'the same sign in at least 3 independent runs, and', why: 'the user\'s bar for a distortion (interview 2026-10-24): every tool reads one run, no module decides it' },
    { file: SHELLS, find: 'the same sign in at least 3 independent runs and', why: 'the user\'s bar for a distortion (interview 2026-10-24), in the checklist: no module decides it' },
    { file: ORDER_SIGNATURES, find: '9.1× isotropic', why: 'the two-site dev case\'s measured peak enhancement, not a threshold' },
    { file: ORDER_SIGNATURES, find: '498× the median, > 8σ', why: 'measured on the ordered-antiphase dev case (a σ above 8 as analyze_modulations writes it): not a threshold' },
    { file: ORDER_SIGNATURES, find: 'every kurtosis within ±0.15', why: 'the range of the thermal-only dev case\'s measured kurtoses, not a gate' },
    {
        file: FIT,
        find: 'at least 3 independent runs or configurations and at least 10 accepted moves per atom',
        why: 'the user\'s sampling bar (interview 2026-10-24): no module decides the 3 runs or the 10 moves; the answer check "single-configuration" flags a distribution quoted without saying it is one configuration\'s'
    },
    { file: FIT, find: 'a partial under 5 % of the signal', why: 'the user\'s line for a weakly weighted partial (interview 2026-10-24): the context gives the weights, no module decides at 5 %' },
    {
        file: LOCAL_ORDER,
        find: 'A partial under 5 % of the signal is',
        why: 'the user\'s line for a weakly weighted partial (interview 2026-10-24, answer 3; rmc-fit-assessment\'s too): the context gives the weights, no module decides at 5 %'
    },
    { file: 'rmc-fit-assessment/run-numbers.md', find: 'with a doubled cell below 29 K', why: 'a published transition temperature' }
];

// --- FIRM_RULES -----------------------------------------------------------------------

// `find` is a phrase of the firm-rule sentence; `enforcedBy` names the code;
// `probe` exercises it (it throws when the code no longer enforces the rule).
export const FIRM_RULES = [
    {
        file: ANISOTROPY,
        find: 'claims a preferred displacement direction when every directions analysis of the question found none',
        enforcedBy: 'the answer check "direction" (the verification\'s FACTS), from analyze_displacement_directions\' facts.direction',
        probe: () => {
            const { facts } = summarizeOrientation(orientationPayload({ peak: 1.2, map: 0.5 }), null);
            findsRule('Site 1 has a preferred displacement direction along [1 1 1].', { facts: [facts] }, 'direction');
            const significant = summarizeOrientation(orientationPayload({ peak: 4, map: 0.5 }), null).facts;
            assert.deepEqual(judgeAnswer('Site 1 has a preferred displacement direction along [1 1 1].', { facts: [significant] }).audit, []);
            lowersOutcome('Site 1 has a preferred displacement direction.', { facts: [facts] });
        }
    },
    {
        file: ANISOTROPY,
        find: 'only the peak test makes a direction significant; the map test is reported as anisotropy',
        enforcedBy: 'analyze_displacement_directions\' facts.direction (summarizeOrientation: significant on the peak test, anisotropic on the map test) and its checks',
        probe: () => {
            // The thermal-only O of the skill: map 7.2σ, peak 2.3σ.
            const thermal = summarizeOrientation(orientationPayload({ peak: 2.3, map: 7.2 }), null);
            assert.equal(thermal.facts.direction.significant, false, 'a map test alone made a direction significant');
            assert.equal(thermal.facts.direction.anisotropic, true, 'the map test is not reported as anisotropy');
            assert.ok(thermal.checks.some((check) => /anisotropic/.test(check.text) && /not a preferred direction/.test(check.text)),
                'the checks do not call a map test alone anisotropy, not a preferred direction');
            assert.ok(!thermal.checks.some((check) => /^Preferred direction/.test(check.text)), 'a map test alone named a preferred direction');
            assert.equal(summarizeOrientation(orientationPayload({ peak: SIGNIFICANCE_SIGMA, map: 0.5 }), null).facts.direction.significant, true,
                'the peak test at 3σ did not make a direction significant');
        }
    },
    {
        file: SHELLS,
        find: '`find_gr_peaks` joins the split shell (CN 6) and offers the first peak alone as `first_peak_only`',
        enforcedBy: 'find_gr_peaks\' join (shells.js readShell: SPLIT_SHELL_REACH_A, IDEAL_COORDINATION), the Workbench\'s window question '
            + '(clarify.js beforeCall, on the CN margin) and the answer check "coordination" (the verification\'s FACTS)',
        probe: async () => {
            const run = splitShellRun();
            const tool = buildRunTools(run).find((item) => item.name === 'find_gr_peaks');
            assert.ok(tool, 'no find_gr_peaks for a run with partials');
            const { data, facts } = await tool.run({ pair: 'Cu-F' });
            assert.ok(cnMatches(data.first_shell?.cn?.F_per_Cu, 6), 'find_gr_peaks did not join the CuF6 shell (CN 6)');
            const alone = data.first_shell.first_peak_only;
            assert.ok(alone && cnMatches(alone.cn?.F_per_Cu, 4), 'find_gr_peaks does not offer the first peak alone (CN 4)');
            // The Workbench asks before an analysis runs in the first peak alone; the whole shell runs unasked.
            const clarifier = buildClarifier(run);
            assert.equal(clarifier.beforeCall(windowedAnalysis, { lo: alone.window_A[0], hi: alone.window_A[1] })?.kind, 'cut',
                'the Workbench did not ask before an analysis in the first peak alone');
            assert.equal(clarifier.beforeCall(windowedAnalysis, { lo: data.first_shell.window_A[0], hi: data.first_shell.window_A[1] }), null,
                'the Workbench asked about the whole shell');
            // An answer that makes Cu four-coordinate gets a warning under it.
            findsRule('Cu is four-coordinate.', { facts: [facts] }, 'coordination');
            lowersOutcome('Cu is four-coordinate.', { facts: [facts] });
        }
    },
    {
        file: SHELLS,
        find: '`compute_bond_angles` refuses, unless `window_confirmed` is set, a window that holds fewer partner atoms than the whole first shell',
        enforcedBy: 'compute_bond_angles\' refusal (summarizeBondAngles, by the one cut rule, shells.js windowCut: the CN margin, cnMarginCut, where '
            + 'a window meets its shell; the coordination counts, countsCut, where there is no shell to compare) and its window_confirmed parameter',
        probe: async ({ tools, bondAngles, catalogue } = {}) => {
            // A window that cuts its shell (4 of 6), its counts uniform: refused,
            // with the counts and the shell reading that show the cut.
            const cut = bondAnglePayload([0, 0, 0, 0, 100]);
            const result = summarizeBondAngles(cut, { shells: [shellReading(4)] });
            assert.ok(refused(result), 'a window holding 4 of 6 returned its angles');
            assert.deepEqual(result.data.coordination.counts, cut.coordination, 'the refusal does not return the counts');
            assert.equal(result.data.shell_ab?.cn_in_window, 4, 'the refusal does not return the shell reading that shows the cut');
            assert.ok(result.checks.some((check) => check.level === 'fail' && /window_confirmed/.test(check.text)));
            assert.ok(!refused(summarizeBondAngles(cut, { shells: [shellReading(4)], windowConfirmed: true })), 'window_confirmed did not keep the user\'s window');
            // With no shell to compare, the counts refuse (the fallback).
            assert.ok(refused(summarizeBondAngles(bondAnglePayload([0, 2, 30, 40, 20, 8, 0]))), 'with no shell, a 40 % modal share returned its angles');
            const tool = tools?.get('compute_bond_angles');
            if (tools) assert.ok(tool?.parameters?.properties?.window_confirmed, 'compute_bond_angles has no window_confirmed parameter');
            if (!tool || !bondAngles || !catalogue) return;
            // On the Demo, as the app runs it: the Ta–Se first peak (3.25 of 6) and the spot check's
            // 2.3–2.7 Å window are cut on the CN margin and refused; the user's window is kept; the
            // whole shell runs.
            assert.equal(typeof (await bondAngles({ ...SE_TA_SE, r12_min: 2.21, r12_max: 2.65 })).refused, 'string', 'the Demo\'s Ta–Se first peak returned its angles');
            assert.equal(catalogue.window('Ta', 'Se', [2.3, 2.7])?.cnMargin.cut, true, 'the 2.3–2.7 Å Ta–Se window is not cut on the CN margin');
            assert.equal(typeof (await bondAngles({ ...SE_TA_SE, r12_min: 2.3, r12_max: 2.7 })).refused, 'string', 'the Demo\'s 2.3–2.7 Å window returned its angles');
            const kept = (await tool.run({ ...SE_TA_SE, r12_min: 2.3, r12_max: 2.7, window_confirmed: true })).data;
            assert.ok(!kept.refused && 'mean_deg' in kept, 'window_confirmed did not keep the user\'s window on the Demo');
            const [lo, hi] = catalogue.shell('Ta', 'Se').window;
            assert.ok(!(await bondAngles({ ...SE_TA_SE, r12_min: lo, r12_max: hi })).refused, 'the whole Ta–Se shell was refused');
        }
    },
    {
        file: SYMMETRY,
        find: '`detect_symmetry` and the context\'s ladder state the point groups\' facts',
        enforcedBy: 'detect_symmetry\'s point_group_facts (summarizeSymmetry), the context\'s ladder (buildRunContext), and the answer checks "chirality", "centrosymmetry", "polar", "inversion-lost" (the verification\'s FACTS)',
        probe: () => {
            const detected = summarizeSymmetry({ spaceGroup: 'P-42_1m', spaceGroupNumber: 113, pointGroup: '-42m', nSpace: 8, orbits: [] }, 0.025);
            assert.ok(detected.data.point_group_facts, 'detect_symmetry states no point-group facts');
            assert.equal(detected.facts.groups[0].sohncke, false);
            const context = buildRunContext({ symmetry: demoSymmetry });
            assert.ok(context.symmetry.ladder.every((rung) => rung.pg), 'the context\'s ladder has no point groups');
            assert.ok(context.symmetry.point_group_facts?.length, 'the context states no point-group facts');
            const audit = { facts: [detected.facts], context };
            findsRule('P-42₁m is chiral.', audit, 'chirality');
            findsRule('F-43m is centrosymmetric.', audit, 'centrosymmetry');
            findsRule('P-42₁m is polar.', audit, 'polar');
            findsRule('Going from F-43m to P-42₁m, the inversion centre is lost.', audit, 'inversion-lost');
            lowersOutcome('P-42₁m is chiral.', audit);
        }
    },
    {
        file: FIT,
        find: 'says the run has converged while the watchdog or `assess_convergence` says improving or diverging',
        enforcedBy: 'the answer check "convergence" (the verification\'s FACTS), from assess_convergence\'s facts and the context\'s watchdog_status',
        probe: async () => {
            const tool = buildRunTools({ rValueFile: fallingHistory }).find((item) => item.name === 'assess_convergence');
            assert.ok(tool, 'no assess_convergence for a χ² history');
            const { facts } = await tool.run({});
            assert.equal(facts.convergence, 'improving');
            findsRule('The run has converged.', { facts: [facts] }, 'convergence');
            const context = buildRunContext({ rValueFile: fallingHistory });
            assert.equal(context.convergence?.watchdog_status, 'improving', 'the context\'s watchdog_status is not improving');
            findsRule('The run has converged.', { context }, 'convergence');
            lowersOutcome('The run has converged.', { context });
        }
    },
    {
        file: FIT,
        find: 'an answer that claims order at a box-scale wavevector',
        enforcedBy: 'the answer check "long-wave" (the verification\'s FACTS, the modulation row), from analyze_modulations\' facts.modulation.long_wave_q (the engine\'s LONG_WAVE_Q, isLongWave)',
        probe: async ({ tools } = {}) => {
            // The facts as analyze_modulations returns them: the engine's bound, a long wave among the strongest modes.
            const q = [0.1, 0.2, -0.1];
            assert.equal(isLongWave(q), true);
            const facts = { modulation: { elements: ['Se'], ordered: [], general: [], long_wave_q: MODULATION_RULES.longWaveQ, long_waves: [{ element: 'Se', q }] } };
            findsRule('Se is ordered at q = (0.1, 0.2, -0.1).', { facts: [facts] }, 'long-wave');
            findsRule('The Se displacements show a long-wavelength modulation.', { facts: [facts] }, 'long-wave');
            lowersOutcome('Se is ordered at q = (0.1, 0.2, -0.1).', { facts: [facts] });
            // The user's own wording, and a point that passed its own test at another q, raise no finding.
            const said = judgeAnswer('Se\'s strongest mode, q = (0.1, 0.2, -0.1), is a long wave: a correlation the finite box can\'t represent, so it says nothing about order at its own q.', { facts: [facts] });
            assert.deepEqual(said.audit, []);
            const ordered = { modulation: { ...facts.modulation, ordered: [{ element: 'Se', q: [0, 0, 0.5] }] } };
            assert.deepEqual(judgeAnswer('Se is ordered at q = (0, 0, ½).', { facts: [ordered] }).audit, []);
            // On the Demo: the tool's facts carry the engine's bound.
            const modulations = tools?.get('analyze_modulations');
            if (modulations) {
                const result = await modulations.run({});
                assert.equal(result.facts.modulation.long_wave_q, LONG_WAVE_Q, 'analyze_modulations\' facts do not carry the long-wave bound');
                assert.ok(result.facts.modulation.long_waves.every((wave) => isLongWave(wave.q)));
            }
        }
    },
    {
        file: FIT,
        find: 'an answer that quotes a distribution from a single configuration without saying so',
        enforcedBy: 'the answer check "single-configuration" (the verification\'s FACTS, the distribution row), from compute_bond_angles\', analyze_displacement_directions\' and analyze_clusters\' facts.distribution',
        probe: async ({ tools } = {}) => {
            // A bond-angle distribution the tool returned (the counts uniform: not refused).
            const { facts } = summarizeBondAngles(bondAnglePayload([0, 0, 0, 0, 0, 0, 1000]));
            assert.equal(facts.distribution?.configurations, 1, 'compute_bond_angles returns no distribution fact');
            findsRule('Se–Ta–Se angles: mean 100°, sd 20°.', { facts: [facts] }, 'single-configuration');
            lowersOutcome('Se–Ta–Se angles: mean 100°, sd 20°.', { facts: [facts] });
            assert.deepEqual(judgeAnswer('In this one configuration the Se–Ta–Se angles have mean 100°, sd 20°.', { facts: [facts] }).audit, []);
            const directions = summarizeOrientation(orientationPayload({ peak: 1.2, map: 0.5 }), null).facts;
            assert.equal(directions.distribution?.kind, 'directions', 'analyze_displacement_directions returns no distribution fact');
            // On the Demo: the whole GaSe4 shell's angles.
            const angles = tools?.get('compute_bond_angles');
            if (angles) {
                const result = await angles.run({ end1: 'Se', apex: 'Ga', end2: 'Se', r12_min: 2.18, r12_max: 2.79 });
                findsRule('Se–Ga–Se angles: mean 109.4°, sd 4.1°.', { facts: [result.facts] }, 'single-configuration');
            }
        }
    },
    {
        file: FIT,
        find: 'the Outcome is computed from the checks; the model\'s line only explains it',
        enforcedBy: 'the computed verdict (verification/verdict.js computeOutcome; the model\'s Outcome line is read by outcome.js as `claimed` only)',
        probe: () => {
            // A check contradicts the answer: not verified, whatever grade the model's line gives.
            const { facts } = summarizeOrientation(orientationPayload({ peak: 1.2, map: 0.5 }), null);
            const contradicted = judgeAnswer('Site 1 has a preferred displacement direction.\n\nOutcome: achieved — the direction is clear.', { facts: [facts] });
            assert.equal(contradicted.outcome.status, 'not verified', 'the model\'s "achieved" set the Outcome');
            assert.equal(contradicted.outcome.claimed, 'achieved');
            assert.equal(contradicted.outcome.computed, true);
            // Nothing contradicts it: the computed Outcome, which the model's line may lower but never raise.
            assert.equal(judgeAnswer('The run is still improving.\n\nOutcome: partial — more steps would help.').outcome.status, 'partial');
        }
    },
    {
        file: LOCAL_ORDER,
        find: 'an answer that claims a superstructure, a cluster distortion, breathing, rigidity or a preferred direction',
        enforcedBy: 'the answer checks "superstructure", "cluster" and "direction" (the verification\'s FACTS) and the evidence gate\'s "no order" row '
            + '(claims.js: analyze_modulations, analyze_displacement_directions with site_ref, analyze_clusters for each element clusterElementsOf names; THY-18)',
        probe: () => {
            findsRule('This structure has a superstructure: a cell doubling.', { facts: [{ modulation: { elements: ['Ti'], ordered: [], general: [] } }] }, 'superstructure');
            const cluster = {
                cluster: { element: 'Mo', members: 3, window: [2.5, 2.9], distorted: false, static: false, breathing: false, rigid: false, split: false, families: [] }
            };
            findsRule('The Mo3 clusters are distorted.', { facts: [cluster] }, 'cluster');
            findsRule('The Mo3 clusters breathe.', { facts: [cluster] }, 'cluster');
            findsRule('The Mo3 clusters are rigid.', { facts: [cluster] }, 'cluster');
            findsRule('Site 1 has a preferred displacement direction.', { facts: [summarizeOrientation(orientationPayload({ peak: 1.2, map: 0.5 }), null).facts] }, 'direction');
            // "No order" before its analyses ran is sent back to run them.
            const catalog = { elements: ['Mo', 'O'], shell: (a, b) => (a === b ? { cn: a === 'Mo' ? 2 : 8 } : null) };
            const judged = createVerification({ catalog }).judge('There is no local order beyond thermal motion.', { canSendBack: true });
            assert.ok(judged?.sendBack, 'a "no order" answer before its analyses was not sent back');
            assert.equal(judged.sendBack.id, 'no-order');
            assert.deepEqual(judged.sendBack.missing.map((item) => item.split(' ')[0]), ['analyze_modulations', 'analyze_displacement_directions', 'analyze_clusters']);
            assert.ok(judged.sendBack.missing[2].startsWith('analyze_clusters for Mo'), 'the gate does not ask for the finite-cluster element');
        }
    }
];

// --- VALUES ---------------------------------------------------------------------------

// A value a skill quotes as measured on the Demo, bound to what the tool shows
// there now: `find` locates the sentence, `number` reads the value as written
// (group n, and m for a window), `measure(probeContext)` gives the tool's
// value (a number, or [lo, hi]) and `tool` names where it comes from. The
// probe context (skillLint.test.js) runs the tools on the Demo: gr(pair) is
// find_gr_peaks' data, bondAngles(args) compute_bond_angles' data,
// triplets(params) the bond-angle engine's payload, catalogue the run's shell
// catalogue (shells.js). Bound so far: the numbers the shell reading gives
// (architecture candidate #6) and the counts summary beside them.
const SE_TA_SE = { end1: 'Se', apex: 'Ta', end2: 'Se' };
const firstShell = async (ctx, pair) => (await ctx.gr(pair)).first_shell;
export const VALUES = [
    {
        file: SHELLS,
        what: 'the two Ta–Se peaks of the TaSe₆ shell',
        find: 'the TaSe₆ shell has Ta–Se peaks at',
        number: /peaks at (?<n>\d+\.\d+) and (?<m>\d+\.\d+) Å/,
        tool: 'find_gr_peaks first_shell.peaks_in_shell_A',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Se')).peaks_in_shell_A
    },
    {
        file: SHELLS,
        what: 'the Ta–Se first-peak window',
        find: 'the first-peak window',
        number: /first-peak window (?<n>\d+\.\d+)–(?<m>\d+\.\d+) Å/,
        tool: 'find_gr_peaks first_shell.first_peak_only.window_A',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Se')).first_peak_only?.window_A ?? null
    },
    {
        file: SHELLS,
        what: 'the Ta–Se first peak\'s CN',
        find: 'Se per Ta, the whole shell',
        before: 24,
        number: /holds (?<n>\d+(?:\.\d+)?) Se per Ta, the whole shell/,
        tool: 'find_gr_peaks first_shell.first_peak_only.cn.Se_per_Ta',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Se')).first_peak_only?.cn?.Se_per_Ta ?? null
    },
    {
        file: SHELLS,
        what: 'the whole TaSe₆ shell\'s CN',
        find: 'Se per Ta, the whole shell',
        number: /the whole shell (?<n>\d+(?:\.\d+)?)\./,
        tool: 'find_gr_peaks first_shell.cn.Se_per_Ta',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Se')).cn.Se_per_Ta
    },
    {
        file: SHELLS,
        what: 'the mean coordination in the spot check\'s 2.3–2.7 Å window',
        find: 'gemma4 used 2.3–2.7 Å, got CN',
        number: /got CN (?<n>\d+\.\d+)/,
        tool: 'compute_bond_angles coordination.mean (Se–Ta–Se, 2.3–2.7 Å)',
        measure: async (ctx) => (await ctx.bondAngles({ ...SE_TA_SE, r12_min: 2.3, r12_max: 2.7 })).coordination.mean
    },
    {
        file: SHELLS,
        what: 'the mean coordination in 2.33–2.71 Å',
        find: 'qwen3.8 started at 2.33–2.71 Å (CN',
        number: /\(CN (?<n>\d+\.\d+)\)/,
        tool: 'compute_bond_angles coordination.mean (Se–Ta–Se, 2.33–2.71 Å)',
        measure: async (ctx) => (await ctx.bondAngles({ ...SE_TA_SE, r12_min: 2.33, r12_max: 2.71 })).coordination.mean
    },
    {
        file: SHELLS,
        what: 'the Ta–Ta g(r) peak',
        find: '"Peak shifted to',
        number: /Peak shifted to (?<n>\d+\.\d+) Å/,
        tool: 'find_gr_peaks peaks_A[0].r (Ta-Ta)',
        measure: async (ctx) => (await ctx.gr('Ta-Ta')).peaks_A[0].r
    },
    {
        file: SHELLS,
        what: 'the Ta–Ta shell mean',
        find: 'sits below the shell mean',
        number: /shell mean \((?<n>\d+\.\d+) Å from r²g\)/,
        tool: 'find_gr_peaks first_shell.mean_d_A (Ta-Ta)',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Ta')).mean_d_A
    },
    {
        file: SHELL_TABLE,
        what: 'the TaSe₆ shell window',
        find: '(find_gr_peaks, 2026-10-11)',
        before: 20,
        number: /≈ (?<n>\d+\.\d+)–(?<m>\d+\.\d+) Å \(find_gr_peaks/,
        tool: 'find_gr_peaks first_shell.window_A (Ta-Se)',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Se')).window_A
    },
    {
        file: SHELL_TABLE,
        what: 'the two peaks of the TaSe₆ shell',
        find: 'Å in the review |',
        number: /\| (?<n>\d+\.\d+), (?<m>\d+\.\d+) Å \|/,
        tool: 'find_gr_peaks first_shell.peaks_in_shell_A (Ta-Se)',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Se')).peaks_in_shell_A
    },
    {
        file: SHELL_TABLE,
        what: 'the TaSe₆ shell\'s CN',
        find: 'in 2.33–2.96 Å) |',
        before: 12,
        number: /\| (?<n>\d+) \(\d+\.\d+ in/,
        tool: 'find_gr_peaks first_shell.cn.Se_per_Ta',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Se')).cn.Se_per_Ta
    },
    {
        file: SHELL_TABLE,
        what: 'the CN in the review\'s 2.33–2.96 Å window',
        find: 'in 2.33–2.96 Å) |',
        before: 12,
        number: /\((?<n>\d+\.\d+) in 2\.33–2\.96 Å\)/,
        tool: 'the shell catalogue\'s window(Ta, Se, [2.33, 2.96]).cn (what the Workbench reads)',
        measure: (ctx) => ctx.catalogue.window('Ta', 'Se', [2.33, 2.96]).cn
    },
    {
        file: SHELL_TABLE,
        what: 'the Ta–Se first-peak window',
        find: 'first-peak window',
        number: /first-peak window (?<n>\d+\.\d+)–(?<m>\d+\.\d+) Å/,
        tool: 'find_gr_peaks first_shell.first_peak_only.window_A',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Se')).first_peak_only?.window_A ?? null
    },
    {
        file: SHELL_TABLE,
        what: 'the Ta–Se first peak\'s CN',
        find: 'Å holds only',
        number: /holds only (?<n>\d+(?:\.\d+)?):/,
        tool: 'find_gr_peaks first_shell.first_peak_only.cn.Se_per_Ta',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Se')).first_peak_only?.cn?.Se_per_Ta ?? null
    },
    {
        file: SHELL_TABLE,
        what: 'the GaSe₄ shell window',
        find: '| GaSe₄ | Ga–Se |',
        number: /Ga–Se \| (?<n>\d+\.\d+)–(?<m>\d+\.\d+) Å/,
        tool: 'find_gr_peaks first_shell.window_A (Ga-Se)',
        measure: async (ctx) => (await firstShell(ctx, 'Ga-Se')).window_A
    },
    {
        file: SHELL_TABLE,
        what: 'the GaSe₄ shell\'s CN',
        find: '| GaSe₄ | Ga–Se |',
        number: /— \| (?<n>\d+\.\d+) \| mean/,
        tool: 'find_gr_peaks first_shell.cn.Se_per_Ga',
        measure: async (ctx) => (await firstShell(ctx, 'Ga-Se')).cn.Se_per_Ga
    },
    {
        file: SHELL_TABLE,
        what: 'the GaSe₄ shell\'s mean distance',
        find: '| GaSe₄ | Ga–Se |',
        number: /mean (?<n>\d+\.\d+) Å, sd/,
        tool: 'find_gr_peaks first_shell.mean_d_A (Ga-Se)',
        measure: async (ctx) => (await firstShell(ctx, 'Ga-Se')).mean_d_A
    },
    {
        file: SHELL_TABLE,
        what: 'the GaSe₄ shell\'s spread',
        find: '| GaSe₄ | Ga–Se |',
        number: /, sd (?<n>\d+\.\d+) Å;/,
        tool: 'find_gr_peaks first_shell.sd_d_A (Ga-Se)',
        measure: async (ctx) => (await firstShell(ctx, 'Ga-Se')).sd_d_A
    },
    {
        file: SHELL_TABLE,
        what: 'the constraint that shapes the GaSe₄ shell',
        find: 'low-r side shaped by the',
        number: /shaped by the (?<n>\d+\.\d+) Å constraint/,
        tool: 'find_gr_peaks min_distance_constraint_A (Ga-Se), when g_at_constraint_fraction flags it',
        measure: async (ctx) => {
            const data = await ctx.gr('Ga-Se');
            return data.g_at_constraint_fraction != null ? data.min_distance_constraint_A : null;
        }
    },
    {
        file: SHELL_TABLE,
        what: 'the Ta₄ shell mean',
        find: '| Ta₄ | Ta–Ta |',
        number: /shell mean (?<n>\d+\.\d+) Å \(r²g\)/,
        tool: 'find_gr_peaks first_shell.mean_d_A (Ta-Ta)',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Ta')).mean_d_A
    },
    {
        file: SHELL_TABLE,
        what: 'the Ta–Ta g(r) peak',
        find: '| Ta₄ | Ta–Ta |',
        number: /g peak (?<n>\d+\.\d+) Å/,
        tool: 'find_gr_peaks peaks_A[0].r (Ta-Ta)',
        measure: async (ctx) => (await ctx.gr('Ta-Ta')).peaks_A[0].r
    },
    {
        file: SHELL_TABLE,
        what: 'the Ta₄ shell\'s CN',
        find: '| Ta₄ | Ta–Ta |',
        number: /g peak \d+\.\d+ Å \| (?<n>\d+) \|/,
        tool: 'find_gr_peaks first_shell.cn.Ta_per_Ta',
        measure: async (ctx) => (await firstShell(ctx, 'Ta-Ta')).cn.Ta_per_Ta
    },
    {
        file: SHELL_TABLE,
        what: 'the modal share in the 2.3–2.7 Å window',
        find: 'of the Ta atoms share the most common Se count',
        before: 12,
        number: /only (?<n>\d+) % of the Ta atoms/,
        unit: '%',
        tool: 'compute_bond_angles coordination.modal_share (Se–Ta–Se, 2.3–2.7 Å: the engine\'s coordinationSummary, kept as information)',
        measure: async (ctx) => (await ctx.bondAngles({ ...SE_TA_SE, r12_min: 2.3, r12_max: 2.7 })).coordination.modal_share
    },
    {
        file: SHELL_TABLE,
        what: 'the modal Se count in the 2.3–2.7 Å window',
        find: 'share the most common Se count',
        number: /most common Se count \((?<n>\d+)\)/,
        tool: 'compute_bond_angles coordination.modal (Se–Ta–Se, 2.3–2.7 Å)',
        measure: async (ctx) => (await ctx.bondAngles({ ...SE_TA_SE, r12_min: 2.3, r12_max: 2.7 })).coordination.modal
    },
    {
        file: RUN_NUMBERS,
        what: 'g(r) just above the Ga–Se constraint',
        find: 'g(r) just above it is',
        number: /just above it is (?<n>\d+(?:\.\d+)?) % of the first peak/,
        unit: '%',
        tool: 'find_gr_peaks g_at_constraint_fraction (Ga-Se; the shell catalogue\'s constraint.gFraction)',
        measure: (ctx) => ctx.catalogue.shell('Ga', 'Se').constraint.gFraction
    }
];

// --- RULES ----------------------------------------------------------------------------

// A rule a skill states in words that the app and arm A+ apply, one reading
// since architecture candidate 8 (validation/armA.js SWEEP_RULES holds the
// app's functions): `find` locates the statement (`before` / `after` widen
// the sentence read), `states` matches a sentence that states the rule as the
// app applies it, `skill` and `app` say each rule in brief (an allow-list
// entry repeats them), and `verify` checks the app and A+ apply `app`.
const SPLIT_RULE = `PC1 resolved, its excess kurtosis ≤ ${String(SPLIT_SITE_KURTOSIS).replace('-', '−')} and ≥ ${SPLIT_SITE_SE} sampling errors below 0 (isSplitSite)`;
const CLUSTER_RULE = `analyze_clusters for each element whose own first shell holds 0.5–${CLUSTER_CN_MAX} neighbours (clusterElementsOf, CLUSTER_CN_MAX)`;
const CONVERGENCE_RULE = `every χ² term of the .log and their total, none still falling beyond ${TREND_NOISE_FACTOR} × its noise (the scatter about `
    + `its least-squares line) over the last tenth (≥ ${TREND_MIN_ROWS} rows, from the last restart): assess_convergence, the watchdog and the `
    + 'context\'s watchdog_status (watchdog/heuristics.js judgeConvergence)';
const TOLERANCE_RULE = `a group counts where it holds to ${NOISE_FLOOR_MARGIN} × the noise floor (lowerGroupTolerance, NOISE_FLOOR_MARGIN)`;
// analyze_displacement_directions today: one test, one tier, one resolution.
const DIRECTION_TIERS_RULE = `a preferred direction at a peak σ ≥ ${SIGNIFICANCE_SIGMA} (unrounded) whatever its enhancement, else "no significant peak direction": no hint tier (summarizeOrientation)`;
const DIRECTION_RESOLUTION_RULE = `any one directions analysis at ≥ ${SIGNIFICANCE_SIGMA}σ backs a preferred direction, at whichever frequency it ran (facts.direction, the answer check "direction")`;
// Order between cells as the answer checks and arm A+ read it today (the
// diffuse verdict, plan item 4.4, was left as the user's call).
const BETWEEN_CELLS_RULE = 'the answer checks and arm A+ read a zone-boundary point and a general mode only, not the diffuse verdict '
    + '(facts.modulation; validation/armA.js MODULATION_FINDINGS_READ)';

const verifySplitRule = () => {
    assert.equal(SWEEP_RULES.splitSite, isSplitSite, 'arm A+ does not read the app\'s split rule');
    const site = (kappa, count = 4000, resolved = true) => ({ referenceNumber: 1, count, axisResolved: [resolved, true, true], excessKurtosis: [kappa, 0, 0] });
    assert.equal(isSplitSite(site(-0.6)), true);
    assert.equal(isSplitSite(site(-0.3)), false, 'a negative kurtosis alone made a site split');
    assert.equal(isSplitSite(site(-1.5, 4000, false)), false, 'an unresolved PC1 made a site split');
    assert.equal(isSplitSite(site(-0.93, 216)), false, 'κ₁ −0.93 at 216 copies (2.8 errors) made a site split');
};
const verifyClusterRule = () => {
    assert.equal(SWEEP_RULES.clusterElements, clusterElementsOf, 'arm A+ does not read the app\'s cluster-element rule');
    const cns = { Mo: 2, F: 4, O: 8, Fe: 6, Mg: 12, V: 1, K: 0.3 };
    const catalog = { elements: Object.keys(cns).sort(), shell: (a, b) => (a === b ? { cn: cns[a] } : null) };
    assert.deepEqual(clusterElementsOf(catalog), ['F', 'Mo', 'V']);
    // The evidence gate (the verification) asks analyze_clusters for exactly these before "no order".
    const missing = createVerification({ catalog }).judge('There is no local order beyond thermal motion.').done.audit[0].missing;
    assert.deepEqual(missing.filter((item) => item.startsWith('analyze_clusters for ')).map((item) => item.split(' ')[2]), clusterElementsOf(catalog),
        'the gate reads another rule');
};
const verifyToleranceRule = () => {
    assert.equal(SWEEP_RULES.lowerGroupTolerance, lowerGroupTolerance, 'arm A+ does not read the app\'s lower-group tolerance');
    // A floor of 3√2 · 0.1/√450 = 0.02 Å: the lower group at 0.025 Å.
    const basis = [{ dispA: 0.1, count: 450 }];
    assert.equal(lowerGroupTolerance(basis), 0.025);
    // The Workbench asks between two groups only when both hold to 1.25 × the floor.
    const ask = (to) => buildClarifier({ symmetry: { ladder: [{ spaceGroup: 'P4mm', from: 0, to }, { spaceGroup: 'Pm-3m', from: to, to: 1 }] } }, { noiseFloor: 0.02 })
        .beforeTurn('What is the space group?');
    assert.equal(ask(0.024), null, 'a group holding to 1.2 × the floor was offered');
    assert.equal(ask(0.026)?.kind, 'tolerance');
};

export const RULES = [
    {
        // The user's convergence rule (interview 2026-10-24, answer 1), the fit-assessment code PR's.
        file: FIT,
        rule: 'convergence',
        find: 'Judge every fit term and the total',
        after: 160,
        states: /every fit term and the total[\s\S]*\bstill falling beyond its noise over the\s+last tenth\b/,
        skill: 'every fit term and the total, none still falling beyond its noise over the last tenth',
        app: CONVERGENCE_RULE,
        verify: async () => {
            // The window is the last tenth (≥ 10 rows), "beyond its noise" 3 × the residual sd.
            assert.deepEqual([TREND_WINDOW_FRACTION, TREND_MIN_ROWS, TREND_NOISE_FACTOR, MIN_TOTAL_DROP], [0.1, 10, 3, 0.1]);
            // One term still falling keeps the run improving, though the charted column alone reads converged.
            const rValueFile = termStillFalling();
            assert.equal(classifyConvergence(rValueFile.plotData.series[0].y), 'converged', 'the probe\'s charted column no longer reads flat');
            const tool = buildRunTools({ rValueFile }).find((item) => item.name === 'assess_convergence');
            const { data, facts } = await tool.run({});
            assert.equal(data.status, 'improving', 'assess_convergence judges the charted column only');
            assert.deepEqual(data.terms.map((term) => [term.term, term.trend]), [['F(Q)_1', 'falling'], ['Curvature', 'flat']]);
            assert.equal(data.total.trend, 'flat');
            assert.equal(facts.convergence, 'improving');
            const context = buildRunContext({ rValueFile });
            assert.equal(context.convergence.watchdog_status, 'improving', 'the context\'s watchdog_status judges the charted column only');
            assert.equal(context.convergence.last_10pct.terms['F(Q)_1'].trend, 'falling');
            findsRule('The run has converged.', { context }, 'convergence');
        }
    },
    {
        file: ANISOTROPY,
        rule: 'a split site',
        find: 'A site is split (two wells) only by the app\'s rule',
        after: 90,
        states: /\bPC1\b/,
        skill: 'a negative Mardia kurtosis',
        app: SPLIT_RULE,
        verify: verifySplitRule
    },
    // The user's tiers and two resolutions (interview 2026-10-24, answers 3 and
    // 4): the displacement-anisotropy code PR builds them (ALLOW_LIST).
    {
        file: ANISOTROPY,
        rule: 'a direction\'s tiers',
        find: 'Three tiers for a direction',
        after: 260,
        states: /whatever (?:its|the) enhancement|no hint tier/,
        skill: 'a finding at ≥ 3σ (unrounded) and ≥ 1.2× isotropic; a hint at ≥ 3σ under 1.2×, or 2.5–3σ; nothing below 2.5σ (the user, 2026-10-24)',
        app: DIRECTION_TIERS_RULE,
        verify: () => {
            // A 3σ peak at 1.1× isotropic is a preferred direction today, and a 2.7σ peak no direction at all.
            const small = summarizeOrientation({ ...orientationPayload({ peak: SIGNIFICANCE_SIGMA, map: 0.5 }), peakEnhancement: 1.1 }, null);
            assert.equal(small.facts.direction.significant, true, 'a 3σ direction under 1.2× isotropic is no longer preferred: update this rule');
            assert.ok(small.checks.some((check) => /^Preferred direction/.test(check.text)));
            const near = summarizeOrientation(orientationPayload({ peak: 2.7, map: 0.5 }), null);
            assert.ok(near.checks.some((check) => /^No significant peak direction/.test(check.text)));
            assert.doesNotMatch(JSON.stringify(near), /\bhint\b/i, 'the tool labels a hint: update this rule');
        }
    },
    {
        file: ANISOTROPY,
        rule: 'a finding at two resolutions',
        find: 'Two resolutions for a finding.',
        after: 360,
        states: /\bone resolution\b|\bwhichever frequency\b/,
        skill: 'a finding holds at the Auto ν and at a neighbouring ν (ν ± 1); one at the Auto ν alone is a hint (the user, 2026-10-24)',
        app: DIRECTION_RESOLUTION_RULE,
        verify: () => {
            // One analysis at one frequency backs a preferred direction today, and its fact names no frequency.
            const { facts } = summarizeOrientation(orientationPayload({ peak: SIGNIFICANCE_SIGMA + 1, map: 0.5 }), null);
            assert.deepEqual(judgeAnswer('Site 1 has a preferred displacement direction along [1 1 1].', { facts: [facts] }).audit, [],
                'one resolution no longer backs a preferred direction: update this rule');
            assert.ok(!('frequency' in facts.direction), 'the direction fact names its resolution: update this rule');
        }
    },
    {
        file: SIGNIFICANCE,
        rule: 'a split site',
        find: 'A site is split (two positions) only by the app\'s rule on its PC1',
        states: /\bPC1\b/,
        skill: 'a negative Mardia kurtosis',
        app: SPLIT_RULE,
        verify: verifySplitRule
    },
    {
        file: SHELLS,
        rule: 'the cluster elements',
        find: 'For a cluster or polyhedron distortion (Jahn–Teller, dimers, trimers), run `analyze_clusters`',
        after: 160,
        states: /own first shell|finite cluster|\b4\.5\b/,
        skill: 'analyze_clusters for any cluster or polyhedron distortion (Jahn–Teller, dimers, trimers)',
        app: CLUSTER_RULE,
        verify: verifyClusterRule
    },
    {
        file: SYMMETRY,
        rule: 'the lower group\'s tolerance',
        find: 'the app\'s margin) is a distortion of the average structure',
        before: 60,
        states: /1\.25|clearly above/,
        skill: 'a rung above noise_floor_A is a distortion',
        app: TOLERANCE_RULE,
        verify: verifyToleranceRule
    },
    {
        file: SYMMETRY,
        rule: 'the lower group\'s tolerance',
        find: 'A rung between the floor and',
        after: 60,
        states: /1\.25/,
        skill: 'a rung between the floor and 1.25 × it: "at the limit of the sampling", never a lowering (the user, 2026-10-24)',
        app: TOLERANCE_RULE,
        verify: verifyToleranceRule
    },
    // The user's interview of local-order (2026-10-24, answer 6): the diffuse
    // verdict counts as order between cells; its code is a later PR (ALLOW_LIST).
    {
        file: LOCAL_ORDER,
        rule: 'order between cells',
        find: 'no diffuse short-range order: the diffuse verdict',
        after: 80,
        states: /not the diffuse|only (?:a|the) zone-boundary point and (?:a|the) general/,
        skill: 'the diffuse verdict past its calibrated bar counts as order between cells, read by arm A+ and the answer checks (the user, 2026-10-24)',
        app: BETWEEN_CELLS_RULE,
        verify: () => {
            // The modulation fact the answer checks read holds no diffuse verdict today.
            // (long_wave_q and long_waves, the long-wave check's, since the fit-assessment code PR.)
            const modulation = FACTS.find((fact) => fact.key === 'modulation');
            assert.deepEqual(Object.keys(modulation.shape), ['elements', 'ordered', 'general', 'long_wave_q', 'long_waves'],
                'the modulation fact reads the diffuse verdict: update this rule');
        }
    },
    // local-order (batch 1) states the same three rules for "no order".
    {
        file: LOCAL_ORDER,
        rule: 'a split site',
        find: 'Note the split sites (two positions, by the app\'s rule',
        after: 60,
        states: /\bPC1\b/,
        skill: 'a negative Mardia kurtosis (flat or split: two positions)',
        app: SPLIT_RULE,
        verify: verifySplitRule
    },
    {
        file: LOCAL_ORDER,
        rule: 'the cluster elements',
        find: 'for every element whose own first shell is closed and holds',
        after: 40,
        // "own first shell" alone was in the draft: the app's rule is its count.
        states: /0\.5–4\.5|\b4\.5\b/,
        skill: 'analyze_clusters for every element whose own first shell is short and closed (dimers, trimers, tetrahedra)',
        app: CLUSTER_RULE,
        verify: verifyClusterRule
    },
    {
        file: LOCAL_ORDER,
        rule: 'the lower group\'s tolerance',
        find: 'Empty means the group of the lowest rung clearly above the noise floor',
        states: /1\.25|clearly above/,
        skill: 'the group of the lowest rung above the noise floor',
        app: TOLERANCE_RULE,
        verify: verifyToleranceRule
    }
];

// --- ALLOW_LIST -----------------------------------------------------------------------

// Today's violations (docs/AGENT_PLAN.md section E). `fragment` (≤ 60
// characters) is skill text at the violation; an entry that no finding needs
// fails the test, so a fix removes its entry. Phase 5's items 5.2–5.4 empty it.
// A value entry also names the skill's line, its value and the tool's, and a
// rule entry the line, the skill's rule and the app's, and a gate entry may (the
// skill's number as written, the constant the app decides with): the test
// fails when any of the three moves. Each `why` opens with its kind. FACT: the
// skill states a number or a rule the tool decides otherwise — the user
// decided on 2026-10-19 that the skills' facts are fixed now, the interviews
// after v1, so none is left (the skill fact-sync PR cleared them). JUDGMENT: a
// method, an order or a judgment call the skill leaves to the user — it stays,
// listed for the interviews (docs/SKILL_INTERVIEWS.md); none is left since the
// user's interviews of 2026-10-24. AWAITS CODE: an interviewed skill states the
// user's rule, which a code PR still has to build (each names its PR); the
// entry leaves when that code lands.
// The user's interview of rmc-fit-assessment (2026-10-24): its four entries (the every-term
// convergence rule, the 5 % constraint and two firm rules) left with the fit-assessment code PR,
// which bound them above (RULES, BINDINGS, FIRM_RULES).
// The user's interview of coordination-shells (2026-10-24): the firm rule's code is a later PR.
const SHELLS_CODE_PR = 'awaits the coordination code PR';
// The user's interview of local-symmetry (2026-10-24): the firm rules' code is a later PR.
const SYMMETRY_CODE_PR = 'awaits the local-symmetry code PR';
// The user's interview of displacement-anisotropy (2026-10-24): the tiers, the two resolutions
// and the firm rule's code are a later PR.
const ANISOTROPY_CODE_PR = 'awaits the displacement-anisotropy code PR';
// The user's interview of local-order (2026-10-24): the extended evidence gate and the diffuse
// verdict's reading are a later answer-check PR; the extended networks' tool is its own PR.
const ORDER_CODE_PR = 'awaits the local-order answer-check code PR';
const POLYHEDRAL_PR = 'awaits the polyhedral_distortion PR';
export const ALLOW_LIST = [
    {
        check: 'tool',
        file: LOCAL_ORDER,
        fragment: 'polyhedral_distortion',
        why: `AWAITS CODE: the user's test for an extended network's polyhedra (interview 2026-10-24, answer 1b); the tool is being built and is not on main: ${POLYHEDRAL_PR}`
    },
    {
        check: 'firm-rule',
        file: LOCAL_ORDER,
        fragment: 'an answer that concludes "no order" is refused until every',
        why: `AWAITS CODE: the user's firm rule (interview 2026-10-24, answer 5); today the evidence gate asks only for analyze_modulations, one per-site directions map and analyze_clusters for each cluster element: ${ORDER_CODE_PR} (the extended gate), then mapped in FIRM_RULES`
    },
    {
        check: 'rule',
        file: LOCAL_ORDER,
        fragment: 'no diffuse short-range order: the diffuse verdict',
        line: 34,
        skill: 'the diffuse verdict past its calibrated bar counts as order between cells, read by arm A+ and the answer checks (the user, 2026-10-24)',
        app: BETWEEN_CELLS_RULE,
        why: `AWAITS CODE: the user's rule for the diffuse verdict (interview 2026-10-24, answer 6); arm A+ and the answer checks read no diffuse verdict today: ${ORDER_CODE_PR} (the diffuse verdict read as order between cells)`
    },
    {
        check: 'gate',
        file: ANISOTROPY,
        fragment: '1.2× isotropic',
        why: `AWAITS CODE: the user's enhancement line for a finding (interview 2026-10-24, answer 3); no module applies an enhancement threshold to a direction today (summarizeOrientation calls one preferred at ≥ 3σ whatever its enhancement): ${ANISOTROPY_CODE_PR} (the 1.2× constant)`
    },
    {
        check: 'gate',
        file: ANISOTROPY,
        fragment: 'below 2.5σ',
        why: `AWAITS CODE: the user's floor for a hint (interview 2026-10-24, answer 3); summarizeOrientation has no hint tier today: ${ANISOTROPY_CODE_PR} (the tier labels and their 2.5σ constant)`
    },
    {
        check: 'rule',
        file: ANISOTROPY,
        fragment: 'Three tiers for a direction',
        line: 52,
        skill: 'a finding at ≥ 3σ (unrounded) and ≥ 1.2× isotropic; a hint at ≥ 3σ under 1.2×, or 2.5–3σ; nothing below 2.5σ (the user, 2026-10-24)',
        app: DIRECTION_TIERS_RULE,
        why: `AWAITS CODE: the user's three tiers (interview 2026-10-24, answer 3); analyze_displacement_directions labels a preferred direction or none, with no hint: ${ANISOTROPY_CODE_PR} (the tier labels)`
    },
    {
        check: 'rule',
        file: ANISOTROPY,
        fragment: 'Two resolutions for a finding.',
        line: 58,
        skill: 'a finding holds at the Auto ν and at a neighbouring ν (ν ± 1); one at the Auto ν alone is a hint (the user, 2026-10-24)',
        app: DIRECTION_RESOLUTION_RULE,
        why: `AWAITS CODE: the user's neighbouring-ν rule (interview 2026-10-24, answer 4); one analysis at one frequency backs a preferred direction today: ${ANISOTROPY_CODE_PR} (a finding at two resolutions)`
    },
    {
        check: 'firm-rule',
        file: ANISOTROPY,
        fragment: 'an answer that calls a split or flat-topped site "disorder"',
        why: `AWAITS CODE: the user's firm rule (interview 2026-10-24, answer 8); today the evidence gate's split-site claim sends back a split site called random disorder until analyze_modulations ran in the question, for any element, and reads no flat-topped site: ${ANISOTROPY_CODE_PR} (its answer check), then mapped in FIRM_RULES`
    },
    {
        check: 'firm-rule',
        file: SYMMETRY,
        fragment: 'an answer that claims a group–subgroup relation that',
        why: `AWAITS CODE: the user's firm rule (interview 2026-10-24, answer 4a); no answer check reads a group–subgroup claim today: ${SYMMETRY_CODE_PR} (its answer check), then mapped in FIRM_RULES`
    },
    {
        check: 'firm-rule',
        file: SYMMETRY,
        fragment: 'an answer that claims a cell doubling the box\'s average',
        why: `AWAITS CODE: the user's firm rule (interview 2026-10-24, answer 4b); the superstructure check reads analyze_modulations, not what the box's average structure can hold: ${SYMMETRY_CODE_PR} (its answer check), then mapped in FIRM_RULES`
    },
    {
        check: 'firm-rule',
        file: SHELLS,
        fragment: 'every CN or angle in an answer carries its window in Å',
        why: `AWAITS CODE: the user's firm rule (interview 2026-10-24, answer 8); no answer check flags a CN or angle without its window today: ${SHELLS_CODE_PR} (its answer check), then mapped in FIRM_RULES`
    }
];
