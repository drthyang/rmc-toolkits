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
//   ALLOW_LIST  today's violations, each with why and what clears it.
//
// When a skill's text changes (the user's review of a draft), its entries here
// change in the same PR: a binding or a firm rule whose sentence moved is
// "stale", a violation the change fixes leaves the allow-list (the test fails
// on an entry nothing needs), and a new threshold or firm rule gets its entry.
// Files are named by their path under .claude/skills/.

import assert from 'node:assert/strict';
import {
    LOW_INDEX_CONE_DEG, SIGNIFICANCE_SIGMA, summarizeBondAngles, summarizeOrientation, summarizeSymmetry
} from '../../copilotTools';
import { buildClarifier } from '../../llm/copilot/clarify';
import { buildRunTools } from '../../llm/copilot/runTools';
import { createVerification } from '../../llm/copilot/verification';
import { buildRunContext, lowerGroupTolerance, NOISE_FLOOR_MARGIN } from '../../llm/context/runContext';
import { CLUSTER_CN_MAX, clusterElementsOf, MODAL_COORDINATION_SHARE, UNBONDED_APEX_FRACTION } from '../../llm/context/shells';
import { SWEEP_RULES } from '../../../validation/armA';
import {
    AXIS_CONE_DEG, AXIS_ENHANCEMENT, CLUSTER_EFFECT, CLUSTER_Z, SPLIT_SITE_KURTOSIS, SPLIT_SITE_SE, TRIGONALITY_EFFECT, isSplitSite
} from '../../workers/clusters';
import { ASYMMETRY_FLAG_SIGMA, MAP_TEST_MIN_PAIRS } from '../../workers/orientation';
import { AXIS_RESOLUTION_SIGMAS } from '../../workers/pcaKde';
import { coordinationSummaryOf } from '../../workers/triplets';

const ANISOTROPY = 'displacement-anisotropy/SKILL.md';
const SIGNIFICANCE = 'displacement-anisotropy/significance-reading.md';
const SHELLS = 'coordination-shells/SKILL.md';
const SHELL_TABLE = 'coordination-shells/lacunar-spinel-shells.md';
const SYMMETRY = 'local-symmetry/SKILL.md';
const FIT = 'rmc-fit-assessment/SKILL.md';
const RUN_NUMBERS = 'rmc-fit-assessment/run-numbers.md';

const constant = (name, value, module) => ({ name, value, module });
const SIGMA = constant('SIGNIFICANCE_SIGMA', SIGNIFICANCE_SIGMA, 'copilotTools.js');
const ASYMMETRY_SIGMA = constant('ASYMMETRY_FLAG_SIGMA', ASYMMETRY_FLAG_SIGMA, 'workers/orientation.js');
// compute_bond_angles' counts rule lives in the shell module (countsCut), pinned as it was.
const MODAL_SHARE = constant('MODAL_COORDINATION_SHARE', MODAL_COORDINATION_SHARE, 'llm/context/shells.js, countsCut');

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

// The answer judged by the verification (llm/copilot/verification/): the
// tools' `facts` recorded as steps, the run `context` given, no claim offered
// (only the checks on the answer speak). Its `done`: the findings and the
// Outcome they compute.
const judgeAnswer = (text, { facts = [], context = null } = {}) => {
    const verification = createVerification({ context, offered: [] });
    facts.forEach((item) => verification.record({ name: 'analysis', args: {}, status: 'ok', error: null, executed: true, facts: item }, () => ({ data: null })));
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
        find: '"Preferred direction" only with a',
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
        constants: [
            constant('LOW_INDEX_CONE_DEG', LOW_INDEX_CONE_DEG, 'copilotTools.js, measure_symmetry_lowering'),
            constant('AXIS_CONE_DEG', AXIS_CONE_DEG, 'workers/clusters.js, analyze_clusters\' axis families')
        ]
    },
    {
        file: ANISOTROPY,
        gate: 'a split site',
        find: 'flat or split, e.g. two wells',
        before: 12,
        number: /< ?(?<n>[−–-]?\d+(?:\.\d+)?) flat or split/,
        constants: [constant('SPLIT_SITE_KURTOSIS', SPLIT_SITE_KURTOSIS, 'workers/clusters.js, isSplitSite (list_sites; the cluster verdict\'s split member sites)')],
        statistic: { pattern: /\bPC1\b|per-axis/, says: `the PC1 excess kurtosis on a resolved axis, at least ${SPLIT_SITE_SE} sampling errors below zero` },
        verify: () => {
            const site = (kappa) => ({ referenceNumber: 1, count: 4000, axisResolved: [true, true, true], excessKurtosis: [kappa, 0, 0] });
            assert.equal(isSplitSite(site(SPLIT_SITE_KURTOSIS - 0.01)), true);
            assert.equal(isSplitSite(site(SPLIT_SITE_KURTOSIS + 0.01)), false, 'a kurtosis above the bound counted as split');
        }
    },
    {
        file: ANISOTROPY,
        gate: 'the direction firm rule',
        find: 'every directions analysis of the question found none at',
        number: /found none at ≥ (?<n>\d+(?:\.\d+)?)σ/,
        constants: [SIGMA]
    },
    {
        file: ANISOTROPY,
        gate: 'the reported direction',
        find: 'the strongest direction with its enhancement and σ',
        number: /"none at ≥ (?<n>\d+(?:\.\d+)?)σ"/,
        constants: [SIGMA]
    },
    {
        file: SIGNIFICANCE,
        gate: 'the significance threshold',
        find: 'counts as significant',
        before: 14,
        number: /≥ (?<n>\d+(?:\.\d+)?)σ counts as significant/,
        constants: [SIGMA, ASYMMETRY_SIGMA, constant('CLUSTER_Z', CLUSTER_Z, 'workers/clusters.js')]
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
        gate: 'a cut window (judgment call)',
        find: 'modal-coordination threshold the app uses to call a window cut',
        before: 12,
        number: /The (?<n>\d+(?:\.\d+)?) % modal-coordination/,
        unit: '%',
        constants: [MODAL_SHARE]
    },
    {
        file: SHELLS,
        gate: 'a cut window (firm rule)',
        find: 'of the central atoms share the most common bond count',
        before: 24,
        number: /fewer than (?<n>\d+(?:\.\d+)?) % of the central atoms/,
        unit: '%',
        constants: [MODAL_SHARE],
        verify: () => {
            // 59 of 100 central atoms share the modal count: refused; 61: kept.
            assert.ok(refused(summarizeBondAngles(bondAnglePayload([0, 0, 0, 0, 41, 0, 59]))), 'a 59 % modal share was not refused');
            assert.ok(!refused(summarizeBondAngles(bondAnglePayload([0, 0, 0, 0, 39, 0, 61]))), 'a 61 % modal share was refused');
        }
    },
    {
        file: SHELLS,
        gate: 'a window that misses the shell (firm rule)',
        find: 'or most have none',
        number: /or (?<n>most) have none/,
        constants: [constant('UNBONDED_APEX_FRACTION', UNBONDED_APEX_FRACTION, 'llm/context/shells.js, countsCut')],
        verify: () => {
            // "Most" is more than half: 51 of 100 central atoms with no bond → refused.
            assert.ok(refused(summarizeBondAngles(bondAnglePayload([51, 0, 49]))), 'a window that misses 51 % was not refused');
        }
    }
];

// --- NOT_GATES ------------------------------------------------------------------------

// A number after a comparator that no module decides with: `find` holds it.
export const NOT_GATES = [
    { file: ANISOTROPY, find: '> 0 heavy-tailed or peaked', why: 'the sign convention of Mardia\'s kurtosis (0 for a Gaussian); the app decides nothing at it' },
    { file: ANISOTROPY, find: 'under 1.2× isotropic', why: 'a judgment call (to confirm with the user); no module applies an enhancement threshold to a direction' },
    { file: SHELLS, find: 'a broad shoulder below 180°', why: 'geometry (the trans angle), not a threshold' },
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
        file: SHELLS,
        find: '`compute_bond_angles` refuses a window in which fewer than',
        enforcedBy: 'compute_bond_angles\' refusal (summarizeBondAngles) and its window_confirmed parameter',
        probe: ({ tools } = {}) => {
            const cut = bondAnglePayload([0, 2, 30, 40, 20, 8, 0]);
            const result = summarizeBondAngles(cut);
            assert.ok(refused(result), 'a window with a 40 % modal share returned its angles');
            assert.deepEqual(result.data.coordination.counts, cut.coordination, 'the refusal does not return the counts that show the cut');
            assert.ok(result.checks.some((check) => check.level === 'fail' && /window_confirmed/.test(check.text)));
            assert.ok(!refused(summarizeBondAngles(cut, { windowConfirmed: true })), 'window_confirmed did not keep the user\'s window');
            const tool = tools?.get('compute_bond_angles');
            if (tools) assert.ok(tool?.parameters?.properties?.window_confirmed, 'compute_bond_angles has no window_confirmed parameter');
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
        tool: 'the bond-angle engine\'s coordinationSummary.modalShare (Se–Ta–Se, 2.3–2.7 Å; compute_bond_angles\' refusal quotes it)',
        measure: async (ctx) => (await ctx.triplets({ end1: 'Se', apex: 'Ta', end2: 'Se', r12Min: 2.3, r12Max: 2.7 })).coordinationSummary.modalShare
    },
    {
        file: SHELL_TABLE,
        what: 'the modal Se count in the 2.3–2.7 Å window',
        find: 'share the most common Se count',
        number: /most common Se count \((?<n>\d+)\)/,
        tool: 'the bond-angle engine\'s coordinationSummary.modal (Se–Ta–Se, 2.3–2.7 Å)',
        measure: async (ctx) => (await ctx.triplets({ end1: 'Se', apex: 'Ta', end2: 'Se', r12Min: 2.3, r12Max: 2.7 })).coordinationSummary.modal
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
const TOLERANCE_RULE = `a group counts where it holds to ${NOISE_FLOOR_MARGIN} × the noise floor (lowerGroupTolerance, NOISE_FLOOR_MARGIN)`;

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
        file: ANISOTROPY,
        rule: 'a split site',
        find: '< 0 flat or split, e.g. two wells',
        before: 70,
        states: /\bPC1\b/,
        skill: 'a negative Mardia kurtosis',
        app: SPLIT_RULE,
        verify: verifySplitRule
    },
    {
        file: SIGNIFICANCE,
        rule: 'a split site',
        find: 'negative for a flat or split site',
        before: 80,
        after: 0,
        states: /\bPC1\b/,
        skill: 'a negative Mardia kurtosis',
        app: SPLIT_RULE,
        verify: verifySplitRule
    },
    {
        file: SHELLS,
        rule: 'the cluster elements',
        find: 'For a cluster or polyhedron distortion (Jahn–Teller, dimers, trimers), run `analyze_clusters`',
        states: /own first shell|finite cluster|\b4\.5\b/,
        skill: 'analyze_clusters for any cluster or polyhedron distortion (Jahn–Teller, dimers, trimers)',
        app: CLUSTER_RULE,
        verify: verifyClusterRule
    },
    {
        file: SYMMETRY,
        rule: 'the lower group\'s tolerance',
        find: 'A rung above it is a distortion of the average structure',
        before: 60,
        states: /1\.25|clearly above/,
        skill: 'a rung above noise_floor_A is a distortion',
        app: TOLERANCE_RULE,
        verify: verifyToleranceRule
    },
    {
        file: SYMMETRY,
        rule: 'the lower group\'s tolerance',
        find: 'Whether a rung just above the noise floor is reported as a distortion',
        after: 60,
        states: /1\.25|clearly above/,
        skill: 'left to the user: a rung just above the noise floor reported as a distortion or as "at the limit of the sampling"',
        app: TOLERANCE_RULE,
        verify: verifyToleranceRule
    }
];

// --- ALLOW_LIST -----------------------------------------------------------------------

// Today's violations (docs/AGENT_PLAN.md section E). `fragment` (≤ 60
// characters) is skill text at the violation; an entry that no finding needs
// fails the test, so a fix removes its entry. Phase 5's items 5.2–5.4 empty it.
// A value entry also names the skill's line, its value and the tool's: the
// test fails when any of the three moves.
const FIRST_PEAK_WHY = 'the first-peak window moved from r ± FWHM, which reached into the second Ta–Se peak, to the inter-peak minimum '
    + '(shells.js; approved by the manager on the user\'s delegation, 2026-10-19: validation/findings/decisions/2026-10-19-shell-catalogue.json). '
    + 'The skill text is the user\'s and is not edited: awaiting the user\'s skill review, a precondition for the freeze';
// Candidate 8: arm A+ reads the app's rules; the skills that state them otherwise wait for the user.
const APP_RULES_WHY = 'arm A+ reads the app\'s rule since candidate 8 (approved by the manager on the user\'s delegation, 2026-10-19: '
    + 'validation/findings/decisions/2026-10-19-aplus-app-rules.json), so the skill, the app and A+ must say one rule. The skill text is the '
    + 'user\'s and is not edited: awaiting the user\'s skill review (plan 5.3)';
export const ALLOW_LIST = [
    {
        check: 'value',
        file: SHELLS,
        fragment: 'the first-peak window 2.33–2.71 Å',
        line: 32,
        skill: '2.33–2.71',
        tool: '2.21–2.65',
        why: `step 2: ${FIRST_PEAK_WHY}`
    },
    {
        check: 'value',
        file: SHELLS,
        fragment: 'holds 3.9 Se per Ta, the whole shell',
        line: 33,
        skill: '3.9',
        tool: '3.25',
        why: `step 2, the first peak's CN (≈ 3, the three short Ta–Se bonds of TaSe₆): ${FIRST_PEAK_WHY}`
    },
    {
        check: 'value',
        file: SHELL_TABLE,
        fragment: 'first-peak window 2.33–2.71 Å',
        line: 17,
        skill: '2.33–2.71',
        tool: '2.21–2.65',
        why: `the TaSe₆ row of the Demo table: ${FIRST_PEAK_WHY}`
    },
    {
        check: 'value',
        file: SHELL_TABLE,
        fragment: 'Å holds only 3.9:',
        line: 17,
        skill: '3.9',
        tool: '3.25',
        why: `the TaSe₆ row of the Demo table, the first peak's CN (≈ 3): ${FIRST_PEAK_WHY}`
    },
    {
        check: 'gate',
        file: ANISOTROPY,
        fragment: 'only when the angle to it is small',
        why: 'no number; the app names a family within 15° (LOW_INDEX_CONE_DEG, AXIS_CONE_DEG). PR #71\'s draft says "under 15°" — waiting for the user (plan 5.2)'
    },
    {
        check: 'gate',
        file: ANISOTROPY,
        fragment: '"Preferred direction" only with a significance ≥ 3σ',
        why: 'any significance counts here; the tool calls a direction preferred on its peak test only since batch 1 L3. PR #71\'s draft says "a peak significance" — waiting for the user (plan 5.2)'
    },
    {
        check: 'gate',
        file: ANISOTROPY,
        fragment: '< 0 flat or split',
        why: 'the skill\'s split site is a negative Mardia kurtosis; the tools\' and arm A+\'s (since candidate 8) is a PC1 kurtosis ≤ −0.5, 3 errors below zero. Left after PR #71: plan 5.3 (the rule entry below holds the same sentence)'
    },
    {
        check: 'rule',
        file: ANISOTROPY,
        fragment: '< 0 flat or split, e.g. two wells',
        line: 27,
        skill: 'a negative Mardia kurtosis',
        app: SPLIT_RULE,
        why: `step 1 calls a negative Mardia kurtosis "flat or split"; a negative kurtosis is also a flat-topped thermal cloud, and at 100–216 copies arm A's prolate rule called 3.6–14 % of prolate Gaussian sites split where isSplitSite called none: ${APP_RULES_WHY}`
    },
    {
        check: 'rule',
        file: SIGNIFICANCE,
        fragment: 'negative for a flat or split site',
        line: 17,
        skill: 'a negative Mardia kurtosis',
        app: SPLIT_RULE,
        why: `the Mardia kurtosis paragraph reads a negative value as a flat or split site; the app calls a site split on its PC1 kurtosis only: ${APP_RULES_WHY}`
    },
    {
        check: 'rule',
        file: SHELLS,
        fragment: 'For a cluster or polyhedron distortion (Jahn–Teller, dimers',
        line: 51,
        skill: 'analyze_clusters for any cluster or polyhedron distortion (Jahn–Teller, dimers, trimers)',
        app: CLUSTER_RULE,
        why: `step 7 sends any polyhedron distortion to analyze_clusters; the app tries clusters of an element only where its own first shell is a finite cluster's (a corner-sharing octahedron's ligands, perovskite O–O at CN 8, are an extended network the engine refuses). Dimers and trimers agree: ${APP_RULES_WHY}`
    },
    {
        check: 'rule',
        file: SYMMETRY,
        fragment: 'A rung above it is a distortion of the average structure',
        line: 38,
        skill: 'a rung above noise_floor_A is a distortion',
        app: TOLERANCE_RULE,
        why: `step 4 takes the noise floor itself as the line; the app reads the lower group at 1.25 × the floor (measure_symmetry_lowering's default) and the Workbench offers a group only where it holds to 1.25 × the floor: ${APP_RULES_WHY}`
    },
    {
        check: 'rule',
        file: SYMMETRY,
        fragment: 'Whether a rung just above the noise floor is reported as a',
        line: 59,
        skill: 'left to the user: a rung just above the noise floor reported as a distortion or as "at the limit of the sampling"',
        app: TOLERANCE_RULE,
        why: `a judgment call the skill leaves to the user, which the app decides (a rung holding only to 1–1.25 × the floor is not read as a group): ${APP_RULES_WHY}`
    }
];
