// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The skill lint (docs/AGENT_PLAN.md item 5.1) on the skills the app bundles
// (copilotSkills.js: every .claude/skills/<name>/SKILL.md and the reference
// files beside it). It reads the skills and never edits them:
//
//   tool       every tool a skill names is one the Copilot can offer — the
//              set built by the tool builders on the Demo run, as the dock
//              builds it (run tools, highlight_cards, the engine tools, the
//              view tools, ask_user, read_skill);
//   gate       every gate number a skill states equals the constant of the
//              module that decides it (skillLint/tables.js BINDINGS), and every
//              number after a comparator is bound or listed as no gate;
//   firm-rule  every "*Firm rule in the app:*" is mapped to the check or
//              refusal that enforces it, and that code is exercised;
//   value      every value a skill quotes as measured on the Demo that
//              tables.js VALUES binds equals what the tool shows there now
//              (the shell readings: find_gr_peaks, the shell catalogue, the
//              bond-angle counts summary).
//
// Today's violations are in ALLOW_LIST (tables.js), each with why: any new
// one fails, and so does an entry no violation needs any more. A value entry
// also names the line, the skill's value and the tool's, and fails when any
// of them moves. The lint's logic is lint.js; its own tests are at the end of
// this file.

import { readdirSync, readFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { beforeAll, describe, expect, it } from 'vitest';
import { copilotSkills } from '../copilotSkills';
import { buildCopilotCards } from '../copilotCards';
import { buildCopilotTools } from '../copilotTools';
import { buildViewTools } from '../copilotViewTools';
import { buildAskUserTool } from '../llm/copilot/askUser';
import { buildReadSkillTool } from '../llm/copilot/skills';
import { buildRunTools } from '../llm/copilot/runTools';
import { buildFocusTools } from '../llm/copilot/viewTools';
import { shellCatalogue } from '../llm/context/shells';
import { handlePcaMessage } from '../workers/pcaKdeWorker';
import { bundledSkills, COPILOT_PAGES, loadDemoAssistantRun } from './modelVisible/helpers';
import {
    classifyIdentifier, compareAllowList, describeFinding, evaluateBinding, evaluateValue, firmRuleSentences, flatten, gateNumberCandidates,
    identifierMentions, lintSkills
} from './skillLint/lint';
import { ALLOW_LIST, BINDINGS, FIRM_RULES, NOT_GATES, VALUES } from './skillLint/tables';

const HERE = dirname(fileURLToPath(import.meta.url));
const SRC = join(HERE, '..');
const SNAPSHOTS = join(HERE, 'modelVisible/__snapshots__');

// --- what the Copilot offers and what the model sees ----------------------------------

/**
 * Every tool the Copilot can offer, from the builders, on the Demo run (all
 * of them are built there), and the probe context the value check measures
 * with: find_gr_peaks' and compute_bond_angles' data (each call once), the
 * bond-angle engine's payload and the run's shell catalogue.
 */
const offeredTools = async () => {
    const { run, rmc6fText } = await loadDemoAssistantRun();
    const request = (kind, params) => handlePcaMessage({ kind, ...params }, async () => rmc6fText);
    const cards = buildCopilotCards({ plotFiles: run.plotFiles, rValueFile: run.rValueFile, symmetry: run.symmetry });
    const tools = [
        ...buildRunTools(run, { cards }),
        ...buildFocusTools({ cards, pages: COPILOT_PAGES }),
        ...buildCopilotTools({ request, structure: run.structure, elements: Object.keys(run.structure.elementCounts || {}).sort() }),
        ...buildViewTools({ cards, pages: COPILOT_PAGES, activePage: () => 'dashboard' }),
        buildAskUserTool(),
        buildReadSkillTool(bundledSkills(), { load: copilotSkills.load })
    ].filter(Boolean);
    const byName = new Map(tools.map((tool) => [tool.name, tool]));
    const once = (compute) => {
        const seen = new Map();
        return (args) => {
            const key = JSON.stringify(args);
            if (!seen.has(key)) seen.set(key, compute(args));
            return seen.get(key);
        };
    };
    const probes = {
        tools: byName,
        gr: once(async (pair) => (await byName.get('find_gr_peaks').run({ pair })).data),
        bondAngles: once(async (args) => (await byName.get('compute_bond_angles').run(args)).data),
        triplets: once((params) => request('triplets', { binWidth: 1, ...params })),
        catalogue: shellCatalogue(run)
    };
    return { tools, probes };
};

const addSchemaKeys = (schema, into) => {
    for (const [key, value] of Object.entries(schema?.properties || {})) {
        into.add(key);
        addSchemaKeys(value, into);
        addSchemaKeys(value?.items, into);
    }
};
const addObjectKeys = (value, into) => {
    if (Array.isArray(value)) value.forEach((item) => addObjectKeys(item, into));
    else if (value && typeof value === 'object') {
        for (const [key, item] of Object.entries(value)) {
            into.add(key);
            addObjectKeys(item, into);
        }
    }
};
// A snake_case key in an object literal (`window_A: …`) or set on one (`block.first_peak_only = …`).
const SNAKE_KEYS = [/['"]?\b([a-z][a-z0-9]*(?:_[A-Za-z0-9]+)+)['"]?\s*:/g, /\.([a-z][a-z0-9]*(?:_[A-Za-z0-9]+)+)\s*=(?!=)/g];

/**
 * The names that are not tools but a skill may name: every tool parameter;
 * every key of the tool results and the run context the model saw on the
 * Demo (the model-visible snapshots, read_skill's own text left out); and the
 * snake_case keys of the modules that build them, for fields only some runs
 * have (first_shell.first_peak_only).
 */
const seenFields = (tools) => {
    const fields = new Set();
    tools.forEach((tool) => addSchemaKeys(tool.parameters, fields));
    for (const file of readdirSync(join(SNAPSHOTS, 'tool-results')).filter((name) => name.endsWith('.json') && !name.startsWith('read_skill'))) {
        addObjectKeys(JSON.parse(readFileSync(join(SNAPSHOTS, 'tool-results', file), 'utf8')).content, fields);
    }
    for (const file of ['dock/first-request.txt', 'harness/first-request.arm-B.txt']) {
        for (const match of readFileSync(join(SNAPSHOTS, file), 'utf8').matchAll(/"([A-Za-z_]\w*)"\s*:/g)) fields.add(match[1]);
    }
    const sources = ['copilotTools.js', 'llm/copilot/runTools.js', ...readdirSync(join(SRC, 'llm/context')).filter((name) => name.endsWith('.js')).map((name) => `llm/context/${name}`)];
    for (const source of sources) {
        const text = readFileSync(join(SRC, source), 'utf8');
        for (const pattern of SNAKE_KEYS) for (const match of text.matchAll(pattern)) fields.add(match[1]);
    }
    return fields;
};

/** The bundled skill files: SKILL.md (eager) and the references (loaded), by repository path. */
const bundledFiles = async () => {
    const files = { ...copilotSkills.files };
    for (const path of copilotSkills.references) files[path] = await copilotSkills.load(path);
    return files;
};

// --- the lint on the bundled skills ------------------------------------------------------

describe('the skill lint (AGENT_PLAN 5.1)', () => {
    let tools;
    let files;
    let result;

    beforeAll(async () => {
        const offered = await offeredTools();
        tools = offered.tools;
        files = await bundledFiles();
        result = await lintSkills({
            files,
            tools: new Set(tools.map((tool) => tool.name)),
            fields: seenFields(tools),
            bindings: BINDINGS,
            notGates: NOT_GATES,
            firmRules: FIRM_RULES,
            values: VALUES,
            probeContext: offered.probes
        });
    }, 120000);

    it('checks against every tool the dock and the harness send', () => {
        const offered = new Set(tools.map((tool) => tool.name));
        for (const file of ['dock/tools.json', 'harness/tools.arm-B.json']) {
            const sent = JSON.parse(readFileSync(join(SNAPSHOTS, file), 'utf8')).map((tool) => tool.function.name);
            expect(sent.filter((name) => !offered.has(name)), file).toEqual([]);
        }
        expect(offered.size).toBe(tools.length);
    });

    it('reads every skill and reference file the app bundles', () => {
        const paths = Object.keys(files);
        expect(paths.length).toBe(Object.keys(copilotSkills.files).length + copilotSkills.references.length);
        expect(paths.filter((path) => path.endsWith('/SKILL.md')).length).toBeGreaterThanOrEqual(4);
        // Each skill names at least one tool.
        for (const path of paths.filter((item) => item.endsWith('/SKILL.md'))) {
            expect(result.named.get(path.replace(/^\.claude\/skills\//, ''))?.size ?? 0, path).toBeGreaterThan(0);
        }
    });

    it.each([
        ['tool', 'every tool a skill names is one the Copilot offers'],
        ['gate', 'every gate number a skill states equals the constant that decides it'],
        ['firm-rule', 'every firm rule names a check or refusal that exists'],
        ['value', 'every value a skill quotes as measured on the Demo equals the tool\'s']
    ])('%s: %s (beyond the allow-list)', (check) => {
        const { unexpected } = compareAllowList(result.findings.filter((finding) => finding.check === check), ALLOW_LIST);
        expect(unexpected.map(describeFinding), `new ${check} violations: fix the skill or the code, or add an entry with why`).toEqual([]);
    });

    it('holds only today\'s violations in its allow-list', () => {
        const { stale } = compareAllowList(result.findings, ALLOW_LIST);
        expect(stale.map((entry) => `${entry.file} [${entry.check}] "${entry.fragment}"`), 'fixed: remove these allow-list entries').toEqual([]);
        for (const entry of ALLOW_LIST) {
            expect(entry.fragment.length, entry.fragment).toBeLessThanOrEqual(60);
            expect(['tool', 'gate', 'firm-rule', 'value']).toContain(entry.check);
            expect(entry.why.length, entry.fragment).toBeGreaterThan(20);
        }
    });

    it('names each value mismatch it holds: the line, the skill\'s value and the tool\'s', () => {
        for (const entry of ALLOW_LIST.filter((item) => item.check === 'value')) {
            const finding = result.findings.find((item) => item.check === 'value' && item.file === entry.file && item.context.includes(entry.fragment));
            expect(finding, entry.fragment).toBeDefined();
            expect({ line: finding.line, skill: finding.stated, tool: finding.measured }, `${entry.file}: ${entry.fragment}`)
                .toEqual({ line: entry.line, skill: entry.skill, tool: entry.tool });
        }
    });

    it.each(BINDINGS.filter((binding) => binding.verify).map((binding) => [`${binding.file}: ${binding.gate}`, binding]))(
        'the code applies the constant as bound — %s',
        (_, binding) => binding.verify()
    );
});

// --- the lint itself -----------------------------------------------------------------------

describe('the lint itself', () => {
    const tools = new Set(['find_gr_peaks', 'analyze_clusters']);
    const fields = new Set(['first_shell', 'window_A', 'window_confirmed', 'largest_abs_kurtosis_refs']);
    const lint = (text, tables = {}) => lintSkills({ files: { '.claude/skills/demo/SKILL.md': text }, tools, fields, ...tables });

    it('finds tool names in code spans and in the prose, and no mathematics', () => {
        const { flat } = flatten('Call `find_gr_peaks`, read `first_shell.window_A` and\n`largest_*_refs`; then analyze_clusters. |q_i| ≤ 0.25; see `data/5K_try1` and `*_FT_*`.');
        expect(identifierMentions(flat).map((mention) => mention.token))
            .toEqual(['find_gr_peaks', 'first_shell.window_A', 'largest_*_refs', 'analyze_clusters']);
        expect(classifyIdentifier('first_shell.window_A', { tools, fields }).kind).toBe('field');
        expect(classifyIdentifier('largest_*_refs', { tools, fields }).kind).toBe('field');
        expect(classifyIdentifier('find_gr_peaks.first_shell', { tools, fields })).toEqual({ kind: 'tool', tool: 'find_gr_peaks' });
    });

    it('reports a tool that does not exist', async () => {
        const { findings, named } = await lint('Run `analyze_local_order` after `find_gr_peaks` with `window_confirmed: true`.');
        expect(findings.map((finding) => [finding.check, finding.kind, finding.fragment])).toEqual([['tool', 'unknown', 'analyze_local_order']]);
        expect([...named.get('demo/SKILL.md')]).toEqual(['find_gr_peaks']);
    });

    it('reads comparator numbers, not direction families or plain values', () => {
        const { flat } = flatten('Peaks at ≥ 3σ, below ~0.1 pairs, fewer than 60 %, under 1.2× isotropic, from 1.25 × the floor; <111> and <100>; z = 11; 19 % vs 14 %.');
        expect(gateNumberCandidates(flat).map((item) => flat.slice(item.at, item.end))).toEqual(['3', '0.1', '60', '1.2', '1.25']);
    });

    it('binds a number to its constant: equal, different, missing, stale, wrong statistic', () => {
        const binding = { gate: 'a direction', find: 'Preferred direction only with a', number: /significance ≥ (?<n>\d+)σ/, constants: [{ name: 'S', value: 3, module: 'm.js' }] };
        const read = (text, extra = {}) => evaluateBinding({ ...binding, ...extra }, flatten(text).flat);
        expect(read('Preferred direction only with a significance ≥ 3σ.').problem).toBeNull();
        expect(read('Preferred direction only with a significance ≥ 4σ.')).toMatchObject({ problem: 'differs', stated: 4 });
        expect(read('Preferred direction only with a significant peak.').problem).toBe('no-number');
        expect(read('A direction needs a peak.').problem).toBe('stale');
        expect(read('Preferred direction only with a significance ≥ 3σ.', { statistic: { pattern: /peak/, says: 'the peak test' } }).problem).toBe('statistic');
        expect(read('Preferred direction only with a peak\nsignificance ≥ 3σ.', { statistic: { pattern: /peak/, says: 'the peak test' } }).problem).toBeNull();
        const percent = { gate: 'a share', find: 'of the atoms', before: 20, number: /fewer than (?<n>\d+) % of the atoms/, unit: '%', constants: [{ name: 'F', value: 0.6, module: 'm.js' }] };
        expect(evaluateBinding(percent, flatten('Refused when fewer than 60 % of the atoms agree.').flat).problem).toBeNull();
    });

    it('requires every comparator number to be bound or listed as no gate', async () => {
        const text = 'Join a peak within 0.5 Å. The trans angle is a shoulder below 180°.';
        const { findings } = await lint(text, { notGates: [{ file: 'demo/SKILL.md', find: 'a shoulder below 180°', why: 'geometry' }] });
        expect(findings.map((finding) => [finding.check, finding.kind])).toEqual([['gate', 'unbound']]);
        expect(findings[0].context).toContain('within 0.5 Å');
        const { findings: stale } = await lint(text, { notGates: [{ file: 'demo/SKILL.md', find: 'below 170°', why: 'geometry' }] });
        expect(stale.map((finding) => finding.kind)).toEqual(['stale', 'unbound', 'unbound']);
    });

    it('maps firm rules to probes, and reports an unmapped, a stale and a broken one', async () => {
        const text = '- *Firm rule in the\napp:* an answer that says X gets a warning.\n- *Firm rule in the app (THY-18):* the Outcome is computed.\n';
        expect(firmRuleSentences(flatten(text).flat).map((rule) => rule.text)).toEqual(['an answer that says X gets a warning.', 'the Outcome is computed.']);
        const firmRules = [
            { file: 'demo/SKILL.md', find: 'says X gets a warning', enforcedBy: 'rule x', probe: () => { throw new Error('rule x is gone'); } },
            { file: 'demo/SKILL.md', find: 'says Y', enforcedBy: 'rule y', probe: () => {} }
        ];
        const { findings } = await lint(text, { firmRules });
        expect(findings.map((finding) => [finding.check, finding.kind])).toEqual([['firm-rule', 'unmapped'], ['firm-rule', 'broken'], ['firm-rule', 'stale']]);
        expect(findings[0].context).toBe('the Outcome is computed.');
        expect(findings[1].detail).toContain('rule x is gone');
    });

    it('reads a measured value as written, and compares it at the skill\'s own precision', async () => {
        const entry = { what: 'the first peak\'s CN', find: 'the first peak holds', number: /holds (?<n>\d+(?:\.\d+)?) atoms/, tool: 't' };
        const read = (text, measured, extra = {}) => evaluateValue({ ...entry, ...extra }, flatten(text).flat, measured);
        expect(read('Here the first peak holds 3.9 atoms.', 3.25)).toMatchObject({ problem: 'differs', stated: '3.9', measured: '3.25' });
        expect(read('Here the first peak holds 3.3 atoms.', 3.25).problem).toBeNull();
        // A count written whole agrees within half an atom.
        expect(read('Here the first peak holds 3 atoms.', 2.75).problem).toBeNull();
        expect(read('Here the first peak holds 3 atoms.', 2.4).problem).toBe('differs');
        expect(read('No such sentence.', 3).problem).toBe('stale');
        expect(read('Here the first peak holds many atoms.', 3).problem).toBe('no-number');
        // A window, and a share in %.
        const window = { what: 'w', find: 'window', number: /window (?<n>\d+\.\d+)–(?<m>\d+\.\d+) Å/, tool: 't' };
        expect(evaluateValue(window, 'the window 2.33–2.71 Å', [2.21, 2.65])).toMatchObject({ problem: 'differs', stated: '2.33–2.71', measured: '2.21–2.65' });
        expect(evaluateValue(window, 'the window 2.21–2.65 Å', [2.21, 2.65]).problem).toBeNull();
        const share = { what: 's', find: 'share', number: /share (?<n>\d+) %/, unit: '%', tool: 't' };
        expect(evaluateValue(share, 'a share 47 %', 0.4748).problem).toBeNull();
        expect(evaluateValue(share, 'a share 47 %', 0.4)).toMatchObject({ problem: 'differs', stated: '47 %', measured: '40 %' });
        // In the lint: a finding with the line of the number, and a measurement that throws is one too.
        const text = 'Line one.\nHere the first peak holds\n3.9 atoms.\n';
        const { findings } = await lintSkills({
            files: { '.claude/skills/demo/SKILL.md': text },
            tools: new Set(),
            fields: new Set(),
            values: [
                { ...entry, file: 'demo/SKILL.md', measure: () => 3.25 },
                { ...entry, file: 'demo/SKILL.md', what: 'x', measure: () => { throw new Error('no run'); } }
            ]
        });
        expect(findings.map((finding) => [finding.check, finding.kind, finding.line, finding.stated, finding.measured]))
            .toEqual([['value', 'differs', 3, '3.9', '3.25'], ['value', 'broken', 2, undefined, undefined]]);
    });

    it('compares findings with the allow-list both ways', () => {
        const findings = [{ check: 'gate', file: 'a/SKILL.md', context: 'only when the angle is small' }];
        const entry = { check: 'gate', file: 'a/SKILL.md', fragment: 'angle is small', why: 'x' };
        const other = { check: 'gate', file: 'a/SKILL.md', fragment: 'not in the text', why: 'x' };
        expect(compareAllowList(findings, [entry])).toEqual({ unexpected: [], stale: [] });
        expect(compareAllowList(findings, [other])).toEqual({ unexpected: findings, stale: [other] });
    });
});
