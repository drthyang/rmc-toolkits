// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Plan item 1.6, part (i): the report beside the case score — its parts (the
// feature-only score), every component against A+, and every disagreement
// with the planted truth — with the case score unchanged. Dev set and
// synthetic records only: the held-out test cases and their records are never
// read here.

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES } from '../cases';
import { restrictedScore, scoreComponents, scoreVerdict } from '../scoring';
import { normalizeVerdict, truthVerdict } from '../verdict';
import { compareWithReference, disagreementsOf, fmtDifference, renderComparison } from '../comparison';
import { renderReport } from '../report';

// Every dev record re-scores as recorded and its components are its own
// (ALL); the report's pins below are on the ten dev cases outside the
// judgment subset (BASELINE), the table copilot-validation.md quotes: the
// subset's results stay out of that document, and are pinned apart
// (arms.test.js, judgment.test.js; plan item 1.6 ii).
const ALL = JSON.parse(readFileSync(fileURLToPath(new URL('../baseline/arm-a-dev.json', import.meta.url)), 'utf8'));
const BASELINE = ALL.filter((record) => !record.truth?.judgment);
const PAIRED_CASES = DEV_CASES.filter((item) => !item.judgment);
const A = 'no model · A';
const APLUS = 'no model · A+';
const B = 'scripted · B';

const recordOf = (id, arm) => ALL.find((record) => record.case === id && record.arm === arm);
const truthOf = (id) => recordOf(id, 'D').truth;
const verdict = (raw) => normalizeVerdict({ features: [], space_group: null, chiral: null, centrosymmetric: null, polar: null, ...raw });

// A live record as the harness writes it, with only what the report reads.
const live = (id, raw, { run = 1, arm = 'B', model = 'scripted', error = null } = {}) => {
    const truth = truthOf(id);
    const parsed = raw == null ? null : verdict(raw);
    return {
        id: `${id}__${arm}__${model}__r${run}`, case: id, split: 'dev', arm, model: { name: model },
        verdict: parsed, verdict_error: parsed ? null : 'no VERDICT line', error, score: scoreVerdict(truth, parsed), truth
    };
};
// B giving A+'s own verdict on a case (n runs).
const asAplus = (id, runs = 1) => Array.from({ length: runs }, (_, k) => live(id, recordOf(id, 'A+').verdict, { run: k + 1 }));
const symmetryOf = (id) => {
    const { space_group: group, chiral, centrosymmetric, polar } = truthOf(id).symmetry;
    return { space_group: group, chiral, centrosymmetric, polar };
};

const meanOf = (components) => Number((components.reduce((sum, c) => sum + c.value, 0) / components.length).toFixed(4));
const part = (result, column) => result.parts.find((item) => item.column === column);
const row = (result, id, component) => result.components.find((item) => item.case === id && item.id === component);
const disagreementsIn = (result, id) => result.disagreements.filter((item) => item.case === id)
    .map((item) => `${item.where} | ${item.answers.join('; ')} | ${item.made.map((m) => `${m.column} ${m.runs}/${m.of}`).join(', ')}`);

describe('the case score is unchanged, and its components are its own', () => {
    it('the dev baselines re-score exactly as recorded, against the dev truths as cases.js builds them today', () => {
        // Ten cases and the judgment subset's five, each with A, A+ and D.
        expect(BASELINE).toHaveLength(30);
        expect(ALL).toHaveLength(45);
        for (const record of ALL) {
            expect(record.split, record.id).toBe('dev');
            expect(JSON.parse(JSON.stringify(scoreVerdict(record.truth, record.verdict))), record.id).toEqual(record.score);
        }
        for (const item of DEV_CASES) {
            expect(JSON.parse(JSON.stringify(buildCase(item, { partials: false }).truth)), item.id).toEqual(truthOf(item.id));
        }
    });

    it('scoreComponents: their mean is the case score, for every arm, the truth, no verdict and the synthetic verdicts', () => {
        const scores = [
            ...ALL.map((record) => [record.id, record.score, record.truth]),
            ...DEV_CASES.map((item) => [`${item.id} truth`, scoreVerdict(truthOf(item.id), truthVerdict(truthOf(item.id))), truthOf(item.id)]),
            ...DEV_CASES.map((item) => [`${item.id} none`, scoreVerdict(truthOf(item.id), null), truthOf(item.id)]),
            ['offcentre magnitude off', live('dev-offcentre', { features: [{ kind: 'off_centre', element: 'Nb', direction_uvw: [0, 0, 1], magnitude_A: 0.3 }] }).score, truthOf('dev-offcentre')],
            ['null over-claimed', live('dev-thermal-null', { features: [{ kind: 'off_centre', element: 'Ti' }], ...symmetryOf('dev-thermal-null') }).score, truthOf('dev-thermal-null')]
        ];
        for (const [what, score, truth] of scores) {
            const components = scoreComponents(score, truth);
            expect(meanOf(components), what).toBe(score.score);
            // Without the truth the features are named from the score: the same components.
            expect(scoreComponents(score).map((c) => c.value), what).toEqual(components.map((c) => c.value));
        }
        expect(scoreComponents(recordOf('dev-jahn-teller', 'A+').score, truthOf('dev-jahn-teller'))).toEqual([
            { id: 'feature jahn_teller Cu', group: 'feature', kind: 'jahn_teller', element: 'Cu', value: 0.5 },
            { id: 'no spurious claim', group: 'spurious', value: 1 },
            { id: 'symmetry space_group', group: 'symmetry', field: 'space_group', value: 1 },
            { id: 'symmetry chiral', group: 'symmetry', field: 'chiral', value: 1 },
            { id: 'symmetry centrosymmetric', group: 'symmetry', field: 'centrosymmetric', value: 1 },
            { id: 'symmetry polar', group: 'symmetry', field: 'polar', value: 1 }
        ]);
        // No verdict: every component 0, the spurious-claim test included.
        expect(scoreComponents(scoreVerdict(truthOf('dev-thermal-null'), null)).every((c) => c.value === 0)).toBe(true);
    });

    it('restrictedScore: the feature-only score is the planted features alone; null for a case with none', () => {
        const components = scoreComponents(recordOf('dev-zone-boundary', 'A').score, truthOf('dev-zone-boundary'));
        expect(restrictedScore(components, 'feature')).toBe(0);
        expect(restrictedScore(components, 'spurious')).toBe(0);
        expect(restrictedScore(components, 'symmetry')).toBe(1);
        expect(restrictedScore(scoreComponents(recordOf('dev-thermal-null', 'A').score), 'feature')).toBeNull();
        expect(restrictedScore(scoreComponents(recordOf('dev-breathing', 'A+').score), 'symmetry')).toBeNull();
    });
});

describe('the dev baselines through the report (A, A+ and D; no B runs yet)', () => {
    const result = compareWithReference(BASELINE);

    it('the parts: the feature-only score separates A from A+ where the case score barely does', () => {
        expect(result.reference).toBe(APLUS);
        expect(result.live).toEqual([]);
        expect(result.differences).toEqual([]);
        const a = part(result, A);
        const sweep = part(result, APLUS);
        expect(a.cases).toBe(10);
        expect(a.score).toBeCloseTo(0.825, 4);
        expect(sweep.score).toBeCloseTo(0.99167, 4);
        expect(a.groups.feature).toEqual({ value: 0.25, cases: 6 });
        expect(sweep.groups.feature.value).toBeCloseTo(5.5 / 6, 10);
        expect(a.groups.spurious).toEqual({ value: 0.9, cases: 10 });
        expect(sweep.groups.spurious).toEqual({ value: 1, cases: 10 });
        expect(a.groups.symmetry.value).toBeCloseTo(8.25 / 9, 10);
        expect(part(result, 'no model · D').groups).toEqual({ feature: { value: 1, cases: 6 }, spurious: { value: 1, cases: 10 }, symmetry: { value: 1, cases: 9 } });
    });

    it('every component, A beside A+ (D left out), 52 in all', () => {
        expect(result.shown).toEqual([A, APLUS]);
        expect(result.components).toHaveLength(52);
        expect(row(result, 'dev-jahn-teller', 'feature jahn_teller Cu').values).toEqual({ [A]: 0.5, [APLUS]: 0.5 });
        expect(row(result, 'dev-zone-boundary', 'no spurious claim').values).toEqual({ [A]: 0, [APLUS]: 1 });
        expect(row(result, 'dev-breathing', 'symmetry space_group')).toBeUndefined();
    });

    it('every disagreement: A+ has one (the JT long axis), A the misleading pages', () => {
        expect(result.disagreements.filter((item) => item.made.some((m) => m.column === APLUS)).map((item) => `${item.case} ${item.where}`))
            .toEqual(['dev-jahn-teller direction of jahn_teller Cu']);
        expect(result.disagreements.filter((item) => item.made.some((m) => m.column === 'no model · D'))).toEqual([]);
        expect(disagreementsIn(result, 'dev-zone-boundary')).toEqual([
            `missed: modulation Ti | read as off_centre Ti (wrong kind) | ${A} 1/1`,
            `spurious claim: off_centre Ti | [-0.03 0.1 1], 0.0987 Å | ${A} 1/1`
        ]);
        expect(disagreementsIn(result, 'dev-symmetry-lowering')).toEqual([
            `missed: symmetry_lowering Ta | not claimed | ${A} 1/1`,
            `space group | Pm-3m | ${A} 1/1`,
            `centrosymmetric | true | ${A} 1/1`,
            `polar | false | ${A} 1/1`
        ]);
        expect(disagreementsIn(result, 'dev-jahn-teller')).toEqual([`direction of jahn_teller Cu | not given | ${A} 1/1, ${APLUS} 1/1`]);
        expect(result.clean).toEqual(['dev-chiral', 'dev-chirality-trap', 'dev-fit-unconverged', 'dev-offcentre', 'dev-thermal-null']);
    });

    it('each disagreement costs a component, and each component below 1 has one', () => {
        for (const record of ALL) {
            const lost = scoreComponents(record.score, record.truth).filter((c) => c.value < 1).length;
            const listed = disagreementsOf(record);
            expect(listed.length >= lost && (lost > 0) === (listed.length > 0), record.id).toBe(true);
        }
    });

    it('the rendered report: the parts, the table without differences, the listing, and no live columns', () => {
        const report = renderReport(BASELINE);
        expect(report).toContain('### The case score\'s parts');
        expect(report).toContain(`| ${A} | 10 | 0.82 | 0.25 (6 cases) | 0.90 | 0.92 (9 cases) |`);
        expect(report).toContain(`| ${APLUS} | 10 | 0.99 | 0.92 (6 cases) | 1.00 | 1.00 (9 cases) |`);
        expect(report).toContain('No live run (arm B or C) in this split yet');
        expect(report).toContain(`| Case | Component | ${A} | ${APLUS} |`);
        expect(report).toContain('| dev-symmetry-lowering | symmetry polar | 0.00 | 1.00 |');
        expect(report).toContain(`| dev-jahn-teller | direction of jahn_teller Cu | not given | [0 0 1] or any cubic image, within 15° | ${A} (1/1), ${APLUS} (1/1) |`);
        expect(report).toContain('| dev-breathing | missed: cluster_breathing Mo | not claimed | planted | no model · A (1/1) |');
        expect(report).toContain('No disagreement in any run: dev-chiral, dev-chirality-trap, dev-fit-unconverged, dev-offcentre, dev-thermal-null.');
        expect(report).not.toMatch(/ − A\+ \|/);
        // The existing tables are still there, unchanged in form.
        expect(report).toContain(`| ${A} | 0.82 | 0/4 | – | 0 | – | 0/0 |`);
    });
});

describe('B against A+ on synthetic records', () => {
    const withBaseline = (records) => compareWithReference([...BASELINE, ...records]);

    it('B matching A+ (3 runs a case): every difference 0, and B shares A+\'s one disagreement', () => {
        const result = withBaseline(PAIRED_CASES.flatMap((item) => asAplus(item.id, 3)));
        expect(result.live).toEqual([B]);
        expect(result.differences).toEqual([{ column: B, cases: 10, score: expect.any(Number), groups: { feature: 0, spurious: 0, symmetry: 0 } }]);
        expect(result.differences[0].score).toBeCloseTo(0, 10);
        expect(result.components.every((item) => item.differences[B] === 0)).toBe(true);
        expect(disagreementsIn(result, 'dev-jahn-teller')).toEqual([`direction of jahn_teller Cu | not given | ${A} 1/1, ${APLUS} 1/1, ${B} 3/3`]);
        const lines = renderComparison([...BASELINE, ...PAIRED_CASES.flatMap((item) => asAplus(item.id, 3))]);
        expect(lines).toContain(`| ${B} − A+ | 10 paired | 0.00 | 0.00 | 0.00 | 0.00 |`);
        expect(lines).toContain(`- ${B} against A+: 0 components higher, 0 lower, 52 the same.`);
        expect(lines).toContain(`| Case | Component | ${A} | ${APLUS} | ${B} | ${B} − A+ |`);
    });

    it('B over-claiming on the no-order case: the spurious-claim test −1, a spurious claim listed with B', () => {
        const over = live('dev-thermal-null', { features: [{ kind: 'off_centre', element: 'Ti', direction_uvw: [1, 1, 2], magnitude_A: 0.05 }], ...symmetryOf('dev-thermal-null') });
        const result = withBaseline([over]);
        expect(result.differences[0]).toMatchObject({ column: B, cases: 1, groups: { feature: null, spurious: -1, symmetry: 0 } });
        expect(result.differences[0].score).toBeCloseTo(-0.2, 10);
        expect(row(result, 'dev-thermal-null', 'no spurious claim').differences[B]).toBe(-1);
        expect(disagreementsIn(result, 'dev-thermal-null')).toEqual([`spurious claim: off_centre Ti | [1 1 2], 0.05 Å | ${B} 1/1`]);
        expect(result.disagreements.find((item) => item.case === 'dev-thermal-null').truth).toBe('no order planted');
        const lines = renderComparison([...BASELINE, over]);
        expect(lines).toContain(`| ${B} − A+ | 1 paired | −0.20 | – | −1.00 | 0.00 |`);
        expect(lines).toContain(`| dev-thermal-null | no spurious claim | 1.00 | 1.00 | 0.00 | **−1.00** |`);
        expect(lines).toContain(`- ${B} against A+: 0 components higher, 1 lower, 4 the same.`);
    });

    it('B missing one feature: the feature-only score −1 where the case score moves by half', () => {
        const missing = live('dev-breathing', { features: [] });
        const result = withBaseline([missing]);
        expect(result.differences[0]).toMatchObject({ cases: 1, score: -0.5, groups: { feature: -1, spurious: 0, symmetry: null } });
        expect(disagreementsIn(result, 'dev-breathing')).toEqual([`missed: cluster_breathing Mo | not claimed | ${A} 1/1, ${B} 1/1`]);
    });

    it('B better than A+ on the one component A+ misses (the JT long axis): +0.50 there, A+ alone keeps the disagreement', () => {
        const aplus = recordOf('dev-jahn-teller', 'A+').verdict;
        const better = live('dev-jahn-teller', {
            ...aplus, features: aplus.features.map((f) => (f.kind === 'jahn_teller' ? { ...f, direction_uvw: [1, 0, 0] } : f))
        });
        expect(better.score.score).toBe(1);
        const result = withBaseline([better]);
        expect(row(result, 'dev-jahn-teller', 'feature jahn_teller Cu')).toMatchObject({ values: { [A]: 0.5, [APLUS]: 0.5, [B]: 1 }, differences: { [B]: 0.5 } });
        expect(result.differences[0].groups.feature).toBe(0.5);
        expect(result.differences[0].score).toBeCloseTo(1 - 0.9167, 10);
        expect(disagreementsIn(result, 'dev-jahn-teller')).toEqual([`direction of jahn_teller Cu | not given | ${A} 1/1, ${APLUS} 1/1`]);
        const lines = renderComparison([...BASELINE, better]);
        expect(lines).toContain(`| dev-jahn-teller | feature jahn_teller Cu | 0.50 | 0.50 | 1.00 | **+0.50** |`);
        expect(lines).toContain(`| ${B} − A+ | 1 paired | +0.08 | +0.50 | 0.00 | 0.00 |`);
    });

    it('runs that differ: each component is the mean over the runs, each disagreement counts its runs', () => {
        const runs = [
            live('dev-zone-boundary', { features: [{ kind: 'modulation', element: 'Ti', q: [0, 0, 0.5], direction_uvw: [0, 0, 1] }], ...symmetryOf('dev-zone-boundary') }, { run: 1 }),
            live('dev-zone-boundary', { features: [{ kind: 'modulation', element: 'Ti', q: [0, 0, 0.375], direction_uvw: [0, 0, 1] }], ...symmetryOf('dev-zone-boundary') }, { run: 2 }),
            live('dev-zone-boundary', { features: [{ kind: 'modulation', element: 'Sr', q: [0.125, 0, 0] }], ...symmetryOf('dev-zone-boundary') }, { run: 3 })
        ];
        const result = withBaseline(runs);
        expect(row(result, 'dev-zone-boundary', 'feature modulation Ti').values[B]).toBeCloseTo(0.5, 10);
        expect(row(result, 'dev-zone-boundary', 'no spurious claim').values[B]).toBeCloseTo(2 / 3, 10);
        expect(disagreementsIn(result, 'dev-zone-boundary')).toEqual([
            `missed: modulation Ti | read as off_centre Ti (wrong kind); read as modulation Sr (wrong element) | ${A} 1/1, ${B} 1/3`,
            `q of modulation Ti | (0, 0, 0.375) | ${B} 1/3`,
            `spurious claim: modulation Sr | q (0.125, 0, 0) | ${B} 1/3`,
            `spurious claim: off_centre Ti | [-0.03 0.1 1], 0.0987 Å | ${A} 1/1`
        ]);
    });

    it('a magnitude off, a symmetry field unanswered, a run with no verdict', () => {
        const runs = [
            live('dev-offcentre', { features: [{ kind: 'off_centre', element: 'Nb', direction_uvw: [0, 0, 1], magnitude_A: 0.3 }], ...symmetryOf('dev-offcentre'), polar: null }),
            live('dev-chiral', null, { error: { name: 'ContextLimitError', message: 'refused' } })
        ];
        const result = withBaseline(runs);
        expect(disagreementsIn(result, 'dev-offcentre')).toEqual([`magnitude of off_centre Nb | 0.3 Å | ${B} 1/1`, `polar | not answered | ${B} 1/1`]);
        expect(result.disagreements.find((item) => item.where === 'magnitude of off_centre Nb').truth).toBe('0.2 ± 0.05 Å');
        expect(disagreementsIn(result, 'dev-chiral')).toEqual([`no parseable verdict | no VERDICT line | ${B} 1/1`]);
        expect(result.components.filter((item) => item.case === 'dev-chiral').map((item) => item.differences[B])).toEqual([-1, -1, -1, -1, -1]);
    });

    it('a split with no A+ record says so, and gives no differences', () => {
        const lines = renderComparison([live('dev-offcentre', recordOf('dev-offcentre', 'A+').verdict)]);
        expect(lines).toContain('No A+ record in this split: the comparison with the engines needs `run.js arm-a` records beside the runs.');
        expect(lines.join('\n')).not.toContain('### Component by component');
        expect(lines).toContain(`| ${B} | 1 | 1.00 | 1.00 | 1.00 | 1.00 |`);
        expect(lines).toContain('No run disagrees with the planted truth.');
        expect(renderComparison([])).toEqual([]);
    });

    it('a difference prints signed, with a true minus, or – when a side is missing', () => {
        expect([0.5, -1, 0, 1e-9, null].map(fmtDifference)).toEqual(['+0.50', '−1.00', '0.00', '0.00', '–']);
    });
});
