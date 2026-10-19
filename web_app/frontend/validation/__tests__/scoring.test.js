// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import { describe, expect, it } from 'vitest';
import { buildCase, DEV_CASES, TEST_CASE_IDS } from '../cases';
import { holdoutEnabled, loadTestCases } from '../holdoutLoader';
import {
    directionError, reproducibility, sameWavevector, scoreVerdict, TOLERANCES, verdictSignature
} from '../scoring';
import { normalizeVerdict, parseVerdict, truthVerdict } from '../verdict';

// Every case the suite may build: the dev set, and the sealed test set where
// the loader may open it (CI's RMC_VALIDATION_HOLDOUT=1).
const HOLDOUT = holdoutEnabled();
const CASES = [...DEV_CASES, ...(HOLDOUT ? await loadTestCases() : [])];
const truths = Object.fromEntries(CASES.map((item) => [item.id, buildCase(item, { partials: false }).truth]));
const verdict = (raw) => normalizeVerdict({ features: [], space_group: null, chiral: null, centrosymmetric: null, polar: null, ...raw });

describe('scoring against the planted truth', () => {
    it.skipIf(!HOLDOUT)('arm D reaches every held-out case', () => {
        expect(TEST_CASE_IDS.every((id) => id in truths)).toBe(true);
    });

    it('the truth itself (arm D) scores 1 on every case, dev and test', () => {
        for (const [id, truth] of Object.entries(truths)) {
            const score = scoreVerdict(truth, truthVerdict(truth));
            expect(score.score, id).toBe(1);
            expect(score.spurious, id).toEqual([]);
        }
    });

    it('a direction: ±u are one answer for a split site, a cubic image for a JT axis, 15° the limit', () => {
        expect(directionError([0, 0, -1], [0, 0, 1], 'axis')).toBeCloseTo(0, 6);
        expect(directionError([1, 0, 0], [0, 0, 1], 'axis')).toBeCloseTo(90, 6);
        expect(directionError([1, 0, 0], [0, 0, 1], 'cubic')).toBeCloseTo(0, 6);
        expect(directionError([1, 1, 1], [0, 0, 1], 'cubic')).toBeCloseTo(54.74, 1);
        const truth = truths['dev-offcentre'];
        const tilted = (deg) => [Math.sin((deg * Math.PI) / 180), 0, Math.cos((deg * Math.PI) / 180)];
        const at = (deg) => scoreVerdict(truth, verdict({ features: [{ kind: 'off_centre', element: 'Nb', direction_uvw: tilted(deg), magnitude_A: 0.2 }] })).features[0];
        expect(at(TOLERANCES.DIRECTION_DEG - 1).direction_ok).toBe(true);
        expect(at(TOLERANCES.DIRECTION_DEG + 1).direction_ok).toBe(false);
    });

    it('a magnitude: 25 % or 0.03 Å, whichever is larger', () => {
        const truth = truths['dev-offcentre'];
        const with_ = (magnitude) => scoreVerdict(truth, verdict({ features: [{ kind: 'off_centre', element: 'Nb', direction_uvw: [0, 0, 1], magnitude_A: magnitude }] }));
        expect(with_(0.24).features[0].magnitude_ok).toBe(true);
        expect(with_(0.26).features[0].magnitude_ok).toBe(false);
        expect(with_(0.26).score).toBeLessThan(1);
        const small = truths['dev-symmetry-lowering'];
        const shifted = (magnitude) => scoreVerdict(small, verdict({ features: [{ kind: 'symmetry_lowering', element: 'Ta', direction_uvw: [0, 0, 1], magnitude_A: magnitude }] })).features[0];
        expect(shifted(0.099).magnitude_ok).toBe(true);
        expect(shifted(0.105).magnitude_ok).toBe(false);
    });

    it('a wavevector: q ≡ −q ≡ q + G, nothing else', () => {
        expect(sameWavevector([0, 0, -0.5], [0, 0, 0.5])).toBe(true);
        expect(sameWavevector([1, 0, 0.5], [0, 0, 0.5])).toBe(true);
        expect(sameWavevector([0, 0, 0.49], [0, 0, 0.5])).toBe(true);
        expect(sameWavevector([0.5, 0, 0], [0, 0, 0.5])).toBe(false);
        expect(sameWavevector([0, 0, 0.375], [0, 0, 0.5])).toBe(false);
    });

    it('the misread zone-boundary case: an ordered mode called random off-centring misses and over-claims', () => {
        const truth = truths['dev-zone-boundary'];
        const misread = scoreVerdict(truth, verdict({ features: [{ kind: 'off_centre', element: 'Ti', direction_uvw: [0, 0, 1], magnitude_A: 0.1 }], ...truth.symmetry }));
        expect(misread.features[0].found).toBe(false);
        expect(misread.spurious).toHaveLength(1);
        expect(misread.contradictions).toEqual([{ what: 'spurious off_centre (Ti)', rule: 'direction' }]);
        const right = scoreVerdict(truth, verdict({ features: [{ kind: 'modulation', element: 'Ti', q: [0, 0, 0.5], direction_uvw: [0, 0, 1] }], ...truth.symmetry }));
        expect(right.score).toBe(1);
        // The decoy long wave named as order is spurious.
        const decoy = scoreVerdict(truth, verdict({ features: [{ kind: 'modulation', element: 'Ti', q: [0, 0, 0.5], direction_uvw: [0, 0, 1] }, { kind: 'modulation', element: 'Sr', q: [0.125, 0, 0] }], ...truth.symmetry }));
        expect(decoy.spurious.map((c) => c.element)).toEqual(['Sr']);
    });

    it('no order: any claim over-claims; a JT case tolerates true partial descriptions', () => {
        const nul = truths['dev-thermal-null'];
        expect(scoreVerdict(nul, verdict({ ...nul.symmetry })).score).toBe(1);
        const over = scoreVerdict(nul, verdict({ features: [{ kind: 'off_centre', element: 'Ti', direction_uvw: [1, 1, 2] }], ...nul.symmetry }));
        expect(over.overclaim).toBe(true);
        expect(over.score).toBeLessThan(1);
        const jt = truths['dev-jahn-teller'];
        const partial = scoreVerdict(jt, verdict({ features: [
            { kind: 'jahn_teller', element: 'Cu', magnitude_A: 0.3, direction_uvw: [1, 0, 0] },
            { kind: 'off_centre', element: 'F', direction_uvw: [1, 0, 0] }
        ], ...jt.symmetry }));
        expect(partial.score).toBe(1);
        expect(partial.tolerated.map((c) => c.kind)).toEqual(['off_centre']);
    });

    it('symmetry fields: wrong and unanswered both cost; wrong positive claims name their answer check', () => {
        const trap = truths['dev-chirality-trap'];
        const wrong = scoreVerdict(trap, verdict({ space_group: 'F-43m', chiral: true, centrosymmetric: false, polar: true }));
        expect(wrong.symmetry).toEqual({ space_group: 'right', chiral: 'wrong', centrosymmetric: 'right', polar: 'wrong' });
        expect(wrong.contradictions.map((c) => c.rule)).toEqual(['chirality', 'polar']);
        const silent = scoreVerdict(trap, verdict({ space_group: 'F -4 3 m' }));
        expect(silent.symmetry.space_group).toBe('right');
        expect(silent.symmetry.chiral).toBe('unanswered');
        expect(scoreVerdict(truths['dev-chiral'], verdict({ space_group: 'P2₁3', chiral: true, centrosymmetric: false, polar: false })).score).toBe(1);
    });

    it('no verdict scores 0; reproducibility compares verdict signatures', () => {
        expect(scoreVerdict(truths['dev-offcentre'], null)).toMatchObject({ format_ok: false, score: 0 });
        const a = verdict({ features: [{ kind: 'modulation', element: 'Ti', q: [0, 0, 0.5] }], space_group: 'Pm-3m' });
        const b = verdict({ features: [{ kind: 'modulation', element: 'ti', q: [0, 0, -0.5] }], space_group: 'Pm−3m' });
        expect(verdictSignature(a)).toBe(verdictSignature(b));
        const repro = reproducibility([{ verdict: a, score: { score: 1 } }, { verdict: b, score: { score: 1 } }, { verdict: null, score: { score: 0 } }]);
        expect(repro).toMatchObject({ runs: 3, identical: false, pair_agreement: 0.333 });
    });
});

describe('the verdict line', () => {
    it('reads the last VERDICT object, in a fence or bold, and normalises its fields', () => {
        const text = 'VERDICT {"features": []}\n\nOutcome: achieved\n\n**VERDICT** ```json\n{"features": [{"kind": "Off-center", "element": "Nb", '
            + '"direction_uvw": [0, 0, 1], "magnitude_A": "0.2", "q": null}], "space_group": "Pm-3m", "chiral": false, '
            + '"centrosymmetric": true, "polar": "no"}\n```';
        const { verdict: parsed, error } = parseVerdict(text);
        expect(error).toBeNull();
        expect(parsed).toEqual({
            features: [{ kind: 'off_centre', element: 'Nb', direction_uvw: [0, 0, 1], magnitude_A: 0.2 }],
            space_group: 'Pm-3m', chiral: false, centrosymmetric: true, polar: null
        });
    });

    it('a missing, unclosed or malformed verdict is an error, never a guess', () => {
        expect(parseVerdict('No order here.').error).toBe('no VERDICT line');
        expect(parseVerdict('VERDICT {"features": [').error).toMatch(/not closed/);
        expect(parseVerdict('VERDICT {"features": [{"kind": "twinning"}]}').error).toMatch(/unknown kind 'twinning'/);
        expect(parseVerdict('VERDICT {"order": "none"}').error).toMatch(/no features array/);
    });
});

describe('the computed verdict against the score (THY-18)', () => {
    it('earned, over, under, held', async () => {
        const { outcomeAgreement } = await import('../scoring');
        expect(outcomeAgreement('achieved', 1)).toEqual({ status: 'achieved', right: true, agreement: 'earned' });
        expect(outcomeAgreement('achieved', 0.75)).toMatchObject({ agreement: 'over' });
        expect(outcomeAgreement('partial', 1)).toMatchObject({ agreement: 'under' });
        expect(outcomeAgreement('not verified', 0.5)).toMatchObject({ agreement: 'held' });
        expect(outcomeAgreement(null, 0)).toMatchObject({ agreement: 'none' });
    });
});
