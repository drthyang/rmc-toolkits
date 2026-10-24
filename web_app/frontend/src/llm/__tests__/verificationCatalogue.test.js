// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The verification's two catalogues (copilot/verification/): CLAIMS, the
// claims that need an analysis, in send-back order; FACTS, the tools' facts
// the answer must not contradict, in check order. The binding test ties
// CLAIMS to the prompt's EVIDENCE_RULES, which stays written by hand
// (prompts/copilot.js): the sentence the rows render equals it byte for byte,
// each row names its first analysis, and what the sentence does not say of
// the rows is an explicit table. Generating the sentence from the catalogue
// is the user's decision: render it, import it in prompts/copilot.js, delete
// the gap table — the model-visible snapshots would show the diff.

import { describe, expect, it } from 'vitest';
import { CLAIMS, createVerification, evidenceSentence, FACTS, promptGaps } from '../copilot/verification';
import { EVIDENCE_RULES } from '../prompts/copilot';

describe('the binding test: EVIDENCE_RULES and the claim catalogue', () => {
    it('the rows render EVIDENCE_RULES byte for byte (the turn\'s claims with every tool offered, and the catalogue itself)', () => {
        expect(evidenceSentence(createVerification().claims)).toBe(EVIDENCE_RULES);
        expect(evidenceSentence(CLAIMS)).toBe(EVIDENCE_RULES);
    });

    it('each row\'s prompt names the first tool that backs it', () => {
        for (const claim of CLAIMS) {
            expect(claim.prompt, claim.id).not.toBeNull();
            expect(claim.prompt.analysis, claim.id).toBe(claim.backedBy.requirements[0].tool);
        }
        expect(createVerification().claims.map((claim) => [claim.id, claim.analyses[0]]))
            .toEqual(CLAIMS.map((claim) => [claim.id, claim.prompt.analysis]));
    });

    // The drift, explicit: what the sentence does not say of each row's
    // backing (read from the code). A row added, an alternative accepted or
    // a requirement changed moves this table, and so does a prompt edit.
    it('what the sentence does not say is the gap table', () => {
        expect(promptGaps(CLAIMS)).toEqual({
            // "distorted polyhedra, polyhedral_distortion": the gate takes analyze_clusters too, and asks for it
            // where the run does not offer polyhedral_distortion (the coordination-shells interview, 2026-10-24).
            polyhedra: ['analyze_clusters'],
            // "a coordination number or split shell, find_gr_peaks": the gate takes compute_bond_angles too.
            coordination: ['compute_bond_angles'],
            // "a lowering …, measure_symmetry_lowering": where it is not offered, the gate asks for detect_symmetry.
            lowering: ['detect_symmetry'],
            // '"no order", every one of them': the gate asks for a per-site directions map and analyze_clusters
            // on each element whose own first shell is a finite cluster's, and for neither polyhedral_distortion,
            // find_gr_peaks nor measure_symmetry_lowering.
            'no-order': [
                'analyze_displacement_directions with site_ref',
                'analyze_clusters with element for each clusterElement',
                'not polyhedral_distortion',
                'not find_gr_peaks',
                'not measure_symmetry_lowering'
            ]
        });
    });
});

describe('the catalogues\' order and shape', () => {
    it('CLAIMS: one row per claim, in send-back order, each with what it reads and what backs it', () => {
        expect(CLAIMS.map((claim) => claim.id)).toEqual(['order-between-cells', 'split-site', 'direction', 'clusters', 'polyhedra', 'coordination', 'lowering', 'no-order']);
        for (const claim of CLAIMS) {
            expect(Object.keys(claim).sort(), claim.id).toEqual(['backedBy', 'id', 'phrasings', 'prompt', 'scope', 'what']);
            expect(['run', 'always'], claim.id).toContain(claim.scope);
            expect(['any', 'first', 'all'], claim.id).toContain(claim.backedBy.need);
            expect(claim.phrasings.asserted.length, claim.id).toBeGreaterThan(0);
            claim.phrasings.asserted.forEach((pattern) => expect(pattern.global, `${claim.id} ${pattern}`).toBe(true));
            expect(Object.isFrozen(claim), claim.id).toBe(true);
        }
        // "No order" is about the run in any sentence; every other claim only in a sentence about the run.
        expect(CLAIMS.filter((claim) => claim.scope === 'always').map((claim) => claim.id)).toEqual(['no-order']);
    });

    it('the send-back order: the first unmet claim goes back first', () => {
        const text = 'Beyond that there is no local order. The Ta4 clusters breathe. Site 2 has a preferred direction. Ti forms a superstructure.';
        const verification = createVerification();
        const sent = [verification.judge(text, { canSendBack: true }).sendBack.id, verification.judge(text, { canSendBack: true }).sendBack.id];
        expect(sent).toEqual(['order-between-cells', 'direction']);
        expect(createVerification().judge(text).done.audit.map((finding) => finding.id)).toEqual(['order-between-cells', 'direction', 'clusters', 'no-order']);
    });

    it('FACTS: one row per facts key, in check order, the producing tools named; a lowering has no check yet', () => {
        expect(FACTS.map((fact) => fact.key)).toEqual(['convergence', 'groups', 'direction', 'modulation', 'cluster', 'coordination', 'distribution', 'lowering']);
        expect(Object.fromEntries(FACTS.map((fact) => [fact.key, fact.from]))).toEqual({
            convergence: ['assess_convergence'],
            groups: ['detect_symmetry'],
            direction: ['analyze_displacement_directions'],
            modulation: ['analyze_modulations'],
            cluster: ['analyze_clusters'],
            coordination: ['find_gr_peaks', 'compute_bond_angles'],
            distribution: ['compute_bond_angles', 'analyze_displacement_directions', 'analyze_clusters'],
            lowering: ['measure_symmetry_lowering']
        });
        expect(FACTS.filter((fact) => !fact.contradicts).map((fact) => fact.key)).toEqual(['lowering']);
        FACTS.forEach((fact) => fact.reads.forEach((pattern) => expect(pattern).toBeInstanceOf(RegExp)));
    });

    it('the turn\'s claims as offered: the rows a tool it offers backs, each with its offered analyses', () => {
        expect(createVerification({ offered: ['compute_bond_angles', 'detect_symmetry'] }).claims).toEqual([
            { id: 'coordination', what: 'a coordination number or a split shell', prompt: { phrase: 'a coordination number or split shell', analysis: 'find_gr_peaks' }, analyses: ['compute_bond_angles'], need: 'any' },
            { id: 'lowering', what: 'a lowering of the average structure', prompt: { phrase: 'a lowering of the average structure', analysis: 'measure_symmetry_lowering' }, analyses: ['detect_symmetry'], need: 'first' }
        ]);
    });
});
