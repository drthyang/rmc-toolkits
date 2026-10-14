// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Scoring a verdict against the planted truth: one function for every arm
// (A, A+, B, C, and the truth itself, which must score 1). Fixed before any
// run; a change here is a change of protocol and is dated in the validation
// doc. Pure, no I/O.

import { normalizeSymbol } from '../src/llm/copilot/answerAudit';

/**
 * The tolerances, each from the physics of the cases:
 *
 * - DIRECTION_DEG 15°: half the smallest angle between two distinct
 *   low-index families (<110> to <111> is 35.3°), so a direction within it
 *   names the right family unambiguously; it also spans more than one cell of
 *   the ν = 10 sphere tiling (~6° across) plus the peak's sampling scatter.
 * - MAGNITUDE_REL 25 % or MAGNITUDE_ABS 0.03 Å, whichever is larger: 0.03 Å
 *   is half the smallest thermal rms of the cases (0.05–0.06 Å), below which
 *   a static shift cannot be told from the thermal width; 25 % accepts
 *   sqrt(rms1² − rms⊥²) and similar estimators (sampling error ~3 % at
 *   500 atoms per site) but fails quoting the thermal-broadened rms1 for a
 *   shift comparable to σ, or a wrong shell.
 * - Q_ABS 0.02 r.l.u.: the supercells are 8 cells long, so allowed q are
 *   multiples of 0.125; 0.02 accepts rounding (0.5 written 0.49) and nothing
 *   else.
 */
export const TOLERANCES = Object.freeze({ DIRECTION_DEG: 15, MAGNITUDE_REL: 0.25, MAGNITUDE_ABS: 0.03, Q_ABS: 0.02 });

const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const norm = (a) => Math.hypot(...a);
const angleDeg = (a, b) => {
    const c = dot(a, b) / (norm(a) * norm(b));
    return (Math.acos(Math.max(-1, Math.min(1, c))) * 180) / Math.PI;
};

// The 48 operations of m-3m on a vector: every permutation with every sign.
const PERMUTATIONS = [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];
const cubicImages = (v) => PERMUTATIONS.flatMap((p) => [-1, 1].flatMap((sx) => [-1, 1].flatMap((sy) => [-1, 1].map((sz) => [
    sx * v[p[0]], sy * v[p[1]], sz * v[p[2]]
]))));

/**
 * The angle (degrees) between a claimed and a planted direction, both as
 * crystal [u v w] of the cubic cells the cases use (so Cartesian-parallel).
 * mode 'axis': ±u are one answer (a split site, a polar axis's domain);
 * 'cubic': any cubic-equivalent of u (a JT long axis disordered over <100>).
 */
export const directionError = (claimed, planted, mode = 'axis') => {
    if (!claimed || norm(claimed) === 0) return null;
    const images = mode === 'cubic' ? cubicImages(planted) : [planted, planted.map((x) => -x)];
    return Math.min(...images.map((image) => angleDeg(claimed, image)));
};

/** Whether two wavevectors are the same order: q ≡ −q and q ≡ q + G. */
export const sameWavevector = (claimed, planted) => {
    if (!claimed) return false;
    const close = (sign) => claimed.every((value, k) => {
        const d = sign * value - planted[k];
        return Math.abs(d - Math.round(d)) <= TOLERANCES.Q_ABS;
    });
    return close(1) || close(-1);
};

export const magnitudeOk = (claimed, planted) => (
    Number.isFinite(claimed) && Math.abs(claimed - planted) <= Math.max(TOLERANCES.MAGNITUDE_REL * planted, TOLERANCES.MAGNITUDE_ABS)
);

const sameElement = (a, b) => typeof a === 'string' && typeof b === 'string' && a.trim().toLowerCase() === b.trim().toLowerCase();

// One claim against one planted feature: what matches and what is accurate.
const compareFeature = (feature, claim) => {
    const out = { kind: claim.kind, element: claim.element };
    const checks = [];
    if (feature.direction_uvw) {
        out.direction_error_deg = directionError(claim.direction_uvw, feature.direction_uvw, feature.direction_mode);
        out.direction_ok = out.direction_error_deg != null && out.direction_error_deg <= TOLERANCES.DIRECTION_DEG;
        checks.push(out.direction_ok);
    }
    if (Number.isFinite(feature.magnitude_A)) {
        out.magnitude_claimed_A = claim.magnitude_A ?? null;
        out.magnitude_ok = magnitudeOk(claim.magnitude_A, feature.magnitude_A);
        checks.push(out.magnitude_ok);
    }
    if (feature.q) {
        out.q_claimed = claim.q ?? null;
        out.q_ok = sameWavevector(claim.q, feature.q);
        checks.push(out.q_ok);
    }
    out.accurate = checks.every(Boolean);
    out.quality = checks.filter(Boolean).length;
    return out;
};

const SYMMETRY_FIELDS = ['space_group', 'chiral', 'centrosymmetric', 'polar'];

const symmetryResult = (truth, verdict) => {
    if (!truth.symmetry) return null;
    const out = {};
    SYMMETRY_FIELDS.forEach((field) => {
        const planted = truth.symmetry[field];
        if (planted == null) return;
        const claimed = verdict?.[field] ?? null;
        if (claimed == null) out[field] = 'unanswered';
        else if (field === 'space_group') out[field] = normalizeSymbol(claimed) === normalizeSymbol(planted) ? 'right' : 'wrong';
        else out[field] = claimed === planted ? 'right' : 'wrong';
    });
    return out;
};

/**
 * The score of one verdict. Components (each 0–1, the case score their mean):
 * every planted feature (1 found and accurate, 0.5 found but a quantity off,
 * 0 missed); "no spurious claim" (1 when every claim is a planted or
 * tolerated feature: for a no-order case this is the over-claim test); and
 * each symmetry field the truth defines (1 right, 0 wrong or unanswered).
 * A missing or malformed verdict scores 0 everywhere.
 */
export const scoreVerdict = (truth, verdict) => {
    const nullCase = truth.features.length === 0;
    if (!verdict) {
        return {
            format_ok: false, null_case: nullCase, features: truth.features.map((f) => ({ kind: f.kind, element: f.element, found: false })),
            spurious: [], tolerated: [], overclaim: false, symmetry: symmetryResult(truth, null), contradictions: [], score: 0
        };
    }
    const claims = verdict.features.map((claim, index) => ({ ...claim, index }));
    const used = new Set();
    const features = truth.features.map((feature) => {
        const candidates = claims
            .filter((claim) => !used.has(claim.index) && feature.accept.includes(claim.kind) && sameElement(claim.element, feature.element))
            .map((claim) => ({ claim, result: compareFeature(feature, claim) }))
            .sort((a, b) => b.result.quality - a.result.quality);
        if (!candidates.length) return { kind: feature.kind, element: feature.element, found: false };
        used.add(candidates[0].claim.index);
        const { quality, ...result } = candidates[0].result;
        return { planted: feature.kind, found: true, claim_index: candidates[0].claim.index, ...result };
    });
    const rest = claims.filter((claim) => !used.has(claim.index));
    const isTolerated = (claim) => truth.tolerated.some((t) => t.kind === claim.kind && sameElement(t.element, claim.element));
    const tolerated = rest.filter(isTolerated).map(({ index, ...claim }) => ({ index, ...claim }));
    const spurious = rest.filter((claim) => !isTolerated(claim)).map(({ index, ...claim }) => ({ index, ...claim }));
    const symmetry = symmetryResult(truth, verdict);

    const components = features.map((f) => (f.found ? (f.accurate ? 1 : 0.5) : 0));
    components.push(spurious.length ? 0 : 1);
    Object.values(symmetry || {}).forEach((value) => components.push(value === 'right' ? 1 : 0));

    // Wrong positive claims, by the answer check that should catch them
    // (answerAudit rules); 'none' = no check covers it.
    const contradictions = [];
    spurious.forEach((claim) => contradictions.push({
        what: `spurious ${claim.kind}${claim.element ? ` (${claim.element})` : ''}`,
        rule: { modulation: 'superstructure', cluster_distortion: 'cluster', cluster_breathing: 'cluster', off_centre: 'direction' }[claim.kind] || 'none'
    }));
    if (symmetry?.chiral === 'wrong' && verdict.chiral === true) contradictions.push({ what: 'called chiral', rule: 'chirality' });
    if (symmetry?.centrosymmetric === 'wrong' && verdict.centrosymmetric === true) contradictions.push({ what: 'called centrosymmetric', rule: 'centrosymmetry' });
    if (symmetry?.polar === 'wrong' && verdict.polar === true) contradictions.push({ what: 'called polar', rule: 'polar' });

    return {
        format_ok: true,
        null_case: nullCase,
        features,
        spurious,
        tolerated,
        overclaim: nullCase && spurious.length > 0,
        symmetry,
        contradictions,
        score: Number((components.reduce((sum, value) => sum + value, 0) / components.length).toFixed(4))
    };
};

/** What must agree across repeated runs for a case to count as reproducible. */
export const verdictSignature = (verdict) => {
    if (!verdict) return 'none';
    const features = verdict.features
        .map((f) => `${f.kind}:${(f.element || '').toLowerCase()}${f.q ? `@${f.q.map((x) => Math.round(Math.abs(x) * 8) / 8).join(',')}` : ''}`)
        .sort();
    const symmetry = ['space_group', 'chiral', 'centrosymmetric', 'polar']
        .map((field) => (field === 'space_group' ? (verdict.space_group ? normalizeSymbol(verdict.space_group) : '-') : String(verdict[field])));
    return `${features.join('|')}#${symmetry.join(',')}`;
};

/**
 * Repeated runs of one case and arm: reproducible when every run gives the
 * same verdict signature; also the share of agreeing pairs and the scores'
 * spread.
 */
export const reproducibility = (runs) => {
    const signatures = runs.map((run) => verdictSignature(run.verdict));
    let pairs = 0;
    let agree = 0;
    for (let i = 0; i < signatures.length; i += 1) {
        for (let j = i + 1; j < signatures.length; j += 1) {
            pairs += 1;
            if (signatures[i] === signatures[j]) agree += 1;
        }
    }
    const scores = runs.map((run) => run.score?.score ?? 0);
    const mean = scores.reduce((sum, value) => sum + value, 0) / Math.max(scores.length, 1);
    const sd = Math.sqrt(scores.reduce((sum, value) => sum + (value - mean) ** 2, 0) / Math.max(scores.length - 1, 1));
    return {
        runs: runs.length,
        identical: runs.length > 1 && new Set(signatures).size === 1,
        pair_agreement: pairs ? Number((agree / pairs).toFixed(3)) : null,
        score_mean: Number(mean.toFixed(4)),
        score_sd: Number(sd.toFixed(4))
    };
};
