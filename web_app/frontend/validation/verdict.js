// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The verdict every arm ends with, in one fixed machine-readable form, so one
// scoring function judges a Copilot answer, the default pages' reading (arm A)
// and the planted truth (arm D) alike. The live arms are asked to end their
// answer with it (VERDICT_INSTRUCTIONS, appended to the case question); the
// parser takes the last VERDICT line and nothing else: no guessing from prose.

export const FEATURE_KINDS = [
    // Atoms displaced from their site at random (split site); the average site stays.
    'off_centre',
    // Ordered displacements with a wavevector q: a superstructure.
    'modulation',
    // Clusters whose shape is distorted beyond independent atoms.
    'cluster_distortion',
    // Clusters whose size fluctuates collectively (breathing).
    'cluster_breathing',
    // Atoms of a cluster moving together: less distortion than independent atoms.
    'rigid_cluster',
    // A split coordination shell: short and long bonds of one polyhedron.
    'jahn_teller',
    // The average structure itself has lower symmetry than its parent.
    'symmetry_lowering'
];

export const VERDICT_INSTRUCTIONS = 'Finish with one line that starts with VERDICT and holds a JSON object, '
    + 'after everything else (the Outcome line included):\n'
    + 'VERDICT {"features": [{"kind": "...", "element": "...", "direction_uvw": [u, v, w] or null, '
    + '"magnitude_A": number or null, "q": [h, k, l] or null}], "space_group": "..." or null, '
    + '"chiral": true/false/null, "centrosymmetric": true/false/null, "polar": true/false/null}\n'
    + `kind is one of ${FEATURE_KINDS.join(', ')}. features is [] when there is no order beyond thermal motion. `
    + 'direction_uvw is a crystal direction; magnitude_A the static displacement (off_centre, symmetry_lowering) '
    + 'or the long minus short bond (jahn_teller); q in reciprocal-lattice units. Use null for what you did not determine.';

const isVector = (value) => Array.isArray(value) && value.length === 3 && value.every((x) => Number.isFinite(Number(x)));
const bool = (value) => (value === true || value === false ? value : null);

/** A verdict object checked and normalised; throws on a malformed one. */
export const normalizeVerdict = (raw) => {
    if (!raw || typeof raw !== 'object' || Array.isArray(raw)) throw new Error('the verdict is not an object');
    if (!Array.isArray(raw.features)) throw new Error('the verdict has no features array');
    const features = raw.features.map((feature, index) => {
        if (!feature || typeof feature !== 'object') throw new Error(`feature ${index} is not an object`);
        const kind = String(feature.kind || '').trim().toLowerCase().replace(/[\s-]+/g, '_').replace('off_center', 'off_centre');
        if (!FEATURE_KINDS.includes(kind)) throw new Error(`feature ${index}: unknown kind '${feature.kind}'`);
        const out = { kind, element: feature.element == null ? null : String(feature.element).trim() };
        if (isVector(feature.direction_uvw)) out.direction_uvw = feature.direction_uvw.map(Number);
        if (Number.isFinite(Number(feature.magnitude_A)) && feature.magnitude_A !== null) out.magnitude_A = Number(feature.magnitude_A);
        if (isVector(feature.q)) out.q = feature.q.map(Number);
        return out;
    });
    return {
        features,
        space_group: raw.space_group == null || raw.space_group === '' ? null : String(raw.space_group).trim(),
        chiral: bool(raw.chiral),
        centrosymmetric: bool(raw.centrosymmetric),
        polar: bool(raw.polar)
    };
};

// The JSON object starting at `start`: braces counted outside strings.
const jsonObjectAt = (text, start) => {
    let depth = 0;
    let inString = false;
    for (let i = start; i < text.length; i += 1) {
        const char = text[i];
        if (inString) {
            if (char === '\\') i += 1;
            else if (char === '"') inString = false;
        } else if (char === '"') inString = true;
        else if (char === '{') depth += 1;
        else if (char === '}') {
            depth -= 1;
            if (depth === 0) return text.slice(start, i + 1);
        }
    }
    return null;
};

/**
 * The verdict in an answer: the LAST "VERDICT" followed by a JSON object
 * (code fences, bold or a colon around it allowed). `{ verdict, error }`:
 * verdict null when there is none or it does not parse.
 */
// "VERDICT", then bold, a colon or a code fence (```json) before the object.
export const VERDICT_MARK = /VERDICT[\s*:_`]*(?:json)?[\s`]*\{/g;

export const parseVerdict = (text) => {
    const source = String(text || '');
    const marks = [...source.matchAll(VERDICT_MARK)];
    if (!marks.length) return { verdict: null, error: 'no VERDICT line' };
    const mark = marks[marks.length - 1];
    const json = jsonObjectAt(source, mark.index + mark[0].length - 1);
    if (!json) return { verdict: null, error: 'the VERDICT JSON is not closed' };
    try {
        return { verdict: normalizeVerdict(JSON.parse(json)), error: null };
    } catch (failure) {
        return { verdict: null, error: `the VERDICT does not parse: ${failure.message}` };
    }
};

/** The planted truth as a verdict (arm D): what a perfect answer says. */
export const truthVerdict = (truth) => normalizeVerdict({
    features: truth.features.map((feature) => ({
        kind: feature.kind,
        element: feature.element,
        direction_uvw: feature.direction_uvw ?? null,
        magnitude_A: feature.magnitude_A ?? null,
        q: feature.q ?? null
    })),
    space_group: truth.symmetry?.space_group ?? null,
    chiral: truth.symmetry?.chiral ?? null,
    centrosymmetric: truth.symmetry?.centrosymmetric ?? null,
    polar: truth.symmetry?.polar ?? null
});
