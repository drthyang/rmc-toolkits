// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The holdout scan of split.test.js: what identifies a held-out test case,
// derived at run time from its parameters in cases.js, and a scanner that
// finds it in text with few false positives. No test-case value is written
// here, and a hit is reported with what matched masked, so a failure never
// prints one either.
//
// A case's tokens:
//   id           its id, anywhere;
//   composition  all its elements together (two or more): a formula with exactly
//                that element set ("AbCdX3", with or without subscripts), or, for
//                three or more elements, an element-dash chain of exactly that
//                set ("Ab–Cd–X"). A chain that repeats an element ("Ab–X–Ab") is a
//                bond angle, which every study of that pair writes: not a hit;
//   cluster      its cluster written as a formula of its own ("Ab4" for a
//                tetrahedron, "Ab2" for a dimer), not inside a larger formula;
//   structure    a word of its structure name that no dev case uses and that is
//                not a formula (a structure type);
//   magnitude    a planted length in Å within NEAR_WORDS words of its element
//                (alone, in a dash chain, a cluster or a formula);
//   direction    a [uvw] or <uvw> of the planted family within NEAR_WORDS words
//                of its element;
//   wavevector   an (h, k, l) equivalent to the planted q within NEAR_WORDS words
//                of its element.
// Magnitudes, directions and wavevectors count only where they differ from the
// dev cases of the same kind (otherwise they are the kind's physics, not the
// test case). Thermal widths are never tokens. A bare element symbol is never
// a hit: element symbols are everywhere.

export const NEAR_WORDS = 6;

const ELEMENTS = new Set((
    'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr '
    + 'Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb '
    + 'Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr '
    + 'Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og'
).split(' '));

const SUBSCRIPTS = '₀₁₂₃₄₅₆₇₈₉';
const SUBTOKEN = /[A-Za-z0-9₀-₉\-‐‑‒–—−]+/g;
const FORMULA = /^(?:[A-Z][a-z]?(?:[0-9]+|[₀-₉]+)?)+$/;
const UNIT = /([A-Z][a-z]?)([0-9]+|[₀-₉]+)?/g;

const plainDigits = (text) => text.replace(/[₀-₉]/g, (c) => String(SUBSCRIPTS.indexOf(c)));
const sameSet = (a, b) => a.size === b.size && [...a].every((item) => b.has(item));
const unique = (items) => [...new Set(items)];

// --- tokens from a case ----------------------------------------------------------

const cubicFamily = (uvw) => uvw.map((value) => Math.abs(value)).sort((x, y) => x - y).join(',');
// q ≡ −q ≡ q + G, and a cubic case's q is read up to permutation.
const reducedQ = (q) => q.map((value) => {
    const frac = ((value % 1) + 1) % 1;
    return Math.min(frac, 1 - frac);
}).sort((x, y) => x - y).map((value) => value.toFixed(6)).join(',');
const round = (value, digits) => Number(value.toFixed(digits));

// What each kind plants, read from its parameters: the elements, the cluster
// (if any), the lengths, directions and wavevectors with the elements they
// belong to, and the structure name. Each item has a name, so a test case's
// item is compared with the same item of the dev cases of its kind.
const PLANTED = {
    offcentre: (p) => ({
        elements: p.elements,
        magnitudes: [{ name: 'delta', value: p.delta, near: [p.elements[1]] }, { name: 'a', value: p.a, near: p.elements }],
        directions: [{ name: 'axis', uvw: p.axis, near: [p.elements[1]] }]
    }),
    zone_boundary: (p) => ({
        elements: p.elements,
        magnitudes: [
            { name: 'delta', value: p.delta, near: [p.elements[1]] }, { name: 'a', value: p.a, near: p.elements },
            { name: 'decoy.amp', value: p.decoy.amp, near: [p.elements[0]] }
        ],
        directions: [{ name: 'axis', uvw: p.axis, near: [p.elements[1]] }, { name: 'decoy.axis', uvw: p.decoy.axis, near: [p.elements[0]] }],
        wavevectors: [{ name: 'q', q: p.q, near: [p.elements[1]] }, { name: 'decoy.q', q: p.decoy.q, near: [p.elements[0]] }]
    }),
    breathing: (p) => ({
        elements: [p.element, p.anion],
        cluster: { element: p.element, count: p.shape === 'triangle' ? 3 : 4 },
        magnitudes: [{ name: 'edge', value: p.edge, near: [p.element] }, { name: 'a', value: p.a, near: [p.element, p.anion] }]
    }),
    jahn_teller: (p) => {
        const bond = p.elements.slice(1);
        const split = p.long - p.short;
        return {
            elements: p.elements,
            magnitudes: [
                { name: 'short', value: p.short, near: bond }, { name: 'long', value: p.long, near: bond },
                { name: 'long-short', value: round(split, 2), near: bond }, { name: 'long-short (3)', value: round(split, 3), near: bond },
                { name: 'mean bond', value: round((p.long + 2 * p.short) / 3, 3), near: bond },
                { name: 'a', value: p.a, near: p.elements }
            ]
        };
    },
    rigid_dimer: (p) => ({
        elements: [p.element],
        cluster: { element: p.element, count: 2 },
        magnitudes: [
            { name: 'bond', value: p.bond, near: [p.element] }, { name: 'sigmaCommon', value: p.sigmaCommon, near: [p.element] },
            { name: 'a', value: p.a, near: [p.element] }
        ],
        directions: [{ name: 'axis', uvw: [0, 1, 2].map((k) => (k === p.axis ? 1 : 0)), near: [p.element] }]
    }),
    thermal_null: (p) => ({
        elements: p.elements,
        magnitudes: [{ name: 'a', value: p.a, near: p.elements }]
    }),
    symmetry_lowering: (p) => ({
        elements: p.elements,
        magnitudes: [{ name: 'delta', value: p.delta, near: [p.elements[1]] }, { name: 'a', value: p.a, near: p.elements }],
        directions: [{ name: 'axis', uvw: p.axis, near: [p.elements[1]] }]
    }),
    sites_only: (p) => {
        const elements = unique(p.sites.map((site) => site.el));
        return { elements, magnitudes: [{ name: 'a', value: p.a, near: elements }], structure: p.structure };
    }
};

const planted = (item) => {
    const read = PLANTED[item.kind];
    if (!read) throw new Error(`holdout scan: no tokens for case kind '${item.kind}' (add it to PLANTED)`);
    return { magnitudes: [], directions: [], wavevectors: [], structure: '', ...read(item.params) };
};

const structureWords = (name) => name.split(/\s+/).filter((word) => word && !FORMULA.test(word));

/**
 * The tokens that identify one test case, derived from its parameters. A
 * magnitude, direction or wavevector is kept only when no dev case of the
 * same kind plants the same item with the same value.
 */
export const holdoutTokens = (item, devCases) => {
    const mine = planted(item);
    const devs = devCases.filter((dev) => dev.kind === item.kind).map(planted);
    const devWords = new Set(devCases.flatMap((dev) => structureWords(dev.params.structure ?? '')).map((word) => word.toLowerCase()));
    const sameAsDev = (field, entry, key) => devs.some((dev) => dev[field].some((other) => other.name === entry.name && key(other) === key(entry)));
    const tokens = [{ case: item.id, kind: 'id', name: 'id', text: item.id }];
    if (unique(mine.elements).length >= 2) tokens.push({ case: item.id, kind: 'composition', name: 'elements', elements: new Set(mine.elements) });
    if (mine.cluster) tokens.push({ case: item.id, kind: 'cluster', name: 'cluster', ...mine.cluster });
    structureWords(mine.structure).filter((word) => !devWords.has(word.toLowerCase()))
        .forEach((word) => tokens.push({ case: item.id, kind: 'structure', name: 'structure', word }));
    mine.magnitudes.filter((m) => !sameAsDev('magnitudes', m, (x) => x.value))
        .forEach((m) => tokens.push({ case: item.id, kind: 'magnitude', name: m.name, value: m.value, near: new Set(m.near) }));
    mine.directions.filter((d) => !sameAsDev('directions', d, (x) => cubicFamily(x.uvw)))
        .forEach((d) => tokens.push({ case: item.id, kind: 'direction', name: d.name, family: cubicFamily(d.uvw), near: new Set(d.near) }));
    mine.wavevectors.filter((w) => !sameAsDev('wavevectors', w, (x) => reducedQ(x.q)))
        .forEach((w) => tokens.push({ case: item.id, kind: 'wavevector', name: w.name, q: reducedQ(w.q), near: new Set(w.near) }));
    // Two magnitudes that round to one value are one token.
    return tokens.filter((token, index) => token.kind !== 'magnitude'
        || tokens.findIndex((other) => other.kind === 'magnitude' && other.value === token.value && sameSet(other.near, token.near)) === index);
};

// --- the scanner -----------------------------------------------------------------

const parseFormula = (text) => {
    if (!FORMULA.test(text)) return null;
    const units = [...text.matchAll(UNIT)];
    if (!units.every(([, symbol]) => ELEMENTS.has(symbol))) return null;
    const elements = unique(units.map(([, symbol]) => symbol));
    if (units.length > 1) return { type: 'formula', elements };
    return units[0][2] ? { type: 'cluster', elements, count: Number(plainDigits(units[0][2])) } : { type: 'element', elements };
};

// A sub-token is a run of letters, digits, subscripts and dashes: "Ab–Cd–X"
// is one chain, "Ab4-cluster" a cluster and a word. Offsets are relative.
const readSubtoken = (text) => {
    const parts = [...text.matchAll(/[^-‐‑‒–—−]+/g)];
    if (parts.length >= 2 && parts.every(([part]) => ELEMENTS.has(part))) {
        return [{ type: 'chain', elements: parts.map(([part]) => part), start: 0, end: text.length }];
    }
    return parts.map((part) => {
        const piece = parseFormula(part[0]);
        return piece && { ...piece, start: part.index, end: part.index + part[0].length };
    }).filter(Boolean);
};

const NUMBER_FRACTIONS = { '½': 0.5, '¼': 0.25, '¾': 0.75 };
const readComponent = (text) => {
    const value = text.trim().replace(/^\+/, '').replace(/[−–]/, '-');
    const sign = value.startsWith('-') ? -1 : 1;
    const body = value.replace(/^-/, '');
    if (body in NUMBER_FRACTIONS) return sign * NUMBER_FRACTIONS[body];
    const fraction = /^(\d+)\/(\d+)$/.exec(body);
    if (fraction) return (sign * Number(fraction[1])) / Number(fraction[2]);
    return /^\d+(?:\.\d+)?$/.test(body) ? sign * Number(body) : null;
};

const escape = (text) => text.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
// The number as written, with trailing zeros or not: 0.5 matches "0.50" but
// not "0.51", "10.5" or the "0.5" of "0.5.1".
const numberPattern = (value) => new RegExp(`(?<![\\d.])${escape(String(value))}0*(?![\\d]|\\.\\d)`, 'g');

const lowerBound = (sorted, value) => {
    let lo = 0;
    let hi = sorted.length;
    while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (sorted[mid] < value) lo = mid + 1; else hi = mid;
    }
    return lo;
};
const indexAt = (starts, offset) => lowerBound(starts, offset + 1) - 1;

const DIRECTION = /[[<⟨]\s*([-−]?\d)\s*,?\s*([-−]?\d)\s*,?\s*([-−]?\d)\s*[\]>⟩]/g;
const WAVEVECTOR = /\(\s*([^,()\s]+)\s*,\s*([^,()\s]+)\s*,\s*([^,()\s]+)\s*\)/g;

/** What a text says, read once: its words, lines and the element-bearing sub-tokens. */
export const readText = (text) => {
    // A bracketed vector ("[1 1 0]", "(½, ½, ½)") is one word, however it is spaced.
    const vectors = [...text.matchAll(DIRECTION), ...text.matchAll(WAVEVECTOR)]
        .map((match) => [match.index, match.index + match[0].length]).sort((a, b) => a[0] - b[0]);
    const wordStarts = [];
    let next = 0;
    let reach = -1;
    for (const match of text.matchAll(/\S+/g)) {
        while (next < vectors.length && vectors[next][0] < match.index) reach = Math.max(reach, vectors[next++][1]);
        if (match.index >= reach) wordStarts.push(match.index);
    }
    const lineStarts = [0];
    for (let i = 0; i < text.length; i += 1) if (text[i] === '\n') lineStarts.push(i + 1);
    const pieces = [];
    for (const match of text.matchAll(SUBTOKEN)) {
        readSubtoken(match[0]).forEach((piece) => pieces.push({ ...piece, start: match.index + piece.start, end: match.index + piece.end }));
    }
    const wordsOf = new Map();
    pieces.forEach((piece) => piece.elements.forEach((el) => {
        if (!wordsOf.has(el)) wordsOf.set(el, []);
        wordsOf.get(el).push({ word: indexAt(wordStarts, piece.start), start: piece.start, end: piece.end });
    }));
    return { text, wordStarts, lineStarts, pieces, wordsOf };
};

// The occurrences of the token's elements within NEAR_WORDS words of an offset.
const nearElement = (read, offset, near) => {
    const word = indexAt(read.wordStarts, offset);
    return [...near].flatMap((el) => (read.wordsOf.get(el) ?? []).filter((at) => Math.abs(at.word - word) <= NEAR_WORDS));
};

/** Every hit of the tokens in a read text: { case, kind, name, line, spans }. */
export const scanText = (read, tokens) => {
    const { text } = read;
    const hits = [];
    const hit = (token, spans) => hits.push({
        case: token.case, kind: token.kind, name: token.name, line: indexAt(read.lineStarts, spans[0][0]) + 1, spans
    });
    const nearHits = (token, start, end) => {
        const close = nearElement(read, start, token.near);
        if (close.length) hit(token, [[start, end], ...close.map((at) => [at.start, at.end])]);
    };
    tokens.forEach((token) => {
        if (token.kind === 'id') {
            for (let at = text.indexOf(token.text); at >= 0; at = text.indexOf(token.text, at + 1)) hit(token, [[at, at + token.text.length]]);
        } else if (token.kind === 'composition') {
            const named = (piece) => piece.type === 'formula'
                || (piece.type === 'chain' && token.elements.size >= 3 && piece.elements.length === token.elements.size);
            read.pieces.filter((piece) => named(piece) && sameSet(new Set(piece.elements), token.elements))
                .forEach((piece) => hit(token, [[piece.start, piece.end]]));
        } else if (token.kind === 'cluster') {
            read.pieces.filter((piece) => piece.type === 'cluster' && piece.elements[0] === token.element && piece.count === token.count)
                .forEach((piece) => hit(token, [[piece.start, piece.end]]));
        } else if (token.kind === 'structure') {
            for (const match of text.matchAll(new RegExp(`(?<![A-Za-z])${escape(token.word)}(?![A-Za-z])`, 'gi'))) {
                hit(token, [[match.index, match.index + match[0].length]]);
            }
        } else if (token.kind === 'magnitude') {
            if (!text.includes(String(token.value))) return;
            for (const match of text.matchAll(numberPattern(token.value))) nearHits(token, match.index, match.index + match[0].length);
        } else if (token.kind === 'direction') {
            for (const match of text.matchAll(DIRECTION)) {
                const uvw = match.slice(1, 4).map((value) => Number(value.replace('−', '-')));
                if (cubicFamily(uvw) === token.family) nearHits(token, match.index, match.index + match[0].length);
            }
        } else if (token.kind === 'wavevector') {
            for (const match of text.matchAll(WAVEVECTOR)) {
                const q = match.slice(1, 4).map(readComponent);
                if (q.every((value) => value !== null) && reducedQ(q) === token.q) nearHits(token, match.index, match.index + match[0].length);
            }
        } else {
            throw new Error(`holdout scan: unknown token kind '${token.kind}'`);
        }
    });
    return hits;
};

/**
 * A hit as one line for a failure message: where, which case, which kind of
 * token, and its line with every matched span in the file masked (■), so the
 * message never prints a test-case value.
 */
export const describeHit = (read, hit, allHits, file) => {
    const lineStart = read.lineStarts[hit.line - 1];
    const lineEnd = hit.line < read.lineStarts.length ? read.lineStarts[hit.line] - 1 : read.text.length;
    const masked = new Array(lineEnd - lineStart).fill(false);
    allHits.flatMap((other) => other.spans).forEach(([start, end]) => {
        for (let i = Math.max(start, lineStart); i < Math.min(end, lineEnd); i += 1) masked[i - lineStart] = true;
    });
    const line = read.text.slice(lineStart, lineEnd).split('').map((char, i) => (masked[i] ? '■' : char)).join('')
        .replace(/■+/g, '■').trim();
    return `${file}:${hit.line} ${hit.case} ${hit.kind} (${hit.name}): ${line.length > 140 ? `${line.slice(0, 137)}...` : line}`;
};
