// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Distinct, recognizable atom colors for the structure views.
//
// One element -> color map per model, used by every page that draws atoms
// (Atomic Density, Bond Geometry, PCA Ellipsoid, Displacement Directions), so
// an element looks the same everywhere. Every page passes the model's full
// species list (minority species of mixed sites included), and the map
// depends only on that set, never on atom order or on the page.
//
// Common elements use a CPK/Jmol-style table; any element outside it takes an
// unused color from a qualitative palette. Two table colors can still look
// alike (Ti/V, Fe/Co, Mn/Sn, Cr/Mn), so the map is then separated
// automatically: an element whose color sits too close to one already placed
// moves to the nearest color that is far enough from all of them. "Too close"
// is measured in OKLab (a perceptually uniform space) for normal vision and
// again under simulated deuteranopia and protanopia, so a pair that collapses
// for a red–green colour-blind reader is separated too.

export const ELEMENT_COLORS = {
    H: '#cfd3d8', He: '#d9ffff', Li: '#cc80ff', Be: '#c2ff00', B: '#ffb5b5',
    C: '#4b4f57', N: '#3050f8', O: '#e6443b', F: '#5bd35b', Ne: '#b3e3f5',
    Na: '#ab5cf2', Mg: '#69d100', Al: '#bfa6a6', Si: '#f0c8a0', P: '#ff8000',
    S: '#e6c84b', Cl: '#1fd01f', Ar: '#80d1e3', K: '#8f40d4', Ca: '#3dff00',
    Ti: '#bfc2c7', V: '#a6a6ab', Cr: '#8a99c7', Mn: '#9c7ac7', Fe: '#e06633',
    Co: '#f090a0', Ni: '#50d050', Cu: '#c88033', Zn: '#7d80b0', Ga: '#3C5488',
    Ge: '#668f8f', As: '#bd80e3', Se: '#00A087', Br: '#a62929', Rb: '#702eb0',
    Sr: '#43d100', Y: '#94ffff', Zr: '#94e0e0', Nb: '#E64B35', Mo: '#54b5b5',
    Ag: '#c0c0c0', Cd: '#ffd98f', In: '#a67573', Sn: '#668080', Sb: '#9e63b5',
    Te: '#d47a00', I: '#940094', Cs: '#57178f', Ba: '#00c900', Ta: '#F39B7F', W: '#2194d6',
    Pt: '#d0d0e0', Au: '#ffd123', Hg: '#b8b8d0', Pb: '#575961', Bi: '#9e4fb5'
};

export const FALLBACK_PALETTE = [
    '#4DBBD5', '#F39B7F', '#8491B4', '#91D1C2', '#7E6148', '#B09C85',
    '#E6A0C4', '#7AA457', '#C77CFF', '#FFB400', '#00B9E3', '#DC6866',
    '#5B8FF9', '#9FB40F', '#D87C7C', '#6DC8EC'
];

export const DEFAULT_ELEMENT_COLOR = '#8A8F98';

// Minimum OKLab distance between two atoms' colors for normal vision. Pairs
// that were hard to tell apart on the demo sat below it (Ta's old cyan beside
// Se's teal: 0.133); clearly distinct pairs sit well above (Ga/Se 0.23,
// Nb/Se 0.30).
export const MIN_COLOR_DISTANCE = 0.15;
// The same under simulated deuteranopia and protanopia. Lower on purpose: it
// only catches pairs that collapse for a colour-blind reader (Cr/Mn 0.05),
// not every red/green pair, which stays readable by lightness.
export const MIN_CVD_COLOR_DISTANCE = 0.06;

// ── Color math (sRGB ⇄ OKLab, colour-vision simulation) ─────────────────────

const srgbToLinear = (channel) => {
    const c = channel / 255;
    return c <= 0.04045 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
};

const linearToSrgb = (value) => {
    const c = value <= 0.0031308 ? 12.92 * value : 1.055 * value ** (1 / 2.4) - 0.055;
    return Math.round(Math.min(1, Math.max(0, c)) * 255);
};

const hexToLinear = (hex) => {
    const n = parseInt(hex.slice(1), 16);
    return [srgbToLinear((n >> 16) & 255), srgbToLinear((n >> 8) & 255), srgbToLinear(n & 255)];
};

const linearToHex = (rgb) => `#${rgb.map((v) => linearToSrgb(v).toString(16).padStart(2, '0')).join('')}`;

const linearToOklab = ([r, g, b]) => {
    const l = Math.cbrt(0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b);
    const m = Math.cbrt(0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b);
    const s = Math.cbrt(0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b);
    return [
        0.2104542553 * l + 0.793617785 * m - 0.0040720468 * s,
        1.9779984951 * l - 2.428592205 * m + 0.4505937099 * s,
        0.0259040371 * l + 0.7827717662 * m - 0.808675766 * s
    ];
};

const oklabToLinear = ([L, a, b]) => {
    const l = (L + 0.3963377774 * a + 0.2158037573 * b) ** 3;
    const m = (L - 0.1055613458 * a - 0.0638541728 * b) ** 3;
    const s = (L - 0.0894841775 * a - 1.291485548 * b) ** 3;
    return [
        4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s,
        -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s,
        -0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s
    ];
};

// Machado, Oliveira & Fernandes (2009), severity 1.0, applied in linear RGB.
const DEUTERANOPIA = [
    [0.367322, 0.860646, -0.227968],
    [0.280085, 0.672501, 0.047413],
    [-0.01182, 0.04294, 0.968881]
];
const PROTANOPIA = [
    [0.152286, 1.052583, -0.204868],
    [0.114503, 0.786281, 0.099216],
    [-0.003882, -0.048116, 1.051998]
];

const simulate = (matrix, rgb) => matrix.map((row) => (
    Math.min(1, Math.max(0, row[0] * rgb[0] + row[1] * rgb[1] + row[2] * rgb[2]))
));

// A color as the three OKLab points compared: normal, deutan, protan.
const perceive = (rgb) => [
    linearToOklab(rgb),
    linearToOklab(simulate(DEUTERANOPIA, rgb)),
    linearToOklab(simulate(PROTANOPIA, rgb))
];

const distance = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]);

/** OKLab distances between two hex colors: { normal, deuteranopia, protanopia }. */
export const colorDistances = (hexA, hexB) => {
    const a = perceive(hexToLinear(hexA));
    const b = perceive(hexToLinear(hexB));
    return { normal: distance(a[0], b[0]), deuteranopia: distance(a[1], b[1]), protanopia: distance(a[2], b[2]) };
};

// How well a color clears every placed one: the worst of the three margins,
// each distance divided by its threshold (≥ 1 means far enough everywhere).
const clearance = (candidate, placed) => {
    let worst = Infinity;
    placed.forEach((other) => {
        worst = Math.min(
            worst,
            distance(candidate[0], other[0]) / MIN_COLOR_DISTANCE,
            distance(candidate[1], other[1]) / MIN_CVD_COLOR_DISTANCE,
            distance(candidate[2], other[2]) / MIN_CVD_COLOR_DISTANCE
        );
    });
    return worst;
};

// Replacement colors: an OKLCH grid of mid lightness and moderate chroma (no
// near-white or near-black, which vanish on the canvas or read as shadow),
// kept inside the sRGB gamut. Built once, in a fixed order, so every tie
// breaks the same way.
let candidateCache = null;
const candidates = () => {
    if (candidateCache) return candidateCache;
    const list = [];
    [0.5, 0.58, 0.66, 0.74, 0.82].forEach((L) => {
        [0.07, 0.11, 0.15, 0.19].forEach((C) => {
            for (let step = 0; step < 48; step += 1) {
                const h = (step * 7.5 * Math.PI) / 180;
                const rgb = oklabToLinear([L, C * Math.cos(h), C * Math.sin(h)]);
                if (rgb.some((v) => v < -1e-6 || v > 1 + 1e-6)) continue;
                const hex = linearToHex(rgb);
                list.push({ hex, seen: perceive(hexToLinear(hex)) });
            }
        });
    });
    candidateCache = list;
    return list;
};

// ── The per-model map ───────────────────────────────────────────────────────

/**
 * Atom counts per species from a site table (`sites.sites` of /api/pca/sites
 * or the worker): the sum of every site's `elementCounts`, so the minority
 * species of a mixed site count too.
 */
export const speciesCounts = (sites) => {
    const counts = {};
    (sites || []).forEach((site) => {
        const composition = site?.elementCounts ?? (site?.element ? { [site.element]: site.count ?? 0 } : {});
        Object.entries(composition).forEach(([element, count]) => {
            counts[element] = (counts[element] ?? 0) + (Number(count) || 0);
        });
    });
    return counts;
};

// Placement order: the more abundant species first (its share of the atoms,
// rounded to a whole percent, so two payloads that differ by a few skipped
// atoms still agree), then those with a table color, then by name.
const placementOrder = (elements, counts) => {
    const unique = [...new Set(elements)].filter(Boolean).sort();
    const total = unique.reduce((sum, element) => sum + (Number(counts?.[element]) || 0), 0);
    const share = (element) => (total > 0 ? Math.round((100 * (Number(counts?.[element]) || 0)) / total) : 0);
    return unique
        .map((element) => ({ element, share: share(element), table: ELEMENT_COLORS[element] ? 0 : 1 }))
        .sort((a, b) => (b.share - a.share) || (a.table - b.table) || (a.element < b.element ? -1 : a.element > b.element ? 1 : 0));
};

/**
 * The element -> color map for a model, and what was moved.
 *
 * `counts` ({ element: atoms }, optional) sets who keeps a table color when
 * two collide: elements are placed in the order above, and each keeps its
 * table (or palette) color when that color clears every color already placed;
 * otherwise it takes the replacement nearest its own color (in OKLab) that
 * clears them all, or, when none does, the one that clears them best. So the
 * majority species keeps its familiar color and a minority species moves.
 * `adjusted` lists each move as { element, from, to, near } (`near`: the
 * element it was too close to).
 */
export const resolveElementColors = (elements, counts = null) => {
    const ordered = placementOrder(elements, counts);
    const colors = {};
    const adjusted = [];
    const placed = [];       // perceived triples, in placement order
    const placedBy = [];     // element per placed entry
    const used = new Set();
    let cursor = 0;

    ordered.forEach(({ element }) => {
        let base = ELEMENT_COLORS[element];
        if (!base || used.has(base.toLowerCase())) {
            while (cursor < FALLBACK_PALETTE.length && used.has(FALLBACK_PALETTE[cursor].toLowerCase())) {
                cursor += 1;
            }
            base = FALLBACK_PALETTE[cursor] || DEFAULT_ELEMENT_COLOR;
            cursor += 1;
        }
        const seen = perceive(hexToLinear(base));
        let color = base;
        let chosen = seen;
        if (placed.length && clearance(seen, placed) < 1) {
            // The element it is closest to, for the record.
            let near = null;
            let nearest = Infinity;
            placed.forEach((other, index) => {
                const d = distance(seen[0], other[0]);
                if (d < nearest) {
                    nearest = d;
                    near = placedBy[index];
                }
            });
            let best = null;
            let bestScore = Infinity;
            let fallback = null;
            let fallbackClearance = -Infinity;
            candidates().forEach((candidate) => {
                if (used.has(candidate.hex)) return;
                const clear = clearance(candidate.seen, placed);
                if (clear >= 1) {
                    const shift = distance(candidate.seen[0], seen[0]);
                    if (shift < bestScore) {
                        bestScore = shift;
                        best = candidate;
                    }
                } else if (clear > fallbackClearance) {
                    fallbackClearance = clear;
                    fallback = candidate;
                }
            });
            const pick = best || fallback;
            if (pick) {
                color = pick.hex;
                chosen = pick.seen;
                adjusted.push({ element, from: base, to: color, near });
            }
        }
        colors[element] = color;
        placed.push(chosen);
        placedBy.push(element);
        used.add(color.toLowerCase());
    });
    // Keys in name order, so legends list elements alphabetically as before.
    const sorted = {};
    Object.keys(colors).sort().forEach((element) => { sorted[element] = colors[element]; });
    return { colors: sorted, adjusted };
};

const withNotes = (colors, adjusted) => {
    const map = { ...colors };
    const notes = {};
    adjusted.forEach(({ element, from, near }) => {
        notes[element] = `${element}'s usual color ${from} was too close to ${near}'s; moved to tell them apart`;
    });
    // Not enumerable: pages iterate the map as { element: color }.
    Object.defineProperty(map, 'notes', { value: notes, enumerable: false });
    return map;
};

/** Why an element's color differs from its table color (for a legend hover), or undefined. */
export const elementColorNote = (elementColors, element) => elementColors?.notes?.[element];

// Build the element -> color map every page draws with (see
// resolveElementColors). Memoised on the placement order, since every page
// asks for the same model's map.
let lastKey = null;
let lastResult = null;
export const buildElementColors = (elements, counts = null) => {
    const key = placementOrder(elements, counts).map(({ element, share }) => `${element}:${share}`).join(',');
    if (key !== lastKey) {
        lastKey = key;
        lastResult = resolveElementColors(elements, counts);
    }
    return withNotes(lastResult.colors, lastResult.adjusted);
};
