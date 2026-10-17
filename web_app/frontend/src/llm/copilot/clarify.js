// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Questions the app asks itself. In the 2026-10-11 checks gemma4:26b never
// called ask_user: on a request that leaves a choice open it picked one and
// moved on. Where the run data show that the choice changes the result, the
// Workbench asks the user before the model guesses, with the options it
// computed and a default. It asks only then:
//
//   (a) a bond-angle or coordination request naming only an element: which
//       neighbours (one option per partner element with a first shell);
//   (a) a request naming a pair or triplet whose first shell is split into
//       peaks, or not separated from the next, when the window changes the
//       coordination number by ≥ max(0.5, 15 %): the whole shell (default), or
//       the first peak, or a cut 0.2 Å either side of the minimum;
//   (b) an analysis about to run in a window that cuts a shell (it holds less
//       than the whole shell by the same margin): the whole shell (default),
//       or the window as asked;
//   (c) displacement directions of an element with several sites whose site
//       symmetry is not cubic: pooling mixes their local frames and can cancel
//       a per-site preference (the review's pooled Ta map), so one site of an
//       orbit (default) or all sites pooled;
//   (c) the space group with no tolerance named, when more than one group
//       holds above the sampling noise floor: which tolerance.
//
// A question that settles the choice itself asks nothing: an ideal polyhedron
// it names (an octahedron: 6) that exactly one option's coordination number
// matches, the shell it names ("only the first Ta–Se peak", "the whole first
// shell", "both peaks"), an element pair or triplet, a site number, a window
// or a tolerance in Å. A window question is asked once per turn: the user's
// choice for a (centre, partner, window) is remembered in the turn's
// `session` and applied to a later call with that window without asking
// again. Pure functions of the run data; the agent loop shows the card
// (agentLoop.js) and hands the model the user's choice.

import { elementDensity, extractPeaks, firstShell } from '../context/pairCorrelations';

const roundTo = (value, digits = 2) => (Number.isFinite(value) ? Number(value.toFixed(digits)) + 0 : null);
const cnText = (value) => (Number.isFinite(value) ? (Math.abs(value - Math.round(value)) < 0.05 ? String(Math.round(value)) : value.toFixed(1)) : '?');
const windowText = ([lo, hi]) => `${roundTo(lo)}–${roundTo(hi)} Å`;

// A window changes the answer when its coordination number differs from the
// whole shell's by at least this much (absolute, or relative to the shell's).
export const CN_CHANGE_ABS = 0.5;
export const CN_CHANGE_REL = 0.15;
// A neighbour shell counts as bonded when it starts within this factor of the
// nearest shell's start.
export const BONDED_REACH = 1.4;
// A ladder rung above the sampling noise floor must reach this factor of it.
export const NOISE_MARGIN = 1.25;
const changes = (cn, reference) => Number.isFinite(cn) && Number.isFinite(reference)
    && Math.abs(cn - reference) >= Math.max(CN_CHANGE_ABS, CN_CHANGE_REL * reference);

const trapezoidRange = (x, f, lo, hi) => {
    let sum = 0;
    for (let i = 0; i < x.length - 1; i += 1) {
        const a = Math.max(lo, x[i]);
        const b = Math.min(hi, x[i + 1]);
        if (b <= a) continue;
        const slope = (f[i + 1] - f[i]) / (x[i + 1] - x[i]);
        const fa = f[i] + slope * (a - x[i]);
        const fb = f[i] + slope * (b - x[i]);
        sum += 0.5 * (fa + fb) * (b - a);
    }
    return sum;
};

// --- the shells the run's partial g(r) hold ------------------------------------

/**
 * The first coordination shell of every element pair in the run's partial-PDF
 * file: `shell(centre, partner)` → `{ window, peaks, cn, firstPeak: { window,
 * cn }, cut: { earlier, later }, troughRatio }` (CN = partner atoms per centre
 * atom), `cnIn(centre, partner, [lo, hi])`, `partners(centre)`, nearest first.
 * Null without partials or densities.
 */
export const buildShellCatalog = ({ plotFiles, structure } = {}) => {
    const partials = (plotFiles || []).find((file) => file.plotKind === 'pdf_partials' && file.plotData?.series?.length);
    if (!partials || !structure?.elementCounts) return null;
    const series = new Map();
    for (const entry of partials.plotData.series) {
        const match = /^([A-Z][a-z]?)\s*-\s*([A-Z][a-z]?)$/.exec((entry.label || '').trim());
        if (!match || !Array.isArray(entry.x) || entry.x.length !== entry.y?.length) continue;
        const radial = entry.y.map((value, i) => Math.max(0, value) * entry.x[i] * entry.x[i]);
        const value = { x: entry.x, y: entry.y, radial };
        series.set(`${match[1]}|${match[2]}`, value);
        series.set(`${match[2]}|${match[1]}`, value);
    }
    if (!series.size) return null;
    const cache = new Map();
    const cnIn = (centre, partner, [lo, hi]) => {
        const data = series.get(`${centre}|${partner}`);
        const density = elementDensity(structure, partner);
        if (!data || !density || !(hi > lo)) return null;
        return 4 * Math.PI * density * trapezoidRange(data.x, data.radial, lo, hi);
    };
    const shell = (centre, partner) => {
        const key = `${centre}|${partner}`;
        if (cache.has(key)) return cache.get(key);
        const data = series.get(key);
        let result = null;
        if (data) {
            const peaks = extractPeaks(data.x, data.y, { rMax: 8, maxPeaks: 3 });
            // With both densities a split shell cut at a deep trough is joined
            // when its summed CN is an ideal polyhedral count (firstShell).
            const found = peaks.length ? firstShell(data.x, data.y, {
                rMax: Math.max(6, peaks[0].r + 3),
                partnerDensity: elementDensity(structure, partner),
                centreDensity: elementDensity(structure, centre)
            }) : null;
            if (found?.boundaryR != null) {
                const step = data.x.length > 1 ? Math.abs(data.x[1] - data.x[0]) : 0;
                // Rounded as shown to the user (0.01 Å), so a window applied
                // from a card is the one its option names.
                const window = [roundTo(Math.max(0, found.onsetR - step)), roundTo(found.boundaryR)];
                const first = peaks[0];
                // A joined split shell: its first peak alone ends at the deep
                // trough the shell was cut at before.
                const firstWindow = found.split
                    ? [window[0], roundTo(found.split.firstBoundaryR)]
                    : [roundTo(Math.max(0, first.r - first.fwhm)), roundTo(first.r + first.fwhm)];
                result = {
                    centre,
                    partner,
                    window,
                    peaks: found.peaksInShell,
                    cn: cnIn(centre, partner, window),
                    firstPeak: { window: firstWindow, cn: cnIn(centre, partner, firstWindow) },
                    cut: {
                        earlier: { window: [window[0], roundTo(window[1] - 0.2)], cn: cnIn(centre, partner, [window[0], roundTo(window[1] - 0.2)]) },
                        later: { window: [window[0], roundTo(window[1] + 0.2)], cn: cnIn(centre, partner, [window[0], roundTo(window[1] + 0.2)]) }
                    },
                    troughRatio: found.troughRatio
                };
            }
        }
        cache.set(key, result);
        return result;
    };
    const elements = [...new Set([...series.keys()].flatMap((key) => key.split('|')))].sort();
    // The bonded neighbours: shells starting within BONDED_REACH of the
    // nearest one (the demo's Ta: Se from 2.21 Å, Ta from 2.65 Å — the Ta4
    // cluster — but not Ga from 3.77 Å).
    const partners = (centre) => {
        const shells = elements
            .map((partner) => shell(centre, partner))
            .filter((item) => item && Number.isFinite(item.cn) && item.cn >= 0.5)
            .sort((a, b) => a.window[0] - b.window[0]);
        const nearest = shells[0]?.window[0];
        return shells.filter((item) => item.window[0] <= BONDED_REACH * nearest);
    };
    return { elements, shell, cnIn, partners };
};

/**
 * The windows the data offer for one shell when the choice changes the
 * answer: [{ label, description, window, cn }], the whole shell first; null
 * when the shell is clear.
 */
export const ambiguousWindows = (shell) => {
    if (!shell || !Number.isFinite(shell.cn)) return null;
    const pair = `${shell.partner} per ${shell.centre}`;
    const whole = {
        label: 'Whole first shell',
        description: `${windowText(shell.window)}${shell.peaks.length > 1 ? `, both peaks (${shell.peaks.map((r) => roundTo(r)).join(' and ')} Å)` : ''}: CN ${cnText(shell.cn)} ${pair}`,
        window: shell.window,
        cn: shell.cn
    };
    if (shell.peaks.length > 1 && changes(shell.firstPeak.cn, shell.cn)) {
        return [whole, {
            label: 'First peak only',
            description: `${windowText(shell.firstPeak.window)}, the shortest bonds: CN ${cnText(shell.firstPeak.cn)} ${pair}`,
            window: shell.firstPeak.window,
            cn: shell.firstPeak.cn
        }];
    }
    const { earlier, later } = shell.cut;
    if (shell.troughRatio > 0.02 && (changes(earlier.cn, shell.cn) || changes(later.cn, shell.cn))) {
        return [
            { ...whole, label: 'Cut at the minimum', description: `${windowText(shell.window)}: CN ${cnText(shell.cn)} ${pair}` },
            { label: 'Cut 0.2 Å earlier', description: `${windowText(earlier.window)}: CN ${cnText(earlier.cn)} ${pair}`, window: earlier.window, cn: earlier.cn },
            { label: 'Cut 0.2 Å later', description: `${windowText(later.window)}: CN ${cnText(later.cn)} ${pair}`, window: later.window, cn: later.cn }
        ];
    }
    return null;
};

// --- what a question names ------------------------------------------------------

const ANGLES = /\b(bond[\s-]?angles?|angles?|triplets?)\b/i;
const COORDINATION = /\b(coordination|coordinated|polyhedr\w*|octahedr\w*|tetrahedr\w*|neighbou?rs?|bond\s+counts?)\b|\bCN\b/i;
const DIRECTIONS = /\b(preferred|preferential)\b|\b(displacement|displacements)\s+(direction|orientation)s?\b|\bdirection(s|al)?\b/i;
const SYMMETRY = /\bspace[\s-]?groups?\b/i;
const IN_ANGSTROM = /\d+(?:\.\d+)?\s*(?:Å|A\b|ångström|angstrom)/i;
const SITE_NUMBER = /\bsites?\s*(?:#|no\.?|ref\.?|number)?\s*\d+\b|\bref\s*\d+\b/i;
// A shell the question itself names: the first peak only, or the whole shell.
const FIRST_PEAK = /\bfirst\s+(?:[A-Z][a-z]?\s*[–—−-]\s*[A-Z][a-z]?\s+)?(?:g\(r\)\s+)?peak\b/i;
const WHOLE_SHELL = /\b(?:whole|entire|full|complete)\s+(?:first\s+)?(?:coordination\s+)?shell\b|\bboth\s+(?:[A-Z][a-z]?\s*[–—−-]\s*[A-Z][a-z]?\s+)?peaks\b/i;
const POLYHEDRA = [[/\boctahedr/i, 6], [/\btetrahedr/i, 4], [/\btrigonal\s+prism/i, 6], [/\bsquare[\s-]planar/i, 4], [/\bcuboctahedr/i, 12]];
// Element symbols that are also English words: named only in an element-like
// spot ("around In", "In atoms", "In–Se").
const WORDLIKE = new Set(['I', 'In', 'As', 'At', 'Be', 'He', 'No', 'S', 'B', 'C', 'O', 'P', 'K', 'U', 'V', 'W', 'Y', 'Am', 'Es', 'Ho', 'Pa']);

const DASH = '[–—−-]';

/** What a question names: its intent and the run's elements, pairs and triplets in it. */
export const readQuestion = (text, elements) => {
    const source = String(text || '');
    const known = new Set(elements || []);
    const triplets = [];
    const pairs = [];
    const tripletPattern = new RegExp(`\\b([A-Z][a-z]?)\\s*${DASH}\\s*([A-Z][a-z]?)\\s*${DASH}\\s*([A-Z][a-z]?)\\b`, 'g');
    let rest = source;
    for (const match of source.matchAll(tripletPattern)) {
        if ([match[1], match[2], match[3]].every((el) => known.has(el))) {
            triplets.push([match[1], match[2], match[3]]);
            rest = rest.replace(match[0], ' ');
        }
    }
    const pairPattern = new RegExp(`\\b([A-Z][a-z]?)\\s*${DASH}\\s*([A-Z][a-z]?)\\b`, 'g');
    for (const match of rest.matchAll(pairPattern)) {
        if (known.has(match[1]) && known.has(match[2])) {
            pairs.push([match[1], match[2]]);
            rest = rest.replace(match[0], ' ');
        }
    }
    const singles = [];
    for (const element of known) {
        const pattern = new RegExp(`(^|[^A-Za-z])${element}(?![A-Za-z])`, 'g');
        for (const match of rest.matchAll(pattern)) {
            const at = match.index + match[1].length;
            if (WORDLIKE.has(element)) {
                const before = rest.slice(Math.max(0, at - 14), at);
                const after = rest.slice(at + element.length, at + element.length + 8);
                if (!/(around|of|for|on|about|central|centred on|centered on)\s*$/i.test(before) && !/^\s*(atoms?|sites?)\b/i.test(after)) continue;
            }
            if (!singles.includes(element)) singles.push(element);
        }
    }
    const polyhedron = POLYHEDRA.find(([pattern]) => pattern.test(source));
    return {
        angles: ANGLES.test(source),
        coordination: COORDINATION.test(source),
        directions: DIRECTIONS.test(source),
        symmetry: SYMMETRY.test(source),
        namesAngstrom: IN_ANGSTROM.test(source),
        namesSite: SITE_NUMBER.test(source),
        polyhedronCn: polyhedron ? polyhedron[1] : null,
        // 'first-peak', 'whole' or null: a window choice the question makes.
        shellChoice: FIRST_PEAK.test(source) ? 'first-peak' : (WHOLE_SHELL.test(source) ? 'whole' : null),
        triplets,
        pairs,
        elements: singles
    };
};

// --- the questions --------------------------------------------------------------

const MAX_OPTIONS = 4;
const CUBIC_SITE = new Set(['23', 'm-3', '432', '-43m', 'm-3m']);

const choice = (label, description, apply, isDefault = false) => ({ label, description, apply, isDefault });

const asRequest = (question, choices, kind) => ({
    kind,
    request: {
        question,
        options: choices.map(({ label, description }) => ({ label, ...(description ? { description } : {}) })),
        allowOther: true,
        defaultIndex: Math.max(0, choices.findIndex((item) => item.isDefault))
    },
    choices
});

// The polyhedron the question names settles a window when exactly one option
// holds its coordination number.
const settledByPolyhedron = (options, polyhedronCn) => {
    if (!polyhedronCn || !options) return null;
    const matching = options.filter((option) => Math.abs(option.cn - polyhedronCn) <= 0.5);
    return matching.length === 1 ? matching[0] : null;
};

// (a) Which neighbours, for an angle or coordination request naming one element.
const neighbourQuestion = (intent, catalog) => {
    const [centre] = intent.elements;
    const partners = catalog.partners(centre);
    if (partners.length < 2) return null;
    if (intent.polyhedronCn) {
        const match = partners.filter((shell) => Math.abs(shell.cn - intent.polyhedronCn) <= 0.5);
        if (match.length === 1) return null;
    }
    const shown = partners.slice(0, intent.angles ? MAX_OPTIONS : MAX_OPTIONS - 1);
    const choices = shown.map((shell, index) => {
        const label = intent.angles ? `${shell.partner}–${centre}–${shell.partner}` : `${shell.partner} around ${centre}`;
        const description = `${shell.partner} in the first ${centre}–${shell.partner} shell, ${windowText(shell.window)}: `
            + `CN ${cnText(shell.cn)} ${shell.partner} per ${centre}`;
        return choice(label, description,
            () => `the user chose ${label}: ${shell.partner} around ${centre} in the whole first ${centre}–${shell.partner} shell, `
                + `${windowText(shell.window)} (CN ${cnText(shell.cn)}). Use exactly this.`,
            index === 0);
    });
    if (!intent.angles) {
        choices.push(choice('Every neighbour type', shown.map((shell) => `${shell.partner} ${cnText(shell.cn)}`).join(', '),
            () => `the user chose every neighbour type of ${centre}: ${shown.map((shell) => `${shell.partner} in ${windowText(shell.window)} `
                + `(CN ${cnText(shell.cn)})`).join('; ')}. Report each.`));
    }
    return asRequest(
        intent.angles ? `Which bond angles around ${centre}?` : `The coordination of ${centre} by which neighbours?`,
        choices,
        'neighbours'
    );
};

// (a) Which window, for a named pair or triplet whose first shell is ambiguous.
const windowQuestion = (intent, catalog) => {
    if (intent.shellChoice) return null;
    const candidates = [];
    for (const [end1, apex, end2] of intent.triplets) {
        candidates.push([apex, end1]);
        if (end2 !== end1) candidates.push([apex, end2]);
    }
    for (const [a, b] of intent.pairs) candidates.push([a, b]);
    for (const [centre, partner] of candidates) {
        const shell = catalog.shell(centre, partner);
        const options = ambiguousWindows(shell);
        if (!options) continue;
        if (settledByPolyhedron(options, intent.polyhedronCn)) return null;
        const choices = options.map((option, index) => choice(option.label, option.description,
            () => `the user chose "${option.label}" for the ${centre}–${partner} bonds: the window ${windowText(option.window)} `
                + `(CN ${cnText(option.cn)} ${partner} per ${centre}). Use exactly this window.`,
            index === 0));
        return asRequest(`Which ${centre}–${partner} window should the analysis use?`, choices, 'window');
    }
    return null;
};

// (c) One site or all, for the displacement directions of a multi-site element.
const directionQuestion = (intent, run) => {
    const [element] = intent.elements;
    const orbits = (run.symmetry?.orbits || []).filter((orbit) => orbit.element === element);
    const sites = (run.structure?.basis || []).filter((site) => site.el === element);
    if (sites.length < 2 || !orbits.length) return null;
    const anisotropic = orbits.filter((orbit) => !CUBIC_SITE.has(String(orbit.site || '').replace(/\./g, '')));
    if (!anisotropic.length) return null;
    const representative = (orbit) => run.structure.basis[orbit.members?.[0]]?.referenceNumber ?? null;
    const choices = orbits.slice(0, MAX_OPTIONS - 1).map((orbit, index) => {
        const ref = representative(orbit);
        const name = `${orbit.wyckoffMultiplicity ?? orbit.size}${orbit.wyckoff || ''}`;
        return choice(
            orbits.length > 1 ? `One site of the ${name} orbit` : `One ${element} site`,
            `site ${ref ?? '?'}; its ${orbit.size} sites (site symmetry ${orbit.site}) are equivalent by symmetry`,
            () => `the user chose one ${element} site: analyse site_ref ${ref} (the ${name} orbit; its sites are equivalent by symmetry). `
                + 'Do not pool the element.',
            index === 0
        );
    });
    choices.push(choice(`All ${sites.length} ${element} sites pooled`,
        'their local axes differ, so a direction each site prefers can cancel',
        () => `the user chose all ${element} sites pooled (element ${element}); say that pooling mixes the sites' local frames.`));
    return asRequest(`Displacement directions of ${element}: one site, or all ${sites.length} pooled?`, choices, 'sites');
};

// (c) Which tolerance, for the space group when several groups hold above the noise floor.
const toleranceQuestion = (run, noiseFloor) => {
    const ladder = run.symmetry?.ladder;
    if (!Array.isArray(ladder) || ladder.length < 2) return null;
    const floor = Number.isFinite(noiseFloor) ? noiseFloor : 0;
    // A group counts when it holds clearly above the floor, not only up to it.
    const above = ladder.filter((brick) => brick.to > NOISE_MARGIN * floor && brick.spaceGroup);
    if (above.length < 2) return null;
    const atDefault = above.find((brick) => brick.from <= 0.2 && 0.2 < brick.to) || above.at(-1);
    const choices = above.slice(-MAX_OPTIONS).reverse().map((brick) => choice(
        brick.spaceGroup,
        `holds from ${roundTo(Math.max(brick.from, floor), 3)} to ${roundTo(brick.to, 3)} Å`,
        () => `the user chose ${brick.spaceGroup}, which holds from ${roundTo(Math.max(brick.from, floor), 3)} to ${roundTo(brick.to, 3)} Å; `
            + 'report it with that tolerance range.',
        brick === atDefault
    ));
    return asRequest('Which tolerance? Several groups hold above the sampling noise floor.', choices, 'tolerance');
};

/**
 * The Workbench's clarifier for one run: `beforeTurn(questionText)` and
 * `beforeCall(tool, args, { questionText })`, each null or
 * `{ kind, request, choices }` (the card's request; each choice's `apply`
 * gives the note for the model, or the new arguments). `run` holds the host
 * props (plotFiles, structure, symmetry); `noiseFloor` the context's.
 */
export const buildClarifier = (run = {}, { noiseFloor = null } = {}) => {
    const catalog = buildShellCatalog(run);
    const elements = Object.keys(run.structure?.elementCounts || {});
    return {
        catalog,
        beforeTurn: (questionText) => {
            const intent = readQuestion(questionText, elements);
            if ((intent.angles || intent.coordination) && catalog && !intent.namesAngstrom) {
                if (!intent.triplets.length && !intent.pairs.length && intent.elements.length === 1) {
                    return neighbourQuestion(intent, catalog);
                }
                if (intent.triplets.length || intent.pairs.length) return windowQuestion(intent, catalog);
            }
            if (intent.directions && !intent.angles && !intent.namesSite && intent.elements.length === 1
                && !intent.pairs.length && !intent.triplets.length) {
                return directionQuestion(intent, run);
            }
            if (intent.symmetry && !intent.namesAngstrom) return toleranceQuestion(run, noiseFloor);
            return null;
        },
        beforeCall: (tool, args, { questionText = '', session = null } = {}) => {
            if (!catalog || typeof tool?.shellWindows !== 'function') return null;
            let windows;
            try {
                windows = tool.shellWindows(args) || [];
            } catch {
                return null;
            }
            const shellChoice = readQuestion(questionText, elements).shellChoice;
            for (const [windowIndex, item] of windows.entries()) {
                const [lo, hi] = item.window || [];
                if (!(hi > lo)) continue;
                // A window the user wrote into the question is theirs.
                const named = [lo, hi].every((value) => new RegExp(`(^|[^\\d.])${String(value).replace('.', '\\.')}0*(?![\\d])`).test(questionText));
                if (named) continue;
                const shell = catalog.shell(item.centre, item.partner);
                if (!shell || !Number.isFinite(shell.cn)) continue;
                const cn = catalog.cnIn(item.centre, item.partner, [lo, hi]);
                const overlaps = hi > shell.window[0] && lo < shell.window[1];
                if (!overlaps || !(cn < shell.cn) || !changes(cn, shell.cn)) continue;
                const pair = `${item.partner} per ${item.centre}`;
                // Answered already in this turn: the same choice, unasked.
                const key = windowKey(item.centre, item.partner, [lo, hi]);
                const remembered = session?.windows?.get(key);
                if (remembered) {
                    return settled(item.apply(remembered.window, { confirmed: remembered.confirmed }),
                        `the user already chose ${windowText(remembered.window)} for the ${item.centre}–${item.partner} window `
                        + `${windowText([lo, hi])} in this question: run in it.`);
                }
                // The question names the shell: the first peak only keeps a
                // partial window (as the user's), the whole shell replaces it.
                if (shellChoice === 'first-peak') {
                    return settled(item.apply([lo, hi], { confirmed: true }),
                        `the question asks for the first peak only, so ${windowText([lo, hi])} (CN ${cnText(cn)} ${pair}) is kept as the user's window.`);
                }
                if (shellChoice === 'whole') {
                    return settled(item.apply(shell.window, { confirmed: false }),
                        `the question asks for the whole first shell, so the run uses ${windowText(shell.window)} (CN ${cnText(shell.cn)} ${pair}) `
                        + `instead of ${windowText([lo, hi])}.`);
                }
                const asked = asRequest(
                    `The window ${windowText([lo, hi])} cuts the first ${item.centre}–${item.partner} shell: it holds ${cnText(cn)} ${pair} `
                        + `of the shell's ${cnText(shell.cn)}. Which window?`,
                    [
                        choice('Whole first shell', `${windowText(shell.window)}: CN ${cnText(shell.cn)} ${pair}`,
                            () => ({ args: item.apply(shell.window, { confirmed: false }), window: shell.window }), true),
                        choice(`Keep ${windowText([lo, hi])}`, `part of the shell: CN ${cnText(cn)} ${pair}`,
                            () => ({ args: item.apply([lo, hi], { confirmed: true }), window: [lo, hi] }))
                    ],
                    'cut'
                );
                return { ...asked, key, windowIndex };
            }
            return null;
        }
    };
};

// A window settled without a question: the arguments to run with, and why.
const settled = (args, note) => ({ kind: 'cut', settled: true, args, note });
const windowKey = (centre, partner, [lo, hi]) => `${centre}|${partner}|${roundTo(lo)}|${roundTo(hi)}`;

// Two numbers in a typed answer: a window.
const typedWindow = (text) => {
    const match = /(\d+(?:\.\d+)?)\s*(?:Å\s*)?(?:[–—−-]|to|and)\s*(\d+(?:\.\d+)?)/i.exec(String(text || ''));
    if (!match) return null;
    const lo = Number(match[1]);
    const hi = Number(match[2]);
    return hi > lo ? [lo, hi] : null;
};

/**
 * What a reply means: for a question before the turn, the note the model
 * reads with the user's question; for a cut window, the arguments to run with
 * and the note, and the choice is remembered in the turn's `session` (a later
 * call with the same window gets it unasked). `item` is the
 * beforeCall/beforeTurn result.
 */
export const resolveClarification = (item, reply, { args = null, tool = null, session = null } = {}) => {
    const remember = (window, confirmed) => {
        if (item.kind === 'cut' && item.key && session?.windows) session.windows.set(item.key, { window, confirmed });
    };
    if (reply?.source === 'other') {
        const answer = String(reply.answer || '').trim();
        if (item.kind === 'cut' && tool) {
            const window = typedWindow(answer);
            // The window the question was about (r12 or r23), not always the first.
            const target = (tool.shellWindows(args) || [])[item.windowIndex ?? 0];
            if (window && target) {
                remember(window, true);
                return { args: target.apply(window, { confirmed: true }), note: `the user typed "${answer}": run in ${windowText(window)}.` };
            }
            return { args, note: `the user answered "${answer}"; the window was kept as asked.` };
        }
        return { note: `asked, the user answered in their own words: "${answer}". Follow it exactly.` };
    }
    const index = reply?.source === 'skipped' ? item.request.defaultIndex : reply?.index;
    const picked = item.choices[Number.isInteger(index) ? index : item.request.defaultIndex] || item.choices[0];
    const result = picked.apply();
    const how = reply?.source === 'skipped' ? 'the user skipped, so the default stands' : 'the user chose it';
    if (item.kind === 'cut') {
        remember(result.window, Boolean(result.args?.window_confirmed) || picked.label.startsWith('Keep'));
        return { args: result.args, note: `${how}: run in ${windowText(result.window)} (${picked.label}).` };
    }
    return { note: reply?.source === 'skipped' ? `${result} (The user skipped the question; this is the default.)` : result };
};
