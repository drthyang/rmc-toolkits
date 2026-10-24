// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// What follows from a crystallographic point group, so the Copilot never has
// to recall it: whether the group has an inversion centre, whether it holds
// proper rotations only (a Sohncke group: the only ones a chiral structure can
// have), whether it is polar, and its crystal system. The spot check of
// 2026-10-10 called P-42₁m chiral and said F-43m lost its centrosymmetry;
// neither -42m nor -43m has an inversion centre or is a Sohncke group.
//
// The symbols are the symmetry finder's (symmetry.js pointGroupOf): '1', '-1',
// '2', 'm', '2/m', '222', 'mm2', 'mmm', '4', '-4', '4/m', '422', '4mm',
// '-42m', '4/mmm', '3', '-3', '32', '3m', '-3m', '6', '-6', '6/m', '622',
// '6mm', '-6m2', '6/mmm', '23', 'm-3', '432', '-43m', 'm-3m'.

const CENTROSYMMETRIC = new Set(['-1', '2/m', 'mmm', '4/m', '4/mmm', '-3', '-3m', '6/m', '6/mmm', 'm-3', 'm-3m']);
const SOHNCKE = new Set(['1', '2', '222', '4', '422', '3', '32', '6', '622', '23', '432']);
const POLAR = new Set(['1', '2', 'm', 'mm2', '4', '4mm', '3', '3m', '6', '6mm']);
const SYSTEMS = [
    ['triclinic', ['1', '-1']],
    ['monoclinic', ['2', 'm', '2/m']],
    ['orthorhombic', ['222', 'mm2', 'mmm']],
    ['tetragonal', ['4', '-4', '4/m', '422', '4mm', '-42m', '4/mmm']],
    ['trigonal', ['3', '-3', '32', '3m', '-3m']],
    ['hexagonal', ['6', '-6', '6/m', '622', '6mm', '-6m2', '6/mmm']],
    ['cubic', ['23', 'm-3', '432', '-43m', 'm-3m']]
];
const SYSTEM_OF = new Map(SYSTEMS.flatMap(([system, groups]) => groups.map((group) => [group, system])));

/** The point group's facts, or null for a symbol the table does not hold. */
export const pointGroupFacts = (pointGroup) => {
    const key = String(pointGroup ?? '').trim();
    if (!SYSTEM_OF.has(key)) return null;
    return {
        point_group: key,
        centrosymmetric: CENTROSYMMETRIC.has(key),
        sohncke: SOHNCKE.has(key),
        polar: POLAR.has(key),
        crystal_system: SYSTEM_OF.get(key)
    };
};

// --- orders and subgroups (the local-symmetry interview, 2026-10-24) -----------------

// Each point group's order (symmetry.js POINT_GROUP_ORDER; duplicated: this
// module may not import the host).
const ORDER = {
    1: 1, '-1': 2, 2: 2, m: 2, '2/m': 4, 222: 4, mm2: 4, mmm: 8,
    4: 4, '-4': 4, '4/m': 8, 422: 8, '4mm': 8, '-42m': 8, '4/mmm': 16,
    3: 3, '-3': 6, 32: 6, '3m': 6, '-3m': 12, 6: 6, '-6': 6, '6/m': 12,
    622: 12, '6mm': 12, '-6m2': 12, '6/mmm': 24,
    23: 12, 'm-3': 24, 432: 24, '-43m': 24, 'm-3m': 48
};
// The maximal subgroups of each crystallographic point group (ITA Vol. A,
// the subgroup diagram of the 32 point groups): an abstract inclusion, in
// whatever orientation the subgroup takes (-42m in -43m, 3m in -43m along
// <111>). Every subgroup is reached through these.
const MAXIMAL = {
    'm-3m': ['m-3', '432', '-43m', '4/mmm', '-3m'],
    'm-3': ['23', 'mmm', '-3'],
    432: ['23', '422', '32'],
    '-43m': ['23', '-42m', '3m'],
    23: ['222', '3'],
    '6/mmm': ['6/m', '622', '6mm', '-6m2', '-3m', 'mmm'],
    '6/m': ['6', '-6', '-3', '2/m'],
    622: ['6', '32', '222'],
    '6mm': ['6', '3m', 'mm2'],
    '-6m2': ['-6', '32', '3m', 'mm2'],
    6: ['3', '2'],
    '-6': ['3', 'm'],
    '-3m': ['-3', '32', '3m', '2/m'],
    '-3': ['3', '-1'],
    32: ['3', '2'],
    '3m': ['3', 'm'],
    3: ['1'],
    '4/mmm': ['4/m', '422', '4mm', '-42m', 'mmm'],
    '4/m': ['4', '-4', '2/m'],
    422: ['4', '222'],
    '4mm': ['4', 'mm2'],
    '-42m': ['-4', '222', 'mm2'],
    4: ['2'],
    '-4': ['2'],
    mmm: ['222', 'mm2', '2/m'],
    222: ['2'],
    mm2: ['2', 'm'],
    '2/m': ['2', 'm', '-1'],
    2: ['1'],
    m: ['1'],
    '-1': ['1'],
    1: []
};
const SUBGROUPS = new Map(Object.keys(MAXIMAL).map((group) => {
    const all = new Set([group]);
    const visit = (g) => (MAXIMAL[g] || []).forEach((sub) => { if (!all.has(sub)) { all.add(sub); visit(sub); } });
    visit(group);
    return [group, all];
}));

/** A point group's order, or null for a symbol the table does not hold. */
export const pointGroupOrder = (pointGroup) => ORDER[String(pointGroup ?? '').trim()] ?? null;

/** Whether `sub` is (isomorphic to) a subgroup of `parent`, itself included; null when either is unknown. */
export const isPointSubgroup = (sub, parent) => {
    const all = SUBGROUPS.get(String(parent ?? '').trim());
    const key = String(sub ?? '').trim();
    if (!all || !(key in ORDER)) return null;
    return all.has(key);
};

// A symbol's positions as the point group reads them: a screw axis's
// rotation, a glide plane's mirror; a slash belongs to its position (4/m).
const SAME_CLASS = { '-4m2': '-42m', '-62m': '-6m2', '2mm': 'mm2', 'm2m': 'mm2', 312: '32', 321: '32', '31m': '3m', '3m1': '3m', '-31m': '-3m', '-3m1': '-3m' };
const pointGroupOfPositions = (positions) => {
    if (!positions || !positions.length || positions.length > 3) return null;
    let parts = positions;
    // A monoclinic full symbol (P 1 2/m 1) names one position; a trigonal or
    // hexagonal one may hold a 1 (P3m1, P321).
    if (parts.length === 3 && parts.filter((p) => p === '1').length === 2) parts = parts.filter((p) => p !== '1');
    if (parts.length === 3 && /^-?[36]/.test(parts[0])) parts = parts.filter((p, i) => i === 0 || p !== '1');
    const symbol = SAME_CLASS[parts.join('')] ?? parts.join('');
    return symbol in ORDER ? symbol : null;
};

// Every reading of a symbol's axes (what follows its lattice letter) as
// positions: each an optional minus, a rotation digit with its screw and a
// slash part, or a mirror or glide letter. A screw marked by the writer
// (4_2) is read as one; an unmarked digit after a proper rotation n and
// smaller than n may be its screw (42 → 4_2) or the next position (422), and
// both are read (`screws` false: never, as for an R lattice, whose groups
// have no screw axes: R32 is 32). [] when the axes read no way.
const ROTATIONS = new Set(['2', '3', '4', '6']);
const isDigit = (char) => char !== undefined && char >= '0' && char <= '9';
const isPlane = (char) => char !== undefined && 'mnabcde'.includes(char);
const readingsOf = (axes, screws, i = 0) => {
    if (i === axes.length) return [[]];
    if (isPlane(axes[i])) return readingsOf(axes, screws, i + 1).map((rest) => ['m', ...rest]);
    const minus = axes[i] === '-' ? '-' : '';
    const at = i + minus.length;
    const digit = axes[at];
    if (!isDigit(digit)) return [];
    // Where the position's rotation part ends: after a marked screw, or after
    // an unmarked one (a proper rotation only), or right after the digit.
    const ends = [];
    if (axes[at + 1] === '_' && isDigit(axes[at + 2])) ends.push(at + 3);
    else {
        ends.push(at + 1);
        if (screws && !minus && ROTATIONS.has(digit) && isDigit(axes[at + 1]) && Number(axes[at + 1]) < Number(digit)) ends.push(at + 2);
    }
    return ends.flatMap((end) => {
        const slash = axes[end] === '/' && isPlane(axes[end + 1]);
        const position = `${minus}${digit}${slash ? '/m' : ''}`;
        return readingsOf(axes, screws, end + (slash ? 2 : 0)).map((rest) => [position, ...rest]);
    });
};

const SUBSCRIPT_DIGITS = '₀₁₂₃₄₅₆₇₈₉';

/**
 * The point group of a space-group symbol as written (P-42₁m, P4_2/mnm, Fm-3m,
 * P2₁2₁2₁, a TeX overbar, P 1 2/m 1), or null for anything that is not one of
 * the 32 point groups' symbols with a lattice letter in front. Screws the
 * writer marked (4_2, 4₂) are read as such; an unmarked digit after a proper
 * rotation is read both as its screw and as the next position (P42/mnm,
 * P4212, P422, P321), and the symbol's point group is the one its readings
 * agree on (null when none reads as a point group, or two disagree).
 */
export const pointGroupOfSymbol = (symbol) => {
    const text = String(symbol ?? '')
        .replace(/\\bar\{?(\d)\}?/g, '-$1')
        .replace(/(\d)̅/g, '-$1')
        .replace(/[₀-₉]/g, (digit) => `_${SUBSCRIPT_DIGITS.indexOf(digit)}`)
        .replace(/[−–—]/g, '-')
        .replace(/\(no\.?\s*\d+\)/gi, '')
        .replace(/[\s${}\\]/g, '')
        .replace(/_+/g, '_');
    const match = /^([PIFCABR])([-\d_/mnabcde]+)$/.exec(text);
    if (!match) return null;
    const groups = new Set(readingsOf(match[2], match[1] !== 'R').map(pointGroupOfPositions).filter(Boolean));
    return groups.size === 1 ? [...groups][0] : null;
};

/** One line: "-43m: cubic; no inversion centre; achiral (not a Sohncke group); non-polar". */
export const describePointGroup = (pointGroup) => {
    const facts = pointGroupFacts(pointGroup);
    if (!facts) return null;
    return `${facts.point_group}: ${facts.crystal_system}; `
        + `${facts.centrosymmetric ? 'centrosymmetric' : 'no inversion centre'}; `
        + `${facts.sohncke ? 'a Sohncke group (a chiral structure is possible)' : 'achiral (not a Sohncke group)'}; `
        + `${facts.polar ? 'polar' : 'non-polar'}`;
};
