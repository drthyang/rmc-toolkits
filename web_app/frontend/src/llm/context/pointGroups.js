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

/** One line: "-43m: cubic; no inversion centre; achiral (not a Sohncke group); non-polar". */
export const describePointGroup = (pointGroup) => {
    const facts = pointGroupFacts(pointGroup);
    if (!facts) return null;
    return `${facts.point_group}: ${facts.crystal_system}; `
        + `${facts.centrosymmetric ? 'centrosymmetric' : 'no inversion centre'}; `
        + `${facts.sohncke ? 'a Sohncke group (a chiral structure is possible)' : 'achiral (not a Sohncke group)'}; `
        + `${facts.polar ? 'polar' : 'non-polar'}`;
};
