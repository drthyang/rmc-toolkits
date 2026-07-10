// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Distinct, recognizable atom colors for the structure views.
//
// Common elements use a CPK/Jmol-style table; any element outside it is assigned
// an unused color from a qualitative palette, so that no two elements in a
// structure ever share a color. Shared by the slab-in-cell canvas and the
// folded-unit-cell 3D view so an element looks the same in both.

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
    Te: '#d47a00', I: '#940094', Cs: '#57178f', Ba: '#00c900', W: '#2194d6',
    Pt: '#d0d0e0', Au: '#ffd123', Hg: '#b8b8d0', Pb: '#575961', Bi: '#9e4fb5'
};

export const FALLBACK_PALETTE = [
    '#4DBBD5', '#F39B7F', '#8491B4', '#91D1C2', '#7E6148', '#B09C85',
    '#E6A0C4', '#7AA457', '#C77CFF', '#FFB400', '#00B9E3', '#DC6866',
    '#5B8FF9', '#9FB40F', '#D87C7C', '#6DC8EC'
];

export const DEFAULT_ELEMENT_COLOR = '#8A8F98';

// Build a stable element -> color map for the structure. Each element keeps its
// CPK-style color when known; the rest draw distinct colors from the palette
// (and, if that is ever exhausted, evenly spaced HSL hues) so every element is
// visually separable. Elements are sorted first so the assignment is
// deterministic regardless of atom order in the file.
export const buildElementColors = (elements) => {
    const map = {};
    const used = new Set();
    let cursor = 0;
    [...new Set(elements)].sort().forEach((element) => {
        let color = ELEMENT_COLORS[element];
        if (!color || used.has(color)) {
            while (cursor < FALLBACK_PALETTE.length && used.has(FALLBACK_PALETTE[cursor])) {
                cursor += 1;
            }
            color = FALLBACK_PALETTE[cursor] || `hsl(${(Object.keys(map).length * 47) % 360}, 70%, 60%)`;
            cursor += 1;
        }
        map[element] = color;
        used.add(color);
    });
    return map;
};
