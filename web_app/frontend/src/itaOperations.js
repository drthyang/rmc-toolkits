// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/itaOperations.js
//
// The symmetry operations of the 230 space groups in their International Tables (ITA)
// standard settings, as the CIF export needs them to move a structure onto ITA's own origin
// (averageStructure.js → itaOriginShift). Each group is its lattice letter and a few
// generator triplets; itaOperations() closes them, with the centring, into the full set.
//
// The settings are those of the Wyckoff table (wyckoffTable.js): unique axis b and cell
// choice 1 for monoclinic groups, ORIGIN CHOICE 2 (an inversion centre at the origin) for the
// 24 groups ITA gives two origins, and HEXAGONAL axes (obverse) for the 7 R groups. Every set
// is pinned to spglib's operations for the same setting by itaOperations.test.js
// (tests/generate_ita_operations_fixture.py), and every Wyckoff row to its multiplicity under
// them by wyckoff.test.js.

import { parseCoordinateForm } from './wyckoff.js';

/** The 24 groups ITA describes on two origins; these are origin choice 2 here. */
export const ORIGIN_CHOICE_2 = new Set([48, 50, 59, 68, 70, 85, 86, 88, 125, 126, 129, 130, 133, 134, 137, 138, 141, 142, 201, 203, 222, 224, 227, 228]);

/** Centring translations of each lattice letter (R: obverse, hexagonal axes). */
export const CENTRING = {
  P: [[0, 0, 0]],
  A: [[0, 0, 0], [0, 1 / 2, 1 / 2]],
  B: [[0, 0, 0], [1 / 2, 0, 1 / 2]],
  C: [[0, 0, 0], [1 / 2, 1 / 2, 0]],
  I: [[0, 0, 0], [1 / 2, 1 / 2, 1 / 2]],
  F: [[0, 0, 0], [0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]],
  R: [[0, 0, 0], [2 / 3, 1 / 3, 1 / 3], [1 / 3, 2 / 3, 2 / 3]],
};

// number: 'lattice letter|generator;generator;…'
const ITA_GENERATORS = {
  1: 'P|',
  2: 'P|-x,-y,-z',
  3: 'P|-x,y,-z',
  4: 'P|-x,y+1/2,-z',
  5: 'C|-x,y,-z',
  6: 'P|x,-y,z',
  7: 'P|x,-y,z+1/2',
  8: 'C|x,-y,z',
  9: 'C|x,-y,z+1/2',
  10: 'P|-x,y,-z;-x,-y,-z',
  11: 'P|-x,y+1/2,-z;-x,-y,-z',
  12: 'C|-x,y,-z;-x,-y,-z',
  13: 'P|-x,y,-z+1/2;-x,-y,-z',
  14: 'P|-x,y+1/2,-z+1/2;-x,-y,-z',
  15: 'C|-x,y,-z+1/2;-x,-y,-z',
  16: 'P|-x,-y,z;x,-y,-z',
  17: 'P|-x,-y,z+1/2;-x,y,-z+1/2',
  18: 'P|-x,-y,z;-x+1/2,y+1/2,-z',
  19: 'P|-x+1/2,-y,z+1/2;-x,y+1/2,-z+1/2',
  20: 'C|-x,-y,z+1/2;x,-y,-z',
  21: 'C|-x,-y,z;x,-y,-z',
  22: 'F|-x,-y,z;x,-y,-z',
  23: 'I|-x,-y,z;x,-y,-z',
  24: 'I|-x+1/2,-y,z+1/2;-x,y+1/2,-z+1/2',
  25: 'P|-x,-y,z;x,-y,z',
  26: 'P|-x,-y,z+1/2;x,-y,z+1/2',
  27: 'P|-x,-y,z;x,-y,z+1/2',
  28: 'P|-x,-y,z;x+1/2,-y,z',
  29: 'P|-x,-y,z+1/2;x+1/2,-y,z',
  30: 'P|-x,-y,z;x,-y+1/2,z+1/2',
  31: 'P|-x+1/2,-y,z+1/2;-x,y,z',
  32: 'P|-x,-y,z;x+1/2,-y+1/2,z',
  33: 'P|-x,-y,z+1/2;x+1/2,-y+1/2,z',
  34: 'P|-x,-y,z;x+1/2,-y+1/2,z+1/2',
  35: 'C|-x,-y,z;x,-y,z',
  36: 'C|-x,-y,z+1/2;x,-y,z+1/2',
  37: 'C|-x,-y,z;x,-y,z+1/2',
  38: 'A|-x,-y,z;x,-y,z',
  39: 'A|-x,-y,z;x,-y+1/2,z',
  40: 'A|-x,-y,z;x+1/2,-y,z',
  41: 'A|-x,-y,z;-x+1/2,y+1/2,z',
  42: 'F|-x,-y,z;-x,y,z',
  43: 'F|-x,-y,z;-x+1/4,y+1/4,z+1/4',
  44: 'I|-x,-y,z;-x,y,z',
  45: 'I|-x,-y,z;-x,y,z+1/2',
  46: 'I|-x,-y,z;-x+1/2,y,z',
  47: 'P|-x,-y,z;x,-y,-z;-x,-y,-z',
  48: 'P|-x+1/2,-y+1/2,z;x,-y+1/2,-z+1/2;-x,-y,-z',
  49: 'P|-x,-y,z;x,-y,-z+1/2;-x,-y,-z',
  50: 'P|-x+1/2,-y+1/2,z;x,-y+1/2,-z;-x,-y,-z',
  51: 'P|-x+1/2,-y,z;x+1/2,-y,-z;-x,-y,-z',
  52: 'P|-x+1/2,-y,z;x,-y+1/2,-z+1/2;-x,-y,-z',
  53: 'P|-x+1/2,-y,z+1/2;x,-y,-z;-x,-y,-z',
  54: 'P|-x+1/2,-y,z;x+1/2,-y,-z+1/2;-x,-y,-z',
  55: 'P|-x,-y,z;x+1/2,-y+1/2,-z;-x,-y,-z',
  56: 'P|-x+1/2,-y+1/2,z;x+1/2,-y,-z+1/2;-x,-y,-z',
  57: 'P|-x,-y,z+1/2;x,-y+1/2,-z;-x,-y,-z',
  58: 'P|-x,-y,z;x+1/2,-y+1/2,-z+1/2;-x,-y,-z',
  59: 'P|-x+1/2,-y+1/2,z;x+1/2,-y,-z;-x,-y,-z',
  60: 'P|-x+1/2,-y+1/2,z+1/2;x+1/2,-y+1/2,-z;-x,-y,-z',
  61: 'P|-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,-y,-z',
  62: 'P|-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z+1/2;-x,-y,-z',
  63: 'C|-x,-y,z+1/2;x,-y,-z;-x,-y,-z',
  64: 'C|-x+1/2,-y,z+1/2;x,-y,-z;-x,-y,-z',
  65: 'C|-x,-y,z;x,-y,-z;-x,-y,-z',
  66: 'C|-x,-y,z;x,-y,-z+1/2;-x,-y,-z',
  67: 'C|-x+1/2,-y,z;x,-y,-z;-x,-y,-z',
  68: 'C|-x+1/2,-y,z;x+1/2,-y,-z+1/2;-x,-y,-z',
  69: 'F|-x,-y,z;x,-y,-z;-x,-y,-z',
  70: 'F|-x+1/4,-y+1/4,z;x,-y+1/4,-z+1/4;-x,-y,-z',
  71: 'I|-x,-y,z;x,-y,-z;-x,-y,-z',
  72: 'I|-x,-y,z;x,-y,-z+1/2;-x,-y,-z',
  73: 'I|-x,-y+1/2,z;x,-y,-z+1/2;-x,-y,-z',
  74: 'I|-x,-y+1/2,z;x,-y,-z;-x,-y,-z',
  75: 'P|-y,x,z',
  76: 'P|-y,x,z+1/4',
  77: 'P|-y,x,z+1/2',
  78: 'P|-y,x,z+3/4',
  79: 'I|-y,x,z',
  80: 'I|-y,x+1/2,z+1/4',
  81: 'P|y,-x,-z',
  82: 'I|y,-x,-z',
  83: 'P|-y,x,z;-x,-y,-z',
  84: 'P|-y,x,z+1/2;-x,-y,-z',
  85: 'P|-y+1/2,x,z;-x,-y,-z',
  86: 'P|-y,x+1/2,z+1/2;-x,-y,-z',
  87: 'I|-y,x,z;-x,-y,-z',
  88: 'I|-y+3/4,x+1/4,z+1/4;-x,-y,-z',
  89: 'P|-y,x,z;x,-y,-z',
  90: 'P|-y+1/2,x+1/2,z;y,x,-z',
  91: 'P|-y,x,z+1/4;-x,y,-z',
  92: 'P|-y+1/2,x+1/2,z+1/4;y,x,-z',
  93: 'P|-y,x,z+1/2;x,-y,-z',
  94: 'P|-y+1/2,x+1/2,z+1/2;y,x,-z',
  95: 'P|-y,x,z+3/4;-x,y,-z',
  96: 'P|-y+1/2,x+1/2,z+3/4;y,x,-z',
  97: 'I|-y,x,z;x,-y,-z',
  98: 'I|-y,x+1/2,z+1/4;-y,-x,-z',
  99: 'P|-y,x,z;x,-y,z',
  100: 'P|-y,x,z;x+1/2,-y+1/2,z',
  101: 'P|-y,x,z+1/2;x,-y,z+1/2',
  102: 'P|-y+1/2,x+1/2,z+1/2;y,x,z',
  103: 'P|-y,x,z;x,-y,z+1/2',
  104: 'P|-y,x,z;x+1/2,-y+1/2,z+1/2',
  105: 'P|-y,x,z+1/2;x,-y,z',
  106: 'P|-y,x,z+1/2;x+1/2,-y+1/2,z',
  107: 'I|-y,x,z;x,-y,z',
  108: 'I|-y,x,z;x,-y,z+1/2',
  109: 'I|-y,x+1/2,z+1/4;-x,y,z',
  110: 'I|-y,x+1/2,z+1/4;-x,y,z+1/2',
  111: 'P|y,-x,-z;x,-y,-z',
  112: 'P|y,-x,-z;x,-y,-z+1/2',
  113: 'P|y,-x,-z;x+1/2,-y+1/2,-z',
  114: 'P|y,-x,-z;x+1/2,-y+1/2,-z+1/2',
  115: 'P|y,-x,-z;y,x,-z',
  116: 'P|y,-x,-z;y,x,-z+1/2',
  117: 'P|y,-x,-z;y+1/2,x+1/2,-z',
  118: 'P|y,-x,-z;y+1/2,x+1/2,-z+1/2',
  119: 'I|y,-x,-z;y,x,-z',
  120: 'I|y,-x,-z;y,x,-z+1/2',
  121: 'I|y,-x,-z;x,-y,-z',
  122: 'I|y,-x,-z;x,-y+1/2,-z+1/4',
  123: 'P|-y,x,z;x,-y,-z;-x,-y,-z',
  124: 'P|-y,x,z;x,-y,-z+1/2;-x,-y,-z',
  125: 'P|-y+1/2,x,z;x,-y+1/2,-z;-x,-y,-z',
  126: 'P|-y+1/2,x,z;x,-y+1/2,-z+1/2;-x,-y,-z',
  127: 'P|-y,x,z;x+1/2,-y+1/2,-z;-x,-y,-z',
  128: 'P|-y,x,z;x+1/2,-y+1/2,-z+1/2;-x,-y,-z',
  129: 'P|-y+1/2,x,z;x+1/2,-y,-z;-x,-y,-z',
  130: 'P|-y+1/2,x,z;x+1/2,-y,-z+1/2;-x,-y,-z',
  131: 'P|-y,x,z+1/2;x,-y,-z;-x,-y,-z',
  132: 'P|-y,x,z+1/2;x,-y,-z+1/2;-x,-y,-z',
  133: 'P|-y+1/2,x,z+1/2;x,-y+1/2,-z;-x,-y,-z',
  134: 'P|-y+1/2,x,z+1/2;x,-y+1/2,-z+1/2;-x,-y,-z',
  135: 'P|-y,x,z+1/2;x+1/2,-y+1/2,-z;-x,-y,-z',
  136: 'P|-y+1/2,x+1/2,z+1/2;x+1/2,-y+1/2,-z+1/2;-x,-y,-z',
  137: 'P|-y+1/2,x,z+1/2;x+1/2,-y,-z;-x,-y,-z',
  138: 'P|-y+1/2,x,z+1/2;x+1/2,-y,-z+1/2;-x,-y,-z',
  139: 'I|-y,x,z;x,-y,-z;-x,-y,-z',
  140: 'I|-y,x,z;x,-y,-z+1/2;-x,-y,-z',
  141: 'I|-y+1/4,x+3/4,z+1/4;x,-y,-z;-x,-y,-z',
  142: 'I|-y+1/4,x+3/4,z+1/4;x,-y,-z+1/2;-x,-y,-z',
  143: 'P|-y,x-y,z',
  144: 'P|-y,x-y,z+1/3',
  145: 'P|-y,x-y,z+2/3',
  146: 'R|-y,x-y,z',
  147: 'P|-y,x-y,z;-x,-y,-z',
  148: 'R|-y,x-y,z;-x,-y,-z',
  149: 'P|-y,x-y,z;-y,-x,-z',
  150: 'P|-y,x-y,z;x-y,-y,-z',
  151: 'P|-y,x-y,z+1/3;x,x-y,-z',
  152: 'P|-y,x-y,z+1/3;y,x,-z',
  153: 'P|-y,x-y,z+2/3;x,x-y,-z',
  154: 'P|-y,x-y,z+2/3;y,x,-z',
  155: 'R|-y,x-y,z;y,x,-z',
  156: 'P|-y,x-y,z;-y,-x,z',
  157: 'P|-y,x-y,z;y,x,z',
  158: 'P|-y,x-y,z;-y,-x,z+1/2',
  159: 'P|-y,x-y,z;y,x,z+1/2',
  160: 'R|-y,x-y,z;-y,-x,z',
  161: 'R|-y,x-y,z;-y,-x,z+1/2',
  162: 'P|-y,x-y,z;-y,-x,-z;-x,-y,-z',
  163: 'P|-y,x-y,z;-y,-x,-z+1/2;-x,-y,-z',
  164: 'P|-y,x-y,z;y,x,-z;-x,-y,-z',
  165: 'P|-y,x-y,z;y,x,-z+1/2;-x,-y,-z',
  166: 'R|-y,x-y,z;y,x,-z;-x,-y,-z',
  167: 'R|-y,x-y,z;y,x,-z+1/2;-x,-y,-z',
  168: 'P|x-y,x,z',
  169: 'P|x-y,x,z+1/6',
  170: 'P|x-y,x,z+5/6',
  171: 'P|x-y,x,z+1/3',
  172: 'P|x-y,x,z+2/3',
  173: 'P|x-y,x,z+1/2',
  174: 'P|-x+y,-x,-z',
  175: 'P|x-y,x,z;-x,-y,-z',
  176: 'P|x-y,x,z+1/2;-x,-y,-z',
  177: 'P|x-y,x,z;x-y,-y,-z',
  178: 'P|x-y,x,z+1/6;x-y,-y,-z',
  179: 'P|x-y,x,z+5/6;x-y,-y,-z',
  180: 'P|x-y,x,z+1/3;x-y,-y,-z',
  181: 'P|x-y,x,z+2/3;x-y,-y,-z',
  182: 'P|x-y,x,z+1/2;x-y,-y,-z',
  183: 'P|x-y,x,z;x-y,-y,z',
  184: 'P|x-y,x,z;x-y,-y,z+1/2',
  185: 'P|x-y,x,z+1/2;x-y,-y,z',
  186: 'P|x-y,x,z+1/2;x-y,-y,z+1/2',
  187: 'P|-x+y,-x,-z;-x+y,y,z',
  188: 'P|-x+y,-x,-z+1/2;-x+y,y,z+1/2',
  189: 'P|-x+y,-x,-z;x-y,-y,z',
  190: 'P|-x+y,-x,-z+1/2;x-y,-y,z+1/2',
  191: 'P|x-y,x,z;x-y,-y,-z;-x,-y,-z',
  192: 'P|x-y,x,z;x-y,-y,-z+1/2;-x,-y,-z',
  193: 'P|x-y,x,z+1/2;x-y,-y,-z+1/2;-x,-y,-z',
  194: 'P|x-y,x,z+1/2;x-y,-y,-z;-x,-y,-z',
  195: 'P|-x,-y,z;z,x,y',
  196: 'F|-x,-y,z;z,x,y',
  197: 'I|-x,-y,z;z,x,y',
  198: 'P|-x+1/2,-y,z+1/2;z,x,y',
  199: 'I|-x+1/2,-y,z+1/2;z,x,y',
  200: 'P|-x,-y,z;z,x,y;-x,-y,-z',
  201: 'P|-x+1/2,-y+1/2,z;z,x,y;-x,-y,-z',
  202: 'F|-x,-y,z;z,x,y;-x,-y,-z',
  203: 'F|-x+1/4,-y+1/4,z;z,x,y;-x,-y,-z',
  204: 'I|-x,-y,z;z,x,y;-x,-y,-z',
  205: 'P|-x+1/2,-y,z+1/2;z,x,y;-x,-y,-z',
  206: 'I|-x+1/2,-y,z+1/2;z,x,y;-x,-y,-z',
  207: 'P|-y,x,z;z,x,y',
  208: 'P|-y+1/2,x+1/2,z+1/2;z,x,y',
  209: 'F|-y,x,z;z,x,y',
  210: 'F|-y+1/4,x+1/4,z+1/4;z,x,y',
  211: 'I|-y,x,z;z,x,y',
  212: 'P|-y+3/4,x+1/4,z+3/4;z,x,y',
  213: 'P|-y+1/4,x+3/4,z+1/4;z,x,y',
  214: 'I|-y+1/4,x+3/4,z+1/4;z,x,y',
  215: 'P|y,-x,-z;z,x,y',
  216: 'F|y,-x,-z;z,x,y',
  217: 'I|y,-x,-z;z,x,y',
  218: 'P|y+1/2,-x+1/2,-z+1/2;z,x,y',
  219: 'F|y+1/2,-x,-z;z,x,y',
  220: 'I|y+1/4,-x+3/4,-z+1/4;z,x,y',
  221: 'P|-y,x,z;z,x,y;-x,-y,-z',
  222: 'P|-y+1/2,x,z;z,x,y;-x,-y,-z',
  223: 'P|-y+1/2,x+1/2,z+1/2;z,x,y;-x,-y,-z',
  224: 'P|-y,x+1/2,z+1/2;z,x,y;-x,-y,-z',
  225: 'F|-y,x,z;z,x,y;-x,-y,-z',
  226: 'F|-y+1/2,x,z;z,x,y;-x,-y,-z',
  227: 'F|-y,x+1/4,z+1/4;z,x,y;-x,-y,-z',
  228: 'F|-y+1/2,x+1/4,z+1/4;z,x,y;-x,-y,-z',
  229: 'I|-y,x,z;z,x,y;-x,-y,-z',
  230: 'I|-y+1/4,x+3/4,z+1/4;z,x,y;-x,-y,-z',
};

const wrap = (v) => {
  const y = v - Math.floor(v);
  return y < 1e-9 || y > 1 - 1e-9 ? 0 : y;
};
const key = ({ R, t }) => `${R.flat().join(',')}|${t.map((v) => Math.round(wrap(v) * 48) % 48).join(',')}`;
const compose = (a, b) => ({
  R: a.R.map((row) => [0, 1, 2].map((j) => row[0] * b.R[0][j] + row[1] * b.R[1][j] + row[2] * b.R[2][j])),
  t: [0, 1, 2].map((i) => wrap(a.R[i][0] * b.t[0] + a.R[i][1] * b.t[1] + a.R[i][2] * b.t[2] + a.t[i])),
});

/** The lattice letter and generator triplets of a group (null for a number outside 1–230). */
export function itaGenerators(number) {
  const packed = ITA_GENERATORS[number];
  if (!packed) return null;
  const [letter, list] = packed.split('|');
  return { letter, centring: CENTRING[letter], generators: list ? list.split(';') : [] };
}

const cache = new Map();

/**
 * Every operation {R, t} of the group in its ITA standard setting, centring included, t in
 * [0, 1): the generators closed under composition. Cached; null for an unknown number.
 */
export function itaOperations(number) {
  if (cache.has(number)) return cache.get(number);
  const spec = itaGenerators(number);
  if (!spec) return null;
  const ops = new Map();
  const add = (o) => {
    const k = key(o);
    if (ops.has(k)) return false;
    ops.set(k, { R: o.R, t: o.t.map(wrap) });
    return true;
  };
  add({ R: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: [0, 0, 0] });
  for (const g of spec.generators) add(parseCoordinateForm(g));
  for (const c of spec.centring) add({ R: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], t: c });
  for (let grew = true; grew;) {
    grew = false;
    const snapshot = [...ops.values()];
    for (const a of snapshot) for (const b of snapshot) if (add(compose(a, b))) grew = true;
  }
  const out = [...ops.values()];
  cache.set(number, out);
  return out;
}
