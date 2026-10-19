// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The app's own phrasings (private to verification/): the words that make a
// kind of claim, as global regexes whose group 1 is the character before the
// claim (reading.js positiveMatches). The claim catalogue (claims.js) and the
// facts catalogue (facts.js) name the ones they read; a phrasing both read
// (a preferred direction, a superstructure, distorted clusters, a
// coordination number) is one pattern here.
//
// The words of the validation's verdict schema (off_centre, cluster_breathing,
// jahn_teller, symmetry_lowering, "features": [] …) are not the app's: the
// harness adds them to its turns (validation/verdict.js VERDICT_PHRASINGS,
// the one seam of createVerification's `phrasings`).

// --- the answer checks' claims (facts.js) -----------------------------------------

export const CONVERGED = /(^|[^\w-])(converged|reached convergence|at convergence|plateau(?:ed|s)?)\b/gi;
export const STILL_MOVING = /(^|\W)(improving|diverg\w*|still (?:falling|decreasing|dropping|declining|descending|going down|moving))\b/gi;

export const CHIRAL = /(^|[^\w-])chiral\b/gi;
export const CENTROSYMMETRIC = /(^|[^\w-])centrosymmetric\b/gi;
// "polar" as a property ("non-polar", "nonpolar" and "polarisation" do not match).
export const POLAR = /(^|[^\w-])polar\b(?!\s+(?:plot|coordinates?|angle))/gi;
export const LOST_INVERSION = [
    /\b(los(?:e|es|t|ing)|break(?:s|ing)?|broke|broken|remov(?:e|es|ed|ing)|destroy(?:s|ed|ing)?)\b[^.;]{0,50}\b(inversion|centrosymmetr\w*|cent(?:re|er) of symmetry)/gi,
    /\b(inversion|centrosymmetr\w*|cent(?:re|er) of symmetry)\b[^.;]{0,30}\b(is|was|are|were|gets|got)\s+(lost|broken|removed|destroyed)/gi
];

// --- read by both catalogues ------------------------------------------------------

export const PREFERRED = /(^|[^\w-])(preferred|preferential)\s+(?:displacement\s+)?(direction|orientation|axis)\b/gi;

// A claimed superstructure: cell doubling or ordered zone-boundary
// displacements (needs a zone-boundary point), or order between cells in
// general (a zone-boundary point or a general wavevector will do).
export const CELL_DOUBLING = /(^|[^\w-])(cell[- ]doubling|doubl(?:es|ed|ing) (?:of )?the (?:unit )?cell|zone[- ]boundary (?:order|modulation|mode|distortion)|commensurate modulation)\b/gi;
export const SUPERSTRUCTURE = /(^|[^\w-])(superstructure|supercell order|ordered (?:\w+ )?displacements? at q|incommensurate (?:order|modulation))\b/gi;

const CLUSTER_NOUN = '(?:clusters?|tetrahedr(?:on|a)|trimers?|dimers?|triangles?|octahedr(?:on|a)|polyhedr(?:on|a))';
// An adverb may stand between the verb and "distorted", but not a negation.
export const DISTORTED = new RegExp(`(^|[^\\w-])(distorted\\s+(?:\\w+\\s+)?${CLUSTER_NOUN}|${CLUSTER_NOUN}\\s+(?:are|is|appear|look|seem)\\s+`
    + '(?:(?!not\\b|never\\b|hardly\\b|barely\\b)\\w+\\s+)?(?:distorted|deformed)|cluster\\s+distortions?)\\b', 'gi');
export const BREATHING = /(^|[^\w-])(breath(?:e|es|ing))\b/gi;
export const RIGID = /(^|[^\w-])(rigid(?:ly)?\b|move(?:s|d)?\s+(?:together\s+)?as\s+(?:a\s+)?(?:rigid\s+)?(?:units?|bod(?:y|ies))\b)/gi;

export const WORD_NUMBER = { two: 2, three: 3, four: 4, five: 5, six: 6, seven: 7, eight: 8, nine: 9, ten: 10, eleven: 11, twelve: 12 };
const COUNT = '(\\d+(?:\\.\\d+)?|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve)';
// A coordination number, its count in group 2.
export const CN_CLAIMS = [
    new RegExp(`(^|[^\\w.-])${COUNT}[- ]?(?:fold[- ]?)?coordinat(?:e|ed|ion)\\b`, 'gi'),
    new RegExp(`(^|[^\\w-])(?:coordination(?:\\s+number)?|CN)\\s*(?:of|is|=|equals?|:)?\\s*${COUNT}(?!\\d|\\.\\d)`, 'gi')
];
// A coordination named by its polyhedron: [pattern, count].
export const POLYHEDRON_CN = [
    [/(^|[^\w-])(tetrahedral(?:ly)?)\s+(?:\w+\s+)?coordinat/gi, 4],
    [/(^|[^\w-])(octahedral(?:ly)?)\s+(?:\w+\s+)?coordinat/gi, 6],
    [/(^|[^\w-])(square[- ]planar)\s+(?:\w+\s+)?coordinat/gi, 4],
    [/(^|[^\w-])(cub(?:ic|ically))\s+(?:\w+\s+)?coordinat/gi, 8]
];
// A clause about part of the shell (the first peak, the short bonds, a window
// in Å) is not a claim about the whole shell's count.
export const PARTIAL_SHELL = /\bfirst[- ]peak\b|\bshort(?:est|er)?\s+bonds?\b|\bwithin\s+[\d.]+\s*(?:[–—−-]|to)\s*[\d.]+\s*(?:Å|A\b)|\d+(?:\.\d+)?\s*[–—−-]\s*\d+(?:\.\d+)?\s*(?:Å|A\b)/i;

// --- the evidence gate's claims (claims.js) ---------------------------------------

export const MODULATION = /(^|[^\w-])(modulat(?:ion|ions|ed)|incommensurate|order(?:ed)? between cells|correlated displacements between cells)\b/gi;
// A split site: named, random off-centring either way round, a ± displacement
// (not a ± sampling error; not random thermal displacements, which are no order).
export const SPLIT_SITE = /(^|[^\w-])(split[- ]sites?|two[- ]site(?:\s+disorder)?|disordered over two|random(?:ly)?\s+(?:(?!thermal\b)\w+\s+)?(?:off[- ]cent\w*|displac\w*)|(?:off[- ]cent\w*|displac\w*)[^.;]{0,40}?\bat random|(?:displac\w*|off[- ]cent\w*|shift\w*|split)\s+(?:by\s+)?±\s*\d+(?:\.\d+)?\s*(?:Å|A\b)|±\s*\d+(?:\.\d+)?\s*(?:Å|A\b)\s+along|off[- ]cent(?:re|er|ring|ering)\s+(?:dis)?order)/gi;
export const SHELL_CLAIMS = /(^|[^\w-])(jahn[-– ]teller|split (?:first )?(?:coordination )?shell|short and long bonds)\b/gi;
// A lowering of the average structure (lowers / lowered its symmetry), or a
// quantified shift every copy shares ("shifted by 0.07 Å", "(coherent)"); a ±
// shift is a split site, "coherent scattering" none, and "lower symmetry
// locally" or a local symmetry lowering is about the clusters, not the average.
export const LOWERING = /(^|[^\w-])((?<!local\s)symmetry[- ]lowering|lower(?:s|ed|ing) (?:the )?(?:average )?symmetry|(?:average structure|space group)\s+(?:is\s+)?lowered|coherent(?:ly)? (?:shift\w*|off[- ]cent\w*|displace\w*)|polar shift|shift(?:ed|s)? (?:by )?\d+(?:\.\d+)?\s*(?:Å|A\b)|\(coherent\))/gi;
export const NO_ORDER = /(^|[^\w-])(no (?:local |short-range |additional )?order(?!\s+(?:that|which)\s+doubles)|no (?:ordered|correlated) (?:displacements?|distortions?)|(?:only|purely) thermal motion|thermal motion only|nothing beyond (?:the )?(?:average structure|thermal motion)|no (?:sign|evidence) of (?:local )?order)\b/gi;
// A sentence about this run names the run itself (this structure, the
// configuration, the clusters, site 3), or one of its elements (claims.js).
export const RUN_NOUN = /\b(?:th(?:is|e) (?:run|configuration|structure|model|box|supercell|sample|compound|crystal|material|system)|the (?:sites?|atoms|clusters|octahedra|tetrahedra|polyhedra|displacements)|sites? \d+)\b/i;
