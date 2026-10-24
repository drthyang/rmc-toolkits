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

// A group–subgroup relation, for the subgroup check (facts.js; the
// local-symmetry interview, 2026-10-24, answer 4a), read in symbolSpelling's
// lower-case text: "H is a (maximal) (t-/k-) subgroup of G", "G is a
// supergroup of H", "G has H as a subgroup", and a lowering from G to H
// ("lowers G to H", "from G to H", "G → H") in a clause that says lower,
// reduce, descend or break. The claimed index and type beside it.
export const SUBGROUP_OF = /(^|\W)(subgroups?\s+of)\b/g;
export const SUPERGROUP_OF = /(^|\W)(supergroups?\s+of)\b/g;
export const HAS_AS_SUBGROUP = /(^|[^\w-])(has)\b[^.;]{0,40}?\bas\s+(?:an?\s+|its\s+)?(?:\w+[-\s]+){0,2}subgroups?\b/g;
export const LOWERING_TO = /(^|[^\w-])(lower(?:s|ed|ing)?|reduc(?:e|es|ed|ing|tion)|descen(?:t|ds|ding)|break(?:s|ing)?\s+down)\b/g;
// What stands between the two symbols of a lowering step: only the
// connector, with a lowering verb at most ("F-43m to P-42₁m", "F-43m →
// P-42₁m", "F-43m lowers to P-42₁m"; a bracket between them, a tolerance,
// is read as a space).
export const ARROW = /^\s*(?:(?:is\s+)?(?:lower(?:s|ed)?|reduc(?:es|ed)|descends?|goes|drops|breaks\s+down)\s+)?(?:down\s+)?(?:→|->|⟶|⇒|to|into)\s*$/;
export const INDEX_CLAIM = /\bindex(?:\s+of)?\s*(?:=|:|is)?\s*(\d+)\b|\b[tk](\d+)\b/;
export const T_TYPE = /\bt-?\s?subgroups?\b|\btranslationengleiche\b|\btranslation-equivalent\b|\bt\d+\b/;
export const K_TYPE = /\bk-?\s?subgroups?\b|\bklassengleiche\b|\bclass-equivalent\b|\bk\d+\b|\bisomorphic\s+subgroups?\b/;
export const MAXIMAL = /\bmaximal\b/;

// A cell doubling, for the box check (facts.js; the local-symmetry
// interview, 2026-10-24, answer 4b): CELL_DOUBLING, or a doubled cell or axis
// in other words ("a doubled cell", "doubles c", "the c axis is doubled", "2a
// × 2a × c"); the axes it names (along a, b, c; a ½ wavevector component; a
// factor 2 in a cell); and an answer that says the box cannot hold it.
// An axis named by its letter: "the c axis", "c-axis", b or c alone, and a
// alone only where no word follows it (not the article in "along a direction").
const AXIS_WORD = '(?:[abc][- ]?axis|[bc](?![a-z])|a(?=\\s*(?:[,.;:)]|and\\b|$)))';
export const DOUBLED_CELL = new RegExp('(^|[^\\w-])(doubled\\s+(?:unit\\s+)?cell|cell\\s+(?:is\\s+)?doubled|doubl(?:e|es|ed|ing)\\s+(?:the\\s+)?'
    + `(?:unit\\s+cell|cell|${AXIS_WORD}|lattice parameter)|doubl(?:es|ing)\\s+along|[abc][- ]?axis\\s+(?:is\\s+)?doubled`
    + '|2\\s*[abc]\\s*[×x*]\\s*\\d?\\s*[abc]\\s*[×x*]\\s*\\d?\\s*[abc])', 'gi');
export const ALONG_AXES = new RegExp(`\\balong\\s+(?:the\\s+)?(${AXIS_WORD})(?:\\s*(?:,|and)\\s*(?:the\\s+)?(${AXIS_WORD}))?`
    + `(?:\\s*(?:,|and)\\s*(?:the\\s+)?(${AXIS_WORD}))?`, 'i');
export const DOUBLED_AXIS = new RegExp(`\\bdoubl(?:es|ed|ing)\\s+(?:of\\s+)?(?:the\\s+)?(${AXIS_WORD})|\\b([abc])[- ]?axis\\s+(?:is\\s+)?doubled`
    + '|\\b2\\s*([abc])(?![a-z])', 'gi');
export const ALL_AXES = /\balong\s+(?:all\s+)?(?:three|every|each)\s+(?:axes|axis)\b|\b2a\s*[×x*]\s*2b\s*[×x*]\s*2c\b/i;
export const HALF_TUPLE = /\(\s*(½|1\/2|0?\.5|0)\s*[, ]\s*(½|1\/2|0?\.5|0)\s*[, ]\s*(½|1\/2|0?\.5|0)\s*\)/g;
export const CANNOT_HOLD = /\b(?:can(?:not|'t)|could(?:n't|\s+not)|does(?:n't|\s+not)|do(?:n't|\s+not)|too\s+small\s+to)\b[^.;]{0,40}\b(?:hold|fit|represent|accommodate|contain|show|support)\b|\bneither\s+confirm\s+nor\s+rule\s+out\b|\bcan\s+neither\b/i;

// --- read by both catalogues ------------------------------------------------------

export const PREFERRED = /(^|[^\w-])(preferred|preferential)\s+(?:displacement\s+)?(direction|orientation|axis)\b/gi;

// A claimed superstructure: cell doubling or ordered zone-boundary
// displacements (needs a zone-boundary point), or order between cells in
// general (a zone-boundary point or a general wavevector will do).
export const CELL_DOUBLING = /(^|[^\w-])(cell[- ]doubling|doubl(?:es|ed|ing) (?:of )?the (?:unit )?cell|zone[- ]boundary (?:order|modulation|mode|distortion)|commensurate modulation)\b/gi;
export const SUPERSTRUCTURE = /(^|[^\w-])(superstructure|supercell order|ordered (?:\w+ )?displacements? at q|incommensurate (?:order|modulation))\b/gi;

// Order at a wavevector, for the long-wave check (facts.js; the user's rule,
// 2026-10-24): any word that calls displacements ordered between cells, or a
// modulation ("in order to" is no order, "disorder" no word of it). A long
// wave named in words: a long-wavelength or box-scale wave or mode.
export const ORDER_CLAIM = /(^|[^\w-])((?<!\bin\s)order(?:ed|ing)?|superstructure|modulat(?:ion|ions|ed)|incommensurate|cell[- ]doubling|correlated (?:\w+ )?displacements?)\b/gi;
export const LONG_WAVE_WORDS = /\b(?:long[- ]wave(?:length)?s?|longest[- ]wave(?:length)?s?|box[- ](?:scale|size[d]?))\b/i;
// "Nothing" negates a claim of order for the long-wave check ("it says
// nothing about order at its own q", the user's own wording).
export const SAYS_NOTHING = /\bnothing\b/i;

// A distribution's numbers, for the single-configuration check (facts.js;
// the user's rule, 2026-10-24): a statistic of a distribution followed by a
// number (mean 109.4°, sd 4.1°, mode 109.5°, peaks at 100–110°), a share of
// it (80 % cis, 20 % of the angles lie …), or a coordination histogram
// (3 948 of 4 000 Ga atoms having 4 neighbours) — in a sentence that names a
// distribution (angles, a histogram, neighbours, cluster shapes, edges or the
// edge ladder, breathing, a direction map) or gives degrees. A test's σ is no number of the
// distribution.
export const DISTRIBUTION_STAT = new RegExp([
    '\\b(?:mean|average|median|sd|s\\.d\\.|standard deviation|spread|width|mode|modal(?: share)?|share)\\b\\s*(?:of|is|=|:|at)?\\s*[~≈]?\\s*[-−]?\\d',
    '\\b(?:peaks?|centred|centered)\\s+at\\s+[-−]?\\d',
    '\\d+(?:\\.\\d+)?\\s?%\\s+(?:of\\s+(?:the|these|all)\\s+|cis\\b|trans\\b|lie\\b|lies\\b|fall\\b|falls\\b|are\\b|have\\b|has\\b)',
    '\\d[\\d\\s,]*\\s+of\\s+[\\d\\s,]+\\s*(?:[A-Z][a-z]?\\s+)?atoms?\\s+(?:have|having|with)\\s+\\d'
].join('|'), 'i');
export const DISTRIBUTION_NOUN = /\b(?:bond[- ])?angles?\b|\bangular\b|\bdistributions?\b|\bhistograms?\b|\bneighbou?rs\b|\bcoordination counts?\b|\b(?:cluster|shape)s?\s+(?:shapes?|rms)\b|\bedge(?:[- ]lengths?|s|\s+ladder)\b|\bbreathing\b|\bdirection(?:al)?\s+(?:map|distribution|histogram)s?\b|\d\s?°/i;
// The answer says the numbers come from one configuration.
export const SINGLE_CONFIGURATION = /\b(?:one|a single|single|only one|just one|this one|the one)\s+(?:RMC\s+|RMCProfile\s+)?(?:configuration|snapshot|box|run|realisation|realization)s?\b|\bsingle[- ](?:configuration|snapshot)\b|\b(?:configuration|snapshot)\s+alone\b/i;

// A distorted cluster or polyhedron: "distorted clusters", "the octahedra are
// (strongly) distorted", "a cluster distortion". An adverb may stand between
// the verb and "distorted", but not a negation. The facts catalogue reads
// every noun (DISTORTED, against analyze_clusters' facts); the claim
// catalogue reads the clusters' (DISTORTED_CLUSTERS: analyze_clusters) apart
// from the coordination polyhedra's (DISTORTED_POLYHEDRA: polyhedral_distortion
// or analyze_clusters; the coordination-shells interview, 2026-10-24, answer
// 7: a corner-sharing polyhedron or another extended network is read with a
// distortion index, never the cluster null).
const CLUSTER_ONLY_NOUN = '(?:clusters?|trimers?|dimers?|triangles?)';
const POLYHEDRON_NOUN = '(?:tetrahedr(?:on|a)|octahedr(?:on|a)|polyhedr(?:on|a))';
const distortedOf = (noun, extra = '') => new RegExp(`(^|[^\\w-])(distorted\\s+(?:\\w+\\s+)?${noun}|${noun}\\s+(?:are|is|appear|look|seem)\\s+`
    + `(?:(?!not\\b|never\\b|hardly\\b|barely\\b)\\w+\\s+)?(?:distorted|deformed)${extra})\\b`, 'gi');
export const DISTORTED = distortedOf(`(?:${CLUSTER_ONLY_NOUN}|${POLYHEDRON_NOUN})`, '|cluster\\s+distortions?');
export const DISTORTED_CLUSTERS = distortedOf(CLUSTER_ONLY_NOUN, '|cluster\\s+distortions?');
export const DISTORTED_POLYHEDRA = distortedOf(POLYHEDRON_NOUN, '|(?:polyhedral|octahedral|tetrahedral)\\s+distortions?');
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
// A clause about part of the shell (the first peak, the short bonds) is not a
// claim about the whole shell's count; nor is one that states a window in Å
// other than the whole shell's (reading.js windowsIn; facts.js).
export const PARTIAL_SHELL = /\bfirst[- ]peak\b|\bshort(?:est|er)?\s+bonds?\b/i;

// A window in Å, for the window check (facts.js; the coordination-shells
// interview, 2026-10-24, answer 8: every CN or angle carries its window in
// Å): a range "2.18–2.79 Å", "2.18 Å–2.79 Å", "from 2.18 to 2.79 Å", "between
// 2.18 and 2.79 Å", a filled reference "[2.18, 2.79] Å", a cutoff "within
// 3.2 Å" / "up to 3.2 Å" / "r < 3.2 Å", or a range after the word "window"
// (its unit understood). Each pattern's groups are the window's edges (one
// for a cutoff). An angle (°) or a share (%) is never a window.
const DISTANCE = '(\\d+(?:\\.\\d+)?)';
const ANGSTROM = '(?:Å|A\\b|angstroms?\\b)';
export const WINDOW_RANGES = Object.freeze([
    new RegExp(`${DISTANCE}\\s*(?:${ANGSTROM}\\s*)?(?:[–—−-]|to|and)\\s*${DISTANCE}\\s*${ANGSTROM}`, 'gi'),
    new RegExp(`\\[\\s*${DISTANCE}\\s*,\\s*${DISTANCE}\\s*\\]\\s*${ANGSTROM}`, 'gi'),
    new RegExp(`\\bwindows?\\b[^.;\\d°%]{0,40}?${DISTANCE}\\s*(?:[–—−-]|to)\\s*${DISTANCE}(?!\\s*[°%\\d])`, 'gi')
]);
export const WINDOW_CUTOFF = new RegExp(`(?:\\bwithin|\\bup to|\\bbelow|\\bunder|[<≤]|<=)\\s*${DISTANCE}\\s*${ANGSTROM}`, 'gi');
// A coordination number as a count of neighbours: "3.99 Se per Ga", "six Se
// around each Ta", "4 neighbours", "12 nearest Ga neighbours" (the count in
// group 2), beside CN_CLAIMS' "6-coordinated" and "CN 6". A polyhedron named
// in words ("octahedrally coordinated") names no number.
// "n A per B" with both elements: B's count of A (group 2 the count, 3 the
// partner A, 4 the centre B), which the coordination rule reads as such.
export const PER_COUNT = new RegExp(`(^|[^\\w.-])${COUNT}\\s+([A-Z][a-z]?)\\s+(?:atoms?\\s+)?(?:per|around(?:\\s+each)?|for each)\\s+([A-Z][a-z]?)(?![a-z])`, 'g');
export const NEIGHBOUR_COUNTS = Object.freeze([
    new RegExp(`(^|[^\\w.-])${COUNT}\\s+(?:[A-Z][a-z]?\\s+)?(?:atoms?\\s+)?(?:per|around(?:\\s+each)?|for each)\\s+[A-Z][a-z]?(?![a-z])`, 'g'),
    new RegExp(`(^|[^\\w.-])${COUNT}\\s+(?:(?:nearest|closest|first[- ]shell)\\s+)?(?:[A-Z][a-z]?\\s+)?neighbou?rs?\\b`, 'gi')
]);
// A bond angle: degrees in a clause that names an angle (bond angle, cis,
// trans, an A–B–C triplet), and not a direction's (a cell radius, a cone, an
// axis, a [u v w]) or a tilt; an ideal or reference value claims nothing.
export const DEGREES = /(^|[^\w.])(\d+(?:\.\d+)?)\s?°(?![CF]\b)/g;
export const BOND_ANGLE_WORDS = /\bangles?\b|\bcis\b|\btrans\b|(?<![A-Za-z])[A-Z][a-z]?\s*[–—−-]\s*[A-Z][a-z]?\s*[–—−-]\s*[A-Z][a-z]?(?![a-z])/;
export const NOT_A_BOND_ANGLE = /\bdirections?\b|\bradius\b|\bcone\b|\bax(?:is|es)\b|\bPC[123]\b|<\s*-?\d\s*-?\d\s*-?\d\s*>|\[\s*-?[\d.]+\s+-?[\d.]+\s+-?[\d.]+\s*\]|\btilts?\b|\bfamil(?:y|ies)\b|\bdeg\b/i;
export const REFERENCE_VALUE = /\bideal\b|\bperfect\b|\btheoretical\b|\breference\b/i;

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
