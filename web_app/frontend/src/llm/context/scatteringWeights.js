// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// How much of a total-scattering signal each partial pair carries — the
// Faber–Ziman weights w_ab = (2 − δ_ab) c_a c_b f_a f_b / (Σ c f)², which sum
// to 1 over the pairs a ≤ b. A partial with a small weight is barely seen by
// that dataset, so the RMC configuration's version of it is weakly constrained
// (Ga correlations under x-rays in GaTa4Se8: ~5 %). X-ray weights use f ≈ Z,
// the Q → 0 limit (heavier atoms gain relative weight at high Q); neutron
// weights use the bound coherent scattering lengths. The tables are duplicated
// from rmc_toolkits/scattering.py on purpose: this module may not import host
// code (see the llm README's import boundary).

// Bound coherent neutron scattering lengths (fm), NIST NCNR after Sears,
// Neutron News 3, 26 (1992); the real part for B, Cd, Dy, Eu, Gd, In, Sm.
// Kept equal to COHERENT_B_FM in rmc_toolkits/scattering.py.
export const COHERENT_B_FM = {
    Ag: 5.922, Al: 3.449, Am: 8.3, Ar: 1.909, As: 6.58, Au: 7.63,
    B: 5.3, Ba: 5.07, Be: 7.79, Bi: 8.532, Br: 6.795, C: 6.646,
    Ca: 4.7, Cd: 4.87, Ce: 4.84, Cl: 9.577, Co: 2.49, Cr: 3.635,
    Cs: 5.42, Cu: 7.718, Dy: 16.9, Er: 7.79, Eu: 7.22, F: 5.654,
    Fe: 9.45, Ga: 7.288, Gd: 6.5, Ge: 8.185, H: -3.739, He: 3.26,
    Hf: 7.7, Hg: 12.692, Ho: 8.01, I: 5.28, In: 4.065, Ir: 10.6,
    K: 3.67, Kr: 7.81, La: 8.24, Li: -1.9, Lu: 7.21, Mg: 5.375,
    Mn: -3.73, Mo: 6.715, N: 9.36, Na: 3.63, Nb: 7.054, Nd: 7.69,
    Ne: 4.566, Ni: 10.3, Np: 10.55, O: 5.803, Os: 10.7, P: 5.13,
    Pa: 9.1, Pb: 9.405, Pd: 5.91, Pm: 12.6, Pr: 4.58, Pt: 9.6,
    Ra: 10.0, Rb: 7.09, Re: 9.2, Rh: 5.88, Ru: 7.03, S: 2.847,
    Sb: 5.57, Sc: 12.29, Se: 7.97, Si: 4.1491, Sm: 0.8, Sn: 6.225,
    Sr: 7.02, Ta: 6.91, Tb: 7.38, Tc: 6.8, Te: 5.8, Th: 10.31,
    Ti: -3.438, Tl: 8.776, Tm: 7.07, U: 8.417, V: -0.3824, W: 4.86,
    Xe: 4.92, Y: 7.75, Yb: 12.43, Zn: 5.68, Zr: 7.16
};

// Element symbols in order of atomic number (Z = index + 1).
const SYMBOLS = (
    'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr '
    + 'Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb '
    + 'Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf'
).split(' ');

export const atomicNumber = (symbol) => {
    const index = SYMBOLS.indexOf(symbol);
    return index >= 0 ? index + 1 : null;
};

const amplitude = (radiation, element) => (
    radiation === 'xray' ? atomicNumber(element) : COHERENT_B_FM[element] ?? null
);

/**
 * Faber–Ziman weight of each partial (pairs a ≤ b in the given element
 * order) for 'xray' or 'neutron', from the element counts; null when an
 * element has no tabulated amplitude or the mean amplitude vanishes.
 * Returns [{ pair: 'Ta-Se', weight }] sorted by |weight|, largest first.
 */
export const partialWeights = (elementCounts, radiation) => {
    const elements = Object.keys(elementCounts || {}).filter((element) => elementCounts[element] > 0);
    const total = elements.reduce((sum, element) => sum + elementCounts[element], 0);
    if (!elements.length || !(total > 0)) return null;
    const amplitudes = elements.map((element) => amplitude(radiation, element));
    if (amplitudes.some((value) => !Number.isFinite(value))) return null;
    const fractions = elements.map((element) => elementCounts[element] / total);
    const mean = fractions.reduce((sum, c, i) => sum + c * amplitudes[i], 0);
    if (!(Math.abs(mean) > 1e-9)) return null;
    const weights = [];
    for (let i = 0; i < elements.length; i++) {
        for (let j = i; j < elements.length; j++) {
            const factor = i === j ? 1 : 2;
            weights.push({
                pair: `${elements[i]}-${elements[j]}`,
                weight: factor * fractions[i] * fractions[j] * amplitudes[i] * amplitudes[j] / (mean * mean)
            });
        }
    }
    return weights.sort((a, b) => Math.abs(b.weight) - Math.abs(a.weight));
};
