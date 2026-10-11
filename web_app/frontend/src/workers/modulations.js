// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Displacement modulations — static-mode port of rmc_toolkits/modulations.py
// (the source of truth; keep the two in sync, the parity fixture pins them).
//
// The average structure folds every cell into one, so it cannot show order
// that enlarges the cell. The configuration can: an ordered pattern with
// wavevector q puts power into the Fourier component, at q, of every site's
// displacement field over the supercell's cell grid. Per site the power
// |U_s(q)|² / (n_s tr C_s) is 1 on average for independent displacements;
// an element's sites are averaged; the zone-boundary points and the
// strongest modes are reported against the empirical distribution over all
// q ≠ 0 (ratio to its median, percentile). numpy's fftn is replaced by a
// separable DFT — the same numbers to floating-point noise.

export const TOP_MODES = 5;
export const PERCENTILE_RTOL = 1e-9;

/** The half-integer points (0/1 flags per axis) an even supercell can hold. */
export const zoneBoundaryPoints = (supercell) => {
    const points = [];
    for (const a of [0, 1]) {
        for (const b of [0, 1]) {
            for (const c of [0, 1]) {
                const flags = [a, b, c];
                if (!flags.some(Boolean)) continue;
                if (flags.every((flag, axis) => !flag || supercell[axis] % 2 === 0)) points.push(flags);
            }
        }
    }
    return points;
};

const folded = (index, size) => {
    const value = index / size;
    return value > 0.5 ? value - 1 : value;
};

// In-place DFT of `re`/`im` (length n·stride… strided lines) along one axis
// of a C-ordered [N0][N1][N2] grid: X_k = Σ_n x_n e^{−2πi kn/N}.
const dftAxis = (re, im, dims, axis) => {
    const n = dims[axis];
    const strides = [dims[1] * dims[2], dims[2], 1];
    const stride = strides[axis];
    const cos = Float64Array.from({ length: n * n }, (_, i) => Math.cos(-2 * Math.PI * ((Math.floor(i / n) * (i % n)) % n) / n));
    const sin = Float64Array.from({ length: n * n }, (_, i) => Math.sin(-2 * Math.PI * ((Math.floor(i / n) * (i % n)) % n) / n));
    const lineRe = new Float64Array(n);
    const lineIm = new Float64Array(n);
    const total = dims[0] * dims[1] * dims[2];
    for (let start = 0; start < total; start += 1) {
        // A line starts where this axis's index is 0.
        if (Math.floor(start / stride) % n !== 0) continue;
        for (let k = 0; k < n; k += 1) {
            lineRe[k] = re[start + k * stride];
            lineIm[k] = im[start + k * stride];
        }
        for (let k = 0; k < n; k += 1) {
            let sumRe = 0;
            let sumIm = 0;
            for (let j = 0; j < n; j += 1) {
                const c = cos[k * n + j];
                const s = sin[k * n + j];
                sumRe += lineRe[j] * c - lineIm[j] * s;
                sumIm += lineRe[j] * s + lineIm[j] * c;
            }
            re[start + k * stride] = sumRe;
            im[start + k * stride] = sumIm;
        }
    }
};

const capitalize = (symbol) => (symbol ? symbol[0].toUpperCase() + symbol.slice(1).toLowerCase() : symbol);

/**
 * Per-element displacement power at every wavevector of the supercell grid,
 * from a parsed configuration (siteDisplacementsFromRmc6f). Mirrors
 * rmc_toolkits.modulations.displacement_modulations; camelCase payload.
 */
export const displacementModulations = (parsed, element = null, { top = TOP_MODES } = {}) => {
    const sites = parsed.sites;
    if (parsed.reconstructed || sites.some((site) => !Array.isArray(site.cells))) {
        throw new Error('modulation analysis needs the reference-site and cell columns of the .rmc6f file');
    }
    const supercell = parsed.supercell.map(Number);
    const nQ = supercell[0] * supercell[1] * supercell[2];
    if (nQ < 8) throw new Error(`a [${supercell.join(', ')}] supercell has too few wavevectors to compare`);
    const labels = sites.map((site) => site.element);
    let elements;
    if (element != null) {
        const symbol = capitalize(String(element).trim());
        if (!labels.includes(symbol)) {
            const species = [...new Set(sites.flatMap((site) => Object.keys(site.elementCounts || { [site.element]: 1 })))].sort();
            throw new Error(`no site holds ${symbol}; elements: ${species.join(', ')}`);
        }
        elements = [symbol];
    } else {
        elements = [...new Set(labels)].sort();
    }
    const boundary = zoneBoundaryPoints(supercell);
    const untestable = [];
    for (const a of [0, 1]) for (const b of [0, 1]) for (const c of [0, 1]) {
        const flags = [a, b, c];
        if (flags.some(Boolean) && !boundary.some((p) => p.every((v, i) => v === flags[i]))) untestable.push(flags);
    }
    const result = {
        supercell,
        wavevectors: nQ - 1,
        untestable,
        elements: [],
        parseWarning: parsed.parseWarning ?? null
    };
    const linear = ([i, j, k]) => {
        const w = [i, j, k].map((v, axis) => ((v % supercell[axis]) + supercell[axis]) % supercell[axis]);
        return (w[0] * supercell[1] + w[1]) * supercell[2] + w[2];
    };
    for (const name of elements) {
        const power = new Float64Array(nQ);
        let missing = 0;
        let used = 0;
        sites.forEach((site) => {
            if (site.element !== name) return;
            const u = site.displacements;
            const n = u.length;
            const mean = [0, 1, 2].map((k) => u.reduce((sum, row) => sum + row[k], 0) / n);
            const variance = [0, 1, 2].reduce((sum, k) => sum + u.reduce((acc, row) => acc + (row[k] - mean[k]) ** 2, 0) / n, 0);
            if (!(variance > 0)) return;
            missing += nQ - n;
            for (let component = 0; component < 3; component += 1) {
                const re = new Float64Array(nQ);
                const im = new Float64Array(nQ);
                site.cells.forEach((cell, row) => { re[linear(cell)] = u[row][component]; });
                for (let axis = 0; axis < 3; axis += 1) dftAxis(re, im, supercell, axis);
                for (let q = 0; q < nQ; q += 1) power[q] += (re[q] * re[q] + im[q] * im[q]) / (n * variance);
            }
            used += 1;
        });
        if (!used) {
            result.elements.push({ element: name, sites: 0 });
            continue;
        }
        for (let q = 0; q < nQ; q += 1) power[q] /= used;
        const others = Array.from(power.subarray(1));
        const sorted = [...others].sort((a, b) => a - b);
        const mid = sorted.length >> 1;
        const median = sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
        const unravel = (q) => [Math.floor(q / (supercell[1] * supercell[2])), Math.floor(q / supercell[2]) % supercell[1], q % supercell[2]];
        const describe = (q) => {
            const value = power[q];
            const index = unravel(q);
            let lower = 0;
            others.forEach((other) => { if (other < value * (1 - PERCENTILE_RTOL)) lower += 1; });
            return {
                q: index.map((v, axis) => folded(v, supercell[axis])),
                ratio: median > 0 ? value / median : null,
                percentile: (lower / others.length) * 100
            };
        };
        const zoneBoundary = boundary.map((flags) => {
            const entry = describe(linear(flags.map((flag, axis) => (flag ? supercell[axis] / 2 : 0))));
            entry.q = flags.map((flag) => 0.5 * flag);
            return entry;
        });
        const order = Array.from({ length: nQ }, (_, q) => q).sort((a, b) => (power[b] - power[a]) || (a - b));
        const strongest = [];
        const seen = new Set();
        for (const q of order) {
            if (q === 0) continue;
            const index = unravel(q);
            const partner = index.map((v, axis) => ((-v % supercell[axis]) + supercell[axis]) % supercell[axis]);
            // +q and −q carry the same power (to rounding): report the pair
            // once, as its member with the smaller grid index, in both engines.
            const key = linear(index) < linear(partner) ? linear(index) : linear(partner);
            if (seen.has(key)) continue;
            seen.add(key);
            strongest.push(describe(key));
            if (strongest.length === top) break;
        }
        result.elements.push({ element: name, sites: used, missingCopies: missing, medianPower: median, zoneBoundary, strongest });
    }
    return result;
};
