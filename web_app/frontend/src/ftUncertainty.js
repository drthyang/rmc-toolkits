// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// RMCProfile's Fourier transform of an X-ray F(Q), and its propagated
// uncertainty: a line-for-line port of rmc_toolkits/ft_uncertainty.py (the
// source of truth; its docstring has the derivation). Pure: no I/O, no DOM —
// ftUncertaintyRun.js reads a run's files and calls assessFtUncertainty.
//
//   G(r) = Σ_k ΔQ · Q_k · F(Q_k) · sin(Q_k r)
//
// over the run's F(Q) grid: a rectangle rule with no 2/π and no ρ₀ (Keen's
// D(r) times π/2), of the run's F(Q)_Expt plus its fitted offset. That column
// is the input data file renormalised by the fitted scale s and, with
// CONVOLVE, convolved with the box function of radius R (the last r of the
// run's PDFpartials). Two reproduction tests establish the chain to
// RMC_FT_TOLERANCE; only then is the input file's own σ propagated:
//
//   Cov[G(r), G(r')] = Σ_k (ΔQ Q_k s σ_k)² sin(Q_k r) sin(Q_k r') = ½[P(r − r') − P(r + r')]
//
// The parity golden is src/__tests__/fixtures/ft_uncertainty_fixture.json
// (tests/generate_ft_uncertainty_fixture.py → ftUncertaintyParity.test.js).

export const RMC_FT_TOLERANCE = 1e-6;
export const CONVOLUTION_CHECK_LIMIT = 0.01;
export const UNIFORM_GRID_RTOL = 1e-6;
export const CONVOLUTION_CHECK_POINTS = 8;
export const FT_FILE_RE = /^(.+)_FT_XFQ(\d+)\.csv$/;

const LINE_BREAK = /\r\n|\r|\n/;
const DAT_KEY_RE = /^([A-Z_\d]+)\s*::\s*(.*)$/;
const DAT_SUB_RE = /^([A-Z_]+)\s*(?:::\s*(.*))?$/;

// --- the transform and its uncertainty -------------------------------------

/** ΔQ of a uniform, strictly increasing grid, or a thrown Error (uniform_step). */
export const uniformStep = (q, name = 'Q') => {
    const n = q.length;
    if (n < 2) throw new Error(`the ${name} grid needs at least two points`);
    for (let k = 0; k < n; k += 1) {
        if (!Number.isFinite(q[k]) || (k > 0 && !(q[k] > q[k - 1]))) {
            throw new Error(`the ${name} grid is not strictly increasing`);
        }
    }
    const step = (q[n - 1] - q[0]) / (n - 1);
    let worst = 0;
    for (let k = 0; k < n; k += 1) worst = Math.max(worst, Math.abs(q[k] - (q[0] + step * k)));
    if (worst > UNIFORM_GRID_RTOL * step) throw new Error(`the ${name} grid is not uniform`);
    return step;
};

/** c_k = ΔQ · Q_k. */
export const rmcFtCoefficients = (q) => {
    const step = uniformStep(q);
    return Float64Array.from(q, (value) => step * value);
};

/** RMCProfile's transform of F(Q) onto `r`: Σ_k ΔQ Q_k F_k sin(Q_k r). */
export const rmcFt = (q, f, r) => {
    const c = rmcFtCoefficients(q);
    return Float64Array.from(r, (at) => {
        let sum = 0;
        for (let k = 0; k < q.length; k += 1) sum += Math.sin(at * q[k]) * c[k] * f[k];
        return sum;
    });
};

const powerOf = (q, sigma) => {
    const c = rmcFtCoefficients(q);
    return Float64Array.from(c, (value, k) => (value * sigma[k]) ** 2);
};

/** σ of the transform at each r for independent σ_k: √Σ_k (c_k σ_k sin Q_k r)². */
export const rmcFtSigma = (q, sigma, r) => {
    const power = powerOf(q, sigma);
    return Float64Array.from(r, (at) => {
        let sum = 0;
        for (let k = 0; k < q.length; k += 1) {
            const s = Math.sin(at * q[k]);
            sum += s * s * power[k];
        }
        return Math.sqrt(sum);
    });
};

/** P(x) = Σ_k (c_k σ_k)² cos(Q_k x), the covariance's Toeplitz kernel. */
export const noiseAutocorrelation = (q, sigma, x) => {
    const power = powerOf(q, sigma);
    return Float64Array.from(x, (at) => {
        let sum = 0;
        for (let k = 0; k < q.length; k += 1) sum += Math.cos(at * q[k]) * power[k];
        return sum;
    });
};

const uniformRStep = (r) => {
    const n = r.length;
    if (n < 2) return null;
    const step = (r[n - 1] - r[0]) / (n - 1);
    if (!(step > 0)) return null;
    let worst = 0;
    for (let i = 0; i < n; i += 1) worst = Math.max(worst, Math.abs(r[i] - (r[0] + step * i)));
    return worst <= 1e-7 * step ? step : null;
};

/**
 * Cov[G(r_i), G(r_j)] as an entry accessor (i, j) → value, without the n × n
 * matrix: uniform r from the O(n) Toeplitz/Hankel values of P, otherwise the
 * direct sum. The diagonal is rmcFtSigma squared, the same numbers.
 */
export const covarianceEntries = (q, sigma, r) => {
    const n = r.length;
    const diagonal = Float64Array.from(rmcFtSigma(q, sigma, r), (value) => value * value);
    const step = uniformRStep(r);
    if (step !== null) {
        const toeplitz = noiseAutocorrelation(q, sigma, Float64Array.from({ length: n }, (_, d) => step * d));
        const hankel = noiseAutocorrelation(q, sigma, Float64Array.from({ length: Math.max(0, 2 * n - 1) }, (_, s) => 2 * r[0] + step * s));
        return (i, j) => (i === j ? diagonal[i] : 0.5 * (toeplitz[Math.abs(i - j)] - hankel[i + j]));
    }
    const power = powerOf(q, sigma);
    const nq = q.length;
    const sines = new Float64Array(n * nq);
    for (let i = 0; i < n; i += 1) for (let k = 0; k < nq; k += 1) sines[i * nq + k] = Math.sin(r[i] * q[k]);
    return (i, j) => {
        if (i === j) return diagonal[i];
        let sum = 0;
        for (let k = 0; k < nq; k += 1) sum += sines[i * nq + k] * power[k] * sines[j * nq + k];
        return sum;
    };
};

/** The full covariance, row-major n × n (rmc_ft_covariance). */
export const rmcFtCovariance = (q, sigma, r) => {
    const n = r.length;
    const entry = covarianceEntries(q, sigma, r);
    const out = new Float64Array(n * n);
    for (let i = 0; i < n; i += 1) {
        out[i * n + i] = entry(i, i);
        for (let j = i + 1; j < n; j += 1) {
            const value = entry(i, j);
            out[i * n + j] = value;
            out[j * n + i] = value;
        }
    }
    return out;
};

/** Neighbour and Nyquist-step correlations, effective and Nyquist point counts (correlation_summary). */
export const correlationSummary = (q, sigma, r) => {
    const n = r.length;
    const gridStep = n > 1 ? (r[n - 1] - r[0]) / (n - 1) : 0;
    const qmax = q[q.length - 1];
    const nyquistStep = Math.PI / qmax;
    const [p0, pGrid, pNyquist] = noiseAutocorrelation(q, sigma, [0, gridStep, nyquistStep]);
    const entry = covarianceEntries(q, sigma, r);
    const sd = Float64Array.from({ length: n }, (_, i) => Math.sqrt(entry(i, i)));
    let m = 0;
    let frobenius = 0;
    for (let i = 0; i < n; i += 1) {
        if (!(sd[i] > 0)) continue;
        m += 1;
        for (let j = 0; j < n; j += 1) {
            if (!(sd[j] > 0)) continue;
            const rho = i === j ? entry(i, i) / (sd[i] * sd[i]) : entry(i, j) / (sd[i] * sd[j]);
            frobenius += rho * rho;
        }
    }
    return {
        points: n,
        gridStep,
        nyquistStep,
        neighbour: p0 > 0 ? pGrid / p0 : 0,
        nyquist: p0 > 0 ? pNyquist / p0 : 0,
        effectivePoints: m && frobenius > 0 ? (m * m) / frobenius : 0,
        nyquistPoints: n ? (Math.max(0, r[n - 1] - r[0]) * qmax) / Math.PI : 0
    };
};

const boxKernelTerms = (q, radius) => {
    const step = uniformStep(q);
    if (!(q[0] > 0)) throw new Error('the box convolution needs Q > 0 (it divides by Q)');
    const n = q.length;
    const toeplitz = new Float64Array(n);
    toeplitz[0] = radius;
    for (let m = 1; m < n; m += 1) {
        const lag = step * m;
        toeplitz[m] = Math.sin(lag * radius) / lag;
    }
    const hankel = new Float64Array(2 * n - 1);
    for (let s = 0; s < 2 * n - 1; s += 1) {
        const total = 2 * q[0] + step * s;
        hankel[s] = Math.sin(total * radius) / total;
    }
    return { toeplitz, hankel };
};

// Σ_k v_k [T(|k − j|) − H(k + j)] for every j.
const boxKernelApply = (v, { toeplitz, hankel }) => {
    const n = v.length;
    const out = new Float64Array(n);
    for (let j = 0; j < n; j += 1) {
        let t = 0;
        let h = 0;
        for (let k = 0; k < n; k += 1) {
            t += v[k] * toeplitz[Math.abs(k - j)];
            h += hankel[k + j] * v[k];
        }
        out[j] = t - h;
    }
    return out;
};

/** RMCProfile's CONVOLVE: (K F)(Q_k) = 1/(π Q_k) Σ_j ΔQ Q_j F_j [T(|k − j|) − H(k + j)]. */
export const boxConvolution = (q, f, radius) => {
    const terms = boxKernelTerms(q, radius);
    const c = rmcFtCoefficients(q);
    const applied = boxKernelApply(Float64Array.from(c, (value, k) => value * f[k]), terms);
    return Float64Array.from(applied, (value, k) => value / (Math.PI * q[k]));
};

/** r_min, and CONVOLUTION_CHECK_POINTS points over the last π/Q_max up to r_max. */
export const convolutionCheckPoints = (r, qmax) => {
    const rMin = r[0];
    const rMax = r[r.length - 1];
    const span = Math.min(Math.PI / qmax, rMax - rMin);
    const count = CONVOLUTION_CHECK_POINTS;
    const tail = [];
    for (let k = count - 1; k >= 0; k -= 1) tail.push(rMax - (span * k) / (count - 1));
    return Float64Array.from([rMin, ...tail]);
};

/** σ(r) with the box convolution over σ(r) without it, minus 1, at each point (convolution_sigma_change). */
export const convolutionSigmaChange = (q, sigma, rPoints, radius) => {
    const step = uniformStep(q);
    const terms = boxKernelTerms(q, radius);
    return Float64Array.from(rPoints, (at) => {
        const sines = Float64Array.from(q, (value) => Math.sin(value * at));
        const applied = boxKernelApply(sines, terms);
        let withBox = 0;
        let without = 0;
        for (let j = 0; j < q.length; j += 1) {
            const exact = ((step * step * q[j]) / Math.PI) * applied[j];
            withBox += (exact * sigma[j]) ** 2;
            without += (step * q[j] * sines[j] * sigma[j]) ** 2;
        }
        return without > 0 ? Math.sqrt(withBox / without) - 1 : 0;
    });
};

// --- the two reproduction tests --------------------------------------------

const deviation = (target, model, tolerance) => {
    let scale = 0;
    let largest = 0;
    let squares = 0;
    for (let i = 0; i < target.length; i += 1) {
        const residual = target[i] - model[i];
        scale = Math.max(scale, Math.abs(target[i]));
        largest = Math.max(largest, Math.abs(residual));
        squares += residual * residual;
    }
    const relative = scale > 0 ? largest / scale : (largest === 0 ? 0 : Infinity);
    return { residual: relative, rms: Math.sqrt(squares / target.length), pass: relative <= tolerance };
};

/** G ≈ A(F + c), c by least squares (reproduce_transform): { residual, rms, pass, offset }. */
export const reproduceTransform = (q, f, r, g, tolerance = RMC_FT_TOLERANCE) => {
    const c = rmcFtCoefficients(q);
    const base = new Float64Array(r.length);
    const unit = new Float64Array(r.length);
    for (let i = 0; i < r.length; i += 1) {
        let b = 0;
        let u = 0;
        for (let k = 0; k < q.length; k += 1) {
            const s = Math.sin(r[i] * q[k]);
            b += s * c[k] * f[k];
            u += s * c[k];
        }
        base[i] = b;
        unit[i] = u;
    }
    let num = 0;
    let norm = 0;
    for (let i = 0; i < r.length; i += 1) {
        num += (g[i] - base[i]) * unit[i];
        norm += unit[i] * unit[i];
    }
    const offset = norm > 0 ? num / norm : 0;
    const model = Float64Array.from(base, (value, i) => value + offset * unit[i]);
    return { ...deviation(g, model, tolerance), offset };
};

/** F_run ≈ s·(K F_input) + c' (reproduce_input): { residual, rms, pass, scale, offset }. */
export const reproduceInput = (q, fRun, fInput, { radius = null, tolerance = RMC_FT_TOLERANCE } = {}) => {
    const base = radius !== null ? boxConvolution(q, fInput, radius) : Float64Array.from(fInput);
    const n = base.length;
    let sx = 0;
    let sy = 0;
    let sxx = 0;
    let sxy = 0;
    for (let i = 0; i < n; i += 1) {
        sx += base[i];
        sy += fRun[i];
        sxx += base[i] * base[i];
        sxy += base[i] * fRun[i];
    }
    const det = n * sxx - sx * sx;
    const scale = det === 0 ? 0 : (n * sxy - sx * sy) / det;
    const offset = det === 0 ? sy / n : (sxx * sy - sx * sxy) / det;
    const model = Float64Array.from(base, (value) => scale * value + offset);
    return { ...deviation(fRun, model, tolerance), scale, offset };
};

// --- the run-control file --------------------------------------------------

const NUMBER_RE = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;

// The plain decimal numbers a value starts with (`1 2619 1 Angstrom` → [1, 2619, 1]).
const leadingNumbers = (value) => {
    const numbers = [];
    for (const token of value.trim().split(/\s+/)) {
        if (!NUMBER_RE.test(token)) break;
        numbers.push(Number(token));
    }
    return numbers;
};

/** The .dat's data blocks, as the transform needs them (run_datasets). */
export const runDatasets = (text) => {
    const datasets = [];
    let block = null;
    let rSpacing = null;
    let numberDensity = null;
    for (const raw of text.split(LINE_BREAK)) {
        const line = raw.trim();
        if (!line) continue;
        if (line.startsWith('>')) {
            if (!block) continue;
            const match = DAT_SUB_RE.exec(line.slice(1).trim());
            if (!match) continue;
            const key = match[1];
            const value = (match[2] || '').trim();
            if (key === 'FILENAME') block.file = value;
            else if (key === 'DATA_TYPE') block.dataType = value;
            else if (key === 'FIT_TYPE') block.fitType = value;
            else if (key === 'CONVOLVE') block.convolve = true;
            else if (key === 'FITTED_SCALE') block.fittedScale = true;
            else if (key === 'FITTED_OFFSET') block.fittedOffset = true;
            else if (key === 'RECIPROCAL_SPACE_PARAMETERS') block.reciprocalParameters = leadingNumbers(value);
            else if (key === 'REAL_SPACE_PARAMETERS') block.realParameters = leadingNumbers(value);
            continue;
        }
        const match = DAT_KEY_RE.exec(line);
        block = null;
        if (!match) continue;
        const [, key, value] = match;
        if (key === 'END') break;
        if (key.endsWith('_DATA')) {
            block = {
                block: key, file: null, dataType: null, fitType: null, convolve: false,
                fittedScale: false, fittedOffset: false, reciprocalParameters: null, realParameters: null
            };
            datasets.push(block);
        } else if (key === 'R_SPACING') {
            const numbers = leadingNumbers(value);
            rSpacing = numbers.length ? numbers[0] : null;
        } else if (key === 'NUMBER_DENSITY') {
            const numbers = leadingNumbers(value);
            numberDensity = numbers.length ? numbers[0] : null;
        }
    }
    return { rSpacing, numberDensity, datasets };
};

/** The `number`-th (1-based) XRAY_RECIPROCAL_SPACE_DATA block, or null. */
export const xrayDataset = (settings, number) => {
    const blocks = (settings?.datasets || []).filter((entry) => entry.block === 'XRAY_RECIPROCAL_SPACE_DATA');
    return number > 0 && number <= blocks.length ? blocks[number - 1] : null;
};

/** A header that names an uncertainty: sigma, σ, uncert…, error, or the word err / esd. */
export const isSigmaLabel = (label) => {
    const text = String(label ?? '').toLowerCase();
    if (/sigma|σ|uncert|error/.test(text)) return true;
    return text.split(/[^a-z0-9σ]+/).some((word) => word === 'err' || word === 'esd');
};

/** The index of a column whose header names an uncertainty (never the x column), or null. */
export const fileSigmaColumn = (labels) => {
    for (let index = 1; index < labels.length; index += 1) {
        if (isSigmaLabel(labels[index])) return index;
    }
    return null;
};

// --- the assessment --------------------------------------------------------

const usableSigma = (q, f, sigma) => {
    for (let i = 0; i < q.length; i += 1) {
        if (Number.isFinite(q[i]) && Number.isFinite(f[i]) && !(Number.isFinite(sigma[i]) && sigma[i] > 0)) return null;
    }
    return sigma;
};

const matchingBlock = (qInput, q, step) => {
    const tolerance = UNIFORM_GRID_RTOL * step;
    for (let start = 0; start < qInput.length; start += 1) {
        if (!(Math.abs(qInput[start] - q[0]) <= tolerance)) continue;
        if (start + q.length > qInput.length) continue;
        let ok = true;
        for (let k = 0; k < q.length && ok; k += 1) ok = Math.abs(qInput[start + k] - q[k]) <= tolerance;
        if (ok) return [start, start + q.length];
    }
    return null;
};

const allFinite = (values) => {
    for (let i = 0; i < values.length; i += 1) if (!Number.isFinite(values[i])) return false;
    return true;
};

/**
 * Whether, and how, the G(r) file's error ribbon follows from the run's files
 * (assess_ft_uncertainty, the same arguments camel-cased: { r, g, q, f,
 * settings, datasetNumber, inputColumns, radius, names, tolerance }).
 * Returns { status, reason, sigma, transform, reproduction, convolutionCheck,
 * correlation, sources }.
 */
export const assessFtUncertainty = ({
    r, g, q = null, f = null, settings = null, datasetNumber = 1, inputColumns = null, radius = null,
    names = {}, tolerance = RMC_FT_TOLERANCE
}) => {
    const transformName = names.transform || 'the transform file';
    const fqName = names.fq;
    const settingsName = names.settings || 'the run-control .dat';
    const missingSettings = names.settings
        ? `no run-control file ${names.settings} beside ${transformName}`
        : `no run-control .dat beside ${transformName}`;
    const inputName = names.input;
    const partialsName = names.partials || "the run's PDFpartials file";
    const fqCandidates = names.fqCandidates || "the run's F(Q) file";
    const used = {
        transform: true, fq: q !== null, settings: settings !== null, input: inputColumns !== null, partials: radius !== null
    };
    const sources = {};
    Object.entries(used).forEach(([key, present]) => { if (present && names[key]) sources[key] = names[key]; });

    const payload = {
        status: 'unavailable',
        reason: null,
        sigma: null,
        transform: null,
        reproduction: { tolerance, gr: null, fq: null },
        convolutionCheck: null,
        correlation: null,
        sources
    };

    const dataset = settings !== null ? xrayDataset(settings, datasetNumber) : null;
    let sigmaInput = null;
    let noSigma;
    if (settings === null) noSigma = missingSettings;
    else if (dataset === null) noSigma = `${settingsName} names no X-ray reciprocal-space dataset ${datasetNumber}`;
    else if (!dataset.file) noSigma = `${settingsName} names no data file for its X-ray dataset ${datasetNumber}`;
    else if (inputColumns === null) {
        noSigma = `F(Q) carries no uncertainties: ${dataset.file} (named in ${settingsName}) is not in the run folder`;
    } else if (inputColumns.length < 3) noSigma = `F(Q) carries no uncertainties: ${inputName} has no error column`;
    else {
        sigmaInput = usableSigma(inputColumns[0], inputColumns[1], inputColumns[2]);
        noSigma = sigmaInput !== null ? null
            : `F(Q) carries no uncertainties: the error column of ${inputName} is not usable (a zero, negative or non-finite σ)`;
    }
    const convolve = Boolean(dataset && dataset.convolve);

    const problems = [];
    let increasing = r.length >= 2;
    for (let i = 0; i < r.length && increasing; i += 1) {
        increasing = Number.isFinite(r[i]) && (i === 0 || r[i] > r[i - 1]);
    }
    if (!increasing) problems.push(`the r column of ${transformName} is not strictly increasing`);
    if (!allFinite(g)) problems.push(`${transformName} has non-finite values`);

    let step = null;
    if (q !== null) {
        try {
            step = uniformStep(q);
        } catch (error) {
            problems.push(`${fqName}: ${error.message}`);
        }
        if (!allFinite(f)) problems.push(`${fqName} has non-finite F(Q) values`);
        if (step !== null && !problems.length) {
            payload.transform = {
                qMin: q[0], qMax: q[q.length - 1], dQ: step, nQ: q.length,
                rMin: r[0], rMax: r[r.length - 1], nR: r.length,
                convolve, radius: convolve && radius !== null ? radius : null
            };
            payload.reproduction.gr = reproduceTransform(q, f, r, g, tolerance);
        }
    }

    if (payload.reproduction.gr !== null && inputColumns !== null && inputColumns.length >= 2) {
        const block = matchingBlock(inputColumns[0], q, step);
        if (block === null) {
            problems.push(`the Q grid of ${inputName} does not hold the Q grid of ${fqName}`);
        } else {
            const fInput = inputColumns[1].slice(block[0], block[1]);
            if (!allFinite(fInput)) {
                problems.push(`${inputName} has non-finite F(Q) values on the run's Q grid`);
            } else if (!(convolve && radius === null)) {
                payload.reproduction.fq = {
                    ...reproduceInput(q, f, fInput, { radius: convolve ? radius : null, tolerance }),
                    convolve
                };
            }
            if (sigmaInput !== null) sigmaInput = sigmaInput.slice(block[0], block[1]);
        }
    }

    const { gr, fq } = payload.reproduction;
    if (noSigma !== null) return { ...payload, reason: noSigma };
    if (problems.length) return { ...payload, status: 'failed', reason: problems[0] };
    if (q === null) return { ...payload, reason: `no F(Q) file (${fqCandidates}) to reproduce the transform from` };
    if (convolve && radius === null) {
        return {
            ...payload,
            reason: `${partialsName} is not in the run folder: the radius of the box convolution RMCProfile applied is unknown`
        };
    }
    if (!gr.pass) {
        return {
            ...payload, status: 'failed',
            reason: `RMCProfile's transform of ${fqName} does not reproduce ${transformName} within the tolerance`
        };
    }
    if (!fq.pass) {
        const how = convolve ? ' and box convolution' : '';
        return {
            ...payload, status: 'failed',
            reason: `${inputName} does not reproduce ${fqName} through the run's renormalisation${how} within the tolerance`
        };
    }
    if (!(fq.scale > 0)) {
        return { ...payload, status: 'failed', reason: `the renormalisation of ${inputName} the run fitted is not positive` };
    }

    const sigmaQ = Float64Array.from(sigmaInput, (value) => fq.scale * value);
    if (convolve) {
        const points = convolutionCheckPoints(r, q[q.length - 1]);
        const changes = convolutionSigmaChange(q, sigmaQ, points, radius);
        let worst = 0;
        for (let i = 1; i < changes.length; i += 1) if (Math.abs(changes[i]) > Math.abs(changes[worst])) worst = i;
        payload.convolutionCheck = {
            radius,
            limit: CONVOLUTION_CHECK_LIMIT,
            maxChange: Math.abs(changes[worst]),
            at: points[worst],
            points: points.length
        };
        if (Math.abs(changes[worst]) > CONVOLUTION_CHECK_LIMIT) {
            return {
                ...payload, status: 'failed',
                reason: `the r range of ${transformName} reaches too close to the box-convolution radius: left out of the `
                    + 'propagation, the convolution would change σ(r) beyond the limit'
            };
        }
    }

    return {
        ...payload,
        status: 'propagated',
        reason: null,
        sigma: Array.from(rmcFtSigma(q, sigmaQ, r)),
        correlation: correlationSummary(q, sigmaQ, r)
    };
};
