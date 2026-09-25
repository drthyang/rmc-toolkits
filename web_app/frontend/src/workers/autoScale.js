// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

/**
 * Static-mode Auto StoG engine: a straight JS port of rmc_toolkits/transforms.py
 * + rmc_toolkits/scaling.py (+ the stog parsers and the Faber-Ziman calculator),
 * so the browser can auto-scale S(Q) without a backend. Keep the two in sync —
 * src/__tests__/autoScale.test.js asserts parity against Python-generated golden
 * numbers (tests/generate_autoscale_fixture.py).
 *
 * All math is float64 on plain arrays / Float64Array; discretization (trapezoid
 * sine transforms on the data grids) is identical to the Python engine's.
 */

// ---------------------------------------------------------------------------
// Small numerics
// ---------------------------------------------------------------------------

const isNum = (value) => Number.isFinite(value);

const median = (values) => {
  const sorted = Array.from(values).sort((a, b) => a - b);
  const n = sorted.length;
  if (!n) return NaN;
  const mid = Math.floor(n / 2);
  return n % 2 ? sorted[mid] : 0.5 * (sorted[mid - 1] + sorted[mid]);
};

const mean = (values) => {
  let total = 0;
  for (let i = 0; i < values.length; i += 1) total += values[i];
  return total / values.length;
};

/** Least squares for a small (<= 3 column) design via pivoted normal equations. */
const solveLeastSquares = (columns, rhs) => {
  const k = columns.length;
  const ata = Array.from({ length: k }, () => new Float64Array(k));
  const atb = new Float64Array(k);
  for (let i = 0; i < k; i += 1) {
    for (let j = i; j < k; j += 1) {
      let dot = 0;
      for (let row = 0; row < rhs.length; row += 1) dot += columns[i][row] * columns[j][row];
      ata[i][j] = dot;
      ata[j][i] = dot;
    }
    let dot = 0;
    for (let row = 0; row < rhs.length; row += 1) dot += columns[i][row] * rhs[row];
    atb[i] = dot;
  }
  // Gaussian elimination with partial pivoting on the k x k system.
  const a = ata.map((rowArr, i) => [...rowArr, atb[i]]);
  for (let col = 0; col < k; col += 1) {
    let pivot = col;
    for (let row = col + 1; row < k; row += 1) {
      if (Math.abs(a[row][col]) > Math.abs(a[pivot][col])) pivot = row;
    }
    [a[col], a[pivot]] = [a[pivot], a[col]];
    const diag = a[col][col];
    if (Math.abs(diag) < 1e-300) throw new Error('singular least-squares system');
    for (let row = 0; row < k; row += 1) {
      if (row === col) continue;
      const factor = a[row][col] / diag;
      for (let j = col; j <= k; j += 1) a[row][j] -= factor * a[col][j];
    }
  }
  return Array.from({ length: k }, (_, i) => a[i][k] / a[i][i]);
};

// ---------------------------------------------------------------------------
// Transforms (port of rmc_toolkits/transforms.py)
// ---------------------------------------------------------------------------

/**
 * Throw unless x is strictly increasing and finite (transforms._increasing_grid):
 * the trapezoid rule integrates with signed panel widths, so a descending grid
 * silently negates every integral. Fewer than two points (an empty filter
 * section) integrate to 0 and pass.
 */
const requireIncreasing = (x, name) => {
  if (x.length < 2) return;
  let ok = Number.isFinite(x[0]);
  for (let i = 1; ok && i < x.length; i += 1) ok = Number.isFinite(x[i]) && x[i] > x[i - 1];
  if (!ok) {
    throw new Error(
      `the ${name} grid must be strictly increasing (and finite); sort it `
      + 'and merge or remove duplicate points first'
    );
  }
};

export const sineTransform = (x, y, xout) => {
  requireIncreasing(x, 'integration (x)');
  const out = new Float64Array(xout.length);
  const n = x.length;
  for (let k = 0; k < xout.length; k += 1) {
    const w = xout[k];
    let total = 0;
    let prev = y[0] * Math.sin(x[0] * w);
    for (let i = 1; i < n; i += 1) {
      const current = y[i] * Math.sin(x[i] * w);
      total += (x[i] - x[i - 1]) * (current + prev) * 0.5;
      prev = current;
    }
    out[k] = total;
  }
  return out;
};

/** sin(v)/v with the limit 1 at v = 0 (transforms._sinc). */
const sinc = (v) => (v === 0 ? 1 : Math.sin(v) / v);
/** (v sin v + cos v - 1)/v^2 = sinc(v) - sinc(v/2)^2 / 2, cancellation-free (transforms._sinc_head). */
const sincHead = (v) => sinc(v) - 0.5 * sinc(0.5 * v) ** 2;

export const lorchWindow = (q, qmax) => {
  const out = new Float64Array(q.length);
  for (let i = 0; i < q.length; i += 1) {
    const x = (Math.PI * q[i]) / qmax;
    out[i] = x === 0 ? 1 : Math.sin(x) / x;
  }
  return out;
};

export const lowQCorrectionBasis = (q, r, { lorch = false, s0Target = 0 } = {}) => {
  requireIncreasing(q, 'Q');
  const coef = new Float64Array(r.length);
  const constant = new Float64Array(r.length);
  const q0 = q[0];
  if (q0 === 0) return { coef, constant };
  const scale = 2 / Math.PI;
  if (lorch) {
    // Cancellation-free form (transforms.low_q_correction_basis): with
    // v = q0 (r -/+ a), sin(v)/(r-a) = q0 sinc(v) and
    // (v sin v + cos v - 1)/(r-a)^2 = q0^2 [sinc(v) - sinc(v/2)^2 / 2]
    // (cos v - 1 = -2 sin^2(v/2)); no removable singularity at r = a, no patch.
    const a = Math.PI / q[q.length - 1];
    for (let i = 0; i < r.length; i += 1) {
      const vm = q0 * (r[i] - a);
      const vp = q0 * (r[i] + a);
      const f1 = (q0 * q0 * (sincHead(vm) - sincHead(vp))) / (2 * a);
      const f2 = (q0 * (sinc(vm) - sinc(vp))) / (2 * a);
      coef[i] = scale * (f1 / q0);
      constant[i] = scale * f2;
    }
  } else {
    for (let i = 0; i < r.length; i += 1) {
      const ri = r[i];
      if (ri === 0) continue;
      const v = q0 * ri;
      const f1 = (2 * v * Math.sin(v) - (v * v - 2) * Math.cos(v) - 2) / (ri * ri * ri);
      const f2 = (Math.sin(v) - v * Math.cos(v)) / (ri * ri);
      coef[i] = scale * (f1 / q0);
      constant[i] = scale * f2;
    }
  }
  if (s0Target !== 0) {
    // S(Q) extrapolates linearly to S(0) = s0Target instead of pystog's 0:
    // only the constant term changes (const' = (1-s0)*const + s0*coef).
    for (let i = 0; i < r.length; i += 1) {
      constant[i] = (1 - s0Target) * constant[i] + s0Target * coef[i];
    }
  }
  return { coef, constant };
};

export const fqToGpdf = (q, fq, r, { lorch = false, lowQCorrection = false, s0Target = 0 } = {}) => {
  requireIncreasing(q, 'Q');
  let weighted = fq;
  if (lorch) {
    const window = lorchWindow(q, q[q.length - 1]);
    weighted = new Float64Array(fq.length);
    for (let i = 0; i < fq.length; i += 1) weighted[i] = fq[i] * window[i];
  }
  const gpdf = sineTransform(q, weighted, r);
  const scale = 2 / Math.PI;
  for (let i = 0; i < gpdf.length; i += 1) gpdf[i] *= scale;
  if (lowQCorrection && q[0] !== 0) {
    const s0 = fq[0] / q[0] + 1;
    const { coef, constant } = lowQCorrectionBasis(q, r, { lorch, s0Target });
    for (let i = 0; i < gpdf.length; i += 1) gpdf[i] += coef[i] * s0 - constant[i];
  }
  return gpdf;
};

export const fourierFilter = (q, sq, r, { rho0, cutoff, lorch = false, lowQCorrection = false, s0Target = 0 }) => {
  requireIncreasing(q, 'Q');
  requireIncreasing(r, 'r');
  if (q[0] <= 0) throw new Error('fourierFilter requires a strictly positive Q grid');
  const n = q.length;
  const fq = new Float64Array(n);
  for (let i = 0; i < n; i += 1) fq[i] = q[i] * (sq[i] - 1);
  const gpdf = fqToGpdf(q, fq, r, { lorch, lowQCorrection, s0Target });
  const g = new Float64Array(r.length);
  for (let i = 0; i < r.length; i += 1) g[i] = gpdf[i] / (4 * Math.PI * rho0 * r[i]) + 1;

  let sectionEnd = 0;
  while (sectionEnd < r.length && r[sectionEnd] <= cutoff) sectionEnd += 1;
  const rSection = r.subarray ? r.subarray(0, sectionEnd) : r.slice(0, sectionEnd);
  const ySection = new Float64Array(sectionEnd);
  for (let i = 0; i < sectionEnd; i += 1) ySection[i] = 4 * Math.PI * rho0 * r[i] * g[i];
  const fqFt = sineTransform(rSection, ySection, q);

  const sqFt = new Float64Array(n);
  const sqFiltered = new Float64Array(n);
  const fqFiltered = new Float64Array(n);
  for (let i = 0; i < n; i += 1) {
    sqFt[i] = fqFt[i] / q[i] + 1;
    fqFiltered[i] = fq[i] - fqFt[i];
    sqFiltered[i] = fqFiltered[i] / q[i] + 1;
  }
  const gpdfFiltered = fqToGpdf(q, fqFiltered, r, { lorch, lowQCorrection, s0Target });
  const gFiltered = new Float64Array(r.length);
  for (let i = 0; i < r.length; i += 1) {
    gFiltered[i] = gpdfFiltered[i] / (4 * Math.PI * rho0 * r[i]) + 1;
  }
  return { sqFiltered, sqFt, gFiltered };
};

export const firstPeakZero = (r, g, { cutoff, peakRmin, peakRmax }) => {
  const out = Float64Array.from(g);
  for (let i = 0; i < r.length; i += 1) {
    if (r[i] <= cutoff && (r[i] >= peakRmax || r[i] <= peakRmin)) out[i] = 0;
  }
  return out;
};

// ---------------------------------------------------------------------------
// Level sweep (port of scaling.level_sweep)
// ---------------------------------------------------------------------------

export const levelSweep = (qIn, sqIn, { minWidth = 3.0, nGrid = 80, slopeNsigma = 2.0 } = {}) => {
  const q = [];
  const sq = [];
  for (let i = 0; i < qIn.length; i += 1) {
    if (isNum(qIn[i]) && isNum(sqIn[i])) {
      q.push(qIn[i]);
      sq.push(sqIn[i]);
    }
  }
  const n = q.length;
  if (n < 32) throw new Error('levelSweep needs at least 32 finite points');

  const ones = new Float64Array(n + 1);
  const sumQ = new Float64Array(n + 1);
  const sumQ2 = new Float64Array(n + 1);
  const sumY = new Float64Array(n + 1);
  const sumQY = new Float64Array(n + 1);
  const sumY2 = new Float64Array(n + 1);
  for (let i = 0; i < n; i += 1) {
    ones[i + 1] = ones[i] + 1;
    sumQ[i + 1] = sumQ[i] + q[i];
    sumQ2[i + 1] = sumQ2[i] + q[i] * q[i];
    sumY[i + 1] = sumY[i] + sq[i];
    sumQY[i + 1] = sumQY[i] + q[i] * sq[i];
    sumY2[i + 1] = sumY2[i] + sq[i] * sq[i];
  }

  // np.linspace(0, n-1, nGrid).astype(int) parity: truncation, not rounding.
  const edgeSet = new Set();
  for (let k = 0; k < nGrid; k += 1) {
    edgeSet.add(Math.trunc((k * (n - 1)) / (nGrid - 1)));
  }
  const edges = [...edgeSet].sort((a, b) => a - b);

  let best = null;
  const admissible = [];
  for (let ii = 0; ii < edges.length - 1; ii += 1) {
    const i = edges[ii];
    for (let jj = ii + 1; jj < edges.length; jj += 1) {
      const j = edges[jj];
      if (q[j] - q[i] < minWidth) continue;
      const count = ones[j + 1] - ones[i];
      const sQ = sumQ[j + 1] - sumQ[i];
      const sQ2 = sumQ2[j + 1] - sumQ2[i];
      const sY = sumY[j + 1] - sumY[i];
      const sQY = sumQY[j + 1] - sumQY[i];
      const sY2 = sumY2[j + 1] - sumY2[i];
      const det = count * sQ2 - sQ * sQ;
      if (det <= 0 || count < 24) continue;
      const slope = (count * sQY - sQ * sY) / det;
      const intercept = (sQ2 * sY - sQ * sQY) / det;
      const rss = sY2 - intercept * sY - slope * sQY;
      const sigma2 = Math.max(rss / (count - 2), 0);
      const varSlope = (sigma2 * count) / det;
      const level = intercept + slope * (sQ / count);
      const varLevel = sigma2 / count;
      const slopeSigma = Math.sqrt(Math.max(varSlope, 1e-30));
      if (Math.abs(slope) < slopeNsigma * slopeSigma) {
        admissible.push(level);
        if (best === null || varLevel < best.score) {
          best = { score: varLevel, level, qLo: q[i], qHi: q[j], slope, slopeSigma };
        }
      }
    }
  }

  if (best === null) {
    const qMax = q[n - 1];
    const tail = [];
    for (let i = 0; i < n; i += 1) if (q[i] >= qMax - minWidth) tail.push(sq[i]);
    return {
      level: median(tail),
      levelUncertainty: NaN,
      qLo: qMax - minWidth,
      qHi: qMax,
      slope: NaN,
      slopeSigma: NaN,
      asymptoteFound: false,
      nAdmissible: 0,
    };
  }
  const levelMean = mean(admissible);
  let variance = 0;
  for (let i = 0; i < admissible.length; i += 1) {
    variance += (admissible[i] - levelMean) ** 2;
  }
  return {
    level: best.level,
    levelUncertainty: Math.sqrt(variance / admissible.length),
    qLo: best.qLo,
    qHi: best.qHi,
    slope: best.slope,
    slopeSigma: best.slopeSigma,
    asymptoteFound: true,
    nAdmissible: admissible.length,
  };
};

// ---------------------------------------------------------------------------
// The auto-scaler (port of scaling.autoscale / scale_pipeline)
// ---------------------------------------------------------------------------

const HUBER_C = 1.345;

const huberWeights = (residuals) => {
  const med = median(residuals);
  const deviations = residuals.map((value) => Math.abs(value - med));
  const scale = 1.4826 * median(deviations);
  const weights = new Float64Array(residuals.length);
  if (scale <= 1e-14) {
    weights.fill(1);
    return weights;
  }
  for (let i = 0; i < residuals.length; i += 1) {
    const z = Math.max(Math.abs(residuals[i]) / scale, 1e-14);
    weights[i] = Math.min(1, HUBER_C / z);
  }
  return weights;
};

/** Gap between the first-shell onset and the automatic window top (scaling.R0_WINDOW_MARGIN). */
export const R0_WINDOW_MARGIN = 0.25;
/** Narrowest low-r fit window autoscale places on its own (scaling.MIN_AUTO_WINDOW). */
export const MIN_AUTO_WINDOW = 0.1;
/** Trial window widths above lo used to locate the first shell (scaling.START_WINDOW_WIDTHS). */
export const START_WINDOW_WIDTHS = [0.3, 1.0];
/** Onsets closer than this (Å) are the same first shell (scaling.ONSET_TOLERANCE). */
export const ONSET_TOLERANCE = 0.15;
/** Most refits autoscale spends confirming first-shell candidates (scaling.MAX_WINDOW_REFITS). */
export const MAX_WINDOW_REFITS = 4;
/** Number densities (atoms/Å³) estimateRho0 may return (scaling.RHO0_PHYSICAL_RANGE). */
export const RHO0_PHYSICAL_RANGE = [0.005, 0.25];
/** Width (Å⁻¹) of the low-Q head the Faber-Ziman extrapolation is fitted on (scaling.FZ_FIT_WIDTH). */
export const FZ_FIT_WIDTH = 1.0;
/** Rounding slack of the <b^2> >= <b>^2 check (scaling.B_SQ_RTOL). */
export const B_SQ_RTOL = 1e-9;
/** <b>^2 values closer than this (relative) are the same scattering-length set (scaling_cli.COEFFICIENT_RTOL). */
export const COEFFICIENT_RTOL = 0.02;

/**
 * <b>^2 and <b^2> from ONE consistent source (port of
 * scaling_cli.resolve_coefficients — same rule): explicit values win; a
 * formula fills what is missing, but its <b^2> is paired with a <b>^2 from
 * elsewhere (a stog.inp or an override) only when the two <b>^2 agree within
 * COEFFICIENT_RTOL — then the formula's ratio <b^2>/<b>^2 is kept on the
 * configured scale. Otherwise <b^2> stays unset (`dropped` carries the
 * formula's pair) rather than fabricating an S(0) target from two sources.
 */
export const resolveCoefficients = ({ bAvgSq, bSqAvg, formula }) => {
  let outBAvgSq = bAvgSq ?? undefined;
  let outBSqAvg = bSqAvg ?? undefined;
  let fromFormula = { bAvgSq: false, bSqAvg: false };
  let dropped = null;
  let disagree = false;
  const name = (formula || '').trim();
  let fz = null;
  if (name) {
    fz = faberZiman(name);
    if (outBAvgSq === undefined) {
      outBAvgSq = fz.bAvgSqBarn;
      fromFormula = { ...fromFormula, bAvgSq: true };
    }
    disagree = Math.abs(fz.bAvgSqBarn - outBAvgSq) > COEFFICIENT_RTOL * Math.abs(outBAvgSq);
    if (outBSqAvg === undefined) {
      if (!disagree) {
        outBSqAvg = outBAvgSq * (fz.bSqAvgBarn / fz.bAvgSqBarn);
        fromFormula = { ...fromFormula, bSqAvg: true };
      } else {
        dropped = { bAvgSq: fz.bAvgSqBarn, bSqAvg: fz.bSqAvgBarn };
      }
    }
  }
  return {
    bAvgSq: outBAvgSq, bSqAvg: outBSqAvg, fz, fromFormula, dropped, disagree,
  };
};

export const defaultConfig = {
  qmin: NaN,
  qmax: NaN,
  rho0: NaN,
  bAvgSq: NaN,
  bSqAvg: null,
  rCutoff: 1.0,
  r0: null,
  fitOffset: true,
  qTailFrac: 0.15,
  rFitMin: null,
  rFitMax: null,
  rmax: 50.0,
  nr: 5000,
  lorch: false,
  lowQCorrection: true,
  c2Weight: 1.0,
  robust: true,
  c1Mode: 'sweep',
  amplitudeCriterion: 'density',
  s0Target: null,
  despike: false,
  despikeWindow: 7,
  despikeNsigma: 6.0,
  maxIter: 50,
  tol: 1e-6,
};

export const makeConfig = (options) => {
  const config = { ...defaultConfig, ...options };
  if (!isNum(config.rho0) || config.rho0 <= 0) throw new Error(`rho0 must be finite and positive, got ${config.rho0}`);
  if (!isNum(config.bAvgSq) || config.bAvgSq <= 0) throw new Error(`bAvgSq must be finite and positive, got ${config.bAvgSq}`);
  if (config.bSqAvg != null) {
    if (!isNum(config.bSqAvg) || config.bSqAvg <= 0) throw new Error(`bSqAvg must be finite and positive, got ${config.bSqAvg}`);
    if (config.bSqAvg < config.bAvgSq * (1 - B_SQ_RTOL)) {
      // scaling.ScalingConfig parity: S(0) = 1 - <b^2>/<b>^2 > 0 is impossible.
      throw new Error(
        `<b^2> = ${fmtP(config.bSqAvg, 6)} barn is smaller than <b>^2 = `
        + `${fmtP(config.bAvgSq, 6)} barn, so S(0) = 1 - <b^2>/<b>^2 = `
        + `${fmtP(1 - config.bSqAvg / config.bAvgSq, 4)} > 0, which is impossible `
        + '(<b^2> >= <b>^2, Cauchy-Schwarz): the two coefficients must come '
        + 'from the same source (same radiation and units) — e.g. normalized '
        + 'x-ray data need <b>^2 = 1 with <b^2> = <Z^2>/<Z>^2, not a neutron '
        + "composition's <b^2>"
      );
    }
  }
  if (!(config.qmax > config.qmin)) throw new Error('qmax must exceed qmin');
  if (!Number.isInteger(config.nr) || config.nr <= 0) throw new Error(`nr must be a positive integer, got ${config.nr}`);
  if (!isNum(config.rmax) || config.rmax <= 0) throw new Error(`rmax must be finite and positive, got ${config.rmax}`);
  if (config.c1Mode !== 'sweep' && config.c1Mode !== 'joint') throw new Error(`c1Mode must be 'sweep' or 'joint', got ${config.c1Mode}`);
  if (config.amplitudeCriterion !== 'density' && config.amplitudeCriterion !== 'fz') {
    throw new Error(`amplitudeCriterion must be 'density' or 'fz', got ${config.amplitudeCriterion}`);
  }
  if (config.amplitudeCriterion === 'fz') {
    if (config.bSqAvg == null) throw new Error("amplitudeCriterion='fz' requires bSqAvg (<b^2>)");
    if (config.c1Mode !== 'sweep') throw new Error("amplitudeCriterion='fz' requires c1Mode='sweep'");
  }
  rFitWindow(config); // validate eagerly
  return config;
};

/** S(0) used by the omitted-low-Q extrapolation (composition-aware default). */
export const effectiveS0Target = (config) => {
  if (config.s0Target != null) return config.s0Target;
  if (config.bSqAvg != null) return 1 - config.bSqAvg / config.bAvgSq;
  return 0;
};

export const rGrid = (config) => {
  const step = config.rmax / config.nr;
  const r = new Float64Array(config.nr);
  for (let i = 0; i < config.nr; i += 1) r[i] = (i + 1) * step;
  return r;
};

export const rFitWindow = (config) => {
  const lo = config.rFitMin != null ? config.rFitMin : config.rCutoff + 0.2;
  let hi;
  if (config.rFitMax != null) hi = config.rFitMax;
  else if (config.r0 != null) hi = config.r0 - R0_WINDOW_MARGIN;
  else hi = lo + 1.0;
  if (hi <= lo) throw new Error(`empty low-r fit window [${lo}, ${hi}]`);
  return [lo, hi];
};

const qTailMin = (config) => config.qmax - config.qTailFrac * (config.qmax - config.qmin);

const despikeKeepMask = (sq, window, nsigma) => {
  const n = sq.length;
  const pad = Math.floor(window / 2);
  const rolled = new Float64Array(n);
  const buffer = new Float64Array(window);
  for (let i = 0; i < n; i += 1) {
    for (let k = 0; k < window; k += 1) {
      let index = i - pad + k;
      if (index < 0) index = 0;
      if (index >= n) index = n - 1;
      buffer[k] = sq[index];
    }
    rolled[i] = median(buffer);
  }
  const residual = new Float64Array(n);
  for (let i = 0; i < n; i += 1) residual[i] = sq[i] - rolled[i];
  const mad = 1.4826 * median(residual.map(Math.abs));
  const bound = nsigma * Math.max(mad, 1e-12);
  const keep = new Uint8Array(n);
  for (let i = 0; i < n; i += 1) keep[i] = Math.abs(residual[i]) <= bound ? 1 : 0;
  return keep;
};

/**
 * The σ column only when it is clean (scaling_cli.usable_sigma — the CLI/API
 * guard): if any row with finite Q and S has a non-finite or non-positive σ,
 * the whole column is dropped (a zero σ would get a 1e12 weight and one NaN σ
 * turns every weight NaN). Returns { sigma: sigma | null, nBad }.
 */
export const usableSigma = (q, sq, sigma) => {
  if (!sigma) return { sigma: null, nBad: 0 };
  let nBad = 0;
  for (let i = 0; i < q.length; i += 1) {
    if (isNum(q[i]) && isNum(sq[i]) && !(isNum(sigma[i]) && sigma[i] > 0)) nBad += 1;
  }
  return { sigma: nBad ? null : sigma, nBad };
};

const fmtQ = (value) => String(Number(value.toPrecision(6)));

/**
 * Index order that sorts the cropped q ascending — or throw (port of
 * scaling._ascending_order, same rules and messages): a descending file is read
 * backwards and non-overlapping segments in any order are sorted; duplicate Q or
 * segments whose Q ranges overlap (concatenated banks) throw.
 */
const ascendingOrder = (q) => {
  const n = q.length;
  const index = Array.from({ length: n }, (_, i) => i);
  let ascending = true;
  let up = 0;
  let down = 0;
  for (let i = 1; i < n; i += 1) {
    if (q[i] > q[i - 1]) up += 1;
    else {
      ascending = false;
      if (q[i] < q[i - 1]) down += 1;
    }
  }
  if (ascending) return index;
  const order = [...index].sort((i, j) => q[i] - q[j] || i - j);
  for (let k = 1; k < n; k += 1) {
    if (q[order[k]] === q[order[k - 1]]) {
      throw new Error(
        `S(Q) has duplicate Q values (e.g. Q = ${fmtQ(q[order[k]])} Å⁻¹ appears more `
        + 'than once) inside [qmin, qmax]: merge or rebin the data (several '
        + 'detector banks?) into one S(Q) before scaling'
      );
    }
  }
  const read = down > up ? [...index].reverse() : index;
  const spans = [];
  let start = 0;
  for (let k = 1; k <= n; k += 1) {
    if (k === n || q[read[k]] < q[read[k - 1]]) {
      spans.push([q[read[start]], q[read[k - 1]]]);
      start = k;
    }
  }
  spans.sort((left, right) => left[0] - right[0]);
  for (let k = 1; k < spans.length; k += 1) {
    const [lo1, hi1] = spans[k - 1];
    const [lo2, hi2] = spans[k];
    if (lo2 < hi1) {
      throw new Error(
        'S(Q) rows are not in one monotonic Q order: segments overlap in Q '
        + `([${fmtQ(lo1)}, ${fmtQ(hi1)}] and [${fmtQ(lo2)}, ${fmtQ(hi2)}] Å⁻¹) — typically `
        + 'several detector banks concatenated, or rows out of order; merge '
        + 'them into one monotonic S(Q) before scaling'
      );
    }
  }
  return order;
};

export const cropSq = (qIn, sqIn, config, sigmaIn = null) => {
  const qKept = [];
  const sqKept = [];
  const sigmaKept = sigmaIn ? [] : null;
  for (let i = 0; i < qIn.length; i += 1) {
    const qi = qIn[i];
    const si = sqIn[i];
    if (!isNum(qi) || !isNum(si)) continue;
    if (qi <= 0 || qi < config.qmin - 1e-12 || qi > config.qmax + 1e-12) continue;
    qKept.push(qi);
    sqKept.push(si);
    if (sigmaKept) sigmaKept.push(sigmaIn[i]);
  }
  if (qKept.length < 16) throw new Error('fewer than 16 usable S(Q) points after cropping');
  // Ascending Q (scaling.crop_sq): the transforms need an increasing grid, and
  // the despike's rolling median runs on the sorted data.
  const order = ascendingOrder(qKept);
  const q = order.map((i) => qKept[i]);
  const sq = order.map((i) => sqKept[i]);
  const sigma = sigmaKept ? order.map((i) => sigmaKept[i]) : null;
  if (!config.despike) {
    return { q: Float64Array.from(q), sq: Float64Array.from(sq), sigma: sigma ? Float64Array.from(sigma) : null, nDespiked: 0 };
  }
  const keep = despikeKeepMask(sq, config.despikeWindow, config.despikeNsigma);
  const q2 = [];
  const sq2 = [];
  const sigma2 = sigma ? [] : null;
  for (let i = 0; i < q.length; i += 1) {
    if (!keep[i]) continue;
    q2.push(q[i]);
    sq2.push(sq[i]);
    if (sigma2) sigma2.push(sigma[i]);
  }
  return {
    q: Float64Array.from(q2),
    sq: Float64Array.from(sq2),
    sigma: sigma2 ? Float64Array.from(sigma2) : null,
    nDespiked: q.length - q2.length,
  };
};

const fitWindows = (q, r, config) => {
  const tailMin = qTailMin(config);
  const tail = [];
  for (let i = 0; i < q.length; i += 1) if (q[i] >= tailMin) tail.push(i);
  const [lo, hi] = rFitWindow(config);
  const window = [];
  for (let i = 0; i < r.length; i += 1) if (r[i] >= lo && r[i] <= hi) window.push(i);
  if (tail.length < 2 || window.length < 2) throw new Error('fit windows contain fewer than 2 points');
  return { tail, window };
};

const solveAffine = (q, sq, deltaSq, r, tailIdx, windowIdx, config, sigma, level) => {
  const rw = Float64Array.from(windowIdx, (index) => r[index]);
  const n = q.length;
  const fqData = new Float64Array(n);
  const fqOne = new Float64Array(n);
  const fqDelta = new Float64Array(n);
  for (let i = 0; i < n; i += 1) {
    fqData[i] = q[i] * (sq[i] - 1);
    fqOne[i] = q[i];
    fqDelta[i] = q[i] * deltaSq[i];
  }
  const gData = fqToGpdf(q, fqData, rw, { lorch: config.lorch });
  const gOne = fqToGpdf(q, fqOne, rw, { lorch: config.lorch });
  const gDelta = fqToGpdf(q, fqDelta, rw, { lorch: config.lorch });
  let coef = new Float64Array(rw.length);
  let constant = new Float64Array(rw.length);
  if (config.lowQCorrection) {
    ({ coef, constant } = lowQCorrectionBasis(q, rw, {
      lorch: config.lorch,
      s0Target: effectiveS0Target(config),
    }));
  }
  const w2 = Math.sqrt(config.c2Weight);

  const n1 = tailIdx.length;
  const n2 = rw.length;
  const colA = new Float64Array(n1 + n2);
  const colB = new Float64Array(n1 + n2);
  const rhs = new Float64Array(n1 + n2);
  for (let row = 0; row < n1; row += 1) {
    const i = tailIdx[row];
    let weight = 1;
    if (sigma) {
      weight = 1 / Math.max(sigma[i], 1e-12);
    }
    colA[row] = sq[i] * weight;
    colB[row] = weight;
    rhs[row] = (1 + deltaSq[i]) * weight;
  }
  if (sigma) {
    // Normalize the sigma weights to unit mean over the tail block (python parity).
    let total = 0;
    for (let row = 0; row < n1; row += 1) total += colB[row];
    const norm = total / n1;
    for (let row = 0; row < n1; row += 1) {
      colA[row] /= norm;
      colB[row] /= norm;
      rhs[row] /= norm;
    }
  }
  const delta0 = deltaSq[0];
  for (let row = 0; row < n2; row += 1) {
    const denom = 4 * Math.PI * config.rho0 * rw[row];
    colA[n1 + row] = (w2 * (gData[row] + gOne[row] + coef[row] * sq[0])) / denom;
    colB[n1 + row] = (w2 * (gOne[row] + coef[row])) / denom;
    rhs[n1 + row] = w2 * ((gDelta[row] + gOne[row] + coef[row] * delta0 + constant[row]) / denom - 1);
  }

  let columns;
  let effRhs = rhs;
  if (level != null) {
    const anchored = new Float64Array(n1 + n2);
    effRhs = new Float64Array(n1 + n2);
    for (let row = 0; row < n1 + n2; row += 1) {
      anchored[row] = colA[row] - level * colB[row];
      effRhs[row] = rhs[row] - colB[row];
    }
    columns = [anchored];
  } else if (config.fitOffset) {
    columns = [colA, colB];
  } else {
    columns = [colA];
  }

  let solution = solveLeastSquares(columns, effRhs);
  if (config.robust) {
    for (let pass = 0; pass < 3; pass += 1) {
      const residuals = new Float64Array(effRhs.length);
      for (let row = 0; row < effRhs.length; row += 1) {
        let predicted = 0;
        for (let c = 0; c < columns.length; c += 1) predicted += columns[c][row] * solution[c];
        residuals[row] = predicted - effRhs[row];
      }
      const weights = new Float64Array(effRhs.length).fill(1);
      const w1 = huberWeights(Array.from(residuals.slice(0, n1)));
      for (let row = 0; row < n1; row += 1) weights[row] = w1[row];
      if (effRhs.length - n1 >= 4) {
        const wc2 = huberWeights(Array.from(residuals.slice(n1)));
        for (let row = 0; row < wc2.length; row += 1) weights[n1 + row] = wc2[row];
      }
      const weightedCols = columns.map((column) => {
        const out = new Float64Array(column.length);
        for (let row = 0; row < column.length; row += 1) out[row] = column[row] * weights[row];
        return out;
      });
      const weightedRhs = new Float64Array(effRhs.length);
      for (let row = 0; row < effRhs.length; row += 1) weightedRhs[row] = effRhs[row] * weights[row];
      solution = solveLeastSquares(weightedCols, weightedRhs);
    }
  }
  const a = solution[0];
  let b;
  if (level != null) b = 1 - a * level;
  else if (config.fitOffset) b = solution[1];
  else b = 0;
  return { a, b };
};

const lowRRmsOf = (r, gFiltered, config) => {
  const [lo, hi] = rFitWindow(config);
  let total = 0;
  let count = 0;
  for (let i = 0; i < r.length; i += 1) {
    if (r[i] >= lo && r[i] <= hi) {
      total += gFiltered[i] * gFiltered[i];
      count += 1;
    }
  }
  return Math.sqrt(total / count);
};

export const amplitudeFromFzLimit = (q, sq, level, config, { fitWidth = FZ_FIT_WIDTH } = {}) => {
  if (config.bSqAvg == null) return null;
  const s0Target = 1 - config.bSqAvg / config.bAvgSq;
  const headIdx = [];
  for (let i = 0; i < q.length; i += 1) if (q[i] <= q[0] + fitWidth) headIdx.push(i);
  if (headIdx.length < 8) return null;
  const qHead = headIdx.map((i) => q[i]);
  const yHead = headIdx.map((i) => sq[i]);
  const qMean = mean(qHead);
  const qc = qHead.map((value) => value - qMean);
  const onesCol = new Float64Array(qc.length).fill(1);
  const qcCol = Float64Array.from(qc);
  let weights = new Float64Array(qc.length).fill(1);
  let solution = [median(yHead), 0];
  for (let pass = 0; pass < 4; pass += 1) {
    const wCols = [
      Float64Array.from(onesCol, (value, i) => value * weights[i]),
      Float64Array.from(qcCol, (value, i) => value * weights[i]),
    ];
    const wRhs = Float64Array.from(yHead, (value, i) => value * weights[i]);
    solution = solveLeastSquares(wCols, wRhs);
    const residuals = yHead.map((value, i) => solution[0] + solution[1] * qc[i] - value);
    weights = huberWeights(residuals);
  }
  const sMeas0 = solution[0] - solution[1] * qMean;
  const denom = sMeas0 - level;
  if (Math.abs(denom) < 1e-9) return null;
  return (s0Target - 1) / denom;
};

export const scalePipeline = (qIn, sqIn, config, a, b, extras = {}) => {
  const { q, sq, nDespiked } = cropSq(qIn, sqIn, config);
  const r = rGrid(config);
  const sqScaled = new Float64Array(q.length);
  for (let i = 0; i < q.length; i += 1) sqScaled[i] = a * sq[i] + b;
  const { sqFiltered, sqFt, gFiltered } = fourierFilter(q, sqScaled, r, {
    rho0: config.rho0,
    cutoff: config.rCutoff,
    lorch: config.lorch,
    lowQCorrection: config.lowQCorrection,
    s0Target: effectiveS0Target(config),
  });
  const gk = new Float64Array(r.length);
  const dr = new Float64Array(r.length);
  for (let i = 0; i < r.length; i += 1) {
    gk[i] = config.bAvgSq * (gFiltered[i] - 1);
    dr[i] = 4 * Math.PI * config.rho0 * r[i] * gk[i];
  }
  const fk = new Float64Array(q.length);
  for (let i = 0; i < q.length; i += 1) fk[i] = config.bAvgSq * (sqFiltered[i] - 1);

  const { tail } = fitWindows(q, r, config);
  let tailTotal = 0;
  for (let i = 0; i < tail.length; i += 1) tailTotal += sqFiltered[tail[i]];

  return {
    a,
    b,
    converged: extras.converged !== undefined ? extras.converged : true,
    iterations: extras.iterations || 0,
    history: extras.history || [],
    q,
    sqRaw: sq,
    sqScaled,
    sqFiltered,
    sqFt,
    r,
    gFiltered,
    gk,
    dr,
    fk,
    lowRRms: lowRRmsOf(r, gFiltered, config),
    c1TailMean: tailTotal / tail.length,
    sweep: extras.sweep || null,
    aFz: extras.aFz !== undefined ? extras.aFz : null,
    nDespiked,
    mode: extras.mode || 'manual',
    c1ModeEffective: extras.c1ModeEffective || 'manual',
    rFitWindowUsed: rFitWindow(config),
    r0Detected: extras.r0Detected !== undefined ? extras.r0Detected : null,
    windowRefined: Boolean(extras.windowRefined),
    fitFailure: extras.fitFailure || null,
  };
};

/**
 * Onsets of the shell-like features of g(r), increasing (port of
 * scaling.first_shell_candidates — keep in sync). Local maxima of |g| beyond a
 * one-ripple-period reference zone (2π/qmax above searchMin) that stand out of
 * the ripple field below them (max |g| from searchMin to the lobe start) by
 * strongProminence, or by prominence while being a major feature (>= major of
 * the range maximum); either sign, so inverted shells count. Each onset is on
 * its shell's own flank (|g| down to max(floor, fraction·peak)); a shell whose
 * flank reaches below searchMin is not separable and ends the scan.
 */
export const firstShellCandidates = (
  r, g, {
    searchMin = 1.0, searchMax = 6.0, fraction = 0.35, floor = 0.5,
    prominence = 2.0, major = 0.5, strongProminence = 4.0, qmax = 0,
  } = {}
) => {
  let first = -1;
  let last = -1;
  for (let i = 0; i < r.length; i += 1) {
    if (r[i] >= searchMin && r[i] <= searchMax) {
      if (first < 0) first = i;
      last = i;
    }
  }
  if (first < 0 || last - first + 1 < 3) return [];
  const a = new Float64Array(r.length);
  for (let i = 0; i < r.length; i += 1) a[i] = Math.abs(g[i]);
  let globalMax = -Infinity;
  for (let i = first; i <= last; i += 1) if (a[i] > globalMax) globalMax = a[i];
  if (globalMax < floor) return [];
  const zoneEnd = searchMin + (qmax > 0 ? (2 * Math.PI) / qmax : 0);
  const onsets = [];
  for (let index = first + 1; index < last; index += 1) {
    const peak = a[index];
    if (r[index] < zoneEnd || peak < floor) continue;
    if (!(peak >= a[index - 1] && peak > a[index + 1])) continue;
    let start = index;
    while (start > first && a[start - 1] < a[start]) start -= 1;
    let ripple = -Infinity;
    for (let i = first; i <= start; i += 1) if (a[i] > ripple) ripple = a[i];
    if (!(peak >= strongProminence * ripple
      || (peak >= prominence * ripple && peak >= major * globalMax))) continue;
    const level = Math.max(floor, fraction * peak);
    let onset = index;
    while (onset > first && a[onset] > level) onset -= 1;
    if (a[onset] > level) break; // the flank reaches below the search range: not separable
    const value = r[onset + 1];
    if (!onsets.length || value > onsets[onsets.length - 1]) onsets.push(value);
  }
  return onsets;
};

/**
 * Data-derived closest approach: the rising flank of the FIRST shell (port of
 * scaling.detect_first_peak_onset — keep in sync): the first entry of
 * firstShellCandidates, or null.
 */
export const detectFirstPeakOnset = (r, g, options = {}) => {
  const onsets = firstShellCandidates(r, g, options);
  return onsets.length ? onsets[0] : null;
};

/**
 * Foot of the first shell (port of scaling.first_shell_foot): walk left from
 * the onset while |g| keeps decreasing and g keeps the shell's sign; stop at
 * the first local minimum of |g|, or return the point just across a sign
 * change.
 */
export const firstShellFoot = (r, g, onset) => {
  let index = 0;
  let best = Infinity;
  for (let i = 0; i < r.length; i += 1) {
    const distance = Math.abs(r[i] - onset);
    if (distance < best) {
      best = distance;
      index = i;
    }
  }
  const sign = Math.sign(g[index]);
  while (index > 0) {
    if (Math.sign(g[index - 1]) !== sign) return r[index - 1];
    if (Math.abs(g[index - 1]) >= Math.abs(g[index])) break;
    index -= 1;
  }
  return r[index];
};

/**
 * Automatic low-r enforcement cutoff (port of scaling.auto_enforcement_cutoff):
 * anchor = the first-shell onset capped at a pinned config.r0 (never
 * overridden upward); cutoff = min(foot of the first shell below the anchor,
 * anchor - R0_WINDOW_MARGIN) — below the rising flank, never on it. onset
 * defaults to detectFirstPeakOnset on g; null when there is neither a
 * detected shell nor a pinned r0.
 */
export const autoEnforcementCutoff = (r, g, config, onset = null) => {
  let anchor = onset;
  if (anchor == null) {
    anchor = detectFirstPeakOnset(r, g, { searchMin: config.rCutoff + 0.3, qmax: config.qmax });
  }
  if (config.r0 != null) anchor = anchor == null ? config.r0 : Math.min(anchor, config.r0);
  if (anchor == null) return null;
  return Math.min(firstShellFoot(r, g, anchor), anchor - R0_WINDOW_MARGIN);
};

const detectOnset = (result, config) => detectFirstPeakOnset(result.r, result.gFiltered, {
  searchMin: config.rCutoff + 0.3,
  qmax: config.qmax,
});

const shellCandidates = (result, config) => firstShellCandidates(result.r, result.gFiltered, {
  searchMin: config.rCutoff + 0.3,
  qmax: config.qmax,
});

const near = (value, onsets) => onsets.some((other) => Math.abs(value - other) <= ONSET_TOLERANCE);

const fmtG = (value) => String(Number(value.toPrecision(6)));
const fmtP = (value, digits) => String(Number(value.toPrecision(digits)));

/** How to make room for a low-r window below a first shell at onset (scaling._cutoff_advice). */
const cutoffAdvice = (onset, config) => (config.rFitMin != null
  ? `lower the fit-window minimum below ${(onset - R0_WINDOW_MARGIN - MIN_AUTO_WINDOW).toFixed(2)} Å`
  : `lower the filter r-cut to <= ${(Math.floor((onset - R0_WINDOW_MARGIN - MIN_AUTO_WINDOW - 0.2) / 0.05) * 0.05).toFixed(2)} Å`);

/**
 * Auto-scale with data-located low-r window (port of scaling.autoscale — keep
 * in sync). With r0 / rFitMax pinned, or in fz mode, one pass runs and the
 * detection only annotates it (fz: the diagnostic window is still refined when
 * the shell leaves room). Otherwise two trial windows [lo, lo + w] propose
 * shell onsets (whatever the sign of their scale); the smallest untried onset
 * c is refitted on [lo, c - 0.25]: a <= 0 throws, the refit re-detecting c
 * (within ONSET_TOLERANCE) confirms it, a lower re-detected shell is tried
 * next, and otherwise c is dropped as a feature of the trial's scale. A
 * confirmed shell too close to lo, no confirmable onset, or an exhausted refit
 * budget throw — a fit across the first shell, or with a <= 0, is never
 * returned. r0Detected is the confirmed onset the window was built from.
 */
export const autoscale = (qIn, sqIn, config, sigmaIn = null) => {
  const lo = rFitWindow(config)[0];
  if (config.r0 != null || config.rFitMax != null) {
    const result = autoscalePass(qIn, sqIn, config, sigmaIn);
    const onset = detectOnset(result, config);
    if (onset != null) result.r0Detected = onset;
    return result;
  }

  if (config.amplitudeCriterion === 'fz') {
    const result = autoscalePass(qIn, sqIn, config, sigmaIn);
    const onset = detectOnset(result, config);
    if (onset == null) return result;
    if (onset - R0_WINDOW_MARGIN - lo < MIN_AUTO_WINDOW) {
      result.r0Detected = onset;
      return result;
    }
    const refined = autoscalePass(qIn, sqIn, { ...config, r0: onset }, sigmaIn);
    refined.r0Detected = onset;
    refined.windowRefined = true;
    return refined;
  }

  return placeLowRWindow((trialConfig) => autoscalePass(qIn, sqIn, trialConfig, sigmaIn), config);
};

/**
 * Locate the first shell and fit below it (port of scaling._place_low_r_window
 * — keep in sync): the trial / confirm loop of autoscale, with the fit pass
 * injected as runPass(config) so the placement logic can be replayed on
 * scripted passes (the parity fixture).
 */
export const placeLowRWindow = (runPass, config) => {
  const lo = rFitWindow(config)[0];
  const pool = [];
  const propose = (onsets) => {
    onsets.forEach((value) => {
      if (!near(value, pool)) pool.push(value);
    });
  };

  const failures = [];
  const trialScales = [];
  for (const width of START_WINDOW_WIDTHS) {
    let trial;
    try {
      trial = runPass({ ...config, rFitMax: lo + width });
    } catch (error) { // e.g. a trial window with < 2 r points
      failures.push(error);
      continue;
    }
    trialScales.push(trial.a);
    propose(shellCandidates(trial, config));
  }
  if (failures.length === START_WINDOW_WIDTHS.length) {
    throw failures[failures.length - 1]; // the data cannot be fitted at all: report why
  }

  const fixes = 'set r₀ (the closest interatomic approach) or the fit-window maximum';
  const refits = new Map();
  const dropped = []; // onsets their own refit no longer shows
  for (;;) {
    const live = pool.filter((value) => !near(value, dropped));
    if (!live.length) break;
    const onset = Math.min(...live);
    const room = onset - R0_WINDOW_MARGIN - lo;
    if (!refits.has(onset)) {
      if (refits.size >= MAX_WINDOW_REFITS) {
        throw new Error(
          'autoscale: no first-shell onset was confirmed within '
          + `${MAX_WINDOW_REFITS} refits of the low-r window (tried `
          + `${[...refits.keys()].sort((x, y) => x - y).map((value) => value.toFixed(2)).join(', ')} Å). `
          + `${fixes[0].toUpperCase()}${fixes.slice(1)}`
        );
      }
      let refined;
      try {
        refined = runPass({ ...config, r0: onset });
      } catch (error) {
        if (room >= MIN_AUTO_WINDOW) throw error;
        throw new Error(
          `autoscale: a shell-like feature starts at ${onset.toFixed(2)} Å, too close to `
          + `the fit-window start ${fmtG(lo)} Å to place or verify a low-r window below it `
          + `(${error.message}). If it is the first coordination shell (a bond that short), `
          + `${cutoffAdvice(onset, config)}; otherwise ${fixes}`
        );
      }
      if (!(refined.a > 0)) {
        const where = `${room < MIN_AUTO_WINDOW ? 'the narrow window ' : ''}`
          + `[${fmtG(lo)}, ${(onset - R0_WINDOW_MARGIN).toFixed(2)}] Å`;
        throw new Error(
          'autoscale: the density-limit fit below the first-shell candidate at '
          + `${onset.toFixed(2)} Å gives a non-physical scale (a = ${fmtP(refined.a, 4)} on `
          + `${where}): the low-r region cannot be modelled as g = 0 there, so no `
          + 'automatic window is trustworthy (typical of data missing structure below '
          + `Qmin, where the density limit is degenerate). ${fixes[0].toUpperCase()}${fixes.slice(1)}, `
          + 'or use the Faber-Ziman Q→0 amplitude criterion when the composition is known'
        );
      }
      const found = shellCandidates(refined, config);
      refits.set(onset, [refined, found]);
      propose(found);
    }
    const [refined, foundAll] = refits.get(onset);
    const found = foundAll.filter((value) => !near(value, dropped));
    if (found.length && Math.abs(found[0] - onset) <= ONSET_TOLERANCE) {
      if (room < MIN_AUTO_WINDOW) {
        const advice = cutoffAdvice(onset, config);
        throw new Error(
          `autoscale: the first coordination shell starts at ${onset.toFixed(2)} Å, `
          + `leaving no low-r fit window between ${fmtG(lo)} Å and `
          + `${(onset - R0_WINDOW_MARGIN).toFixed(2)} Å (onset - ${R0_WINDOW_MARGIN}); a window `
          + 'across the shell would force it to zero and bias the scale. '
          + `${advice[0].toUpperCase()}${advice.slice(1)}, or set r₀ / the fit-window maximum`
        );
      }
      refined.r0Detected = onset;
      refined.windowRefined = true;
      return refined;
    }
    if (found.length && found[0] < onset - ONSET_TOLERANCE) continue; // a lower shell: tried next
    dropped.push(onset);
  }

  const tried = dropped.length
    ? `; onsets not confirmed by their own refit: ${dropped.map((value) => value.toFixed(2)).join(', ')} Å`
    : '';
  throw new Error(
    'autoscale: could not locate the first coordination shell in the data '
    + `(trial low-r windows [${fmtG(lo)}, ${fmtG(lo + START_WINDOW_WIDTHS[0])}] `
    + `and [${fmtG(lo)}, ${fmtG(lo + START_WINDOW_WIDTHS[START_WINDOW_WIDTHS.length - 1])}] Å, `
    + `scales a = ${trialScales.map((value) => fmtP(value, 3)).join(', ')}${tried}; no shell `
    + 'stands out of the ripples and survives a refit below it), so the '
    + `density-limit window cannot be placed below it. ${fixes[0].toUpperCase()}${fixes.slice(1)}; `
    + `if the first bond is shorter than ~${(lo + 0.55).toFixed(1)} Å, also lower the filter `
    + `r-cut (now ${fmtG(config.rCutoff)} Å)`
  );
};

const autoscalePass = (qIn, sqIn, config, sigmaIn = null) => {
  const { q, sq, sigma } = cropSq(qIn, sqIn, config, sigmaIn);
  const r = rGrid(config);
  const { tail, window } = fitWindows(q, r, config);

  let sweep = null;
  let level = null;
  if (config.c1Mode === 'sweep') {
    sweep = levelSweep(q, sq);
    if (sweep.asymptoteFound) level = sweep.level;
  }

  if (config.amplitudeCriterion === 'fz') {
    if (level == null) {
      throw new Error(
        "amplitudeCriterion='fz': no statistically flat high-Q window found; "
        + 'there is no measured level to anchor'
      );
    }
    const aFz = amplitudeFromFzLimit(q, sq, level, config);
    if (aFz == null || !isNum(aFz) || aFz <= 0) {
      throw new Error(`amplitudeCriterion='fz': degenerate Q->0 extrapolation (aFz=${aFz})`);
    }
    return scalePipeline(qIn, sqIn, config, aFz, 1 - aFz * level, {
      converged: true,
      iterations: 0,
      sweep,
      aFz,
      mode: 'auto',
      c1ModeEffective: 'sweep',
    });
  }

  let deltaSq = new Float64Array(q.length);
  let aPrev = Infinity;
  let bPrev = Infinity;
  let a = 1;
  let b = 0;
  const history = [];
  let converged = false;
  let iterations = 0;

  // The loop variable is kept separate so an unconverged run reports maxIter
  // (Python's `for iterations in range(1, max_iter + 1)`), not maxIter + 1.
  for (let iteration = 1; iteration <= config.maxIter; iteration += 1) {
    iterations = iteration;
    ({ a, b } = solveAffine(q, sq, deltaSq, r, tail, window, config, sigma, level));
    const sqScaled = new Float64Array(q.length);
    for (let i = 0; i < q.length; i += 1) sqScaled[i] = a * sq[i] + b;
    // Same filter as scaling._pipeline, including the composition-aware S(0)
    // target: the fed-back deltaSq must come from the omitted-low-Q model the
    // affine solve (solveAffine) and the final scalePipeline use.
    const { sqFt, gFiltered } = fourierFilter(q, sqScaled, r, {
      rho0: config.rho0,
      cutoff: config.rCutoff,
      lorch: config.lorch,
      lowQCorrection: config.lowQCorrection,
      s0Target: effectiveS0Target(config),
    });
    deltaSq = new Float64Array(q.length);
    for (let i = 0; i < q.length; i += 1) deltaSq[i] = sqFt[i] - 1;
    history.push([a, b, lowRRmsOf(r, gFiltered, config)]);
    if (Math.abs(a - aPrev) <= config.tol * Math.max(1, Math.abs(a))
      && Math.abs(b - bPrev) <= config.tol * Math.max(1, Math.abs(b))) {
      converged = true;
      break;
    }
    aPrev = a;
    bPrev = b;
  }

  // A non-physical scale is a failed fit however stable the iteration was
  // (scaling._autoscale_pass): never converged; the worker refuses to export it.
  let fitFailure = null;
  if (!(isNum(a) && a > 0)) {
    converged = false;
    const [lo, hi] = rFitWindow(config);
    fitFailure = `non-physical scale a = ${fmtP(a, 6)} (a <= 0): the density-limit fit `
      + 'cannot describe these data as a*S + b with a positive scale on the '
      + `low-r window [${fmtP(lo, 4)}, ${fmtP(hi, 4)}] Å`;
  }

  let aFz = null;
  if (level != null) aFz = amplitudeFromFzLimit(q, sq, level, config);

  return scalePipeline(qIn, sqIn, config, a, b, {
    converged,
    fitFailure,
    iterations,
    history,
    sweep,
    aFz,
    mode: 'auto',
    c1ModeEffective: level != null ? 'sweep' : 'joint',
  });
};

/**
 * Self-consistent number density from amplitude-criteria concordance
 * (straight port of scaling.py estimate_rho0 — keep in sync).
 *
 * The density-limit amplitude depends on rho0 (the C2 rows scale with the
 * density line -4*pi*rho0*r); the Q->0 Faber-Ziman amplitude does not (it
 * needs only the composition <b^2> and the measured high-Q level). rho0 is
 * therefore the root of concordance(rho0) = aFz / aDensity(rho0) = 1; since
 * aDensity grows ~linearly with rho0 the fixed-point update
 * rho *= concordance converges in a few autoscale passes. Requires
 * config.bSqAvg; `extrapolated` flags data whose first measured Q (qFirst,
 * after cropping) exceeds the FZ fit width (the estimate is then a starting
 * point, not a measurement). The
 * iterate stays in [rhoMin, rhoMax] (RHO0_PHYSICAL_RANGE) and a concordant
 * root counts only where the density limit holds; every non-converged exit
 * sets `reason`.
 */
export const estimateRho0 = (qIn, sqIn, config, sigmaIn = null, {
  rtol = 1e-3, maxIter = 8, rhoMin = RHO0_PHYSICAL_RANGE[0], rhoMax = RHO0_PHYSICAL_RANGE[1],
} = {}) => {
  if (config.bSqAvg == null) {
    throw new Error(
      'estimateRho0 requires bSqAvg (<b^2>): without a composition the '
      + 'density and the amplitude are degenerate'
    );
  }
  if (!(rhoMin > 0 && rhoMin < rhoMax)) {
    throw new Error(`need 0 < rhoMin < rhoMax, got [${rhoMin}, ${rhoMax}]`);
  }
  let work = { ...config, amplitudeCriterion: 'density', c1Mode: 'sweep' };
  // Where the measured data actually start: the Q->0 extrapolation spans [0, qFirst].
  const qFirst = cropSq(qIn, sqIn, work).q[0];
  let rho = Math.min(Math.max(work.rho0, rhoMin), rhoMax);
  const history = [];
  let converged = false;
  let stopped = null;
  let reason = null;
  let exhausted = true;
  for (let iteration = 0; iteration < maxIter; iteration += 1) {
    work = { ...work, rho0: rho };
    let result;
    try {
      result = autoscale(qIn, sqIn, work, sigmaIn);
    } catch (error) {
      if (!history.length) throw error; // the seed density itself cannot be fitted
      // The fixed-point step left the physical range (the density limit
      // cannot even be placed below the first shell there): stop.
      stopped = `autoscale failed at rho0 = ${fmtG(rho)}: ${error.message}`;
      reason = stopped;
      exhausted = false;
      break;
    }
    if (result.aFz == null || !isNum(result.aFz) || result.aFz <= 0) {
      throw new Error(
        'estimateRho0: no usable Faber-Ziman amplitude (no flat high-Q '
        + 'level, or a degenerate Q->0 extrapolation) — the density cannot '
        + 'be anchored on this data'
      );
    }
    const concordance = result.aFz / result.a;
    history.push([rho, result.a, result.aFz, concordance]);
    if (Math.abs(concordance - 1) <= rtol) {
      // The amplitudes agree; the root is the density only if the fit it
      // rests on satisfies the density limit there (scaling.estimate_rho0).
      const summary = diagnosticsSummary(result, work);
      if (summary.density_limit_satisfied) {
        converged = true;
      } else {
        const [lo, hi] = summary.r_fit_window;
        reason = `the amplitudes agree at rho0 = ${fmtG(rho)} A^-3, but the `
          + 'density-limit criterion fails there (mean g on the low-r '
          + `window [${fmtP(lo, 3)}, ${fmtP(hi, 3)}] A = ${fmtP(summary.g_window_mean, 3)}, `
          + 'not ~0): a spurious root, not the sample density';
      }
      exhausted = false;
      break;
    }
    if (result.a <= 0 || concordance <= 0) {
      // Non-physical density-limit amplitude: the criteria are discordant on
      // this data (typically missing low-Q structure) and no density can
      // reconcile them. Stop instead of iterating garbage.
      reason = `the density-limit amplitude is non-physical at rho0 = ${fmtG(rho)} `
        + `A^-3 (a = ${fmtP(result.a, 4)}, concordance ${fmtP(concordance, 4)}): the `
        + 'two amplitude criteria cannot be reconciled on these data';
      exhausted = false;
      break;
    }
    const rhoNext = Math.min(Math.max(rho * concordance, rhoMin), rhoMax);
    if (rhoNext === rho) { // pinned at a bound: no progress possible
      reason = 'the fixed-point step leaves the physical density range '
        + `[${fmtG(rhoMin)}, ${fmtG(rhoMax)}] A^-3 (next rho0 = ${fmtP(rho * concordance, 4)} A^-3)`;
      exhausted = false;
      break;
    }
    rho = rhoNext;
  }
  if (exhausted) {
    reason = `no concordant density within ${maxIter} passes (last concordance `
      + `${fmtP(history[history.length - 1][3], 4)})`;
  }
  const last = history[history.length - 1];
  return {
    rho0: last[0],
    converged,
    iterations: history.length,
    concordance: last[3],
    aDensity: last[1],
    aFz: last[2],
    // Judged on the first measured Q, not config.qmin (scaling.estimate_rho0).
    extrapolated: qFirst > FZ_FIT_WIDTH,
    qFirst,
    history,
    stopped,
    reason: converged ? null : reason,
  };
};

/**
 * Why a density self-consistency (estimateRho0) did not converge, for the
 * worker and the page to throw: the iteration either stopped at a trial
 * density the auto-scale could not fit (estimate.stopped), or ran out with
 * the two amplitude criteria still disagreeing. scaling_cli.py words the same
 * two cases for the CLI.
 */
export const rho0NonConvergenceMessage = (estimate) => {
  const concordance = Number(estimate.concordance.toPrecision(3));
  let reason;
  if (estimate.stopped) {
    reason = `the iteration stopped at a density the auto-scale cannot fit (${estimate.stopped})`;
  } else if (estimate.reason) {
    reason = estimate.reason;
  } else {
    reason = 'the density-limit and Q→0 Faber-Ziman amplitudes disagree at every density — '
      + 'typically data missing structure below Qmin';
  }
  return `ρ₀ self-consistency did not converge (final concordance ${concordance}): ${reason}. `
    + 'Set ρ₀ explicitly (value, data header, or mass density) and consider the '
    + 'Faber-Ziman Q→0 amplitude criterion for the scale.';
};

export const diagnosticsSummary = (result, config) => {
  const [lo, hi] = result.rFitWindowUsed || rFitWindow(config);
  let total = 0;
  let count = 0;
  for (let i = 0; i < result.r.length; i += 1) {
    if (result.r[i] >= lo && result.r[i] <= hi) {
      total += result.gFiltered[i];
      count += 1;
    }
  }
  const gWindowMean = total / count;
  const summary = {
    a: result.a,
    b: result.b,
    converged: result.converged,
    iterations: result.iterations,
    c1_tail_mean: result.c1TailMean,
    low_r_rms_pre_enforcement: result.lowRRms,
    g_window_mean: gWindowMean,
    r_fit_window: [lo, hi],
    gk_low_r_theory: -config.bAvgSq,
    d_r_low_r_slope_theory: -4 * Math.PI * config.rho0 * config.bAvgSq,
    density_limit_satisfied: Math.abs(gWindowMean) < 0.1,
  };
  if (result.fitFailure) summary.fit_failure = result.fitFailure;
  if (result.r0Detected != null) {
    summary.r0_detected = result.r0Detected;
    summary.window_refined = Boolean(result.windowRefined);
    if (!summary.window_refined && config.r0 != null) {
      // A given r0 is respected, but a first shell detected clearly below it
      // means the low-r window may cut into that shell (Python parity).
      summary.first_shell_below_r0 = result.r0Detected < config.r0 - 0.1;
    }
  }
  if (result.sweep) {
    summary.level = result.sweep.level;
    summary.level_uncertainty = result.sweep.levelUncertainty;
    summary.level_window = [result.sweep.qLo, result.sweep.qHi];
    summary.asymptote_found = result.sweep.asymptoteFound;
  }
  if (result.aFz != null) {
    summary.a_fz = result.aFz;
    if (config.amplitudeCriterion !== 'fz') {
      summary.amplitude_concordance = result.aFz / result.a;
      summary.amplitudes_concordant = Math.abs(result.aFz / result.a - 1) < 0.1;
    }
  }
  if (config.bSqAvg != null) {
    summary.fk_qmin = result.fk[0];
    summary.fk_q0_theory = -config.bSqAvg;
  }
  return summary;
};

// ---------------------------------------------------------------------------
// STOG file parsing / writing (port of the rmc_toolkits.parsers additions)
// ---------------------------------------------------------------------------

const NUMERIC_TOKEN = /^[+-]?((\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?|nan|inf(inity)?)$/i;

const tokenToFloat = (token) => {
  if (/^[+-]?nan$/i.test(token)) return NaN;
  if (/^[+-]?inf(inity)?$/i.test(token)) return token.startsWith('-') ? -Infinity : Infinity;
  return Number(token.replace(/[dD]/, 'e'));
};

export const readStogXy = (text) => {
  const groups = new Map();
  for (const line of text.split(/\r?\n/)) {
    const parts = line.trim().split(/\s+/).filter(Boolean);
    if (parts.length < 2) continue;
    if (!parts.every((token) => NUMERIC_TOKEN.test(token))) continue;
    const values = parts.map(tokenToFloat);
    if (!groups.has(values.length)) groups.set(values.length, []);
    groups.get(values.length).push(values);
  }
  if (!groups.size) throw new Error('file does not contain STOG numeric rows');
  let rows = null;
  for (const list of groups.values()) {
    if (!rows || list.length > rows.length) rows = list;
  }
  const nCols = rows[0].length;
  return Array.from({ length: nCols }, (_, col) => Float64Array.from(rows, (row) => row[col]));
};

const stogFlag = (token) => token.trim().toUpperCase().startsWith('Y');

export const readStogInp = (text) => {
  const lines = text.split(/\r?\n/).map((line) => line.trim()).filter(Boolean);
  if (lines.length < 22) throw new Error(`stog input has ${lines.length} non-empty lines; expected >= 22`);
  const nFiles = parseInt(lines[0].split(/\s+/)[0], 10);
  if (nFiles !== 1) throw new Error('only single-dataset stog inputs supported');
  const [qmin, qmax] = lines[2].split(/\s+/).slice(0, 2).map(Number);
  const [yoffset, yscale] = lines[3].split(/\s+/).slice(0, 2).map(Number);
  if (yscale === 0 || !isNum(yscale) || !isNum(yoffset)) {
    throw new Error('invalid yoffset/yscale line: yscale must be finite and nonzero');
  }
  if (Number(lines[4].split(/\s+/)[0]) !== 0) throw new Error('nonzero Q offset not supported');
  if (Number(lines[11].split(/\s+/)[0]) !== 0) throw new Error('nonzero second y offset not supported');
  if (stogFlag(lines[12])) throw new Error("interactive 'try again' loops not supported");
  if (!stogFlag(lines[13])) throw new Error('only filter-enabled stog inputs supported');
  const [peakCutoff, peakRmin, peakRmax] = lines[21].split(/\s+/).slice(0, 3).map(Number);
  return {
    dataFile: lines[1],
    qmin,
    qmax,
    yoffset,
    yscale,
    a: 1 / yscale,
    b: yoffset,
    outSq: lines[5],
    outGr: lines[6],
    rmax: Number(lines[7].split(/\s+/)[0]),
    nr: parseInt(lines[8].split(/\s+/)[0], 10),
    lorch: stogFlag(lines[9]),
    rho0: Number(lines[10].split(/\s+/)[0]),
    rCutoff: Number(lines[14].split(/\s+/)[0]),
    outFtSq: lines[15],
    outFtGr: lines[16],
    bAvgSq: Number(lines[17].split(/\s+/)[0]),
    outRmcFq: lines[18],
    outRmcGr: lines[19],
    outRmcDr: lines[20],
    peakCutoff,
    peakRmin,
    peakRmax,
  };
};

/**
 * Closest-approach proxy from a classic stog.inp first-peak line (port of
 * scaling_cli.stog_inp_closest_approach — keep in sync): the line zeroes g for
 * r <= peakCutoff except inside [peakRmin, peakRmax], so the asserted g = 0
 * region ends at peakRmin when a genuine window starts below the cutoff,
 * else at peakCutoff. null when it leaves no default fit window above rCutoff.
 */
export const stogInpClosestApproach = (inp, rCutoff) => {
  let candidate = inp.peakCutoff;
  if (inp.peakRmin > 0 && inp.peakRmin < inp.peakCutoff && inp.peakRmax > inp.peakRmin) {
    candidate = inp.peakRmin;
  }
  return candidate - R0_WINDOW_MARGIN > rCutoff + 0.2 ? candidate : null;
};

/**
 * Enforcement descriptor for the Auto StoG page, with the CLI's precedence
 * (scaling_cli._resolve_enforcement): null (off) | 'auto' (the worker
 * enforces at the first-shell foot) | {cutoff, peakRmin, peakRmax}. The page
 * pre-fills its Cutoff field with the loaded stog.inp's peak cutoff, so a
 * cutoff equal to inp.peakCutoff keeps the inp's first-peak window (classic
 * first_peak_zero semantics, as the CLI does without --enforce-cutoff); any
 * other typed cutoff is a flat replacement below it.
 */
export const resolveEnforcementDescriptor = ({ enforce, cutoff }, inp) => {
  if (!enforce) return null;
  const value = cutoff ?? (inp ? inp.peakCutoff : undefined);
  if (value === undefined || value === null) return 'auto';
  const usingInpWindow = Boolean(inp) && value === inp.peakCutoff;
  return {
    cutoff: value,
    peakRmin: usingInpWindow ? inp.peakRmin : value,
    peakRmax: usingInpWindow ? inp.peakRmax : value,
  };
};

export const readDatHeader = (text) => {
  const raw = {};
  for (const line of text.split(/\r?\n/)) {
    const idx = line.indexOf('::');
    if (idx < 0) continue;
    raw[line.slice(0, idx).trim().toUpperCase()] = line.slice(idx + 2).trim();
  }
  const result = { raw };
  if (raw.TITLE) result.title = raw.TITLE;
  if (raw.NUMBER_DENSITY) {
    for (const token of raw.NUMBER_DENSITY.split(/\s+/)) {
      const value = Number(token);
      if (isNum(value)) { result.numberDensity = value; break; }
    }
  }
  if (raw.MINIMUM_DISTANCES) {
    const values = raw.MINIMUM_DISTANCES.split(/\s+/).map(Number).filter(isNum);
    if (values.length) result.minDistance = Math.min(...values);
  }
  return result;
};

/** Fortran-style %.16E used by write_stog_xy (matches the Python writer). */
const fortranE = (value) => {
  if (!isNum(value)) return String(value).toUpperCase();
  let text = value.toExponential(16).toUpperCase();
  return text.replace(/E([+-])(\d)$/, 'E$10$2');
};

export const writeStogXy = (x, y, { title = '', extra = null } = {}) => {
  const lines = [String(x.length).padStart(12), title];
  for (let i = 0; i < x.length; i += 1) {
    let row = `  ${fortranE(x[i])}  ${fortranE(y[i])}`;
    if (extra) row += `  ${fortranE(extra[i])}`;
    lines.push(row);
  }
  return `${lines.join('\n')}\n`;
};

// ---------------------------------------------------------------------------
// Faber-Ziman coefficients (port of rmc_toolkits/scattering.py)
// ---------------------------------------------------------------------------

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
  Xe: 4.92, Y: 7.75, Yb: 12.43, Zn: 5.68, Zr: 7.16,
};

//: Standard atomic weights (u), generated from `periodictable` (CIAAW) for
//: exactly the COHERENT_B_FM element set (kept in sync with scattering.py).
export const ATOMIC_MASS_U = {
  Ag: 107.8682, Al: 26.98154, Am: 243.0, Ar: 39.95, As: 74.92159,
  Au: 196.96657, B: 10.81, Ba: 137.327, Be: 9.01218, Bi: 208.9804,
  Br: 79.904, C: 12.011, Ca: 40.078, Cd: 112.414, Ce: 140.116,
  Cl: 35.45, Co: 58.93319, Cr: 51.9961, Cs: 132.90545, Cu: 63.546,
  Dy: 162.5, Er: 167.259, Eu: 151.964, F: 18.9984, Fe: 55.845,
  Ga: 69.723, Gd: 157.25, Ge: 72.63, H: 1.008, He: 4.0026,
  Hf: 178.486, Hg: 200.592, Ho: 164.93033, I: 126.90447, In: 114.818,
  Ir: 192.217, K: 39.0983, Kr: 83.798, La: 138.90547, Li: 6.94,
  Lu: 174.9668, Mg: 24.305, Mn: 54.93804, Mo: 95.95, N: 14.007,
  Na: 22.98977, Nb: 92.90637, Nd: 144.242, Ne: 20.1797, Ni: 58.6934,
  Np: 237.0, O: 15.999, Os: 190.23, P: 30.97376, Pa: 231.03588,
  Pb: 207.2, Pd: 106.42, Pm: 145.0, Pr: 140.90766, Pt: 195.084,
  Ra: 226.0, Rb: 85.4678, Re: 186.207, Rh: 102.90549, Ru: 101.07,
  S: 32.06, Sb: 121.76, Sc: 44.95591, Se: 78.971, Si: 28.085,
  Sm: 150.36, Sn: 118.71, Sr: 87.62, Ta: 180.94788, Tb: 158.92535,
  Tc: 98.0, Te: 127.6, Th: 232.0377, Ti: 47.867, Tl: 204.38,
  Tm: 168.93422, U: 238.02891, V: 50.9415, W: 183.84, Xe: 131.293,
  Y: 88.90584, Yb: 173.045, Zn: 65.38, Zr: 91.224,
};

const AVOGADRO_PER_A3 = 6.02214076e23 / 1.0e24;

export const molarMass = (composition) => {
  const counts = typeof composition === 'string' ? parseFormula(composition) : { ...composition };
  let mass = 0;
  let atoms = 0;
  for (const [element, count] of Object.entries(counts)) {
    if (!(element in ATOMIC_MASS_U)) throw new Error(`no atomic mass for element '${element}'`);
    mass += ATOMIC_MASS_U[element] * count;
    atoms += count;
  }
  if (!(atoms > 0)) throw new Error('composition counts must sum to a positive number');
  return { mass, atoms };
};

/** ADDIE convention: g/cm^3 -> atoms/A^3 (rho0 = rho_m * N_A/1e24 * n/M). */
export const numberDensityFromMassDensity = (composition, massDensityGCm3) => {
  if (!(massDensityGCm3 > 0)) throw new Error('mass density must be positive');
  const { mass, atoms } = molarMass(composition);
  return massDensityGCm3 * AVOGADRO_PER_A3 * atoms / mass;
};

export const massDensityFromNumberDensity = (composition, numberDensityPerA3) => {
  if (!(numberDensityPerA3 > 0)) throw new Error('number density must be positive');
  const { mass, atoms } = molarMass(composition);
  return numberDensityPerA3 * mass / atoms / AVOGADRO_PER_A3;
};

export const parseFormula = (formula) => {
  const compact = formula.replace(/\s+/g, '');
  if (!compact) throw new Error('empty formula');
  const parseGroup = (start) => {
    const counts = {};
    let i = start;
    while (i < compact.length) {
      const char = compact[i];
      if (char === ')') return { counts, end: i };
      if (char === '(') {
        const inner = parseGroup(i + 1);
        if (inner.end >= compact.length || compact[inner.end] !== ')') {
          throw new Error(`unbalanced parentheses in '${formula}'`);
        }
        let j = inner.end + 1;
        const multMatch = compact.slice(j).match(/^\d*\.?\d*/);
        const mult = multMatch[0] ? Number(multMatch[0]) : 1;
        j += multMatch[0].length;
        for (const [el, count] of Object.entries(inner.counts)) {
          counts[el] = (counts[el] || 0) + count * mult;
        }
        i = j;
        continue;
      }
      const match = compact.slice(i).match(/^([A-Z][a-z]?)(\d*\.?\d*)/);
      if (!match || !match[1]) throw new Error(`cannot parse formula '${formula}' at '${compact.slice(i)}'`);
      const count = match[2] ? Number(match[2]) : 1;
      counts[match[1]] = (counts[match[1]] || 0) + count;
      i += match[0].length;
    }
    return { counts, end: i };
  };
  const { counts, end } = parseGroup(0);
  if (end !== compact.length) throw new Error(`unbalanced parentheses in '${formula}'`);
  if (!Object.keys(counts).length) throw new Error(`no elements found in '${formula}'`);
  return counts;
};

export const faberZiman = (composition, bOverridesFm = null) => {
  const counts = typeof composition === 'string' ? parseFormula(composition) : { ...composition };
  const total = Object.values(counts).reduce((sum, value) => sum + value, 0);
  if (!(total > 0)) throw new Error('composition counts must sum to a positive number');
  const bValues = {};
  for (const el of Object.keys(counts)) {
    if (bOverridesFm && el in bOverridesFm) bValues[el] = bOverridesFm[el];
    else if (el in COHERENT_B_FM) bValues[el] = COHERENT_B_FM[el];
    else throw new Error(`no coherent scattering length for element '${el}'`);
  }
  let bAvg = 0;
  let bSqAvg = 0;
  for (const [el, count] of Object.entries(counts)) {
    const fraction = count / total;
    bAvg += fraction * bValues[el];
    bSqAvg += fraction * bValues[el] * bValues[el];
  }
  const bAvgSq = bAvg * bAvg;
  if (bAvgSq < 1e-4 * bSqAvg) {
    throw new Error('near-null-matrix composition: <b>^2 is negligible relative to <b^2>');
  }
  return {
    bAvgSqFm2: bAvgSq,
    bSqAvgFm2: bSqAvg,
    bAvgSqBarn: bAvgSq / 100,
    bSqAvgBarn: bSqAvg / 100,
  };
};
