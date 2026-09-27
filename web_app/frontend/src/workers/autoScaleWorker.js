// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Static-mode Auto StoG worker: runs the JS engine off-thread. One message per
// job: { id, config, mode, a, b, q, sq, sigma, enforcement } -> { id, ok,
// result | error }. Large arrays travel as transferable Float64Array buffers.
// kind: 'estimateRho0' runs only the density self-consistency and returns the
// small estimate dict; a scaling job with estimateRho0: true runs it first and
// adopts the estimated density for the fit.

import {
  autoEnforcementCutoff,
  autoscale,
  detectFirstPeakOnset,
  diagnosticsSummary,
  effectiveS0Target,
  estimateRho0,
  firstPeakZero,
  fqToGpdf,
  makeConfig,
  rho0NonConvergenceError,
  scalePipeline,
  usableSigma,
  validateEnforcement,
} from './autoScale';
import { isBlankRequestValue, requestNumber, requestObject } from './requestGuards.js';

// app._resolve_scaling_mode's rules: mode 'auto' | 'manual'; manual needs a
// finite, non-zero scale a (S_corr = a*S_meas + b; a = 0 discards the data)
// and a finite offset b (default 0).
const resolveMode = (mode, a, b) => {
  const name = mode ?? 'auto';
  if (name !== 'auto' && name !== 'manual') {
    throw new Error(`mode must be 'auto' or 'manual', got '${name}'`);
  }
  if (name !== 'manual') return { mode: name, a: 0, b: 0 };
  if (isBlankRequestValue(a)) throw new Error("manual mode requires a scale 'a'");
  const scale = requestNumber(a, 'a');
  if (scale === 0) throw new Error(`manual mode requires a finite, non-zero scale 'a', got ${scale}`);
  return { mode: name, a: scale, b: requestNumber(b, 'b', { fallback: 0 }) };
};

// app._require_finite_scaling (+ the preview's raw-S(Q) check): a result whose
// a, b or curves hold NaN/Infinity means float64 overflowed on extreme input.
// It is an error, never ok: true with NaN curves for the page to plot or export.
const requireFiniteScaling = (result) => {
  const arrays = [result.sqScaled, result.sqFiltered, result.sqFt, result.gFiltered,
    result.gk, result.dr, result.fk, result.sqRaw];
  const finite = (array) => {
    if (!array) return true;
    for (let i = 0; i < array.length; i += 1) if (!Number.isFinite(array[i])) return false;
    return true;
  };
  if (!(Number.isFinite(result.a) && Number.isFinite(result.b)) || !arrays.every(finite)) {
    throw new Error(
      `the scaling result contains NaN or Infinity (a = ${result.a}, b = ${result.b}): `
      + 'the scale or offset is too extreme for float64 arithmetic'
    );
  }
};

/**
 * Run one job; returns `{ message, transfers }` to post, or throws. Exported so
 * tests can drive it without a Worker.
 */
export const runAutoScaleJob = (data) => {
  const {
    id, kind, config: rawConfig, mode: rawMode, a: rawA, b: rawB, q, sq, sigma, enforcement,
    estimateRho0: wantEstimate,
  } = requestObject(data);
  if (kind != null && kind !== 'estimateRho0') {
    throw new Error(`unknown request kind '${kind}' (expected estimateRho0, or none for a scaling job)`);
  }
  const qArr = new Float64Array(q);
  const sqArr = new Float64Array(sq);
  // CLI/API parity: a σ column with any zero / negative / non-finite value on
  // a usable row is dropped as a whole (the page also warns about it).
  const sigmaArr = sigma ? usableSigma(qArr, sqArr, new Float64Array(sigma)).sigma : null;
  if (kind === 'estimateRho0') {
    const estimate = estimateRho0(qArr, sqArr, makeConfig(rawConfig), sigmaArr);
    return { message: { id, ok: true, result: { estimate } }, transfers: [] };
  }
  const { mode, a, b } = resolveMode(rawMode, rawA, rawB);
  let config = makeConfig(rawConfig);
  // An explicit enforcement (not 'auto' / off) is checked before any work,
  // as the CLI and the API do (scaling.validate_enforcement).
  if (enforcement && typeof enforcement === 'object') validateEnforcement(enforcement, config.rmax);
  let rho0Estimate = null;
  if (wantEstimate) {
    rho0Estimate = estimateRho0(qArr, sqArr, config, sigmaArr);
    if (!rho0Estimate.converged) {
      // Never fit with a garbage density: surface the physics (or the
      // trial density that could not be fitted) instead.
      throw rho0NonConvergenceError(rho0Estimate);
    }
    config = { ...config, rho0: rho0Estimate.rho0 };
  }
  const result = mode === 'manual'
    ? scalePipeline(qArr, sqArr, config, a, b, { mode: 'manual' })
    : autoscale(qArr, sqArr, config, sigmaArr);
  requireFiniteScaling(result);
  if (mode !== 'manual' && result.fitFailure) {
    // A non-physical scale never becomes RMCProfile input (scaling_cli
    // refuse_failed_fit parity).
    throw new Error(
      `Auto-fit failed: ${result.fitFailure}. Check that the S(Q) is not `
      + 'sign-inverted or corrupted, set r₀ / the fit-window maximum, or use the '
      + 'Faber-Ziman Q→0 amplitude criterion when the composition is known.'
    );
  }
  // 'auto' enforcement is anchored on the data-derived first shell. Manual
  // runs skip the detection inside autoscale, so recover the onset here
  // exactly like the CLI does (scaling_cli.py post-run detection) —
  // otherwise a checked "Enforce low-r" would silently become a no-op.
  if (enforcement === 'auto' && result.r0Detected == null) {
    const onset = detectFirstPeakOnset(result.r, result.gFiltered, {
      searchMin: config.rCutoff + 0.3,
      qmax: config.qmax,
    });
    if (onset != null) result.r0Detected = onset;
  }
  const summary = diagnosticsSummary(result, config);

  // 'auto' enforces at the FOOT of the first shell (below its rising flank,
  // so no first-shell signal leaves the RMC files) — scaling_cli parity.
  let effectiveEnforcement = enforcement;
  if (enforcement === 'auto') {
    const cutoff = autoEnforcementCutoff(result.r, result.gFiltered, config, result.r0Detected);
    effectiveEnforcement = cutoff != null
      ? {
        cutoff,
        peakRmin: cutoff,
        peakRmax: cutoff,
        // Name the anchor like the CLI: a given r0 caps the onset, and
        // anchors the cutoff alone when no shell was detected.
        source: config.r0 != null && (result.r0Detected == null || config.r0 < result.r0Detected)
          ? 'auto (given r0)'
          : 'auto (first-shell foot)',
        firstShellOnset: result.r0Detected ?? null,
      }
      : null;
  }

  // Unfiltered g(r) (the classic scale.gr): recomputed here exactly like the
  // CLI writer does — it is not part of the engine result.
  const fqScaled = new Float64Array(result.q.length);
  for (let i = 0; i < result.q.length; i += 1) {
    fqScaled[i] = result.q[i] * (result.sqScaled[i] - 1);
  }
  const gpdfUnfiltered = fqToGpdf(result.q, fqScaled, result.r, {
    lorch: config.lorch,
    lowQCorrection: config.lowQCorrection,
    s0Target: effectiveS0Target(config),
  });
  const gUnfiltered = new Float64Array(result.r.length);
  for (let i = 0; i < result.r.length; i += 1) {
    gUnfiltered[i] = gpdfUnfiltered[i] / (4 * Math.PI * config.rho0 * result.r[i]) + 1;
  }

  let gkEnforced = null;
  let drEnforced = null;
  if (effectiveEnforcement) {
    const gFinal = firstPeakZero(result.r, result.gFiltered, effectiveEnforcement);
    gkEnforced = new Float64Array(result.r.length);
    drEnforced = new Float64Array(result.r.length);
    for (let i = 0; i < result.r.length; i += 1) {
      gkEnforced[i] = config.bAvgSq * (gFinal[i] - 1);
      drEnforced[i] = 4 * Math.PI * config.rho0 * result.r[i] * gkEnforced[i];
    }
  }

  const payload = {
    id,
    ok: true,
    result: {
      a: result.a,
      b: result.b,
      converged: result.converged,
      iterations: result.iterations,
      history: result.history,
      lowRRms: result.lowRRms,
      c1TailMean: result.c1TailMean,
      mode: result.mode,
      c1ModeEffective: result.c1ModeEffective,
      nDespiked: result.nDespiked,
      sweep: result.sweep,
      aFz: result.aFz,
      r0Detected: result.r0Detected,
      windowRefined: result.windowRefined,
      rFitWindowUsed: result.rFitWindowUsed,
      enforcement: effectiveEnforcement,
      rho0Estimate,
      rho0Used: config.rho0,
      summary,
      q: result.q.buffer,
      sqRaw: result.sqRaw.buffer,
      sqScaled: result.sqScaled.buffer,
      sqFiltered: result.sqFiltered.buffer,
      sqFt: result.sqFt.buffer,
      r: result.r.buffer,
      gk: result.gk.buffer,
      dr: result.dr.buffer,
      fk: result.fk.buffer,
      gUnfiltered: gUnfiltered.buffer,
      gkEnforced: gkEnforced ? gkEnforced.buffer : null,
      drEnforced: drEnforced ? drEnforced.buffer : null,
    },
  };
  const transfers = [
    payload.result.q, payload.result.sqRaw, payload.result.sqScaled,
    payload.result.sqFiltered, payload.result.sqFt, payload.result.r,
    payload.result.gk, payload.result.dr, payload.result.fk,
    payload.result.gUnfiltered,
  ];
  if (payload.result.gkEnforced) transfers.push(payload.result.gkEnforced, payload.result.drEnforced);
  return { message: payload, transfers };
};

// Guarded so the module can be imported by tests outside a worker context.
if (typeof self !== 'undefined' && typeof self.postMessage === 'function') {
  self.onmessage = (event) => {
    // The id is read defensively and the job validated inside the try, so a
    // null or malformed message still gets an answer and the page never hangs.
    const data = event?.data;
    const id = data !== null && typeof data === 'object' ? data.id : undefined;
    try {
      const { message, transfers } = runAutoScaleJob(data);
      self.postMessage(message, transfers);
    } catch (error) {
      // `summary` (additive, engine errors only): the page's one-line form.
      self.postMessage({ id, ok: false, error: error?.message || String(error), summary: error?.summary });
    }
  };
}
