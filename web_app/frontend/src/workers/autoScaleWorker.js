// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Static-mode Auto StoG worker: runs the JS engine off-thread. One message per
// job: { id, config, mode, a, b, q, sq, sigma, enforcement } -> { id, ok,
// result | error }. Large arrays travel as transferable Float64Array buffers.

import {
  autoscale,
  diagnosticsSummary,
  firstPeakZero,
  fqToGpdf,
  makeConfig,
  scalePipeline,
} from './autoScale';

self.onmessage = (event) => {
  const { id, config: rawConfig, mode, a, b, q, sq, sigma, enforcement } = event.data;
  try {
    const config = makeConfig(rawConfig);
    const qArr = new Float64Array(q);
    const sqArr = new Float64Array(sq);
    const sigmaArr = sigma ? new Float64Array(sigma) : null;
    const result = mode === 'manual'
      ? scalePipeline(qArr, sqArr, config, a, b, { mode: 'manual' })
      : autoscale(qArr, sqArr, config, sigmaArr);
    const summary = diagnosticsSummary(result, config);

    // Unfiltered g(r) - 1 (the classic scale.gr): recomputed here exactly like
    // the CLI writer does — it is not part of the engine result.
    const fqScaled = new Float64Array(result.q.length);
    for (let i = 0; i < result.q.length; i += 1) {
      fqScaled[i] = result.q[i] * (result.sqScaled[i] - 1);
    }
    const gpdfUnfiltered = fqToGpdf(result.q, fqScaled, result.r, {
      lorch: config.lorch,
      lowQCorrection: config.lowQCorrection,
    });
    const gm1Unfiltered = new Float64Array(result.r.length);
    for (let i = 0; i < result.r.length; i += 1) {
      gm1Unfiltered[i] = gpdfUnfiltered[i] / (4 * Math.PI * config.rho0 * result.r[i]);
    }

    let gkEnforced = null;
    let drEnforced = null;
    if (enforcement) {
      const gFinal = firstPeakZero(result.r, result.gFiltered, enforcement);
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
        gm1Unfiltered: gm1Unfiltered.buffer,
        gkEnforced: gkEnforced ? gkEnforced.buffer : null,
        drEnforced: drEnforced ? drEnforced.buffer : null,
      },
    };
    const transfers = [
      payload.result.q, payload.result.sqRaw, payload.result.sqScaled,
      payload.result.sqFiltered, payload.result.sqFt, payload.result.r,
      payload.result.gk, payload.result.dr, payload.result.fk,
      payload.result.gm1Unfiltered,
    ];
    if (payload.result.gkEnforced) transfers.push(payload.result.gkEnforced, payload.result.drEnforced);
    self.postMessage(payload, transfers);
  } catch (error) {
    self.postMessage({ id, ok: false, error: error?.message || String(error) });
  }
};
