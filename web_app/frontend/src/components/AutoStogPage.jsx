// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Auto StoG is *pre*-processing: it turns a measured S(Q) into the
// RMCProfile-ready PDF file family. Unlike the post-processing pages it is
// fully independent of the run folder — users upload the S(Q) here, the
// parity-tested engine runs in a Web Worker in both runtimes, and the outputs
// download as a zip. Files never leave the browser.

import React, { useEffect, useMemo, useRef, useState } from 'react';
import InteractivePlot from './InteractivePlot';
import { downloadBlob, sanitizeFilename } from '../figureExport';
import { buildZip } from '../zipArchive';
import {
  faberZiman,
  makeConfig,
  numberDensityFromMassDensity,
  readDatHeader,
  readStogInp,
  readStogXy,
  writeStogXy,
} from '../workers/autoScale';
import './AutoStogPage.css';

const STORAGE_KEY = 'autostog-session';

// Files that can seed a scaling session: classic stog inputs, plus data files
// that are not themselves stog *outputs* (so dropping a whole run folder does
// not offer scale.fq / ft.dat as sources).
const isInpCandidate = (name) => name.endsWith('.inp') || name === 'stog_input.dat';
const isDataCandidate = (name) => (
  /\.(sq|fq|dat)$/i.test(name)
  && !/^(scale[._]|ft\.dat)/i.test(name)
  && !/_rmc\.(fq|gr|dr)$/i.test(name)
);

const numberOr = (value) => {
  if (value === '' || value === null || value === undefined) return undefined;
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : undefined;
};

const fmt = (value, digits = 4) => (
  Number.isFinite(value) ? Number(value).toPrecision(digits) : '—'
);

const EMPTY_FORM = {
  qmin: '', qmax: '', rho0: '', massDensity: '', bAvgSq: '', bSqAvg: '', formula: '',
  rCutoff: '', r0: '', rFitMin: '', rFitMax: '', rmax: '', nr: '',
  c1Mode: 'sweep', amplitude: 'density',
  lorch: false, despike: false, robust: true, lowQCorrection: true, useSigma: true,
  enforce: true, enforceCutoff: '', manualA: '', manualB: '',
};

const OUTPUT_LIST = [
  ['sq_scaled', 'scaled S(Q)'],
  ['gr_unfiltered', 'unfiltered g(r)−1'],
  ['sq_filtered', 'filtered S(Q)'],
  ['gr_filtered', 'filtered g(r)−1 + D(r)'],
  ['rmc_fq', 'FK(Q) → RMCProfile'],
  ['rmc_gr', 'GK(r) → RMCProfile'],
  ['rmc_dr', 'D(r) → RMCProfile'],
  ['ft_correction', 'ft.dat correction'],
  ['provenance', 'provenance JSON'],
];

// Q extent / point count of the finite rows (NaN-padded rebin files are common).
const dataExtent = (q, sq, sigma) => {
  let qlo = null;
  let qhi = null;
  let count = 0;
  for (let i = 0; i < q.length; i += 1) {
    if (Number.isFinite(q[i]) && Number.isFinite(sq[i])) {
      if (qlo === null) qlo = q[i];
      qhi = q[i];
      count += 1;
    }
  }
  return { qlo, qhi, count, hasSigma: Boolean(sigma) };
};

// ---------------------------------------------------------------------------
// Config resolution: form > stog.inp > data header, exactly like the CLI.
// When rho0 cannot be resolved but a composition (⟨b²⟩) is available, the run
// proceeds with a seed and the worker replaces it by the self-consistent
// estimate (density-limit vs Q→0 amplitude concordance).
// ---------------------------------------------------------------------------

const resolveConfig = (form, inp, header) => {
  const pick = (formValue, fallback) => {
    const value = numberOr(formValue);
    return value === undefined ? fallback : value;
  };
  let bAvgSq = numberOr(form.bAvgSq) ?? (inp ? inp.bAvgSq : undefined);
  let bSqAvg = numberOr(form.bSqAvg);
  const formula = form.formula.trim();
  if (formula) {
    const coefficients = faberZiman(formula);
    if (bSqAvg === undefined) bSqAvg = coefficients.bSqAvgBarn;
    if (bAvgSq === undefined) bAvgSq = coefficients.bAvgSqBarn;
  }
  if (bAvgSq === undefined) throw new Error('⟨b⟩² unknown: give a composition, or set it under Advanced → Coefficients');
  const rCutoff = pick(form.rCutoff, inp ? inp.rCutoff : 1.0);
  let r0 = numberOr(form.r0);
  if (r0 === undefined && header?.minDistance != null) r0 = header.minDistance;
  if (r0 === undefined && inp) {
    const candidate = Math.max(inp.peakCutoff, inp.peakRmin);
    if (candidate - 0.25 > rCutoff + 0.2) r0 = candidate;
  }
  const qmin = pick(form.qmin, inp ? inp.qmin : undefined);
  const qmax = pick(form.qmax, inp ? inp.qmax : undefined);
  if (qmin === undefined || qmax === undefined) throw new Error('set the Q window (Qmin and Qmax)');
  let rho0 = pick(form.rho0, inp ? inp.rho0 : undefined);
  if (rho0 === undefined && header?.numberDensity != null) rho0 = header.numberDensity;
  if (rho0 === undefined) {
    const massDensity = numberOr(form.massDensity);
    if (massDensity !== undefined && formula) {
      rho0 = numberDensityFromMassDensity(formula, massDensity);
    }
  }
  let wantEstimate = false;
  if (rho0 === undefined) {
    if (bSqAvg !== undefined) {
      rho0 = 0.05; // seed only — the worker adopts the self-consistent estimate
      wantEstimate = true;
    } else {
      throw new Error(
        'number density unknown: set ρ₀ or mass density — or give a '
        + 'composition and Auto StoG estimates ρ₀ self-consistently'
      );
    }
  }
  const config = makeConfig({
    qmin, qmax, rho0, bAvgSq,
    bSqAvg: bSqAvg === undefined ? null : bSqAvg,
    rCutoff,
    r0: r0 === undefined ? null : r0,
    rFitMin: numberOr(form.rFitMin) ?? null,
    rFitMax: numberOr(form.rFitMax) ?? null,
    rmax: pick(form.rmax, inp ? inp.rmax : 50.0),
    nr: Math.round(pick(form.nr, inp ? inp.nr : 5000)),
    lorch: Boolean(form.lorch),
    lowQCorrection: Boolean(form.lowQCorrection),
    robust: Boolean(form.robust),
    c1Mode: form.c1Mode,
    amplitudeCriterion: form.amplitude,
    despike: Boolean(form.despike),
  });
  return { config, wantEstimate };
};

const resolveEnforcement = (form, inp) => {
  if (!form.enforce) return null;
  const cutoff = numberOr(form.enforceCutoff) ?? (inp ? inp.peakCutoff : undefined);
  if (cutoff === undefined) return 'auto'; // worker enforces at the detected r0
  const usingInpWindow = inp && numberOr(form.enforceCutoff) === undefined;
  return {
    cutoff,
    peakRmin: usingInpWindow ? inp.peakRmin : cutoff,
    peakRmax: usingInpWindow ? inp.peakRmax : cutoff,
  };
};

const AutoStogPage = () => {
  const [sources, setSources] = useState([]); // uploaded {name, text}
  const [selectedName, setSelectedName] = useState('');
  const [inspect, setInspect] = useState(null);
  const [form, setForm] = useState(EMPTY_FORM);
  const [advancedOpen, setAdvancedOpen] = useState(false);
  const [preview, setPreview] = useState(null);
  const [running, setRunning] = useState(false);
  const [estimating, setEstimating] = useState(false);
  const [rho0Info, setRho0Info] = useState(null); // last self-consistency estimate
  const [error, setError] = useState(null);
  const [exportStem, setExportStem] = useState('');
  const [exportResult, setExportResult] = useState(null);
  const [dragActive, setDragActive] = useState(false);
  const workerRef = useRef(null);
  const jobRef = useRef(0);
  const dataRef = useRef(null); // { q, sq, sigma, inp, header, name }
  const fileInputRef = useRef(null);

  // ── session restore (parameters only; uploads cannot persist) ────────────
  useEffect(() => {
    try {
      const saved = JSON.parse(sessionStorage.getItem(STORAGE_KEY) || 'null');
      if (saved?.form) {
        setForm((current) => ({ ...current, ...saved.form }));
        setAdvancedOpen(Boolean(saved.advancedOpen));
      }
    } catch { /* corrupted state: start clean */ }
  }, []);
  useEffect(() => {
    try {
      sessionStorage.setItem(STORAGE_KEY, JSON.stringify({ form, advancedOpen }));
    } catch { /* storage full/blocked: persistence is best-effort */ }
  }, [form, advancedOpen]);

  useEffect(() => () => workerRef.current?.terminate(), []);

  // ── worker plumbing ──────────────────────────────────────────────────────
  const postJob = (message, transfers) => new Promise((resolve, reject) => {
    if (!workerRef.current) {
      workerRef.current = new Worker(
        new URL('../workers/autoScaleWorker.js', import.meta.url),
        { type: 'module' }
      );
    }
    const worker = workerRef.current;
    jobRef.current += 1;
    const id = jobRef.current;
    const onMessage = (event) => {
      if (event.data.id !== id) return;
      worker.removeEventListener('message', onMessage);
      if (!event.data.ok) reject(new Error(event.data.error));
      else resolve(event.data.result);
    };
    worker.addEventListener('message', onMessage);
    worker.postMessage({ id, ...message }, transfers);
  });

  const packedData = () => {
    const data = dataRef.current;
    const q = Float64Array.from(data.q);
    const sq = Float64Array.from(data.sq);
    const sigma = form.useSigma && data.sigma ? Float64Array.from(data.sigma) : null;
    return {
      payload: { q: q.buffer, sq: sq.buffer, sigma: sigma ? sigma.buffer : null },
      transfers: sigma ? [q.buffer, sq.buffer, sigma.buffer] : [q.buffer, sq.buffer],
    };
  };

  // ── uploads ──────────────────────────────────────────────────────────────
  const selectSource = (name, list) => {
    setSelectedName(name);
    setPreview(null);
    setExportResult(null);
    setError(null);
    setRho0Info(null);
    dataRef.current = null;
    if (!name) { setInspect(null); return; }
    try {
      const entry = list.find((file) => file.name === name);
      if (!entry) throw new Error('File is no longer available — upload it again');
      let inp = null;
      let dataText = entry.text;
      let dataName = entry.name;
      if (isInpCandidate(entry.name)) {
        inp = readStogInp(entry.text);
        const dataEntry = list.find((file) => file.name === inp.dataFile);
        if (!dataEntry) {
          throw new Error(`stog input references '${inp.dataFile}' — upload that file too`);
        }
        dataText = dataEntry.text;
        dataName = inp.dataFile;
      }
      const columns = readStogXy(dataText);
      const header = readDatHeader(dataText);
      const sigma = columns.length >= 3 ? columns[2] : null;
      dataRef.current = { q: columns[0], sq: columns[1], sigma, inp, header, name: dataName };
      const extent = dataExtent(columns[0], columns[1], sigma);
      setInspect({ kind: inp ? 'inp' : 'data', inp, header, dataFile: dataName, extent });
      setForm((current) => ({
        ...current,
        // Q window: stog.inp > kept form value > the data's finite extent.
        qmin: inp ? String(inp.qmin) : (current.qmin || (extent.qlo != null ? String(Number(extent.qlo.toPrecision(4))) : '')),
        qmax: inp ? String(inp.qmax) : (current.qmax || (extent.qhi != null ? String(Number(extent.qhi.toPrecision(4))) : '')),
        rho0: inp ? String(inp.rho0) : (header?.numberDensity != null ? String(header.numberDensity) : current.rho0),
        bAvgSq: inp ? String(inp.bAvgSq) : current.bAvgSq,
        rCutoff: inp ? String(inp.rCutoff) : current.rCutoff,
        r0: header?.minDistance != null ? String(header.minDistance) : current.r0,
        rmax: inp ? String(inp.rmax) : current.rmax,
        nr: inp ? String(inp.nr) : current.nr,
        lorch: inp ? Boolean(inp.lorch) : current.lorch,
        enforce: inp ? true : current.enforce,
        enforceCutoff: inp ? String(inp.peakCutoff) : current.enforceCutoff,
        manualA: inp ? String(inp.a) : current.manualA,
        manualB: inp ? String(inp.b) : current.manualB,
      }));
    } catch (err) {
      setInspect(null);
      dataRef.current = null;
      setError(err.message || 'Could not read the uploaded file');
    }
  };

  const ingestFiles = async (fileList) => {
    setError(null);
    const accepted = [];
    for (const file of Array.from(fileList || [])) {
      if (isInpCandidate(file.name) || isDataCandidate(file.name)) {
        accepted.push({ name: file.name, text: await file.text() });
      }
    }
    if (!accepted.length) {
      setError('No usable files: upload a rebinned S(Q) (.sq / .fq / .dat, 2–3 columns), optionally with its classic stog.inp.');
      return;
    }
    const merged = [
      ...sources.filter((existing) => !accepted.some((file) => file.name === existing.name)),
      ...accepted,
    ].sort((left, right) => {
      const leftInp = isInpCandidate(left.name) ? 0 : 1;
      const rightInp = isInpCandidate(right.name) ? 0 : 1;
      return leftInp - rightInp || left.name.localeCompare(right.name);
    });
    setSources(merged);
    const preferred = accepted.find((file) => isInpCandidate(file.name)) || accepted[0];
    selectSource(preferred.name, merged);
  };

  const onDrop = (event) => {
    event.preventDefault();
    setDragActive(false);
    ingestFiles(event.dataTransfer?.files);
  };

  // ── scaling runs ─────────────────────────────────────────────────────────
  const runScaling = async (mode) => {
    const data = dataRef.current;
    if (!data) return;
    setRunning(true);
    setError(null);
    setExportResult(null);
    try {
      if (mode === 'manual' && numberOr(form.manualA) === undefined && !data.inp) {
        throw new Error('fixed scaling needs a (and optionally b)');
      }
      const { config, wantEstimate } = resolveConfig(form, data.inp, data.header);
      const enforcement = resolveEnforcement(form, data.inp);
      const { payload, transfers } = packedData();
      const raw = await postJob({
        config,
        mode,
        estimateRho0: wantEstimate,
        a: mode === 'manual' ? numberOr(form.manualA) ?? (data.inp ? data.inp.a : undefined) : undefined,
        b: mode === 'manual' ? numberOr(form.manualB) ?? (data.inp ? data.inp.b : 0) : undefined,
        enforcement,
        ...payload,
      }, transfers);
      // The worker may have replaced the rho0 seed by the self-consistent
      // estimate — everything downstream (guides, provenance) uses the
      // effective config.
      const effective = { ...config, rho0: raw.rho0Used ?? config.rho0 };
      if (raw.rho0Estimate) {
        setRho0Info(raw.rho0Estimate);
        setForm((current) => ({
          ...current,
          rho0: String(Number(raw.rho0Estimate.rho0.toPrecision(5))),
        }));
      }
      const arr = (buffer) => (buffer ? Array.from(new Float64Array(buffer)) : null);
      setPreview({
        kind: data.inp ? 'inp' : 'data',
        mode,
        inp: data.inp,
        header: data.header,
        result: {
          a: raw.a, b: raw.b, converged: raw.converged, iterations: raw.iterations,
          lowRRms: raw.lowRRms, c1TailMean: raw.c1TailMean, history: raw.history,
        },
        diagnostics: raw.summary,
        enforcement: raw.enforcement,
        rho0Estimate: raw.rho0Estimate || null,
        config: effective,
        guides: {
          asymptote: 1.0,
          gkLowR: -effective.bAvgSq,
          drSlope: -4 * Math.PI * effective.rho0 * effective.bAvgSq,
          s0Target: effective.bSqAvg == null ? null : 1 - effective.bSqAvg / effective.bAvgSq,
          level: raw.sweep ? raw.sweep.level : null,
          levelWindow: raw.sweep ? [raw.sweep.qLo, raw.sweep.qHi] : null,
          rFitWindow: raw.rFitWindowUsed,
          r0Detected: raw.r0Detected,
        },
        series: {
          q: arr(raw.q), sqRaw: arr(raw.sqRaw), sqScaled: arr(raw.sqScaled),
          sqFiltered: arr(raw.sqFiltered), sqFt: arr(raw.sqFt), r: arr(raw.r),
          gk: arr(raw.gk), dr: arr(raw.dr), fk: arr(raw.fk),
          gm1Unfiltered: arr(raw.gm1Unfiltered),
          gkEnforced: arr(raw.gkEnforced), drEnforced: arr(raw.drEnforced),
        },
      });
    } catch (err) {
      setPreview(null);
      setError(err.message || 'Scaling failed');
    } finally {
      setRunning(false);
    }
  };

  const runEstimate = async () => {
    const data = dataRef.current;
    if (!data) return;
    setEstimating(true);
    setError(null);
    try {
      const { config } = resolveConfig(form, data.inp, data.header);
      const { payload, transfers } = packedData();
      const result = await postJob({ kind: 'estimateRho0', config, ...payload }, transfers);
      setRho0Info(result.estimate);
      setForm((current) => ({
        ...current,
        rho0: String(Number(result.estimate.rho0.toPrecision(5))),
      }));
    } catch (err) {
      setError(err.message || 'Could not estimate the density');
    } finally {
      setEstimating(false);
    }
  };

  // ── export (zip of the classic stog output family) ───────────────────────
  const writeFiles = () => {
    if (!preview) return;
    setError(null);
    try {
      const { series, result, diagnostics, config, enforcement, rho0Estimate } = preview;
      const label = `rmc-autoscale (browser): a=${result.a.toPrecision(8)} b=${result.b.toPrecision(8)}`;
      const stem = sanitizeFilename(
        exportStem.trim() || dataRef.current?.name?.replace(/\.[^.]+$/, '') || 'autoscale'
      );
      const gm1 = series.gk.map((value) => value / config.bAvgSq);
      const encoder = new TextEncoder();
      const entries = [
        [`${stem}.sq`, writeStogXy(series.q, series.sqScaled, { title: label })],
        [`${stem}.gr`, writeStogXy(series.r, series.gm1Unfiltered || gm1, { title: label })],
        [`${stem}_ft.sq`, writeStogXy(series.q, series.sqFiltered, { title: label })],
        [`${stem}_ft.gr`, writeStogXy(series.r, gm1, {
          title: label,
          extra: series.r.map((radius, i) => 4 * Math.PI * config.rho0 * radius * gm1[i]),
        })],
        [`${stem}_rmc.fq`, writeStogXy(series.q, series.fk, { title: label })],
        [`${stem}_rmc.gr`, writeStogXy(series.r, series.gkEnforced || series.gk, { title: label })],
        [`${stem}_rmc.dr`, writeStogXy(series.r, series.drEnforced || series.dr, { title: label })],
        ['ft.dat', writeStogXy(series.q, series.sqFt, { title: label })],
        [`${stem}_provenance.json`, JSON.stringify({
          tool: 'rmc-autoscale (browser engine)',
          source: dataRef.current?.name,
          a: result.a,
          b: result.b,
          mode: preview.mode,
          enforcement,
          rho0Estimate,
          config,
          diagnostics,
        }, null, 2)],
      ].map(([name, text]) => ({ name, data: encoder.encode(text) }));
      downloadBlob(buildZip(entries), `${stem}_autoscale.zip`);
      setExportResult({ zip: `${stem}_autoscale.zip`, count: entries.length });
    } catch (err) {
      setError(err.message || 'Could not build the output zip');
    }
  };

  const setField = (key) => (event) => {
    const value = event.target.type === 'checkbox' ? event.target.checked : event.target.value;
    setForm((current) => ({ ...current, [key]: value }));
  };

  const fzMissing = form.amplitude === 'fz'
    && numberOr(form.bSqAvg) === undefined
    && !form.formula.trim();

  const canEstimate = inspect != null
    && (Boolean(form.formula.trim()) || numberOr(form.bSqAvg) !== undefined);

  // Live ⟨b⟩²/⟨b²⟩ + S(0) chip so the composition's role is visible up front.
  const compositionChip = useMemo(() => {
    const formula = form.formula.trim();
    if (!formula) return null;
    try {
      const fz = faberZiman(formula);
      const s0 = 1 - fz.bSqAvgBarn / fz.bAvgSqBarn;
      return `⟨b⟩² ${fz.bAvgSqBarn.toPrecision(4)} · ⟨b²⟩ ${fz.bSqAvgBarn.toPrecision(4)} barn · S(0) ${s0.toPrecision(3)}`;
    } catch {
      return null; // incomplete formula while typing
    }
  }, [form.formula]);

  // ── plot data ────────────────────────────────────────────────────────────
  const RMAX_DISPLAY = 8;

  const sqPlot = useMemo(() => {
    if (!preview) return null;
    const { series, guides } = preview;
    const qEnd = series.q[series.q.length - 1];
    const plots = [
      { label: 'Auto-scaled a·S + b', x: series.q, y: series.sqScaled },
      { label: 'Filtered S(Q)', x: series.q, y: series.sqFiltered },
      { label: 'Measured (unscaled)', x: series.q, y: series.sqRaw, defaultHidden: true },
      { label: 'S → 1', x: [series.q[0], qEnd], y: [1, 1], role: 'guide' },
    ];
    if (guides.level != null && guides.levelWindow) {
      plots.push({
        label: `Level ${fmt(guides.level, 5)}`,
        x: guides.levelWindow,
        y: [guides.level, guides.level],
        role: 'guide',
        color: '#4c7df0',
      });
    }
    if (guides.s0Target != null) {
      plots.push({
        label: `S(0) FZ target ${fmt(guides.s0Target, 3)}`,
        x: [series.q[0], Math.min(series.q[0] + 1.5, qEnd)],
        y: [guides.s0Target, guides.s0Target],
        role: 'guide',
        color: '#0c8599',
      });
    }
    return {
      title: 'S(Q) — scaled and filtered',
      xLabel: 'Q (Å^{-1})',
      yLabel: 'S(Q)',
      series: plots,
    };
  }, [preview]);

  // Full-range real-space plots; the theory guide lines only extend over the
  // low-r region where they mean something (below/around the first shell).
  const gkPlot = useMemo(() => {
    if (!preview) return null;
    const { series, guides } = preview;
    const level = guides.gkLowR;
    const guideEnd = Math.min(RMAX_DISPLAY, series.r[series.r.length - 1]);
    // The enforced curve (flat -<b>^2 below r0) IS the output that reaches
    // RMCProfile — show it as the primary series; the ripply pre-enforcement
    // fit stays one legend click away as the honesty diagnostic.
    const plots = series.gkEnforced
      ? [
        { label: 'G_K(r) output (RMC file)', x: series.r, y: series.gkEnforced },
        { label: 'Pre-enforcement fit', x: series.r, y: series.gk, defaultHidden: true },
      ]
      : [{ label: 'G_K(r) filtered', x: series.r, y: series.gk }];
    plots.push({ label: '−⟨b⟩² theory', x: [0, guideEnd], y: [level, level], role: 'guide' });
    return {
      title: 'G_K(r) — full range',
      xLabel: 'r (Å)',
      yLabel: 'G_K(r)',
      series: plots,
      // Default zoom keeps the theory level readable; box-zoom for detail.
      initialYDomain: [level * 2.1, -level * 3.2],
    };
  }, [preview]);

  const drPlot = useMemo(() => {
    if (!preview) return null;
    const { series, guides } = preview;
    const guideEnd = Math.min(RMAX_DISPLAY, series.r[series.r.length - 1]);
    const plots = series.drEnforced
      ? [
        { label: 'D(r) output (RMC file)', x: series.r, y: series.drEnforced },
        { label: 'Pre-enforcement fit', x: series.r, y: series.dr, defaultHidden: true },
      ]
      : [{ label: 'D(r)', x: series.r, y: series.dr }];
    plots.push({
      label: '−4πρ₀⟨b⟩²r theory',
      x: [0, guideEnd],
      y: [0, guides.drSlope * guideEnd],
      role: 'guide',
    });
    return {
      title: 'D(r) — full range',
      xLabel: 'r (Å)',
      yLabel: 'D(r)',
      series: plots,
    };
  }, [preview]);

  const diagnostics = preview?.diagnostics;
  const reference = preview?.inp;
  const trajectory = useMemo(() => {
    const history = preview?.result?.history;
    if (!history?.length || history.length < 2) return null;
    return history.slice(-6).map((entry) => fmt(entry[0], 5)).join(' → ');
  }, [preview]);

  return (
    <div className="autostog-page">
      <div className="autostog-controls">
        <div className="autostog-cluster autostog-cluster--source">
          <span className="autostog-label">DATA</span>
          <div
            className={`autostog-dropzone${dragActive ? ' is-drag' : ''}`}
            onDragOver={(event) => { event.preventDefault(); setDragActive(true); }}
            onDragLeave={() => setDragActive(false)}
            onDrop={onDrop}
          >
            <input
              ref={fileInputRef}
              type="file"
              multiple
              accept=".sq,.fq,.dat,.inp"
              className="visually-hidden"
              onChange={(event) => { ingestFiles(event.target.files); event.target.value = ''; }}
            />
            <button type="button" className="autostog-pill" onClick={() => fileInputRef.current?.click()}>
              Upload S(Q)…
            </button>
            <span className="autostog-dropzone-hint">or drop files (S(Q) ± stog.inp)</span>
          </div>
          {sources.length > 1 && (
            <select
              value={selectedName}
              onChange={(event) => selectSource(event.target.value, sources)}
              aria-label="Scaling source file"
            >
              {sources.map((file) => (
                <option key={file.name} value={file.name}>
                  {isInpCandidate(file.name) ? `⚙ ${file.name}` : file.name}
                </option>
              ))}
            </select>
          )}
          {inspect?.extent && (
            <span
              className="autostog-chip autostog-chip--file"
              title={`${inspect.dataFile}: ${inspect.extent.count} points, Q ${fmt(inspect.extent.qlo, 3)}–${fmt(inspect.extent.qhi, 4)} Å⁻¹${inspect.extent.hasSigma ? ', σ column present' : ''}`}
            >
              {inspect.dataFile}: {inspect.extent.count} pts
              · Q {fmt(inspect.extent.qlo, 3)}–{fmt(inspect.extent.qhi, 4)} Å⁻¹
              {inspect.extent.hasSigma ? ' · σ' : ''}
            </span>
          )}
          {inspect?.kind === 'inp' && (
            <span className="autostog-chip" title="Hand scaling recorded in the stog input">
              hand a = {fmt(inspect.inp.a, 4)}, b = {fmt(inspect.inp.b, 4)}
            </span>
          )}
        </div>

        <div className="autostog-cluster">
          <span className="autostog-label">SAMPLE</span>
          <label
            className="autostog-field autostog-field--formula"
            title="Chemical composition (neutron Sears table drives ⟨b⟩², ⟨b²⟩, and the S(0) target; for x-ray data set ⟨b²⟩ in Advanced instead)"
          >
            <span>Composition</span>
            <input value={form.formula} onChange={setField('formula')} placeholder="e.g. Mn3Sn" spellCheck="false" />
          </label>
          <label
            className="autostog-field"
            title="Number density in atoms/Å³. Resolved from: your value → NUMBER_DENSITY :: data header → mass density + composition → left empty with a composition, it is estimated self-consistently"
          >
            <span>ρ₀ Å⁻³</span>
            <input value={form.rho0} onChange={setField('rho0')} inputMode="decimal" placeholder="auto / estimate" />
          </label>
          <button
            type="button"
            className="autostog-pill"
            disabled={!canEstimate || estimating || running}
            onClick={runEstimate}
            title="Self-consistent density: iterate until the density-limit amplitude agrees with the composition's ρ₀-independent Q→0 Faber-Ziman amplitude"
          >
            {estimating ? 'Estimating…' : 'Estimate ρ₀'}
          </button>
          <label className="autostog-field" title="Alternative to ρ₀: mass density + composition convert via N_A (ADDIE convention)">
            <span>or ρ g/cm³</span>
            <input value={form.massDensity} onChange={setField('massDensity')} inputMode="decimal" />
          </label>
          {compositionChip && <span className="autostog-chip" title="Computed from the composition (Sears neutron lengths)">{compositionChip}</span>}
          {rho0Info && (
            <span
              className={`autostog-chip ${rho0Info.extrapolated ? 'autostog-chip--warn' : 'autostog-chip--good'}`}
              title={`Concordance a_fz/a_density = ${fmt(rho0Info.concordance, 5)} after ${rho0Info.iterations} passes${rho0Info.extrapolated ? ' — long Q→0 extrapolation: treat as a starting point' : ''}`}
            >
              ρ₀ ≈ {fmt(rho0Info.rho0, 5)} (self-consistent{rho0Info.extrapolated ? ', extrapolated' : ''})
            </span>
          )}
        </div>

        <div className="autostog-cluster">
          <span className="autostog-label">Q WINDOW</span>
          <label className="autostog-field" title="Fit/transform lower bound — data below is replaced by the analytic low-Q correction">
            <span>Qmin Å⁻¹</span>
            <input value={form.qmin} onChange={setField('qmin')} inputMode="decimal" />
          </label>
          <label className="autostog-field" title="Fit/transform upper bound — cut before the detector edge noise">
            <span>Qmax Å⁻¹</span>
            <input value={form.qmax} onChange={setField('qmax')} inputMode="decimal" />
          </label>
        </div>

        <div className="autostog-cluster autostog-cluster--actions">
          <button
            type="button"
            className="autostog-primary"
            disabled={!inspect || running || estimating || fzMissing}
            onClick={() => runScaling('auto')}
          >
            {running ? 'Fitting…' : 'Auto-scale'}
          </button>
          <button
            type="button"
            className={`autostog-pill${advancedOpen ? ' is-active' : ''}`}
            aria-expanded={advancedOpen}
            onClick={() => setAdvancedOpen((open) => !open)}
          >
            Advanced
          </button>
        </div>
      </div>

      {advancedOpen && (
        <div className="autostog-controls autostog-controls--advanced">
          <fieldset className="autostog-group">
            <legend>Amplitude &amp; offset</legend>
            <p className="autostog-group-desc">How the correction S′ = a·S + b is determined.</p>
            <div className="autostog-group-fields">
              <label className="autostog-field autostog-field--select" title="Level sweep (default): measure the flat high-Q level and tie b to it, leaving one amplitude dof. Joint: original 2-dof (a, b) fit">
                <span>High-Q</span>
                <select value={form.c1Mode} onChange={setField('c1Mode')}>
                  <option value="sweep">Level sweep</option>
                  <option value="joint">Joint 2-dof</option>
                </select>
              </label>
              <label className="autostog-field autostog-field--select" title="What pins the amplitude: the low-r density limit (default), or the composition's Q→0 Faber-Ziman limit S(0) = 1 − ⟨b²⟩/⟨b⟩²">
                <span>Amplitude</span>
                <select value={form.amplitude} onChange={setField('amplitude')}>
                  <option value="density">Density limit</option>
                  <option value="fz">Faber-Ziman Q→0</option>
                </select>
              </label>
              <label className={`autostog-toggle${form.robust ? ' is-on' : ''}`} title="Huber re-weighting so isolated outliers cannot drag the closed-form fit (recommended)">
                <input type="checkbox" checked={form.robust} onChange={setField('robust')} />
                Robust
              </label>
              <label className={`autostog-toggle${form.useSigma ? ' is-on' : ''}`} title="1/σ-weight the high-Q fit with the data file's third column, when present">
                <input type="checkbox" checked={form.useSigma} onChange={setField('useSigma')} />
                σ column
              </label>
              <label className={`autostog-toggle${form.despike ? ' is-on' : ''}`} title="Drop narrow rolling-median outliers before fitting. For detector glitches only — it also flags real Bragg maxima on crystalline data">
                <input type="checkbox" checked={form.despike} onChange={setField('despike')} />
                Despike
              </label>
              {fzMissing && (
                <span className="autostog-chip autostog-chip--danger">FZ needs ⟨b²⟩ or a composition</span>
              )}
            </div>
          </fieldset>

          <fieldset className="autostog-group">
            <legend>Coefficients</legend>
            <p className="autostog-group-desc">Override the composition-derived neutron values (required for x-ray data).</p>
            <div className="autostog-group-fields">
              <label className="autostog-field" title="⟨b⟩² = (Σ cᵢbᵢ)² in barns — the classic stog 'Faber-Ziman coefficient' (x-ray normalized data: 1)">
                <span>⟨b⟩² barn</span>
                <input value={form.bAvgSq} onChange={setField('bAvgSq')} inputMode="decimal" placeholder="from composition" />
              </label>
              <label className="autostog-field" title="⟨b²⟩ = Σ cᵢbᵢ² in barns — sets the S(0) limit; enables the FZ amplitude, ρ₀ estimation, and the Q→0 diagnostic (x-ray: ⟨Z²⟩/⟨Z⟩²)">
                <span>⟨b²⟩ barn</span>
                <input value={form.bSqAvg} onChange={setField('bSqAvg')} inputMode="decimal" placeholder="from composition" />
              </label>
            </div>
          </fieldset>

          <fieldset className="autostog-group">
            <legend>Transform</legend>
            <p className="autostog-group-desc">Fourier grid and filter — the defaults suit most data.</p>
            <div className="autostog-group-fields">
              <label className="autostog-field" title="Classic stog Fourier filter: below this radius g(r) is unphysical and the corresponding S(Q) correction is subtracted">
                <span>Filter r-cut Å</span>
                <input value={form.rCutoff} onChange={setField('rCutoff')} inputMode="decimal" placeholder="1.0" />
              </label>
              <label className="autostog-field" title="Real-space grid extent">
                <span>rmax Å</span>
                <input value={form.rmax} onChange={setField('rmax')} inputMode="decimal" placeholder="50" />
              </label>
              <label className="autostog-field" title="Number of r grid points">
                <span>r points</span>
                <input value={form.nr} onChange={setField('nr')} inputMode="numeric" placeholder="5000" />
              </label>
              <label className={`autostog-toggle${form.lorch ? ' is-on' : ''}`} title="Lorch window: damps termination ripples at the cost of real-space resolution">
                <input type="checkbox" checked={form.lorch} onChange={setField('lorch')} />
                Lorch
              </label>
              <label className={`autostog-toggle${form.lowQCorrection ? ' is-on' : ''}`} title="Analytic correction for the unmeasured [0, Qmin] range — keeps the fitted scale unbiased (disable only for strict classic-stog parity)">
                <input type="checkbox" checked={form.lowQCorrection} onChange={setField('lowQCorrection')} />
                Low-Q corr.
              </label>
            </div>
          </fieldset>

          <fieldset className="autostog-group">
            <legend>Low-r region</legend>
            <p className="autostog-group-desc">Below the first shell g(r) must vanish — where the density limit is read and outputs are cleaned.</p>
            <div className="autostog-group-fields">
              <label className="autostog-field" title="Closest interatomic approach. Empty: detected from the data (first-shell flank), or a MINIMUM_DISTANCES :: header">
                <span>r₀ approach Å</span>
                <input value={form.r0} onChange={setField('r0')} inputMode="decimal" placeholder="auto" />
              </label>
              <label className="autostog-field" title="Low-r fit window minimum (default: filter r-cut + 0.2 Å)">
                <span>Fit win min</span>
                <input value={form.rFitMin} onChange={setField('rFitMin')} inputMode="decimal" placeholder="auto" />
              </label>
              <label className="autostog-field" title="Low-r fit window maximum (default: r₀ − 0.25 Å)">
                <span>Fit win max</span>
                <input value={form.rFitMax} onChange={setField('rFitMax')} inputMode="decimal" placeholder="auto" />
              </label>
              <label className={`autostog-toggle${form.enforce ? ' is-on' : ''}`} title="Classic stog final step: replace the RMC outputs below the cutoff by the exact theoretical low-r values (ripple removal)">
                <input type="checkbox" checked={form.enforce} onChange={setField('enforce')} />
                Enforce low-r
              </label>
              {form.enforce && (
                <label className="autostog-field" title="Enforcement cutoff (empty: the detected r₀)">
                  <span>Cutoff Å</span>
                  <input value={form.enforceCutoff} onChange={setField('enforceCutoff')} inputMode="decimal" placeholder="auto" />
                </label>
              )}
            </div>
          </fieldset>

          <fieldset className="autostog-group">
            <legend>Fixed scaling</legend>
            <p className="autostog-group-desc">Skip the auto-fit: apply a hand (a, b), e.g. to reproduce a classic stog run.</p>
            <div className="autostog-group-fields">
              <label className="autostog-field" title="Fixed a in S′ = a·S + b">
                <span>a</span>
                <input value={form.manualA} onChange={setField('manualA')} inputMode="decimal" />
              </label>
              <label className="autostog-field" title="Fixed b in S′ = a·S + b">
                <span>b</span>
                <input value={form.manualB} onChange={setField('manualB')} inputMode="decimal" />
              </label>
              <button
                type="button"
                className="autostog-pill"
                disabled={!inspect || running}
                onClick={() => runScaling('manual')}
              >
                Run fixed (a, b)
              </button>
            </div>
          </fieldset>
        </div>
      )}

      {error && <div className="autostog-banner autostog-banner--danger">{error}</div>}

      <details className="autostog-explainer">
        <summary>How Auto StoG works</summary>
        <ol>
          <li><b>Inputs:</b> an uploaded S(Q), the composition, and the [Qmin, Qmax] window
            are all that is required. ρ₀ comes from the data header, mass density, your
            value — or is estimated self-consistently; every other parameter has a
            physics-based default.</li>
          <li><b>Composition →</b> ⟨b⟩², ⟨b²⟩ (Sears neutron lengths) and the Keen Eq. 21
            target S(0) = 1 − ⟨b²⟩/⟨b⟩², which also anchors the analytic correction for the
            unmeasured [0, Qmin] range.</li>
          <li><b>Level sweep:</b> every high-Q window is tested for statistical flatness;
            the measured level pins the offset (b = 1 − a·level).</li>
          <li><b>Amplitude:</b> the low-r density limit (g → 0 below the first shell,
            self-consistent with the Fourier filter), with the Q→0 Faber-Ziman limit as an
            independent criterion — their concordance is the absolute-scale trust metric.</li>
          <li><b>ρ₀ self-consistency:</b> the density-limit amplitude depends on ρ₀, the
            FZ amplitude does not — iterating ρ₀ until they agree recovers the density
            from the data (needs a composition; long Q→0 extrapolations are flagged).</li>
          <li><b>r₀ detection:</b> the first coordination shell is located from the data
            (|g| flank, robust to inverted negative-b shells) and refines the fit window and
            the classic low-r enforcement cutoff.</li>
          <li><b>Outputs:</b> the classic stog file family (S(Q), g(r)−1, filtered pair,
            F<sub>K</sub>(Q), G<sub>K</sub>(r), D(r), ft.dat) + a provenance JSON. Read the
            flags: a violated density limit means the absolute scale needs the composition
            (FZ) route or external validation.</li>
        </ol>
      </details>

      {!preview && !error && (
        <div className="autostog-empty">
          <h2>Automatic total-scattering scaling</h2>
          <p>
            Upload a rebinned S(Q) file (optionally with its classic <code>stog.inp</code>),
            enter the composition and Q window, then hit <b>Auto-scale</b>. The engine
            determines the scale and offset from the physics — the statistically flat
            high-Q level, the low-r density limit, and the composition&apos;s Faber-Ziman
            constraints — replacing the classic stog “try again” loop, and reports the
            honest fit quality before any enforcement. This page is independent of the
            run folder used by the other tabs: everything runs in your browser and files
            never leave your device.
          </p>
        </div>
      )}

      {preview && diagnostics && (
        <div className="autostog-readout">
          <div className="autostog-stat">
            <span className="autostog-stat-label">Correction</span>
            <span className="autostog-stat-value">a = {fmt(preview.result.a, 5)} · b = {fmt(preview.result.b, 5)}</span>
            <span className="autostog-stat-sub">
              {reference
                ? `hand: a = ${fmt(reference.a, 5)}, b = ${fmt(reference.b, 5)}`
                : preview.mode === 'manual' ? 'fixed by you' : 'auto-fit'}
            </span>
          </div>
          <div className="autostog-stat">
            <span className="autostog-stat-label">{preview.mode === 'manual' ? 'Mode' : 'Convergence'}</span>
            <span className="autostog-stat-value">
              {preview.mode === 'manual'
                ? 'fixed (a, b)'
                : `${preview.result.converged ? '✓' : '✗'} ${preview.result.iterations} iterations`}
            </span>
            <span className="autostog-stat-sub">
              {trajectory ? `a: ${trajectory}` : diagnostics.level != null
                ? `level ${fmt(diagnostics.level, 5)}${Number.isFinite(diagnostics.level_uncertainty) ? ` ± ${fmt(diagnostics.level_uncertainty, 2)}` : ''}`
                : ''}
            </span>
          </div>
          {preview.rho0Estimate && (
            <div className={`autostog-stat ${preview.rho0Estimate.extrapolated ? 'is-warn' : 'is-good'}`}>
              <span className="autostog-stat-label">ρ₀ self-consistency</span>
              <span className="autostog-stat-value">{fmt(preview.rho0Estimate.rho0, 5)} Å⁻³ (adopted)</span>
              <span className="autostog-stat-sub">
                concordance {fmt(preview.rho0Estimate.concordance, 4)}
                {preview.rho0Estimate.extrapolated ? ' · long Q→0 extrapolation' : ` · ${preview.rho0Estimate.iterations} passes`}
              </span>
            </div>
          )}
          {diagnostics.level != null && trajectory && (
            <div className="autostog-stat">
              <span className="autostog-stat-label">High-Q level</span>
              <span className="autostog-stat-value">{fmt(diagnostics.level, 5)}</span>
              <span className="autostog-stat-sub">
                {Number.isFinite(diagnostics.level_uncertainty) ? `± ${fmt(diagnostics.level_uncertainty, 2)} · ` : ''}
                Q ∈ [{fmt(diagnostics.level_window?.[0], 4)}, {fmt(diagnostics.level_window?.[1], 4)}]
              </span>
            </div>
          )}
          <div className="autostog-stat">
            <span className="autostog-stat-label">Fit quality</span>
            <span className="autostog-stat-value">low-r rms {fmt(diagnostics.low_r_rms_pre_enforcement, 3)}</span>
            <span className="autostog-stat-sub">C1 tail mean {fmt(diagnostics.c1_tail_mean, 5)}</span>
          </div>
          <div className={`autostog-stat ${diagnostics.density_limit_satisfied ? 'is-good' : 'is-bad'}`}>
            <span className="autostog-stat-label">Density limit</span>
            <span className="autostog-stat-value">{diagnostics.density_limit_satisfied ? 'satisfied' : 'NOT satisfiable'}</span>
            <span className="autostog-stat-sub">
              {diagnostics.density_limit_satisfied
                ? 'necessary, not sufficient'
                : 'absolute scale needs external validation'}
            </span>
          </div>
          {diagnostics.r0_detected != null && (
            <div className="autostog-stat">
              <span className="autostog-stat-label">First shell r₀</span>
              <span className="autostog-stat-value">{fmt(diagnostics.r0_detected, 4)} Å (detected)</span>
              <span className="autostog-stat-sub">
                {diagnostics.window_refined ? 'fit window refined to it' : 'window unchanged'}
                {preview.enforcement ? ` · enforced below ${fmt(preview.enforcement.cutoff ?? preview.enforcement[0], 3)} Å` : ''}
              </span>
            </div>
          )}
          {diagnostics.amplitude_concordance != null && (
            <div className={`autostog-stat ${diagnostics.amplitudes_concordant ? 'is-good' : 'is-warn'}`}>
              <span className="autostog-stat-label">Concordance</span>
              <span className="autostog-stat-value">a_fz / a = {fmt(diagnostics.amplitude_concordance, 3)}</span>
              <span className="autostog-stat-sub">
                {diagnostics.amplitudes_concordant ? 'independent criteria agree' : 'disagree — check ρ₀ / low-Q'}
              </span>
            </div>
          )}
        </div>
      )}

      {preview && (
        <div className="autostog-plots">
          {sqPlot && (
            <section className="autostog-card autostog-card--span">
              <header><h3>{sqPlot.title}</h3></header>
              <InteractivePlot file={{ path: 'autostog-sq', name: 'autostog-sq' }} variant="wide" plotData={sqPlot} />
            </section>
          )}
          {gkPlot && (
            <section className="autostog-card">
              <header>
                <h3>{gkPlot.title}</h3>
                <span>
                  fit window {fmt(preview.guides.rFitWindow?.[0], 3)}–{fmt(preview.guides.rFitWindow?.[1], 3)} Å
                  {' · double-click for full amplitude'}
                </span>
              </header>
              <InteractivePlot file={{ path: 'autostog-gk', name: 'autostog-gk' }} plotData={gkPlot} />
            </section>
          )}
          {drPlot && (
            <section className="autostog-card">
              <header><h3>{drPlot.title}</h3></header>
              <InteractivePlot file={{ path: 'autostog-dr', name: 'autostog-dr' }} plotData={drPlot} />
            </section>
          )}
        </div>
      )}

      {preview && (
        <div className="autostog-export">
          <span className="autostog-label">EXPORT</span>
          <label className="autostog-field autostog-field--wide">
            <span>Name stem</span>
            <input
              value={exportStem}
              onChange={(event) => setExportStem(event.target.value)}
              placeholder="default: data file stem"
              spellCheck="false"
            />
          </label>
          <button type="button" className="autostog-primary" onClick={writeFiles}>
            Download .zip (RMCProfile files)
          </button>
          {exportResult?.zip && (
            <span className="autostog-chip autostog-chip--good">
              {exportResult.zip} · {exportResult.count} files
            </span>
          )}
          <span className="autostog-export-note">
            {OUTPUT_LIST.length} files: {OUTPUT_LIST.map(([, label]) => label).join(' · ')}
          </span>
        </div>
      )}
    </div>
  );
};

export default AutoStogPage;
