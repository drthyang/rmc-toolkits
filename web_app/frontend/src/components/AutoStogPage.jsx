// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useEffect, useMemo, useRef, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import InteractivePlot from './InteractivePlot';
import { isStaticMode } from '../browserData';
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

// Files that can seed a scaling session: classic stog inputs first, then data.
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

// ---------------------------------------------------------------------------
// Static-mode helpers: resolve config exactly like the Flask backend does.
// ---------------------------------------------------------------------------

const resolveStaticConfig = (form, inp, header) => {
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
  if (bAvgSq === undefined) throw new Error('data mode requires ⟨b⟩²: set it or give a formula');
  const rCutoff = pick(form.rCutoff, inp ? inp.rCutoff : 1.0);
  let r0 = numberOr(form.r0);
  if (r0 === undefined && header?.minDistance != null) r0 = header.minDistance;
  if (r0 === undefined && inp) {
    const candidate = Math.max(inp.peakCutoff, inp.peakRmin);
    if (candidate - 0.25 > rCutoff + 0.2) r0 = candidate;
  }
  const qmin = pick(form.qmin, inp ? inp.qmin : undefined);
  const qmax = pick(form.qmax, inp ? inp.qmax : undefined);
  if (qmin === undefined || qmax === undefined) throw new Error('data mode requires Qmin and Qmax');
  let rho0 = pick(form.rho0, inp ? inp.rho0 : undefined);
  if (rho0 === undefined && header?.numberDensity != null) rho0 = header.numberDensity;
  if (rho0 === undefined) {
    const massDensity = numberOr(form.massDensity);
    if (massDensity !== undefined && formula) {
      rho0 = numberDensityFromMassDensity(formula, massDensity);
    }
  }
  if (rho0 === undefined) throw new Error('number density unknown: set ρ₀, or mass density with a composition');
  return makeConfig({
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
};

const resolveStaticEnforcement = (form, inp) => {
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

const AutoStogPage = ({ directory, localRun }) => {
  const staticMode = isStaticMode();
  const [files, setFiles] = useState([]);
  const [filesError, setFilesError] = useState(null);
  const [selectedKey, setSelectedKey] = useState('');
  const [inspect, setInspect] = useState(null);
  const [form, setForm] = useState(EMPTY_FORM);
  const [advancedOpen, setAdvancedOpen] = useState(false);
  const [preview, setPreview] = useState(null);
  const [running, setRunning] = useState(false);
  const [error, setError] = useState(null);
  const [exportForm, setExportForm] = useState({ outDir: '', outStem: '', force: false });
  const [exportResult, setExportResult] = useState(null);
  const [exporting, setExporting] = useState(false);
  const workerRef = useRef(null);
  const jobRef = useRef(0);
  const staticDataRef = useRef(null); // { q, sq, sigma, inp, header, name }
  const restoredRef = useRef(false);

  // ── session restore ──────────────────────────────────────────────────────
  useEffect(() => {
    try {
      const saved = JSON.parse(sessionStorage.getItem(STORAGE_KEY) || 'null');
      if (saved?.form) {
        setForm((current) => ({ ...current, ...saved.form }));
        setAdvancedOpen(Boolean(saved.advancedOpen));
        restoredRef.current = saved.selectedKey || false;
      }
    } catch { /* corrupted state: start clean */ }
  }, []);
  useEffect(() => {
    try {
      sessionStorage.setItem(
        STORAGE_KEY,
        JSON.stringify({ selectedKey, form, advancedOpen })
      );
    } catch { /* storage full/blocked: persistence is best-effort */ }
  }, [selectedKey, form, advancedOpen]);

  useEffect(() => () => workerRef.current?.terminate(), []);

  // ── candidate files ──────────────────────────────────────────────────────
  useEffect(() => {
    if (staticMode) {
      const candidates = (localRun?.files || [])
        .map((file) => ({ name: file.name, key: file.name, sourceFile: file.sourceFile }))
        .filter((file) => isInpCandidate(file.name) || isDataCandidate(file.name))
        .sort((left, right) => {
          const leftInp = isInpCandidate(left.name) ? 0 : 1;
          const rightInp = isInpCandidate(right.name) ? 0 : 1;
          return leftInp - rightInp || left.name.localeCompare(right.name);
        });
      setFiles(candidates);
      setFilesError(null);
      return undefined;
    }
    let cancelled = false;
    setFilesError(null);
    axios.get(`${API_BASE_URL}/api/files`, { params: { dir: directory || '.' } })
      .then((response) => {
        if (cancelled) return;
        const candidates = (response.data.files || [])
          .filter((item) => item.type === 'file')
          .filter((item) => isInpCandidate(item.name) || isDataCandidate(item.name))
          .map((item) => ({ name: item.name, key: item.path }))
          .sort((left, right) => {
            const leftInp = isInpCandidate(left.name) ? 0 : 1;
            const rightInp = isInpCandidate(right.name) ? 0 : 1;
            return leftInp - rightInp || left.name.localeCompare(right.name);
          });
        setFiles(candidates);
      })
      .catch((err) => {
        if (!cancelled) setFilesError(err.response?.data?.error || 'Could not list the run folder');
      });
    return () => { cancelled = true; };
  }, [directory, staticMode, localRun]);

  // Re-select the restored source once the candidate list knows it.
  useEffect(() => {
    if (!restoredRef.current) return;
    const wanted = restoredRef.current;
    if (files.some((file) => file.key === wanted)) {
      restoredRef.current = false;
      handleSelect(wanted);
    }
  }, [files]); // eslint-disable-line react-hooks/exhaustive-deps

  const applyInspect = (info) => {
    setInspect(info);
    const { inp, header } = info;
    setForm((current) => ({
      ...current,
      qmin: inp ? String(inp.qmin) : current.qmin,
      qmax: inp ? String(inp.qmax) : current.qmax,
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
  };

  const handleSelect = async (key) => {
    setSelectedKey(key);
    setPreview(null);
    setExportResult(null);
    setError(null);
    staticDataRef.current = null;
    if (!key) { setInspect(null); return; }
    try {
      if (staticMode) {
        const entry = files.find((file) => file.key === key);
        if (!entry?.sourceFile) throw new Error('File is no longer available — re-select the run folder');
        const text = await entry.sourceFile.text();
        let inp = null;
        let dataText = text;
        let dataName = entry.name;
        if (isInpCandidate(entry.name)) {
          inp = readStogInp(text);
          const dataEntry = (localRun?.files || []).find((file) => file.name === inp.dataFile);
          if (!dataEntry?.sourceFile) {
            throw new Error(`stog input references '${inp.dataFile}' — include it in the selected folder`);
          }
          dataText = await dataEntry.sourceFile.text();
          dataName = inp.dataFile;
        }
        const columns = readStogXy(dataText);
        const header = readDatHeader(dataText);
        staticDataRef.current = {
          q: columns[0],
          sq: columns[1],
          sigma: columns.length >= 3 ? columns[2] : null,
          inp,
          header,
          name: dataName,
        };
        applyInspect({ kind: inp ? 'inp' : 'data', inp, header, dataFile: dataName });
      } else {
        const response = await axios.post(`${API_BASE_URL}/api/scaling/preview`, { path: key, inspect: true });
        applyInspect(response.data);
      }
    } catch (err) {
      setInspect(null);
      setError(err.response?.data?.error || err.message || 'Could not inspect the selected file');
    }
  };

  const buildApiParams = () => ({
    path: selectedKey,
    qmin: numberOr(form.qmin),
    qmax: numberOr(form.qmax),
    rho0: numberOr(form.rho0),
    massDensity: numberOr(form.massDensity),
    bAvgSq: numberOr(form.bAvgSq),
    bSqAvg: numberOr(form.bSqAvg),
    formula: form.formula.trim() || undefined,
    rCutoff: numberOr(form.rCutoff),
    r0: numberOr(form.r0),
    rFitMin: numberOr(form.rFitMin),
    rFitMax: numberOr(form.rFitMax),
    rmax: numberOr(form.rmax),
    nr: numberOr(form.nr),
    lorch: form.lorch,
    despike: form.despike,
    robust: form.robust,
    lowQCorrection: form.lowQCorrection,
    useSigma: form.useSigma,
    c1Mode: form.c1Mode,
    amplitude: form.amplitude,
    enforce: form.enforce,
    enforceCutoff: form.enforce ? numberOr(form.enforceCutoff) : undefined,
  });

  const runStatic = (mode) => new Promise((resolve, reject) => {
    const data = staticDataRef.current;
    if (!data) { reject(new Error('Select a source file first')); return; }
    let config;
    let enforcement;
    try {
      config = resolveStaticConfig(form, data.inp, data.header);
      enforcement = resolveStaticEnforcement(form, data.inp);
    } catch (resolveError) { reject(resolveError); return; }
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
      if (!event.data.ok) { reject(new Error(event.data.error)); return; }
      const raw = event.data.result;
      const arr = (buffer) => (buffer ? Array.from(new Float64Array(buffer)) : null);
      resolve({
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
        config,
        guides: {
          asymptote: 1.0,
          gkLowR: -config.bAvgSq,
          drSlope: -4 * Math.PI * config.rho0 * config.bAvgSq,
          s0Target: config.bSqAvg == null ? null : 1 - config.bSqAvg / config.bAvgSq,
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
    };
    worker.addEventListener('message', onMessage);
    const q = Float64Array.from(data.q);
    const sq = Float64Array.from(data.sq);
    const sigma = form.useSigma && data.sigma ? Float64Array.from(data.sigma) : null;
    worker.postMessage(
      {
        id,
        config,
        mode,
        a: mode === 'manual' ? numberOr(form.manualA) ?? (data.inp ? data.inp.a : undefined) : undefined,
        b: mode === 'manual' ? numberOr(form.manualB) ?? (data.inp ? data.inp.b : 0) : undefined,
        q: q.buffer,
        sq: sq.buffer,
        sigma: sigma ? sigma.buffer : null,
        enforcement,
      },
      sigma ? [q.buffer, sq.buffer, sigma.buffer] : [q.buffer, sq.buffer]
    );
  });

  const runScaling = async (mode) => {
    if (!selectedKey) return;
    setRunning(true);
    setError(null);
    setExportResult(null);
    try {
      if (staticMode) {
        if (mode === 'manual' && numberOr(form.manualA) === undefined && !staticDataRef.current?.inp) {
          throw new Error('fixed scaling needs a (and optionally b)');
        }
        setPreview(await runStatic(mode));
      } else {
        const params = buildApiParams();
        if (mode === 'manual') {
          params.mode = 'manual';
          params.a = numberOr(form.manualA);
          params.b = numberOr(form.manualB);
        }
        const response = await axios.post(`${API_BASE_URL}/api/scaling/preview`, params);
        setPreview(response.data);
      }
    } catch (err) {
      setPreview(null);
      setError(err.response?.data?.error || err.message || 'Scaling failed');
    } finally {
      setRunning(false);
    }
  };

  const writeFiles = async () => {
    if (!selectedKey || !preview) return;
    setExporting(true);
    setError(null);
    try {
      if (staticMode) {
        const { series, result, diagnostics, config, enforcement } = preview;
        const label = `rmc-autoscale (browser): a=${result.a.toPrecision(8)} b=${result.b.toPrecision(8)}`;
        const stem = sanitizeFilename(
          exportForm.outStem.trim() || staticDataRef.current?.name?.replace(/\.[^.]+$/, '') || 'autoscale'
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
            source: staticDataRef.current?.name,
            a: result.a,
            b: result.b,
            mode: preview.mode,
            enforcement,
            config,
            diagnostics,
          }, null, 2)],
        ].map(([name, text]) => ({ name, data: encoder.encode(text) }));
        downloadBlob(buildZip(entries), `${stem}_autoscale.zip`);
        setExportResult({ zip: `${stem}_autoscale.zip`, count: entries.length });
      } else {
        const params = { ...buildApiParams(), mode: preview.mode };
        if (preview.mode === 'manual') {
          params.a = preview.result.a;
          params.b = preview.result.b;
        }
        if (exportForm.outDir.trim()) params.outDir = exportForm.outDir.trim();
        if (exportForm.outStem.trim()) params.outStem = exportForm.outStem.trim();
        params.force = exportForm.force;
        const response = await axios.post(`${API_BASE_URL}/api/scaling/run`, params);
        setExportResult(response.data);
      }
    } catch (err) {
      const message = err.response?.data?.error || err.message || 'Could not write the output files';
      if (err.response?.status === 409) setExportResult({ conflict: true, message });
      else setError(message);
    } finally {
      setExporting(false);
    }
  };

  const setField = (key) => (event) => {
    const value = event.target.type === 'checkbox' ? event.target.checked : event.target.value;
    setForm((current) => ({ ...current, [key]: value }));
  };

  const fzMissing = form.amplitude === 'fz'
    && numberOr(form.bSqAvg) === undefined
    && !form.formula.trim();

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

  const sourceStatus = staticMode && !localRun
    ? 'Open a run folder (header) to list scaling sources.'
    : null;

  return (
    <div className="autostog-page">
      <div className="autostog-controls">
        <div className="autostog-cluster autostog-cluster--source">
          <span className="autostog-label">SOURCE</span>
          <select
            value={selectedKey}
            onChange={(event) => handleSelect(event.target.value)}
            aria-label="Scaling source file"
          >
            <option value="">Select a file…</option>
            {files.map((file) => (
              <option key={file.key} value={file.key}>
                {isInpCandidate(file.name) ? `⚙ ${file.name}` : file.name}
              </option>
            ))}
          </select>
          {inspect?.kind === 'inp' && (
            <span className="autostog-chip" title="Hand scaling recorded in the stog input">
              hand a = {fmt(inspect.inp.a, 4)}, b = {fmt(inspect.inp.b, 4)}
            </span>
          )}
          {sourceStatus && <span className="autostog-chip autostog-chip--muted">{sourceStatus}</span>}
          {filesError && <span className="autostog-chip autostog-chip--danger">{filesError}</span>}
        </div>

        <div className="autostog-cluster">
          <label
            className="autostog-field autostog-field--formula"
            title="Chemical composition (neutron Sears table drives ⟨b⟩², ⟨b²⟩, and the S(0) target; for x-ray data set ⟨b²⟩ in Advanced instead)"
          >
            <span>Composition</span>
            <input value={form.formula} onChange={setField('formula')} placeholder="e.g. Mn3Sn" spellCheck="false" />
          </label>
          <label className="autostog-field">
            <span>Qmin Å⁻¹</span>
            <input value={form.qmin} onChange={setField('qmin')} inputMode="decimal" />
          </label>
          <label className="autostog-field">
            <span>Qmax Å⁻¹</span>
            <input value={form.qmax} onChange={setField('qmax')} inputMode="decimal" />
          </label>
          <label className="autostog-field" title="Number density in atoms/Å³ (from a NUMBER_DENSITY :: header when present)">
            <span>ρ₀ Å⁻³</span>
            <input value={form.rho0} onChange={setField('rho0')} inputMode="decimal" placeholder="auto" />
          </label>
          <label className="autostog-field" title="Alternative to ρ₀: mass density + composition convert via N_A (ADDIE convention)">
            <span>or ρ g/cm³</span>
            <input value={form.massDensity} onChange={setField('massDensity')} inputMode="decimal" />
          </label>
          {compositionChip && <span className="autostog-chip" title="Computed from the composition (Sears neutron lengths)">{compositionChip}</span>}
        </div>

        <div className="autostog-cluster autostog-cluster--actions">
          <button
            type="button"
            className="autostog-primary"
            disabled={!selectedKey || running || fzMissing}
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
          <div className="autostog-cluster">
            <label className="autostog-field">
              <span>Filter r-cut Å</span>
              <input value={form.rCutoff} onChange={setField('rCutoff')} inputMode="decimal" />
            </label>
            <label className="autostog-field">
              <span>r₀ approach Å</span>
              <input value={form.r0} onChange={setField('r0')} inputMode="decimal" placeholder="auto" />
            </label>
            <label className="autostog-field">
              <span>Fit win min</span>
              <input value={form.rFitMin} onChange={setField('rFitMin')} inputMode="decimal" placeholder="auto" />
            </label>
            <label className="autostog-field">
              <span>Fit win max</span>
              <input value={form.rFitMax} onChange={setField('rFitMax')} inputMode="decimal" placeholder="auto" />
            </label>
            <label className="autostog-field">
              <span>rmax Å</span>
              <input value={form.rmax} onChange={setField('rmax')} inputMode="decimal" />
            </label>
            <label className="autostog-field">
              <span>r points</span>
              <input value={form.nr} onChange={setField('nr')} inputMode="numeric" />
            </label>
          </div>
          <div className="autostog-cluster">
            <label className="autostog-field" title="Override the composition-derived ⟨b⟩² (needed for x-ray data)">
              <span>⟨b⟩² barn</span>
              <input value={form.bAvgSq} onChange={setField('bAvgSq')} inputMode="decimal" placeholder="from composition" />
            </label>
            <label className="autostog-field" title="Override the composition-derived ⟨b²⟩">
              <span>⟨b²⟩ barn</span>
              <input value={form.bSqAvg} onChange={setField('bSqAvg')} inputMode="decimal" placeholder="from composition" />
            </label>
            <label className="autostog-field autostog-field--select">
              <span>High-Q</span>
              <select value={form.c1Mode} onChange={setField('c1Mode')}>
                <option value="sweep">Level sweep</option>
                <option value="joint">Joint 2-dof</option>
              </select>
            </label>
            <label className="autostog-field autostog-field--select">
              <span>Amplitude</span>
              <select value={form.amplitude} onChange={setField('amplitude')}>
                <option value="density">Density limit</option>
                <option value="fz">Faber-Ziman Q→0</option>
              </select>
            </label>
            {fzMissing && (
              <span className="autostog-chip autostog-chip--danger">FZ needs ⟨b²⟩ or a formula</span>
            )}
          </div>
          <div className="autostog-cluster autostog-cluster--toggles">
            {[
              ['lorch', 'Lorch'],
              ['lowQCorrection', 'Low-Q corr.'],
              ['robust', 'Robust'],
              ['despike', 'Despike'],
              ['useSigma', 'σ column'],
              ['enforce', 'Enforce low-r'],
            ].map(([key, label]) => (
              <label key={key} className={`autostog-toggle${form[key] ? ' is-on' : ''}`}>
                <input type="checkbox" checked={form[key]} onChange={setField(key)} />
                {label}
              </label>
            ))}
            {form.enforce && (
              <label className="autostog-field">
                <span>Cutoff Å</span>
                <input value={form.enforceCutoff} onChange={setField('enforceCutoff')} inputMode="decimal" />
              </label>
            )}
          </div>
          <div className="autostog-cluster autostog-cluster--manual">
            <span className="autostog-label">FIXED SCALING</span>
            <label className="autostog-field">
              <span>a</span>
              <input value={form.manualA} onChange={setField('manualA')} inputMode="decimal" />
            </label>
            <label className="autostog-field">
              <span>b</span>
              <input value={form.manualB} onChange={setField('manualB')} inputMode="decimal" />
            </label>
            <button
              type="button"
              className="autostog-pill"
              disabled={!selectedKey || running}
              onClick={() => runScaling('manual')}
            >
              Run fixed (a, b)
            </button>
          </div>
        </div>
      )}

      {error && <div className="autostog-banner autostog-banner--danger">{error}</div>}

      <details className="autostog-explainer">
        <summary>How Auto StoG works</summary>
        <ol>
          <li><b>Inputs:</b> the composition and your [Qmin, Qmax] window are all that is
            required. ρ₀ comes from the data header, from mass density (converted via N_A),
            or your value; every other parameter has a physics-based default.</li>
          <li><b>Composition →</b> ⟨b⟩², ⟨b²⟩ (Sears neutron lengths) and the Keen Eq. 21
            target S(0) = 1 − ⟨b²⟩/⟨b⟩², which also anchors the analytic correction for the
            unmeasured [0, Qmin] range.</li>
          <li><b>Level sweep:</b> every high-Q window is tested for statistical flatness;
            the measured level pins the offset (b = 1 − a·level).</li>
          <li><b>Amplitude:</b> the low-r density limit (g → 0 below the first shell,
            self-consistent with the Fourier filter), with the Q→0 Faber-Ziman limit as an
            independent criterion — their concordance is the absolute-scale trust metric.</li>
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
            Pick a classic <code>stog.inp</code> or a rebinned S(Q) file, enter the
            composition and Q window, then hit <b>Auto-scale</b>. The engine determines the
            scale and offset from the physics — the statistically flat high-Q level, the
            low-r density limit, and the composition's Faber-Ziman constraints — replacing
            the classic stog “try again” loop, and reports the honest fit quality before
            any enforcement.
            {staticMode && ' Everything runs in your browser; files never leave your device.'}
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
          {staticMode ? (
            <>
              <label className="autostog-field autostog-field--wide">
                <span>Name stem</span>
                <input
                  value={exportForm.outStem}
                  onChange={(event) => setExportForm((current) => ({ ...current, outStem: event.target.value }))}
                  placeholder="default: data file stem"
                  spellCheck="false"
                />
              </label>
              <button type="button" className="autostog-primary" disabled={exporting} onClick={writeFiles}>
                {exporting ? 'Packing…' : 'Download .zip (RMCProfile files)'}
              </button>
              {exportResult?.zip && (
                <span className="autostog-chip autostog-chip--good">
                  {exportResult.zip} · {exportResult.count} files
                </span>
              )}
            </>
          ) : (
            <>
              <label className="autostog-field autostog-field--wide">
                <span>Output folder</span>
                <input
                  value={exportForm.outDir}
                  onChange={(event) => setExportForm((current) => ({ ...current, outDir: event.target.value }))}
                  placeholder="default: autoscale/ beside the input"
                  spellCheck="false"
                />
              </label>
              <label className="autostog-field">
                <span>Name stem</span>
                <input
                  value={exportForm.outStem}
                  onChange={(event) => setExportForm((current) => ({ ...current, outStem: event.target.value }))}
                  placeholder="classic names"
                  spellCheck="false"
                />
              </label>
              <label className={`autostog-toggle${exportForm.force ? ' is-on' : ''}`}>
                <input
                  type="checkbox"
                  checked={exportForm.force}
                  onChange={(event) => setExportForm((current) => ({ ...current, force: event.target.checked }))}
                />
                Force overwrite
              </label>
              <button type="button" className="autostog-primary" disabled={exporting} onClick={writeFiles}>
                {exporting ? 'Writing…' : 'Write RMCProfile files'}
              </button>
              {exportResult?.conflict && (
                <span className="autostog-chip autostog-chip--warn">
                  Files exist — enable Force overwrite or change the folder/stem.
                </span>
              )}
              {exportResult?.outputs && (
                <span className="autostog-chip autostog-chip--good" title={Object.values(exportResult.outputs).join('\n')}>
                  {Object.keys(exportResult.outputs).length} files → {exportResult.outDir}
                </span>
              )}
            </>
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
