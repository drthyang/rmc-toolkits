// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

import React, { useEffect, useMemo, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from '../api';
import InteractivePlot from './InteractivePlot';
import { isStaticMode } from '../browserData';
import './AutoStogPage.css';

// Files that can seed a scaling session: classic stog inputs first, then data.
const isInpCandidate = (name) => name.endsWith('.inp') || name === 'stog_input.dat';
const isDataCandidate = (name) => /\.(sq|fq|dat)$/i.test(name) && !/^(scale|ft\.dat|.*_rmc)/i.test(name);

const numberOr = (value) => {
  if (value === '' || value === null || value === undefined) return undefined;
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : undefined;
};

const fmt = (value, digits = 4) => (
  Number.isFinite(value) ? Number(value).toPrecision(digits) : '—'
);

const EMPTY_FORM = {
  qmin: '', qmax: '', rho0: '', bAvgSq: '', bSqAvg: '', formula: '',
  rCutoff: '', r0: '', rFitMin: '', rFitMax: '', rmax: '', nr: '',
  c1Mode: 'sweep', amplitude: 'density',
  lorch: false, despike: false, robust: true, lowQCorrection: true, useSigma: true,
  enforce: true, enforceCutoff: '', manualA: '', manualB: '',
};

const AutoStogPage = ({ directory }) => {
  const [files, setFiles] = useState([]);
  const [filesError, setFilesError] = useState(null);
  const [selectedPath, setSelectedPath] = useState('');
  const [inspect, setInspect] = useState(null);
  const [form, setForm] = useState(EMPTY_FORM);
  const [preview, setPreview] = useState(null);
  const [running, setRunning] = useState(false);
  const [error, setError] = useState(null);
  const [exportForm, setExportForm] = useState({ outDir: '', outStem: '', force: false });
  const [exportResult, setExportResult] = useState(null);
  const [exporting, setExporting] = useState(false);

  const staticMode = isStaticMode();

  // Candidate files in the selected run folder (header controls the folder).
  useEffect(() => {
    if (staticMode) return;
    let cancelled = false;
    setFilesError(null);
    axios.get(`${API_BASE_URL}/api/files`, { params: { dir: directory || '.' } })
      .then((response) => {
        if (cancelled) return;
        const candidates = (response.data.files || [])
          .filter((item) => item.type === 'file')
          .filter((item) => isInpCandidate(item.name) || isDataCandidate(item.name))
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
  }, [directory, staticMode]);

  const applyInspect = (info) => {
    setInspect(info);
    const { inp, header } = info;
    setForm((current) => ({
      ...current,
      qmin: inp ? String(inp.qmin) : current.qmin,
      qmax: inp ? String(inp.qmax) : current.qmax,
      rho0: inp ? String(inp.rho0) : (header?.numberDensity != null ? String(header.numberDensity) : ''),
      bAvgSq: inp ? String(inp.bAvgSq) : current.bAvgSq,
      rCutoff: inp ? String(inp.rCutoff) : current.rCutoff,
      r0: header?.minDistance != null ? String(header.minDistance) : '',
      rmax: inp ? String(inp.rmax) : current.rmax,
      nr: inp ? String(inp.nr) : current.nr,
      lorch: inp ? Boolean(inp.lorch) : current.lorch,
      enforce: Boolean(inp),
      enforceCutoff: inp ? String(inp.peakCutoff) : '',
      manualA: inp ? String(inp.a) : '',
      manualB: inp ? String(inp.b) : '',
    }));
  };

  const handleSelect = async (path) => {
    setSelectedPath(path);
    setPreview(null);
    setExportResult(null);
    setError(null);
    if (!path) { setInspect(null); return; }
    try {
      const response = await axios.post(`${API_BASE_URL}/api/scaling/preview`, { path, inspect: true });
      applyInspect(response.data);
    } catch (err) {
      setInspect(null);
      setError(err.response?.data?.error || 'Could not inspect the selected file');
    }
  };

  const buildParams = () => ({
    path: selectedPath,
    qmin: numberOr(form.qmin),
    qmax: numberOr(form.qmax),
    rho0: numberOr(form.rho0),
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

  const runScaling = async (mode) => {
    if (!selectedPath) return;
    setRunning(true);
    setError(null);
    setExportResult(null);
    try {
      const params = buildParams();
      if (mode === 'manual') {
        params.mode = 'manual';
        params.a = numberOr(form.manualA);
        params.b = numberOr(form.manualB);
      }
      const response = await axios.post(`${API_BASE_URL}/api/scaling/preview`, params);
      setPreview(response.data);
    } catch (err) {
      setPreview(null);
      setError(err.response?.data?.error || 'Scaling failed');
    } finally {
      setRunning(false);
    }
  };

  const writeFiles = async (force = exportForm.force) => {
    if (!selectedPath || !preview) return;
    setExporting(true);
    setError(null);
    try {
      const params = { ...buildParams(), mode: preview.mode };
      if (preview.mode === 'manual') {
        params.a = preview.result.a;
        params.b = preview.result.b;
      }
      if (exportForm.outDir.trim()) params.outDir = exportForm.outDir.trim();
      if (exportForm.outStem.trim()) params.outStem = exportForm.outStem.trim();
      params.force = force;
      const response = await axios.post(`${API_BASE_URL}/api/scaling/run`, params);
      setExportResult(response.data);
    } catch (err) {
      const message = err.response?.data?.error || 'Could not write the output files';
      if (err.response?.status === 409) {
        setExportResult({ conflict: true, message });
      } else {
        setError(message);
      }
    } finally {
      setExporting(false);
    }
  };

  const setField = (key) => (event) => {
    const value = event.target.type === 'checkbox' ? event.target.checked : event.target.value;
    setForm((current) => ({ ...current, [key]: value }));
  };

  const sqPlot = useMemo(() => {
    if (!preview) return null;
    const { series, guides } = preview;
    const plots = [
      { label: 'Measured S(Q)', x: series.q, y: series.sqRaw },
      { label: 'Auto-scaled a·S + b', x: series.q, y: series.sqScaled },
      { label: 'Filtered S(Q)', x: series.q, y: series.sqFiltered },
      {
        label: 'S → 1 asymptote',
        x: [series.q[0], series.q[series.q.length - 1]],
        y: [guides.asymptote, guides.asymptote],
      },
    ];
    if (guides.level != null && guides.levelWindow) {
      plots.push({
        label: `Measured level ${fmt(guides.level)}`,
        x: guides.levelWindow,
        y: [guides.level, guides.level],
      });
    }
    return {
      title: 'S(Q): measured → scaled → filtered',
      xLabel: 'Q (Å^{-1})',
      yLabel: 'S(Q)',
      series: plots,
    };
  }, [preview]);

  const RMAX_DISPLAY = 8;

  const slicedR = useMemo(() => {
    if (!preview) return null;
    const { r } = preview.series;
    let end = r.length;
    for (let index = 0; index < r.length; index += 1) {
      if (r[index] > RMAX_DISPLAY) { end = index; break; }
    }
    return { end, r: r.slice(0, end) };
  }, [preview]);

  const gkPlot = useMemo(() => {
    if (!preview || !slicedR) return null;
    const { series, guides } = preview;
    const plots = [
      { label: 'G_K(r) filtered', x: slicedR.r, y: series.gk.slice(0, slicedR.end) },
      {
        label: 'Theory −⟨b⟩²',
        x: [0, RMAX_DISPLAY],
        y: [guides.gkLowR, guides.gkLowR],
      },
    ];
    if (series.gkEnforced) {
      plots.push({
        label: 'Enforced (written to RMC files)',
        x: slicedR.r,
        y: series.gkEnforced.slice(0, slicedR.end),
      });
    }
    return {
      title: `G_K(r) low-r region (fit window ${fmt(guides.rFitWindow?.[0], 3)}–${fmt(guides.rFitWindow?.[1], 3)} Å)`,
      xLabel: 'r (Å)',
      yLabel: 'G_K(r)',
      series: plots,
    };
  }, [preview, slicedR]);

  const drPlot = useMemo(() => {
    if (!preview || !slicedR) return null;
    const { series, guides } = preview;
    const plots = [
      { label: 'D(r)', x: slicedR.r, y: series.dr.slice(0, slicedR.end) },
      {
        label: 'Theory −4πρ₀⟨b⟩²r',
        x: [0, RMAX_DISPLAY],
        y: [0, guides.drSlope * RMAX_DISPLAY],
      },
    ];
    if (series.drEnforced) {
      plots.push({
        label: 'Enforced (written to RMC files)',
        x: slicedR.r,
        y: series.drEnforced.slice(0, slicedR.end),
      });
    }
    return {
      title: 'D(r) low-r region',
      xLabel: 'r (Å)',
      yLabel: 'D(r)',
      series: plots,
    };
  }, [preview, slicedR]);

  if (staticMode) {
    return (
      <div className="autostog-page">
        <div className="autostog-static-note">
          <h2>Auto StoG needs the local backend</h2>
          <p>
            Automatic S(Q) scaling runs server-side. Start the Flask app
            (<code>python web_app/backend/app.py</code>) and open this page from there —
            the hosted static dashboard cannot run it yet.
          </p>
        </div>
      </div>
    );
  }

  const diagnostics = preview?.diagnostics;
  const reference = preview?.inp;

  return (
    <div className="autostog-page">
      <aside className="autostog-sidebar">
        <section className="autostog-card">
          <h3>1 · Source</h3>
          <p className="autostog-hint">
            Pick a classic <code>stog.inp</code> (everything pre-fills) or a rebinned
            S(Q) data file from the run folder selected in the header.
          </p>
          <select
            value={selectedPath}
            onChange={(event) => handleSelect(event.target.value)}
            aria-label="Scaling source file"
          >
            <option value="">Select a file…</option>
            {files.map((file) => (
              <option key={file.path} value={file.path}>
                {isInpCandidate(file.name) ? `⚙ ${file.name}` : file.name}
              </option>
            ))}
          </select>
          {filesError && <div className="autostog-error">{filesError}</div>}
          {inspect && (
            <div className="autostog-inspect">
              <span>{inspect.kind === 'inp' ? 'Classic stog input' : 'Data file'}</span>
              {inspect.kind === 'inp' && (
                <span>hand a = {fmt(inspect.inp.a)}, b = {fmt(inspect.inp.b)}</span>
              )}
              {inspect.header?.title && <span>{inspect.header.title}</span>}
            </div>
          )}
        </section>

        <section className="autostog-card">
          <h3>2 · Parameters</h3>
          <div className="autostog-grid">
            <label>Qmin (Å⁻¹)<input value={form.qmin} onChange={setField('qmin')} inputMode="decimal" /></label>
            <label>Qmax (Å⁻¹)<input value={form.qmax} onChange={setField('qmax')} inputMode="decimal" /></label>
            <label>ρ₀ (Å⁻³)<input value={form.rho0} onChange={setField('rho0')} inputMode="decimal" /></label>
            <label>⟨b⟩² (barn)<input value={form.bAvgSq} onChange={setField('bAvgSq')} inputMode="decimal" /></label>
            <label>Formula<input value={form.formula} onChange={setField('formula')} placeholder="e.g. SrTiO3" spellCheck="false" /></label>
            <label>⟨b²⟩ (barn)<input value={form.bSqAvg} onChange={setField('bSqAvg')} inputMode="decimal" placeholder="from formula" /></label>
          </div>
          <details className="autostog-advanced">
            <summary>Advanced</summary>
            <div className="autostog-grid">
              <label>Filter r-cutoff<input value={form.rCutoff} onChange={setField('rCutoff')} inputMode="decimal" /></label>
              <label>r₀ closest approach<input value={form.r0} onChange={setField('r0')} inputMode="decimal" /></label>
              <label>Fit window min<input value={form.rFitMin} onChange={setField('rFitMin')} inputMode="decimal" placeholder="auto" /></label>
              <label>Fit window max<input value={form.rFitMax} onChange={setField('rFitMax')} inputMode="decimal" placeholder="auto" /></label>
              <label>rmax (Å)<input value={form.rmax} onChange={setField('rmax')} inputMode="decimal" /></label>
              <label>r points<input value={form.nr} onChange={setField('nr')} inputMode="numeric" /></label>
              <label>High-Q architecture
                <select value={form.c1Mode} onChange={setField('c1Mode')}>
                  <option value="sweep">Level sweep (recommended)</option>
                  <option value="joint">Joint 2-dof fit</option>
                </select>
              </label>
              <label>Amplitude criterion
                <select value={form.amplitude} onChange={setField('amplitude')}>
                  <option value="density">Low-r density limit</option>
                  <option value="fz">Faber-Ziman Q→0 limit</option>
                </select>
              </label>
            </div>
            <div className="autostog-checks">
              <label><input type="checkbox" checked={form.lorch} onChange={setField('lorch')} />Lorch window</label>
              <label><input type="checkbox" checked={form.lowQCorrection} onChange={setField('lowQCorrection')} />Low-Q correction</label>
              <label><input type="checkbox" checked={form.robust} onChange={setField('robust')} />Robust (Huber)</label>
              <label><input type="checkbox" checked={form.despike} onChange={setField('despike')} />Despike glitches</label>
              <label><input type="checkbox" checked={form.useSigma} onChange={setField('useSigma')} />Use σ column</label>
              <label><input type="checkbox" checked={form.enforce} onChange={setField('enforce')} />Classic low-r enforcement</label>
              {form.enforce && (
                <label className="autostog-inline">Enforce cutoff (Å)
                  <input value={form.enforceCutoff} onChange={setField('enforceCutoff')} inputMode="decimal" />
                </label>
              )}
            </div>
            <div className="autostog-manual">
              <span>Fixed scaling (expert):</span>
              <label>a<input value={form.manualA} onChange={setField('manualA')} inputMode="decimal" /></label>
              <label>b<input value={form.manualB} onChange={setField('manualB')} inputMode="decimal" /></label>
              <button
                type="button"
                className="autostog-secondary"
                disabled={!selectedPath || running}
                onClick={() => runScaling('manual')}
              >
                Run fixed (a, b)
              </button>
            </div>
          </details>
          <button
            type="button"
            className="autostog-primary"
            disabled={!selectedPath || running}
            onClick={() => runScaling('auto')}
          >
            {running ? 'Fitting…' : 'Auto-scale'}
          </button>
        </section>

        {preview && (
          <section className="autostog-card">
            <h3>3 · Export</h3>
            <p className="autostog-hint">
              Writes the classic stog file family (incl. <code>scale_ft_rmc.*</code> for
              RMCProfile) plus a provenance JSON. Nothing is overwritten without Force.
            </p>
            <label className="autostog-block">Output folder
              <input
                value={exportForm.outDir}
                onChange={(event) => setExportForm((current) => ({ ...current, outDir: event.target.value }))}
                placeholder="default: autoscale/ beside the input"
                spellCheck="false"
              />
            </label>
            <label className="autostog-block">Name stem (optional)
              <input
                value={exportForm.outStem}
                onChange={(event) => setExportForm((current) => ({ ...current, outStem: event.target.value }))}
                placeholder="default: classic names"
                spellCheck="false"
              />
            </label>
            <div className="autostog-export-row">
              <label>
                <input
                  type="checkbox"
                  checked={exportForm.force}
                  onChange={(event) => setExportForm((current) => ({ ...current, force: event.target.checked }))}
                />
                Force overwrite
              </label>
              <button
                type="button"
                className="autostog-primary"
                disabled={exporting}
                onClick={() => writeFiles()}
              >
                {exporting ? 'Writing…' : 'Write RMCProfile files'}
              </button>
            </div>
            {exportResult?.conflict && (
              <div className="autostog-warn">
                Files already exist. Enable Force overwrite (or set a different folder/stem) and write again.
              </div>
            )}
            {exportResult?.outputs && (
              <div className="autostog-outputs">
                <b>Written to {exportResult.outDir}</b>
                <ul>
                  {Object.entries(exportResult.outputs).map(([key, path]) => (
                    <li key={key}><code>{path.split('/').pop()}</code></li>
                  ))}
                </ul>
              </div>
            )}
          </section>
        )}
      </aside>

      <div className="autostog-main">
        {error && <div className="autostog-error autostog-banner">{error}</div>}
        {!preview && !error && (
          <div className="autostog-empty">
            <h2>Automatic total-scattering scaling</h2>
            <p>
              Pick a source file, then hit <b>Auto-scale</b>. The engine finds the scale and
              offset from the physics — the high-Q asymptote (level sweep) and the low-r
              density limit — instead of the classic stog “try again” loop, and shows the
              honest fit quality before any enforcement.
            </p>
          </div>
        )}
        {preview && diagnostics && (
          <div className="autostog-readout">
            <div className="autostog-stat">
              <span className="autostog-stat-label">Correction</span>
              <span className="autostog-stat-value">
                a = {fmt(preview.result.a, 5)} · b = {fmt(preview.result.b, 5)}
              </span>
              {reference && (
                <span className="autostog-stat-sub">
                  stog.inp hand: a = {fmt(reference.a, 5)}, b = {fmt(reference.b, 5)}
                </span>
              )}
            </div>
            <div className="autostog-stat">
              <span className="autostog-stat-label">{preview.mode === 'manual' ? 'Mode' : 'Convergence'}</span>
              <span className="autostog-stat-value">
                {preview.mode === 'manual'
                  ? 'fixed (a, b)'
                  : `${preview.result.converged ? '✓' : '✗'} ${preview.result.iterations} iterations`}
              </span>
              {diagnostics.level != null && (
                <span className="autostog-stat-sub">
                  level {fmt(diagnostics.level, 5)}
                  {Number.isFinite(diagnostics.level_uncertainty) ? ` ± ${fmt(diagnostics.level_uncertainty, 2)}` : ''}
                </span>
              )}
            </div>
            <div className="autostog-stat">
              <span className="autostog-stat-label">Fit quality</span>
              <span className="autostog-stat-value">low-r rms {fmt(diagnostics.low_r_rms_pre_enforcement, 3)}</span>
              <span className="autostog-stat-sub">C1 tail mean {fmt(diagnostics.c1_tail_mean, 5)}</span>
            </div>
            <div className={`autostog-stat ${diagnostics.density_limit_satisfied ? 'is-good' : 'is-bad'}`}>
              <span className="autostog-stat-label">Density limit</span>
              <span className="autostog-stat-value">
                {diagnostics.density_limit_satisfied ? 'satisfied' : 'NOT satisfiable'}
              </span>
              <span className="autostog-stat-sub">
                {diagnostics.density_limit_satisfied
                  ? 'necessary, not sufficient — validate the absolute scale'
                  : 'absolute scale unrecoverable from this data alone'}
              </span>
            </div>
            {diagnostics.amplitude_concordance != null && (
              <div className={`autostog-stat ${diagnostics.amplitudes_concordant ? 'is-good' : 'is-warn'}`}>
                <span className="autostog-stat-label">Amplitude concordance</span>
                <span className="autostog-stat-value">a_fz / a = {fmt(diagnostics.amplitude_concordance, 3)}</span>
                <span className="autostog-stat-sub">
                  {diagnostics.amplitudes_concordant ? 'independent criteria agree' : 'criteria disagree — check ρ₀ / low-Q'}
                </span>
              </div>
            )}
          </div>
        )}
        {sqPlot && (
          <div className="autostog-plot">
            <InteractivePlot file={{ path: 'autostog-sq', name: 'autostog-sq' }} variant="wide" plotData={sqPlot} />
          </div>
        )}
        {gkPlot && (
          <div className="autostog-plot">
            <InteractivePlot file={{ path: 'autostog-gk', name: 'autostog-gk' }} variant="wide" plotData={gkPlot} />
          </div>
        )}
        {drPlot && (
          <div className="autostog-plot">
            <InteractivePlot file={{ path: 'autostog-dr', name: 'autostog-dr' }} variant="wide" plotData={drPlot} />
          </div>
        )}
      </div>
    </div>
  );
};

export default AutoStogPage;
