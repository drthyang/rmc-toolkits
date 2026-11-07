// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The termination-ripple assessment of a browser-loaded run (a picked folder,
// the Demo, the Copilot's in-process run): the files beside the run's fit
// files read as rmc_toolkits.ripples.ripples_for_run reads them on the
// server, and assessed by the port of its engine (ripples.js). Flask mode gets
// the same result from /api/ripples. ripplesParity.test.js holds this module
// to the Python golden.
//
// Per fit file (<stem>_FT_XFQ<n>.csv, _Xray_PDF<n>.csv, _PDF<n>.csv,
// _FQ<n>.csv, _SQ<n>.csv): <stem>.dat (ripple_settings), <stem>_PDFpartials.csv
// (the closest approach; a pair's g(r)), for the transform of an X-ray F(Q)
// <stem>_XFQ<n>.csv else <stem>_FQ<n>.csv, for a real-space block the input
// file it names and the Auto StoG provenance (*_provenance.json, in the
// folder or its autoscale/) or classic stog input (*.inp, stog_input.dat)
// that made it. One result is kept per set of file versions and arguments.

import { readRmcCsv, rwpColumns } from './browserData';
import {
    assessDataset, assessRun, autoStogTitle, BLOCKS, closestApproach, datasetFile, provenanceTransform, rippleSettings,
    stogTransform
} from './ripples';
import { readStogInp, readStogXy } from './workers/autoScale';

const dirname = (path) => (path.includes('/') ? path.slice(0, path.lastIndexOf('/')) : '');
const joinPath = (dir, name) => (dir ? `${dir}/${name}` : name);
const versionOf = (file) => (file ? `${file.path}:${file.modified ?? ''}:${file.size ?? ''}` : 'absent');
const byCodePoint = (a, b) => (a < b ? -1 : a > b ? 1 : 0);

// `name` (a .dat FILENAME) against the fit file's folder, `.`/`..` resolved;
// null for an absolute path (never inside a picked folder).
const resolveRelative = (dir, name) => {
    const clean = String(name).replace(/\\/g, '/');
    if (clean.startsWith('/') || /^[A-Za-z]:\//.test(clean)) return null;
    const parts = dir ? dir.split('/') : [];
    for (const part of clean.split('/')) {
        if (!part || part === '.') continue;
        if (part === '..') {
            if (!parts.length) return null;
            parts.pop();
        } else parts.push(part);
    }
    return parts.join('/');
};

/** A fit CSV's x, model and experiment columns (rwpColumns' roles). */
const readColumns = async (file) => {
    const series = readRmcCsv(await file.sourceFile.text(), file.name);
    const roles = rwpColumns(series.labels, series.data.length);
    if (!roles) throw new Error(`${file.name} needs an experimental column beside the calculated one`);
    return [series.data[0], series.data[roles[0]], series.data[roles[1]]];
};

/** The data file's own transform and processing, from the files beside it (data_processing). */
const dataProcessing = async (byPath, runFiles, dir, dataName, inputText) => {
    let transform = null;
    for (const directory of [dir, joinPath(dir, 'autoscale')]) {
        if (transform !== null) continue;
        const inside = (runFiles || []).filter((file) => dirname(file.path || file.name) === directory);
        const provenance = inside.filter((file) => /_provenance\.json$/.test(file.name))
            .sort((a, b) => byCodePoint(a.path || a.name, b.path || b.name));
        for (const file of provenance) {
            let found = null;
            try {
                found = provenanceTransform(JSON.parse(await file.sourceFile.text()), dataName);
            } catch {
                found = null;
            }
            if (found !== null) { transform = { ...found, file: file.name }; break; }
        }
        if (transform !== null) continue;
        const inputs = inside.filter((file) => /\.inp$/.test(file.name) || file.name === 'stog_input.dat')
            .sort((a, b) => byCodePoint(a.path || a.name, b.path || b.name));
        for (const file of inputs) {
            let found = null;
            try {
                found = stogTransform(readStogInp(await file.sourceFile.text()), dataName);
            } catch {
                found = null;
            }
            if (found !== null) { transform = { ...found, file: file.name }; break; }
        }
    }
    const auto = Boolean(inputText) && autoStogTitle(inputText);
    let processing = null;
    if (transform !== null || auto) {
        processing = {
            autoStog: Boolean(auto || (transform !== null && transform.source === 'provenance')),
            filterCutoff: transform ? transform.filterCutoff ?? null : null,
            enforceCutoff: transform ? transform.enforceCutoff ?? null : null,
            source: transform ? transform.file ?? null : null
        };
    }
    return [transform, processing];
};

const failed = (file, info, error) => ({
    status: 'failed', reason: `could not read the run's files: ${error?.message || error}`, file: file.name, kind: info.kind,
    number: info.number, total: true, block: null, fit: null, transform: null, resolution: null, peaks: [], lowR: null,
    processing: null, space: null, rFirst: null, _analysis: null
});

// Results by the versions of the files they were computed from and the arguments.
const MAX_CACHED = 8;
const cache = new Map();
const remember = (key, value) => {
    cache.delete(key);
    cache.set(key, value);
    while (cache.size > MAX_CACHED) cache.delete(cache.keys().next().value);
    return value;
};

/** Forget every kept result (tests). */
export const clearRipplesCache = () => cache.clear();

/**
 * assessRun for a run's fit files (`fitFiles`, each { name, path, sourceFile,
 * modified, size }) from `runFiles`, every file of the run; `args` the
 * caller's { features, rRange, split, pair, qmax, qmin, window }. A file that
 * cannot be read gives a failed dataset.
 */
export const ripplesForRun = async (fitFiles, runFiles, args = {}) => {
    const {
        features = null, rRange = null, split = null, pair = null, qmax = null, qmin = null, window = null
    } = args;
    const byPath = new Map((runFiles || []).map((file) => [file.path || file.name, file]));
    const override = qmax !== null && qmax !== undefined ? { qMax: qmax, qMin: qmin, window } : null;
    const relevant = (fitFiles || []).filter((file) => datasetFile(file.name));
    const companions = new Set();
    for (const file of relevant) {
        const info = datasetFile(file.name);
        const dir = dirname(file.path || file.name);
        for (const name of [`${info.stem}.dat`, `${info.stem}_PDFpartials.csv`, `${info.stem}_XFQ${info.number}.csv`,
            `${info.stem}_FQ${info.number}.csv`]) companions.add(joinPath(dir, name));
    }
    const versions = [...relevant.map(versionOf), ...[...companions].sort().map((path) => versionOf(byPath.get(path))),
        ...(runFiles || []).filter((file) => /_provenance\.json$|\.inp$|^stog_input\.dat$|\.(gr|sq|fq|dat)$/i.test(file.name))
            .map(versionOf).sort()];
    const key = JSON.stringify([versions, features, rRange, split, pair, qmax, qmin, window]);
    if (cache.has(key)) return cache.get(key);

    const memo = new Map();
    const once = (id, make) => {
        if (!memo.has(id)) memo.set(id, make());
        return memo.get(id);
    };
    const payloads = [];
    let pairSeries = null;
    for (const file of relevant) {
        const info = datasetFile(file.name);
        const dir = dirname(file.path || file.name);
        const find = (name) => byPath.get(joinPath(dir, name)) || null;
        try {
            const settings = await once(`dat:${dir}:${info.stem}`, async () => {
                const dat = find(`${info.stem}.dat`);
                return dat ? rippleSettings(await dat.sourceFile.text()) : null;
            });
            const partials = await once(`partials:${dir}:${info.stem}`, async () => {
                const entry = find(`${info.stem}_PDFpartials.csv`);
                return entry ? readRmcCsv(await entry.sourceFile.text(), entry.name) : null;
            });
            const rFirst = closestApproach(partials ? partials.data : null, settings ? settings.minimumDistances : null);
            if (pair !== null && partials !== null && pairSeries === null) {
                const wanted = String(pair).split('-').map((part) => part.trim()).sort();
                partials.labels.slice(1).forEach((label, k) => {
                    if (pairSeries !== null) return;
                    const parts = label.split('-').map((part) => part.trim()).sort();
                    if (parts.length === wanted.length && parts.every((part, i) => part === wanted[i])) {
                        pairSeries = { label: label.trim(), r: partials.data[0], g: partials.data[k + 1] };
                    }
                });
            }
            const [x, calc, data] = await readColumns(file);
            let fq = null;
            let inputColumns = null;
            let dataTransform = null;
            let processing = null;
            if (info.kind === 'xpdf') {
                const fqFile = find(`${info.stem}_XFQ${info.number}.csv`) || find(`${info.stem}_FQ${info.number}.csv`);
                if (fqFile) fq = await readColumns(fqFile);
            } else if (info.kind === 'npdf' || info.kind === 'xray-real') {
                const blocks = (settings?.datasets || []).filter((entry) => entry.block === BLOCKS[info.kind]);
                const dataset = info.number > 0 && info.number <= blocks.length ? blocks[info.number - 1] : null;
                if (dataset?.file) {
                    const relative = dataset.file.replace(/\\/g, '/');
                    const inputPath = resolveRelative(dir, relative);
                    const inputFile = inputPath === null ? null : byPath.get(inputPath) || null;
                    let text = null;
                    if (inputFile) {
                        text = await inputFile.sourceFile.text();
                        inputColumns = readStogXy(text);
                    }
                    [dataTransform, processing] = await dataProcessing(byPath, runFiles, dir, relative.split('/').pop(), text);
                }
            }
            payloads.push(assessDataset({
                name: file.name, kind: info.kind, number: info.number, settings, x, calc, data, fq, inputColumns, dataTransform,
                override, rFirst, processing
            }));
        } catch (error) {
            payloads.push(failed(file, info, error));
        }
    }
    const result = assessRun(payloads, { features, rRange, split, pair: pair !== null ? pairSeries : null });
    return remember(key, result);
};
