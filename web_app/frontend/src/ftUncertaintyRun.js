// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The G(r) error ribbon of a browser-loaded run (a picked folder, the Demo,
// the Copilot's in-process run): the files beside <stem>_FT_XFQ<n>.csv read
// as ft_uncertainty.ft_uncertainty_for_run reads them on the server, and
// assessed by the port of its engine (ftUncertainty.js). Flask mode gets the
// same payload from /api/plot/data (`uncertainty` on an xPDF plot), so a
// server folder and the same files picked in the browser draw the same
// ribbon; ftUncertaintyParity.test.js holds this module to the Python golden.
//
// The run's files: `<stem>.dat` (the n-th X-ray dataset's FILENAME and
// CONVOLVE), the input file it names (resolved against the transform's folder),
// `<stem>_XFQ<n>.csv`, else `<stem>_FQ<n>.csv` (the run's F(Q)_Expt), and,
// for the box radius, `<stem>_PDFpartials.csv`. One payload per set of file
// versions is kept, so a Live Data round that changed none of them returns the
// same object (and a chart keeps its zoom).

import { readRmcCsv, rwpColumns } from './browserData';
import { assessFtUncertainty, FT_FILE_RE, RMC_FT_TOLERANCE, runDatasets, xrayDataset } from './ftUncertainty';
import { readStogXy } from './workers/autoScale';

/** The payload of a computation that threw: no ribbon, the message as the reason. */
export const ftFailure = (message) => ({
    status: 'failed',
    reason: message,
    sigma: null,
    transform: null,
    reproduction: { tolerance: RMC_FT_TOLERANCE, gr: null, fq: null },
    convolutionCheck: null,
    correlation: null,
    sources: {}
});

const dirname = (path) => (path.includes('/') ? path.slice(0, path.lastIndexOf('/')) : '');
const joinPath = (dir, name) => (dir ? `${dir}/${name}` : name);

// `name` (a .dat FILENAME) against the transform's folder, `.`/`..` resolved;
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

const versionOf = (file) => (file ? `${file.path}:${file.modified ?? ''}:${file.size ?? ''}` : 'absent');

const experimentalColumn = ({ labels, data }) => {
    const roles = rwpColumns(labels, data.length);
    if (!roles) throw new Error('needs an experimental column beside the calculated one');
    return data[roles[1]];
};

const largestFinite = (values) => {
    let best = -Infinity;
    for (const value of values) if (Number.isFinite(value) && value > best) best = value;
    return Number.isFinite(best) ? best : null;
};

// Payloads by the versions of the files they were computed from.
const MAX_CACHED = 8;
const cache = new Map();
const remember = (key, payload) => {
    cache.delete(key);
    cache.set(key, payload);
    while (cache.size > MAX_CACHED) cache.delete(cache.keys().next().value);
    return payload;
};

/** Forget every kept payload (tests). */
export const clearFtUncertaintyCache = () => cache.clear();

/**
 * The `uncertainty` payload of the transform file `ftFile` ({ name, path,
 * sourceFile, modified, size }), from `runFiles`, every file of its run (each
 * { name, path, sourceFile, modified, size, plotData? }). Rejects only when a
 * file cannot be read before the versions are known; a later failure is the
 * payload of ftFailure.
 */
export const ftUncertaintyForRun = async (ftFile, runFiles) => {
    const match = FT_FILE_RE.exec(ftFile.name);
    if (!match) throw new Error(`${ftFile.name} is not an RMCProfile transform file (<stem>_FT_XFQ<n>.csv)`);
    const [, stem, numberText] = match;
    const number = Number(numberText);
    const dir = dirname(ftFile.path || ftFile.name);
    const byPath = new Map((runFiles || []).map((file) => [file.path || file.name, file]));
    const find = (name) => byPath.get(joinPath(dir, name)) || null;
    const candidates = [`${stem}_XFQ${number}.csv`, `${stem}_FQ${number}.csv`];
    const names = {
        transform: ftFile.name,
        settings: `${stem}.dat`,
        fqCandidates: `${candidates[0]} or ${candidates[1]}`
    };

    const datFile = find(`${stem}.dat`);
    const settings = datFile ? runDatasets(await datFile.sourceFile.text()) : null;
    const dataset = settings ? xrayDataset(settings, number) : null;
    let inputFile = null;
    let partialsFile = null;
    if (dataset?.file) {
        names.input = dataset.file.replace(/\\/g, '/').split('/').pop();
        const inputPath = resolveRelative(dir, dataset.file);
        inputFile = inputPath === null ? null : byPath.get(inputPath) || null;
        if (dataset.convolve) {
            names.partials = `${stem}_PDFpartials.csv`;
            partialsFile = find(names.partials);
        }
    }
    const fqFile = find(candidates[0]) || find(candidates[1]);
    if (fqFile) names.fq = fqFile.name;

    const key = [ftFile, datFile, inputFile, fqFile, partialsFile].map(versionOf).join('|');
    if (cache.has(key)) return cache.get(key);

    let payload;
    try {
        const transform = readRmcCsv(await ftFile.sourceFile.text(), ftFile.name);
        let q = null;
        let f = null;
        if (fqFile) {
            const series = readRmcCsv(await fqFile.sourceFile.text(), fqFile.name);
            q = series.data[0];
            f = experimentalColumn(series);
        }
        const inputColumns = inputFile ? readStogXy(await inputFile.sourceFile.text()) : null;
        let radius = null;
        if (partialsFile) {
            const x = partialsFile.plotData?.series?.[0]?.x
                || readRmcCsv(await partialsFile.sourceFile.text(), partialsFile.name).data[0];
            radius = largestFinite(x);
        }
        payload = assessFtUncertainty({
            r: transform.data[0], g: experimentalColumn(transform), q, f, settings, datasetNumber: number,
            inputColumns, radius, names
        });
    } catch (error) {
        payload = ftFailure(`could not compute the error ribbon: ${error?.message || error}`);
    }
    return remember(key, payload);
};

/**
 * `files` with every parsed xPDF transform's plotData carrying its
 * `uncertainty` (as Flask's /api/plot/data sends it), from `runFiles`, every
 * file of the run (makeRunFromEntries' allFiles; `files` when absent). A file
 * whose payload did not change is returned as it was.
 */
export const attachFtUncertainty = async (files, runFiles = files) => Promise.all(files.map(async (file) => {
    if (file.plotKind !== 'xpdf' || !file.plotData || !FT_FILE_RE.test(file.name) || !file.sourceFile) return file;
    let uncertainty;
    try {
        uncertainty = await ftUncertaintyForRun(file, runFiles);
    } catch (error) {
        uncertainty = ftFailure(`could not compute the error ribbon: ${error?.message || error}`);
    }
    if (file.plotData.uncertainty === uncertainty) return file;
    return { ...file, plotData: { ...file.plotData, uncertainty } };
}));
