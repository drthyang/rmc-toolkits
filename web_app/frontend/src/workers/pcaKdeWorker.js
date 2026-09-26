// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Static-mode worker for the PCA / thermal-ellipsoid KDE. It parses the loaded
// `.rmc6f` file into per-site displacement clouds once (cached by content), then
// answers three request kinds off the main thread:
//   { kind: 'sites' }  -> anisotropic displacement tensors for every site
//   { kind: 'kde', referenceNumber | element, ...options } -> one PCA-KDE volume
//   { kind: 'orientation', referenceNumber | element, ...options } -> hex-binned
//     solid-angle histogram of the displacement directions
//   { kind: 'triplets', triplet, bond12, ...options } -> bond-angle summary
// mirroring the Flask /api/pca/sites, /api/pca/kde, /api/pca/orientation and
// /api/triplets endpoints so both runtime modes present the identical shape
// to the UI.

import {
    DEFAULT_CLUSTER_THRESHOLD,
    siteDisplacementsFromRmc6f,
    siteEllipsoids,
    sitePcaKde,
    siteSpecies
} from './pcaKde.js';
import { siteOrientationHistogram } from './orientation.js';
import { APP_MAX_ANGLES, bondAngleSummary, isBlankValue } from './triplets.js';
import { assertFiniteResult, requestNumber, requestObject } from './requestGuards.js';

// Parsing the whole configuration into clouds is the expensive part, so cache it
// and re-parse only when the .rmc6f text itself changes. The key MUST come from
// the text content, not a file path: two runs can share a filename, and a
// path-based key would hand back the previous model's clouds for a new dataset.
let cache = { key: null, parsed: null };

// Cheap content signature (FNV-1a over a stride sample + full length). Distinct
// structure files differ in length and/or sampled bytes, so this re-parses on a
// real dataset change while staying a cache hit for the same text.
const textSignature = (text) => {
    let hash = 0x811c9dc5;
    for (let i = 0; i < text.length; i += 64) {
        hash = Math.imul(hash ^ text.charCodeAt(i), 0x01000193);
    }
    return `${text.length}:${(hash >>> 0).toString(36)}`;
};

// The clustering threshold changes the reconstructed sites of an old (coords-only)
// file, so it is part of the cache key; for files with real site columns it is inert
// and every threshold hits the same cached parse.
const parseCached = (text, clusterThreshold, name) => {
    const key = `${textSignature(text)}@${clusterThreshold}`;
    if (cache.key !== key || !cache.parsed) {
        cache = { key, parsed: siteDisplacementsFromRmc6f(text, { clusterThreshold, name }) };
    }
    return cache.parsed;
};

// As app.MAX_ORIENTATION_SMOOTHING (the UI offers 0-12).
const MAX_ORIENTATION_SMOOTHING = 64;
const DECIMAL = /^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$/;

// As app._bw_argument: blank -> 'scott', the two names in any case, and a
// numeric string -> its number (anything else reaches sitePcaKde's own error).
const bwArgument = (raw) => {
    if (raw == null || (typeof raw === 'string' && !raw.trim())) return 'scott';
    if (typeof raw !== 'string') return raw;
    const text = raw.trim();
    const name = text.toLowerCase();
    if (name === 'scott' || name === 'silverman') return name;
    return DECIMAL.test(text) ? Number(text) : raw;
};

const summarizeSites = (parsed, probability) => {
    const ellipsoids = siteEllipsoids(parsed.sites, probability);
    return {
        referenceNumbers: parsed.referenceNumbers,
        // Every species present, including the minority species of a mixed-
        // occupancy site (site labels carry only the majority) -- as /api/pca/sites.
        elements: siteSpecies(parsed.sites),
        totalAtoms: parsed.sites.reduce((sum, site) => sum + site.count, 0),
        latticeVectors: parsed.latticeVectors,
        supercell: parsed.supercell,
        // True when the sites were reconstructed by folding + clustering because the
        // file lacked site/cell columns; the UI shows the threshold knob + count/N.
        reconstructed: Boolean(parsed.reconstructed),
        probability,
        sites: ellipsoids,
        // Atom lines the shared .rmc6f grammar skipped (non-finite coordinates,
        // unparsed lines, a header count mismatch), or null -- the same text as
        // /api/pca/sites' parseWarning.
        parseWarning: parsed.parseWarning ?? null
    };
};

const KINDS = ['sites', 'kde', 'orientation', 'triplets'];

export const handlePcaMessage = async (data, getText) => {
    const { kind = 'kde', probability = 0.5, clusterThreshold = DEFAULT_CLUSTER_THRESHOLD } = requestObject(data);
    if (!KINDS.includes(kind)) {
        throw new Error(`unknown request kind '${kind}' (expected ${KINDS.join(', ')})`);
    }
    const text = await getText();
    const name = data.file?.name || data.file?.sourceFile?.name || 'structure file';
    const parsed = parseCached(text, clusterThreshold, name);

    if (kind === 'sites') {
        return summarizeSites(parsed, probability);
    }

    if (kind === 'triplets') {
        // Flat scalar params (end1/apex/end2, r12Min...) so the identical
        // request shape works as an HTTP query string against /api/triplets.
        // Same request caps as the Flask route (the engine is unrestricted):
        // rmax <= 15 A bounds the neighbour search, binWidth >= 0.05 deg the
        // response, and the shared APP_MAX_ANGLES budget the pairing work —
        // the angle count grows ~rmax^6, so without it a window well inside
        // the rmax cap would freeze the shared worker.
        // A null/''/blank bound is *missing* — the same errors as the Flask
        // route — never Number()'d into 0, which silently widened a window
        // with a cleared rmin box down to 0 A. Bounds then go to the engine
        // unconverted; its validateWindow rejects anything non-numeric.
        const window = (minKey, maxKey, fallback) => {
            const missing = [minKey, maxKey].filter((key) => isBlankValue(data[key]));
            if (missing.length === 2 && fallback !== undefined) return fallback;
            if (missing.length) {
                throw new Error(`${minKey}/${maxKey} are required together; missing ${missing[0]}`);
            }
            return [data[minKey], data[maxKey]];
        };
        const bond12 = window('r12Min', 'r12Max');
        const bond23 = window('r23Min', 'r23Max', null);
        // Absent → the 1° default (as the route's query default); present
        // but blank → an error, never a 0-degree bin width.
        if (typeof data.binWidth === 'string' && isBlankValue(data.binWidth)) {
            throw new Error(`binWidth must be a number, got '${data.binWidth}'`);
        }
        const binWidth = data.binWidth == null ? 1.0 : Number(data.binWidth);
        if (binWidth < 0.05) {
            throw new Error(`binWidth is capped at >= 0.05 deg, got ${binWidth}`);
        }
        for (const [key, raw] of [['r12Max', data.r12Max], ['r23Max', data.r23Max]]) {
            if (!isBlankValue(raw) && Number(raw) > 15) {
                throw new Error(`${key} is capped at 15 A, got ${raw}`);
            }
        }
        const summary = bondAngleSummary(
            parsed.atomList.fractional,
            parsed.atomList.elements,
            parsed.latticeVectors,
            {
                triplet: [data.end1, data.apex, data.end2],
                bond12,
                bond23,
                binWidth,
                maxAngles: APP_MAX_ANGLES
            }
        );
        // Atom lines the shared .rmc6f grammar skipped, or null -- the same
        // text as /api/triplets' parseWarning (bond_angle_summary_from_file).
        return { ...summary, parseWarning: parsed.parseWarning ?? null };
    }

    if (kind === 'orientation') {
        // ''/'all' mean "every site pooled", normalised to null exactly as the
        // Flask route does, so both transports return the same payload shape.
        const element = data.element === '' || data.element === 'all' ? null : data.element ?? null;
        // As /api/pca/orientation's _query_number rules: an integer frequency
        // (its [1, 64] range is the tiling's own check) and an integer
        // smoothing in [0, 64] -- 1e9 passes would pin the worker for minutes.
        const frequency = requestNumber(data.frequency, 'frequency', { fallback: null, integer: true });
        const smoothing = requestNumber(data.smoothing, 'smoothing', {
            fallback: 0, integer: true, ge: 0, le: MAX_ORIENTATION_SMOOTHING
        });
        const histogram = siteOrientationHistogram(parsed, {
            referenceNumber: data.referenceNumber ?? null,
            element,
            frequency,
            weight: data.weight ?? 'count',
            minAmplitude: data.minAmplitude ?? 0,
            minAmplitudeQuantile: data.minAmplitudeQuantile ?? 0,
            smoothing,
            frame: data.frame ?? 'cartesian',
            geometry: data.geometry ?? true
        });
        // As /api/pca/orientation: the parse warning rides along, and a
        // result holding NaN/Infinity is an error (_strict_result_response).
        return assertFiniteResult({ ...histogram, parseWarning: parsed.parseWarning ?? null });
    }

    // ''/'all' mean "every site pooled" -> null, as /api/pca/kde normalises
    // them (and the orientation branch above), so both transports agree.
    const element = data.element === '' || data.element === 'all' ? null : data.element ?? null;
    // Extreme but finite bw / extent underflow or overflow the kernels: an
    // error with the Flask message, never a posted NaN volume.
    return assertFiniteResult(sitePcaKde(parsed, {
        referenceNumber: data.referenceNumber ?? null,
        element,
        bw: bwArgument(data.bw),
        bwScale: data.bwScale ?? 1,
        grid: data.grid ?? 48,
        extent: data.extent ?? 3,
        cubicBox: data.cubicBox ?? false,
        probability,
        projections: data.projections ?? true
    }));
};

// Guarded so the module can be imported by tests outside a worker context.
if (typeof self !== 'undefined' && typeof self.postMessage === 'function') {
    self.onmessage = async (event) => {
        const data = event?.data;
        const id = data !== null && typeof data === 'object' ? data.id : undefined;
        try {
            const { file, text } = requestObject(data);
            const getText = async () => {
                if (typeof text === 'string') return text;
                if (file?.sourceFile) return file.sourceFile.text();
                throw new Error('No browser structure file available');
            };
            const result = await handlePcaMessage(data, getText);
            self.postMessage({ id, result });
        } catch (error) {
            self.postMessage({ id, error: error.message || 'Browser PCA-KDE computation failed' });
        }
    };
}
