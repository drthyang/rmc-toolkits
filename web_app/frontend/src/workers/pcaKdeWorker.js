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
    sitePcaKde
} from './pcaKde.js';
import { siteOrientationHistogram } from './orientation.js';
import { APP_MAX_ANGLES, bondAngleSummary, isBlankValue } from './triplets.js';

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
const parseCached = (text, clusterThreshold) => {
    const key = `${textSignature(text)}@${clusterThreshold}`;
    if (cache.key !== key || !cache.parsed) {
        cache = { key, parsed: siteDisplacementsFromRmc6f(text, { clusterThreshold }) };
    }
    return cache.parsed;
};

const summarizeSites = (parsed, probability) => {
    const ellipsoids = siteEllipsoids(parsed.sites, probability);
    return {
        referenceNumbers: parsed.referenceNumbers,
        elements: [...new Set(parsed.sites.map((site) => site.element))].sort(),
        totalAtoms: parsed.sites.reduce((sum, site) => sum + site.count, 0),
        latticeVectors: parsed.latticeVectors,
        supercell: parsed.supercell,
        // True when the sites were reconstructed by folding + clustering because the
        // file lacked site/cell columns; the UI shows the threshold knob + count/N.
        reconstructed: Boolean(parsed.reconstructed),
        probability,
        sites: ellipsoids
    };
};

export const handlePcaMessage = async (data, getText) => {
    const { kind = 'kde', probability = 0.5, clusterThreshold = DEFAULT_CLUSTER_THRESHOLD } = data;
    const text = await getText();
    const parsed = parseCached(text, clusterThreshold);

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
        return bondAngleSummary(
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
    }

    if (kind === 'orientation') {
        return siteOrientationHistogram(parsed, {
            referenceNumber: data.referenceNumber ?? null,
            element: data.element ?? null,
            frequency: data.frequency ?? null,
            weight: data.weight ?? 'count',
            minAmplitude: data.minAmplitude ?? 0,
            minAmplitudeQuantile: data.minAmplitudeQuantile ?? 0,
            smoothing: data.smoothing ?? 0,
            frame: data.frame ?? 'cartesian',
            geometry: data.geometry ?? true
        });
    }

    return sitePcaKde(parsed, {
        referenceNumber: data.referenceNumber ?? null,
        element: data.element ?? null,
        bw: data.bw ?? 'scott',
        bwScale: data.bwScale ?? 1,
        grid: data.grid ?? 48,
        extent: data.extent ?? 3,
        cubicBox: data.cubicBox ?? false,
        probability,
        projections: data.projections ?? true
    });
};

// Guarded so the module can be imported by tests outside a worker context.
if (typeof self !== 'undefined' && typeof self.postMessage === 'function') {
    self.onmessage = async (event) => {
        const { id, file, text } = event.data;
        try {
            const getText = async () => {
                if (typeof text === 'string') return text;
                if (file?.sourceFile) return file.sourceFile.text();
                throw new Error('No browser structure file available');
            };
            const result = await handlePcaMessage(event.data, getText);
            self.postMessage({ id, result });
        } catch (error) {
            self.postMessage({ id, error: error.message || 'Browser PCA-KDE computation failed' });
        }
    };
}
