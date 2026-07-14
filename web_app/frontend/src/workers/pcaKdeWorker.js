// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Static-mode worker for the PCA / thermal-ellipsoid KDE. It parses the loaded
// `.rmc6f` file into per-site displacement clouds once (cached by content), then
// answers two request kinds off the main thread:
//   { kind: 'sites' }  -> anisotropic displacement tensors for every site
//   { kind: 'kde', referenceNumber | element, ...options } -> one PCA-KDE volume
// mirroring the Flask /api/pca/sites and /api/pca/kde endpoints so both runtime
// modes present the identical shape to the UI.

import {
    siteDisplacementsFromRmc6f,
    siteEllipsoids,
    sitePcaKde
} from './pcaKde.js';

// Parsing the whole configuration into clouds is the expensive part, so cache it
// per file text; selecting a site or changing bandwidth then only re-runs the
// (fast) volume solve.
let cache = { key: null, parsed: null };

const parseCached = (key, text) => {
    if (cache.key !== key || !cache.parsed) {
        cache = { key, parsed: siteDisplacementsFromRmc6f(text) };
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
        probability,
        sites: ellipsoids
    };
};

export const handlePcaMessage = async (data, getText) => {
    const { kind = 'kde', cacheKey, probability = 0.5 } = data;
    const text = await getText();
    const parsed = parseCached(cacheKey ?? text.length, text);

    if (kind === 'sites') {
        return summarizeSites(parsed, probability);
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
