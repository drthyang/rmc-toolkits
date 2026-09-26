// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Shared data plumbing for pages built on per-site RMC displacement clouds
// (PCA Ellipsoid, Orientation): reads the local `.rmc6f` text once per file,
// routes requests to the PCA worker (browser-loaded runs, both runtimes) or
// the Flask /api/pca/* endpoints (typed backend directories), loads the
// per-site ellipsoid table, and manages the selected reference site.
//
// Extracted from PcaKdePage when the Orientation view became its own page:
// the two pages are scientifically independent (the orientation histogram is
// not a PCA product) but share exactly this plumbing and the site picker.

import { useCallback, useEffect, useMemo, useState } from 'react';
import axios from 'axios';
import API_BASE_URL from './api';
import { unitCellVectors } from './pcaCrystalFrame';

// One worker for the whole app, never torn down: its parse cache is keyed by
// file content (+ cluster threshold), so pages share the parsed model and
// switching between them never re-parses the configuration.
let sharedWorker = null;
let nextRequestId = 0;
const getWorker = () => {
    if (!sharedWorker) {
        sharedWorker = new Worker(new URL('./workers/pcaKdeWorker.js', import.meta.url), { type: 'module' });
    }
    return sharedWorker;
};

// `dataEpoch` is App.jsx's Flask-mode configuration epoch: it changes when
// the run's .rmc6f changes on disk. A typed backend directory keeps the same
// `directory` across such a save, so the epoch is what tells this hook (and
// every consumer of `requestPca`) to re-read. A browser-loaded run needs no
// epoch — its new text changes `rmc6fText` instead.
export default function useSiteCloud({ directory, localRun, probability = 0.5, clusterThreshold, dataEpoch = 0 }) {
    // The loaded .rmc6f text, tagged with the file it came from, so a
    // just-changed dataset never runs against the previous model's text.
    const [loadedText, setLoadedText] = useState({ file: null, text: null });
    const [sites, setSites] = useState(null);
    const [selectedRef, setSelectedRef] = useState(null);
    const [sitesError, setSitesError] = useState(null);
    const [loadingSites, setLoadingSites] = useState(false);

    const structureFile = localRun?.structureFile || null;
    // A locally-loaded run (the Demo, or a picked folder) carries its .rmc6f as
    // a browser file, parsed in the worker in BOTH runtimes; only a typed
    // backend directory goes through the Flask API.
    const localFile = structureFile?.sourceFile || null;
    // Only treat the text as current when it was read from the file in props
    // right now; on a dataset switch this is null until the new text loads.
    const rmc6fText = loadedText.file === localFile ? loadedText.text : null;
    // Stable per-dataset identity (see PcaKdePage's framing logic): changes on
    // a genuinely different dataset, not on a Live Data refresh of the same run.
    const datasetKey = localRun?.runId ?? directory ?? null;

    // --- Read the raw .rmc6f text once per local file. ------------------------
    useEffect(() => {
        let cancelled = false;
        if (localFile) {
            localFile.text().then((text) => {
                if (!cancelled) setLoadedText({ file: localFile, text });
            }).catch(() => {
                if (!cancelled) { setLoadedText({ file: null, text: null }); setSitesError('Could not read the structure file.'); }
            });
        } else {
            setLoadedText({ file: null, text: null });
        }
        return () => { cancelled = true; };
    }, [localFile]);

    // A single request path for both data sources: the shared worker when a
    // local file is loaded, otherwise the Flask API against the run directory.
    // Its identity is the configuration's identity: it changes with the
    // loaded text, the directory, or (Flask Live Data) the configuration
    // epoch, so every effect keyed on it — the site table below, the PCA
    // page's KDE volume, the orientation histogram — re-reads in place when
    // a new .rmc6f is saved, while the pages keep their picks and settings.
    const requestPca = useCallback((kind, params) => {
        if (rmc6fText) {
            return new Promise((resolve, reject) => {
                const worker = getWorker();
                nextRequestId += 1;
                const id = nextRequestId;
                const handler = (event) => {
                    if (event.data.id !== id) return;
                    worker.removeEventListener('message', handler);
                    if (event.data.error) reject(new Error(event.data.error));
                    else resolve(event.data.result);
                };
                worker.addEventListener('message', handler);
                worker.postMessage({ id, kind, text: rmc6fText, ...params });
            });
        }
        const endpoint = {
            sites: '/api/pca/sites',
            orientation: '/api/pca/orientation',
            triplets: '/api/triplets'
        }[kind] ?? '/api/pca/kde';
        return axios
            .get(`${API_BASE_URL}${endpoint}`, { params: { dir: directory || '.', ...params } })
            .then((response) => response.data)
            .catch((error) => {
                // Surface the server's own message ("displacement cloud has zero
                // spread", a non-finite coordinate, ...) rather than axios's generic
                // "Request failed with status code 400".
                const message = error?.response?.data?.error;
                throw message ? new Error(message) : error;
            });
    // dataEpoch is deliberately a dependency the body never reads: it gives
    // the callback (and so every effect keyed on it) a new identity when the
    // backend's .rmc6f changes under an unchanged directory.
    // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [rmc6fText, directory, dataEpoch]);

    // --- Load the per-site ellipsoid table. -----------------------------------
    useEffect(() => {
        let cancelled = false;
        const loadSites = async () => {
            // Wait for a local file's text before requesting; a backend
            // directory needs no text and proceeds immediately.
            if (localFile && !rmc6fText) { setSites(null); return; }
            setLoadingSites(true);
            setSitesError(null);
            try {
                const data = await requestPca('sites', { probability, clusterThreshold });
                if (cancelled) return;
                setSites(data);
                setSelectedRef((current) => {
                    const refs = data.sites.map((site) => site.referenceNumber);
                    if (current != null && refs.includes(current)) return current;
                    // Prefer a clean reconstructed site (members == one-per-cell) so
                    // the default view is a genuine thermal site, not a disordered
                    // shell. Inert for files with real site columns.
                    const clean = data.sites.find((site) => site.copiesPerCell && site.count === site.copiesPerCell);
                    return clean?.referenceNumber ?? refs[0] ?? null;
                });
            } catch (error) {
                if (!cancelled) { setSites(null); setSitesError(error.message); }
            } finally {
                if (!cancelled) setLoadingSites(false);
            }
        };
        loadSites();
        return () => { cancelled = true; };
    }, [requestPca, localFile, rmc6fText, probability, clusterThreshold]);

    const selectedEllipsoid = useMemo(
        () => sites?.sites.find((site) => site.referenceNumber === selectedRef) || null,
        [sites, selectedRef]
    );

    // Unit-cell edge vectors (a, b, c) in the shared Cartesian frame; null
    // until the sites (and their lattice metadata) have loaded.
    const unitCell = useMemo(
        () => (sites?.latticeVectors && sites?.supercell
            ? unitCellVectors(sites.latticeVectors, sites.supercell)
            : null),
        [sites]
    );

    return {
        sites,
        sitesError,
        loadingSites,
        selectedRef,
        setSelectedRef,
        selectedEllipsoid,
        requestPca,
        localFile,
        rmc6fText,
        // True once requests can be issued (no local file, or its text loaded).
        ready: !localFile || Boolean(rmc6fText),
        unitCell,
        datasetKey
    };
}
