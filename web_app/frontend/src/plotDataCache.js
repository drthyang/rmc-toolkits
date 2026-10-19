// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Flask mode's plot data, fetched once per file and version: the Dashboard's
// charts (InteractivePlot) and the run the Dashboard publishes for the AI
// Copilot (Dashboard.jsx) both read GET /api/plot/data through here, so a
// chart and the Copilot share one request and one payload. Whichever asks
// first fetches; the other gets the same promise.
//
// The key is the file's path and its version (the Dashboard's "modified:size"
// refreshKey; for the χ² history, the version of every log the server
// concatenates). An entry holds one version per path, and a newer version
// replaces it. A failed request is not kept, so the next asker tries again.
// An undefined version is never cached.
//
// The payload is the browser parser's form (plotDataFromText): a masked point,
// which the strict JSON provider sends as null, is NaN again.

import axios from 'axios';
import API_BASE_URL from './api';
import { plotPayloadError } from './plotDomain';

// Paths kept: a folder's plot files, with room for the folder before.
const MAX_ENTRIES = 64;

const entries = new Map();

const nanForNull = (values) => (Array.isArray(values) ? values.map((value) => (value === null ? NaN : value)) : values);

/** A Flask plot payload as plotDataFromText gives it: null in a series' x or y is NaN. */
export const plotDataFromServer = (payload) => ({
    ...payload,
    series: payload.series.map((series) => ({ ...series, x: nanForNull(series.x), y: nanForNull(series.y) }))
});

const failure = (message, invalidPayload = false) => Object.assign(new Error(message), { invalidPayload });

const request = (path) => axios.get(`${API_BASE_URL}/api/plot/data`, { params: { path } }).then(
    (response) => {
        // A body that is not valid JSON arrives as a raw string: never hand
        // that on, say what went wrong instead.
        const payloadError = plotPayloadError(response.data);
        if (payloadError) throw failure(payloadError, true);
        return plotDataFromServer(response.data);
    },
    (error) => { throw failure(error.response?.data?.error || 'Failed to load interactive plot'); }
);

/**
 * The plot data of the file at `path`, at `version`: a promise of the payload
 * (plotDataFromServer's form), or a rejection with the message to show
 * (`invalidPayload` set when the server's answer was not plot data).
 */
export const fetchPlotData = (path, version) => {
    if (version == null) return request(path);
    const cached = entries.get(path);
    if (cached && cached.version === version) return cached.promise;
    const promise = request(path);
    entries.delete(path);
    entries.set(path, { version, promise });
    while (entries.size > MAX_ENTRIES) entries.delete(entries.keys().next().value);
    promise.catch(() => {
        if (entries.get(path)?.promise === promise) entries.delete(path);
    });
    return promise;
};

/** Forget every cached payload (tests). */
export const clearPlotDataCache = () => entries.clear();
