// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// A server folder's termination-ripple source (copilotRun.js LoadedRun
// `ripples`): Flask's /api/ripples, which runs rmc_toolkits.ripples on the
// files where they are. The query repeats `path` (URLSearchParams: axios's own
// serialisation would write `path[]`); the server's message is the error.

import axios from 'axios';
import API_BASE_URL from './api';

const abortError = () => {
    const error = new Error('The request was cancelled');
    error.name = 'AbortError';
    return error;
};

/** The /api/ripples query of a call: its fit files' paths and the arguments. */
export const ripplesQuery = (paths, args = {}) => {
    const query = new URLSearchParams();
    for (const path of paths || []) query.append('path', path);
    for (const key of ['qmax', 'qmin', 'window']) {
        if (args[key] !== null && args[key] !== undefined) query.append(key, String(args[key]));
    }
    return query;
};

/** (paths, args, { signal }) → the server's result. */
export const flaskRipplesSource = () => (paths, args = {}, { signal } = {}) => axios
    .get(`${API_BASE_URL}/api/ripples`, { params: ripplesQuery(paths, args), ...(signal ? { signal } : {}) })
    .then((response) => response.data)
    .catch((error) => {
        if (signal?.aborted) throw abortError();
        const message = error?.response?.data?.error;
        throw message ? new Error(message) : error;
    });
