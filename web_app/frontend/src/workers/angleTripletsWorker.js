// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Bond Geometry page's hover highlight, off the main thread: after a
// Compute, the page posts the same request and the folded cell's drawn
// region, and this worker answers with every bin's triplets folded onto the
// average sites (angleHighlight.js siteAngleTriplets). A worker of its own,
// so it never queues ahead of the next Compute in the shared analysis
// worker; the page terminates it when a newer request makes it stale.
//
//   { id, text, params, region, clusterThreshold } -> { id, result } | { id, error }

import { siteAngleTriplets } from '../angleHighlight.js';
import { DEFAULT_CLUSTER_THRESHOLD, siteDisplacementsFromRmc6f } from './pcaKde.js';
import { requestObject } from './requestGuards.js';
import { tripletOptionsFromRequest } from './triplets.js';

let cache = { text: null, threshold: null, parsed: null };

export const handleAngleTriplets = (data) => {
    const { text, params, region = null, clusterThreshold = DEFAULT_CLUSTER_THRESHOLD } = requestObject(data);
    if (typeof text !== 'string') throw new Error('No structure text to read the triplets from');
    if (cache.text !== text || cache.threshold !== clusterThreshold) {
        cache = { text, threshold: clusterThreshold, parsed: siteDisplacementsFromRmc6f(text, { clusterThreshold }) };
    }
    const { parsed } = cache;
    const siteFractional = new Map(parsed.sites.map((site) => [site.referenceNumber, site.siteFractional]));
    return siteAngleTriplets(parsed.atomList, parsed.latticeVectors, parsed.supercell, siteFractional,
        tripletOptionsFromRequest(params ?? {}), region);
};

if (typeof self !== 'undefined' && typeof self.postMessage === 'function') {
    self.onmessage = (event) => {
        const id = event?.data?.id;
        try {
            self.postMessage({ id, result: handleAngleTriplets(event.data) });
        } catch (error) {
            self.postMessage({ id, error: error?.message || 'Could not find the triplets' });
        }
    };
}
