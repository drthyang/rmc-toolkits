// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Bond Geometry page's hover highlight (angleHighlight.js), wired up:
// after a Compute on a run opened in the browser, the request goes to its own
// worker (angleTripletsWorker.js) with the folded cell's drawn region; the
// answer is kept only if every bin's angle count equals the histogram's, so
// a highlight is never from another configuration or request. The angle
// range dragged on the plot (`range`, degrees, at bin edges) picks its
// triplets; with no range, the hovered bin (onHoverX: the plot's snapped x,
// or null) previews its own.
//
// A typed server directory has no structure text in the browser, so there
// the plot still marks the hovered bin but nothing is highlighted
// (`available` false).

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { binIndexOf, binsOfRange, highlightForBins } from './angleHighlight';

const sameCounts = (data, result) => (
    Array.isArray(data?.bins) && data.bins.length === result.counts.length
    && data.bins.every((bin, index) => bin.total === result.counts[index])
);

export default function useAngleHighlight({ rmc6fText, result, request, network, sites, range = null, createWorker }) {
    const [found, setFound] = useState(null); // { result, data }
    const [hoveredX, setHoveredX] = useState(null);
    const workerRef = useRef(null);
    const busyRef = useRef(false);

    useEffect(() => () => {
        workerRef.current?.terminate();
        workerRef.current = null;
    }, []);

    const region = network?.limited ? { center: network.center, radiusSq: network.radiusSq } : null;
    const regionKey = region ? `${region.center.join()}|${region.radiusSq}` : 'all';
    const hasNetwork = Boolean(network);
    useEffect(() => {
        if (!rmc6fText || !result || !request || !result.angleCount || !Array.isArray(result.counts) || !hasNetwork) return undefined;
        // A request still running is stale now: stop it rather than queue.
        if (busyRef.current) {
            workerRef.current?.terminate();
            workerRef.current = null;
        }
        if (!workerRef.current) workerRef.current = createWorker();
        const worker = workerRef.current;
        let cancelled = false;
        busyRef.current = true;
        const onMessage = (event) => {
            worker.removeEventListener('message', onMessage);
            busyRef.current = false;
            if (cancelled || event.data.error) return;
            if (sameCounts(event.data.result, result)) setFound({ result, data: event.data.result });
        };
        worker.addEventListener('message', onMessage);
        worker.postMessage({ id: 1, text: rmc6fText, params: request, region });
        return () => {
            cancelled = true;
            worker.removeEventListener('message', onMessage);
        };
    // region is rebuilt each render; regionKey is its identity.
    // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [rmc6fText, result, request, regionKey, hasNetwork, createWorker]);

    const data = found && found.result === result ? found.data : null;
    const siteFractional = useMemo(
        () => new Map((sites?.sites ?? []).map((site) => [site.referenceNumber, site.siteFractional])),
        [sites]
    );
    const nbins = Array.isArray(result?.counts) ? result.counts.length : 0;
    const rangeBins = nbins ? binsOfRange(range, result.binWidth, nbins) : null;
    const bin = nbins && Number.isFinite(hoveredX) ? binIndexOf(hoveredX, result.binWidth, nbins) : null;
    const [first, last] = rangeBins ?? (bin != null ? [bin, bin] : [null, null]);
    const isRange = Boolean(rangeBins);
    const highlight = useMemo(
        () => (data && first != null && network?.unit
            ? highlightForBins(data, first, last, siteFractional, network.unit, { range: isRange })
            : null),
        [data, first, last, isRange, siteFractional, network]
    );
    const onHoverX = useCallback((x) => setHoveredX(x), []);
    return { highlight, onHoverX, available: Boolean(rmc6fText), ready: Boolean(data) };
}
