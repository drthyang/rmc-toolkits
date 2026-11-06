// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Bond Geometry page's hover highlight (angleHighlight.js), wired up:
// after a Compute on a run opened in the browser, the request goes to its own
// worker (angleTripletsWorker.js) with the folded cell's drawn region; the
// answer is kept only if every bin's angle count equals the histogram's, so
// a highlight is never from another configuration or request. The hovered
// bin (onHoverX: the plot's snapped x, or null) then picks its triplets.
//
// A typed server directory has no structure text in the browser, so there
// the plot still marks the hovered bin but nothing is highlighted
// (`available` false).

import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { binIndexOf, highlightForBin } from './angleHighlight';

const sameCounts = (data, result) => (
    data.bins.length === result.counts.length
    && data.bins.every((bin, index) => bin.total === result.counts[index])
);

export default function useAngleHighlight({ rmc6fText, result, request, network, sites, createWorker }) {
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
        if (!rmc6fText || !result || !request || !result.angleCount || !hasNetwork) return undefined;
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
    const bin = result && Number.isFinite(hoveredX) ? binIndexOf(hoveredX, result.binWidth, result.counts.length) : null;
    const highlight = useMemo(
        () => (data && bin != null && network?.unit ? highlightForBin(data, bin, siteFractional, network.unit) : null),
        [data, bin, siteFractional, network]
    );
    const onHoverX = useCallback((x) => setHoveredX(x), []);
    return { highlight, onHoverX, available: Boolean(rmc6fText), ready: Boolean(data) };
}
