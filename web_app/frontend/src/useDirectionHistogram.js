// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Displacement Directions page's two engine results for the selected site
// (directionSettings.js builds both requests):
//
// - `readout`: the raw histogram at the page's analysis settings (smoothing 0),
//   what every number under the sphere is read from. It is exactly what
//   analyze_displacement_directions computes at the same settings.
// - `map`: what the sphere draws, the same request with the page's smoothing;
//   the readout itself when smoothing is 0.
//
// Both land together, so the drawn map and its readouts never come from two
// different sites or settings. A smoothing change asks only for the map: the
// last readout answers the same request, so it is reused (nothing in it
// depends on smoothing). A new `requestPca` (another run, or a new
// configuration under Flask Live Data) always asks again.

import { useEffect, useRef, useState } from 'react';
import { directionMapRequest, directionReadoutRequest } from './directionSettings';

const EMPTY = Object.freeze({ readout: null, map: null, ref: null });

export default function useDirectionHistogram({ requestPca, ready, selectedRef, clusterThreshold, settings }) {
    const { frequency, weight, frame, minQuantile, smoothing } = settings;
    const [loaded, setLoaded] = useState(EMPTY);
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState(null);
    // The last readout, with the request it answered and the requestPca it came through.
    const lastReadout = useRef({ key: null, requestPca: null, data: null });

    useEffect(() => {
        let cancelled = false;
        const load = async () => {
            if (!ready || selectedRef == null) {
                setLoaded(EMPTY);
                setError(null);
                setLoading(false);
                return;
            }
            setLoading(true);
            setError(null);
            const site = { referenceNumber: selectedRef, frequency, weight, frame, minQuantile, smoothing };
            // clusterThreshold rides along for the browser engine's site
            // reconstruction; geometry: the readout is drawn when smoothing is 0.
            const extra = { clusterThreshold, geometry: true };
            const readoutParams = { ...directionReadoutRequest(site), ...extra };
            const key = JSON.stringify(readoutParams);
            const last = lastReadout.current;
            try {
                const readoutPromise = last.key === key && last.requestPca === requestPca
                    ? Promise.resolve(last.data)
                    : requestPca('orientation', readoutParams);
                const mapPromise = smoothing > 0
                    ? requestPca('orientation', { ...directionMapRequest(site), ...extra })
                    : readoutPromise;
                const [readout, map] = await Promise.all([readoutPromise, mapPromise]);
                if (cancelled) return;
                lastReadout.current = { key, requestPca, data: readout };
                setLoaded({ readout, map, ref: selectedRef });
            } catch (requestError) {
                if (cancelled) return;
                // A cluster-distance change can transiently request a re-numbered
                // site (the sites effect re-selects a valid one a moment later):
                // keep the current sphere rather than surfacing it.
                if (/Unknown reference number/i.test(requestError.message || '')) return;
                setLoaded(EMPTY);
                setError(requestError.message);
            } finally {
                if (!cancelled) setLoading(false);
            }
        };
        load();
        return () => { cancelled = true; };
    }, [requestPca, ready, selectedRef, frequency, weight, frame, smoothing, minQuantile, clusterThreshold]);

    return { readout: loaded.readout, map: loaded.map, loadedRef: loaded.ref, loading, error };
}
