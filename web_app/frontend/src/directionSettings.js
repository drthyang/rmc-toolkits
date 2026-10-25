// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Displacement Directions page's settings: its defaults, its resolution
// choices and the two engine requests it draws from. Read by the page
// (OrientationPage through useDirectionHistogram); the validation's arm A
// pins its own frozen settings (validation/armA.js, protocol 2026-10-24.6)
// and builds its request with directionReadoutRequest.
//
// - The page opens at ν = 10 (1002 cells) with the map smoothed 4× (four
//   neighbour passes), the user's choice (2026-10-24, decision record
//   2026-10-24-directions-page-defaults). Auto stays a resolution choice:
//   the engine's recommended frequency, the finest ν whose cells still
//   average at least 12 displacements (recommendedFrequency in
//   workers/orientation.js), the tiling analyze_displacement_directions
//   reads, so the page and the AI Copilot read the same tiling only when the
//   user picks Auto.
// - Every readout comes from the raw histogram. The engines pick the peak
//   cell on the smoothed map, so with smoothing on the peak's cell, direction
//   and σ would follow the Smoothing slider. The readout request therefore
//   carries the page's analysis settings at smoothing 0, as the tool's does.
//   The map request is the same request smoothed for display, so smoothing
//   only redraws the map: a display setting only, no readout, and no score,
//   depends on it.

/** The page's opening settings. `smoothing` is a display setting only. */
export const DIRECTION_DEFAULTS = Object.freeze({
    frequency: 10,
    weight: 'count',
    frame: 'cartesian',
    minQuantile: 0,
    smoothing: 4
});

/** Resolution choices: 'auto', or a geodesic frequency ν (10ν² + 2 cells). */
export const FREQUENCY_OPTIONS = Object.freeze(['auto', 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24]);

/** The cell count of the frequency-ν Goldberg tiling. */
export const cellCountOf = (frequency) => 10 * frequency * frequency + 2;

/**
 * The engine parameters of the page's readouts (the peak, its direction and
 * σ, the map, ± asymmetry and anisotropy tests) for one site: its analysis
 * settings at smoothing 0. 'auto' is an absent frequency (null), which the
 * engine resolves to its recommended frequency.
 */
export const directionReadoutRequest = ({ referenceNumber, frequency, weight, frame, minQuantile }) => ({
    referenceNumber,
    frequency: frequency === 'auto' ? null : frequency,
    weight,
    frame,
    minAmplitudeQuantile: minQuantile,
    smoothing: 0
});

/** The map the page draws: the readout request, smoothed for display. */
export const directionMapRequest = (settings) => ({
    ...directionReadoutRequest(settings),
    smoothing: settings.smoothing
});

/**
 * The Auto choice's label: the frequency the engine recommends for the site on
 * screen (its readout), so the user sees which tiling Auto stands for, also
 * while another resolution is picked. Plain "Auto" before a readout.
 */
export const autoResolutionLabel = (readout) => {
    const frequency = readout?.recommendedFrequency;
    return Number.isInteger(frequency) && frequency > 0
        ? `Auto · ν=${frequency} (${cellCountOf(frequency)})`
        : 'Auto';
};
