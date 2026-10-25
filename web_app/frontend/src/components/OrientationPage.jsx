// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Orientation page: the solid-angle distribution of a site's displacement
// directions, hex-binned on a Goldberg sphere. A standalone workspace page —
// the histogram is not a PCA product, so it does not live on the PCA Ellipsoid
// page; the two share only the site picker, provided by the useSiteCloud hook.
//
// Design logic mirrors the PCA Ellipsoid page: all options live in the top
// controls bar, the three viewport panels sit in one equal-height grid (axis
// views : sphere : site picker = 3 : 6.5 : 6.5), and each panel's actions
// (frame toggle, Reset, Save) live in its own header.
//
// The settings open at DIRECTION_DEFAULTS (directionSettings.js): ν = 10
// with the map smoothed 4× (the user's choice, 2026-10-24). Auto, the
// engine's recommended frequency and the tiling the AI Copilot reads, shows
// the ν it chose in its label. Like every option here, a pick
// is page state: it stays when the site changes, the tab is left and the page
// returns, or another run opens (the page is never remounted), until the user
// changes it (Auto included) or a Copilot step's Show in sets it. The
// readouts come from the raw histogram whatever the smoothing
// (useDirectionHistogram).

import React, { useEffect, useMemo, useState } from 'react';
import { isStaticMode } from '../browserData';
import { buildElementColors, speciesCounts } from '../atomColors';
import { COLORMAP_NAMES } from '../colormaps';
import { Control, ControlGroup, ControlsBar, Hint, Page, PageIssues, Switch, useIssue } from '../ui';
import InfoBadge from '../ui/InfoBadge';
import OrientationView from './OrientationView';
import SiteStructurePanel from './SiteStructurePanel';
import useSiteCloud from '../useSiteCloud';
import useDirectionHistogram from '../useDirectionHistogram';
import { autoResolutionLabel, cellCountOf, DIRECTION_DEFAULTS, FREQUENCY_OPTIONS } from '../directionSettings';
import { siteLabel } from '../siteLabel';
import AppFooter from './AppFooter';
import './OrientationPage.css';

const numberFormat = (value, digits = 4) =>
    Number.isFinite(value) ? value.toFixed(digits) : '—';

const DEFAULT_CLUSTER_THRESHOLD = 1.5;

const WEIGHT_OPTIONS = [
    { value: 'count', label: 'Count' },
    { value: 'amplitude', label: '|Δr|' },
    { value: 'amplitude2', label: '|Δr|²' }
];

export default function OrientationPage({ directory, localRun, dataEpoch = 0, request = null }) {
    // Fold-and-cluster distance (Å), used only when the loaded file has no
    // reference-site/cell columns and its sites must be reconstructed.
    const [clusterThreshold, setClusterThreshold] = useState(DEFAULT_CLUSTER_THRESHOLD);

    // Histogram + display options (owned here, in the top controls bar).
    const [frequency, setFrequency] = useState(DIRECTION_DEFAULTS.frequency);
    const [weight, setWeight] = useState(DIRECTION_DEFAULTS.weight);
    const [frame, setFrame] = useState(DIRECTION_DEFAULTS.frame);
    const [smoothing, setSmoothing] = useState(DIRECTION_DEFAULTS.smoothing);
    const [minQuantile, setMinQuantile] = useState(DIRECTION_DEFAULTS.minQuantile);
    const [colormap, setColormap] = useState('viridis');
    const [contrast, setContrast] = useState(1);
    const [relief, setRelief] = useState(0.5);
    const [showOutline, setShowOutline] = useState(true);
    const [showAxes, setShowAxes] = useState(true);

    const {
        sites,
        sitesError,
        loadingSites,
        selectedRef,
        setSelectedRef,
        selectedEllipsoid,
        requestPca,
        localFile,
        ready,
        unitCell
    } = useSiteCloud({ directory, localRun, clusterThreshold, dataEpoch });

    // The selected site's readout (raw histogram) and drawn map.
    const histogram = useDirectionHistogram({
        requestPca,
        ready,
        selectedRef,
        clusterThreshold,
        settings: { frequency, weight, frame, minQuantile, smoothing }
    });
    // Auto's label names the ν the engine chose for the site on screen.
    const autoLabel = autoResolutionLabel(histogram.loadedRef === selectedRef ? histogram.readout : null);

    // "Show in Displacement Directions" from an AI Copilot step (App.jsx
    // `request`, a new `id` per click): select its site and, when it is one
    // of the offered resolutions, its sphere resolution ('auto' when the
    // Copilot let the engine choose). Smoothing, weight and frame stay as the
    // user set them, so the drawn map can differ from the Copilot's (which is
    // unsmoothed); the readouts use the raw counts either way.
    useEffect(() => {
        const params = request?.params;
        if (!params) return;
        if (params.referenceNumber != null) setSelectedRef(params.referenceNumber);
        if (FREQUENCY_OPTIONS.includes(params.frequency)) setFrequency(params.frequency);
    // Applied once per request (its id), never re-applied over the user's picks.
    // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [request?.id]);

    const staticMode = isStaticMode();
    const noRun = staticMode && !localFile;

    const siteRefs = useMemo(() => (sites?.sites || []).map((site) => site.referenceNumber), [sites]);

    const elementColors = useMemo(
        () => buildElementColors(sites?.elements ?? [], speciesCounts(sites?.sites)),
        [sites]
    );

    // The page's Problems section (the sphere card reports its own failures).
    useIssue('sites', sitesError, { source: 'Sites' });
    useIssue('sites-atoms-skipped', sitesError ? null : sites?.parseWarning, { severity: 'warning', source: 'Atoms skipped' });

    return (
        <Page as="div" column>
            <PageIssues />
            <ControlsBar>
                {/* Site & resolution */}
                <ControlGroup label="Site and resolution">
                    <Control
                        label={(
                            <>
                                Site
                                <InfoBadge label="About the site picker">
                                    <p>
                                        Each reference site (an RMCProfile reference number) is one
                                        crystallographic position. Only the <em>directions</em> of its
                                        per-atom displacements are analysed here — measured from the
                                        site&apos;s own mean position in this configuration, so an
                                        off-centring shared by every copy is not visible — and the
                                        amplitude enters solely through the optional weighting and the
                                        amplitude height.
                                    </p>
                                </InfoBadge>
                            </>
                        )}
                    >
                        <select className="ui-select"
                            value={selectedRef ?? ''}
                            onChange={(event) => setSelectedRef(Number(event.target.value))}
                            disabled={!sites}
                            aria-label="Site"
                        >
                            {sites?.sites.map((site) => (
                                <option key={site.referenceNumber} value={site.referenceNumber}>
                                    {`#${site.referenceNumber} ${siteLabel(site)} — U=${numberFormat(site.uIso, 4)} Å²`}
                                    {site.copiesPerCell ? ` (${site.count}/${site.copiesPerCell})` : ''}
                                </option>
                            ))}
                        </select>
                    </Control>
                    <Control
                        label={(
                            <>
                                Resolution
                                <InfoBadge label="About the sphere resolution">
                                    <p>
                                        Geodesic frequency ν of the hex tiling (10ν² + 2 cells — hexagons
                                        plus the 12 pentagons every hexagonal tiling of a sphere must
                                        contain). Auto is the engine&apos;s choice and the
                                        tiling the AI Copilot reads: the finest ν whose cells still
                                        average at least 12 displacements, the guard against reading
                                        Poisson noise as structure. Its label shows the ν it chose for
                                        this site.
                                    </p>
                                    <p>
                                        A tiling names a direction only to within its cell — at ν = 1
                                        (twelve cells) up to 37.4° off, at ν = 2 about 20°; the peak
                                        readout gives its cell&apos;s radius as ± R°. A finer ν places a
                                        lobe more closely, with fewer atoms per cell.
                                    </p>
                                </InfoBadge>
                            </>
                        )}
                    >
                        <select className="ui-select"
                            value={frequency}
                            onChange={(event) => setFrequency(event.target.value === 'auto' ? 'auto' : Number(event.target.value))}
                            aria-label="Sphere resolution"
                        >
                            {FREQUENCY_OPTIONS.map((option) => (
                                <option key={option} value={option}>
                                    {option === 'auto' ? autoLabel : `ν=${option} (${cellCountOf(option)})`}
                                </option>
                            ))}
                        </select>
                    </Control>
                    {sites?.reconstructed && (
                        <Control
                            label={(
                                <>
                                    Cluster
                                    <InfoBadge label="About site clustering">
                                        <p>
                                            This file carries no reference-site or cell columns, so sites
                                            are rebuilt by folding every atom into one unit cell and
                                            grouping atoms of the same element within this distance. Each
                                            site should gather one copy per supercell image; the count
                                            beside a site (e.g. 27/27) is its members against that
                                            expected number. Raise the distance to merge over-split sites,
                                            lower it to separate ones that ran together.
                                        </p>
                                    </InfoBadge>
                                </>
                            )}
                            value={<>{clusterThreshold.toFixed(1)} Å</>}
                        >
                            <input
                                className="ui-range"
                                type="range" min="0.4" max="2.5" step="0.1"
                                value={clusterThreshold}
                                onChange={(event) => setClusterThreshold(Number(event.target.value))}
                                aria-label="Site clustering distance in Angstrom"
                            />
                        </Control>
                    )}
                </ControlGroup>

                {/* Weighting & cutoff */}
                <ControlGroup label="Weighting and cutoff">
                    <Control
                        label={(
                            <>
                                Weight
                                <InfoBadge label="About the weighting">
                                    <p>
                                        Count: every atom votes once — the orientation distribution
                                        proper. |Δr| / |Δr|²: longer displacements vote more; the |Δr|²
                                        map is the angular decomposition of the mean-square displacement
                                        the U tensor summarizes.
                                    </p>
                                </InfoBadge>
                            </>
                        )}
                    >
                        <select className="ui-select" value={weight} onChange={(event) => setWeight(event.target.value)} aria-label="Cell weighting">
                            {WEIGHT_OPTIONS.map((option) => (
                                <option key={option.value} value={option.value}>{option.label}</option>
                            ))}
                        </select>
                    </Control>
                    <Control
                        label={(
                            <>
                                Min |Δr|
                                <InfoBadge label="About the amplitude cutoff">
                                    <p>
                                        Drops the shortest displacements before binning (quantile of
                                        |Δr|). A near-zero displacement has a direction dominated by
                                        noise, which dilutes a real pattern toward uniform.
                                    </p>
                                </InfoBadge>
                            </>
                        )}
                        value={<>{Math.round(minQuantile * 100)}%</>}
                    >
                        <input
                            className="ui-range"
                            type="range" min="0" max="0.5" step="0.05"
                            value={minQuantile}
                            onChange={(event) => setMinQuantile(Number(event.target.value))}
                            aria-label="Minimum displacement quantile"
                        />
                    </Control>
                    <Control
                        label={(
                            <>
                                Smoothing
                                <InfoBadge label="About smoothing">
                                    <p>
                                        Neighbour-diffusion passes over the cell graph (mass-conserving):
                                        each pass shares part of every cell's count with its neighbours,
                                        trading resolution for a smoother, lower-noise map. It only
                                        redraws the map: the peak, its direction and σ, and every test
                                        under the sphere come from the raw, unsmoothed counts.
                                    </p>
                                </InfoBadge>
                            </>
                        )}
                        value={<>{smoothing}×</>}
                    >
                        <input
                            className="ui-range"
                            type="range" min="0" max="12" step="1"
                            value={smoothing}
                            onChange={(event) => setSmoothing(Number(event.target.value))}
                            aria-label="Neighbour smoothing passes"
                        />
                    </Control>
                </ControlGroup>

                {/* Appearance */}
                <ControlGroup label="Appearance">
                    <Control
                        label={(
                            <>
                                Amplitude height
                                <InfoBadge label="About the amplitude height">
                                    <p>
                                        Raises the sphere surface radially by each cell's mean |Δr|
                                        relative to the site average — directions where atoms move
                                        farther stick out, shorter ones dent in, like a relief map.
                                        Color (how often) and height (how far) then carry independent
                                        information. 0% keeps a perfect sphere.
                                    </p>
                                </InfoBadge>
                            </>
                        )}
                        value={<>{Math.round(relief * 100)}%</>}
                    >
                        <input
                            className="ui-range"
                            type="range" min="0" max="1" step="0.05"
                            value={relief}
                            onChange={(event) => setRelief(Number(event.target.value))}
                            aria-label="Amplitude height"
                        />
                    </Control>
                    <Control label="Colormap">
                        <select className="ui-select" value={colormap} onChange={(event) => setColormap(event.target.value)} aria-label="Sphere colormap">
                            {COLORMAP_NAMES.map((name) => <option key={name} value={name}>{name}</option>)}
                        </select>
                    </Control>
                    <Control
                        label={(
                            <>
                                Contrast
                                <InfoBadge label="About the contrast">
                                    <p>
                                        Stretches the color scale about the isotropic level (1× = chance):
                                        higher values push cells away from the isotropic tone so faint
                                        lobes and depletions stand out, lower values flatten toward it.
                                        The colorbar tracks the same transfer. Surface height and
                                        significance are unaffected.
                                    </p>
                                </InfoBadge>
                            </>
                        )}
                        value={<>{contrast.toFixed(1)}×</>}
                    >
                        <input
                            className="ui-range"
                            type="range" min="0.5" max="3" step="0.1"
                            value={contrast}
                            onChange={(event) => setContrast(Number(event.target.value))}
                            aria-label="Color contrast"
                        />
                    </Control>
                    <Switch
                        label="Cell borders"
                        checked={showOutline}
                        onChange={(event) => setShowOutline(event.target.checked)}
                        inputProps={{ 'aria-label': 'Show cell borders' }}
                    />
                    <Switch
                        label="Axes"
                        checked={showAxes}
                        onChange={(event) => setShowAxes(event.target.checked)}
                        inputProps={{ 'aria-label': 'Show axis rods' }}
                    />
                </ControlGroup>
            </ControlsBar>

            {noRun && <Hint>Open a run folder with an <code>.rmc6f</code> file.</Hint>}

            {/* Three equal-height panels: axis views : sphere : site picker = 3 : 6.5 : 6.5. */}
            <div className="orient-layout">
                <OrientationView
                    result={histogram.map}
                    readout={histogram.readout}
                    loading={histogram.loading}
                    error={histogram.error}
                    loadedRef={histogram.loadedRef}
                    selectedRef={selectedRef}
                    selectedEllipsoid={selectedEllipsoid}
                    unitCell={unitCell}
                    frame={frame}
                    onFrameChange={setFrame}
                    colormap={colormap}
                    contrast={contrast}
                    relief={relief}
                    showOutline={showOutline}
                    showAxes={showAxes}
                    loadingSites={loadingSites}
                    onSelectSite={setSelectedRef}
                    siteRefs={siteRefs}
                />
                <SiteStructurePanel
                    cardId="direction-sites"
                    sites={sites}
                    selectedRef={selectedRef}
                    onSelectSite={setSelectedRef}
                    selectedEllipsoid={selectedEllipsoid}
                    elementColors={elementColors}
                    loadingSites={loadingSites}
                />
            </div>

            <AppFooter tight />
        </Page>
    );
}
