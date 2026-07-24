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

import React, { useMemo, useState } from 'react';
import { isStaticMode } from '../browserData';
import { buildElementColors } from '../atomColors';
import { COLORMAP_NAMES } from '../colormaps';
import InfoBadge from './InfoBadge';
import OrientationView from './OrientationView';
import SiteStructurePanel from './SiteStructurePanel';
import useSiteCloud from '../useSiteCloud';
import './PcaKdePage.css';

const numberFormat = (value, digits = 4) =>
    Number.isFinite(value) ? value.toFixed(digits) : '—';

const DEFAULT_CLUSTER_THRESHOLD = 1.5;

// Manual resolution choices (geodesic frequency ν → 10ν²+2 cells). 'auto' asks
// the engine for recommended_frequency, the ~12-points-per-cell over-binning
// guard.
const FREQUENCY_OPTIONS = ['auto', 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24];

const WEIGHT_OPTIONS = [
    { value: 'count', label: 'Count' },
    { value: 'amplitude', label: '|Δr|' },
    { value: 'amplitude2', label: '|Δr|²' }
];

export default function OrientationPage({ directory, localRun }) {
    // Fold-and-cluster distance (Å), used only when the loaded file has no
    // reference-site/cell columns and its sites must be reconstructed.
    const [clusterThreshold, setClusterThreshold] = useState(DEFAULT_CLUSTER_THRESHOLD);

    // Histogram + display options (owned here, in the top controls bar).
    // Default ν=10 (1002 cells) + 2× smoothing gives a legible map on a
    // typical ~1000-copy site out of the box.
    const [frequency, setFrequency] = useState(10);
    const [weight, setWeight] = useState('count');
    const [frame, setFrame] = useState('cartesian');
    const [smoothing, setSmoothing] = useState(2);
    const [minQuantile, setMinQuantile] = useState(0);
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
    } = useSiteCloud({ directory, localRun, clusterThreshold });

    const staticMode = isStaticMode();
    const noRun = staticMode && !localFile;

    const elementColors = useMemo(
        () => buildElementColors(sites?.elements ?? []),
        [sites]
    );

    return (
        <div className="pca-page">
            <div className="pca-controls">
                {/* Site & resolution */}
                <div className="control-group" role="group" aria-label="Site and resolution">
                    <label className="control">
                        <span className="control-name">
                            Site
                            <InfoBadge label="About the site picker">
                                <p>
                                    Each reference site (an RMCProfile reference number) is one
                                    crystallographic position. Only the <em>directions</em> of its
                                    per-atom displacements are analysed here — the amplitude enters
                                    solely through the optional weighting and the amplitude height.
                                </p>
                            </InfoBadge>
                        </span>
                        <select
                            value={selectedRef ?? ''}
                            onChange={(event) => setSelectedRef(Number(event.target.value))}
                            disabled={!sites}
                            aria-label="Site"
                        >
                            {sites?.sites.map((site) => (
                                <option key={site.referenceNumber} value={site.referenceNumber}>
                                    {`#${site.referenceNumber} ${site.element} — U=${numberFormat(site.uIso, 4)} Å²`}
                                    {site.copiesPerCell ? ` (${site.count}/${site.copiesPerCell})` : ''}
                                </option>
                            ))}
                        </select>
                    </label>
                    <label className="control">
                        <span className="control-name">
                            Resolution
                            <InfoBadge label="About the sphere resolution">
                                <p>
                                    Geodesic frequency ν of the hex tiling (10ν² + 2 cells — hexagons
                                    plus the 12 pentagons every hexagonal tiling of a sphere must
                                    contain). Auto targets ~12 displacements per cell, the guard
                                    against reading Poisson noise as structure.
                                </p>
                            </InfoBadge>
                        </span>
                        <select
                            value={frequency}
                            onChange={(event) => setFrequency(event.target.value === 'auto' ? 'auto' : Number(event.target.value))}
                            aria-label="Sphere resolution"
                        >
                            {FREQUENCY_OPTIONS.map((option) => (
                                <option key={option} value={option}>
                                    {option === 'auto' ? 'Auto' : `ν=${option} (${10 * option * option + 2})`}
                                </option>
                            ))}
                        </select>
                    </label>
                    {sites?.reconstructed && (
                        <label className="control">
                            <span className="control-name">
                                Cluster
                                <InfoBadge label="About site clustering">
                                    <p>
                                        This file carries no reference-site or cell columns, so sites
                                        are rebuilt by folding every atom into one unit cell and
                                        grouping atoms of the same element within this distance.
                                    </p>
                                </InfoBadge>
                            </span>
                            <input
                                type="range" min="0.4" max="2.5" step="0.1"
                                value={clusterThreshold}
                                onChange={(event) => setClusterThreshold(Number(event.target.value))}
                                aria-label="Site clustering distance in Angstrom"
                            />
                            <span className="control-value">{clusterThreshold.toFixed(1)} Å</span>
                        </label>
                    )}
                </div>

                {/* Weighting & cutoff */}
                <div className="control-group" role="group" aria-label="Weighting and cutoff">
                    <label className="control">
                        <span className="control-name">
                            Weight
                            <InfoBadge label="About the weighting">
                                <p>
                                    Count: every atom votes once — the orientation distribution
                                    proper. |Δr| / |Δr|²: longer displacements vote more; the |Δr|²
                                    map is the angular decomposition of the mean-square displacement
                                    the U tensor summarizes.
                                </p>
                            </InfoBadge>
                        </span>
                        <select value={weight} onChange={(event) => setWeight(event.target.value)} aria-label="Cell weighting">
                            {WEIGHT_OPTIONS.map((option) => (
                                <option key={option.value} value={option.value}>{option.label}</option>
                            ))}
                        </select>
                    </label>
                    <label className="control">
                        <span className="control-name">
                            Min |Δr|
                            <InfoBadge label="About the amplitude cutoff">
                                <p>
                                    Drops the shortest displacements before binning (quantile of
                                    |Δr|). A near-zero displacement has a direction dominated by
                                    noise, which dilutes a real pattern toward uniform.
                                </p>
                            </InfoBadge>
                        </span>
                        <input
                            type="range" min="0" max="0.5" step="0.05"
                            value={minQuantile}
                            onChange={(event) => setMinQuantile(Number(event.target.value))}
                            aria-label="Minimum displacement quantile"
                        />
                        <span className="control-value">{Math.round(minQuantile * 100)}%</span>
                    </label>
                    <label className="control">
                        <span className="control-name">
                            Smoothing
                            <InfoBadge label="About smoothing">
                                <p>
                                    Neighbour-diffusion passes over the cell graph (mass-conserving):
                                    each pass shares part of every cell's count with its neighbours,
                                    trading resolution for a smoother, lower-noise map. The z-scores
                                    and significance always come from the raw, unsmoothed counts.
                                </p>
                            </InfoBadge>
                        </span>
                        <input
                            type="range" min="0" max="12" step="1"
                            value={smoothing}
                            onChange={(event) => setSmoothing(Number(event.target.value))}
                            aria-label="Neighbour smoothing passes"
                        />
                        <span className="control-value">{smoothing}×</span>
                    </label>
                </div>

                {/* Appearance */}
                <div className="control-group" role="group" aria-label="Appearance">
                    <label className="control">
                        <span className="control-name">
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
                        </span>
                        <input
                            type="range" min="0" max="1" step="0.05"
                            value={relief}
                            onChange={(event) => setRelief(Number(event.target.value))}
                            aria-label="Amplitude height"
                        />
                        <span className="control-value">{Math.round(relief * 100)}%</span>
                    </label>
                    <label className="control">
                        <span className="control-name">Colormap</span>
                        <select value={colormap} onChange={(event) => setColormap(event.target.value)} aria-label="Sphere colormap">
                            {COLORMAP_NAMES.map((name) => <option key={name} value={name}>{name}</option>)}
                        </select>
                    </label>
                    <label className="control">
                        <span className="control-name">
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
                        </span>
                        <input
                            type="range" min="0.5" max="3" step="0.1"
                            value={contrast}
                            onChange={(event) => setContrast(Number(event.target.value))}
                            aria-label="Color contrast"
                        />
                        <span className="control-value">{contrast.toFixed(1)}×</span>
                    </label>
                    <label className="control switch">
                        <span className="control-name">Cell borders</span>
                        <input type="checkbox" checked={showOutline} onChange={(event) => setShowOutline(event.target.checked)} aria-label="Show cell borders" />
                        <i className="switch-track" aria-hidden="true" />
                    </label>
                    <label className="control switch">
                        <span className="control-name">Axes</span>
                        <input type="checkbox" checked={showAxes} onChange={(event) => setShowAxes(event.target.checked)} aria-label="Show axis rods" />
                        <i className="switch-track" aria-hidden="true" />
                    </label>
                </div>
            </div>

            {noRun && (
                <p className="pca-hint">Open a run folder (with an <code>.rmc6f</code> file) to view displacement orientations.</p>
            )}
            {sitesError && <p className="pca-error-banner">{sitesError}</p>}

            {/* Three equal-height panels: axis views : sphere : site picker = 3 : 6.5 : 6.5. */}
            <div className="orient-layout">
                <OrientationView
                    requestPca={requestPca}
                    ready={ready}
                    selectedRef={selectedRef}
                    selectedEllipsoid={selectedEllipsoid}
                    clusterThreshold={clusterThreshold}
                    unitCell={unitCell}
                    frequency={frequency}
                    weight={weight}
                    frame={frame}
                    onFrameChange={setFrame}
                    smoothing={smoothing}
                    minQuantile={minQuantile}
                    colormap={colormap}
                    contrast={contrast}
                    relief={relief}
                    showOutline={showOutline}
                    showAxes={showAxes}
                />
                <SiteStructurePanel
                    sites={sites}
                    selectedRef={selectedRef}
                    onSelectSite={setSelectedRef}
                    selectedEllipsoid={selectedEllipsoid}
                    elementColors={elementColors}
                    loadingSites={loadingSites}
                />
            </div>

            <footer className="app-footer">
                &copy; 2026 Tsung-Han Yang &middot;{' '}
                <a
                    href="https://github.com/drthyang/rmc-toolkits/blob/main/LICENSE"
                    target="_blank"
                    rel="noreferrer"
                >
                    AGPLv3
                </a>
                {' '}&middot;{' '}
                <a
                    href="https://github.com/drthyang/rmc-toolkits#readme"
                    target="_blank"
                    rel="noreferrer"
                >
                    About & documentation
                </a>
            </footer>
        </div>
    );
}
