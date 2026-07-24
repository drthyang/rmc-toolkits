// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Orientation page: the solid-angle distribution of a site's displacement
// directions, hex-binned on a Goldberg sphere. A standalone workspace page —
// the histogram is not a PCA product, so it does not live on the PCA Ellipsoid
// page; the two share only the site picker, provided by the useSiteCloud hook
// (same data plumbing, same worker parse cache).

import React, { useState } from 'react';
import { isStaticMode } from '../browserData';
import InfoBadge from './InfoBadge';
import OrientationView from './OrientationView';
import useSiteCloud from '../useSiteCloud';
import './PcaKdePage.css';

const numberFormat = (value, digits = 4) =>
    Number.isFinite(value) ? value.toFixed(digits) : '—';

const DEFAULT_CLUSTER_THRESHOLD = 1.5;

export default function OrientationPage({ directory, localRun }) {
    // Fold-and-cluster distance (Å), used only when the loaded file has no
    // reference-site/cell columns and its sites must be reconstructed.
    const [clusterThreshold, setClusterThreshold] = useState(DEFAULT_CLUSTER_THRESHOLD);

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

    return (
        <div className="pca-page">
            <div className="pca-controls">
                <div className="control-group" role="group" aria-label="Site">
                    <label className="control">
                        <span className="control-name">
                            Site
                            <InfoBadge label="About the site picker">
                                <p>
                                    Each reference site (an RMCProfile reference number) is one
                                    crystallographic position. Only the <em>directions</em> of its
                                    per-atom displacements are analysed here — the amplitude enters
                                    solely through the optional weighting and the relief.
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
            </div>

            {noRun && (
                <p className="pca-hint">Open a run folder (with an <code>.rmc6f</code> file) to view displacement orientations.</p>
            )}
            {sitesError && <p className="pca-error-banner">{sitesError}</p>}

            <div className="pca-panel pca-viewport orient-page-panel">
                <h3>
                    <span className="panel-title-label">
                        {selectedEllipsoid
                            ? `${selectedEllipsoid.element} site #${selectedEllipsoid.referenceNumber} — displacement directions`
                            : 'Displacement directions'}
                    </span>
                    <span className="panel-title-actions">
                        {selectedEllipsoid && (
                            <span className="panel-title-count">{selectedEllipsoid.count.toLocaleString()} atoms</span>
                        )}
                        {loadingSites && <span className="panel-title-count">Loading sites…</span>}
                    </span>
                </h3>
                <OrientationView
                    requestPca={requestPca}
                    ready={ready}
                    selectedRef={selectedRef}
                    selectedEllipsoid={selectedEllipsoid}
                    clusterThreshold={clusterThreshold}
                    unitCell={unitCell}
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
