// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Workbench cards the AI Copilot can point at: the host's half of the
// llm module's `cards` prop. Each card in the pages carries the same id as
// data-copilot-card; the entry says where it lives (page key), what it shows
// (for the model), and how the module may use it: `kind` / `file` tie a card
// to the run data a tool reads (a plot kind, 'r_value', 'symmetry', 'model'),
// `marks` that its plot can carry x marks (Dashboard.jsx → InteractivePlot),
// `view3d` that it is a 3D view orient_3d_view can turn (viewer3d.js), and
// `sites` that it shows the unit cell's sites.

// The pages an AI Copilot step can open ("Show in Bond Geometry") and a card
// can live on, by page key, with their tab labels: CopilotChat's `pages` and
// the view tools' page names (App.jsx). Auto StoG is not one of them.
export const COPILOT_PAGES = Object.freeze({
    dashboard: 'Dashboard',
    structure: 'Atomic Density',
    geometry: 'Bond Geometry',
    ellipsoids: 'PCA Ellipsoid',
    orientation: 'Displacement Directions'
});

const PLOT_ABOUT = {
    pdf_partials: 'the partial g(r) of every element pair (x: r in Å)',
    xpdf: 'the fitted X-ray G(r), model against data (x: r in Å)',
    npdf: 'the fitted neutron G(r), model against data (x: r in Å)',
    xray_sq: 'the fitted X-ray F(Q), model against data (x: Q in Å⁻¹)',
    neutron_sq: 'the fitted neutron S(Q), model against data (x: Q in Å⁻¹)',
    bragg: 'the fitted Bragg profile (x: TOF or d)',
    exafs_q: 'the fitted EXAFS in k space',
    exafs_r: 'the fitted EXAFS in r space (x: r in Å)',
    stog: 'a STOG data file (x: r in Å or Q in Å⁻¹)'
};

export const PAGE_CARDS = [
    { id: 'model-info', page: 'dashboard', kind: 'model', title: 'Model information', about: 'cell, supercell, atoms per element, move counters' },
    { id: 'detected-sg', page: 'dashboard', kind: 'symmetry', title: 'Detected SG', about: 'space group of the average structure and its tolerance ladder' },
    { id: 'kde-slice', page: 'structure', title: 'KDE slice', about: 'density map of one element in a slab through the folded cell' },
    { id: 'slab-in-cell', page: 'structure', title: 'Slab in cell', about: 'where that slab sits in the unit cell' },
    { id: 'density-cell', page: 'structure', title: 'Folded unit cell (3D)', view3d: true, about: 'every atom folded into one unit cell, by element' },
    { id: 'bond-angles', page: 'geometry', title: 'Bond-angle distribution', about: 'the A–B–C angle histogram computed on the page' },
    { id: 'bond-cell', page: 'geometry', title: 'Folded cell with bonds (3D)', view3d: true, about: 'the folded atom cloud with the detected bonds' },
    { id: 'bond-windows', page: 'geometry', title: 'Partial g(r) windows', about: 'the bond windows on the partial g(r)' },
    { id: 'pca-density', page: 'ellipsoids', title: 'Displacement density (3D)', view3d: true, about: 'the selected site\'s displacement cloud: KDE shell or isosurface and thermal ellipsoid' },
    { id: 'pca-stats', page: 'ellipsoids', title: 'Site statistics', about: 'U_iso, principal RMS amplitudes, anisotropy, kurtosis and the axes in the crystal' },
    { id: 'pca-sites', page: 'ellipsoids', title: 'Site ellipsoids (3D)', view3d: true, sites: true, about: 'every site of the unit cell as its thermal ellipsoid; the selected site in orange' },
    { id: 'direction-sphere', page: 'orientation', title: 'Displacement-direction sphere (3D)', view3d: true, about: 'where the selected site\'s displacements point, as enhancement over isotropic on the sphere' },
    { id: 'direction-views', page: 'orientation', title: 'Axis views', about: 'that sphere seen down a, b and c' },
    { id: 'direction-sites', page: 'orientation', title: 'Site ellipsoids (3D)', view3d: true, sites: true, about: 'the unit-cell site picker of the Displacement Directions page' }
];

/**
 * The cards for the open run: one per plot on the Dashboard (`plot:<file>`),
 * the χ² history when the run has logs, then the fixed cards of every page
 * (Detected SG only with a symmetry analysis: a browser-loaded run).
 */
export const buildCopilotCards = ({ plotFiles = [], rValueFile = null, symmetry = null } = {}) => {
    const plots = (plotFiles || [])
        .filter((file) => file?.name && file.plotKind && file.plotKind !== 'r_value')
        .map((file) => ({
            id: `plot:${file.name}`,
            page: 'dashboard',
            kind: file.plotKind,
            file: file.name,
            marks: true,
            title: file.plotData?.title && file.plotData.title !== file.name
                ? `${file.plotData.title} (${file.name})`
                : file.name,
            about: PLOT_ABOUT[file.plotKind] || 'a fitted dataset'
        }));
    const history = rValueFile
        ? [{
            id: 'chi2-history',
            page: 'dashboard',
            kind: 'r_value',
            file: rValueFile.name,
            marks: true,
            title: 'χ² history',
            about: 'ln χ² of the last .log column against time, with the convergence badge'
        }]
        : [];
    const [model, detected, ...rest] = PAGE_CARDS;
    return [model, ...(symmetry ? [detected] : []), ...history, ...plots, ...rest];
};
