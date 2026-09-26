// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The label of one PCA/orientation site. A mixed-occupancy site (one
// reference number, several species) reads as its composition, majority
// first and ties by name, e.g. Ga0.75In0.25; a pure site as its element.
// Shared by the PCA Ellipsoid and Displacement Directions pages, so the two
// name the same site the same way (headings, pickers, saved files).
export const siteLabel = (site) => {
    if (!site?.mixed || !site.elementCounts) return site?.element ?? '';
    const entries = Object.entries(site.elementCounts);
    const total = entries.reduce((sum, [, count]) => sum + count, 0) || 1;
    return entries
        .sort(([nameA, a], [nameB, b]) => b - a || (nameA < nameB ? -1 : nameA > nameB ? 1 : 0))
        .map(([name, count]) => `${name}${(count / total).toFixed(2)}`)
        .join('');
};
