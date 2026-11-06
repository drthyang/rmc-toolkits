// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// The Bond Geometry page's hover highlight: pointing at a bin of the angle
// histogram picks out, in the folded cell, the A–B–C triplets the histogram
// counts in that bin. Plain numbers, no Three.js, so it is unit-tested.
//
// Two halves. siteAngleTriplets runs in a worker (angleTripletsWorker.js):
// it walks the engine's own triplet loop (workers/triplets.js
// visitAngleTriplets — the windows, periodic images and central atoms of the
// histogram, bin by bin) and folds each supercell triplet onto the average
// sites the folded cell draws: the central atom's site in the cell, each end
// at the image of its site that the bond reaches. Many supercell triplets
// fold onto one site triplet; each keeps its count. highlightForBins runs on
// hover, or as the angle range on the plot is dragged: the site triplets of
// those bins (whole bins: the range's edges are bin edges) as tube
// segments, at most MAX_HIGHLIGHT_TRIPLETS of them, nearest the cell centre.
//
// A large cell draws only the bonds nearest its centre (bondNetwork.js), and
// the highlight keeps to that region: a triplet is kept when both its tubes'
// midpoints lie within the drawn bonds' radius. Every angle of the bin is
// counted (`total` equals the histogram's count); `shown` is how many of
// them the drawn site triplets carry.

import { visitAngleTriplets } from './workers/triplets.js';

// Site triplets drawn per bin: two tubes each, so the highlight never draws
// more than the bond network's own cap (MAX_DRAWN_BONDS = 1,500).
export const MAX_HIGHLIGHT_TRIPLETS = 750;
// Site triplets the worker stores per bin; the angles of any further ones
// are counted, not stored, so memory stays bounded on any cell.
export const MAX_BIN_TRIPLETS = 20000;

const fracToCartesian = (unit, f) => [0, 1, 2].map(
    (axis) => f[0] * unit[0][axis] + f[1] * unit[1][axis] + f[2] * unit[2][axis]
);

const byDistance = (p, q) => p.d - q.d || p.id - q.id;

const distanceSq = (p, q) => (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 + (p[2] - q[2]) ** 2;

export const unitCellOf = (latticeVectors, supercell) => latticeVectors.map(
    (row, i) => row.map((value) => value / (supercell[i] || 1))
);

/**
 * Every bin's triplets folded onto the average sites.
 *
 * @param {{ fractional: number[][], elements: string[], sites: Array<number|null> }} atoms
 *   the configuration (supercell fractions) and each atom's site (reference number)
 * @param {number[][]} latticeVectors supercell vectors (rows, Å)
 * @param {number[]} supercell        unit cells along a, b, c
 * @param {Map<number, number[]>} siteFractional site → in-cell fractional position
 * @param {object} options the engine options (triplets.js tripletOptionsFromRequest)
 * @param {{ center: number[], radiusSq: number|null }|null} region the drawn
 *   region (Cartesian Å; bondNetwork's center / radiusSq), null: the whole cell
 * @returns {{ nbins, binWidth, triplet, sharedEnds, bins: Array<{ total, triplets }> }}
 *   per bin: `total` angles (the histogram's count) and its site triplets
 *   { id, apex, ab, bc, count, d }, nearest the centre first (d: the apex
 *   site's squared distance from it; id: the same site triplet in every
 *   bin). ab / bc are
 *   [end site, image a, image b, image c]: the end at its site's position
 *   plus that whole-cell image; the apex sits at its own site in the cell.
 */
export const siteAngleTriplets = (atoms, latticeVectors, supercell, siteFractional, options,
    region = null, { maxPerBin = MAX_BIN_TRIPLETS } = {}) => {
    const unit = unitCellOf(latticeVectors, supercell);
    const center = region?.center ?? fracToCartesian(unit, [0.5, 0.5, 0.5]);
    const radiusSq = region?.radiusSq ?? null;
    let wrapped = null;
    let sharedEnds = false;

    // One site bond per distinct (apex site, end site, image), folded once
    // per supercell bond, not once per angle.
    const bondIds = new Map();
    const siteBonds = [];
    const siteBondOf = (apexRow, bond) => {
        if (bond.siteBond !== undefined) return bond.siteBond;
        const apex = atoms.sites[apexRow];
        const end = atoms.sites[bond.row];
        const fa = siteFractional.get(apex);
        const fe = siteFractional.get(end);
        if (!fa || !fe) {
            bond.siteBond = null;
            return null;
        }
        // The bond in unit-cell fractions, then the whole-cell image of the
        // end site that it reaches from the apex site.
        const from = wrapped[apexRow];
        const to = wrapped[bond.row];
        const image = [bond.ix, bond.iy, bond.iz].map((shift, axis) => Math.round(
            (to[axis] - from[axis] + shift) * supercell[axis] - (fe[axis] - fa[axis])
        ));
        const key = `${apex}|${end}|${image.join(',')}`;
        let id = bondIds.get(key);
        if (id === undefined) {
            id = siteBonds.length;
            bondIds.set(key, id);
            const mid = fracToCartesian(unit, [0, 1, 2].map((axis) => (fa[axis] + fe[axis] + image[axis]) / 2));
            siteBonds.push({ apex, end, image, inside: radiusSq == null || distanceSq(mid, center) <= radiusSq });
        }
        bond.siteBond = id;
        return id;
    };

    // Per bin: A–B site bond → B–C site bond → angle count.
    const binMaps = [];
    const sizes = [];
    const totals = [];
    const info = visitAngleTriplets(atoms.fractional, atoms.elements, latticeVectors, options, (bin, apexRow, one, two) => {
        totals[bin] += 1;
        const ab = siteBondOf(apexRow, one);
        const bc = siteBondOf(apexRow, two);
        if (ab === null || bc === null || !siteBonds[ab].inside || !siteBonds[bc].inside) return;
        // One shared window: the two bonds are interchangeable, so the
        // unordered pair is the triplet.
        const [x, y] = sharedEnds && bc < ab ? [bc, ab] : [ab, bc];
        let inner = binMaps[bin].get(x);
        const count = inner?.get(y);
        if (count !== undefined) {
            inner.set(y, count + 1);
            return;
        }
        if (sizes[bin] >= maxPerBin) return;
        if (!inner) {
            inner = new Map();
            binMaps[bin].set(x, inner);
        }
        inner.set(y, 1);
        sizes[bin] += 1;
    }, (start) => {
        wrapped = start.wrapped;
        sharedEnds = start.sharedEnds;
        for (let bin = 0; bin < start.nbins; bin += 1) {
            binMaps.push(new Map());
            sizes.push(0);
            totals.push(0);
        }
    });

    // One id per site triplet across the bins, in the order first found.
    const tripletIds = new Map();
    const tripletId = (abId, bcId) => {
        const key = abId * (siteBonds.length + 1) + bcId;
        let id = tripletIds.get(key);
        if (id === undefined) {
            id = tripletIds.size;
            tripletIds.set(key, id);
        }
        return id;
    };
    const apexDistance = new Map();
    const distanceOf = (apex) => {
        if (!apexDistance.has(apex)) apexDistance.set(apex, distanceSq(fracToCartesian(unit, siteFractional.get(apex)), center));
        return apexDistance.get(apex);
    };
    const bins = binMaps.map((map, bin) => {
        const found = [];
        map.forEach((inner, abId) => inner.forEach((count, bcId) => {
            const ab = siteBonds[abId];
            const bc = siteBonds[bcId];
            found.push({
                id: tripletId(abId, bcId), apex: ab.apex, ab: [ab.end, ...ab.image], bc: [bc.end, ...bc.image], count, d: distanceOf(ab.apex)
            });
        }));
        // Nearest the centre first; ties in the order first found.
        found.sort(byDistance);
        return { total: totals[bin], triplets: found };
    });
    return { nbins: info.nbins, binWidth: info.binWidth, triplet: info.triplet, sharedEnds: info.sharedEnds, bins };
};

/** The bin a hovered x (a bin centre, or anywhere in the bin) falls in. */
export const binIndexOf = (x, binWidth, nbins) => {
    if (!Number.isFinite(x) || !(binWidth > 0)) return null;
    const index = Math.floor(x / binWidth);
    return Math.min(Math.max(index, 0), nbins - 1);
};

// A bin edge as the legend prints it: whole degrees, else up to two decimals.
const degrees = (value) => String(Number(value.toFixed(2)));

/** An angle range as one plain line: "85–95°". */
export const rangeText = (lo, hi) => `${degrees(lo)}–${degrees(hi)}°`;

/**
 * The legend line of a highlight: "830 angles at 88–89°" for one hovered bin,
 * "2,310 angles in 85–95°" for a range; "12 of 830 … shown" when the drawn
 * triplets carry fewer than all of them.
 */
export const highlightLabel = ({ shown, total, lo, hi, range = false }) => {
    const where = `${range ? 'in' : 'at'} ${rangeText(lo, hi)}`;
    const noun = total === 1 ? 'angle' : 'angles';
    return shown < total
        ? `${shown.toLocaleString('en-US')} of ${total.toLocaleString('en-US')} ${noun} ${where} shown`
        : `${total.toLocaleString('en-US')} ${noun} ${where}`;
};

/**
 * The highlight of bins first..last (inclusive): their site triplets as tube
 * segments (six numbers each, apex → end, Cartesian Å), the A–B and B–C
 * bonds apart (with one shared window, all in `ab`), at most `max` triplets
 * — those nearest the cell centre. A site triplet in several of the bins is
 * drawn once, with their counts summed. `total` is the bins' angle count
 * (the histogram's); `shown` counts the angles the drawn triplets carry.
 *
 * @param {object} data     siteAngleTriplets' result
 * @param {number} first    first bin
 * @param {number} last     last bin (>= first)
 * @param {Map<number, number[]>} siteFractional site → in-cell fractional position
 * @param {number[][]} unit unit-cell vectors (rows, Å)
 */
export const highlightForBins = (data, first, last, siteFractional, unit, { max = MAX_HIGHLIGHT_TRIPLETS, range = false } = {}) => {
    if (!data?.bins || !(first >= 0) || !(last >= first) || last >= data.bins.length) return null;
    let total = 0;
    let merged;
    if (first === last) {
        total = data.bins[first].total;
        merged = data.bins[first].triplets;
    } else {
        const byId = new Map();
        for (let bin = first; bin <= last; bin += 1) {
            total += data.bins[bin].total;
            data.bins[bin].triplets.forEach((triplet) => {
                const seen = byId.get(triplet.id);
                if (seen) seen.count += triplet.count;
                else byId.set(triplet.id, { ...triplet });
            });
        }
        merged = [...byId.values()].sort(byDistance);
    }
    const drawn = merged.slice(0, max);
    const ab = new Float64Array(drawn.length * 6 * (data.sharedEnds ? 2 : 1));
    const bc = new Float64Array(data.sharedEnds ? 0 : drawn.length * 6);
    let shown = 0;
    let abAt = 0;
    let bcAt = 0;
    const put = (out, at, apex, [end, ...image]) => {
        const fe = siteFractional.get(end);
        out.set(apex, at);
        out.set(fracToCartesian(unit, [fe[0] + image[0], fe[1] + image[1], fe[2] + image[2]]), at + 3);
    };
    drawn.forEach((triplet) => {
        const apex = fracToCartesian(unit, siteFractional.get(triplet.apex));
        put(ab, abAt, apex, triplet.ab);
        abAt += 6;
        if (data.sharedEnds) {
            put(ab, abAt, apex, triplet.bc);
            abAt += 6;
        } else {
            put(bc, bcAt, apex, triplet.bc);
            bcAt += 6;
        }
        shown += triplet.count;
    });
    return {
        first,
        last,
        range,
        lo: first * data.binWidth,
        hi: (last + 1) * data.binWidth,
        total,
        shown,
        ab,
        bc
    };
};

/** One bin's highlight (the hovered one): highlightForBins of that bin alone. */
export const highlightForBin = (data, bin, siteFractional, unit, options) => (
    highlightForBins(data, bin, bin, siteFractional, unit, options)
);

/**
 * An angle range [lo, hi] (degrees) → its bins. The range's edges are bin
 * edges (the plot snaps them), so each bin is wholly in or out; an edge off
 * the grid rounds to the nearest edge. Null when it covers no bin.
 */
export const binsOfRange = (range, binWidth, nbins) => {
    if (!range || !(binWidth > 0)) return null;
    const first = Math.max(0, Math.round(range[0] / binWidth));
    const last = Math.min(nbins - 1, Math.round(range[1] / binWidth) - 1);
    return last >= first ? [first, last] : null;
};
