// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// Bond-angle (triplet) distribution engine — static-mode port of
// rmc_toolkits/triplets.py (the source of truth; keep the two in sync).
//
// Given a triplet (A, B, C) with B the central atom and an inclusive length
// window for each of the A–B and B–C bonds, finds every A–B–C triplet in a
// periodic configuration and histograms the angle at B. Neighbours come from
// a linked-cell search that carries explicit periodic-image shifts instead of
// assuming minimum image, so it is exact for any cell shape — triclinic
// included — and for boxes smaller than the cutoff, where several images of
// one atom are genuine distinct neighbours. A pair at exactly zero length
// (bitwise-coincident atoms under rmin = 0) is never a bond.
//
// Angle counting mirrors the Python engine: equivalent ends (same element,
// same window) count each unordered pair of distinct bonds once; distinct
// windows count ordered (1→2, 2→3) assignments minus the combinations where
// both bonds reach the same atom in the same periodic image.
//
// The summary payload (bondAngleSummary) is the contract shared with the
// Flask /api/triplets route — camelCase keys, plain arrays — so the page
// renders identical numbers in both runtimes. `sinCorrected` divides each
// bin's count fraction by the exact isotropic fraction (cosθ₁ − cosθ₂)/2,
// flat 1.0 for random directions and finite at 0°/180°.
//
// Cross-engine caveat (also noted in triplets.py): libm acos can differ from
// Math.acos by 1 ulp, so a bitwise-ideal geometry whose cosine lands exactly
// on a bin edge (e.g. cos = 0.5 in an undisplaced average configuration) can
// shift one count into the neighbouring bin between the two engines. Element
// symbols are ASCII, where the capitalize below matches str.capitalize().

export const MAX_CELLS_PER_AXIS = 64;
export const LENGTH_BINS = 40;
const REACH_HEADROOM = 1e-9;

const capitalize = (symbol) => {
  const text = String(symbol).trim();
  return text.charAt(0).toUpperCase() + text.slice(1).toLowerCase();
};

const validateWindow = (name, window) => {
  if (!Array.isArray(window) || window.length !== 2) {
    throw new Error(`${name} must be a (rmin, rmax) pair`);
  }
  // Number(null) and Number('') coerce to 0 — reject them explicitly so a
  // missing bound errors like Python's float() instead of becoming rmin = 0.
  if (window.some((value) => value == null || value === '')) {
    throw new Error(`${name} bounds must be finite, got (${window[0]}, ${window[1]})`);
  }
  const rmin = Number(window[0]);
  const rmax = Number(window[1]);
  if (!Number.isFinite(rmin) || !Number.isFinite(rmax)) {
    throw new Error(`${name} bounds must be finite, got (${window[0]}, ${window[1]})`);
  }
  if (rmin < 0 || rmax <= rmin) {
    throw new Error(`${name} needs 0 <= rmin < rmax, got (${rmin}, ${rmax})`);
  }
  return [rmin, rmax];
};

// Bin index exactly as numpy.histogram assigns it for uniform bins: truncate
// (x - lo) * n / (hi - lo), make the right edge of the last bin inclusive,
// then correct against the linspace edge values — float rounding can put the
// truncation one bin off for values sitting exactly on an interior edge.
const histogramIndex = (x, lo, hi, n) => {
  const norm = n / (hi - lo);
  let index = Math.floor((x - lo) * norm);
  if (index >= n) index = n - 1;
  const step = (hi - lo) / n;
  const edge = (i) => (i === n ? hi : i * step + lo);
  if (x < edge(index)) index -= 1;
  else if (index !== n - 1 && x >= edge(index + 1)) index += 1;
  return index;
};

const cross = (a, b) => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0]
];

const norm = (v) => Math.hypot(v[0], v[1], v[2]);

// Perpendicular width of the box along each lattice direction:
// w_i = V / |a_j × a_k| (mirrors _perpendicular_widths).
const perpendicularWidths = (lattice) => {
  const [a, b, c] = lattice;
  const crosses = [cross(b, c), cross(c, a), cross(a, b)];
  const volume = Math.abs(
    a[0] * crosses[0][0] + a[1] * crosses[0][1] + a[2] * crosses[0][2]
  );
  const areas = crosses.map(norm);
  if (!(volume > 0) || areas.some((area) => !(area > 0))) {
    throw new Error('lattice vectors are singular (zero cell volume)');
  }
  return areas.map((area) => volume / area);
};

// All (center, candidate image) bonds whose length falls inside one window,
// grouped per center (mirrors _neighbor_bonds; grouping replaces the numpy
// sort-by-center step). Each bond: { row, ix, iy, iz, vx, vy, vz, length }.
const neighborBonds = (wrapped, centerRows, candidateRows, lattice, rmin, rmax) => {
  const widths = perpendicularWidths(lattice);
  const cells = widths.map((width) =>
    Math.min(Math.max(1, Math.floor(width / rmax)), MAX_CELLS_PER_AXIS)
  );
  const reach = widths.map((width, axis) =>
    Math.ceil((rmax * cells[axis]) / width + REACH_HEADROOM)
  );
  const cellOf = (value, axis) =>
    Math.min(Math.floor(value * cells[axis]), cells[axis] - 1);

  const nCells = cells[0] * cells[1] * cells[2];
  const buckets = new Array(nCells);
  for (const row of candidateRows) {
    const flat =
      (cellOf(wrapped[row][0], 0) * cells[1] + cellOf(wrapped[row][1], 1)) * cells[2] +
      cellOf(wrapped[row][2], 2);
    (buckets[flat] ??= []).push(row);
  }

  const rminSq = rmin * rmin;
  const rmaxSq = rmax * rmax;
  const perCenter = new Array(centerRows.length);
  for (let center = 0; center < centerRows.length; center += 1) {
    const centerRow = centerRows[center];
    const [fx, fy, fz] = wrapped[centerRow];
    const cx = cellOf(fx, 0);
    const cy = cellOf(fy, 1);
    const cz = cellOf(fz, 2);
    const bonds = [];
    for (let dx = -reach[0]; dx <= reach[0]; dx += 1) {
      const sx = cx + dx;
      const wx = ((sx % cells[0]) + cells[0]) % cells[0];
      const ix = (sx - wx) / cells[0];
      for (let dy = -reach[1]; dy <= reach[1]; dy += 1) {
        const sy = cy + dy;
        const wy = ((sy % cells[1]) + cells[1]) % cells[1];
        const iy = (sy - wy) / cells[1];
        for (let dz = -reach[2]; dz <= reach[2]; dz += 1) {
          const sz = cz + dz;
          const wz = ((sz % cells[2]) + cells[2]) % cells[2];
          const iz = (sz - wz) / cells[2];
          const bucket = buckets[(wx * cells[1] + wy) * cells[2] + wz];
          if (!bucket) continue;
          for (const row of bucket) {
            // A center is never its own neighbour in the unshifted image;
            // other images of the same atom are genuine neighbours and stay.
            if (row === centerRow && ix === 0 && iy === 0 && iz === 0) continue;
            const dxf = wrapped[row][0] + ix - fx;
            const dyf = wrapped[row][1] + iy - fy;
            const dzf = wrapped[row][2] + iz - fz;
            const vx = dxf * lattice[0][0] + dyf * lattice[1][0] + dzf * lattice[2][0];
            const vy = dxf * lattice[0][1] + dyf * lattice[1][1] + dzf * lattice[2][1];
            const vz = dxf * lattice[0][2] + dyf * lattice[1][2] + dzf * lattice[2][2];
            const distSq = vx * vx + vy * vy + vz * vz;
            // Lower bound exclusive at exactly zero even when rmin == 0: a
            // zero-length pair has no direction (mirrors the Python engine).
            if (distSq < rminSq || distSq > rmaxSq || !(distSq > 0)) continue;
            bonds.push({ row, ix, iy, iz, vx, vy, vz, length: Math.sqrt(distSq) });
          }
        }
      }
    }
    perCenter[center] = bonds;
  }
  return perCenter;
};

const toDegrees = 180 / Math.PI;

// Validate + select + bond search + angle formation: the shared core behind
// bondAngleSummary (mirrors _triplet_core + _pair_angles).
const tripletCore = (fractional, elements, latticeVectors, { triplet, bond12, bond23 = null }) => {
  if (!Array.isArray(fractional) || fractional.some((row) => !Array.isArray(row) || row.length !== 3)) {
    throw new Error('fractional must be an (N, 3) array');
  }
  if (fractional.some((row) => row.some((value) => !Number.isFinite(value)))) {
    throw new Error('fractional coordinates contain non-finite values');
  }
  const lattice = latticeVectors;
  if (!Array.isArray(lattice) || lattice.length !== 3
    || lattice.some((row) => !Array.isArray(row) || row.length !== 3 || row.some((v) => !Number.isFinite(v)))) {
    throw new Error('latticeVectors must be a finite (3, 3) matrix');
  }
  const symbols = elements.map(capitalize);
  if (symbols.length !== fractional.length) {
    throw new Error(`${symbols.length} elements for ${fractional.length} coordinates`);
  }
  if (!Array.isArray(triplet) || triplet.length !== 3) {
    throw new Error('triplet must name three atom types (A, B, C)');
  }
  const [end1, apex, end2] = triplet.map(capitalize);
  const window12 = validateWindow('bond12', bond12);
  const window23 = bond23 == null ? window12 : validateWindow('bond23', bond23);

  const selections = new Map();
  for (const symbol of new Set([end1, apex, end2])) {
    const rows = [];
    symbols.forEach((candidate, row) => {
      if (candidate === symbol) rows.push(row);
    });
    if (!rows.length) {
      const available = [...new Set(symbols)].sort().join(', ');
      throw new Error(`No '${symbol}' atoms in the configuration; available: ${available}`);
    }
    selections.set(symbol, rows);
  }

  // Fold every coordinate into [0, 1); the image bookkeeping restores the
  // true relative geometry, so pre-wrapped and drifted inputs agree.
  const wrapped = fractional.map((row) => row.map((value) => value - Math.floor(value)));

  const apexRows = selections.get(apex);
  const sharedEnds = end1 === end2
    && window12[0] === window23[0] && window12[1] === window23[1];
  const bonds12 = neighborBonds(wrapped, apexRows, selections.get(end1), lattice, ...window12);
  const bonds23 = sharedEnds
    ? bonds12
    : neighborBonds(wrapped, apexRows, selections.get(end2), lattice, ...window23);

  const angles = [];
  for (let center = 0; center < apexRows.length; center += 1) {
    const first = bonds12[center];
    const second = bonds23[center];
    for (let i = 0; i < first.length; i += 1) {
      const one = first[i];
      // Shared ends: strict upper triangle of one list, so each unordered
      // pair of distinct bonds counts once. Distinct windows: every (1→2,
      // 2→3) combination minus a bond paired with its own atom-image.
      for (let j = sharedEnds ? i + 1 : 0; j < second.length; j += 1) {
        const two = second[j];
        if (!sharedEnds
          && one.row === two.row && one.ix === two.ix && one.iy === two.iy && one.iz === two.iz) {
          continue;
        }
        const cosine = (one.vx * two.vx + one.vy * two.vy + one.vz * two.vz)
          / (one.length * two.length);
        angles.push(Math.acos(Math.min(1, Math.max(-1, cosine))) * toDegrees);
      }
    }
  }

  return {
    triplet: [end1, apex, end2],
    window12,
    window23,
    sharedEnds,
    apexCount: apexRows.length,
    bonds12,
    bonds23,
    angles
  };
};

const lengthHistogram = (perCenter, [lo, hi]) => {
  const counts = new Array(LENGTH_BINS).fill(0);
  const width = (hi - lo) / LENGTH_BINS;
  let total = 0;
  let sum = 0;
  for (const bonds of perCenter) {
    for (const bond of bonds) {
      total += 1;
      sum += bond.length;
      counts[histogramIndex(bond.length, lo, hi, LENGTH_BINS)] += 1;
    }
  }
  const binCenters = Array.from(
    { length: LENGTH_BINS },
    (_, index) => lo + (index + 0.5) * width
  );
  return {
    binCenters,
    counts,
    count: total,
    meanLength: total ? sum / total : null
  };
};

/**
 * The Bond Geometry payload: angle histogram (counts / per-degree density /
 * exact sin-corrected), bond-length histograms per window, and the
 * coordination distribution — mirrors rmc_toolkits.triplets.bond_angle_summary.
 * `collectAngles: true` adds `sortedAngles` (degrees) for tests.
 */
export const bondAngleSummary = (fractional, elements, latticeVectors, options = {}) => {
  const { binWidth = 1.0, collectAngles = false } = options;
  if (!Number.isFinite(binWidth) || binWidth <= 0 || binWidth > 180) {
    throw new Error(`binWidth must be in (0, 180], got ${binWidth}`);
  }
  const core = tripletCore(fractional, elements, latticeVectors, options);

  // Math.round is floor(x + 0.5) for positives — the same half-up rule the
  // Python engine's _bin_count uses, so both build identical bin counts.
  const nbins = Math.max(1, Math.round(180 / binWidth));
  const width = 180 / nbins;
  const counts = new Array(nbins).fill(0);
  for (const angle of core.angles) {
    counts[histogramIndex(angle, 0, 180, nbins)] += 1;
  }
  const total = core.angles.length;
  const binCenters = Array.from({ length: nbins }, (_, index) => (index + 0.5) * width);
  const density = counts.map((count) => (total ? count / (total * width) : 0));
  const sinCorrected = counts.map((count, index) => {
    if (!total) return 0;
    const lo = (index * width * Math.PI) / 180;
    const hi = ((index + 1) * width * Math.PI) / 180;
    // Exact isotropic reference fraction per bin: integral of sin(theta)/2.
    return count / total / ((Math.cos(lo) - Math.cos(hi)) / 2);
  });

  let meanAngle = null;
  let stdAngle = null;
  if (total) {
    const mean = core.angles.reduce((acc, value) => acc + value, 0) / total;
    const variance = core.angles.reduce((acc, value) => acc + (value - mean) ** 2, 0) / total;
    meanAngle = mean;
    stdAngle = Math.sqrt(variance);
  }

  // coordination[n] = how many central atoms have exactly n window-1 bonds.
  const maxBonds = core.bonds12.reduce((acc, bonds) => Math.max(acc, bonds.length), 0);
  const coordination = new Array(maxBonds + 1).fill(0);
  for (const bonds of core.bonds12) coordination[bonds.length] += 1;

  const summary = {
    triplet: core.triplet,
    bond12: core.window12,
    bond23: core.window23,
    sharedEnds: core.sharedEnds,
    binWidth: width,
    binCenters,
    counts,
    density,
    sinCorrected,
    angleCount: total,
    meanAngle,
    stdAngle,
    apexCount: core.apexCount,
    lengths12: lengthHistogram(core.bonds12, core.window12),
    lengths23: core.sharedEnds ? null : lengthHistogram(core.bonds23, core.window23),
    coordination
  };
  if (collectAngles) {
    summary.sortedAngles = [...core.angles].sort((a, b) => a - b);
  }
  return summary;
};
