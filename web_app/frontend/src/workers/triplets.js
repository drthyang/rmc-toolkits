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
// Angle counting mirrors the Python engine: when A and C are the same
// element, each unordered pair of distinct bond images {x, y} is one physical
// triplet, counted once when either assignment puts one bond in the A–B
// window and the other in the B–C window (with equal windows: every unordered
// pair) — continuous in the windows, so nudging a bound never doubles the
// count. Different end elements count every (A-bond, C-bond) combination.
//
// The summary payload (bondAngleSummary) is the contract shared with the
// Flask /api/triplets route — camelCase keys, plain arrays — so the page
// renders identical numbers in both runtimes. `sinCorrected` divides each
// bin's count fraction by the exact isotropic fraction (cosθ₁ − cosθ₂)/2 =
// sin θc · sin(w/2): the bin-centre 1/sin θc correction times the constant
// 1/sin(w/2), flat 1.0 for random directions (RMCProfile's TRIPLETS
// norm/sin(theta) is the same curve × sin(w/2)/w_deg ≈ π/360).
//
// Bin edges (also documented in triplets.py): an undisplaced configuration
// puts every symmetry angle (60/90/120° …) exactly on a bin edge up to float
// noise, and libm acos and Math.acos differ by 1 ulp on ~17% of inputs. So an
// angle within EDGE_SNAP_DEG of an edge is binned as exactly on it — into the
// bin that edge starts — identically in both engines: a symmetry class lands
// whole in one bin whatever the rounding. Only binning snaps; statistics and
// raw angles keep the computed values. Element symbols are ASCII, where the
// capitalize below matches str.capitalize().

export const MAX_CELLS_PER_AXIS = 64;
export const LENGTH_BINS = 40;
const REACH_HEADROOM = 1e-9;

// Work budget for one app request (mirrors APP_MAX_ANGLES in triplets.py —
// keep the two equal; the parity fixture pins it). The exact angle count is
// taken from the bond lists before any angle exists and a spec above it is
// refused, because the count grows ~rmax^6 and the 15 Å rmax cap bounds only
// the neighbour search. The engine itself is unrestricted unless the caller
// passes `maxAngles` (the worker's 'triplets' handler does).
export const APP_MAX_ANGLES = 50_000_000;

// Distances within this many Å of a window bound count as on it — inside,
// since windows are inclusive (mirrors WINDOW_TOL in triplets.py): an ideal
// shell at exactly a typed bound is kept whole, not by rounding luck.
export const WINDOW_TOL = 1e-9;

// Angles closer than this to a bin edge bin as exactly on it (mirrors
// EDGE_SNAP_DEG in triplets.py): far above float noise in an angle (~1e-13°),
// far below any bin width or real displacement.
export const EDGE_SNAP_DEG = 1e-9;

const capitalize = (symbol) => {
  const text = String(symbol).trim();
  return text.charAt(0).toUpperCase() + text.slice(1).toLowerCase();
};

// A missing request value: null/undefined, '' or whitespace only. Number()
// turns every one of them into 0 — so a cleared "rmin" box would silently
// become rmin = 0 — where Python's float() raises. Callers test this first.
export const isBlankValue = (value) =>
  value == null || (typeof value === 'string' && value.trim() === '');

const validateWindow = (name, window) => {
  if (!Array.isArray(window) || window.length !== 2) {
    throw new Error(`${name} must be a (rmin, rmax) pair`);
  }
  // Reject missing bounds explicitly so they error like Python's float()
  // instead of becoming rmin = 0.
  if (window.some(isBlankValue)) {
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

// Angle bin with the edge snap (mirrors _angle_bins): an angle within
// EDGE_SNAP_DEG of an edge goes to the bin that edge starts (180° to the
// last bin); any other angle bins exactly as numpy.histogram would.
const angleBin = (angle, nbins, width) => {
  const nearest = Math.floor(angle / width + 0.5);
  if (Math.abs(angle - nearest * width) < EDGE_SNAP_DEG) return Math.min(nearest, nbins - 1);
  return histogramIndex(angle, 0, 180, nbins);
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

// Linked-cell periodic neighbour search with explicit image shifts (mirrors
// _neighbor_bonds). Calls visit(center, row, ix, iy, iz, vx, vy, vz, distSq,
// member) for every (center, candidate image) pair whose length falls inside
// at least one of `windows`, with bit w of `member` set when window w holds
// it. Centers are visited in order, each one's pairs in stencil order.
const visitNeighbors = (wrapped, centerRows, candidateRows, lattice, windows, visit) => {
  const rmax = Math.max(...windows.map((window) => window[1]));
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

  // Inclusive bounds, widened by WINDOW_TOL against float noise.
  const loSq = windows.map(([lo]) => Math.max(lo - WINDOW_TOL, 0) ** 2);
  const hiSq = windows.map(([, hi]) => (hi + WINDOW_TOL) ** 2);
  const nWindows = windows.length;
  for (let center = 0; center < centerRows.length; center += 1) {
    const centerRow = centerRows[center];
    const [fx, fy, fz] = wrapped[centerRow];
    const cx = cellOf(fx, 0);
    const cy = cellOf(fy, 1);
    const cz = cellOf(fz, 2);
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
            // Difference first, image shift second (as in the Python engine):
            // the other end's vector of the same bond is its exact negative,
            // so a bond is in or out of a window from both ends alike.
            const dxf = (wrapped[row][0] - fx) + ix;
            const dyf = (wrapped[row][1] - fy) + iy;
            const dzf = (wrapped[row][2] - fz) + iz;
            const vx = dxf * lattice[0][0] + dyf * lattice[1][0] + dzf * lattice[2][0];
            const vy = dxf * lattice[0][1] + dyf * lattice[1][1] + dzf * lattice[2][1];
            const vz = dxf * lattice[0][2] + dyf * lattice[1][2] + dzf * lattice[2][2];
            const distSq = vx * vx + vy * vy + vz * vz;
            // Lower bound exclusive at exactly zero even when rmin == 0: a
            // zero-length pair has no direction (mirrors the Python engine).
            if (!(distSq > 0)) continue;
            let member = 0;
            for (let window = 0; window < nWindows; window += 1) {
              if (distSq >= loSq[window] && distSq <= hiSq[window]) member |= 1 << window;
            }
            if (!member) continue;
            visit(center, row, ix, iy, iz, vx, vy, vz, distSq, member);
          }
        }
      }
    }
  }
};

// All bonds grouped per center (grouping replaces the numpy sort-by-center
// step). Each bond: { row, ix, iy, iz, vx, vy, vz, length, member }.
const neighborBonds = (wrapped, centerRows, candidateRows, lattice, windows) => {
  const perCenter = Array.from({ length: centerRows.length }, () => []);
  visitNeighbors(wrapped, centerRows, candidateRows, lattice, windows,
    (center, row, ix, iy, iz, vx, vy, vz, distSq, member) => {
      perCenter[center].push({ row, ix, iy, iz, vx, vy, vz, length: Math.sqrt(distSq), member });
    });
  return perCenter;
};

// Per-center bond counts by window-membership pattern, storing no bond: the
// counting pass of a budgeted request (mirrors _neighbor_bonds(count_only)).
// patterns[center * 2^W + member] = number of that center's bonds.
const countPatterns = (wrapped, centerRows, candidateRows, lattice, windows) => {
  const nPatterns = 1 << windows.length;
  const patterns = new Float64Array(centerRows.length * nPatterns);
  visitNeighbors(wrapped, centerRows, candidateRows, lattice, windows,
    (center, row, ix, iy, iz, vx, vy, vz, distSq, member) => {
      patterns[center * nPatterns + member] += 1;
    });
  return patterns;
};

// Membership-pattern counts from stored bonds (same layout as countPatterns).
const patternsOf = (perCenter, nWindows) => {
  const nPatterns = 1 << nWindows;
  const patterns = new Float64Array(perCenter.length * nPatterns);
  perCenter.forEach((bonds, center) => {
    for (const bond of bonds) patterns[center * nPatterns + bond.member] += 1;
  });
  return patterns;
};

const pairs2 = (n) => (n * (n - 1)) / 2;

// Exact angle count summed over centers, from membership-pattern counts
// (mirrors _angles_per_center). Different ends: Σ n1·n2. Same end, one table
// over both windows (bit 0 = A–B, bit W−1 = B–C): the combined list's
// unordered pairs minus those with no bond in one window,
// C(a+b+c, 2) − C(a, 2) − C(b, 2) (a only A–B, b only B–C, c in both).
const angleTotal = (sameEnd, nWindows, nCenters, patterns1, patterns2 = null) => {
  const nPatterns = 1 << nWindows;
  let total = 0;
  for (let center = 0; center < nCenters; center += 1) {
    const base = center * nPatterns;
    if (!sameEnd) {
      total += patterns1[base + 1] * patterns2[base + 1];
      continue;
    }
    let all = 0;
    let only12 = 0;
    let only23 = 0;
    for (let pattern = 1; pattern < nPatterns; pattern += 1) {
      const count = patterns1[base + pattern];
      const in12 = (pattern & 1) !== 0;
      const in23 = ((pattern >> (nWindows - 1)) & 1) !== 0;
      all += count;
      if (in12 && !in23) only12 += count;
      else if (in23 && !in12) only23 += count;
    }
    total += pairs2(all) - pairs2(only12) - pairs2(only23);
  }
  return total;
};

const toDegrees = 180 / Math.PI;

// Validate + select + bond search + exact angle count: the shared core behind
// bondAngleSummary (mirrors _triplet_core + _Pairing.angle_counts). No angle
// is formed here, so an over-budget spec is refused before any pairing work.
const tripletCore = (fractional, elements, latticeVectors,
  { triplet, bond12, bond23 = null, maxAngles = null }) => {
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
  const sameEnd = end1 === end2;
  const sharedEnds = sameEnd
    && window12[0] === window23[0] && window12[1] === window23[1];

  const windowsSameEnd = sharedEnds ? [window12] : [window12, window23];
  const refuse = (count) => {
    throw new Error(
      `${end1}-${apex}-${end2} with bond windows ${window12[0]}-${window12[1]} / `
      + `${window23[0]}-${window23[1]} A would form ${count.toLocaleString('en-US')} angles, `
      + `over the limit of ${Number(maxAngles).toLocaleString('en-US')} for one request; `
      + 'narrow the bond windows (the angle count grows roughly as rmax^6)'
    );
  };
  if (maxAngles != null) {
    // Budgeted (app) request: an exact count from a search that stores
    // nothing, so an oversized spec is refused before either the bond lists
    // or a single angle take up memory.
    const counted = sameEnd
      ? angleTotal(true, windowsSameEnd.length, apexRows.length,
        countPatterns(wrapped, apexRows, selections.get(end1), lattice, windowsSameEnd))
      : angleTotal(false, 1, apexRows.length,
        countPatterns(wrapped, apexRows, selections.get(end1), lattice, [window12]),
        countPatterns(wrapped, apexRows, selections.get(end2), lattice, [window23]));
    if (counted > maxAngles) refuse(counted);
  }

  // Same end element: one search over both windows, each bond image once,
  // tagged with its window(s) (bit 1 = A–B, bit 2 = B–C; one window = bit 1).
  // Different end elements: one search per end.
  const both = sameEnd
    ? neighborBonds(wrapped, apexRows, selections.get(end1), lattice, windowsSameEnd)
    : null;
  const bit23 = sharedEnds ? 1 : 2;
  const bonds12 = sameEnd
    ? both.map((bonds) => bonds.filter((bond) => bond.member & 1))
    : neighborBonds(wrapped, apexRows, selections.get(end1), lattice, [window12]);
  const bonds23 = sameEnd
    ? (sharedEnds ? bonds12 : both.map((bonds) => bonds.filter((bond) => bond.member & bit23)))
    : neighborBonds(wrapped, apexRows, selections.get(end2), lattice, [window23]);

  // Exact angle count (equal to the counting pass's when one ran).
  const angleCount = sameEnd
    ? angleTotal(true, windowsSameEnd.length, apexRows.length, patternsOf(both, windowsSameEnd.length))
    : angleTotal(false, 1, apexRows.length, patternsOf(bonds12, 1), patternsOf(bonds23, 1));

  return {
    triplet: [end1, apex, end2],
    window12,
    window23,
    sharedEnds,
    sameEnd,
    bit23,
    apexCount: apexRows.length,
    both,
    bonds12,
    bonds23,
    angleCount
  };
};

// Visit every angle (degrees) of a core in the engine's deterministic order,
// without storing any: the histogram and moments accumulate as it streams.
const forEachAngle = (core, visit) => {
  const { sameEnd, bit23, both, bonds12, bonds23 } = core;
  for (let center = 0; center < core.apexCount; center += 1) {
    const first = sameEnd ? both[center] : bonds12[center];
    const second = sameEnd ? both[center] : bonds23[center];
    for (let i = 0; i < first.length; i += 1) {
      const one = first[i];
      // Same end element: strict upper triangle of the one list, so each
      // unordered pair of distinct bond images counts once, kept when either
      // assignment puts one bond in each window. Otherwise every (A-bond,
      // C-bond) combination.
      for (let j = sameEnd ? i + 1 : 0; j < second.length; j += 1) {
        const two = second[j];
        if (sameEnd
          && !(((one.member & 1) && (two.member & bit23)) || ((one.member & bit23) && (two.member & 1)))) {
          continue;
        }
        const cosine = (one.vx * two.vx + one.vy * two.vy + one.vz * two.vz)
          / (one.length * two.length);
        visit(Math.acos(Math.min(1, Math.max(-1, cosine))) * toDegrees);
      }
    }
  }
};

// Bond-length histogram of one window. `count` is the number of B-centred
// bond vectors (the histogram total); `uniqueBonds` the physical bonds, each
// once — half of `count` when the end element is the central one, since each
// such bond is found from both of its ends (mirrors _unique_bonds).
const lengthHistogram = (perCenter, [lo, hi], homonuclear) => {
  const counts = new Array(LENGTH_BINS).fill(0);
  const width = (hi - lo) / LENGTH_BINS;
  let total = 0;
  let sum = 0;
  for (const bonds of perCenter) {
    for (const bond of bonds) {
      total += 1;
      sum += bond.length;
      // Clipped into the window (a bond admitted by WINDOW_TOL just outside
      // a bound goes to the edge bin), exactly as the Python engine.
      counts[histogramIndex(Math.min(Math.max(bond.length, lo), hi), lo, hi, LENGTH_BINS)] += 1;
    }
  }
  const binCenters = Array.from(
    { length: LENGTH_BINS },
    (_, index) => lo + (index + 0.5) * width
  );
  if (homonuclear && total % 2) {
    throw new Error(`internal error: ${total} directed homonuclear bonds do not pair up`);
  }
  return {
    binCenters,
    counts,
    count: total,
    uniqueBonds: homonuclear ? total / 2 : total,
    meanLength: total ? sum / total : null
  };
};

/**
 * The Bond Geometry page's input boxes → the flat `triplets` request (the
 * shape both the worker and /api/triplets take). Every box is a string; a
 * cleared or non-numeric one throws an Error naming the field — it is never
 * sent as 0. The B–C window is included only when `split23` is on.
 */
export const tripletRequestFromInputs = ({
  end1, apex, end2, r12Min, r12Max, split23, r23Min, r23Max, binWidth
}) => {
  const number = (value, label) => {
    if (isBlankValue(value)) throw new Error(`${label} is empty — enter a number.`);
    const parsed = Number(value);
    if (!Number.isFinite(parsed)) throw new Error(`${label} is not a number: "${value}".`);
    return parsed;
  };
  const params = {
    end1,
    apex,
    end2,
    r12Min: number(r12Min, 'A–B window minimum'),
    r12Max: number(r12Max, 'A–B window maximum'),
    binWidth: number(binWidth, 'Bin width')
  };
  if (split23) {
    params.r23Min = number(r23Min, 'B–C window minimum');
    params.r23Max = number(r23Max, 'B–C window maximum');
  }
  return params;
};

/**
 * The Bond Geometry payload: angle histogram (counts / per-degree density /
 * exact sin-corrected), bond-length histograms per window, and the
 * coordination distribution — mirrors rmc_toolkits.triplets.bond_angle_summary.
 * Angles stream straight into the histogram and running moments (Welford),
 * so memory does not grow with the angle count; `collectAngles: true` also
 * keeps them and adds `sortedAngles` (degrees) for tests. `maxAngles` refuses
 * a spec whose exact angle count exceeds it, before any angle is formed.
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
  const collected = collectAngles ? [] : null;
  let total = 0;
  let mean = 0;
  let m2 = 0;
  forEachAngle(core, (angle) => {
    counts[angleBin(angle, nbins, width)] += 1;
    total += 1;
    const delta = angle - mean;
    mean += delta / total;
    m2 += delta * (angle - mean);
    if (collected) collected.push(angle);
  });
  const binCenters = Array.from({ length: nbins }, (_, index) => (index + 0.5) * width);
  const density = counts.map((count) => (total ? count / (total * width) : 0));
  const sinCorrected = counts.map((count, index) => {
    if (!total) return 0;
    const lo = (index * width * Math.PI) / 180;
    const hi = ((index + 1) * width * Math.PI) / 180;
    // Exact isotropic reference fraction per bin: integral of sin(theta)/2.
    return count / total / ((Math.cos(lo) - Math.cos(hi)) / 2);
  });

  const meanAngle = total ? mean : null;
  const stdAngle = total ? Math.sqrt(m2 / total) : null;

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
    lengths12: lengthHistogram(core.bonds12, core.window12, core.triplet[0] === core.triplet[1]),
    lengths23: core.sharedEnds
      ? null
      : lengthHistogram(core.bonds23, core.window23, core.triplet[2] === core.triplet[1]),
    coordination
  };
  if (collected) {
    summary.sortedAngles = collected.sort((a, b) => a - b);
  }
  return summary;
};
