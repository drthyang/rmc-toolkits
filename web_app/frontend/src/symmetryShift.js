// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Tsung-Han Yang

// web_app/frontend/src/symmetryShift.js
//
// The coherent shift of a symmetry lowering: how far each average site sits from the position
// a higher-symmetry ("parent") group puts it at. The displacement analyses cannot see it — each
// site's cloud is centred on its own mean, so a shift every copy shares moves the mean and
// leaves the cloud alone — and the tolerance ladder only says at which tolerance the parent
// stops holding. The 2026-10-14 review's L4 (approved for batch 1): the shift had to be worked
// out by hand from detect_symmetry's orbit positions.
//
//   1. The parent group at `parentTol` (analyseSymmetry: the finder's closed group, its
//      least-squares operations and orbits).
//   2. Each parent orbit averaged over its operations (averageStructure.js averageOrbit: the
//      representative projected onto its site symmetry), and every member's ideal position
//      the image of that average nearest the member.
//   3. A member's shift = its mean − its ideal position (Cartesian Å). The finder places the
//      parent's operations by least squares over every atom, so the shifts share an arbitrary
//      common part (along a polar axis the origin is free): it is removed by the geometric
//      median of the members' shifts (one atom each), which leaves the parent on the majority of
//      the atoms — the description that moves the fewest (a P4mm Ta shift reads as Ta moving,
//      not as Ta, K and O all moving by parts).
//
// Browser-only (the finder has no Python source of truth), as detect_symmetry.

import { analyseSymmetry } from './symmetryModel.js';
import { averageOrbit } from './averageStructure.js';

const ZERO3 = () => [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
const add = (u, v) => [u[0] + v[0], u[1] + v[1], u[2] + v[2]];
const sub = (u, v) => [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
const mulV = (M, v) => M.map((row) => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
const cyc = (x) => x - Math.round(x);
// A fractional row vector in Cartesian Å (cell rows A).
const toCartesian = (A, d) => [0, 1, 2].map((k) => d[0] * A[0][k] + d[1] * A[1][k] + d[2] * A[2][k]);
const norm = (v) => Math.hypot(v[0], v[1], v[2]);

/**
 * Weiszfeld's geometric median of 3-vectors (equal weights): the point with the least sum of
 * distances. Starts from the componentwise median; a point it lands on exactly is kept.
 */
export function geometricMedian(points, { iterations = 200, tolerance = 1e-12 } = {}) {
  if (!points.length) return [0, 0, 0];
  const median = (values) => {
    const sorted = [...values].sort((a, b) => a - b);
    const mid = sorted.length >> 1;
    return sorted.length % 2 ? sorted[mid] : 0.5 * (sorted[mid - 1] + sorted[mid]);
  };
  let y = [0, 1, 2].map((k) => median(points.map((p) => p[k])));
  for (let step = 0; step < iterations; step += 1) {
    let weight = 0;
    let sum = [0, 0, 0];
    let onPoint = false;
    for (const p of points) {
      const d = norm(sub(p, y));
      if (d < tolerance) { onPoint = true; continue; }
      weight += 1 / d;
      sum = add(sum, p.map((x) => x / d));
    }
    if (!(weight > 0)) break;
    const next = sum.map((x) => x / weight);
    // On a data point the update can only stay or leave towards the others;
    // Weiszfeld's rule: keep it unless the pull of the others exceeds 1.
    if (onPoint) {
      const pull = norm(sub(next, y)) * weight;
      if (pull <= 1) break;
    }
    if (norm(sub(next, y)) < tolerance) { y = next; break; }
    y = next;
  }
  return y;
}

/**
 * The shift of every basis site from the parent group's positions at `parentTol` (Å):
 *   { parent: { spaceGroup, pointGroup, tolerance }, sites: [{ basisIndex, orbit, element,
 *     referenceNumber, shift (Cartesian Å, origin removed), length, standardError }], origin }
 * or null without a basis, or the finder's 'not analysed' description (skipped: true).
 * `standardError` is the sampling error (Å) of the site's mean along its shift.
 */
export function symmetryShifts(structure, parentTol = 0.2) {
  const analysis = analyseSymmetry(structure, parentTol);
  if (!analysis || analysis.skipped) return analysis;
  const { A, sg, found } = analysis;
  const parent = { spaceGroup: sg.spaceGroup, spaceGroupNumber: sg.spaceGroupNumber ?? null, pointGroup: sg.pointGroup, tolerance: parentTol };
  if (!sg.ops?.length || !found.length) return { parent, sites: [], origin: [0, 0, 0] };
  const basis = structure.basis;
  const raw = [];
  found.forEach((orbit, orbitIndex) => {
    const members = orbit.index.map((k) => ({
      pos: basis[k].mean ?? basis[k].frac,
      cov: basis[k].covFrac ?? ZERO3(),
      weight: basis[k].count ?? 1,
    }));
    const { x } = averageOrbit(sg.ops, members, A, parentTol);
    members.forEach((member, j) => {
      let best = null;
      let bestLength = Infinity;
      for (const g of sg.ops) {
        const image = add(mulV(g.R, x), g.t);
        const d = [0, 1, 2].map((q) => cyc(member.pos[q] - image[q]));
        const length = norm(toCartesian(A, d));
        if (length < bestLength) { bestLength = length; best = d; }
      }
      raw.push({ basisIndex: orbit.index[j], orbit: orbitIndex, cartesian: toCartesian(A, best) });
    });
  });
  const origin = geometricMedian(raw.map((entry) => entry.cartesian));
  const sites = raw.map(({ basisIndex, orbit, cartesian }) => {
    const site = basis[basisIndex];
    const shift = sub(cartesian, origin);
    const length = norm(shift);
    // The mean's sampling error along the shift: Cartesian covariance AᵀCA over the copies.
    let standardError = null;
    if (site.covFrac && site.count > 0 && length > 0) {
      const u = shift.map((v) => v / length);
      // u·(Aᵀ C A)·u = (A u)ᵀ C (A u), with A u the fractional-basis projection.
      const a = mulV(A, u);
      const variance = [0, 1, 2].reduce((sum, i) => sum + [0, 1, 2].reduce((inner, j) => inner + a[i] * site.covFrac[i][j] * a[j], 0), 0);
      standardError = Math.sqrt(Math.max(variance, 0) / site.count);
    }
    return { basisIndex, orbit, element: site.el, referenceNumber: site.referenceNumber ?? basisIndex + 1, shift, length, standardError };
  });
  return { parent, sites, origin };
}
