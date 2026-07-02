// web_app/frontend/src/symmetry.js
//
// Pure, fully-offline crystal symmetry finder (a lightweight FINDSYM-like tool),
// ported from the RMC-phonon-dynamics app (web/src/math/symmetry.js). Given the
// conventional cell A_conv (rows = lattice vectors, Å) + a basis and a tolerance,
// it returns the space-group operations {R|t} that map the basis onto itself,
// plus the residual of that fit, so the symmetry can be traced as a function of
// tolerance on a disordered RMC average structure — no external dependency, no
// WASM, runs client-side in the static dashboard.
//
// Method (spglib-lite, but table-free and bounded):
//   1. Point operations = integer matrices R (entries in {-1,0,1}) that preserve
//      the lattice metric G = A·Aᵀ, i.e. RᵀGR = G. These are the lattice
//      automorphisms in the CONVENTIONAL direct basis; |det R| = 1.
//   2. For each R, the space-group translations t are found from atom images:
//      t = frac(b) − R·frac(a0) for candidate partners b (same element), kept when
//      {R|t} maps every atom onto a same-element atom within `tol` (cartesian Å).
//
// Fractional coords are COLUMN vectors here: x' = R·x + t. Lattice rows: A = [a1,a2,a3].
// NOTE: entries in {-1,0,1} cover the standard conventional settings (cubic,
// tetragonal, orthorhombic, hexagonal, rhombohedral-in-hex, monoclinic, triclinic).

/** Determinant of a 3×3 matrix (rows). */
export function det3(m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

const wrap01 = (x) => x - Math.floor(x);
const nearestInt = (x) => x - Math.round(x);
const meanEdge = (A) => (Math.hypot(A[0][0], A[0][1], A[0][2]) + Math.hypot(A[1][0], A[1][1], A[1][2]) + Math.hypot(A[2][0], A[2][1], A[2][2])) / 3;

/** Metric tensor G = A·Aᵀ (A rows = lattice vectors, cartesian). */
export function metricTensor(A) {
  const G = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++)
    G[i][j] = A[i][0] * A[j][0] + A[i][1] * A[j][1] + A[i][2] * A[j][2];
  return G;
}

// Rᵀ · G · R for a 3×3 integer R and symmetric G.
function conjugate(R, G) {
  // M = Rᵀ G R.  (Rᵀ)_{ij} = R_{ji}
  const GR = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++)
    GR[i][j] = G[i][0] * R[0][j] + G[i][1] * R[1][j] + G[i][2] * R[2][j];
  const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++)
    M[i][j] = R[0][i] * GR[0][j] + R[1][i] * GR[1][j] + R[2][i] * GR[2][j];
  return M;
}

/**
 * Lattice point operations: integer R (entries in {-1,0,1}, |det|=1) with RᵀGR=G.
 * `tol` is a RELATIVE metric tolerance (fraction of the mean squared edge).
 */
export function latticePointOps(A, tol = 1e-3) {
  const G = metricTensor(A);
  const scale = (G[0][0] + G[1][1] + G[2][2]) / 3 || 1;
  const eps = tol * scale;
  const ops = [];
  const v = [-1, 0, 1];
  const R = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  // Iterate all 3^9 sign/zero patterns; keep metric-preserving unimodular ones.
  for (let code = 0; code < 19683; code++) {
    let c = code;
    for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) { R[a][b] = v[c % 3]; c = (c / 3) | 0; }
    const d = det3(R);
    if (d !== 1 && d !== -1) continue;
    const M = conjugate(R, G);
    let ok = true;
    for (let i = 0; i < 3 && ok; i++) for (let j = 0; j < 3 && ok; j++)
      if (Math.abs(M[i][j] - G[i][j]) > eps) ok = false;
    if (ok) ops.push(R.map(r => r.slice()));
  }
  return ops;
}

// Apply R (integer) to a fractional column vector.
function applyR(R, x) {
  return [
    R[0][0] * x[0] + R[0][1] * x[1] + R[0][2] * x[2],
    R[1][0] * x[0] + R[1][1] * x[1] + R[1][2] * x[2],
    R[2][0] * x[0] + R[2][1] * x[1] + R[2][2] * x[2],
  ];
}

// Cartesian distance between two fractional points (minimal image), rows A.
function cartDist(fa, fb, A) {
  const d = [nearestInt(fa[0] - fb[0]), nearestInt(fa[1] - fb[1]), nearestInt(fa[2] - fb[2])];
  const c = [
    d[0] * A[0][0] + d[1] * A[1][0] + d[2] * A[2][0],
    d[0] * A[0][1] + d[1] * A[1][1] + d[2] * A[2][1],
    d[0] * A[0][2] + d[1] * A[1][2] + d[2] * A[2][2],
  ];
  return Math.hypot(c[0], c[1], c[2]);
}

// For op {R|t}, the worst cartesian mapping error over all atoms (∞ if some atom
// has no same-element image within `tol`).
function mappingResidual(R, t, basis, byEl, A, tol) {
  let worst = 0;
  for (const s of basis) {
    const img = applyR(R, s.frac);
    img[0] = wrap01(img[0] + t[0]); img[1] = wrap01(img[1] + t[1]); img[2] = wrap01(img[2] + t[2]);
    let best = Infinity;
    for (const o of byEl.get(s.el)) { const dd = cartDist(img, o.frac, A); if (dd < best) best = dd; }
    if (best > tol) return Infinity;
    if (best > worst) worst = best;
  }
  return worst;
}

/**
 * Space-group operations of (A, basis) within a cartesian tolerance `tol` (Å).
 * basis: [{ el, frac:[x,y,z] }]. Returns operations + summary counts + residual.
 *
 * @returns {{ ops:{R,t,residual}[], nSpace, nPoint, order, maxResidual }}
 *   nPoint : distinct rotation parts present in the space group (its point group).
 *   nSpace : total {R|t} (= point-group order × #centering-type cosets for the cell).
 */
export function findSpaceGroupOps(A, basis, tol = 0.1, metricTol = 1e-2) {
  if (!basis || basis.length === 0) return { ops: [], nSpace: 0, nPoint: 0, order: 0, maxResidual: 0, centering: 'P', pointGroup: '1', spaceGroup: 'P1', spaceGroupNumber: 1 };
  const tolFrac = tol / meanEdge(A); // ~ tol in cell fractions (mean edge, not just |a1|)
  const pointOps = latticePointOps(A, metricTol);
  const byEl = new Map();
  for (const s of basis) { if (!byEl.has(s.el)) byEl.set(s.el, []); byEl.get(s.el).push(s); }

  const ops = [];
  const rotSeen = new Set();
  let maxResidual = 0;
  // Use the rarest element for candidate translations (fewest partners → fastest).
  let refEl = basis[0].el;
  for (const [el, arr] of byEl) if (arr.length < byEl.get(refEl).length) refEl = el;
  const refAtom = byEl.get(refEl)[0];

  for (const R of pointOps) {
    const Ra0 = applyR(R, refAtom.frac);
    const tSeen = [];
    for (const cand of byEl.get(refEl)) {
      const t = [wrap01(cand.frac[0] - Ra0[0]), wrap01(cand.frac[1] - Ra0[1]), wrap01(cand.frac[2] - Ra0[2])];
      if (tSeen.some(u => cartDist(u, t, A) < tol)) continue;   // dedupe translations mod lattice
      const res = mappingResidual(R, t, basis, byEl, A, tol);
      if (res === Infinity) continue;
      tSeen.push(t);
      ops.push({ R, t, residual: res });
      if (res > maxResidual) maxResidual = res;
    }
    if (tSeen.length) rotSeen.add(R.flat().join(','));
  }
  return { ops, order: ops.length, maxResidual, ...classifyOperations(ops, tolFrac) };
}

const isIdentityR = (R) => R[0][0] === 1 && R[1][1] === 1 && R[2][2] === 1
  && R[0][1] === 0 && R[0][2] === 0 && R[1][0] === 0 && R[1][2] === 0 && R[2][0] === 0 && R[2][1] === 0;

/** Classify a set of operations {R,t} into point group + centering → H–M symbol. */
export function classifyOperations(ops, tolFrac = 0.02) {
  const rotMap = new Map();
  const centerings = [];
  const transSeen = new Set();          // distinct pure (identity-rotation) translations
  for (const { R, t } of ops) {
    const key = R.flat().join(',');
    if (!rotMap.has(key)) rotMap.set(key, R);
    if (isIdentityR(R)) {
      transSeen.add(t.map(x => Math.round((((x % 1) + 1) % 1) * 1000)).join(','));
      if (t[0] > tolFrac || t[1] > tolFrac || t[2] > tolFrac) centerings.push(t);
    }
  }
  const centering = matchCentering(centerings);
  const pointGroup = pointGroupOf([...rotMap.values()]);
  const sg = spaceGroupHM(centering, pointGroup);
  return { centering, pointGroup, spaceGroup: sg.symbol, spaceGroupNumber: sg.number, nSpace: ops.length, nPoint: rotMap.size, nTrans: transSeen.size };
}

const POINT_GROUP_ORDER = {
  '1': 1, '-1': 2, '2': 2, 'm': 2, '2/m': 4, '222': 4, 'mm2': 4, 'mmm': 8,
  '4': 4, '-4': 4, '4/m': 8, '422': 8, '4mm': 8, '-42m': 8, '4/mmm': 16,
  '3': 3, '-3': 6, '32': 6, '3m': 6, '-3m': 12, '6': 6, '-6': 6, '6/m': 12,
  '622': 12, '6mm': 12, '-6m2': 12, '6/mmm': 24,
  '23': 12, 'm-3': 24, '432': 24, '-43m': 24, 'm-3m': 48,
};
// point group → crystal system → centerings that system allows.
const PG_SYSTEM = {
  '1': 'tri', '-1': 'tri', '2': 'mono', 'm': 'mono', '2/m': 'mono',
  '222': 'orth', 'mm2': 'orth', 'mmm': 'orth',
  '4': 'tet', '-4': 'tet', '4/m': 'tet', '422': 'tet', '4mm': 'tet', '-42m': 'tet', '4/mmm': 'tet',
  '3': 'trig', '-3': 'trig', '32': 'trig', '3m': 'trig', '-3m': 'trig',
  '6': 'hex', '-6': 'hex', '6/m': 'hex', '622': 'hex', '6mm': 'hex', '-6m2': 'hex', '6/mmm': 'hex',
  '23': 'cub', 'm-3': 'cub', '432': 'cub', '-43m': 'cub', 'm-3m': 'cub',
};
const ALLOWED_CENTERING = { tri: 'P', mono: 'PC', orth: 'PCIFAB', tet: 'PI', trig: 'PR', hex: 'P', cub: 'PFI' };
// A classified op set is a real space group only if (a) it closes: the op count
// equals point-group order × the actual number of pure translations (partial
// mid-transition sets don't — and a tiled supercell has extra translations, so we
// count them rather than assuming a standard centering multiplicity), and (b) the
// centering is compatible with the point group's crystal system.
function isValidGroup(cls) {
  const pg = POINT_GROUP_ORDER[cls.pointGroup] || 0;
  if (pg === 0 || cls.nSpace !== pg * (cls.nTrans || 1)) return false;
  return (ALLOWED_CENTERING[PG_SYSTEM[cls.pointGroup]] || 'P').includes(cls.centering);
}

/**
 * Symmetry-vs-tolerance ladder in ONE detection pass. Detect at the loosest
 * tolerance (all candidate ops with their residuals), then threshold: an op holds
 * at tolerance t iff its residual ≤ t. Distinct thresholds where the qualifying set
 * changes → the rungs, each a space group over a tolerance range (merged when the
 * group repeats). Monotonic: looser tol ⇒ higher symmetry (P1 → … → full group).
 *
 * @returns {{from:number, to:number, spaceGroup:string, spaceGroupNumber:number|null,
 *            pointGroup:string, nSpace:number}[]} bricks, tight→loose.
 */
export function symmetryLadder(A, basis, tolMax = 1.5, metricTol = 1e-2) {
  const full = findSpaceGroupOps(A, basis, tolMax, metricTol);
  if (!full.ops.length) return [];
  const tolFrac = tolMax / meanEdge(A);
  const thresholds = [...new Set(full.ops.map(o => o.residual))].sort((a, b) => a - b);
  const bricks = [];
  for (let i = 0; i < thresholds.length; i++) {
    const r = thresholds[i];
    const to = i + 1 < thresholds.length ? thresholds[i + 1] : tolMax;  // full group holds for all looser tol
    const cls = classifyOperations(full.ops.filter(o => o.residual <= r + 1e-9), tolFrac);
    // Skip partial op sets that aren't a closed group; extend the current rung.
    if (!isValidGroup(cls)) { if (bricks.length) bricks[bricks.length - 1].to = to; continue; }
    const last = bricks[bricks.length - 1];
    if (last && last.spaceGroup === cls.spaceGroup) { last.to = to; last.nSpace = cls.nSpace; }
    else bricks.push({ from: r, to, spaceGroup: cls.spaceGroup, spaceGroupNumber: cls.spaceGroupNumber, pointGroup: cls.pointGroup, nSpace: cls.nSpace });
  }
  if (bricks.length) bricks[0].from = 0;
  return bricks;
}

/**
 * The VALID space group holding at cartesian tolerance `tol` (Å), with its
 * operations. Same rule as symmetryLadder: take every op whose residual ≤ tol,
 * but if that raw set is not a closed group (a partial mid-transition set can
 * classify to a spurious symbol like "B-43m"), fall back to the largest valid
 * closed subgroup at a tighter threshold. Guarantees a real space group.
 *
 * @returns {{ centering, pointGroup, spaceGroup, spaceGroupNumber, nSpace, nPoint,
 *             maxResidual, ops:{R,t,residual}[] }}
 */
export function spaceGroupAtTolerance(A, basis, tol = 0.2, metricTol = 1e-2) {
  const full = findSpaceGroupOps(A, basis, Math.max(tol, 1e-3), metricTol);
  const empty = { centering: 'P', pointGroup: '1', spaceGroup: 'P1', spaceGroupNumber: 1, nSpace: 0, nPoint: 0, maxResidual: 0, ops: [] };
  if (!full.ops.length) return empty;
  const tolFrac = Math.max(tol, 1e-6) / meanEdge(A);
  // Distinct residual thresholds ≤ tol, largest first: the valid group at `tol`
  // is the classification at the loosest threshold ≤ tol that closes.
  const thresholds = [...new Set(full.ops.map(o => o.residual))].filter(r => r <= tol + 1e-9).sort((a, b) => b - a);
  for (const r of thresholds) {
    const ops = full.ops.filter(o => o.residual <= r + 1e-9);
    const cls = classifyOperations(ops, tolFrac);
    if (isValidGroup(cls)) return { ...cls, maxResidual: r, ops };
  }
  // Nothing closes at/below `tol` → the trivial group (identity + pure lattice).
  const idOps = full.ops.filter(o => isIdentityR(o.R));
  return { ...classifyOperations(idOps, tolFrac), maxResidual: 0, ops: idOps };
}

// Match a set of fractional centering translations against the Bravais centerings.
const CENTERING_SETS = {
  F: [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
  I: [[0.5, 0.5, 0.5]],
  A: [[0, 0.5, 0.5]], B: [[0.5, 0, 0.5]], C: [[0.5, 0.5, 0]],
};
function matchCentering(translations, tol = 0.1) {
  const has = (v) => translations.some(t => Math.abs(((t[0] - v[0] + 0.5) % 1) - 0.5) < tol
    && Math.abs(((t[1] - v[1] + 0.5) % 1) - 0.5) < tol && Math.abs(((t[2] - v[2] + 0.5) % 1) - 0.5) < tol);
  for (const [letter, vecs] of Object.entries(CENTERING_SETS)) if (vecs.every(has)) return letter;
  return 'P';
}

/* ── Space-group identification ──────────────────────────────────────────────
 * Classify the detected operations into a point group + centering → Hermann–
 * Mauguin symbol. Rotation TYPE is basis-independent (det & trace are similarity
 * invariants): proper (det +1) trace 3,2,1,0,−1 → 1,6,4,3,2-fold; improper
 * (det −1) trace −3,1,−1,0,−2 → inversion, mirror m, −4, −3, −6. The point group
 * is then fixed by (crystal system, order, has-inversion, which proper folds).
 * Symmorphic H–M = centering letter + point-group symbol (correct for GaTa4Se8
 * F-43m and most cases; non-symmorphic screw/glide refinement is a follow-up). */

export function classifyRotation(R) {
  const d = det3(R), t = R[0][0] + R[1][1] + R[2][2];
  if (d === 1) return t === 3 ? '1' : t === 2 ? '6' : t === 1 ? '4' : t === 0 ? '3' : '2';
  return t === -3 ? '-1' : t === 1 ? 'm' : t === -1 ? '-4' : t === 0 ? '-3' : '-6';
}

// The distinct rotation parts of a space group → its point-group H–M symbol.
// The crystal CLASS is derived from the rotation content itself (not the lattice
// metric), so a structure whose symmetry is a proper subgroup of its lattice — the
// generic case along the tolerance ladder — is classified correctly.
export function pointGroupOf(rotations) {
  const h = { '1': 0, '2': 0, '3': 0, '4': 0, '6': 0, '-1': 0, 'm': 0, '-3': 0, '-4': 0, '-6': 0 };
  for (const R of rotations) h[classifyRotation(R)]++;
  const order = rotations.length, inv = h['-1'] > 0, nm = h.m;
  if (h['3'] >= 8) {                                     // cubic: 4 three-fold axes
    if (order === 12) return '23';
    if (order === 48) return 'm-3m';
    return inv ? 'm-3' : (h['4'] > 0 ? '432' : '-43m');
  }
  if (h['6'] > 0 || h['-6'] > 0) {                       // hexagonal
    if (order === 24) return '6/mmm';
    if (order === 12) return inv ? '6/m' : (h['6'] > 0 ? (nm > 0 ? '6mm' : '622') : '-6m2');
    return h['6'] > 0 ? '6' : '-6';
  }
  if (h['4'] > 0 || h['-4'] > 0) {                       // tetragonal
    if (order === 16) return '4/mmm';
    if (order === 8) return inv ? '4/m' : (h['4'] > 0 ? (nm > 0 ? '4mm' : '422') : '-42m');
    return h['4'] > 0 ? '4' : '-4';
  }
  if (h['3'] > 0 || h['-3'] > 0) {                       // trigonal
    if (order === 12) return '-3m';
    if (order === 6) return inv ? '-3' : (nm > 0 ? '3m' : '32');
    return '3';
  }
  if (h['2'] > 0 || nm > 0) {                            // ortho / mono
    if (order === 8) return 'mmm';
    if (order === 4) return inv ? '2/m' : (h['2'] >= 3 ? '222' : 'mm2');
    return h['2'] > 0 ? '2' : 'm';
  }
  return inv ? '-1' : '1';                               // triclinic
}

// Symmorphic space-group number for centering + point group (common groups).
const SG_NUMBER = {
  'P1': 1, 'P-1': 2, 'P2': 3, 'C2': 5, 'Pm': 6, 'Cm': 8, 'P2/m': 10, 'C2/m': 12,
  'P222': 16, 'C222': 21, 'A222': 21, 'B222': 21, 'F222': 22, 'I222': 23,
  'Pmm2': 25, 'Fmm2': 42, 'Imm2': 44,
  'Pmmm': 47, 'Cmmm': 65, 'Ammm': 65, 'Bmmm': 65, 'Fmmm': 69, 'Immm': 71,
  'P4': 75, 'P-4': 81, 'P4/m': 83, 'P422': 89, 'P4mm': 99, 'P-42m': 111, 'P4/mmm': 123,
  'I4': 79, 'I-4': 82, 'I4/m': 87, 'I422': 97, 'I4mm': 107, 'I-42m': 121, 'I4/mmm': 139,
  'Cmm2': 35, 'Amm2': 38, 'Bmm2': 38,
  'P3': 143, 'P-3': 147, 'P32': 149, 'P3m': 156, 'P-3m': 162, 'R3': 146, 'R-3m': 166,
  'P6': 168, 'P-6': 174, 'P6/m': 175, 'P622': 177, 'P6mm': 183, 'P-6m2': 187, 'P6/mmm': 191,
  'P23': 195, 'F23': 196, 'I23': 197, 'Pm-3': 200, 'Fm-3': 202, 'Im-3': 204,
  'P432': 207, 'F432': 209, 'I432': 211, 'P-43m': 215, 'F-43m': 216, 'I-43m': 217,
  'Pm-3m': 221, 'Fm-3m': 225, 'Im-3m': 229,
};

// Combine centering letter + point group into the (symmorphic) H–M symbol + number.
export function spaceGroupHM(centering, pointGroup) {
  const symbol = (centering || 'P') + pointGroup;
  return { symbol, number: SG_NUMBER[symbol] || null };
}

/**
 * Partition the basis into symmetry orbits under the operations `ops`: two sites
 * are in the same orbit if some {R|t} maps one onto the other (same element,
 * within `tol` Å). Orbits = the Wyckoff structure of the (average) crystal — the
 * sets of sites the detected symmetry says are equivalent.
 *
 * Each orbit also carries its multiplicity (size), site-symmetry point group (the
 * stabilizer of a representative — robust, no tables), and a representative frac.
 *
 * @returns {{ index:number[], element:string, size:number, site:string, rep:number[] }[]}
 *   orbits, largest first.
 */
export function siteOrbits(A, basis, ops, tol = 0.1) {
  const n = basis.length;
  const parent = Array.from({ length: n }, (_, i) => i);
  const find = (x) => { while (parent[x] !== x) { parent[x] = parent[parent[x]]; x = parent[x]; } return x; };
  const union = (a, b) => { const ra = find(a), rb = find(b); if (ra !== rb) parent[ra] = rb; };

  const byEl = new Map();
  basis.forEach((s, i) => { if (!byEl.has(s.el)) byEl.set(s.el, []); byEl.get(s.el).push(i); });

  for (const { R, t } of ops) {
    for (let i = 0; i < n; i++) {
      const img = applyR(R, basis[i].frac);
      img[0] = wrap01(img[0] + t[0]); img[1] = wrap01(img[1] + t[1]); img[2] = wrap01(img[2] + t[2]);
      let best = -1, bestD = tol;
      for (const j of byEl.get(basis[i].el)) { const d = cartDist(img, basis[j].frac, A); if (d < bestD) { bestD = d; best = j; } }
      if (best >= 0) union(i, best);
    }
  }
  const groups = new Map();
  for (let i = 0; i < n; i++) { const r = find(i); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(i); }
  return [...groups.values()]
    .map(index => {
      const rep = basis[index[0]];
      // Site symmetry = rotation parts of the ops that FIX the representative point.
      const stab = []; const seen = new Set();
      for (const { R, t } of ops) {
        const img = applyR(R, rep.frac);
        img[0] = wrap01(img[0] + t[0]); img[1] = wrap01(img[1] + t[1]); img[2] = wrap01(img[2] + t[2]);
        if (cartDist(img, rep.frac, A) < tol) { const k = R.flat().join(','); if (!seen.has(k)) { seen.add(k); stab.push(R); } }
      }
      return { index, element: rep.el, size: index.length, site: pointGroupOf(stab), rep: rep.frac };
    })
    .sort((a, b) => b.size - a.size);
}

// ── Wyckoff letters (common groups; letter omitted when not confident) ───────
// [letter, multiplicity, site symmetry, fixed-coord or null (free position)].
const WYCKOFF = {
  216: [ // F-43m
    ['a', 4, '-43m', [0, 0, 0]], ['b', 4, '-43m', [0.5, 0.5, 0.5]],
    ['c', 4, '-43m', [0.25, 0.25, 0.25]], ['d', 4, '-43m', [0.75, 0.75, 0.75]],
    ['e', 16, '3m', null], ['i', 96, '1', null],
  ],
  225: [ // Fm-3m
    ['a', 4, 'm-3m', [0, 0, 0]], ['b', 4, 'm-3m', [0.5, 0.5, 0.5]],
    ['c', 8, '-43m', [0.25, 0.25, 0.25]], ['f', 32, '3m', null], ['l', 192, '1', null],
  ],
  229: [ // Im-3m
    ['a', 2, 'm-3m', [0, 0, 0]], ['b', 6, '4/mmm', [0, 0.5, 0.5]],
    ['c', 8, '-3m', [0.25, 0.25, 0.25]], ['f', 16, '3m', null], ['l', 96, '1', null],
  ],
  221: [ // Pm-3m
    ['a', 1, 'm-3m', [0, 0, 0]], ['b', 1, 'm-3m', [0.5, 0.5, 0.5]],
    ['c', 3, '4/mmm', [0, 0.5, 0.5]], ['d', 3, '4/mmm', [0.5, 0, 0]],
    ['g', 8, '3m', null], ['n', 48, '1', null],
  ],
};
const CEN_VECS = {
  P: [[0, 0, 0]], I: [[0, 0, 0], [0.5, 0.5, 0.5]],
  F: [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]],
  C: [[0, 0, 0], [0.5, 0.5, 0]], A: [[0, 0, 0], [0, 0.5, 0.5]], B: [[0, 0, 0], [0.5, 0, 0.5]],
};

/**
 * Wyckoff letter for an orbit (multiplicity + site symmetry + representative),
 * or null when not confidently determined (unlisted group, non-standard setting,
 * or ambiguous). Special positions are matched modulo the centering.
 */
export function wyckoffLetter(sgNumber, centering, mult, site, rep) {
  const table = WYCKOFF[sgNumber];
  if (!table) return null;
  const cands = table.filter(p => p[1] === mult && p[2] === site);
  if (cands.length === 1 && !cands[0][3]) return cands[0][0];   // unique free position
  const cens = CEN_VECS[centering] || CEN_VECS.P;
  const hit = cands.filter(p => p[3] && cens.some(c => {
    const near = (a, b) => { const d = ((a - b) % 1 + 1.5) % 1 - 0.5; return Math.abs(d) < 0.15; };
    return near(rep[0], p[3][0] + c[0]) && near(rep[1], p[3][1] + c[1]) && near(rep[2], p[3][2] + c[2]);
  }));
  if (hit.length === 1) return hit[0][0];
  if (cands.length === 1) return cands[0][0];
  return null;
}
