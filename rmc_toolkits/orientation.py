# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Orientation distribution of RMC displacement vectors, hex-binned on a sphere.

The thermal-ellipsoid view in :mod:`rmc_toolkits.pca_kde` describes a site's
displacement cloud by its second moment: the covariance of ``dr = r_i - r_avg``,
which is an ellipsoid and therefore always centrosymmetric and always convex.
That is the harmonic picture. This module asks a different, complementary
question of the *same* cloud -- ignore how far each atom moved and keep only
*which way* it moved:

    u_i = dr_i / |dr_i|,   a point on the unit sphere,

then estimate the probability density of ``u`` over solid angle. A harmonic
isotropic site gives a uniform sphere; a harmonic anisotropic site gives a
smooth two-lobed pattern along the long ellipsoid axis; a site hopping between
discrete off-centre positions gives *discrete spots* that no ellipsoid can
represent, and a site with odd-order anharmonicity gives a sphere that is not
antipodally symmetric. None of those three signatures is visible in U.

Reference point: ``dr`` is measured from the configuration's *own* site mean
(``pca_kde.load_site_displacements`` subtracts it), so ``sum_i dr_i = 0`` for
every site. A coherent (ordered) off-centring -- every copy displaced the same
way from the ideal position -- therefore lives entirely in the site mean and
leaves this map, ``antipodalAsymmetry`` and the orientation tensor exactly as
they would be without it; read it from the site mean against the ideal
(Wyckoff / average-structure) position instead. What the map does show is the
*shape* of the cloud about its mean: a +u/-u imbalance is skewness (odd-order
anharmonicity, or unequal occupation of opposite off-centre wells), and with
partial off-centring the minority of copies left near the ideal site is what
appears as a lobe, pointing *away* from the off-centring direction.

Binning
-------
Solid-angle bins are the faces of a **Goldberg polyhedron**: the dual of a
frequency-``nu`` geodesic subdivision of the icosahedron. It has

    C = 10*nu^2 + 2 cells,  of which exactly 12 are pentagons

and the rest hexagons. The 12 pentagons are not an implementation shortcut and
cannot be removed: a closed surface tiled by hexagons alone would need Euler
characteristic 0 (a torus), so any hexagonal tiling of a *sphere* carries
exactly 12 pentagonal defects, sitting at the vertices of the parent
icosahedron. This is the same construction used by geodesic domes, fullerenes
and H3-style geospatial grids, and it is as close to a uniform hexagonal
tiling as the sphere allows -- far more uniform than the equal-angle
(theta, phi) grid, whose cells collapse to slivers at the poles.

Cells are *not* exactly equal in area (they range over roughly +/-10% at usual
frequencies), so raw counts are not a density. Every count here is divided by
the cell's own exactly-computed solid angle, which removes the pentagon
artefact that otherwise prints the parent icosahedron onto the map.

Cell boundaries are the spherical Voronoi diagram of the cell centres: a
direction is assigned to the cell whose centre it is nearest to, and the
polygon drawn for that cell is the exact region of directions that land in it.
Assignment is O(1) per point, not O(cells): the icosahedral face is found by a
ray/cone test, inverted through the (linear) gnomonic map to a lattice index,
then a greedy walk on the cell-adjacency graph fixes the small distortion the
gnomonic map introduces. Greedy walks on a Delaunay graph are guaranteed to
reach the true nearest site, and the seed is already within a cell or two, so
the walk terminates in a step or two regardless of resolution.

Reading the output
------------------
``density`` integrates to 1 over the sphere. ``enhancement = 4*pi*density`` is
the more useful display quantity: it is dimensionless, **1 everywhere for an
isotropic site**, and reads directly as "this direction is 1.8x more likely
than chance". ``zScore`` is each cell's *local* Gaussian z against the
isotropic null -- not a significance: the map has C cells, so its largest
local z is routinely 3-5 on pure noise, and at the ~0.2-1 expected counts per
cell of a fine tiling the Gaussian reading of a Poisson count is badly
anti-conservative. The calibrated readouts are ``peakSignificance`` (the peak
cell's exact Poisson tail, corrected for the C cells searched) and the other
``...Significance`` fields, each a one-sided normal deviate. Over-binning
itself is guarded by :func:`recommended_frequency`.

Conventions
-----------
Directions are Cartesian, in the same frame as ``pca_kde`` displacements, and
cell centres/polygon vertices are unit vectors in that frame. Orientation-tensor
axes use the same canonical sign and handedness convention as the PCA axes
(shared ``_eigen_decomposition``) so the two are directly comparable.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.special import gammainc, gammaincc, ndtri

# Shared with the PCA engine on purpose: the orientation tensor's axes are meant
# to be compared against a site's PCA axes, so they must use the identical
# sign/handedness canonicalization or the comparison is meaningless.
from .pca_kde import SiteDisplacements, _eigen_decomposition, displacement_cloud

# Resolution bounds. nu = 1 is the dodecahedron (12 cells, all pentagons); 64
# is 40962 cells, far beyond what any RMC configuration can populate and already
# ~1 s to tessellate. Real use sits between 4 and 24.
MIN_FREQUENCY = 1
MAX_FREQUENCY = 64

# Displacement magnitudes below this (Angstrom) carry no usable direction --
# the orientation of a ~0 vector is numerical noise, not physics -- so they are
# dropped even when no explicit cutoff is requested.
NEGLIGIBLE_AMPLITUDE = 1e-9

# One neighbour-averaging pass moves this fraction of each cell's mass out to
# its neighbours. 0.5 is a visible but not destructive smoothing step.
SMOOTHING_ALPHA = 0.5

# Counts per cell that make a map worth looking at: at n per cell the Poisson
# noise on an isotropic map is 1/sqrt(n), so 12 gives ~29% cell-to-cell scatter
# -- enough to see a 1.5x lobe, not enough to invent one. recommended_frequency
# floors to it: the auto resolution averages at least this many per cell.
DEFAULT_TARGET_PER_CELL = 12

# Cells whose enhancement lies within this relative distance of the maximum
# are one tied peak. Symmetry-equivalent cells have mathematically equal solid
# angles that differ only in the last bits (differently in NumPy and JS), so a
# plain argmax would pick the peak by round-off; 1e-9 is far above that
# (~1e-15) and far below any physical difference. The lowest index wins.
PEAK_TIE_RTOL = 1e-9

# Every significance readout is a tail probability p converted to the
# equivalent one-sided standard-normal deviate z (P(Z >= z) = p). Tail
# probabilities are floored here before the conversion, so |z| <= 37.04 and
# the payload never carries an infinity (which JSON cannot represent).
SIGNIFICANCE_TAIL_FLOOR = 1e-300

# The antipodal-asymmetry flag: A exceeds its inversion-symmetric null mean by
# more than this many null standard deviations.
ASYMMETRY_FLAG_SIGMA = 3.0

# Isotropic expectation of orientationAnisotropy = 3 lambda_1 - 1 is
# ISOTROPIC_ANISOTROPY_SCALE / sqrt(N_eff) to leading order: sqrt(N) (T - I/3)
# tends to a traceless GOE matrix with off-diagonal variance 1/15, whose mean
# largest eigenvalue is 3 sqrt(3 / (2 pi)) / sqrt(15), so 3 E[lambda_1] - 1 ->
# 9 / sqrt(10 pi N) = 1.6057 / sqrt(N).
ISOTROPIC_ANISOTROPY_SCALE = 9.0 / np.sqrt(10.0 * np.pi)

WEIGHTS = ("count", "amplitude", "amplitude2")
FRAMES = ("cartesian", "pca")


# --- Goldberg (hex + 12 pentagon) sphere tessellation ------------------------


@dataclass(frozen=True)
class SphereTiling:
    """A frequency-``nu`` Goldberg tiling of the unit sphere.

    ``polygons`` is padded to six vertices; a pentagon repeats its last vertex,
    which contributes a degenerate (zero-area) triangle to the fan and so is
    invisible to every downstream calculation. ``sizes`` gives the true vertex
    count. ``neighbours`` is padded with -1 for the 12 pentagons.
    """

    frequency: int
    centers: np.ndarray  # (C, 3) unit vectors -- the cell centres
    polygons: np.ndarray  # (C, 6, 3) unit vectors, CCW seen from outside
    sizes: np.ndarray  # (C,) 5 or 6
    areas: np.ndarray  # (C,) solid angle in steradian, sums to 4*pi
    neighbors: np.ndarray  # (C, 6) adjacent cells, -1 padded
    # Index of the cell whose centre is exactly -centers[c] (the tiling is
    # centrosymmetric). Powers the antipodal-asymmetry readout: comparing a
    # cell against its antipode is how the u vs -u imbalance -- the odd
    # anharmonicity no ellipsoid can show -- is quantified.
    antipode: np.ndarray  # (C,)
    is_pentagon: np.ndarray  # (C,) bool -- exactly 12 True
    # Lookup tables for O(1) direction -> cell assignment.
    face_inverse: np.ndarray  # (20, 3, 3)
    lattice: np.ndarray  # (20, nu+1, nu+1) -> cell index, -1 outside the triangle

    @property
    def cell_count(self) -> int:
        return int(self.centers.shape[0])


_PHI = (1.0 + 5.0**0.5) / 2.0


def _normalize(vectors: np.ndarray) -> np.ndarray:
    lengths = np.linalg.norm(vectors, axis=-1, keepdims=True)
    return vectors / np.maximum(lengths, 1e-300)


def _icosahedron() -> tuple[np.ndarray, np.ndarray]:
    """12 unit vertices and 20 outward, counter-clockwise faces.

    Faces are *derived* (edges are the vertex pairs at the minimum separation,
    faces are the mutually-adjacent triples) rather than hard-coded, so the
    winding and completeness are properties of the construction instead of a
    table that could be transcribed wrong; the counts are asserted below.
    """
    raw = []
    for s1 in (1.0, -1.0):
        for s2 in (1.0, -1.0):
            raw.append((0.0, s1, s2 * _PHI))
            raw.append((s1, s2 * _PHI, 0.0))
            raw.append((s2 * _PHI, 0.0, s1))
    vertices = np.asarray(raw, dtype=float)

    distances = np.linalg.norm(vertices[:, None, :] - vertices[None, :, :], axis=-1)
    edge_length = distances[distances > 1e-9].min()
    adjacent = np.abs(distances - edge_length) < 1e-9 * max(edge_length, 1.0)

    faces = []
    for a in range(12):
        for b in range(a + 1, 12):
            if not adjacent[a, b]:
                continue
            for c in range(b + 1, 12):
                if adjacent[a, c] and adjacent[b, c]:
                    faces.append((a, b, c))
    if len(faces) != 20:
        raise RuntimeError(f"icosahedron construction produced {len(faces)} faces")

    vertices = _normalize(vertices)
    oriented = []
    for a, b, c in faces:
        normal = np.cross(vertices[b] - vertices[a], vertices[c] - vertices[a])
        # Outward normal must agree with the face centroid's direction.
        if normal @ (vertices[a] + vertices[b] + vertices[c]) < 0:
            oriented.append((a, c, b))
        else:
            oriented.append((a, b, c))
    return vertices, np.asarray(oriented, dtype=int)


def _subdivide(vertices: np.ndarray, faces: np.ndarray, frequency: int):
    """Geodesic vertices, the (20, nu+1, nu+1) lattice map, and the triangles.

    Lattice point ``P(f, i, j) = ((nu-i-j)*A + i*B + j*C) / nu`` normalized, for
    face ``f = (A, B, C)``. Points shared between faces are identified by an
    exact *combinatorial* key (corner id, or edge id plus integer position along
    the edge) rather than by rounding coordinates, so the merge cannot fail on
    a floating-point tie.
    """
    nu = int(frequency)
    lattice = np.full((faces.shape[0], nu + 1, nu + 1), -1, dtype=int)
    index_of: dict[tuple, int] = {}
    points: list[np.ndarray] = []

    def key_for(face_id, corners, i, j, k):
        a, b, c = corners
        if i == nu:
            return ("v", int(b))
        if j == nu:
            return ("v", int(c))
        if k == nu:
            return ("v", int(a))
        if k == 0:  # edge b--c, position j from b
            return _edge_key(b, c, j, nu)
        if i == 0:  # edge a--c, position j from a
            return _edge_key(a, c, j, nu)
        if j == 0:  # edge a--b, position i from a
            return _edge_key(a, b, i, nu)
        return ("f", int(face_id), int(i), int(j))

    for face_id, (ia, ib, ic) in enumerate(faces):
        corner_a, corner_b, corner_c = vertices[ia], vertices[ib], vertices[ic]
        for i in range(nu + 1):
            for j in range(nu + 1 - i):
                k = nu - i - j
                key = key_for(face_id, (ia, ib, ic), i, j, k)
                index = index_of.get(key)
                if index is None:
                    point = (k * corner_a + i * corner_b + j * corner_c) / nu
                    index = len(points)
                    index_of[key] = index
                    points.append(point / np.linalg.norm(point))
                lattice[face_id, i, j] = index

    centers = np.asarray(points, dtype=float)
    expected = 10 * nu * nu + 2
    if centers.shape[0] != expected:
        raise RuntimeError(
            f"geodesic subdivision produced {centers.shape[0]} vertices, expected {expected}"
        )

    # Triangles, all wound counter-clockwise as seen from outside (the (i, j)
    # parameter plane maps A->(0,0), B->(nu,0), C->(0,nu), which preserves the
    # face's own outward orientation).
    upward = [
        (lattice[f, i, j], lattice[f, i + 1, j], lattice[f, i, j + 1])
        for f in range(faces.shape[0])
        for i in range(nu)
        for j in range(nu - i)
    ]
    downward = [
        (lattice[f, i + 1, j], lattice[f, i + 1, j + 1], lattice[f, i, j + 1])
        for f in range(faces.shape[0])
        for i in range(nu - 1)
        for j in range(nu - 1 - i)
    ]
    triangles = np.asarray(upward + downward, dtype=int)
    return centers, lattice, triangles


def _edge_key(u: int, v: int, position: int, nu: int) -> tuple:
    """Orientation-independent key for a point on the icosahedron edge u--v."""
    if int(u) < int(v):
        return ("e", int(u), int(v), int(position))
    return ("e", int(v), int(u), nu - int(position))


def _tangent_basis(normals: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Right-handed (e1, e2, n) frames: e1 x e2 = n, so atan2(.e2, .e1) is CCW."""
    reference = np.zeros_like(normals)
    # Pick the axis least aligned with n, so the projection is well conditioned.
    reference[np.arange(normals.shape[0]), np.argmin(np.abs(normals), axis=1)] = 1.0
    e1 = _normalize(reference - normals * (reference * normals).sum(1, keepdims=True))
    e2 = np.cross(normals, e1)
    return e1, e2


def _angular_order(
    centers: np.ndarray, owners: np.ndarray, points: np.ndarray, keys: np.ndarray
):
    """Sort (owner, point) incidences counter-clockwise about each owner's centre.

    Returns the incidence order plus each entry's rank within its own owner, so
    callers can scatter into a padded (C, 6) table with one vectorized pass.

    The CCW cycle is then rotated to start at the entry with the smallest
    integer ``key`` (unique per owner). The raw atan2 start is not reproducible:
    a neighbour lying on the -e1 ray has an e2 component of +/-1e-17 round-off,
    so one engine sorts it at +pi (last) and the other at -pi (first). The
    combinatorial start makes the exported tables identical across engines.
    """
    e1, e2 = _tangent_basis(centers)
    local = points - centers[owners] * (points * centers[owners]).sum(1, keepdims=True)
    angle = np.arctan2((local * e2[owners]).sum(1), (local * e1[owners]).sum(1))
    order = np.lexsort((angle, owners))
    ordered_owners = owners[order]
    degree = np.bincount(ordered_owners, minlength=centers.shape[0])
    start = np.concatenate([[0], np.cumsum(degree)[:-1]])
    rank = np.arange(ordered_owners.size) - start[ordered_owners]
    ordered_keys = np.asarray(keys)[order]
    smallest = np.full(centers.shape[0], np.iinfo(np.int64).max, dtype=np.int64)
    np.minimum.at(smallest, ordered_owners, ordered_keys)
    first = ordered_keys == smallest[ordered_owners]
    shift = np.zeros(centers.shape[0], dtype=np.int64)
    shift[ordered_owners[first]] = rank[first]
    rank = (rank - shift[ordered_owners]) % np.maximum(degree[ordered_owners], 1)
    return order, rank, degree


def _spherical_polygon_areas(polygons: np.ndarray) -> np.ndarray:
    """Solid angle of each padded polygon, by a signed triangle fan.

    Uses the van Oosterom-Strackee form of Girard's theorem,
    ``Omega = 2*atan2(a.(b x c), 1 + a.b + b.c + c.a)``, which stays accurate for
    the small, near-equilateral cells here where the naive spherical-excess
    formula loses precision. A padded pentagon's repeated vertex makes the last
    fan triangle degenerate: the numerator is exactly 0 and the denominator
    positive, so it contributes exactly 0.
    """
    a = polygons[:, 0, :]
    total = np.zeros(polygons.shape[0])
    for i in range(1, polygons.shape[1] - 1):
        b = polygons[:, i, :]
        c = polygons[:, i + 1, :]
        numerator = (a * np.cross(b, c)).sum(1)
        denominator = 1.0 + (a * b).sum(1) + (b * c).sum(1) + (c * a).sum(1)
        total += 2.0 * np.arctan2(numerator, denominator)
    return total


def _assign(
    centers: np.ndarray,
    neighbors: np.ndarray,
    face_inverse: np.ndarray,
    lattice: np.ndarray,
    frequency: int,
    directions: np.ndarray,
) -> np.ndarray:
    """Nearest cell centre for each unit ``direction``, in O(1) per direction.

    Seeded by inverting the gnomonic map (exact ray/cone test for the
    icosahedral face, then linear barycentric -> lattice index), then corrected
    by a greedy walk on the cell-adjacency graph. The seed and the true nearest
    cell differ only near cell boundaries and never by more than a cell or two,
    and greedy descent on a Delaunay graph provably reaches the nearest site, so
    the walk exits after one or two rounds.
    """
    nu = int(frequency)
    if directions.shape[0] == 0:
        return np.zeros(0, dtype=int)

    # lambda[f, n, :] solves u = l0*A + l1*B + l2*C for face f. The face whose
    # cone contains u is the one with all lambda >= 0; argmax of the minimum
    # picks it without an epsilon and breaks exact edge ties deterministically.
    lam = np.einsum("fab,nb->fna", face_inverse, directions)
    face = np.argmax(lam.min(axis=2), axis=0)
    chosen = lam[face, np.arange(directions.shape[0]), :]
    barycentric = chosen / np.maximum(chosen.sum(axis=1, keepdims=True), 1e-300)

    # Round (k, i, j) = nu * barycentric to integers that still sum to nu, by
    # the largest-remainder rule.
    raw = nu * np.clip(barycentric, 0.0, 1.0)
    floor = np.floor(raw)
    deficit = (nu - floor.sum(axis=1)).astype(int)
    remainder = raw - floor
    priority = np.argsort(-remainder, axis=1, kind="stable")
    counts = np.zeros_like(floor)
    slot = np.arange(3)[None, :]
    np.put_along_axis(counts, priority, (slot < deficit[:, None]).astype(float), axis=1)
    integral = (floor + counts).astype(int)

    current = lattice[face, integral[:, 1], integral[:, 2]]
    # A lattice slot is only -1 outside the triangle, which the sum-to-nu
    # rounding cannot produce; guard anyway so a pathological input degrades to
    # a valid cell rather than an index error.
    current = np.where(current >= 0, current, 0)

    valid = neighbors >= 0
    safe_neighbors = np.where(valid, neighbors, 0)
    for _ in range(8):
        best = (centers[current] * directions).sum(1)
        candidate_dots = np.einsum("nkd,nd->nk", centers[safe_neighbors[current]], directions)
        candidate_dots = np.where(valid[current], candidate_dots, -np.inf)
        pick = np.argmax(candidate_dots, axis=1)
        rows = np.arange(current.size)
        improved = candidate_dots[rows, pick] > best
        if not improved.any():
            break
        current = np.where(improved, safe_neighbors[current, pick], current)
    return current


@lru_cache(maxsize=8)
def goldberg_tiling(frequency: int = 8) -> SphereTiling:
    """Hex-and-12-pentagon tiling of the sphere at the given geodesic frequency.

    Cached: the tiling depends on nothing but ``frequency``, and building one is
    far more expensive than binning into it.
    """
    if not np.isfinite(float(frequency)):
        raise ValueError(f"frequency must lie in [{MIN_FREQUENCY}, {MAX_FREQUENCY}]")
    nu = int(frequency)
    if not MIN_FREQUENCY <= nu <= MAX_FREQUENCY:
        raise ValueError(f"frequency must lie in [{MIN_FREQUENCY}, {MAX_FREQUENCY}]")

    vertices, faces = _icosahedron()
    centers, lattice, triangles = _subdivide(vertices, faces, nu)
    cell_count = centers.shape[0]

    # Dual: each geodesic vertex becomes a cell whose polygon vertices are the
    # circumcentres of its incident triangles. For a CCW triangle of unit
    # vectors the plane normal is equidistant from all three, i.e. exactly the
    # spherical circumcentre, and points outward.
    p0, p1, p2 = centers[triangles[:, 0]], centers[triangles[:, 1]], centers[triangles[:, 2]]
    circumcenters = _normalize(np.cross(p1 - p0, p2 - p0))

    owners = triangles.reshape(-1)  # each triangle contributes to its 3 vertices
    incident = np.repeat(np.arange(triangles.shape[0]), 3)
    order, rank, degree = _angular_order(centers, owners, circumcenters[incident], incident)
    # Exactly 12 pentagons, everything else hexagonal. At nu = 1 the tiling is
    # the dodecahedron: all 12 cells are pentagons and no hexagons exist.
    if int((degree == 5).sum()) != 12 or not np.isin(degree, (5, 6)).all():
        raise RuntimeError("dual tiling is not a Goldberg polyhedron (bad cell degrees)")

    polygons = np.zeros((cell_count, 6, 3))
    polygons[owners[order], rank] = circumcenters[incident[order]]
    # Pad each pentagon by repeating its last real vertex (zero-area fan step).
    pentagons = np.flatnonzero(degree == 5)
    polygons[pentagons, 5] = polygons[pentagons, 4]

    # Directed edges: with every triangle wound CCW, each ordered pair (v, w)
    # occurs exactly once across all triangles, so this is the adjacency list
    # with no dedup needed.
    edge_from = triangles[:, [0, 1, 2]].reshape(-1)
    edge_to = triangles[:, [1, 2, 0]].reshape(-1)
    edge_order, edge_rank, edge_degree = _angular_order(
        centers, edge_from, centers[edge_to], edge_to
    )
    if not np.array_equal(edge_degree, degree):
        raise RuntimeError("cell adjacency disagrees with the dual polygon degrees")
    neighbors = np.full((cell_count, 6), -1, dtype=int)
    neighbors[edge_from[edge_order], edge_rank] = edge_to[edge_order]

    areas = _spherical_polygon_areas(polygons)
    if not np.isclose(areas.sum(), 4.0 * np.pi, rtol=1e-9):
        raise RuntimeError(f"cell areas sum to {areas.sum()}, not 4*pi")

    face_inverse = np.linalg.inv(np.stack([vertices[face].T for face in faces]))

    # The icosahedron is centrosymmetric and subdivision preserves that, so -c
    # is exactly another cell centre; resolving it through the same lookup keeps
    # one code path (and asserts the symmetry).
    antipode = _assign(centers, neighbors, face_inverse, lattice, nu, -centers)
    if np.abs(centers[antipode] + centers).max() > 1e-9:
        raise RuntimeError("tiling is not centrosymmetric; the antipode map would be wrong")

    return SphereTiling(
        frequency=nu,
        centers=centers,
        polygons=polygons,
        sizes=degree.astype(int),
        areas=areas,
        neighbors=neighbors,
        antipode=antipode,
        is_pentagon=(degree == 5),
        face_inverse=face_inverse,
        lattice=lattice,
    )


def _opposite_hemisphere(directions: np.ndarray) -> np.ndarray:
    """True where the first non-zero component is negative.

    For any u != 0 exactly one of u and -u is flagged, so it splits the sphere
    into two antipodal halves (the zero vector counts as unflagged).
    """
    x, y, z = directions[:, 0], directions[:, 1], directions[:, 2]
    return (x < 0) | ((x == 0) & ((y < 0) | ((y == 0) & (z < 0))))


def assign_cells(tiling: SphereTiling, directions: np.ndarray) -> np.ndarray:
    """Cell index for each row of ``directions`` (need not be normalized).

    Exactly inversion-equivariant: ``assign(-u) == antipode[assign(u)]`` for
    every u, including directions that sit exactly on a Voronoi boundary. The
    icosahedron puts 2-fold axes on x, y, z, so the crystal axes (odd nu) and
    <111> (nu not a multiple of 3) are such ties, and a first-maximum
    tie-break applied to +u and -u independently lands them in cells that are
    not antipodes -- an exactly centrosymmetric cloud then reads as maximally
    one-sided. Only one hemisphere is therefore resolved directly; the other
    is mapped through the (exact, involutive) antipode table.
    """
    directions = _normalize(np.atleast_2d(np.asarray(directions, dtype=float)))
    flip = _opposite_hemisphere(directions)
    resolved = _assign(
        tiling.centers,
        tiling.neighbors,
        tiling.face_inverse,
        tiling.lattice,
        tiling.frequency,
        np.where(flip[:, None], -directions, directions),
    )
    return np.where(flip, tiling.antipode[resolved], resolved)


def recommended_frequency(
    n_points: int,
    *,
    target_per_cell: int = DEFAULT_TARGET_PER_CELL,
    max_frequency: int = 24,
) -> int:
    """Largest frequency whose cells still average ``target_per_cell`` points.

    Returns the largest ``nu`` in ``[MIN_FREQUENCY, max_frequency]`` with
    ``n_points / (10*nu**2 + 2) >= target_per_cell`` -- a floor, so the average
    occupancy never drops below the target. The one exception is the bottom
    of the range: fewer than ``12 * target_per_cell`` points still get the
    12-cell dodecahedron (``MIN_FREQUENCY``), the coarsest tiling there is.

    Over-binning is the failure mode of this kind of plot: push the resolution
    past the data and every cell holds 0 or 1 points, the map turns into Poisson
    confetti, and the confetti looks like structure. Callers should use this as
    the default resolution and let the user override it deliberately.
    """
    target = _validated_target(target_per_cell)
    if not np.isfinite(float(max_frequency)) or int(max_frequency) < MIN_FREQUENCY:
        raise ValueError(f"max_frequency must be a finite number >= {MIN_FREQUENCY}")
    if not np.isfinite(float(n_points)):
        raise ValueError("n_points must be a finite number")
    cap = min(int(max_frequency), MAX_FREQUENCY)
    if n_points <= 0:
        return MIN_FREQUENCY
    # Seed from the continuous inverse of C = 10 nu^2 + 2, then settle the
    # boundary with the exact comparison target * C <= N (integer-valued, so
    # it is exact in both engines and immune to sqrt round-off).
    seed = int(np.floor(np.sqrt(max(float(n_points) / target - 2.0, 0.0) / 10.0)))
    frequency = min(max(seed, MIN_FREQUENCY), cap)
    while frequency < cap and target * (10 * (frequency + 1) ** 2 + 2) <= n_points:
        frequency += 1
    while frequency > MIN_FREQUENCY and target * (10 * frequency**2 + 2) > n_points:
        frequency -= 1
    return int(frequency)


def _validated_target(target_per_cell) -> int:
    """``target_per_cell`` as the integer the resolution guard divides by."""
    if not np.isfinite(float(target_per_cell)) or float(target_per_cell) < 1:
        raise ValueError("target_per_cell must be a finite number >= 1")
    return int(target_per_cell)


# --- Orientation histogram ----------------------------------------------------


def _normal_deviate(upper: float, lower: float) -> float:
    """One-sided normal deviate z with P(Z >= z) = ``upper``.

    ``lower`` is ``1 - upper`` computed independently by the caller, so the
    conversion stays accurate in both tails (``1 - upper`` itself would lose
    every digit once ``upper`` is within 1e-16 of 1).
    """
    if upper <= 0.5:
        return float(-ndtri(max(float(upper), SIGNIFICANCE_TAIL_FLOOR)))
    return float(ndtri(max(float(lower), SIGNIFICANCE_TAIL_FLOOR)))


def _peak_significance(count: int, expected: float, trials: int) -> dict:
    """Look-elsewhere-corrected Poisson significance of one cell's raw count.

    The local p-value is the exact Poisson upper tail P(X >= n | e) =
    P(n, e) (regularized lower incomplete gamma), not the Gaussian reading of
    (n - e)/sqrt(e), which is badly anti-conservative at e ~ 0.2-1. The cell
    was chosen as the maximum of ``trials`` cells, so the global p-value is the
    Sidak correction 1 - (1 - p)^C, evaluated in log space.
    """
    if count <= 0:
        local, local_lower = 1.0, 0.0
    else:
        local = float(gammainc(count, expected))
        local_lower = float(gammaincc(count, expected))
    if local < 0.5:
        log_lower = trials * np.log1p(-local)
    elif local_lower > 0.0:
        log_lower = trials * np.log(local_lower)
    else:
        log_lower = -np.inf
    corrected = float(-np.expm1(log_lower))
    corrected_lower = float(np.exp(log_lower))
    return {
        "local": local,
        "corrected": corrected,
        "deviate": _normal_deviate(corrected, corrected_lower),
    }


def _central_binomial(kmax: int) -> np.ndarray:
    """``c_k = C(2k, k) / 4**k`` for k = 0..kmax, by the exact recurrence.

    ``c_k = c_{k-1} (2k - 1) / (2k)``, evaluated as one sequential product so
    the JS port (same loop) reproduces it bit for bit.
    """
    k = np.arange(1, kmax + 1)
    return np.concatenate([[1.0], np.cumprod((2 * k - 1) / (2 * k))])


def _antipodal_asymmetry_test(counts: np.ndarray, antipode: np.ndarray, used: int) -> dict:
    """A = sum_pairs |n+ - n-| / N and its inversion-symmetric null.

    Conditional on each antipodal pair's total T (the data fix it), an
    inversion-symmetric (p(u) = p(-u)) population splits every pair as
    X ~ Bin(T, 1/2), independently across pairs. Then D = |2X - T| has
    E[D] = T * c_{floor(T/2)} with c_k = C(2k, k)/4^k, and E[D^2] = T, so the
    null mean and variance of A are exact sums over pairs -- valid for any
    anisotropic, centrosymmetric site and any count level, and bounded by 1
    like A itself.
    """
    first = np.arange(counts.size) < antipode  # each unordered pair once
    plus = counts[first]
    minus = counts[antipode[first]]
    total = plus + minus
    central = _central_binomial(int(total.max()) // 2 if total.size else 0)
    pair_mean = total * central[total // 2]
    pair_variance = np.maximum(total - pair_mean * pair_mean, 0.0)
    asymmetry = float(np.abs(plus - minus).sum() / used)
    null_mean = float(pair_mean.sum() / used)
    null_sd = float(np.sqrt(pair_variance.sum()) / used)
    if null_sd > 0.0:
        z_value = (asymmetry - null_mean) / null_sd
        significant = bool(z_value > ASYMMETRY_FLAG_SIGMA)
    else:
        # Every occupied pair holds a single atom: A equals its null exactly
        # and carries no information either way.
        z_value, significant = None, False
    return {
        "value": asymmetry,
        "null": null_mean,
        "nullSd": null_sd,
        "z": z_value,
        "significant": significant,
    }


def _smooth(mass: np.ndarray, neighbors: np.ndarray, passes: int) -> np.ndarray:
    """Neighbour diffusion on the cell graph; exactly mass-conserving.

    Each cell keeps ``1 - alpha`` of its mass and splits the rest evenly over
    its own neighbours, so the total is invariant by construction -- the smooth
    map still integrates to the same number of points as the raw one.
    """
    valid = neighbors >= 0
    degree = valid.sum(axis=1)
    source = np.repeat(np.arange(mass.size), degree)
    target = neighbors[valid]
    for _ in range(int(passes)):
        share = SMOOTHING_ALPHA * mass / np.maximum(degree, 1)
        updated = (1.0 - SMOOTHING_ALPHA) * mass
        np.add.at(updated, target, share[source])
        mass = updated
    return mass


def orientation_histogram(
    vectors: np.ndarray,
    *,
    frequency: int | None = None,
    weight: str = "count",
    min_amplitude: float = 0.0,
    min_amplitude_quantile: float = 0.0,
    smoothing: int = 0,
    frame: str = "cartesian",
    geometry: bool = True,
    target_per_cell: int = DEFAULT_TARGET_PER_CELL,
) -> dict:
    """Solid-angle distribution of the directions of ``vectors``.

    ``vectors`` is an ``(N, 3)`` array of displacements (Cartesian Angstrom);
    only their directions are used. The map is never antipodally folded: seeing
    the full, possibly inversion-*asymmetric* distribution is the point of this
    view -- a +u/-u imbalance about the site mean is real physics (skewness:
    odd-order anharmonicity, or unequal occupation of opposite off-centre
    wells) that the ellipsoid's second moment is structurally blind to.
    ``antipodalAsymmetry`` quantifies exactly that imbalance. It is measured
    about the cloud's own mean, so a coherent shift of every copy (ordered
    off-centring) is invisible here by construction -- see the module notes.

    Parameters
    ----------
    frequency:
        Geodesic frequency; the tiling has ``10*frequency**2 + 2`` cells. ``None``
        picks :func:`recommended_frequency` for the surviving point count.
    weight:
        ``"count"`` (the orientation distribution proper -- every atom counts
        once regardless of how far it moved), or ``"amplitude"`` / ``"amplitude2"``
        to weight each direction by ``|dr|`` or ``|dr|**2``. The ``amplitude2``
        map is the angular decomposition of the mean-square displacement, so its
        second moment reproduces the U tensor the ellipsoid is built from.
    min_amplitude, min_amplitude_quantile:
        Drop displacements shorter than an absolute cutoff (Angstrom) and/or
        than the given quantile of ``|dr|``. Near-zero displacements have a
        direction dominated by noise, so including them dilutes any real
        anisotropy toward uniform; a modest quantile cut sharpens a weak
        pattern without inventing one.
    smoothing:
        Neighbour-diffusion passes over the cell graph (mass-conserving).
    frame:
        ``"cartesian"`` keeps the crystal frame. ``"pca"`` rotates directions
        into the cloud's own principal axes first, which puts PC1 on +x and
        makes different sites directly comparable.
    geometry:
        Include cell polygons in the result. They depend only on ``frequency``,
        so a caller that caches the tiling can turn this off.

    Returns
    -------
    dict
        ``density`` integrates to 1 over the sphere; ``enhancement`` is
        ``4*pi*density``, i.e. 1 for an isotropic cloud; ``zScore`` is the
        local (uncorrected, Gaussian-approximation) z of each raw count
        against the isotropic null; ``peakSignificance`` is the calibrated
        peak test (exact Poisson tail of the peak cell's raw count, Sidak-
        corrected over all cells, as a one-sided normal deviate);
        ``antipodalAsymmetry`` is ``sum_pairs |n(u) - n(-u)| / N`` (0 for an
        inversion-symmetric cloud, 1 for a fully one-sided one), with
        ``antipodalAsymmetryNull`` / ``...NullSd`` its exact mean and spread for
        an inversion-symmetric population with the same pair totals, and
        ``antipodalAsymmetrySignificant`` the flag ``A > null + 3 SD``;
        ``cellMeanAmplitude`` is the mean ``|dr|`` (Angstrom) of the atoms that
        moved into each cell (0 for empty cells) -- the radial-relief quantity.

    Raises
    ------
    ValueError
        For a wrong shape, an unknown ``weight``/``frame``, an out-of-range or
        non-finite option, and for any non-finite (NaN/inf) displacement row --
        a corrupt vector is rejected by name, never silently dropped.
    """
    vectors = np.asarray(vectors, dtype=float)
    if vectors.ndim != 2 or vectors.shape[1] != 3:
        raise ValueError("vectors must be a numeric array with shape (N, 3)")
    if weight not in WEIGHTS:
        raise ValueError(f"weight must be one of {WEIGHTS}")
    if frame not in FRAMES:
        raise ValueError(f"frame must be one of {FRAMES}")
    if not 0.0 <= float(min_amplitude_quantile) < 1.0:
        raise ValueError("min_amplitude_quantile must lie in [0, 1)")
    if not np.isfinite(float(min_amplitude)):
        raise ValueError("min_amplitude must be a finite number")
    if not np.isfinite(float(smoothing)) or float(smoothing) < 0:
        raise ValueError("smoothing must be a finite, non-negative number of passes")
    if frequency is not None and not np.isfinite(float(frequency)):
        raise ValueError(f"frequency must lie in [{MIN_FREQUENCY}, {MAX_FREQUENCY}]")
    _validated_target(target_per_cell)
    # A NaN/inf row would otherwise surface as LAPACK's opaque 'Eigenvalues did
    # not converge' from the PCA fit (even in the cartesian frame), and the two
    # engines used to fail differently on it. Reject it by name instead.
    finite_rows = np.isfinite(vectors).all(axis=1)
    if not finite_rows.all():
        bad = np.flatnonzero(~finite_rows)
        raise ValueError(
            f"vectors contain {bad.size} non-finite (NaN/inf) row(s); first at row {int(bad[0])}"
        )

    total_points = int(vectors.shape[0])
    amplitude = np.linalg.norm(vectors, axis=1)

    # The PCA rotation is fitted on every vector, before any amplitude cut, so
    # the frame is the same one the ellipsoid view reports for this site.
    pca_axes = np.eye(3)
    if total_points >= 4:
        centered = vectors - vectors.mean(axis=0)
        _, pca_axes = _eigen_decomposition(np.cov(centered, rowvar=False, bias=False))

    threshold = float(min_amplitude)
    if min_amplitude_quantile > 0 and amplitude.size:
        threshold = max(threshold, float(np.quantile(amplitude, min_amplitude_quantile)))
    keep = amplitude > max(threshold, NEGLIGIBLE_AMPLITUDE)
    used = int(keep.sum())
    if used < 1:
        raise ValueError("no displacement vectors survive the amplitude cutoff")

    directions = vectors[keep] / amplitude[keep, None]
    kept_amplitude = amplitude[keep]
    if frame == "pca":
        directions = directions @ pca_axes.T

    if frequency is None:
        frequency = recommended_frequency(used, target_per_cell=target_per_cell)
    tiling = goldberg_tiling(int(frequency))

    if weight == "amplitude":
        weights = kept_amplitude
    elif weight == "amplitude2":
        weights = kept_amplitude**2
    else:
        weights = np.ones(used)

    cells = assign_cells(tiling, directions)
    cell_count = tiling.cell_count
    counts = np.bincount(cells, minlength=cell_count)
    mass = np.bincount(cells, weights=weights, minlength=cell_count)

    # Mean |dr| of the atoms that moved into each cell -- the radial-relief
    # quantity (independent of the color weighting, so shape and color carry
    # separate information). Smoothing is applied to the numerator and
    # denominator sums, not the ratio, so a smoothed relief stays the mean of
    # the same smoothed population the color shows.
    amplitude_sum = np.bincount(cells, weights=kept_amplitude, minlength=cell_count)
    count_field = counts.astype(float)

    if smoothing > 0:
        mass = _smooth(mass, tiling.neighbors, smoothing)
        amplitude_sum = _smooth(amplitude_sum, tiling.neighbors, smoothing)
        count_field = _smooth(count_field, tiling.neighbors, smoothing)

    cell_mean_amplitude = np.where(
        count_field > 1e-12, amplitude_sum / np.maximum(count_field, 1e-12), 0.0
    )

    total_mass = float(mass.sum())
    if total_mass <= 0:
        raise ValueError("displacement weights sum to zero")
    density = mass / (total_mass * tiling.areas)
    enhancement = density * 4.0 * np.pi

    # Local z against the isotropic null, always from the *raw* counts:
    # smoothing correlates neighbouring cells, so a z computed after it would
    # overstate the evidence. It is a per-cell description, not a test -- the
    # peak's calibrated significance is _peak_significance below.
    expected = used * tiling.areas / (4.0 * np.pi)
    z_score = (counts - expected) / np.sqrt(np.maximum(expected, 1e-12))

    # Antipodal (inversion) asymmetry: the +u vs -u imbalance the ellipsoid is
    # blind to, and the reason this view is never folded --
    # sum_pairs |n(u) - n(-u)| / N, 0 for an inversion-symmetric cloud, 1 for a
    # fully one-sided one -- with its exact null conditional on the observed
    # pair totals (see _antipodal_asymmetry_test).
    asymmetry_test = _antipodal_asymmetry_test(counts, tiling.antipode, used)

    # Orientation tensor T = <u u^T>: I/3 exactly for a uniform sphere, so its
    # eigenvalues (summing to 1) read as the fraction of directional variance on
    # each axis. Computed from the points, not the binned map, so it is
    # resolution independent.
    tensor = (directions * weights[:, None]).T @ directions / weights.sum()
    tensor_eigenvalues, tensor_axes = _eigen_decomposition(tensor)

    # Bingham's test of uniformity on the tensor: S = (15 n / 2) sum_i
    # (lambda_i - 1/3)^2 ~ chi^2_5 for isotropic directions. With weights the
    # sample size is the effective n = (sum w)^2 / sum w^2 (exact for weights
    # independent of direction). Evaluated from the deviatoric part's
    # Frobenius norm, so it does not depend on the eigensolver.
    effective_points = float(weights.sum() ** 2 / np.sum(weights * weights))
    deviator = tensor - np.trace(tensor) / 3.0 * np.eye(3)
    bingham_statistic = float(7.5 * effective_points * np.sum(deviator * deviator))
    bingham_p_value = float(gammaincc(2.5, bingham_statistic / 2.0))
    anisotropy_significance = _normal_deviate(
        bingham_p_value, float(gammainc(2.5, bingham_statistic / 2.0))
    )
    # Tie-tolerant argmax (lowest index within PEAK_TIE_RTOL of the maximum),
    # identical in both engines -- see PEAK_TIE_RTOL.
    tied = np.flatnonzero(enhancement >= enhancement.max() * (1.0 - PEAK_TIE_RTOL))
    peak = int(tied[0])
    peak_test = _peak_significance(int(counts[peak]), float(expected[peak]), cell_count)

    # Whole-map test: Pearson's X^2 = sum z^2 against chi^2 with C - 1 degrees
    # of freedom (the counts are multinomial with N fixed), as a one-sided
    # normal deviate. The RMS of z ("significance") is not a sigma level: its
    # null is 1 +/- 1/sqrt(2C).
    chi_square = float(np.sum(z_score**2))
    degrees_of_freedom = cell_count - 1
    map_p_value = float(gammaincc(degrees_of_freedom / 2.0, chi_square / 2.0))
    map_significance = _normal_deviate(
        map_p_value, float(gammainc(degrees_of_freedom / 2.0, chi_square / 2.0))
    )

    result = {
        "frequency": int(tiling.frequency),
        "cellCount": cell_count,
        "pentagonCount": int(tiling.is_pentagon.sum()),
        "totalPoints": total_points,
        "usedPoints": used,
        "rejectedPoints": total_points - used,
        "amplitudeCutoff": float(max(threshold, 0.0)),
        "weight": weight,
        "smoothing": int(smoothing),
        "frame": frame,
        "pcaAxes": pca_axes.tolist(),
        "centers": tiling.centers.tolist(),
        "areas": tiling.areas.tolist(),
        "sizes": tiling.sizes.tolist(),
        "antipode": tiling.antipode.tolist(),
        "counts": counts.tolist(),
        "mass": mass.tolist(),
        "density": density.tolist(),
        "enhancement": enhancement.tolist(),
        "expected": expected.tolist(),
        "zScore": z_score.tolist(),
        "vmin": float(enhancement.min()),
        "vmax": float(enhancement.max()),
        "meanCount": float(used) / cell_count,
        "emptyFraction": float((counts == 0).mean()),
        "meanAmplitude": float(kept_amplitude.mean()),
        "rmsAmplitude": float(np.sqrt((kept_amplitude**2).mean())),
        "cellMeanAmplitude": cell_mean_amplitude.tolist(),
        "antipodalAsymmetry": asymmetry_test["value"],
        # Exact mean and standard deviation of A for an inversion-symmetric
        # population with the observed antipodal-pair totals; the z-score
        # (None when the spread is 0) and the flag A > null + 3 SD.
        "antipodalAsymmetryNull": asymmetry_test["null"],
        "antipodalAsymmetryNullSd": asymmetry_test["nullSd"],
        "antipodalAsymmetryZ": asymmetry_test["z"],
        "antipodalAsymmetrySignificant": asymmetry_test["significant"],
        "orientationTensor": tensor.tolist(),
        "orientationEigenvalues": tensor_eigenvalues.tolist(),
        "orientationAxes": tensor_axes.tolist(),
        # 3*t1 - 1: 2 for a perfect single axis; it tends to 0 for a uniform
        # sphere only as N grows (finite samples are biased upward, see
        # orientationAnisotropyNull). A scalar summary of how much the cloud
        # prefers one axis at all.
        "orientationAnisotropy": float(3.0 * tensor_eigenvalues[0] - 1.0),
        # What an isotropic site of this (effective) size reads, and Bingham's
        # calibrated test of isotropy (statistic, p-value, normal deviate).
        "orientationEffectivePoints": effective_points,
        "orientationAnisotropyNull": float(ISOTROPIC_ANISOTROPY_SCALE / np.sqrt(effective_points)),
        "orientationBinghamStatistic": bingham_statistic,
        "orientationBinghamPValue": bingham_p_value,
        "orientationAnisotropySignificance": anisotropy_significance,
        "peakCell": peak,
        # Number of cells tied for the maximum (1 = a unique peak).
        "peakTieCount": int(tied.size),
        "peakDirection": tiling.centers[peak].tolist(),
        "peakEnhancement": float(enhancement[peak]),
        # Local, uncorrected Gaussian z of the peak cell's raw count -- kept
        # for API compatibility; NOT a significance (see peakSignificance).
        "peakZScore": float(z_score[peak]),
        # The peak test: raw count and isotropic expectation of the peak cell,
        # its exact Poisson upper-tail p-value, that p Sidak-corrected for the
        # C cells the maximum was searched over, and the corrected p as a
        # one-sided normal deviate (sigma). This is the calibrated readout.
        "peakCount": int(counts[peak]),
        "peakExpected": float(expected[peak]),
        "peakLocalPValue": peak_test["local"],
        "peakPValue": peak_test["corrected"],
        "peakSignificance": peak_test["deviate"],
        # Legacy whole-map summary: the RMS of the local z (1 +/- 1/sqrt(2C)
        # for pure noise). NOT a sigma level -- use mapSignificance.
        "significance": float(np.sqrt(np.mean(z_score**2))),
        # Pearson's chi-square test of isotropy over all cells.
        "mapChiSquare": chi_square,
        "mapDegreesOfFreedom": int(degrees_of_freedom),
        "mapPValue": map_p_value,
        "mapSignificance": map_significance,
        "recommendedFrequency": recommended_frequency(used, target_per_cell=target_per_cell),
    }
    if geometry:
        result["polygons"] = [
            tiling.polygons[index, : tiling.sizes[index]].tolist() for index in range(cell_count)
        ]
        result["neighbors"] = tiling.neighbors.tolist()
    return result


def site_orientation_histogram(
    sites: SiteDisplacements,
    *,
    reference_number: int | None = None,
    element: str | None = None,
    **kwargs,
) -> dict:
    """:func:`orientation_histogram` for one site (or one element's pooled sites)."""
    cloud = displacement_cloud(sites, reference_number=reference_number, element=element)
    result = orientation_histogram(cloud, **kwargs)
    if reference_number is not None:
        position = sites.site_position(reference_number)
        result["referenceNumber"] = int(reference_number)
        result["element"] = sites.elements[position]
        result["siteFractional"] = sites.site_fractional[position].tolist()
    elif element not in (None, "", "all"):
        result["element"] = str(element)
    return result
