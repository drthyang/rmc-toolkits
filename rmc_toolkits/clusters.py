# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Cluster distortion: the shape of finite atom clusters in an RMC configuration.

A *cluster* is a finite group of one element's atoms bonded within a length
window in the AVERAGE structure -- the Ta4 tetrahedra of GaTa4Se8, a V-V
dimer, an Mo3 trimer. The templates come from the average site positions, so a
cluster in the configuration is always the same atoms (one copy of each member
site, in the right box copy), however far its bonds stretch.

For every cluster the engine measures its edges (the template's bonds) in the
configuration and asks one question: are the clusters distorted beyond what the
atoms' own displacements make by chance? The comparison is a **shuffle null**:
each member keeps the displacement distribution of its own site, but takes it
from a random other cluster, which keeps every site's spread, anisotropy and
non-Gaussianity and removes only the correlation between the atoms of one
cluster. RMC moves atoms one at a time, so without data a configuration looks
like this null; a departure is a correlation the data put in. Each statistic is
reported as observed, null mean and sd over ``shuffles`` shuffles, and
z = (observed - null mean) / null sd.

That null also removes collective motion, and a collective size fluctuation
(breathing: every member moving out or in together) then lowers the shape
statistics, because the null spreads the same radial displacements over
independent atoms: breathing triangles came out "28 % LESS distorted", and the
gyration-axis test read breathing tetrahedra as a shape preference (the
2026-10-14 review's L2). So a second null, ``shapeNull``, keeps each cluster's
own breathing -- its least-squares scale about the template centroid -- and
shuffles only the breathing-free displacements, with the same permutations.
The shape statistics are read against ``shapeNull``, the size statistics
(``meanEdge``, ``breathingSd``) against the independent-atom null.

Statistics, pooled over the templates of the dominant kind (members, edges):

* ``meanEdge``, ``ladder`` (mean of each cluster's k-th shortest edge),
  ``ladderWidth`` (mean longest - shortest edge), ``shapeRms`` (rms of the
  edges about their own cluster's mean: shape distortion), ``breathingSd``
  (sd of each cluster's mean edge: size fluctuation);
* tetrahedra (4 members, 6 edges): the symmetry-adapted edge modes of Td --
  ``eRms`` from the three opposite-pair sums (tetragonal-type, D2d) and
  ``t2Rms`` from the opposite-pair differences, whose direction
  ``t2Trigonality`` = mean(sum d^4 / (sum d^2)^2) tells the patterns apart:
  1/3 when all three opposite pairs differ alike (three edges at one vertex
  long or short: trigonal, C3v), 1 when one pair alone does (one edge
  stretched against its opposite);
* 3+ members: ``longAxis`` / ``shortAxis``, the fraction of clusters whose
  gyration-tensor axis lies within 15 deg of a <100>, <110> or <111> lattice
  direction (for a dimer, the bond itself).

The average-structure geometry is reported too (``ideal``): the template's own
edges and each edge's sampling error from the two sites' displacement spreads
along it -- a static distortion shows there, not in the shuffle test, which
keeps the site means fixed.

Outliers: ``brokenEdgeFraction`` (template edges whose length in the
configuration falls outside the window) and ``extraContacts`` (pairs of the
element within the window that are not template bonds: atoms crossing the gap
between clusters).

The engine also decides what these statistics show (``verdict``, always
present): each question -- distorted, less distorted, breathing, rigid, a
static distortion of the average cluster, split member sites, a trigonal
lean, favoured gyration axes, and what the window did -- is a *Finding*
``{answer, because, evidence}`` read by fixed rules (``CLUSTER_RULES``) from
the unrounded statistics, so every caller (the Copilot's ``analyze_clusters``,
the validation's arm A+, a script) reads the same answer. See the verdict
section below.

The browser port is ``workers/clusters.js``; keep the two in sync (the parity
fixture ``tests/generate_clusters_fixture.py`` pins them, verdicts included).
Shuffles come from a counter-based hash, identical in both engines.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from math import ceil

import numpy as np

from rmc_toolkits.pca_kde import SiteDisplacements, load_site_displacements, site_ellipsoids
from rmc_toolkits.triplets import WINDOW_TOL, _neighbor_bonds, _perpendicular_widths

#: Shuffles of the null by default: enough for a null sd within ~15 %.
DEFAULT_SHUFFLES = 24
#: Seed of the counter-based shuffle hash (both engines).
SHUFFLE_SEED = 0x5EED
#: A gyration axis within this angle of a lattice direction counts for its family.
AXIS_CONE_DEG = 15.0
#: Opposite edge pairs of a tetrahedron, as indices into its sorted edge list
#: (0,1) (0,2) (0,3) (1,2) (1,3) (2,3).
TETRAHEDRON_OPPOSITE = ((0, 5), (1, 4), (2, 3))
FAMILIES = ("100", "110", "111")
#: The size statistics: read against the independent-atom null only (the
#: shape null keeps every cluster's own size).
SIZE_STATISTICS = ("meanEdge", "breathingSd")

_U32 = np.uint32


def mix32(x: np.ndarray) -> np.ndarray:
    """The lowbias32 integer hash (C. Wellons), on uint32 arrays (wrapping)."""
    x = np.asarray(x, dtype=_U32)
    with np.errstate(over="ignore"):  # uint32 products wrap modulo 2**32, as intended
        x = (x ^ (x >> _U32(16))) * _U32(0x7FEB352D)
        x = (x ^ (x >> _U32(15))) * _U32(0x846CA68B)
    return x ^ (x >> _U32(16))


def shuffle_permutation(n: int, shuffle: int, template: int, member: int, seed: int = SHUFFLE_SEED) -> np.ndarray:
    """A permutation of range(n) from a counter-based hash (stable on ties).

    key(c) = mix32(mix32(mix32(mix32(seed ^ mix32(shuffle + 1)) ^ (template + 1)) ^ (member + 1)) ^ c)
    -- the same keys and the same stable sort in workers/clusters.js.
    """
    head = mix32(_U32(seed) ^ mix32(_U32(shuffle + 1)))
    head = mix32(head ^ _U32(template + 1))
    head = mix32(head ^ _U32(member + 1))
    keys = mix32(head ^ np.arange(n, dtype=_U32))
    return np.argsort(keys, kind="stable")


def _rows_times(rows: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    """``rows @ matrix`` by einsum: these products are tiny, and Accelerate's
    BLAS (the macOS numpy build) raises spurious floating-point warnings on them."""
    rows = np.asarray(rows, dtype=float)
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim == 1:
        return np.einsum("...i,i->...", rows, matrix)
    return np.einsum("...i,ij->...j", rows, matrix)


def lattice_families(unit_vectors: np.ndarray) -> dict[str, np.ndarray]:
    """Unit vectors of the <100>, <110> and <111> lattice directions (Cartesian)."""
    a, b, c = (np.asarray(v, dtype=float) for v in unit_vectors)
    groups = {
        "100": [a, b, c],
        "110": [a + b, a - b, a + c, a - c, b + c, b - c],
        "111": [a + b + c, a + b - c, a - b + c, -a + b + c],
    }
    return {name: np.array([v / np.linalg.norm(v) for v in vectors]) for name, vectors in groups.items()}


@dataclass
class _Template:
    sites: list[int]  # site positions (indices into the SiteDisplacements arrays)
    offsets: np.ndarray  # (M, 3) int unit-cell offsets of the members
    edges: list[tuple[int, int]]  # member index pairs, sorted


def _cell_positions(sites: SiteDisplacements) -> np.ndarray:
    """Average site positions consistent with the atoms' cell indices
    (``site_cell_fractional``, unwrapped): a cluster pairs the atoms of given
    cells, so a site at fraction 0 must not wrap to 0.9999 and pull its
    partners from the next cell. Falls back to ``site_fractional``."""
    position = sites.site_cell_fractional
    if position is None:
        position = sites.site_fractional
    return np.asarray(position, dtype=float)


def _templates(
    sites: SiteDisplacements, members: list[int], window: tuple[float, float]
) -> tuple[list[_Template], bool, int]:
    """Finite clusters of the given sites under the window, from the average structure.

    Returns the templates, whether the bonds percolate (an extended network:
    a site bonded to its own image, or two paths reaching one site at
    different cell offsets), and how many sites have no bond in the window.
    """
    unit = np.asarray(sites.unit_vectors, dtype=float)
    frac = _cell_positions(sites)
    lo, hi = window
    reach = np.ceil(hi / _perpendicular_widths(unit)).astype(int) + 1
    shifts = np.array(
        [(i, j, k) for i in range(-reach[0], reach[0] + 1)
         for j in range(-reach[1], reach[1] + 1) for k in range(-reach[2], reach[2] + 1)],
        dtype=int,
    )
    adjacency: dict[int, list[tuple[int, np.ndarray]]] = {s: [] for s in members}
    extended = False
    for s in members:
        for t in members:
            lengths = np.linalg.norm(_rows_times((frac[t] + shifts - frac[s]), unit), axis=1)
            for index in np.flatnonzero((lengths >= max(lo - WINDOW_TOL, 0.0)) & (lengths <= hi + WINDOW_TOL)):
                shift = shifts[index]
                if s == t:
                    if shift.any():
                        extended = True  # bonded to its own image: a chain or net
                    continue
                adjacency[s].append((t, shift))
    templates: list[_Template] = []
    seen: set[int] = set()
    isolated = 0
    for anchor in members:
        if anchor in seen:
            continue
        if not adjacency[anchor]:
            seen.add(anchor)
            isolated += 1
            continue
        offsets = {anchor: np.zeros(3, dtype=int)}
        queue = [anchor]
        while queue:
            s = queue.pop(0)
            for t, shift in adjacency[s]:
                reached = offsets[s] + shift
                if t in offsets:
                    if not np.array_equal(offsets[t], reached):
                        extended = True
                    continue
                offsets[t] = reached
                queue.append(t)
        order = sorted(offsets)
        seen.update(order)
        index = {s: i for i, s in enumerate(order)}
        edges = sorted({
            tuple(sorted((index[s], index[t])))
            for s in order for t, shift in adjacency[s]
            if t in index and np.array_equal(offsets[t] - offsets[s], shift)
        })
        templates.append(_Template(order, np.array([offsets[s] for s in order]), edges))
    return templates, extended, isolated


def _cell_linear(cells: np.ndarray, supercell: np.ndarray) -> np.ndarray:
    cells = np.mod(cells, supercell)
    return (cells[:, 0] * supercell[1] + cells[:, 1]) * supercell[2] + cells[:, 2]


def _edge_vectors(positions: np.ndarray, edges) -> np.ndarray:
    """(C, E, 3) edge vectors of (C, M, 3) clusters. ``edges`` holds member
    index pairs: one list for every cluster (E, 2), or one per cluster
    (C, E, 2) -- templates of one kind need not share their edges (two chains
    of three, one bonded 0-1-2 and one 1-0-2)."""
    pairs = np.asarray(edges, dtype=int)
    if pairs.ndim == 2:
        pairs = np.broadcast_to(pairs, (positions.shape[0],) + pairs.shape)
    cluster = np.arange(positions.shape[0])[:, None]
    return positions[cluster, pairs[..., 1], :] - positions[cluster, pairs[..., 0], :]


def _statistics(positions: np.ndarray, edges, families: dict[str, np.ndarray]) -> dict:
    """Shape statistics of clusters given as (C, M, 3) member positions and
    their edges (see :func:`_edge_vectors`)."""
    vectors = _edge_vectors(positions, edges)
    lengths = np.sqrt(vectors[..., 0] ** 2 + vectors[..., 1] ** 2 + vectors[..., 2] ** 2)
    ladder = np.sort(lengths, axis=1)
    means = lengths.mean(axis=1)
    stats = {
        "meanEdge": float(lengths.mean()),
        "ladder": ladder.mean(axis=0).tolist(),
        "ladderWidth": float((ladder[:, -1] - ladder[:, 0]).mean()),
        "shapeRms": float(np.sqrt(((lengths - means[:, None]) ** 2).mean())),
        "breathingSd": float(means.std()),
    }
    members = positions.shape[1]
    if members == 4 and vectors.shape[1] == 6:
        sums = np.stack([lengths[:, a] + lengths[:, b] for a, b in TETRAHEDRON_OPPOSITE], axis=1)
        diffs = np.stack([lengths[:, a] - lengths[:, b] for a, b in TETRAHEDRON_OPPOSITE], axis=1)
        e_modes = sums - sums.mean(axis=1, keepdims=True)
        stats["eRms"] = float(np.sqrt((e_modes ** 2).sum(axis=1).mean() / 2.0))
        stats["t2Rms"] = float(np.sqrt((diffs ** 2).sum(axis=1).mean() / 3.0))
        norm2 = (diffs ** 2).sum(axis=1)
        valid = norm2 > 0
        stats["t2Trigonality"] = (
            float(((diffs[valid] ** 4).sum(axis=1) / norm2[valid] ** 2).mean()) if valid.any() else None
        )
    if members >= 2:
        if members == 2:
            axes = {"longAxis": vectors[:, 0, :]}
        else:
            centred = positions - positions.mean(axis=1, keepdims=True)
            gyration = np.einsum("cmi,cmj->cij", centred, centred) / members
            _, eigenvectors = np.linalg.eigh(gyration)
            axes = {"longAxis": eigenvectors[:, :, 2], "shortAxis": eigenvectors[:, :, 0]}
        cone = np.cos(np.radians(AXIS_CONE_DEG))
        for name, axis in axes.items():
            unit = axis / np.linalg.norm(axis, axis=1, keepdims=True)
            stats[name] = {
                family: float((np.abs(_rows_times(unit, directions.T)).max(axis=1) >= cone).mean())
                for family, directions in families.items()
            }
    return stats


def _null_summary(observed: dict, nulls: list[dict]) -> tuple[dict, dict, dict]:
    """Mean, sd and z of every statistic, elementwise for lists and dicts."""
    def combine(key_path, values, obs):
        array = np.array([v for v in values if v is not None], dtype=float)
        if obs is None or array.size < 2:
            return None, None, None
        mean = float(array.mean())
        sd = float(array.std(ddof=1))
        z = (float(obs) - mean) / sd if sd > 0 else None
        return mean, sd, z

    means, sds, zs = {}, {}, {}
    for key, obs in observed.items():
        if isinstance(obs, list):
            parts = [combine(key, [n[key][k] for n in nulls], obs[k]) for k in range(len(obs))]
            means[key], sds[key], zs[key] = ([p[i] for p in parts] for i in range(3))
        elif isinstance(obs, dict):
            parts = {k: combine(key, [n[key][k] for n in nulls], obs[k]) for k in obs}
            means[key] = {k: p[0] for k, p in parts.items()}
            sds[key] = {k: p[1] for k, p in parts.items()}
            zs[key] = {k: p[2] for k, p in parts.items()}
        else:
            means[key], sds[key], zs[key] = combine(key, [n.get(key) for n in nulls], obs)
    return means, sds, zs


def cluster_distortion(
    sites: SiteDisplacements,
    element: str,
    window: tuple[float, float],
    *,
    shuffles: int = DEFAULT_SHUFFLES,
    seed: int = SHUFFLE_SEED,
) -> dict:
    """Cluster templates, observed shape statistics, their shuffle null, and the verdict.

    ``window`` is the inclusive bond-length window (Angstrom) that links the
    element's atoms into clusters in the average structure. Raises
    ``ValueError`` for an unknown element, a malformed window, or a file
    without cell indices (a coordinates-only configuration). ``verdict`` is
    always returned: for an extended network or no clusters every shape
    finding is None (see :func:`cluster_verdict`).
    """
    lo, hi = (float(window[0]), float(window[1]))
    if not (np.isfinite(lo) and np.isfinite(hi)) or lo < 0 or hi <= lo:
        raise ValueError(f"window must be 0 <= min < max, got [{window[0]}, {window[1]}]")
    if sites.cell_indices is None:
        raise ValueError("cluster analysis needs the reference-site and cell columns of the .rmc6f file")
    if shuffles < 2:
        raise ValueError(f"shuffles must be at least 2, got {shuffles}")
    element = element.strip().capitalize()
    members = [s for s, label in enumerate(sites.elements) if label == element]
    if not members:
        raise ValueError(f"no site holds {element}; elements: {', '.join(sites.species)}")

    unit = np.asarray(sites.unit_vectors, dtype=float)
    supercell = np.asarray(sites.supercell, dtype=int)
    frac = _cell_positions(sites)
    families = lattice_families(unit)
    n_cells = int(np.prod(supercell))
    row_of = -np.ones((len(sites.reference_numbers), n_cells), dtype=int)
    row_of[sites.site_index, _cell_linear(np.asarray(sites.cell_indices, dtype=int), supercell)] = np.arange(
        len(sites.site_index)
    )
    cell_grid = np.stack(np.unravel_index(np.arange(n_cells), tuple(supercell)), axis=1)

    templates, extended, isolated = _templates(sites, members, (lo, hi))
    result: dict = {
        "element": element,
        "window": [lo, hi],
        "shuffles": int(shuffles),
        "sites": [int(sites.reference_numbers[s]) for s in members],
        "extended": bool(extended),
        "isolatedSites": int(isolated),
        "templates": [],
        "kinds": [],
        "parseWarning": sites.parse_warning,
    }
    if extended or not templates:
        result["verdict"] = _no_cluster_verdict("extended_network" if extended else "no_clusters")
        return result

    per_template = []
    rows_all = []  # per template: (M, cells) rows, -1 where the atom is missing
    for t_index, template in enumerate(templates):
        rows = np.stack([
            row_of[s, _cell_linear(cell_grid + template.offsets[m], supercell)]
            for m, s in enumerate(template.sites)
        ])  # (M, C)
        rows_all.append(rows)
        complete = (rows >= 0).all(axis=0)
        rows = rows[:, complete]
        base = np.stack([_rows_times(frac[s] + template.offsets[m], unit) for m, s in enumerate(template.sites)])
        ideal = [float(np.linalg.norm(base[j] - base[i])) for i, j in template.edges]
        edge_se = []
        for i, j in template.edges:
            direction = (base[j] - base[i]) / np.linalg.norm(base[j] - base[i])
            variance = 0.0
            for s in (template.sites[i], template.sites[j]):
                cloud = _rows_times(sites.displacements[sites.site_index == s], direction)
                variance += float((cloud ** 2).mean()) / cloud.size
            edge_se.append(float(np.sqrt(variance)))
        result["templates"].append({
            "referenceNumbers": [int(sites.reference_numbers[s]) for s in template.sites],
            "offsets": template.offsets.tolist(),
            "edges": [list(edge) for edge in template.edges],
            "idealEdges": ideal,
            "edgeSe": edge_se,
            "clusters": int(complete.sum()),
            "incomplete": int((~complete).sum()),
        })
        per_template.append((t_index, template, rows, base))

    kinds: dict[tuple[int, int], int] = {}
    for _, template, rows, _ in per_template:
        key = (len(template.sites), len(template.edges))
        kinds[key] = kinds.get(key, 0) + rows.shape[1]
    result["kinds"] = [
        # Most clusters first, ties to the smaller kind: the dominant kind leads.
        {"members": m, "edges": e, "clusters": n} for (m, e), n in sorted(kinds.items(), key=lambda kv: (-kv[1], kv[0][0]))
    ]
    dominant = max(kinds, key=lambda key: (kinds[key], -key[0]))
    chosen = [entry for entry in per_template if (len(entry[1].sites), len(entry[1].edges)) == dominant]
    result["kind"] = {"members": dominant[0], "edges": dominant[1], "clusters": kinds[dominant], "templates": len(chosen)}

    # Each cluster is measured with its own template's edges.
    edges = np.concatenate([
        np.broadcast_to(np.asarray(template.edges, dtype=int), (rows.shape[1], len(template.edges), 2))
        for _, template, rows, _ in chosen
    ])

    def positions(shuffle: int | None) -> np.ndarray:
        blocks = []
        for t_index, template, rows, base in chosen:
            n = rows.shape[1]
            picked = []
            for m in range(len(template.sites)):
                order = np.arange(n) if shuffle is None else shuffle_permutation(n, shuffle, t_index, m, seed)
                picked.append(base[m] + sites.displacements[rows[m, order]])
            blocks.append(np.stack(picked, axis=1))
        return np.concatenate(blocks, axis=0)

    # Each cluster's breathing: its least-squares scale about the template
    # centroid, eps_c = sum_m d_mc . r_m / sum_m |r_m|^2 (r_m = base_m - centroid),
    # and the displacements with it taken out, d'_mc = d_mc - eps_c r_m. The
    # fit also takes out a share beta_m = |r_m|^2 / sum_k |r_k|^2 of each
    # member's own radial noise: for members of equal radial variance v the
    # residual's radial variance is v (1 - beta_m), whatever the geometry, so
    # its radial component is scaled back by 1 / sqrt(1 - beta_m). Without it
    # the shuffled shapes came out ~6 % too regular for independent atoms.
    breathing = []
    for _, template, rows, base in chosen:
        radial = base - base.mean(axis=0)
        displaced = np.stack([sites.displacements[rows[m]] for m in range(len(template.sites))], axis=1)  # (C, M, 3)
        lengths2 = (radial ** 2).sum(axis=1)
        scale = np.einsum("cmi,mi->c", displaced, radial) / float(lengths2.sum())
        residual = displaced - scale[:, None, None] * radial[None, :, :]
        unit = np.zeros_like(radial)
        np.divide(radial, np.sqrt(lengths2)[:, None], out=unit, where=lengths2[:, None] > 0)
        boost = 1.0 / np.sqrt(1.0 - lengths2 / lengths2.sum()) - 1.0
        along = np.einsum("cmi,mi->cm", residual, unit)
        residual = residual + (boost[None, :] * along)[:, :, None] * unit[None, :, :]
        breathing.append((radial, scale, residual))

    def shape_positions(shuffle: int) -> np.ndarray:
        """The shape null: every cluster keeps exactly its own breathing, and
        its members' breathing-free displacements come from random other
        clusters (the same permutations as the independent-atom null), with
        the breathing their combination happens to carry taken out again."""
        blocks = []
        for (t_index, template, rows, base), (radial, scale, residual) in zip(chosen, breathing):
            n = rows.shape[1]
            mixed = np.stack([
                residual[shuffle_permutation(n, shuffle, t_index, m, seed), m] for m in range(len(template.sites))
            ], axis=1)  # (C, M, 3)
            stray = np.einsum("cmi,mi->c", mixed, radial) / float((radial ** 2).sum())
            blocks.append(base[None, :, :] + mixed + (scale - stray)[:, None, None] * radial[None, :, :])
        return np.concatenate(blocks, axis=0)

    observed_positions = positions(None)
    observed = _statistics(observed_positions, edges, families)
    nulls = [_statistics(positions(k), edges, families) for k in range(shuffles)]
    null_mean, null_sd, z = _null_summary(observed, nulls)
    shape_mean, shape_sd, shape_z = _null_summary(
        observed, [_statistics(shape_positions(k), edges, families) for k in range(shuffles)]
    )
    # The shape null holds each cluster's size by construction: its size
    # statistics are no test, and are left out.
    for summary in (shape_mean, shape_sd, shape_z):
        for key in SIZE_STATISTICS:
            summary.pop(key, None)
    result.update({
        "observed": observed, "nullMean": null_mean, "nullSd": null_sd, "z": z,
        "shapeNull": {"mean": shape_mean, "sd": shape_sd, "z": shape_z},
    })

    ideal_sorted = np.mean([np.sort(result["templates"][t]["idealEdges"]) for t, *_ in chosen], axis=0)
    spreads = [np.ptp(result["templates"][t]["idealEdges"]) for t, *_ in chosen]
    errors = [float(np.median(result["templates"][t]["edgeSe"])) for t, *_ in chosen]
    result["ideal"] = {
        "ladder": ideal_sorted.tolist(),
        "spread": float(np.max(spreads)),
        "edgeSe": float(np.median(errors)),
        "spreadOverSe": float(np.max(np.array(spreads) / np.maximum(errors, 1e-12))),
    }

    # Outliers: template edges outside the window, and contacts in the window
    # that are not template bonds (counted over the whole element).
    lengths = np.linalg.norm(_edge_vectors(observed_positions, edges), axis=2)
    inside = (lengths >= max(lo - WINDOW_TOL, 0.0)) & (lengths <= hi + WINDOW_TOL)
    lattice = np.asarray(sites.lattice_vectors, dtype=float)
    element_rows = np.flatnonzero(np.isin(sites.site_index, members))
    fractional = np.mod(
        (frac[sites.site_index[element_rows]] + sites.cell_indices[element_rows]) / supercell
        + _rows_times(sites.displacements[element_rows], np.linalg.inv(lattice)),
        1.0,
    )
    rows_local = np.arange(element_rows.size)
    directed = _neighbor_bonds(
        fractional, rows_local, fractional, rows_local, lattice, [(lo, hi)], count_only=True
    ).sum()
    template_bonds_inside = 0
    # Every template bond whose two atoms are there, in incomplete clusters
    # too: a cluster missing an atom still has its other bonds, and they are
    # not contacts between clusters.
    for (t_index, template, _, base), rows in zip(per_template, rows_all):
        for i, j in template.edges:
            present = (rows[i] >= 0) & (rows[j] >= 0)
            span = (base[j] + sites.displacements[rows[j, present]]) - (base[i] + sites.displacements[rows[i, present]])
            tl = np.linalg.norm(span, axis=1)
            template_bonds_inside += int(((tl >= max(lo - WINDOW_TOL, 0.0)) & (tl <= hi + WINDOW_TOL)).sum())
    extra = int(directed) // 2 - template_bonds_inside
    result["outliers"] = {
        "brokenEdgeFraction": float(1.0 - inside.mean()),
        "extraContacts": extra,
        "extraContactsPerAtom": float(2 * extra / element_rows.size) if element_rows.size else 0.0,
    }

    # The analysed kind's member sites, from the clouds this engine holds:
    # their site rows decide which are split (one rule, is_split_site).
    member_references = sorted({int(sites.reference_numbers[s]) for _, template, _, _ in chosen for s in template.sites})
    rows = {row["referenceNumber"]: row for row in site_ellipsoids(sites)}
    result["verdict"] = cluster_verdict(result, [rows[reference] for reference in member_references])
    return result


def cluster_distortion_from_file(path: str, element: str, rmin: float, rmax: float, shuffles: int = DEFAULT_SHUFFLES) -> dict:
    """:func:`cluster_distortion` of an ``.rmc6f`` file, read afresh."""
    return cluster_distortion(load_site_displacements(path), element, (rmin, rmax), shuffles=shuffles)


# --- the verdict ------------------------------------------------------------------
#
# The engine decides what its statistics show, so every caller reads one answer:
# the Copilot's analyze_clusters only renders it, and arm A+ maps it to claims.
# Each question is a Finding,
#
#     {answer, because, evidence}
#
# answer True/False, or None where the question does not apply or its statistic
# is undefined (never "not significant"); because, a code from CLUSTER_BECAUSE
# saying why; evidence, every statistic the answer consulted (a negative keeps
# them all), those that clear their bar first, then the rest in a fixed order.
# An Evidence record is
#
#     {statistic, against, observed, null, z, change, changeKind, bar, clears,
#      side?, strict?, site?}
#
# against: the null it is read against ('shapeNull', 'independentNull'), the
# sampling error ('samplingError'), a site's kurtosis ('siteKurtosis') or the
# window ('window'); change = observed / null - 1 ('ratio') or observed - null
# ('difference'); bar maps a field of the record (z, change, observed) to its
# threshold, side says which way it must be cleared ('above': >= bar, 'below':
# <= -bar, absent: |value| >= bar; strict: > instead of >=), and clears is
# evidence_clears(record): recomputable from the record alone. Decisions use the
# unrounded statistics; every number is finite or None (Flask's strict JSON).
# The rules are the Copilot's (they were its summarizeClusters' until
# 2026-10-19), ported rule for rule, and echoed in the verdict (CLUSTER_RULES).

#: A cluster statistic departs from its null only when both its z and its
#: effect size clear these: with thousands of clusters a 1 % effect reaches
#: z > 3 and still means nothing physically.
CLUSTER_Z = 3
#: Relative change (observed / null - 1) of a width or an rms.
CLUSTER_EFFECT = 0.05
#: Change (observed - null) of the tetrahedral T2 trigonality (range 1/3 to 1).
TRIGONALITY_EFFECT = 0.02
#: Observed / null fraction of clusters with a gyration axis near one lattice
#: family (one-sided: a family more clusters favour than by chance).
AXIS_ENHANCEMENT = 1.1
#: The average cluster is statically distorted when its edges spread over at
#: least this many sampling errors (six independent normal edges span ~2.5 on
#: average and exceed 5 in under 1 % of cases).
STATIC_SPREAD_SE = 5
#: Template bonds outside the window, or contacts linking clusters per atom,
#: above this: the window, not the clusters.
CLUSTER_OUTLIER_FRACTION = 0.05
#: Clusters of other kinds (sizes) above this share of all: the window may cut
#: through or merge clusters.
MIXED_SIZE_SHARE = 0.1
#: A site is split (two positions) when its PC1 is resolved and its PC1 excess
#: kurtosis is at most -0.5 -- the value at which two equal wells +-delta of
#: thermal width sigma separate (delta = sigma: kappa = -2 delta^4 /
#: (delta^2 + sigma^2)^2 = -0.5) -- and at least SPLIT_SITE_SE sampling errors
#: (sqrt(24/n)) below zero. The one rule for the app: list_sites
#: (workers/clusters.js isSplitSite) and the split member sites here.
SPLIT_SITE_KURTOSIS = -0.5
SPLIT_SITE_SE = 3

CLUSTER_RULES = {
    "z": CLUSTER_Z,
    "effect": CLUSTER_EFFECT,
    "trigonalityEffect": TRIGONALITY_EFFECT,
    "axisEnhancement": AXIS_ENHANCEMENT,
    "staticSpreadSe": STATIC_SPREAD_SE,
    "outlierFraction": CLUSTER_OUTLIER_FRACTION,
    "mixedSizeShare": MIXED_SIZE_SHARE,
    "splitKurtosis": SPLIT_SITE_KURTOSIS,
    "splitSe": SPLIT_SITE_SE,
}

_NO_CLUSTERS = ("extended_network", "no_clusters")
#: Every ``because`` each finding can give (closed sets; the parity fixture pins
#: the strings, and the Copilot's renderer has a sentence entry for each).
CLUSTER_BECAUSE = {
    "distorted": ("edge_spread", "axis_family", "trigonality", "less_distorted", "within_null",
                  "undefined_statistic") + _NO_CLUSTERS,
    "lessDistorted": ("below_null", "above_null", "within_null", "single_edge", "undefined_statistic") + _NO_CLUSTERS,
    "breathing": ("above_null", "below_null", "within_null", "undefined_statistic") + _NO_CLUSTERS,
    "rigid": ("less_distorted", "stiffer_bonds", "split_sites", "pattern", "within_null",
              "undefined_statistic") + _NO_CLUSTERS,
    "staticDistortion": ("above_sampling_error", "within_sampling_error") + _NO_CLUSTERS,
    "splitMemberSites": ("half_or_more_split", "fewer_than_half_split", "none_split") + _NO_CLUSTERS,
    "trigonalLean": ("leans_trigonal", "leans_one_edge", "within_null", "not_tetrahedra",
                     "undefined_statistic") + _NO_CLUSTERS,
    "favouredAxes": ("enhanced", "within_null", "undefined_statistic") + _NO_CLUSTERS,
    "favouredAxes.families": ("enhanced", "within_null", "undefined_statistic"),
    "window.extendedNetwork": ("extended_network", "finite_clusters", "no_clusters"),
    "window.mixedSizes": ("above_bar", "within_bar") + _NO_CLUSTERS,
    "window.bondsOutsideWindow": ("above_bar", "within_bar") + _NO_CLUSTERS,
    "window.contactsBetweenClusters": ("above_bar", "within_bar") + _NO_CLUSTERS,
}
#: What an Evidence record is read against, and how its change is formed.
EVIDENCE_AGAINST = ("shapeNull", "independentNull", "samplingError", "siteKurtosis", "window")
CHANGE_KINDS = ("ratio", "difference")


def evidence_clears(record: dict) -> bool:
    """Whether an Evidence record clears its bar, from the record alone.

    Every field the bar names must be defined and clear its threshold the
    record's way: ``side`` 'above' value >= bar, 'below' value <= -bar, absent
    |value| >= bar; ``strict`` makes each comparison strict. A record without a
    bar clears nothing.
    """
    bar = record.get("bar") or {}
    if not bar:
        return False
    side = record.get("side")
    strict = bool(record.get("strict"))
    for field, threshold in bar.items():
        value = record.get(field)
        if value is None:
            return False
        if side == "above":
            passed = value > threshold if strict else value >= threshold
        elif side == "below":
            passed = value < -threshold if strict else value <= -threshold
        else:
            passed = abs(value) > threshold if strict else abs(value) >= threshold
        if not passed:
            return False
    return True


def _finite(value) -> float | None:
    """A finite float, or None (undefined, or non-finite)."""
    if value is None:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def _evidence(statistic: str, against: str, observed, null, z, change_kind, bar: dict, side=None, **extra) -> dict:
    """One Evidence record; its change from observed and null by its kind."""
    observed, null, z = _finite(observed), _finite(null), _finite(z)
    change = None
    if change_kind is not None and observed is not None and null is not None:
        if change_kind == "ratio":
            change = _finite(observed / null - 1.0) if null != 0 else None
        else:
            change = _finite(observed - null)
    record = {"statistic": statistic, "against": against, **extra, "observed": observed, "null": null, "z": z,
              "change": change, "changeKind": change_kind, "bar": dict(bar)}
    if side is not None:
        record["side"] = side
    record["clears"] = evidence_clears(record)
    return record


def split_site_evidence(row: dict) -> dict:
    """The split rule on one site row of :func:`rmc_toolkits.pca_kde.site_ellipsoids`.

    observed = the PC1 excess kurtosis (its change from 0, a Gaussian's), z =
    that kurtosis in sampling errors of a Gaussian's sample kurtosis,
    kappa / sqrt(24/n) -- None unless PC1 is resolved (an unresolved axis's
    direction, and so its kurtosis, is noise). It clears at z <= -SPLIT_SITE_SE
    and kappa <= SPLIT_SITE_KURTOSIS.
    """
    kurtosis = row.get("excessKurtosis") or [None]
    kappa = _finite(kurtosis[0])
    count = int(row.get("count") or 0)
    resolved = bool((row.get("axisResolved") or [False])[0])
    z = kappa / math.sqrt(24.0 / count) if kappa is not None and resolved and count > 0 else None
    return _evidence(
        "pc1ExcessKurtosis", "siteKurtosis", kappa, 0.0, z, "difference",
        {"z": SPLIT_SITE_SE, "change": -SPLIT_SITE_KURTOSIS}, "below", site=int(row["referenceNumber"]),
    )


def is_split_site(row: dict) -> bool:
    """Whether one site row is split (two positions): its :func:`split_site_evidence` clears."""
    return split_site_evidence(row)["clears"]


def _defined(record: dict | None) -> bool:
    return record is not None and record["z"] is not None and record["change"] is not None


def _finding(answer, because: str, evidence=(), **extra) -> dict:
    """A Finding; its evidence with the records that clear first, the rest after, each in the given order."""
    records = list(evidence)
    return {"answer": answer, "because": because, **extra,
            "evidence": [r for r in records if r["clears"]] + [r for r in records if not r["clears"]]}


def _no_cluster_verdict(because: str) -> dict:
    """An extended network or no clusters: every shape finding None, ``because`` says which."""
    def empty(**extra):
        return _finding(None, because, **extra)

    return {
        "rules": dict(CLUSTER_RULES),
        "distorted": empty(),
        "lessDistorted": empty(),
        "breathing": empty(),
        "rigid": empty(),
        "staticDistortion": empty(),
        "splitMemberSites": empty(sites=[], memberSites=[]),
        "trigonalLean": empty(),
        "favouredAxes": empty(families=[]),
        "window": {
            "extendedNetwork": _finding(because == "extended_network", because),
            "mixedSizes": empty(),
            "bondsOutsideWindow": empty(),
            "contactsBetweenClusters": empty(),
        },
    }


def cluster_verdict(result: dict, member_rows: list[dict]) -> dict:
    """The findings of a :func:`cluster_distortion` result that has clusters.

    ``member_rows`` are the analysed kind's member sites' rows of
    :func:`rmc_toolkits.pca_kde.site_ellipsoids`, in reference order. The
    shape statistics (ladder width, trigonality, axis families) are read
    against the shape null (each cluster's own breathing kept), the size
    statistic (breathing sd) against independent atoms (the 2026-10-14
    review's L2).
    """
    observed = result["observed"]
    nulls = {"shapeNull": result["shapeNull"], "independentNull": {"mean": result["nullMean"], "z": result["z"]}}
    kind = result["kind"]

    def null_test(path, against, effect, side, change_kind="ratio"):
        """A statistic (a key of ``observed``, or a key path) against one of the nulls."""
        values = [observed, nulls[against]["mean"], nulls[against]["z"]]
        for part in path:
            values = [value.get(part) if isinstance(value, dict) else None for value in values]
        return _evidence(".".join(path), against, *values, change_kind, {"z": CLUSTER_Z, "change": effect}, side)

    width_up = null_test(("ladderWidth",), "shapeNull", CLUSTER_EFFECT, "above")
    width_down = null_test(("ladderWidth",), "shapeNull", CLUSTER_EFFECT, "below")
    breath_up = null_test(("breathingSd",), "independentNull", CLUSTER_EFFECT, "above")
    breath_down = null_test(("breathingSd",), "independentNull", CLUSTER_EFFECT, "below")
    tetrahedra = _finite(observed.get("eRms")) is not None
    trig = null_test(("t2Trigonality",), "shapeNull", TRIGONALITY_EFFECT, None, "difference") if tetrahedra else None
    # obs / null >= AXIS_ENHANCEMENT exactly when obs / null - 1 >= AXIS_ENHANCEMENT - 1:
    # the subtraction is exact for ratios in [0.5, 2] (Sterbenz), and outside it both fail or both pass.
    family_tests = []
    for axis_key, axis in (("longAxis", "long"), ("shortAxis", "short")):
        if not isinstance(observed.get(axis_key), dict):
            continue
        for family in observed[axis_key]:
            record = null_test((axis_key, family), "shapeNull", AXIS_ENHANCEMENT - 1.0, "above")
            family_tests.append((axis, family, record))
    family_records = [record for _, _, record in family_tests]
    favoured = [record for record in family_records if record["clears"]]
    trigonal = trig is not None and trig["clears"]
    pattern = bool(favoured) or trigonal

    # Distorted: the edge spread above the shape null, or a pattern (an axis
    # family more clusters favour, or a trigonal lean) whatever the spread.
    consulted = [width_up] + ([trig] if trig is not None else []) + family_records
    if width_up["clears"]:
        distorted = (True, "edge_spread")
    elif favoured:
        distorted = (True, "axis_family")
    elif trigonal:
        distorted = (True, "trigonality")
    elif not any(_defined(record) for record in consulted):
        distorted = (None, "undefined_statistic")
    elif width_down["clears"]:
        distorted = (False, "less_distorted")
    else:
        distorted = (False, "within_null")

    if kind["edges"] == 1:
        less = (None, "single_edge")  # a dimer: its one edge has no spread
    elif not _defined(width_down):
        less = (None, "undefined_statistic")
    elif width_down["clears"]:
        less = (True, "below_null")
    elif width_up["clears"]:
        less = (False, "above_null")
    else:
        less = (False, "within_null")

    if not _defined(breath_up):
        breathing = (None, "undefined_statistic")
    elif breath_up["clears"]:
        breathing = (True, "above_null")
    elif breath_down["clears"]:
        breathing = (False, "below_null")
    else:
        breathing = (False, "within_null")

    # Split member sites: at least half of the analysed kind's member sites.
    split_records = [split_site_evidence(row) for row in member_rows]
    split_sites = [record["site"] for record in split_records if record["clears"]]
    member_sites = [int(row["referenceNumber"]) for row in member_rows]
    if split_sites and 2 * len(split_sites) >= len(member_sites):
        split = (True, "half_or_more_split")
    elif split_sites:
        split = (False, "fewer_than_half_split")
    else:
        split = (False, "none_split")

    # Rigid: less distorted (and no pattern), or a size stiffer than independent
    # atoms -- never with split sites, where the null mixes the sites' static
    # choices and so already holds the distortion's amplitude.
    is_split = split[0] is True
    if not _defined(width_down) and not _defined(breath_down):
        rigid = (None, "undefined_statistic")
    elif is_split and (width_down["clears"] or breath_down["clears"]):
        rigid = (False, "split_sites")
    elif not is_split and width_down["clears"] and not pattern:
        rigid = (True, "less_distorted")
    elif not is_split and breath_down["clears"]:
        rigid = (True, "stiffer_bonds")
    elif width_down["clears"] and pattern:
        rigid = (False, "pattern")
    else:
        rigid = (False, "within_null")

    # The average cluster: its edge spread (null 0, a regular cluster) in
    # sampling errors -- the largest spread / error over the templates.
    ideal = result["ideal"]
    static_record = _evidence(
        "idealSpread", "samplingError", ideal["spread"], 0.0, ideal["spreadOverSe"], None,
        {"z": STATIC_SPREAD_SE}, "above",
    )
    static = (True, "above_sampling_error") if static_record["clears"] else (False, "within_sampling_error")

    if trig is None:
        trigonal_lean = (None, "not_tetrahedra")
    elif not _defined(trig):
        trigonal_lean = (None, "undefined_statistic")
    elif trig["clears"]:
        # Below the null: three edges at one vertex differ alike (C3v); above: one edge against its opposite.
        trigonal_lean = (True, "leans_trigonal" if trig["change"] < 0 else "leans_one_edge")
    else:
        trigonal_lean = (False, "within_null")

    families = []
    for axis, family, record in family_tests:
        if record["clears"]:
            answer = (True, "enhanced")
        elif _defined(record):
            answer = (False, "within_null")
        else:
            answer = (None, "undefined_statistic")
        families.append({"axis": axis, "family": family, "answer": answer[0], "because": answer[1]})
    if favoured:
        axes = (True, "enhanced")
    elif any(_defined(record) for record in family_records):
        axes = (False, "within_null")
    else:
        axes = (None, "undefined_statistic")

    # The window: other cluster sizes, bonds outside it, contacts between clusters.
    kinds = result["kinds"]
    total = sum(entry["clusters"] for entry in kinds)
    share = sum(entry["clusters"] for entry in kinds[1:]) / total if total else None

    def window_finding(statistic, value, bar):
        record = _evidence(statistic, "window", value, None, None, None, {"observed": bar}, "above", strict=True)
        return _finding(record["clears"], "above_bar" if record["clears"] else "within_bar", [record])

    outliers = result["outliers"]
    return {
        "rules": dict(CLUSTER_RULES),
        "distorted": _finding(*distorted, consulted),
        "lessDistorted": _finding(*less, [width_down]),
        "breathing": _finding(*breathing, [breath_up]),
        "rigid": _finding(*rigid, [width_down, breath_down]),
        "staticDistortion": _finding(*static, [static_record]),
        "splitMemberSites": _finding(*split, split_records, sites=split_sites, memberSites=member_sites),
        "trigonalLean": _finding(*trigonal_lean, [trig] if trig is not None else []),
        "favouredAxes": _finding(*axes, family_records, families=families),
        "window": {
            "extendedNetwork": _finding(False, "finite_clusters"),
            "mixedSizes": window_finding("otherKindsShare", share, MIXED_SIZE_SHARE),
            "bondsOutsideWindow": window_finding(
                "brokenEdgeFraction", outliers["brokenEdgeFraction"], CLUSTER_OUTLIER_FRACTION),
            "contactsBetweenClusters": window_finding(
                "extraContactsPerAtom", outliers["extraContactsPerAtom"], CLUSTER_OUTLIER_FRACTION),
        },
    }
