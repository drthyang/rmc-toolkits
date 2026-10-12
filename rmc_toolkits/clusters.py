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

The browser port is ``workers/clusters.js``; keep the two in sync (the parity
fixture ``tests/generate_clusters_fixture.py`` pins them). Shuffles come from a
counter-based hash, identical in both engines.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil

import numpy as np

from rmc_toolkits.pca_kde import SiteDisplacements, load_site_displacements
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
    """Cluster templates, observed shape statistics and their shuffle null.

    ``window`` is the inclusive bond-length window (Angstrom) that links the
    element's atoms into clusters in the average structure. Raises
    ``ValueError`` for an unknown element, a malformed window, or a file
    without cell indices (a coordinates-only configuration).
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
        {"members": m, "edges": e, "clusters": n} for (m, e), n in sorted(kinds.items(), key=lambda kv: -kv[1])
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

    observed_positions = positions(None)
    observed = _statistics(observed_positions, edges, families)
    nulls = [_statistics(positions(k), edges, families) for k in range(shuffles)]
    null_mean, null_sd, z = _null_summary(observed, nulls)
    result.update({"observed": observed, "nullMean": null_mean, "nullSd": null_sd, "z": z})

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
    return result


def cluster_distortion_from_file(path: str, element: str, rmin: float, rmax: float, shuffles: int = DEFAULT_SHUFFLES) -> dict:
    """:func:`cluster_distortion` of an ``.rmc6f`` file, read afresh."""
    return cluster_distortion(load_site_displacements(path), element, (rmin, rmax), shuffles=shuffles)
