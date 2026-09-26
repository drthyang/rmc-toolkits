# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""PCA + KDE engine for RMC displacement clouds (thermal-ellipsoid analysis).

PCA and kernel density estimation are standard statistical tools. The specific
analysis here -- applying them per crystallographic site to RMC displacement
clouds to get the anisotropic displacement ellipsoid, plus the 3D Gaussian-KDE
isosurface with 2D projections on the box walls -- follows Maksim Eremenko's
PCA_KDE utilities:
https://github.com/MaximEremenko/Utilities/tree/main/RMCProfileUtilities/PCA_KDE

We followed that approach and reimplemented it independently -- this is not a
port of his code. His ``KDE.js`` evaluates a full multivariate Gaussian KDE
(Cholesky whitening, one kernel sum per grid point); the engine here samples on a
grid aligned to the cloud's principal axes so the estimator factorizes and costs
two orders of magnitude less (see below).

RMCProfile tags every atom with the reference site it belongs to and the cell
index of the box copy it sits in. Subtracting a site's average position from
each of its copies leaves one displacement cloud per crystallographic site: the
covariance of that cloud is the anisotropic displacement tensor U (its eigen
decomposition is the thermal ellipsoid), and a Gaussian KDE of the cloud is the
smooth probability density that the ellipsoid only approximates.

The KDE is the expensive half. With SciPy's covariance-scaled bandwidth, an
``M = grid**3`` volume costs ``O(N * M)`` kernel evaluations -- 1000 copies of a
site on a 48^3 grid is 1.1e8 exponentials, and a run like ``data/5K_try1`` has
52 sites.

This module never pays that cost, and never approximates to avoid it. SciPy's
bandwidth matrix is ``H = factor**2 * C`` with ``C`` the data covariance, so in
the eigenbasis of ``C`` -- the PCA frame, which is the frame this analysis wants
anyway -- both ``C`` and ``H`` are diagonal and the 3D Gaussian factorizes into
three independent 1D Gaussians. Evaluated on a grid aligned with those axes the
density becomes a tensor product,

    density[i, j, k] = norm * sum_m  Ax[i, m] * Ay[j, m] * Az[k, m]

which needs ``N * (Gx + Gy + Gz)`` exponentials (2e5 rather than 1.1e8) and
contracts through BLAS matrix products instead of a Python-level loop. The
result is the estimator ``scipy.stats.gaussian_kde`` defines, equal to
floating-point round-off -- ``tests/test_pca_kde.py`` asserts exactly that -- at
roughly two orders of magnitude less work.

Conventions
-----------
Displacements are Cartesian Angstrom. Principal axes are returned as rows of a
right-handed orthonormal matrix, sign-fixed so the largest component of each
axis is positive, which keeps results reproducible across platforms. Volumes are
sampled in the PCA frame; ``axes`` and ``mean`` map them back to Cartesian
(``x_cart = mean + pca_coords @ axes``).
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.stats import chi2

from .parsers import Rmc6fParseReport, iter_rmc6f_atoms, read_cell_vectors

# Cap on the number of cloud points fed to the KDE (the ADP table never
# subsamples). A site in an n x n x n box has n^3 copies, so this binds for
# element-pooled clouds (e.g. 32000 Se atoms) AND for a single site once the box
# has >= 28 cells per edge. The subset is a seeded draw that differs from the
# JS port's (mulberry32 Fisher-Yates), so above the cap the two runtimes agree
# only statistically; the density estimate is stable well below it.
MAX_PCA_FIT_POINTS = 20000

# The KDE needs a positive-definite bandwidth. A site whose cloud is flat along
# one axis (a frozen or perfectly ordered direction) drives an eigenvalue to
# zero and makes H singular; floor the eigenvalues at this fraction of the
# largest one, which bounds the condition number without visibly moving a
# well-conditioned site.
EIGENVALUE_FLOOR_RATIO = 1e-8

# Below this ratio of smallest to largest eigenvalue the cloud is effectively
# planar or linear and the ellipsoid is reported as degenerate. The same ratio
# marks an individual principal axis as collapsed (no spread along it), where
# the excess kurtosis m4/m2^2 - 3 is 0/0 and is reported as undefined (None).
DEGENERATE_RATIO = 1e-6

# Absolute floor on the largest displacement variance (Angstrom^2): below
# (1e-4 A)^2 a site has no displacement at all. An RMCProfile *AVERAGE.rmc6f or an
# ideal starting configuration puts every copy of a site at the same offset, and
# the covariance is float round-off (~1e-28 A^2) whose ratios, axes and kurtosis
# are noise. Such a site is flagged ``zeroSpread`` with anisotropy, kurtosis and
# axes reported as None, and ``pca_kde_volume`` refuses it. Five orders of
# magnitude below any physical ADP (zero-point motion alone gives U >~ 1e-3 A^2)
# and far above the round-off, including the ~1e-5 A quantisation of a 7-decimal
# coordinate in a ~100 A box.
ZERO_SPREAD_VARIANCE = 1e-8

# A principal axis is "resolved" -- its direction, per-axis kurtosis and crystal
# orientation mean something -- only when its eigenvalue is separated from each
# neighbour's by more than this many standard errors of the gap. For a truly
# degenerate pair (a cubic site, the in-plane pair of a uniaxial one) the gap
# statistic exceeds 3 in < 1% of samples at n = 216..8000 (Gaussian, t8 and
# <111>-split clouds, by simulation); its median is ~1, so a 1-SE rule would
# call half of all degenerate pairs resolved.
AXIS_RESOLUTION_SIGMAS = 3.0


@dataclass(frozen=True)
class SiteDisplacements:
    """Per-site displacement clouds for one ``.rmc6f`` configuration.

    ``displacements`` holds every atom's Cartesian (Angstrom) offset from the
    average position of its own reference site, so each site's rows are already
    centered. ``site_index`` maps a row to its site in ``reference_numbers``.

    A reference number may be carried by more than one species (a solid
    solution, or RMCProfile swap moves). ``elements`` then holds the site's
    *majority* species (ties to the alphabetically first), ``element_counts`` its
    full composition, and ``atom_elements`` every row's own species -- element
    pooling selects atoms by the latter, never by the site label.

    ``parse_warning`` is the atom section's problem list
    (:meth:`Rmc6fParseReport.warning`), e.g. atom lines skipped for non-finite
    coordinates, or ``None`` for a clean file; the API and the pages show it,
    so a dropped atom is never silent.
    """

    reference_numbers: np.ndarray  # (S,) RMCProfile reference number per site
    elements: list[str]  # (S,) majority species per site
    counts: np.ndarray  # (S,) box copies per site
    displacements: np.ndarray  # (N, 3) Cartesian Angstrom, centered per site
    site_index: np.ndarray  # (N,) row -> site
    site_fractional: np.ndarray  # (S, 3) average site position, unit-cell fraction
    lattice_vectors: np.ndarray  # (3, 3) supercell vectors (Angstrom)
    supercell: np.ndarray  # (3,)
    unit_vectors: np.ndarray  # (3, 3) single-cell vectors (Angstrom)
    atom_elements: np.ndarray | None = None  # (N,) each row's own species
    element_counts: tuple[dict[str, int], ...] | None = None  # (S,) composition per site
    parse_warning: str | None = None  # Rmc6fParseReport.warning() of the atom section

    @property
    def species(self) -> list[str]:
        """Every species present, sorted -- including minority species of mixed sites."""
        if self.atom_elements is not None:
            return sorted({str(element) for element in self.atom_elements})
        return sorted(set(self.elements))

    def site_position(self, reference_number: int) -> int:
        """Index of ``reference_number`` in the per-site arrays."""
        matches = np.flatnonzero(self.reference_numbers == int(reference_number))
        if matches.size == 0:
            available = ", ".join(str(value) for value in self.reference_numbers)
            raise ValueError(
                f"Unknown reference number {reference_number}; available: {available}"
            )
        return int(matches[0])


def _circular_site_centres(offsets: np.ndarray, site_index: np.ndarray, site_count: int) -> np.ndarray:
    """Per-site, per-axis circular mean of period-1 offsets, shape (S, 3).

    ``atan2(sum sin 2 pi o, sum cos 2 pi o) / 2 pi`` -- the centre of a cloud on a
    circle, which no wrap of an individual offset can move. Only used as the
    reference that the unwrap is taken about; the site mean itself is the plain
    arithmetic mean of the unwrapped offsets.
    """
    angle = 2.0 * np.pi * offsets
    centres = np.empty((site_count, 3))
    for axis in range(3):
        sin_sum = np.bincount(site_index, weights=np.sin(angle[:, axis]), minlength=site_count)
        cos_sum = np.bincount(site_index, weights=np.cos(angle[:, axis]), minlength=site_count)
        centres[:, axis] = np.arctan2(sin_sum, cos_sum) / (2.0 * np.pi)
    return centres


def load_site_displacements(rmc6f_path: str | Path) -> SiteDisplacements:
    """Read an ``.rmc6f`` file into per-site Cartesian displacement clouds.

    An atom's stored coordinate is a fraction of the *supercell*, and its cell
    indices locate the box copy it belongs to, so ``coords - cell_indices /
    supercell`` is the offset from that copy's origin -- the site's own position
    in the cell plus its displacement, ``(s + delta) / N``, known only modulo one
    supercell period because the stored coordinate wraps at the box edge.

    The wrap is undone about each site's own centre, not about zero: per site and
    axis the circular mean ``c = arg(sum exp(2 pi i o)) / 2 pi`` of the offsets
    (period 1) is the reference, and every offset is moved by the whole number of
    periods that brings it within half a period of ``c``. A fold about zero (the
    earlier half-box minimum image) cuts a site in two whenever ``s / N`` reaches
    one half -- a site at x = 1/2 in a one-cell-thick box, or near x = 1 in a
    two-cell one -- and puts the halves a whole box edge apart. About the
    site's centre only a genuine displacement of half a supercell could fold.
    Subtracting each site's mean offset then leaves the displacement about the
    average structure.

    Atoms are read with the shared ``.rmc6f`` grammar (:func:`iter_rmc6f_atoms`,
    the browser's ``parseRmc6fAtoms``): a line with a non-finite coordinate is
    skipped -- one NaN would otherwise poison its site's mean -- and counted,
    and the report's warning naming the first such line is returned as
    ``parse_warning``. No parseable full-layout atom at all is a ``ValueError``
    that says what was found.
    """
    rmc6f_path = Path(rmc6f_path)
    lattice_vectors, supercell = read_cell_vectors(rmc6f_path)
    unit_vectors = lattice_vectors / supercell[:, None]

    coords: list[np.ndarray] = []
    cells: list[np.ndarray] = []
    references: list[int] = []
    elements: list[str] = []
    report = Rmc6fParseReport()
    for atom in iter_rmc6f_atoms(rmc6f_path, report=report):
        coords.append(atom["coords"])
        cells.append(atom["cell_indices"])
        references.append(atom["reference_number"])
        elements.append(atom["element"])

    if not coords:
        detail = report.warning() or (
            f"{report.coords_only_atoms} coordinate-only atom lines carry no reference "
            "site or cell indices"
            if report.coords_only_atoms
            else "the Atoms section is empty" if report.has_atoms_section
            else "there is no Atoms section"
        )
        raise ValueError(f"{rmc6f_path}: no atoms could be parsed — {detail}")

    coords_array = np.asarray(coords, dtype=float)
    cells_array = np.asarray(cells, dtype=float)
    references_array = np.asarray(references, dtype=int)

    reference_numbers, site_index = np.unique(references_array, return_inverse=True)
    site_count = reference_numbers.size
    counts = np.bincount(site_index, minlength=site_count)

    # Offset from the atom's own box copy, unwrapped over the supercell edge about
    # the site's circular mean (see the docstring): o -= round(o - centre).
    offsets = coords_array - cells_array / supercell
    offsets -= np.round(offsets - _circular_site_centres(offsets, site_index, site_count)[site_index])

    site_mean = np.column_stack(
        [
            np.bincount(site_index, weights=offsets[:, axis], minlength=site_count) / counts
            for axis in range(3)
        ]
    )
    centered = offsets - site_mean[site_index]
    displacements = centered @ lattice_vectors

    atom_elements = np.asarray(elements, dtype=object)
    element_counts, site_elements = _site_compositions(atom_elements, site_index, site_count)

    return SiteDisplacements(
        reference_numbers=reference_numbers,
        elements=site_elements,
        counts=counts,
        displacements=displacements,
        site_index=site_index,
        site_fractional=np.mod(site_mean * supercell, 1.0),
        lattice_vectors=lattice_vectors,
        supercell=supercell,
        unit_vectors=unit_vectors,
        atom_elements=atom_elements,
        element_counts=element_counts,
        parse_warning=report.warning(),
    )


def _site_compositions(
    atom_elements: np.ndarray, site_index: np.ndarray, site_count: int
) -> tuple[tuple[dict[str, int], ...], list[str]]:
    """Per-site species counts and the site label: the majority species, ties to
    the alphabetically first -- the same rule as ``workers/pcaKde.js``."""
    counts: list[dict[str, int]] = [{} for _ in range(site_count)]
    for element, index in zip(atom_elements, site_index):
        composition = counts[index]
        composition[element] = composition.get(element, 0) + 1
    ordered = tuple(dict(sorted(composition.items())) for composition in counts)
    labels = [
        min(composition, key=lambda name: (-composition[name], name)) if composition else ""
        for composition in ordered
    ]
    return ordered, labels


@lru_cache(maxsize=8)
def cached_site_displacements(path: str, mtime: float) -> SiteDisplacements:
    """:func:`load_site_displacements` memoized on (path, ``mtime``) for library callers.

    ``mtime`` is whatever the CALLER passes: the file is not re-stat'ed here,
    so a stale value (or a coarse mtime that a quick rewrite does not change)
    returns the previous configuration's clouds. For files that change while
    they are served (Live Data), call the uncached
    :func:`load_site_displacements` under a cache keyed on the full file
    signature, as ``app.py``'s ``_SITES_CACHE`` does.
    """
    return load_site_displacements(path)


def displacement_cloud(
    sites: SiteDisplacements,
    *,
    reference_number: int | None = None,
    element: str | None = None,
) -> np.ndarray:
    """Select one site's cloud, or pool every site of an element.

    Pooling is meaningful because each site is centered on its own average
    position: the union describes how an element moves, independent of where its
    sites sit in the cell.
    """
    if reference_number is not None:
        position = sites.site_position(reference_number)
        return sites.displacements[sites.site_index == position]

    if element not in (None, "", "all"):
        wanted = str(element).lower()
        if sites.atom_elements is not None:
            # Each atom by its OWN species: a mixed site contributes only the
            # matching atoms (centred on the site's all-species mean position).
            rows = np.asarray([str(name).lower() == wanted for name in sites.atom_elements], dtype=bool)
        else:
            positions = np.flatnonzero(np.asarray([name.lower() for name in sites.elements]) == wanted)
            rows = np.isin(sites.site_index, positions)
        if not rows.any():
            available = ", ".join(sites.species)
            raise ValueError(f"Unknown element {element!r}; available: {available}")
        return sites.displacements[rows]

    return sites.displacements


def probability_scale(probability: float) -> float:
    """Ellipsoid scale factor k such that k*sigma encloses ``probability``.

    The squared Mahalanobis radius of a 3D Gaussian is chi-square with three
    degrees of freedom, so k = sqrt(chi2.ppf(p, 3)) -- 1.5382 at the
    crystallographic 50% convention.
    """
    probability = float(probability)
    if not 0.0 < probability < 1.0:
        raise ValueError("probability must lie strictly between 0 and 1")
    return float(np.sqrt(chi2.ppf(probability, df=3)))


def _canonical_axes(vectors: np.ndarray) -> np.ndarray:
    """Fix eigenvector signs and handedness so results are reproducible.

    ``eigh`` fixes each axis only up to a sign, and the sign it picks depends on
    the LAPACK build. Force the largest-magnitude component of every axis
    positive, then flip the last axis if needed to keep the frame right-handed.
    Axes of a near-isotropic cloud stay arbitrary within their degenerate
    subspace -- no convention can fix that. Anything measured *along* such an
    axis (its kurtosis, its crystal orientation) is then frame noise, which is
    why ``_shape_statistics`` reports a rotation-invariant non-Gaussianity and
    flags unresolved axes.
    """
    vectors = np.array(vectors, dtype=float, copy=True)
    lead = np.argmax(np.abs(vectors), axis=-2)
    signs = np.sign(np.take_along_axis(vectors, lead[..., None, :], axis=-2)[..., 0, :])
    signs[signs == 0] = 1.0
    vectors *= signs[..., None, :]
    handedness = np.sign(np.linalg.det(vectors))
    vectors[..., 2] *= handedness[..., None]
    return vectors


def _eigen_decomposition(covariance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Descending eigenvalues and canonical axes (rows) of symmetric matrices."""
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    eigenvalues = eigenvalues[..., ::-1]
    eigenvectors = _canonical_axes(eigenvectors[..., ::-1])
    # Rows are easier to read and serialize than eigh's column layout.
    axes = np.swapaxes(eigenvectors, -1, -2)
    return np.maximum(eigenvalues, 0.0), axes


def site_ellipsoids(
    sites: SiteDisplacements,
    *,
    probability: float = 0.5,
) -> list[dict]:
    """Anisotropic displacement tensor and ellipsoid for every site, in one pass.

    All sites share the same per-atom arrays, so the six independent components
    of every covariance come from six ``bincount`` reductions and the eigen
    problem is solved batched -- the whole table costs milliseconds even for a
    52-site, 52000-atom configuration.
    """
    displacements = sites.displacements
    site_index = sites.site_index
    site_count = sites.reference_numbers.size
    counts = sites.counts.astype(float)
    finite_rows = np.isfinite(displacements).all(axis=1)
    if not finite_rows.all():
        bad = sorted({int(sites.reference_numbers[i]) for i in site_index[~finite_rows]})
        raise ValueError(f"displacement cloud contains non-finite coordinates (reference sites {bad})")

    covariance = np.zeros((site_count, 3, 3))
    for a in range(3):
        for b in range(a, 3):
            products = np.bincount(
                site_index,
                weights=displacements[:, a] * displacements[:, b],
                minlength=site_count,
            )
            covariance[:, a, b] = products
            covariance[:, b, a] = products
    # Clouds are already centered per site, so the sum of outer products over
    # n - 1 is the unbiased covariance -- the same estimator gaussian_kde uses.
    covariance /= np.maximum(counts - 1.0, 1.0)[:, None, None]

    eigenvalues, axes = _eigen_decomposition(covariance)
    scale = probability_scale(probability)
    u_eq = eigenvalues.mean(axis=1)
    zero_spread = eigenvalues[:, 0] < ZERO_SPREAD_VARIANCE
    largest = np.maximum(eigenvalues[:, 0], 1e-30)
    ratio = eigenvalues[:, 2] / largest

    # Shape statistics in each site's own PCA frame. Projecting per atom onto its
    # site axes costs one einsum.
    projected = np.einsum("nij,nj->ni", axes[site_index], displacements)
    excess_kurtosis, non_gaussianity, resolved = _shape_statistics(
        projected, site_index, site_count, eigenvalues
    )

    composition = sites.element_counts
    ellipsoids: list[dict] = []
    for index in range(site_count):
        zero = bool(zero_spread[index])
        ellipsoids.append(
            {
                "referenceNumber": int(sites.reference_numbers[index]),
                "element": sites.elements[index],
                "mixed": composition is not None and len(composition[index]) > 1,
                "elementCounts": dict(composition[index]) if composition is not None
                else {sites.elements[index]: int(sites.counts[index])},
                "count": int(sites.counts[index]),
                "siteFractional": sites.site_fractional[index].tolist(),
                "covariance": covariance[index].tolist(),
                "eigenvalues": eigenvalues[index].tolist(),
                # A zero-spread site's eigenvectors are round-off: no axes.
                "axes": None if zero else axes[index].tolist(),
                "rms": np.sqrt(eigenvalues[index]).tolist(),
                "semiAxes": (scale * np.sqrt(eigenvalues[index])).tolist(),
                "probability": float(probability),
                "uIso": float(u_eq[index]),
                "bIso": float(8.0 * np.pi**2 * u_eq[index]),
                "rmsIso": float(np.sqrt(max(u_eq[index], 0.0))),
                "excessKurtosis": _nullable(excess_kurtosis[index]),
                "axisResolved": resolved[index].tolist(),
                "nonGaussianity": _nullable(non_gaussianity[index:index + 1])[0],
                # sqrt of the eigenvalue ratio: the ellipsoid's long/short axis.
                "anisotropy": None if zero
                else float(np.sqrt(largest[index] / max(eigenvalues[index, 2], 1e-30))),
                "degenerate": bool(zero or ratio[index] < DEGENERATE_RATIO),
                "zeroSpread": zero,
            }
        )
    return ellipsoids


def _shape_statistics(
    projected: np.ndarray,
    site_index: np.ndarray,
    site_count: int,
    eigenvalues: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Kurtosis of every site's cloud: per axis, multivariate, and axis resolution.

    ``projected`` holds each (centred) atom in its own site's PCA frame and
    ``eigenvalues`` the sites' descending covariance eigenvalues, shape (S, 3).
    Moments are population (1/n) moments. Returns, per site:

    * ``kappa`` (S, 3) -- excess kurtosis ``m4 / m2**2 - 3`` along each principal
      axis; NaN where undefined: no spread (``lambda_1 < ZERO_SPREAD_VARIANCE``)
      or a collapsed axis (``lambda_a < DEGENERATE_RATIO * lambda_1``), where m2 is
      round-off and the ratio 0/0. One rule, no floor constant.
    * ``non_gaussianity`` (S,) -- Mardia's multivariate kurtosis
      ``b2 = mean[(u^T S^-1 u)^2]`` over the d defined axes, normalised as
      ``3 (b2 - d(d+2)) / (d(d+2))`` -- ``(b2 - 15)/5`` for d = 3. It is affine
      invariant, so it cannot depend on how sampling noise orients the PCA
      frame of a (near-)isotropic site, and for any elliptical distribution it
      equals the excess kurtosis along every direction. NaN when d = 0.
    * ``resolved`` (S, 3) bool -- the axis is separated from each neighbour by
      more than ``AXIS_RESOLUTION_SIGMAS`` standard errors of the eigenvalue gap,
      ``SE(lambda_a) = sqrt((m4_a - m2_a**2) / n)``; a gap between two collapsed
      (undefined) axes never counts. An unresolved axis lies in a
      (near-)degenerate subspace: its direction, and so its kappa and crystal
      orientation, is set by sampling noise and biased toward the outliers.

    ``workers/pcaKde.js`` (``shapeStatistics``) is the line-for-line port.
    """
    counts = np.maximum(np.bincount(site_index, minlength=site_count).astype(float), 1.0)
    m2 = np.column_stack(
        [np.bincount(site_index, weights=projected[:, a] ** 2, minlength=site_count) for a in range(3)]
    ) / counts[:, None]
    m4 = np.column_stack(
        [np.bincount(site_index, weights=projected[:, a] ** 4, minlength=site_count) for a in range(3)]
    ) / counts[:, None]

    has_spread = eigenvalues[:, :1] >= ZERO_SPREAD_VARIANCE
    defined = has_spread & (eigenvalues >= DEGENERATE_RATIO * eigenvalues[:, :1])
    with np.errstate(divide="ignore", invalid="ignore"):
        kappa = np.where(defined, m4 / (m2 * m2) - 3.0, np.nan)
        inverse = np.where(defined, 1.0 / m2, 0.0)

    radius2 = np.einsum("ni,ni->n", projected * projected, inverse[site_index])
    b2 = np.bincount(site_index, weights=radius2 * radius2, minlength=site_count) / counts
    d = defined.sum(axis=1)
    reference = d * (d + 2.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        non_gaussianity = np.where(d > 0, 3.0 * (b2 - reference) / reference, np.nan)

    error = np.sqrt(np.maximum(m4 - m2 * m2, 0.0) / counts[:, None])
    # A gap between two collapsed axes compares two round-off eigenvalues: it
    # is never a resolution (a linear cloud's PC2/PC3 directions are noise).
    # One collapsed axis beside a defined one is still resolved (plane normal).
    gap = [
        ((m2[:, a] - m2[:, a + 1]) > AXIS_RESOLUTION_SIGMAS * np.hypot(error[:, a], error[:, a + 1]))
        & (defined[:, a] | defined[:, a + 1])
        for a in (0, 1)
    ]
    resolved = np.column_stack([gap[0], gap[0] & gap[1], gap[1]]) & has_spread
    return kappa, non_gaussianity, resolved


def _nullable(values: np.ndarray) -> list:
    """List of floats with NaN (undefined) mapped to None -- JSON null, never a NaN token."""
    return [None if not np.isfinite(value) else float(value) for value in np.asarray(values)]


def _bandwidth_factor(method: str | float, count: int, dimensions: int) -> float:
    """SciPy's ``gaussian_kde`` covariance factor for ``count`` unweighted points."""
    if isinstance(method, (int, float)) and not isinstance(method, bool):
        factor = float(method)
        if not (np.isfinite(factor) and factor > 0):
            raise ValueError("numeric bandwidth must be a positive finite number")
        return factor

    name = str(method).lower()
    if name == "scott":
        return float(count ** (-1.0 / (dimensions + 4)))
    if name == "silverman":
        return float((count * (dimensions + 2) / 4.0) ** (-1.0 / (dimensions + 4)))
    raise ValueError("bw must be 'scott', 'silverman', or a positive number")


def _subsample(points: np.ndarray, limit: int, rng_seed: int) -> np.ndarray:
    if points.shape[0] <= limit:
        return points
    rng = np.random.default_rng(rng_seed)
    choice = rng.choice(points.shape[0], limit, replace=False)
    choice.sort()
    return points[choice]


def _kernel_matrix(grid_axis: np.ndarray, projected: np.ndarray, bandwidth: float) -> np.ndarray:
    """1D Gaussian kernel values, shape (grid, points), for one principal axis."""
    scaled = (grid_axis[:, None] - projected[None, :]) / bandwidth
    return np.exp(-0.5 * scaled * scaled)


def _iso_levels(
    density: np.ndarray,
    cell_volume: float,
    probabilities: np.ndarray,
) -> tuple[list[dict], list[dict], float]:
    """Isosurface thresholds by enclosed probability mass and by raw density.

    The mass levels are the physically meaningful ones: the surface at level
    ``tau(p)`` bounds the smallest region holding fraction ``p`` of the cloud, so
    it is directly comparable with the p% ellipsoid from the covariance.
    """
    flat = np.sort(density.ravel())[::-1]
    cumulative = np.cumsum(flat) * cell_volume
    mass = float(cumulative[-1]) if cumulative.size else 0.0

    mass_levels: list[dict] = []
    if mass > 0:
        targets = probabilities * mass
        indices = np.searchsorted(cumulative, targets, side="left")
        indices = np.clip(indices, 0, flat.size - 1)
        mass_levels = [
            {"p": float(p), "level": float(flat[index])}
            for p, index in zip(probabilities, indices)
        ]

    vmin = float(density.min()) if density.size else 0.0
    vmax = float(density.max()) if density.size else 0.0
    density_levels = [
        {"p": float(p), "level": float(vmin + p * (vmax - vmin))} for p in probabilities
    ]
    return mass_levels, density_levels, mass


def _projection(
    kernels: list[np.ndarray],
    bandwidths: np.ndarray,
    axis_coords: list[np.ndarray],
    first: int,
    second: int,
    count: int,
) -> dict:
    """2D KDE of the cloud projected onto a principal plane.

    The 3D covariance is diagonal in the PCA frame, so every 2x2 marginal is
    diagonal too: the plane's own Scott bandwidth stays axis-aligned and the
    whole projection is a single matrix product of the two 1D kernel matrices.
    """
    density = (kernels[first] @ kernels[second].T) / (
        count * 2.0 * np.pi * bandwidths[first] * bandwidths[second]
    )
    return {
        "density": density.tolist(),
        "extent": [
            float(axis_coords[first][0]),
            float(axis_coords[first][-1]),
            float(axis_coords[second][0]),
            float(axis_coords[second][-1]),
        ],
        "axes": [int(first), int(second)],
        "bandwidth": [float(bandwidths[first]), float(bandwidths[second])],
        "vmax": float(density.max()) if density.size else 0.0,
    }


def pca_kde_volume(
    points: np.ndarray,
    *,
    bw: str | float = "scott",
    bw_scale: float = 1.0,
    grid: int = 48,
    extent: float = 3.0,
    cubic_box: bool = False,
    probability: float = 0.5,
    probabilities: np.ndarray | None = None,
    projections: bool = True,
    max_fit_points: int = MAX_PCA_FIT_POINTS,
    rng_seed: int = 0,
) -> dict:
    """PCA statistics and a 3D Gaussian KDE volume for one displacement cloud.

    The volume is sampled on a grid aligned with the cloud's principal axes,
    which is what makes the estimator separable (see the module docstring): the
    returned ``density`` is exactly what ``scipy.stats.gaussian_kde`` would give
    at the same points, computed with ``N * 3 * grid`` exponentials instead of
    ``N * grid**3``.

    ``extent`` is the half-width of the box in units of the kernel-broadened
    standard deviation along each axis, so the default 3.0 captures >99% of the
    density; the returned ``mass`` reports how much it actually captured.

    The volume, ``mass``, ``massLevels`` and the wall ``projections`` are always
    sampled on that per-axis box (``halfWidths``), whose G nodes resolve every
    axis's kernel. ``cubic_box`` only sizes the *display* box returned as
    ``boxHalfWidths`` -- a cube of half-width ``max(halfWidths)`` for a shadow box
    whose walls are all the same size. (Until 1.0 it also sampled every axis on
    that cube, which puts a thin axis's kernels between nodes: a planar cloud
    read 0% captured mass on an even grid and ~2e5 % on an odd one.)
    """
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must be a numeric array with shape (N, 3)")
    if points.shape[0] < 4:
        raise ValueError("a 3D KDE needs at least four points")
    if not np.isfinite(points).all():
        raise ValueError("displacement cloud contains non-finite coordinates")

    # NaN passes every `<= 0` test, so each parameter is checked for finiteness
    # too -- otherwise a NaN extent returns an all-NaN volume, not an error.
    if not np.isfinite(float(grid)):
        raise ValueError("grid must be a finite number")
    grid = int(max(8, min(int(grid), 128)))
    bw_scale = float(bw_scale)
    if not (np.isfinite(bw_scale) and bw_scale > 0):
        raise ValueError("bw_scale must be a positive finite number")
    extent = float(extent)
    if not (np.isfinite(extent) and extent > 0):
        raise ValueError("extent must be a positive finite number")

    total = int(points.shape[0])
    fit = _subsample(points, int(max_fit_points), rng_seed)
    count = int(fit.shape[0])

    mean = fit.mean(axis=0)
    centered = fit - mean
    covariance = np.cov(centered, rowvar=False, bias=False)
    eigenvalues, axes = _eigen_decomposition(covariance)

    # No spread at all (an average/ideal configuration): nothing to estimate.
    largest = float(eigenvalues[0])
    if not largest >= ZERO_SPREAD_VARIANCE:
        raise ValueError("displacement cloud has zero spread (RMS below 1e-4 A on every axis)")
    # A flat direction would make the bandwidth singular; floor it, and say so.
    ratio = float(eigenvalues[2] / largest)
    raw_eigenvalues = eigenvalues
    eigenvalues = np.maximum(eigenvalues, largest * EIGENVALUE_FLOOR_RATIO)

    factor = _bandwidth_factor(bw, count, 3) * bw_scale
    sigma = np.sqrt(eigenvalues)
    bandwidths = factor * sigma

    # Kernel-broadened spread: the KDE convolves the cloud with the kernel, so
    # its variance along an axis is sigma^2 + h^2 = sigma^2 (1 + factor^2).
    half_widths = extent * sigma * np.sqrt(1.0 + factor * factor)
    box_half_widths = np.full(3, float(half_widths.max())) if cubic_box else half_widths

    projected = centered @ axes.T  # cloud in the PCA frame, shape (N, 3)
    axis_coords = [np.linspace(-half, half, grid) for half in half_widths]
    kernels = [
        _kernel_matrix(axis_coords[axis], projected[:, axis], bandwidths[axis])
        for axis in range(3)
    ]

    # density[i, j, k] = norm * sum_m kx[i, m] ky[j, m] kz[k, m]. Contract one
    # PC3 slice at a time: each slice is a (grid x N) @ (N x grid) BLAS product,
    # so the whole volume never materializes an N x grid^3 temporary.
    norm = 1.0 / (
        count * (2.0 * np.pi) ** 1.5 * bandwidths[0] * bandwidths[1] * bandwidths[2]
    )
    density = np.empty((grid, grid, grid))
    for k in range(grid):
        density[:, :, k] = (kernels[0] * kernels[2][k]) @ kernels[1].T
    density *= norm

    cell_volume = float(np.prod([coords[1] - coords[0] for coords in axis_coords]))
    if probabilities is None:
        probabilities = np.linspace(0.0, 1.0, 101)
    probabilities = np.asarray(probabilities, dtype=float)
    mass_levels, density_levels, mass = _iso_levels(density, cell_volume, probabilities)

    # Kurtosis of the fit cloud -- the same statistics as ``site_ellipsoids``.
    kappa, non_gaussianity, resolved = _shape_statistics(
        projected, np.zeros(count, dtype=int), 1, raw_eigenvalues[None, :]
    )

    scale = probability_scale(probability)
    result = {
        "count": total,
        "fitCount": count,
        "mean": mean.tolist(),
        "covariance": covariance.tolist(),
        "eigenvalues": eigenvalues.tolist(),
        "axes": axes.tolist(),
        "rms": sigma.tolist(),
        "semiAxes": (scale * sigma).tolist(),
        "probability": float(probability),
        "uIso": float(eigenvalues.mean()),
        "bIso": float(8.0 * np.pi**2 * eigenvalues.mean()),
        "anisotropy": float(sigma[0] / sigma[2]),
        "excessKurtosis": _nullable(kappa[0]),
        "axisResolved": resolved[0].tolist(),
        "nonGaussianity": _nullable(non_gaussianity)[0],
        "degenerate": bool(ratio < DEGENERATE_RATIO),
        "zeroSpread": False,
        "bw": bw if isinstance(bw, str) else float(bw),
        "bwScale": bw_scale,
        "factor": float(factor),
        "bandwidth": bandwidths.tolist(),
        "grid": grid,
        "extent": float(extent),
        "cubicBox": bool(cubic_box),
        "halfWidths": half_widths.tolist(),
        "boxHalfWidths": box_half_widths.tolist(),
        "axisCoords": [coords.tolist() for coords in axis_coords],
        "cellVolume": cell_volume,
        # C-order over (PC1, PC2, PC3): index = (i * grid + j) * grid + k.
        "density": density.ravel().tolist(),
        "vmin": float(density.min()),
        "vmax": float(density.max()),
        # Fraction of the KDE's unit mass the box actually holds: the estimator
        # integrates to 1 over all space, so this is the truncation error.
        "mass": mass,
        "massLevels": mass_levels,
        "densityLevels": density_levels,
    }

    if projections:
        result["projections"] = {
            "pc12": _projection(kernels, bandwidths, axis_coords, 0, 1, count),
            "pc13": _projection(kernels, bandwidths, axis_coords, 0, 2, count),
            "pc23": _projection(kernels, bandwidths, axis_coords, 1, 2, count),
        }
    return result


def site_pca_kde(
    sites: SiteDisplacements,
    *,
    reference_number: int | None = None,
    element: str | None = None,
    **kwargs,
) -> dict:
    """``pca_kde_volume`` for one site (or one element's pooled sites)."""
    cloud = displacement_cloud(sites, reference_number=reference_number, element=element)
    result = pca_kde_volume(cloud, **kwargs)
    if reference_number is not None:
        position = sites.site_position(reference_number)
        result["referenceNumber"] = int(reference_number)
        result["element"] = sites.elements[position]
        if sites.element_counts is not None:
            result["elementCounts"] = dict(sites.element_counts[position])
            result["mixed"] = len(sites.element_counts[position]) > 1
        result["siteFractional"] = sites.site_fractional[position].tolist()
    elif element not in (None, "", "all"):
        result["element"] = str(element)
    return result
