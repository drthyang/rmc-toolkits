# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Reusable analysis helpers for RMCProfile post-processing."""

__version__ = "0.3.0"

from .kde import (
    MAX_KDE_FIT_POINTS,
    UnitCellPositions,
    kde_slice,
    load_unit_cell_positions,
    oriented_kde_slice,
)
from .pca_kde import (
    MAX_PCA_FIT_POINTS,
    SiteDisplacements,
    cached_site_displacements,
    displacement_cloud,
    load_site_displacements,
    pca_kde_volume,
    probability_scale,
    site_ellipsoids,
    site_pca_kde,
)
from .parsers import (
    CsvSeries,
    Rmc6fAtom,
    RmcStructure,
    frac_lines_from_rmc6f,
    iter_rmc6f_atoms,
    read_atom_indices,
    read_cell_vectors,
    read_chi,
    read_exafs_csv,
    read_rmc_csv,
    read_stog,
    read_structure,
    rwp,
    write_frac_from_rmc6f,
)
from .plots import PlotResult, close_plot, detect_plot_kind, make_plot, plot_to_png

__all__ = [
    "__version__",
    "CsvSeries",
    "MAX_KDE_FIT_POINTS",
    "MAX_PCA_FIT_POINTS",
    "PlotResult",
    "Rmc6fAtom",
    "RmcStructure",
    "SiteDisplacements",
    "UnitCellPositions",
    "cached_site_displacements",
    "close_plot",
    "detect_plot_kind",
    "displacement_cloud",
    "frac_lines_from_rmc6f",
    "iter_rmc6f_atoms",
    "kde_slice",
    "load_site_displacements",
    "load_unit_cell_positions",
    "make_plot",
    "oriented_kde_slice",
    "pca_kde_volume",
    "plot_to_png",
    "probability_scale",
    "site_ellipsoids",
    "site_pca_kde",
    "read_atom_indices",
    "read_cell_vectors",
    "read_chi",
    "read_exafs_csv",
    "read_rmc_csv",
    "read_stog",
    "read_structure",
    "rwp",
    "write_frac_from_rmc6f",
]
