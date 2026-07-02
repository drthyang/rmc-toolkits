"""Reusable analysis helpers for RMCProfile post-processing."""

__version__ = "0.3.0"

from .kde import (
    MAX_KDE_FIT_POINTS,
    UnitCellPositions,
    kde_slice,
    load_unit_cell_positions,
    oriented_kde_slice,
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
    "PlotResult",
    "Rmc6fAtom",
    "RmcStructure",
    "UnitCellPositions",
    "close_plot",
    "detect_plot_kind",
    "frac_lines_from_rmc6f",
    "iter_rmc6f_atoms",
    "kde_slice",
    "load_unit_cell_positions",
    "make_plot",
    "oriented_kde_slice",
    "plot_to_png",
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
