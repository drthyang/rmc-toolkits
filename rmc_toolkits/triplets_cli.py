# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""``rmc-triplets`` -- bond-angle distribution CLI for RMC configurations.

Offline front end for :mod:`rmc_toolkits.triplets`: point it at an ``.rmc6f``
configuration (or a run folder containing one), name the A-B-C triplet with B
central, bound the two bond lengths, and it writes the angle histogram as a
commented CSV -- optionally with a PNG plot and the raw angle list.

Example
-------
::

    rmc-triplets data/5K_try1 --triplet Se Nb Se --bond12 2.2 2.9 \\
        --output se_nb_se.csv --plot se_nb_se.png
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import uuid
from pathlib import Path

from . import __version__

# The run-folder rule is shared with the web app (both runtimes): one
# configuration per folder everywhere.
from .parsers import find_run_configuration
from .triplets import BondAngleDistribution, bond_angles_from_rmc6f


def resolve_config(target: str | Path) -> Path:
    """Accept an ``.rmc6f`` file or a run folder (see ``find_run_configuration``)."""
    target = Path(target)
    if target.is_dir():
        return find_run_configuration(target)
    if not target.exists():
        raise FileNotFoundError(f"{target} does not exist")
    return target


def default_output_name(config: Path, triplet: tuple[str, str, str]) -> str:
    label = "-".join(triplet)
    return f"triplets_{label}_{config.stem}.csv"


def bond_count_text(unique: int, directed: int, end: str, apex: str) -> str:
    """Physical bond count, with the B-centred count when the two differ.

    An A-B bond with A = B is found from both of its ends, so the count of
    bond vectors seen from the central atoms is twice the number of bonds.
    """
    if end != apex:
        return f"{unique} physical bonds"
    return (
        f"{unique} physical bonds ({directed} bond vectors counted from the central "
        f"atoms: {end} is the central element, so each bond is seen from both ends)"
    )


def rmcprofile_sinth_factor(width_deg: float) -> float:
    """sin_corrected -> RMCProfile TRIPLETS ``norm/sin(theta)`` for ``width_deg`` bins.

    RMCProfile's column is the per-degree density over sin(bin centre); the
    engine's sin_corrected is the count fraction over sin(centre) sin(w/2).
    Their ratio is the constant ``sin(w/2) / w`` (``w/2`` in radians, ``w`` in
    degrees), ~pi/360 for small bins.
    """
    return math.sin(math.radians(width_deg) / 2.0) / width_deg


def write_csv(path: Path, config: Path, result: BondAngleDistribution) -> None:
    end1, apex, end2 = result.triplet
    lines = [
        f"# rmc-triplets bond-angle distribution",
        f"# configuration: {config}",
        f"# triplet (B central): {end1}-{apex}-{end2}",
        f"# bond12 window (Ang): {result.bond12[0]:g} .. {result.bond12[1]:g}",
        f"# bond23 window (Ang): {result.bond23[0]:g} .. {result.bond23[1]:g}",
        f"# central atoms: {result.apex_count}",
        "# bonds in window12: "
        + bond_count_text(result.unique_bonds12, result.bond12_count, end1, apex)
        + (f" (mean {result.mean_length12:.4f} Ang)" if result.mean_length12 else ""),
        "# bonds in window23: "
        + bond_count_text(result.unique_bonds23, result.bond23_count, end2, apex)
        + (f" (mean {result.mean_length23:.4f} Ang)" if result.mean_length23 else ""),
        f"# angles: {result.angle_count}",
        "# density is per degree with unit integral over [0, 180];",
        "# sin_corrected divides by the exact isotropic bin fraction (flat 1 = random).",
        # Same shape as RMCProfile's TRIPLETS norm/sin(theta), another scale:
        # its column is density / sin(centre) = sin_corrected * sin(w/2) / w.
        "# RMCProfile TRIPLETS norm/sin(theta) = sin_corrected * "
        f"{rmcprofile_sinth_factor(float(result.bin_edges[1] - result.bin_edges[0])):.10g}"
        " for this bin width (sin(w/2) / w, w/2 in rad, w in deg; ~pi/360).",
        "angle_deg,counts,density_per_deg,sin_corrected",
    ]
    for center, count, density, corrected in zip(
        result.bin_centers, result.counts, result.density, result.sin_corrected
    ):
        lines.append(f"{center:.4f},{int(count)},{density:.8e},{corrected:.8e}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def supported_plot_formats() -> dict[str, str]:
    """Matplotlib's savefig formats on this install: {extension: description}."""
    from matplotlib.figure import Figure

    return Figure().canvas.get_supported_filetypes()


def write_plot(path: Path, result: BondAngleDistribution, fmt: str | None = None) -> None:
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    figure = Figure(figsize=(7.0, 4.5), dpi=150)
    axes = figure.add_subplot(111)
    label = "-".join(result.triplet)
    axes.plot(
        result.bin_centers, result.sin_corrected, color="#1b6ca8", label="sin-corrected"
    )
    scale = result.sin_corrected.max() / result.density.max() if result.density.max() else 1.0
    axes.plot(
        result.bin_centers,
        result.density * scale,
        color="#c05640",
        linestyle="--",
        linewidth=1.0,
        # One axis: the density is drawn for its shape only, scaled so its
        # peak meets the sin-corrected peak (values are in the CSV).
        label="density, rescaled to the sin-corrected peak",
    )
    axes.set_xlim(0, 180)
    axes.set_xlabel("angle (deg)")
    axes.set_ylabel("sin-corrected distribution")
    axes.set_title(
        f"{label}  |  {result.bond12[0]:g}-{result.bond12[1]:g} A"
        + (
            f" / {result.bond23[0]:g}-{result.bond23[1]:g} A"
            if result.bond23 != result.bond12
            else ""
        )
        + f"  |  {result.angle_count} angles"
    )
    axes.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(path, format=fmt)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rmc-triplets",
        description=(
            "Bond-angle (triplet) distribution of an RMC configuration: pick an "
            "A-B-C triplet with B central, bound the two bond lengths, histogram "
            "the angle at B."
        ),
    )
    parser.add_argument("config", help=".rmc6f file, or a run folder containing one")
    parser.add_argument(
        "--triplet",
        nargs=3,
        metavar=("A", "B", "C"),
        required=True,
        help="the three atom types; the middle one (B) is the central atom",
    )
    parser.add_argument(
        "--bond12",
        nargs=2,
        type=float,
        metavar=("RMIN", "RMAX"),
        required=True,
        help="inclusive A-B bond-length window in Angstrom",
    )
    parser.add_argument(
        "--bond23",
        nargs=2,
        type=float,
        metavar=("RMIN", "RMAX"),
        default=None,
        help="inclusive B-C window in Angstrom (defaults to the A-B window)",
    )
    parser.add_argument(
        "--bin-width", type=float, default=1.0, help="histogram bin width in degrees"
    )
    parser.add_argument(
        "--output",
        default=None,
        help="CSV destination (default: triplets_<A-B-C>_<config>.csv in the "
        "current directory)",
    )
    parser.add_argument("--plot", default=None, help="also save a PNG plot here")
    parser.add_argument(
        "--dump-angles",
        default=None,
        metavar="PATH",
        help="also write every raw angle (degrees, one per line)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite existing output files instead of refusing (never the "
        "configuration, a directory, or two outputs onto one file)",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


class DestinationError(ValueError):
    """A destination that cannot be written: reported on one line, exit 1."""


def _same_destination(first: Path, second: Path) -> bool:
    # Case-folded: one file on the default macOS/Windows filesystems.
    if os.path.normcase(str(first.resolve())).casefold() == os.path.normcase(str(second.resolve())).casefold():
        return True
    try:
        return first.exists() and second.exists() and os.path.samefile(first, second)
    except OSError:
        return False


def check_destinations(destinations: dict[str, Path], config: Path, force: bool) -> None:
    """Refuse, before computing, any destination set that cannot be written whole.

    The histogram, plot and angle list must be different files, none of them
    the configuration or a directory, each with a folder that exists or can be
    created, the plot's extension a format matplotlib can write here, and --
    unless ``force`` -- none may exist yet.
    """
    items = list(destinations.items())
    for index, (flag, path) in enumerate(items):
        for other_flag, other in items[index + 1:]:
            if _same_destination(path, other):
                raise DestinationError(f"{flag} and {other_flag} must be different files (both {path})")
        if _same_destination(path, config):
            raise DestinationError(f"{flag} {path} is the input configuration")
        if path.is_dir():
            raise DestinationError(f"{flag} {path} is a directory")
        ancestor = path.parent
        while not ancestor.exists() and ancestor != ancestor.parent:
            ancestor = ancestor.parent
        if not ancestor.is_dir():
            raise DestinationError(f"{flag} {path}: its folder {ancestor} is not a directory")
    plot = destinations.get("--plot")
    if plot is not None and plot.suffix:
        formats = supported_plot_formats()
        if plot.suffix[1:].lower() not in formats:
            raise DestinationError(
                f"unsupported plot format '{plot.suffix}' for --plot {plot}; "
                f"use one of: {', '.join(sorted(formats))}"
            )
    if not force:
        existing = [str(path) for path in destinations.values() if path.exists()]
        if existing:
            raise DestinationError(
                "refusing to overwrite " + ", ".join(existing) + "; pass --force to replace"
            )


def _temporary_sibling(path: Path) -> Path:
    # Same folder (atomic rename) and same extension (matplotlib reads it).
    return path.with_name(f".{path.stem}.{os.getpid()}.{uuid.uuid4().hex[:12]}.tmp{path.suffix}")


def write_all(writers: list[tuple[Path, object]]) -> None:
    """Write every file through a temporary sibling, then rename them all.

    Nothing is renamed into place until every write has succeeded, so a
    failure leaves no partial set of outputs (and any previous files intact).
    """
    temporaries: list[tuple[Path, Path]] = []
    try:
        for path, write in writers:
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = _temporary_sibling(path)
            temporaries.append((temporary, path))
            write(temporary)
        for temporary, path in temporaries:
            os.replace(temporary, path)
    except BaseException:
        for temporary, _ in temporaries:
            try:
                temporary.unlink()
            except OSError:
                pass
        raise


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        config = resolve_config(args.config)
        # Every destination is checked before any work: an unsupported plot
        # format or a clash used to surface only after the CSV was written.
        triplet = tuple(str(symbol).strip().capitalize() for symbol in args.triplet)
        output = Path(args.output) if args.output else Path.cwd() / default_output_name(config, triplet)
        destinations = {"--output": output}
        if args.dump_angles:
            destinations["--dump-angles"] = Path(args.dump_angles)
        if args.plot:
            destinations["--plot"] = Path(args.plot)
        check_destinations(destinations, config, args.force)
        result = bond_angles_from_rmc6f(
            config,
            triplet=tuple(args.triplet),
            bond12=tuple(args.bond12),
            bond23=tuple(args.bond23) if args.bond23 else None,
            bin_width=args.bin_width,
            collect_angles=args.dump_angles is not None,
        )
    # IndexError covers a truncated header (a "Lattice" line with fewer than
    # three vector lines after it); OSError an unreadable input.
    except (FileNotFoundError, ValueError, IndexError, OSError) as error:
        print(f"rmc-triplets: {error}", file=sys.stderr)
        return 1

    writers: list[tuple[Path, object]] = [(output, lambda path: write_csv(path, config, result))]
    if args.dump_angles and result.angles is not None:
        writers.append((
            destinations["--dump-angles"],
            lambda path: path.write_text(
                "".join(f"{value:.6f}\n" for value in result.angles), encoding="utf-8"
            ),
        ))
    if args.plot:
        plot = destinations["--plot"]
        # The temporary name ends in .tmp + the suffix, so name the format
        # (no suffix: PNG, as --help says).
        writers.append((plot, lambda path: write_plot(path, result, plot.suffix[1:].lower() or "png")))
    try:
        write_all(writers)
    except (OSError, ValueError) as error:
        print(f"rmc-triplets: cannot write the outputs: {error}", file=sys.stderr)
        return 1

    label = "-".join(result.triplet)
    print(f"configuration: {config}")
    if result.parse_warning:
        # Atom lines the shared .rmc6f grammar skipped: the histogram covers
        # the atoms that remain, so say so instead of reporting a short model.
        print(f"rmc-triplets: warning: {result.parse_warning}", file=sys.stderr)
    print(f"triplet:       {label} (central {result.triplet[1]})")
    for name, end, unique, directed, mean in (
        ("bonds 1-2", result.triplet[0], result.unique_bonds12, result.bond12_count, result.mean_length12),
        ("bonds 2-3", result.triplet[2], result.unique_bonds23, result.bond23_count, result.mean_length23),
    ):
        # Physical bonds once each (with the B-centred bond-vector count when
        # the end element is the central one, as in the CSV header); the
        # coordination is the B-centred count per B, which sees a B-B bond
        # from both of its ends.
        print(
            f"{name}:     {bond_count_text(unique, directed, end, result.triplet[1])}"
            + (
                f" (mean length {mean:.4f} Ang; coordination "
                f"{directed / result.apex_count:.2f} per central atom)"
                if directed
                else ""
            )
        )
    print(
        f"angles:        {result.angle_count}"
        + (
            f"  (mean {result.mean_angle:.2f} deg, std {result.std_angle:.2f} deg)"
            if result.angle_count
            else ""
        )
    )
    print(f"histogram:     {output}")
    if args.dump_angles and result.angles is not None:
        print(f"angles list:   {destinations['--dump-angles']}")
    if args.plot:
        print(f"plot:          {destinations['--plot']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
