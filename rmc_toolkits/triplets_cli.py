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
import sys
from pathlib import Path

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


def write_plot(path: Path, result: BondAngleDistribution) -> None:
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
    figure.savefig(path)


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
        help="overwrite existing output files instead of refusing",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        config = resolve_config(args.config)
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

    output = Path(args.output) if args.output else Path.cwd() / default_output_name(
        config, result.triplet
    )
    destinations = [output]
    if args.dump_angles:
        destinations.append(Path(args.dump_angles))
    if args.plot:
        destinations.append(Path(args.plot))
    if not args.force:
        existing = [str(path) for path in destinations if path.exists()]
        if existing:
            print(
                "rmc-triplets: refusing to overwrite "
                + ", ".join(existing)
                + "; pass --force to replace",
                file=sys.stderr,
            )
            return 1
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        write_csv(output, config, result)
    except OSError as error:
        print(f"rmc-triplets: cannot write {output}: {error}", file=sys.stderr)
        return 1

    label = "-".join(result.triplet)
    print(f"configuration: {config}")
    if result.parse_warning:
        # Atom lines the shared .rmc6f grammar skipped: the histogram covers
        # the atoms that remain, so say so instead of reporting a short model.
        print(f"rmc-triplets: warning: {result.parse_warning}", file=sys.stderr)
    print(f"triplet:       {label} (central {result.triplet[1]})")
    for name, unique, directed, mean in (
        ("bonds 1-2", result.unique_bonds12, result.bond12_count, result.mean_length12),
        ("bonds 2-3", result.unique_bonds23, result.bond23_count, result.mean_length23),
    ):
        # Bonds once each; "per central atom" is the B-centred count per B
        # (the coordination), which counts a B-B bond at both of its ends.
        print(
            f"{name}:     {unique} "
            + (
                f" (mean length {mean:.4f} Ang, "
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

    try:
        if args.dump_angles and result.angles is not None:
            angles_path = Path(args.dump_angles)
            angles_path.parent.mkdir(parents=True, exist_ok=True)
            angles_path.write_text(
                "".join(f"{value:.6f}\n" for value in result.angles), encoding="utf-8"
            )
            print(f"angles list:   {angles_path}")

        if args.plot:
            plot_path = Path(args.plot)
            plot_path.parent.mkdir(parents=True, exist_ok=True)
            write_plot(plot_path, result)
            print(f"plot:          {plot_path}")
    except OSError as error:
        print(f"rmc-triplets: cannot write output: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
