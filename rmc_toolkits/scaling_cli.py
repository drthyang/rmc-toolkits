# SPDX-License-Identifier: AGPL-3.0-or-later
# Copyright (C) 2026 Tsung-Han Yang

"""Auto StoG command line: automatic total-scattering scaling to RMC-ready files.

Drop-in replacement for an interactive classic-stog session, with the manual
"try again" loop replaced by :func:`rmc_toolkits.scaling.autoscale`::

    rmc-autoscale stog.inp                      # auto-fit (a, b); classic file names
    rmc-autoscale stog.inp --manual             # reproduce the stog.inp hand scaling
    rmc-autoscale --data sample.dat --qmin 0.6 --qmax 30 --formula SrTiO3

    python -m rmc_toolkits.scaling_cli stog.inp # same tool, module form

Reads the classic ``stog.inp`` (or direct arguments), fits the affine
correction ``S_corr = a*S_meas + b`` unless a fixed scaling is requested, and
writes the classic stog output family in the Fortran's own conventions —
scaled S(Q), unfiltered g(r), filtered S(Q), filtered g(r) (+ an r*[g(r)-1]
column), the ``ft.dat`` correction, and the RMCProfile-ready ``FK(Q)`` /
``GK(r)`` / ``D(r)`` — plus a provenance JSON with the full configuration and
fit diagnostics.

Safety: outputs default into an ``autoscale/`` directory next to the input, and
nothing is ever overwritten without ``--force`` — so the tool cannot silently
clobber the real STOG outputs a ``stog.inp`` typically sits beside. An output
that would land on the input data file or the ``stog.inp`` itself is refused
even with ``--force``, as are two outputs naming one file, an output path that
is a directory and an output folder that cannot be created -- all checked
before any computation. The family is written through temporary files and
renamed into place only once every write has succeeded, so a failure never
leaves a half-written family. Classic
low-r enforcement (the Fortran's final ripple removal) is applied to the RMC
files by default: at the ``stog.inp`` cutoff/first-peak window in ``stog.inp``
mode (parity), at ``--enforce-cutoff`` when given, and otherwise at the foot of
the detected first shell (:func:`rmc_toolkits.scaling.auto_enforcement_cutoff`,
never above a given r0); ``--no-enforce`` disables it. The honest
*pre*-enforcement low-r residual is always reported.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import uuid
from dataclasses import replace
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np

from . import __version__
from .parsers import (
    StogInput,
    read_dat_header,
    read_stog_inp,
    read_stog_xy,
    write_stog_xy,
)
from .scaling import (
    MIN_AUTO_WINDOW,
    R0_WINDOW_MARGIN,
    RHO0_SEED,
    ScalingConfig,
    ScalingResult,
    auto_enforcement_cutoff,
    autoscale,
    detect_first_peak_onset,
    diagnostics_summary,
    estimate_rho0,
    scale_pipeline,
)
from .scattering import faber_ziman, number_density_from_mass_density
from .transforms import (
    first_peak_zero,
    fq_to_gpdf,
    g_to_gk,
    gk_to_dr,
    gpdf_to_g,
    sq_to_fq,
)


class CliError(Exception):
    """User-facing error: rendered without a traceback, exit code 2."""


#: Output family, in write order: (logical key, stem-mode suffix, description).
_OUTPUTS = (
    ("sq_scaled", ".sq", "scaled S(Q), unfiltered"),
    ("gr_unfiltered", ".gr", "g(r), unfiltered transform (classic scale.gr)"),
    ("sq_filtered", "_ft.sq", "Fourier-filtered S(Q)"),
    ("gr_filtered", "_ft.gr", "filtered g(r) (+ r*[g(r)-1] column, classic scale_ft.gr)"),
    ("rmc_fq", "_rmc.fq", "FK(Q), barns (RMCProfile input)"),
    ("rmc_gr", "_rmc.gr", "Keen GK(r), barns (RMCProfile input)"),
    ("rmc_dr", "_rmc.dr", "D(r) (RMCProfile input)"),
)

#: Classic fixed-name Fourier-filter correction (stog and pystog both emit it).
#: Written on the data grid; the Fortran's extra sub-Qmin stub carries no data.
_FT_NAME = "ft.dat"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="rmc-autoscale",
        description=(
            "Auto StoG: automatic scale/offset determination for measured "
            "S(Q), STOG-compatible Fourier filter, and RMCProfile-ready "
            "outputs. Feed it a classic stog.inp, or a data file plus "
            "--qmin/--qmax and a composition source."
        ),
    )
    parser.add_argument(
        "stog_inp",
        nargs="?",
        default=None,
        metavar="stog.inp",
        help="classic stog input file (omit when using --data)",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    data = parser.add_argument_group("direct data mode (instead of a stog.inp)")
    data.add_argument("--data", metavar="FILE", help="S(Q) data file (2-3 columns; NaN padding ok)")
    data.add_argument("--qmin", type=float, help="fit/transform Q minimum (1/A)")
    data.add_argument("--qmax", type=float, help="fit/transform Q maximum (1/A)")
    data.add_argument(
        "--rho0",
        type=float,
        help="number density (1/A^3); default: NUMBER_DENSITY :: from the data "
        "header, else derived from --mass-density + --formula",
    )
    data.add_argument(
        "--mass-density",
        type=float,
        help="sample mass density (g/cm^3); with --formula this converts to "
        "rho0 via N_A (ADDIE convention)",
    )
    data.add_argument(
        "--estimate-rho0",
        action="store_true",
        help="self-consistent number density: iterate the density-limit fit "
        "until its amplitude agrees with the rho0-independent Q->0 "
        "Faber-Ziman amplitude (requires <b^2> via --b-sq-avg or --formula); "
        "the estimate replaces rho0 for the run. rho0 (--rho0 / header / "
        "--mass-density / stog.inp) seeds it; with no density source the seed "
        "is 0.05 1/A^3",
    )
    data.add_argument(
        "--b-avg-sq",
        type=float,
        help="<b>^2 in barns (the stog 'Faber-Ziman coefficient'); default: from --formula",
    )
    data.add_argument(
        "--b-sq-avg",
        type=float,
        help="<b^2> in barns, enables the Q->0 amplitude diagnostic; default: from --formula",
    )

    physics = parser.add_argument_group("physics and fit options")
    physics.add_argument(
        "--formula",
        help="chemical formula (e.g. SrTiO3) for <b>^2/<b^2> via the Sears table",
    )
    physics.add_argument("--r-cutoff", type=float, help="Fourier-filter r cutoff (A); default 1.0 / stog.inp")
    physics.add_argument(
        "--r0",
        type=float,
        help="closest interatomic approach (A); default: MINIMUM_DISTANCES :: header, "
        "then the stog.inp first-peak line (peak_rmin when its window starts below "
        "the cutoff, else peak_cutoff), else detected from the data",
    )
    physics.add_argument("--r-fit-min", type=float, help="low-r fit window minimum (A)")
    physics.add_argument("--r-fit-max", type=float, help="low-r fit window maximum (A)")
    physics.add_argument("--rmax", type=float, help="r-grid extent (A); default 50 / stog.inp")
    physics.add_argument("--nr", type=int, help="number of r points; default 5000 / stog.inp")
    physics.add_argument(
        "--lorch",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Lorch window (default: stog.inp flag, else off)",
    )
    physics.add_argument(
        "--c1-mode",
        choices=("sweep", "joint"),
        default="sweep",
        help="high-Q architecture: level-sweep anchored (default) or joint 2-dof fit",
    )
    physics.add_argument(
        "--amplitude",
        choices=("density", "fz"),
        default="density",
        help="amplitude criterion: low-r density limit (default), or 'fz' — "
        "subtract the measured high-Q level, scale Q->0 onto the Faber-Ziman "
        "limit S(0) = 1 - <b^2>/<b>^2, shift the level back to 1 "
        "(requires --b-sq-avg or --formula)",
    )
    physics.add_argument(
        "--fit-offset",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="fit the additive offset b (joint mode only; sweep mode ties b to the level)",
    )
    physics.add_argument(
        "--robust",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Huber IRLS re-weighting of the fit (default on)",
    )
    physics.add_argument(
        "--despike",
        action="store_true",
        help="drop rolling-median outliers before fitting (detector glitches; "
        "beware: also flags real Bragg maxima on crystalline data)",
    )
    physics.add_argument(
        "--low-q-correction",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="analytic correction for the omitted [0, Qmin] range (default on)",
    )
    physics.add_argument(
        "--sigma",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="1/sigma-weight the high-Q fit with the data file's third column (default on)",
    )

    manual = parser.add_argument_group("fixed scaling (skip the auto-fit)")
    manual.add_argument(
        "--manual",
        action="store_true",
        help="use the stog.inp yscale/yoffset unchanged (classic-stog parity run)",
    )
    manual.add_argument("--scale", type=float, help="fix a in S_corr = a*S + b (implies --manual)")
    manual.add_argument("--offset", type=float, help="fix b in S_corr = a*S + b (implies --manual)")

    enforce = parser.add_argument_group("classic low-r enforcement of the RMC outputs")
    enforce.add_argument(
        "--enforce",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Fortran-stog final ripple removal on the RMC files (default: on — "
        "at --enforce-cutoff when given, else the stog.inp cutoff and "
        "first-peak window, else automatically at the foot of the detected "
        "first shell, below its rising flank and never above a given r0; "
        "--no-enforce to disable)",
    )
    enforce.add_argument(
        "--enforce-cutoff",
        type=float,
        help="enforcement r cutoff (A); overrides the stog.inp / automatic cutoff",
    )
    enforce.add_argument(
        "--peak-window",
        type=float,
        nargs=2,
        metavar=("RMIN", "RMAX"),
        help="first-peak window kept below the cutoff (stog.inp line 22 semantics)",
    )

    out = parser.add_argument_group("output")
    out.add_argument(
        "--out-dir",
        metavar="DIR",
        help="output directory (default: 'autoscale/' next to the input file)",
    )
    out.add_argument(
        "--out-stem",
        metavar="STEM",
        help="name outputs STEM.sq, STEM.gr, STEM_ft.sq, ... instead of the "
        "stog.inp declared names (default in --data mode: the data file's stem)",
    )
    out.add_argument(
        "--force",
        action="store_true",
        help="overwrite existing output files (never the input data or stog.inp, "
        "never a directory, and never two outputs onto one file)",
    )
    return parser


def _json_safe(value: Any) -> Any:
    """Recursively convert numpy scalars/arrays and non-finite floats for JSON."""
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_json_safe(item) for item in value.tolist()]
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def usable_sigma(q: np.ndarray, sq: np.ndarray, sigma: Optional[np.ndarray]):
    """The sigma column only when clean, else None (CLI, API; JS ``usableSigma``).

    Any non-finite or non-positive sigma on a row with finite Q and S drops the
    whole column: a zero sigma would get a 1e12 weight and one NaN sigma turns
    every weight NaN, so a broken uncertainty column must not poison the fit.
    """
    if sigma is None:
        return None
    usable = np.isfinite(q) & np.isfinite(sq)
    if not np.all(np.isfinite(sigma[usable])) or np.any(sigma[usable] <= 0):
        return None
    return sigma


def _load_dataset(data_path: Path, use_sigma: bool):
    """Read (q, sq, sigma) from a STOG-style data file; sigma only when clean."""
    if not data_path.exists():
        raise CliError(f"data file not found: {data_path}")
    columns = read_stog_xy(data_path)
    q, sq = columns[0], columns[1]
    sigma = None
    if use_sigma and columns.shape[0] >= 3:
        sigma = usable_sigma(q, sq, columns[2])
    return q, sq, sigma


#: <b>^2 values closer than this (relative) are the same scattering-length set.
COEFFICIENT_RTOL = 0.02


def resolve_coefficients(
    *,
    b_avg_sq: Optional[float],
    b_avg_sq_source: Optional[str],
    b_sq_avg: Optional[float],
    formula: Optional[str],
) -> "dict[str, Any]":
    """<b>^2 and <b^2> from ONE consistent source (CLI and scaling API).

    Explicit values win. A ``formula`` (Sears neutron table) fills what is
    missing, but its <b^2> is paired with a <b>^2 from elsewhere (``stog.inp``
    or an explicit value) only when the two <b>^2 agree within
    :data:`COEFFICIENT_RTOL` — then the formula's ratio <b^2>/<b>^2 (the S(0)
    target) is kept on the configured <b>^2 scale. Otherwise the configured
    <b>^2 belongs to another radiation or normalization (e.g. <b>^2 = 1 for
    normalized x-ray data) and mixing would fabricate an S(0) target, so
    <b^2> stays unset (pass it explicitly). Returns ``{"b_avg_sq", "b_sq_avg",
    "b_avg_sq_source", "b_sq_avg_source", "warnings"}``.
    """
    warnings: list[str] = []
    b_sq_avg_source = "explicit" if b_sq_avg is not None else None
    formula = (formula or "").strip()
    if formula:
        coefficients = faber_ziman(formula)
        if b_avg_sq is None:
            b_avg_sq, b_avg_sq_source = coefficients.b_avg_sq_barn, f"formula {formula}"
        agree = abs(coefficients.b_avg_sq_barn - b_avg_sq) <= COEFFICIENT_RTOL * abs(b_avg_sq)
        if b_sq_avg is None:
            if agree:
                ratio = coefficients.b_sq_avg_barn / coefficients.b_avg_sq_barn
                b_sq_avg, b_sq_avg_source = b_avg_sq * ratio, f"formula {formula}"
            else:
                warnings.append(
                    f"<b>^2 from formula {formula} = {coefficients.b_avg_sq_barn:.6f} "
                    f"barn differs from the {b_avg_sq_source} value {b_avg_sq:.6f} "
                    "barn; using <b>^2 = "
                    f"{b_avg_sq:.6f} barn and NOT the formula's <b^2> = "
                    f"{coefficients.b_sq_avg_barn:.6f} barn (a pair from two sources "
                    "fabricates the S(0) target) — pass <b^2> explicitly (--b-sq-avg; "
                    "<Z^2>/<Z>^2 for normalized x-ray data) for the Q->0 criteria"
                )
        elif not agree:
            warnings.append(
                f"<b>^2 from formula {formula} = {coefficients.b_avg_sq_barn:.6f} "
                f"barn differs from the {b_avg_sq_source} value {b_avg_sq:.6f} barn; "
                f"using <b>^2 = {b_avg_sq:.6f} barn and the explicit <b^2> = "
                f"{b_sq_avg:.6f} barn"
            )
    return {
        "b_avg_sq": b_avg_sq,
        "b_sq_avg": b_sq_avg,
        "b_avg_sq_source": b_avg_sq_source,
        "b_sq_avg_source": b_sq_avg_source,
        "warnings": warnings,
    }


def refuse_failed_fit(result: ScalingResult) -> None:
    """Raise :class:`CliError` when an auto-fit returned a non-physical scale.

    ``autoscale`` flags ``a <= 0`` (or a non-finite scale) as
    ``provenance["fit_failure"]`` with ``converged=False``; such a result must
    not become RMCProfile input. Shared by the CLI and the scaling API.
    """
    failure = result.provenance.get("fit_failure")
    if failure:
        raise CliError(
            f"auto-fit failed: {failure}. No files were written. Check that the "
            "S(Q) is not sign-inverted or corrupted, pin the closest approach "
            "(--r0 / --r-fit-max), or use '--amplitude fz' when the composition "
            "is known"
        )


def stog_inp_closest_approach(inp: StogInput, r_cutoff: float) -> Optional[float]:
    """Closest-approach proxy from a classic stog.inp first-peak line (line 22).

    ``peak_cutoff peak_rmin peak_rmax`` zeroes g for ``r <= peak_cutoff``
    *except* inside ``[peak_rmin, peak_rmax]`` (Fortran ``first_peak_zero``
    semantics), so the region the expert asserts g = 0 ends at ``peak_rmin``
    when a genuine first-peak window starts below the cutoff
    (``0 < peak_rmin < peak_cutoff < ...``, ``peak_rmax > peak_rmin``) — the
    window exists precisely for first peaks that begin inside the cleanup
    radius — and at ``peak_cutoff`` otherwise (window outside ``[0, cutoff]``,
    or ``'1.0 0 0'``-style lines). Returned only when the default fit window it
    leaves above ``r_cutoff`` is at least :data:`MIN_AUTO_WINDOW` wide — the
    floor the automatic placement uses (else None: r0 is detected). A sliver
    used to pass: ``'1.46 0 0'`` at r_cutoff 1.0 pinned [1.2, 1.21] A and
    FeCoSn's scale came out 15 % low. Shared by the CLI, the API and (ported)
    the Auto StoG page.
    """
    candidate = float(inp.peak_cutoff)
    if 0.0 < inp.peak_rmin < inp.peak_cutoff and inp.peak_rmax > inp.peak_rmin:
        candidate = float(inp.peak_rmin)
    if candidate - R0_WINDOW_MARGIN - (r_cutoff + 0.2) >= MIN_AUTO_WINDOW:
        return candidate
    return None


def _default_r0(
    args: argparse.Namespace,
    header: dict,
    inp: Optional[StogInput],
    r_cutoff: float,
) -> Optional[float]:
    """r0 chain: flag > data-header MINIMUM_DISTANCES > usable stog.inp peak line."""
    if args.r0 is not None:
        return args.r0
    if "min_distance" in header:
        return float(header["min_distance"])
    if inp is not None:
        return stog_inp_closest_approach(inp, r_cutoff)
    return None


def _build_config(
    args: argparse.Namespace,
    inp: Optional[StogInput],
    header: dict,
) -> ScalingConfig:
    if inp is not None:
        qmin = args.qmin if args.qmin is not None else inp.qmin
        qmax = args.qmax if args.qmax is not None else inp.qmax
        rho0 = args.rho0 if args.rho0 is not None else inp.rho0
        b_avg_sq = args.b_avg_sq if args.b_avg_sq is not None else inp.b_avg_sq
        r_cutoff = args.r_cutoff if args.r_cutoff is not None else inp.r_cutoff
        rmax = args.rmax if args.rmax is not None else inp.rmax
        nr = args.nr if args.nr is not None else inp.nr
        lorch = args.lorch if args.lorch is not None else inp.lorch
    else:
        if args.qmin is None or args.qmax is None:
            raise CliError("--data mode requires --qmin and --qmax")
        qmin, qmax = args.qmin, args.qmax
        rho0 = args.rho0
        if rho0 is None:
            rho0 = header.get("number_density")
        if rho0 is None and args.mass_density is not None:
            if not (args.formula or "").strip():
                raise CliError("--mass-density needs --formula to convert to rho0")
            rho0 = number_density_from_mass_density(args.formula, args.mass_density)
        b_avg_sq = args.b_avg_sq
        r_cutoff = args.r_cutoff if args.r_cutoff is not None else 1.0
        rmax = args.rmax if args.rmax is not None else 50.0
        nr = args.nr if args.nr is not None else 5000
        lorch = bool(args.lorch)

    if args.b_avg_sq is not None:
        b_avg_sq_source = "--b-avg-sq"
    elif inp is not None:
        b_avg_sq_source = "stog.inp"
    else:
        b_avg_sq_source = None
    try:
        resolved = resolve_coefficients(
            b_avg_sq=b_avg_sq, b_avg_sq_source=b_avg_sq_source,
            b_sq_avg=args.b_sq_avg, formula=args.formula,
        )
    except ValueError as exc:
        raise CliError(str(exc)) from exc
    for warning in resolved["warnings"]:
        print(f"warning: {warning}", file=sys.stderr)
    b_avg_sq, b_sq_avg = resolved["b_avg_sq"], resolved["b_sq_avg"]
    if b_avg_sq is None:
        raise CliError("--data mode requires <b>^2: pass --b-avg-sq or --formula")
    if rho0 is None:
        if args.estimate_rho0 and b_sq_avg is not None:
            # No density source, but the self-consistency is requested: seed it
            # like the Auto StoG page (the estimate replaces the seed).
            rho0 = RHO0_SEED
            print(
                f"rho0: no density source — seeding the self-consistency at "
                f"{RHO0_SEED:g} 1/A^3",
                file=sys.stderr,
            )
        else:
            raise CliError(
                "number density unknown: pass --rho0, or --mass-density with "
                "--formula, or use a data file with a NUMBER_DENSITY :: header "
                "(or --estimate-rho0 with <b^2> from --formula / --b-sq-avg)"
            )
    if args.amplitude == "fz" and b_sq_avg is None:
        raise CliError(
            "--amplitude fz requires <b^2>: pass --b-sq-avg, or --formula when "
            "<b>^2 also comes from it"
        )

    try:
        config = ScalingConfig(
            qmin=float(qmin),
            qmax=float(qmax),
            rho0=float(rho0),
            b_avg_sq=float(b_avg_sq),
            b_sq_avg=None if b_sq_avg is None else float(b_sq_avg),
            r_cutoff=float(r_cutoff),
            r0=_default_r0(args, header, inp, float(r_cutoff)),
            fit_offset=args.fit_offset,
            r_fit_min=args.r_fit_min,
            r_fit_max=args.r_fit_max,
            rmax=float(rmax),
            nr=int(nr),
            lorch=bool(lorch),
            low_q_correction=args.low_q_correction,
            robust=args.robust,
            c1_mode=args.c1_mode,
            amplitude_criterion=args.amplitude,
            despike=args.despike,
        )
        config.r_fit_window  # validate now, with CLI-error rendering
    except ValueError as exc:
        raise CliError(f"invalid configuration: {exc}") from exc
    return config


def _resolve_enforcement(
    args: argparse.Namespace, inp: Optional[StogInput]
) -> Optional[tuple[float, float, float]]:
    """Return (cutoff, peak_rmin, peak_rmax) or None when enforcement is off."""
    if args.enforce is False:
        if args.enforce_cutoff is not None or args.peak_window is not None:
            raise CliError("--no-enforce contradicts --enforce-cutoff/--peak-window")
        return None
    cutoff = args.enforce_cutoff
    if cutoff is None and inp is not None:
        cutoff = inp.peak_cutoff
    if cutoff is None:
        if args.peak_window is not None:
            raise CliError("--peak-window requires --enforce-cutoff (or a stog input)")
        return None  # data mode: resolved post-run from the detected r0
    if args.peak_window is not None:
        peak_rmin, peak_rmax = args.peak_window
    elif inp is not None and args.enforce_cutoff is None:
        peak_rmin, peak_rmax = inp.peak_rmin, inp.peak_rmax
    else:
        peak_rmin = peak_rmax = cutoff  # degenerate window: flat replacement
    return float(cutoff), float(peak_rmin), float(peak_rmax)


def _same_file(target: Path, source: Path) -> bool:
    """True when ``target`` names ``source`` (same path, symlink or hard link)."""
    try:
        if target.exists() and source.exists():
            return os.path.samefile(target, source)
        return target.resolve() == source.resolve()
    except OSError:
        return False


def _resolve_targets(
    args: argparse.Namespace,
    inp: Optional[StogInput],
    inp_path: Optional[Path],
    data_path: Path,
) -> "dict[str, Path]":
    anchor = inp_path if inp_path is not None else data_path
    out_dir = Path(args.out_dir) if args.out_dir else anchor.parent / "autoscale"
    stem = args.out_stem
    if stem is None and inp is None:
        stem = data_path.stem
    targets: dict[str, Path] = {}
    if stem is None:
        declared = (
            inp.out_sq, inp.out_gr, inp.out_ft_sq, inp.out_ft_gr,
            inp.out_rmc_fq, inp.out_rmc_gr, inp.out_rmc_dr,
        )
        for (key, _, _), name in zip(_OUTPUTS, declared):
            targets[key] = out_dir / name
        json_name = f"{inp_path.stem}_provenance.json"
    else:
        for key, suffix, _ in _OUTPUTS:
            targets[key] = out_dir / f"{stem}{suffix}"
        json_name = f"{stem}_provenance.json"
    targets["ft_correction"] = out_dir / _FT_NAME
    targets["provenance"] = out_dir / json_name

    # The inputs are never an output, --force or not: in --data mode the
    # default stem makes <out-dir>/<stem>.sq the measured file itself when
    # --out-dir is its own folder, and a stog.inp may declare an output name
    # equal to its data file. Overwriting either destroys the measured data
    # (and a rerun would silently re-scale already-scaled data).
    inputs = [path for path in (data_path, inp_path) if path is not None]
    clashes = [
        str(target)
        for target in targets.values()
        if any(_same_file(target, source) for source in inputs)
    ]
    if clashes:
        raise CliError(
            "output would overwrite an input file (never allowed, even with "
            "--force; pick --out-dir/--out-stem):\n  " + "\n  ".join(clashes)
        )

    _check_targets_writable(targets)

    if not args.force:
        existing = [str(path) for path in targets.values() if path.exists()]
        if existing:
            raise CliError(
                "refusing to overwrite existing outputs (use --force, or pick "
                "--out-dir/--out-stem):\n  " + "\n  ".join(existing)
            )
    return targets


def _target_key(path: Path) -> str:
    """Identity of an output path for the distinctness check.

    Case-folded, so two names differing only in case count as one file: they
    are one file on the default macOS and Windows filesystems, and nobody wants
    both as separate outputs.
    """
    return os.path.normcase(str(path.resolve())).casefold()


def _check_targets_writable(targets: "dict[str, Path]") -> None:
    """Refuse, before any computation or write, a family that cannot be written whole.

    Every target must be its own file: no two targets may name the same file
    (a stog.inp declaring FK(Q) as ``ft.dat`` used to exit 0 with the RMCProfile
    input overwritten by the Fourier-filter correction), none may be an existing
    directory, and its nearest existing ancestor must be a directory (the
    writer creates the missing folders). ``--force`` never relaxes these.
    """
    items = list(targets.items())
    duplicates = [
        f"{first_key} and {second_key} -> {second}"
        for index, (first_key, first) in enumerate(items)
        for second_key, second in items[index + 1:]
        if _target_key(first) == _target_key(second) or _same_file(first, second)
    ]
    if duplicates:
        raise CliError(
            "two outputs would be the same file (never allowed, even with --force; "
            "give them distinct names):\n  " + "\n  ".join(duplicates)
        )

    directories = [str(path) for path in targets.values() if path.is_dir()]
    if directories:
        raise CliError(
            "output path is a directory (never replaced, even with --force):\n  "
            + "\n  ".join(directories)
        )

    blocked = []
    for path in targets.values():
        ancestor = path.parent
        while not ancestor.exists() and ancestor != ancestor.parent:
            ancestor = ancestor.parent
        if not ancestor.is_dir():
            blocked.append(f"{path} (its folder {ancestor} is not a directory)")
    if blocked:
        raise CliError(
            "cannot create the output folder, a path component is not a directory:\n  "
            + "\n  ".join(blocked)
        )


def _temporary_sibling(path: Path) -> Path:
    """A fresh hidden name beside ``path``: same folder, so the rename is atomic."""
    return path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex[:12]}.tmp")


def _write_family_atomically(targets: "dict[str, Path]", writers: "dict[str, Any]") -> None:
    """Write every target through a temporary sibling, then rename them all.

    ``writers`` maps a target key to a callable writing that file's content to
    the path it is given. Nothing is renamed into place until every write has
    succeeded, so a failure (disk full, permissions, an I/O error) removes the
    temporaries and leaves the previous files -- or none -- exactly as they were:
    never a half-written or mixed family.
    """
    for path in targets.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    temporaries: dict[str, Path] = {}
    try:
        for key, write in writers.items():
            temporary = _temporary_sibling(targets[key])
            temporaries[key] = temporary
            write(temporary)
        for key, temporary in temporaries.items():
            os.replace(temporary, targets[key])
    except BaseException:
        for temporary in temporaries.values():
            try:
                temporary.unlink()
            except OSError:
                pass
        raise


def _write_outputs(
    result: ScalingResult,
    config: ScalingConfig,
    targets: "dict[str, Path]",
    enforcement: Optional[tuple[float, float, float]],
    payload: "dict[str, Any]",
) -> None:
    # The unfiltered transform (classic scale.gr) is not part of ScalingResult;
    # recompute it with the same discretization the filter used internally.
    gpdf = fq_to_gpdf(
        result.q,
        sq_to_fq(result.q, result.sq_scaled),
        result.r,
        lorch=config.lorch,
        low_q_correction=config.low_q_correction,
        s0_target=config.effective_s0_target,
    )
    g_unfiltered = gpdf_to_g(result.r, gpdf, config.rho0)

    if enforcement is not None:
        cutoff, peak_rmin, peak_rmax = enforcement
        g_final = first_peak_zero(
            result.r,
            result.g_filtered,
            cutoff=cutoff,
            peak_rmin=peak_rmin,
            peak_rmax=peak_rmax,
        )
        gk_out = g_to_gk(g_final, config.b_avg_sq)
        dr_out = gk_to_dr(result.r, gk_out, config.rho0)
    else:
        gk_out, dr_out = result.gk, result.d_r

    label = f"rmc-autoscale {__version__}: a={result.a:.8g} b={result.b:.8g}"

    def xy(x, y, **kwargs):
        return lambda path: write_stog_xy(path, x, y, title=label, **kwargs)

    def provenance(path: Path) -> None:
        with path.open("w", encoding="utf-8") as handle:
            json.dump(_json_safe(payload), handle, indent=2)
            handle.write("\n")

    # Classic stog conventions (verified against the Fortran runs in
    # data/stog_tests: scale.gr column 2 is g(r), oscillating about 1, and
    # scale_ft.gr column 3 is exactly r*[g(r) - 1]).
    _write_family_atomically(
        targets,
        {
            "sq_scaled": xy(result.q, result.sq_scaled),
            "gr_unfiltered": xy(result.r, g_unfiltered),
            "sq_filtered": xy(result.q, result.sq_filtered),
            "gr_filtered": xy(result.r, result.g_filtered, extra=result.r * (result.g_filtered - 1.0)),
            "rmc_fq": xy(result.q, result.fk),
            "rmc_gr": xy(result.r, gk_out),
            "rmc_dr": xy(result.r, dr_out),
            "ft_correction": xy(result.q, result.sq_ft),
            "provenance": provenance,
        },
    )


def _print_report(
    result: ScalingResult,
    summary: "dict[str, Any]",
    targets: "dict[str, Path]",
    reference: Optional[tuple[float, float]],
    manual: bool,
    enforcement: Optional[tuple[float, float, float]],
    n_points: int,
    enforcement_note: Optional[str] = None,
    config: Optional[ScalingConfig] = None,
) -> None:
    mode = "manual (fixed a, b)" if manual else f"auto ({summary.get('c1_mode', 'fit')})"
    print(f"Auto StoG (rmc-toolkits {__version__})")
    print(f"  mode      : {mode}")
    print(f"  data      : {n_points} S(Q) points used")
    if config is not None:
        # The coefficients actually in effect (after --formula / stog.inp /
        # explicit-value resolution), with the S(0) target they imply.
        pair = f"<b>^2 = {config.b_avg_sq:.6g} barn, <b^2> = " + (
            "not set" if config.b_sq_avg is None else f"{config.b_sq_avg:.6g} barn"
        )
        print(f"  coeffs    : {pair}, S(0) target = {config.effective_s0_target:.4g}")
    line = f"  result    : a = {result.a:.6g}, b = {result.b:.6g}"
    if reference is not None and not manual:
        line += f"   [stog.inp hand values: a = {reference[0]:.6g}, b = {reference[1]:.6g}]"
    print(line)
    if not manual:
        status = "yes" if result.converged else "NO"
        print(f"  converged : {status} ({result.iterations} iterations)")
    if "level" in summary and summary["level"] is not None:
        level_sigma = summary.get("level_uncertainty")
        spread = "" if level_sigma is None or not math.isfinite(level_sigma) else f" +/- {level_sigma:.2g}"
        window = summary.get("level_window")
        where = f" over Q = [{window[0]:.2f}, {window[1]:.2f}]" if window else ""
        print(f"  level     : {summary['level']:.6g}{spread}{where}")
    print(f"  C1 tail   : filtered S(Q) mean = {summary['c1_tail_mean']:.6g} (target 1)")
    print(
        f"  low-r rms : {summary['low_r_rms_pre_enforcement']:.4g} "
        "(pre-enforcement, g-space, target 0)"
    )
    if not summary["density_limit_satisfied"]:
        print(
            "  density limit NOT satisfied: the absolute scale is not "
            "recoverable from this data alone (missing low-Q information); "
            "validate the scale externally"
        )
    if summary.get("rmax_beyond_alias_limit"):
        print(
            f"  WARNING   : rmax = {config.rmax if config is not None else float('nan'):g} A "
            f"exceeds the aliasing limit pi/dQ = {summary['r_alias_limit']:.4g} A of the "
            "coarsest S(Q) step: G(r)/D(r) beyond it are folded (a negated mirror image "
            "on a uniform grid) or corrupted by coarse steps (log binning, despike gaps) "
            "— lower --rmax, or use finer, uniformly binned data"
        )
    if summary.get("r0_detected") is not None:
        refined = " (fit window refined)" if summary.get("window_refined") else ""
        print(f"  r0 (data) : first-shell onset detected at {summary['r0_detected']:.2f} A{refined}")
        if summary.get("first_shell_below_r0"):
            print(
                "  WARNING   : the first shell starts below the given r0 "
                f"({summary['r0_detected']:.2f} A); the low-r window "
                f"[{summary['r_fit_window'][0]:g}, {summary['r_fit_window'][1]:g}] A may "
                "cut into it — check r0 (--r0 / MINIMUM_DISTANCES / stog.inp)"
            )
    if "amplitude_concordance" in summary:
        verdict = "concordant" if summary["amplitudes_concordant"] else "DISCORDANT"
        print(
            f"  amplitude concordance: a_fz/a = "
            f"{summary['amplitude_concordance']:.3f} ({verdict})"
        )
    if summary.get("a_fz_reliable") is True:
        # The flag is statistical: a systematic head bias passes it (Mn3Sn
        # 55537 / 54139: reliable a_fz drifting 11 -> 6 / 16 -> 24 with Qmin).
        also = "; see also the concordance above" if "amplitude_concordance" in summary else ""
        print(
            f"  Q->0 amplitude: a_fz = {summary['a_fz']:.4g} (relative error "
            f"{summary['a_fz_rel_se']:.0%}, resolved) — necessary, not sufficient: a biased "
            f"low-Q head passes too; re-run at a few --qmin values to check a_fz is "
            f"stable{also}"
        )
    if summary.get("a_fz_reliable") is False:
        use = "the concordance" if "amplitude_concordance" in summary else "it as the scale"
        print(
            f"  WARNING   : the Q->0 Faber-Ziman amplitude a_fz = {summary['a_fz']:.4g} is "
            f"ill-conditioned (relative error {summary['a_fz_rel_se']:.0%}: S_meas(0) - "
            "level is not resolved from its uncertainty — Bragg-contaminated or long "
            f"low-Q head); do not trust {use}"
        )
    if enforcement is not None:
        cutoff, peak_rmin, peak_rmax = enforcement
        where = (
            enforcement_note
            if enforcement_note
            else f"first-peak window [{peak_rmin:g}, {peak_rmax:g}]"
        )
        print(f"  enforcement: RMC outputs hard-set below r = {cutoff:.4g} A ({where})")
    elif enforcement_note:
        print(f"  enforcement: {enforcement_note}")
    print(f"Outputs -> {targets['provenance'].parent}")
    for key, _, description in _OUTPUTS:
        print(f"  {targets[key].name:<28s} {description}")
    print(f"  {targets['ft_correction'].name:<28s} Fourier-filter correction (classic ft.dat)")
    print(f"  {targets['provenance'].name:<28s} configuration + fit diagnostics (JSON)")


def main(argv: Optional[Sequence[str]] = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        if (args.stog_inp is None) == (args.data is None):
            raise CliError("pass exactly one input: a stog.inp path, or --data FILE")

        inp: Optional[StogInput] = None
        inp_path: Optional[Path] = None
        if args.stog_inp is not None:
            inp_path = Path(args.stog_inp)
            if not inp_path.exists():
                raise CliError(f"stog input file not found: {inp_path}")
            inp = read_stog_inp(inp_path)
            data_path = inp_path.parent / inp.data_file
        else:
            data_path = Path(args.data)

        q, sq, sigma = _load_dataset(data_path, use_sigma=args.sigma)
        header: dict = {}
        try:
            header = read_dat_header(data_path)
        except OSError:
            pass

        config = _build_config(args, inp, header)
        enforcement = _resolve_enforcement(args, inp)
        enforcement_source = None
        if enforcement is not None:
            enforcement_source = "user" if args.enforce_cutoff is not None else "stog.inp"
        targets = _resolve_targets(args, inp, inp_path, data_path)

        manual = args.manual or args.scale is not None or args.offset is not None
        if manual and args.amplitude != "density":
            raise CliError(
                "--amplitude selects the auto-fit criterion; it cannot be "
                "combined with --manual/--scale/--offset"
            )

        rho0_estimate = None
        if args.estimate_rho0:
            if manual:
                raise CliError(
                    "--estimate-rho0 drives the auto-fit; it cannot be "
                    "combined with --manual/--scale/--offset"
                )
            if config.b_sq_avg is None:
                raise CliError(
                    "--estimate-rho0 requires <b^2>: pass --b-sq-avg, or --formula "
                    "when <b>^2 also comes from it"
                )
            rho0_estimate = estimate_rho0(q, sq, config, sigma=sigma)
            if not rho0_estimate["converged"]:
                if rho0_estimate.get("stopped"):
                    reason = (
                        "the iteration stopped at a density the auto-scale cannot "
                        f"fit ({rho0_estimate['stopped']})"
                    )
                elif rho0_estimate.get("reason"):
                    reason = rho0_estimate["reason"]
                else:
                    reason = (
                        "the density-limit and Q->0 Faber-Ziman amplitudes "
                        "disagree at every density — typically data missing "
                        "structure below Qmin"
                    )
                raise CliError(
                    "rho0 self-consistency did not converge (final concordance "
                    f"{rho0_estimate['concordance']:.3g}): {reason}. Set rho0 "
                    "explicitly (--rho0 / --mass-density / data header) and "
                    "consider '--amplitude fz' for the scale"
                )
            config = replace(config, rho0=float(rho0_estimate["rho0"]))
            note = (
                "; long Q->0 extrapolation — treat as a starting point"
                if rho0_estimate["extrapolated"]
                else ""
            )
            if rho0_estimate.get("a_fz_reliable") is False:
                note += (
                    "; WARNING: its Faber-Ziman anchor is ill-conditioned (a_fz "
                    f"relative error {rho0_estimate['a_fz_rel_se']:.0%})"
                )
            print(
                f"rho0 self-consistency: {rho0_estimate['rho0']:.6f} 1/A^3 "
                f"(concordance {rho0_estimate['concordance']:.4f}, "
                f"{rho0_estimate['iterations']} passes{note})"
            )

        if manual:
            if args.scale is not None:
                a = args.scale
            elif inp is not None:
                a = inp.a
            else:
                raise CliError("fixed scaling in --data mode requires --scale")
            if args.offset is not None:
                b = args.offset
            elif inp is not None and args.scale is None:
                b = inp.b
            else:
                b = 0.0
            result = scale_pipeline(q, sq, config, float(a), float(b))
        else:
            result = autoscale(q, sq, config, sigma=sigma)
            refuse_failed_fit(result)

        # No explicit cutoff (data mode) and enforcement not refused: enforce
        # automatically at the FOOT of the first shell, below its rising flank
        # (auto_enforcement_cutoff), so no first-shell signal is removed.
        auto_note: Optional[str] = None
        if enforcement is None and args.enforce is not False:
            r0_detected = result.provenance.get("r0_detected")
            if r0_detected is None:
                r0_detected = detect_first_peak_onset(
                    result.r, result.g_filtered, config.qmax,
                    search_min=config.r_cutoff + 0.3,
                )
                if r0_detected is not None:
                    result.provenance["r0_detected"] = float(r0_detected)
            cutoff = auto_enforcement_cutoff(
                result.r, result.g_filtered, config, onset=r0_detected
            )
            if cutoff is not None:
                enforcement = (cutoff,) * 3
                # Name what anchored the cutoff: the detected onset, or the
                # given r0 when it caps the onset or no shell was detected.
                if config.r0 is not None and (
                    r0_detected is None or float(config.r0) < float(r0_detected)
                ):
                    enforcement_source = "auto (given r0)"
                    auto_note = (
                        f"automatic: anchored on the given r0 {config.r0:g} A "
                        + (
                            "(no shell detected)"
                            if r0_detected is None
                            else f"(below the detected onset {r0_detected:.2f} A)"
                        )
                    )
                else:
                    enforcement_source = "auto (first-shell foot)"
                    auto_note = (
                        f"automatic: foot of the first shell, onset {r0_detected:.2f} A"
                    )
            else:
                auto_note = (
                    "none: no first shell detected to anchor an automatic "
                    "cutoff (pass --enforce-cutoff to enforce)"
                )

        summary = diagnostics_summary(result, config)
        summary["c1_mode"] = result.provenance.get("c1_mode_effective", "manual")
        if not manual and config.amplitude_criterion == "fz":
            summary["c1_mode"] += ", FZ-limit amplitude"
        payload = {
            "tool": "rmc-autoscale",
            "rmc_toolkits_version": __version__,
            "argv": argv,
            "stog_inp": None if inp_path is None else str(inp_path),
            "data_file": str(data_path),
            "stog_inp_reference": None
            if inp is None
            else {"a": inp.a, "b": inp.b, "yscale": inp.yscale, "yoffset": inp.yoffset},
            "enforcement": None
            if enforcement is None
            else {
                **dict(zip(("cutoff", "peak_rmin", "peak_rmax"), enforcement)),
                "source": enforcement_source,
            },
            "outputs": {key: str(path) for key, path in targets.items()},
            "rho0_estimate": rho0_estimate,
            "diagnostics": summary,
            "provenance": result.provenance,
        }

        _write_outputs(result, config, targets, enforcement, payload)
        reference = None if inp is None else (inp.a, inp.b)
        _print_report(
            result, summary, targets, reference, manual, enforcement,
            n_points=int(result.provenance["n_q_points"]),
            enforcement_note=auto_note,
            config=config,
        )
        return 0
    except CliError as exc:
        print(f"rmc-autoscale: error: {exc}", file=sys.stderr)
        return 2
    except (ValueError, NotImplementedError, OSError) as exc:
        print(f"rmc-autoscale: error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover - direct module invocation
    raise SystemExit(main())
