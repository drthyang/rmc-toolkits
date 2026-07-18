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
writes the classic stog output family — scaled S(Q), unfiltered g(r)-1,
filtered S(Q), filtered g(r)-1 (+ D(r) column), and the RMCProfile-ready
``FK(Q)`` / ``GK(r)`` / ``D(r)`` — plus a provenance JSON with the full
configuration and fit diagnostics.

Safety: outputs default into an ``autoscale/`` directory next to the input, and
nothing is ever overwritten without ``--force`` — so the tool cannot silently
clobber the real STOG outputs a ``stog.inp`` typically sits beside. Classic
low-r enforcement (the Fortran's final ripple removal) is applied to the RMC
files by default in ``stog.inp`` mode for parity; the honest *pre*-enforcement
low-r residual is always reported.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
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
    ScalingConfig,
    ScalingResult,
    autoscale,
    diagnostics_summary,
    scale_pipeline,
)
from .scattering import faber_ziman
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
    ("gr_unfiltered", ".gr", "g(r) - 1, unfiltered transform"),
    ("sq_filtered", "_ft.sq", "Fourier-filtered S(Q)"),
    ("gr_filtered", "_ft.gr", "filtered g(r) - 1 (+ 4*pi*rho0*r*[g-1] column)"),
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
        help="number density (1/A^3); default: NUMBER_DENSITY :: from the data header",
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
        "then the stog.inp first-peak line",
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
        help="Fortran-stog final ripple removal on the RMC files "
        "(default: on in stog.inp mode, off in --data mode)",
    )
    enforce.add_argument("--enforce-cutoff", type=float, help="enforcement r cutoff (A)")
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
    out.add_argument("--force", action="store_true", help="overwrite existing output files")
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


def _load_dataset(data_path: Path, use_sigma: bool):
    """Read (q, sq, sigma) from a STOG-style data file; sigma only when clean."""
    if not data_path.exists():
        raise CliError(f"data file not found: {data_path}")
    columns = read_stog_xy(data_path)
    q, sq = columns[0], columns[1]
    sigma = None
    if use_sigma and columns.shape[0] >= 3:
        sigma = columns[2]
        usable = np.isfinite(q) & np.isfinite(sq)
        if not np.all(np.isfinite(sigma[usable])) or np.any(sigma[usable] <= 0):
            sigma = None  # a broken uncertainty column must not poison the fit
    return q, sq, sigma


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
        # peak_cutoff is the ripple-cleanup radius and peak_rmin the first-peak
        # start; the larger is the better closest-approach proxy — but only
        # when it leaves a non-empty default fit window above r_cutoff.
        candidate = max(inp.peak_cutoff, inp.peak_rmin)
        if candidate - 0.25 > r_cutoff + 0.2:
            return candidate
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
        if rho0 is None:
            raise CliError(
                "number density unknown: pass --rho0 or use a data file with a "
                "NUMBER_DENSITY :: header"
            )
        b_avg_sq = args.b_avg_sq
        r_cutoff = args.r_cutoff if args.r_cutoff is not None else 1.0
        rmax = args.rmax if args.rmax is not None else 50.0
        nr = args.nr if args.nr is not None else 5000
        lorch = bool(args.lorch)

    b_sq_avg = args.b_sq_avg
    if args.formula:
        coefficients = faber_ziman(args.formula)
        if b_sq_avg is None:
            b_sq_avg = coefficients.b_sq_avg_barn
        if b_avg_sq is None:
            b_avg_sq = coefficients.b_avg_sq_barn
        elif abs(coefficients.b_avg_sq_barn - b_avg_sq) > 0.02 * abs(b_avg_sq):
            print(
                f"warning: <b>^2 from --formula {args.formula} = "
                f"{coefficients.b_avg_sq_barn:.6f} barn differs from the "
                f"configured {b_avg_sq:.6f} barn; keeping the configured value",
                file=sys.stderr,
            )
    if b_avg_sq is None:
        raise CliError("--data mode requires <b>^2: pass --b-avg-sq or --formula")
    if args.amplitude == "fz" and b_sq_avg is None:
        raise CliError("--amplitude fz requires <b^2>: pass --b-sq-avg or --formula")

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
        if args.enforce or args.peak_window is not None:
            raise CliError("enforcement in --data mode requires --enforce-cutoff")
        return None
    if args.peak_window is not None:
        peak_rmin, peak_rmax = args.peak_window
    elif inp is not None and args.enforce_cutoff is None:
        peak_rmin, peak_rmax = inp.peak_rmin, inp.peak_rmax
    else:
        peak_rmin = peak_rmax = cutoff  # degenerate window: flat replacement
    return float(cutoff), float(peak_rmin), float(peak_rmax)


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

    if not args.force:
        existing = [str(path) for path in targets.values() if path.exists()]
        if existing:
            raise CliError(
                "refusing to overwrite existing outputs (use --force, or pick "
                "--out-dir/--out-stem):\n  " + "\n  ".join(existing)
            )
    return targets


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
    gm1 = result.g_filtered - 1.0
    write_stog_xy(targets["sq_scaled"], result.q, result.sq_scaled, title=label)
    write_stog_xy(targets["gr_unfiltered"], result.r, g_unfiltered - 1.0, title=label)
    write_stog_xy(targets["sq_filtered"], result.q, result.sq_filtered, title=label)
    write_stog_xy(
        targets["gr_filtered"],
        result.r,
        gm1,
        title=label,
        extra=4.0 * np.pi * config.rho0 * result.r * gm1,
    )
    write_stog_xy(targets["rmc_fq"], result.q, result.fk, title=label)
    write_stog_xy(targets["rmc_gr"], result.r, gk_out, title=label)
    write_stog_xy(targets["rmc_dr"], result.r, dr_out, title=label)
    write_stog_xy(targets["ft_correction"], result.q, result.sq_ft, title=label)
    with targets["provenance"].open("w", encoding="utf-8") as handle:
        json.dump(_json_safe(payload), handle, indent=2)
        handle.write("\n")


def _print_report(
    result: ScalingResult,
    summary: "dict[str, Any]",
    targets: "dict[str, Path]",
    reference: Optional[tuple[float, float]],
    manual: bool,
    enforcement: Optional[tuple[float, float, float]],
    n_points: int,
) -> None:
    mode = "manual (fixed a, b)" if manual else f"auto ({summary.get('c1_mode', 'fit')})"
    print(f"Auto StoG (rmc-toolkits {__version__})")
    print(f"  mode      : {mode}")
    print(f"  data      : {n_points} S(Q) points used")
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
    if "amplitude_concordance" in summary:
        verdict = "concordant" if summary["amplitudes_concordant"] else "DISCORDANT"
        print(
            f"  amplitude concordance: a_fz/a = "
            f"{summary['amplitude_concordance']:.3f} ({verdict})"
        )
    if enforcement is not None:
        cutoff, peak_rmin, peak_rmax = enforcement
        print(
            f"  enforcement: RMC outputs hard-set below r = {cutoff:g} A "
            f"(first-peak window [{peak_rmin:g}, {peak_rmax:g}])"
        )
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
        targets = _resolve_targets(args, inp, inp_path, data_path)

        manual = args.manual or args.scale is not None or args.offset is not None
        if manual and args.amplitude != "density":
            raise CliError(
                "--amplitude selects the auto-fit criterion; it cannot be "
                "combined with --manual/--scale/--offset"
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
            else dict(zip(("cutoff", "peak_rmin", "peak_rmax"), enforcement)),
            "outputs": {key: str(path) for key, path in targets.items()},
            "diagnostics": summary,
            "provenance": result.provenance,
        }

        targets["provenance"].parent.mkdir(parents=True, exist_ok=True)
        _write_outputs(result, config, targets, enforcement, payload)
        reference = None if inp is None else (inp.a, inp.b)
        _print_report(
            result, summary, targets, reference, manual, enforcement,
            n_points=int(result.provenance["n_q_points"]),
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
