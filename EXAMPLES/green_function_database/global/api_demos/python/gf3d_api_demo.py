#!/usr/bin/env python
"""A runnable tour of the gf3d Python API, on the global example database.

Four things, in order, each printing what it measured:

  1. open the database and say what is in it
  2. extract seismograms and their ten partial derivatives, in memory
  3. check the identity the moment-tensor partials satisfy exactly
  4. use a centroid partial as a derivative -- predict the seismogram of a
     relocated source, and show the prediction improves quadratically as the
     relocation shrinks

Step 4 is the one worth reading. Everything before it can be had from
`xgf3d` and a pile of SAC files; what the in-memory API is for is the
inversion loop, where the partials are Frechet derivatives to be used, not
traces to be written. A first-order Taylor step is the smallest honest
demonstration that they are derivatives at all: get the sign, the units or
the scaling wrong and the prediction does not converge, whatever the traces
look like.

The script writes one figure and changes nothing else. It is not part of
the Snakemake workflow -- nothing in the Snakefile refers to this directory
-- so run it by hand, after the workflow has built `../../GFDB`.

Usage
-----
  cd EXAMPLES/green_function_database/global/api_demos/python
  ../../../.venv/bin/python gf3d_api_demo.py [--station NET.STA] [--component N|E|Z]
                                          [--step DEGREES] [--png FILE] [--no-plot]

`--station` and `--component` choose what the figure draws; the numbers are
computed for every station either way. `--step` is the relocation used in
step 4, in degrees of latitude.

Requires numpy, matplotlib (for the figure), and a built `lib/libgf3d.so`
-- see `utils/green_function/gf3d/README.md`.
"""

from __future__ import annotations

import argparse
import copy
import sys
import time
from pathlib import Path

import numpy as np

# The package normally comes from an install (`pip install -e
# utils/green_function`) or from PYTHONPATH. The example's own virtual
# environment has neither, so that a reader who has only run `uv sync` can
# still execute this file, fall back to the tree this script lives in.
try:
    import gf3d
except ImportError:
    # .../EXAMPLES/green_function_database/global/api_demos/python/this.py
    _repo = Path(__file__).resolve().parents[5]
    sys.path.insert(0, str(_repo / "utils" / "green_function"))
    try:
        import gf3d
    except ImportError:
        raise SystemExit(
            f"cannot import gf3d, and it is not at {_repo / 'utils' / 'green_function'}\n"
            "either run this script from inside a specfem3d_globe tree, or\n"
            "  pip install -e <repo>/utils/green_function"
        ) from None


HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parents[1]          # .../green_function_database/global
DEFAULT_DB = EXAMPLE / "GFDB"
DEFAULT_CMT = EXAMPLE / "validation_data" / "CMTSOLUTION"

TIMINGS = []


def rule(title):
    print()
    print(title)
    print("-" * len(title))


def timed(label, fn, *args, **kwargs):
    """Call fn, record how long it took, and report it."""
    t0 = time.perf_counter()
    result = fn(*args, **kwargs)
    ms = (time.perf_counter() - t0) * 1e3
    TIMINGS.append((label, ms))
    return result, ms


def _mpl():
    """matplotlib, headless. Imported here rather than at the top so that
    --no-plot needs no display and no matplotlib at all."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


# ---------------------------------------------------------------------------
# 1-3: open, extract, and the identity
# ---------------------------------------------------------------------------


def describe(db):
    rule("1. the database")
    info = db.info
    print(f"  path              {db.path}")
    print(f"  elements          {info['nelem']}")
    print(f"  stored samples    {info['nt_subsampled']} of {info['nstep']} solver steps")
    print(f"  sample spacing    {info['dt']} s x {info['subsample_step']} "
          f"= {info['dt'] * info['subsample_step']:.4g} s")
    print(f"  topography        {info['topography']}   ellipticity {info['ellipticity']}")
    print(f"  stations          {info['nstations']}")
    for s in db.stations:
        print(f"    {s.id:<10s} {s.latitude:9.4f} {s.longitude:10.4f}   "
              f"burial {s.depth_m:6.1f} m")


def extract(db, cmt):
    rule("2. seismograms and partial derivatives")

    # locating is separated out only to time it: db.partials() would do it
    # anyway. The first locate in a process is the expensive one -- it builds
    # the element search tree, and loads the topography grid -- and every
    # later one reuses both.
    loc, t_locate = timed("db.locate(...)  first call", db.locate,
                          cmt.latitude, cmt.longitude, cmt.depth)
    _, t_locate2 = timed("db.locate(...)  again", db.locate,
                         cmt.latitude, cmt.longitude, cmt.depth)
    r, t_part = timed("db.partials(cmt)", db.partials, cmt)

    print(f"  db.locate         {t_locate:7.1f} ms   (first call: builds the search tree)")
    print(f"                    {t_locate2:7.1f} ms   (second call: the tree is already there)")
    print(f"  db.partials       {t_part:7.1f} ms   element {loc.ielem}, "
          f"{loc.morton_hex}")
    print()
    print(f"  seismograms       {r.data.shape}   (stations, components N/E/Z, samples)")
    print(f"  partials          {r.dp.shape}   (stations, parameters, components, samples)")
    print(f"  time axis         {r.t.shape}, {r.t[0]:.2f} s to {r.t[-1]:.2f} s "
          f"relative to the centroid time")
    print(f"  peak displacement {np.abs(r.data).max():.4e} m")
    print()
    print("  parameter  unit")
    for name, unit in zip(r.dp_names, r.dp_units):
        print(f"    {name:<8s} {unit}")

    if cmt.centroid_time is not None:
        print()
        print(f"  note: t = 0 is the centroid time, {cmt.centroid_time.isoformat()},")
        print(f"        which is the origin time plus the file's {cmt.time_shift:g} s shift.")
        print("        The shift is header metadata; it is not in the trace.")

    return r


def identity(r, cmt):
    rule("3. the moment-tensor partials are exact")

    m = np.asarray(cmt.tensor)
    reconstructed = (r.dp[:, :6] * m[None, :, None, None]).sum(axis=1)
    err = np.abs(reconstructed - r.data).max() / np.abs(r.data).max()

    print("  The first six partials are per dyne-cm and the seismogram is linear")
    print("  in the moment tensor, so contracting them with the CMTSOLUTION's own")
    print("  components has to give the seismogram back:")
    print()
    print("    (dp[:, :6] * cmt.tensor).sum(1)  vs  data")
    print(f"    max difference {err:.3e} of the trace peak")
    print()
    print("  So a moment-tensor perturbation needs no re-extraction at all: the")
    print("  six partials *are* the forward operator. That is what makes a linear")
    print("  moment-tensor inversion one matrix solve over these arrays.")

    return err


# ---------------------------------------------------------------------------
# 4: the centroid partials are derivatives
# ---------------------------------------------------------------------------


def frechet(db, cmt, r, step, nhalve=2):
    """Predict a relocated seismogram from dp/dlat, and halve the step."""
    rule("4. the centroid partials are derivatives")

    print("  Relocating the source changes the seismogram non-linearly, so here")
    print("  the partial is a genuine derivative and one Taylor term is only an")
    print("  approximation:")
    print()
    print("    u(lat + h)  ~=  u(lat) + h * du/dlat")
    print()
    print("  The error of that must fall as h^2. Halving h therefore divides it")
    print("  by four -- and a partial with the wrong sign, units or scaling")
    print("  cannot produce that, however plausible its traces look.")
    print()

    ilat = r.dp_names.index("lat")
    peak = np.abs(r.data).max()

    steps, errors, truths, predictions = [], [], [], []
    for k in range(nhalve + 1):
        h = step / (2 ** k)

        moved = copy.deepcopy(cmt)
        moved.latitude = cmt.latitude + h

        try:
            truth, _ = timed(f"db.seismograms(...)  relocated {h:g} deg",
                             lambda s: db.seismograms(s).data, moved)
        except gf3d.GF3DError as exc:
            print(f"    step {h:<9g} could not be extracted: {exc.name}")
            print("    (the relocated source left the database; try a smaller --step)")
            break

        predicted = r.data + r.dp[:, ilat] * h
        err = np.abs(predicted - truth).max() / peak

        steps.append(h)
        errors.append(err)
        truths.append(truth)
        predictions.append(predicted)

    if not steps:
        print("  nothing to report: no relocation could be extracted.")
        return steps, errors, truths, predictions

    print(f"    {'step [deg]':<12s} {'error / peak':<14s} {'ratio':<8s}")
    for k, (h, err) in enumerate(zip(steps, errors)):
        ratio = f"{errors[k - 1] / err:.2f}" if k else ""
        print(f"    {h:<12g} {err:<14.3e} {ratio:<8s}")

    if len(errors) > 1:
        ratios = [errors[k - 1] / errors[k] for k in range(1, len(errors))]
        mean = sum(ratios) / len(ratios)
        verdict = "second order, as a first derivative must be" if 3.0 < mean < 5.0 \
            else "NOT second order -- something is wrong"
        print()
        print(f"  mean ratio {mean:.2f}: {verdict}")

    return steps, errors, truths, predictions


# ---------------------------------------------------------------------------
# the figure
# ---------------------------------------------------------------------------


def plot(r, cmt, ista, icomp, steps, errors, truths, predictions, png):
    plt = _mpl()

    sid = r.station_ids[ista]
    comp = "NEZ"[icomp]
    t = r.t

    fig = plt.figure(figsize=(13, 10.5))
    gs = fig.add_gridspec(nrows=4, ncols=1, height_ratios=[1.0, 1.7, 1.0, 0.85],
                          hspace=0.38)

    # -- the seismogram ----------------------------------------------------
    ax = fig.add_subplot(gs[0])
    for j, c in enumerate("NEZ"):
        ax.plot(t, r.data[ista, j], lw=0.9, label=c)
    ax.axvline(0.0, color="0.5", lw=0.5)
    ax.set_xlim(t[0], t[-1])
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=3, loc="upper right")
    ax.set_ylabel("displacement [m]", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title(f"{sid}   db.partials(cmt).data[{ista}]   "
                 f"(t = 0 is the centroid time)", fontsize=10, loc="left")

    # -- the ten Frechet kernels -------------------------------------------
    # normalised individually: they span thirty orders of magnitude between
    # m/dyne-cm and m/s, so only the shapes can share an axis
    ax = fig.add_subplot(gs[1])
    ndp = r.dp.shape[1]
    peaks = []
    for ip in range(ndp):
        trace = r.dp[ista, ip, icomp]
        scale = np.abs(trace).max()
        peaks.append(scale)
        ax.plot(t, 0.42 * trace / scale + ip, "k-", lw=0.7)
    ax.axvline(0.0, color="0.5", lw=0.5)
    ax.set_xlim(t[0], t[-1])
    ax.set_ylim(-0.7, ndp - 0.3)
    ax.set_yticks(range(ndp))
    ax.set_yticklabels(
        [f"{n}   {p:.1e} {u}" for n, p, u in zip(r.dp_names, peaks, r.dp_units)],
        fontsize=6.5, fontfamily="monospace",
    )
    ax.grid(True, axis="x", alpha=0.3)
    ax.tick_params(labelsize=7)
    ax.set_title(f"the ten partial derivatives, {sid} {comp}, "
                 "each normalised by the peak given with its name",
                 fontsize=10, loc="left")

    # -- the Frechet check, the overlay ------------------------------------
    ax = fig.add_subplot(gs[2])
    if steps:
        h = steps[0]
        ax.plot(t, truths[0][ista, icomp], "k-", lw=1.0,
                label=f"re-extracted at latitude + {h:g}°")
        ax.plot(t, predictions[0][ista, icomp], "r--", lw=1.0,
                label="predicted from the unperturbed run, one Taylor term")
    ax.axvline(0.0, color="0.5", lw=0.5)
    ax.set_xlim(t[0], t[-1])
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc="upper right")
    ax.set_ylabel("displacement [m]", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title("a relocated source, predicted rather than re-extracted",
                 fontsize=10, loc="left")

    # -- and the residuals, on an axis where they are visible ---------------
    # the whole point of the panel: on the scale above, a 0.3 % residual is
    # a line thickness
    ax = fig.add_subplot(gs[3])
    for k, hk in enumerate(steps):
        resid = predictions[k][ista, icomp] - truths[k][ista, icomp]
        ax.plot(t, resid, lw=0.8,
                label=f"{hk:g}°   max {errors[k]:.2e} of the trace peak")
    # The last khalf samples are where the source time function kernel runs
    # off the end of the stored record and convolves zeros, so everything is
    # a little wrong there -- for the same reason in every trace, which is
    # why the ratios above are unaffected. The plan says how many.
    if r.plan.khalf > 0:
        edge = t[-1] - r.plan.khalf * r.plan.dt_sub
        ax.axvspan(edge, t[-1], color="0.88", zorder=0)
        ax.text(edge, ax.get_ylim()[0],
                f" last {r.plan.khalf} samples: the kernel runs off the record",
                fontsize=6.5, va="bottom", ha="left", color="0.35")
    ax.axhline(0.0, color="0.6", lw=0.6, ls=":")
    ax.axvline(0.0, color="0.5", lw=0.5)
    ax.set_xlim(t[0], t[-1])
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, ncol=3, loc="lower left", title="relocation step",
              title_fontsize=7)
    ax.set_xlabel("time relative to the centroid time [s]", fontsize=8)
    ax.set_ylabel("prediction − truth [m]", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title("the residual, which falls fourfold each time the step is halved",
                 fontsize=10, loc="left")

    fig.suptitle(
        f"gf3d Python API — event {cmt.event_name or '(unnamed)'} at "
        f"{cmt.latitude:.3f}°, {cmt.longitude:.3f}°, {cmt.depth:.1f} km",
        fontsize=11,
    )

    fig.savefig(png, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {png}")


# ---------------------------------------------------------------------------


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--db", default=str(DEFAULT_DB), help="the GFDB directory")
    ap.add_argument("--cmt", default=str(DEFAULT_CMT), help="a CMTSOLUTION")
    ap.add_argument("--station", default=None, help="NET.STA to draw (default: the first)")
    ap.add_argument("--component", default="Z", choices=list("NEZ"),
                    help="component to draw (default: Z)")
    ap.add_argument("--step", type=float, default=0.05,
                    help="relocation for step 4, in degrees of latitude (default: 0.05)")
    ap.add_argument("--png", default=str(HERE / "gf3d_api_demo.png"), help="figure to write")
    ap.add_argument("--no-plot", action="store_true", help="numbers only")
    args = ap.parse_args(argv)

    db_path = Path(args.db)
    if not (db_path / "mesh_info.h5").exists():
        raise SystemExit(
            f"no database at {db_path}\n"
            "The example's database is built by the workflow, and is gitignored:\n"
            f"    cd {EXAMPLE}\n"
            "    snakemake -j1\n"
            "Or point --db at another GFDB directory."
        )

    cmt_path = Path(args.cmt)
    if not cmt_path.exists():
        raise SystemExit(f"no CMTSOLUTION at {cmt_path}")

    print("gf3d Python API demonstration")
    print(f"  package  {gf3d.__version__}   library {gf3d.library_version()}")
    print(f"  from     {gf3d.library_path}")

    cmt = gf3d.CMTSource.read(cmt_path)

    db, t_open = timed("gf3d.Database(path)", gf3d.Database, db_path)
    with db:
        print(f"\n  opening the database took {t_open:.1f} ms "
              "(metadata only; elements are read on demand)")
        describe(db)
        r = extract(db, cmt)
        identity(r, cmt)
        steps, errors, truths, predictions = frechet(db, cmt, r, args.step)
        if not steps:
            raise SystemExit(
                f"\nthe relocation of {args.step}° left the database, so step 4 showed "
                "nothing.\nTry a smaller --step (the default, 0.05, works on the shipped "
                "example)."
            )

        if not args.no_plot:
            if args.station and args.station not in r.station_ids:
                raise SystemExit(
                    f"\nno station {args.station!r} in this database.\n"
                    f"It has: {', '.join(r.station_ids)}"
                )
            ista = r.station_ids.index(args.station) if args.station else 0
            icomp = "NEZ".index(args.component)
            plot(r, cmt, ista, icomp, steps, errors, truths, predictions, args.png)

    report_timings(db.info["nstations"], r.plan.nt)

    return 0


def report_timings(nsta, nt):
    rule("timings")

    width = max(len(label) for label, _ in TIMINGS)
    for label, ms in TIMINGS:
        print(f"  {label:<{width}s}  {ms:8.1f} ms")

    total = sum(ms for _, ms in TIMINGS)
    print(f"  {'':<{width}s}  {'-' * 8}")
    print(f"  {'in the library':<{width}s}  {total:8.1f} ms")
    print()
    print(f"  Every extraction above produced {nsta} stations x 3 components x {nt}")
    print("  samples. The database is opened once and stays open; the first locate")
    print("  pays for the element search tree and the topography grid, and every")
    print("  source after it reuses them. That difference -- an open per event")
    print("  against an open per run -- is the whole argument for the in-memory")
    print("  API over a subprocess and a SAC file per iteration.")


if __name__ == "__main__":
    sys.exit(main())
