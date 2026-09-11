#!/usr/bin/env python
"""Seismograms for several events from one open database.

The global example's database was not built for a single earthquake. Its
`db_base/DATA/GF_LOCATIONS` names three hypocentres --

    -35.909  -72.733   35.0 km    2010 Chile earthquake (Maule)
     -5.812  -75.270  122.6 km    2019 Peru earthquake
    -13.0    -67.0     10.0 km    near the centre of the chunk

-- and the reciprocal runs stored the strain around all of them. That is
what a Green function database is for: the expensive simulations are done
per *station*, once, and any source inside the covered volume is then a
cheap look-up.

This script is that use case. It opens the database once and extracts
seismograms for every hypocentre the database declares, timing each one, so
the shape of the cost is visible: a few milliseconds to open, a first
locate that pays for the element search tree and the topography grid, and
then one extraction per event that reads only the elements it needs.

The moment tensor is held fixed at the validation event's, so the three
records differ only by where the source sits. A real catalogue would vary
the mechanism too, and that costs nothing extra: the moment-tensor partials
are exact, so a new mechanism at the same hypocentre is a contraction of
arrays already in memory rather than another extraction. `gf3d_api_demo.py`
shows that identity.

Usage
-----
  cd EXAMPLES/green_function_database/global/api_demos/python
  ../../../.venv/bin/python gf3d_events_demo.py [--component N|E|Z] [--png FILE]
                                             [--no-plot]

Requires numpy, matplotlib (for the figure), and a built `lib/libgf3d.so`.
"""

from __future__ import annotations

import argparse
import copy
import sys
import time
from pathlib import Path

import numpy as np

# See the note in gf3d_api_demo.py: the example's virtual environment does
# not have the package, so fall back to the tree this script lives in.
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
DEFAULT_LOCATIONS = EXAMPLE / "db_base" / "DATA" / "GF_LOCATIONS"


def read_locations(path):
    """The hypocentres a GF_LOCATIONS file declares, with their labels.

    Three columns, latitude longitude depth_km, with `#` comments. The
    comment immediately above a row is taken as that row's name, which is
    how the shipped file labels its events.
    """
    events = []
    label = ""
    for line in Path(path).read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            label = stripped.lstrip("#").strip()
            continue
        fields = stripped.split()
        if len(fields) < 3:
            continue
        events.append((label or f"event {len(events) + 1}",
                       float(fields[0]), float(fields[1]), float(fields[2])))
        label = ""
    return events


def _mpl():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def plot(results, station_ids, comp, png, cmt):
    plt = _mpl()

    icomp = "NEZ".index(comp)
    nsta = len(station_ids)

    fig, axes = plt.subplots(nsta, 1, figsize=(13, 2.1 * nsta + 1.2), sharex=True)
    if nsta == 1:
        axes = [axes]

    colours = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

    for ista, (ax, sid) in enumerate(zip(axes, station_ids)):
        for k, res in enumerate(results):
            r = res["result"]
            ax.plot(r.t, r.data[ista, icomp], lw=0.9, color=colours[k % len(colours)],
                    label=(f"{res['label']}   {res['lat']:.3f}°, {res['lon']:.3f}°, "
                           f"{res['depth']:.1f} km" if ista == 0 else None))
        ax.axvline(0.0, color="0.5", lw=0.5)
        ax.grid(True, alpha=0.3)
        ax.set_ylabel(f"{sid}\n{comp} [m]", fontsize=8)
        ax.tick_params(labelsize=7)
        if ista == 0:
            ax.legend(fontsize=7, loc="upper right")

    axes[-1].set_xlabel("time relative to the centroid time [s]", fontsize=8)
    axes[0].set_xlim(results[0]["result"].t[0], results[0]["result"].t[-1])

    fig.suptitle(
        f"one open database, {len(results)} hypocentres — same moment tensor "
        f"(M0 from {cmt.event_name or 'the validation event'}), same stations",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(png, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {png}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--db", default=str(DEFAULT_DB), help="the GFDB directory")
    ap.add_argument("--cmt", default=str(DEFAULT_CMT),
                    help="a CMTSOLUTION, for the moment tensor and origin time")
    ap.add_argument("--locations", default=str(DEFAULT_LOCATIONS),
                    help="a GF_LOCATIONS file listing the hypocentres")
    ap.add_argument("--component", default="Z", choices=list("NEZ"),
                    help="component to draw (default: Z)")
    ap.add_argument("--png", default=str(HERE / "gf3d_events_demo.png"),
                    help="figure to write")
    ap.add_argument("--no-plot", action="store_true", help="numbers only")
    args = ap.parse_args(argv)

    db_path = Path(args.db)
    if not (db_path / "mesh_info.h5").exists():
        raise SystemExit(
            f"no database at {db_path}\n"
            "The example's database is built by the workflow, and is gitignored:\n"
            f"    cd {EXAMPLE}\n"
            "    snakemake -j1"
        )
    if not Path(args.locations).exists():
        raise SystemExit(f"no GF_LOCATIONS file at {args.locations}")

    events = read_locations(args.locations)
    if not events:
        raise SystemExit(f"{args.locations} declares no hypocentres")

    cmt = gf3d.CMTSource.read(args.cmt)

    print("gf3d: several events, one open database")
    print(f"  library  {gf3d.library_version()} from {gf3d.library_path}")
    print(f"  events   {len(events)}, from {args.locations}")
    print(f"  mechanism held fixed at {cmt.event_name or 'the validation event'}'s, "
          "so only the hypocentre differs")

    t0 = time.perf_counter()
    db = gf3d.Database(db_path)
    t_open = (time.perf_counter() - t0) * 1e3

    with db:
        station_ids = db.station_ids
        print(f"  stations {', '.join(station_ids)}")
        print()
        print(f"  gf3d.Database(...)            {t_open:7.1f} ms")
        print()

        header = (f"  {'event':<28s} {'lat':>9s} {'lon':>10s} {'depth':>7s} "
                  f"{'elem':>5s} {'locate':>9s} {'extract':>9s} {'peak [m]':>11s}")
        print(header)
        print("  " + "-" * (len(header) - 2))

        results = []
        for label, lat, lon, depth in events:
            src = copy.deepcopy(cmt)
            src.latitude, src.longitude, src.depth = lat, lon, depth

            try:
                t0 = time.perf_counter()
                loc = db.locate(lat, lon, depth)
                t_loc = (time.perf_counter() - t0) * 1e3

                t0 = time.perf_counter()
                r = db.seismograms(src)
                t_ext = (time.perf_counter() - t0) * 1e3
            except gf3d.GF3DError as exc:
                print(f"  {label:<28s} {lat:9.3f} {lon:10.3f} {depth:7.1f} "
                      f"   --  {exc.name}")
                continue

            print(f"  {label:<28s} {lat:9.3f} {lon:10.3f} {depth:7.1f} "
                  f"{loc.ielem:5d} {t_loc:7.1f} ms {t_ext:7.1f} ms "
                  f"{np.abs(r.data).max():11.4e}")

            results.append({"label": label, "lat": lat, "lon": lon, "depth": depth,
                            "result": r, "locate_ms": t_loc, "extract_ms": t_ext})

        if not results:
            raise SystemExit("\nnone of the declared hypocentres is inside this database")

        # per station, so the reader can see which event each station favours
        print()
        print("  peak amplitude per station [m]")
        width = max(len(r["label"]) for r in results)
        print(f"    {'event':<{width}s}  " +
              "  ".join(f"{s:>11s}" for s in station_ids))
        for res in results:
            peaks = np.abs(res["result"].data).max(axis=(1, 2))
            print(f"    {res['label']:<{width}s}  " +
                  "  ".join(f"{p:11.4e}" for p in peaks))

        extracts = [r["extract_ms"] for r in results]
        print()
        print(f"  {len(results)} events: {t_open:.1f} ms to open, "
              f"{results[0]['locate_ms']:.1f} ms for the first locate, then "
              f"{np.mean(extracts):.1f} ms per event on average.")
        print("  The reciprocal simulations that built this database took about")
        print("  forty minutes. Every source inside it is now a look-up.")

        if not args.no_plot:
            plot(results, station_ids, args.component, args.png, cmt)

    return 0


if __name__ == "__main__":
    sys.exit(main())
