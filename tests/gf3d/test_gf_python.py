#!/usr/bin/env python
"""test_gf_python -- the ctypes package against xgf3d itself.

Tier 2: needs the shared library and a built example database.

The oracle is the executable. Everything the Python package returns has
come through a chain -- ctypes, the bind(C) facade, the transposing copy --
that the Fortran tests do not exercise, and the way to check that chain end
to end is to ask ``xgf3d --seis --format ascii`` for the same source and
compare columns. Those ASCII files are the same ones the workflow's
comparison harness gates on, so agreement here ties the in-memory API to
the numbers the project is validated by.

The rest is the package's own contract: shapes, station order, the names
and units of the partials, the linearity identity the moment-tensor
partials satisfy, and -- the reason the facade exists -- that every
mistake raises a Python exception instead of ending the interpreter.

Usage: test_gf_python.py <xgf3d> <GFDB> <CMTSOLUTION> [FORCESOLUTION] [second GFDB]
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import TimeoutError as FuturesTimeout
from pathlib import Path

import numpy as np

import gf3d

TOL = 1.0e-12

nfail = 0


def ok(name, cond):
    global nfail
    if cond:
        print(f"  ok   {name}")
    else:
        print(f"  FAIL {name}")
        nfail += 1


def ok_err(name, err, tol=TOL):
    global nfail
    if err <= tol and not np.isnan(err):
        print(f"  ok   {name:<40s} error = {err:12.5e}  (tol {tol:12.5e})")
    else:
        print(f"  FAIL {name:<40s} error = {err:12.5e}  (tol {tol:12.5e})")
        nfail += 1


def ok_raises(name, code, fn, *args, **kwargs):
    """The call must raise GF3DError with this status, and leave us running."""
    global nfail
    try:
        fn(*args, **kwargs)
    except gf3d.GF3DError as exc:
        if code is None or exc.code == code:
            print(f"  ok   {name:<40s} -> {exc.name}")
            return
        print(f"  FAIL {name:<40s} -> {exc.name} (wanted code {code})")
        nfail += 1
        return
    except Exception as exc:  # noqa: BLE001
        print(f"  FAIL {name:<40s} -> {type(exc).__name__}: {exc}")
        nfail += 1
        return
    print(f"  FAIL {name:<40s} -> no exception")
    nfail += 1


def read_ascii(path):
    """One NET.STA.gf3d.txt: the time column and the three components."""
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 1:4]


def main(argv):
    if len(argv) < 4:
        print(__doc__)
        return 2

    xgf3d = Path(argv[1]).resolve()
    dbpath = Path(argv[2]).resolve()
    cmtpath = Path(argv[3]).resolve()
    forcepath = Path(argv[4]).resolve() if len(argv) > 4 and argv[4] else None
    dbpath2 = Path(argv[5]).resolve() if len(argv) > 5 and argv[5] else None

    print()
    print(" ******************************")
    print(" test_gf_python")
    print(" ******************************")
    print()

    # ------------------------------------------------------------------
    print(" 1. the package and the library it found")
    print(f"       library: {gf3d.library_path}")
    print(f"       version: {gf3d.library_version()}")
    ok("a version string came back", len(gf3d.library_version()) > 0)
    # the struct layout check runs at import; getting here means it passed
    ok("the struct layouts agree with the header", True)

    # ------------------------------------------------------------------
    print("\n 2. opening")

    ok_raises("a database that is not there", gf3d.GF_ERR_NO_PATH,
              gf3d.Database, "/nonexistent/gf3d/database")
    try:
        gf3d.Database("/nonexistent/gf3d/database")
    except gf3d.GF3DError as exc:
        ok("the message names the path", "/nonexistent/gf3d/database" in str(exc))

    # Stations before info, on a database that has cached neither.
    #
    # These properties compose -- stations needs info -- and both call the
    # library under its lock, so with a non-reentrant lock this deadlocks
    # unless something happened to populate the cache first. Every other
    # caller here reads info first and hid it. Do it in this order, once,
    # and with a watchdog, because the symptom of a regression is a hang
    # rather than a failure and a hung test tells nobody anything.
    def _first_call_order():
        fresh = gf3d.Database(dbpath)
        try:
            return fresh.station_ids
        finally:
            fresh.close()

    with ThreadPoolExecutor(max_workers=1) as pool:
        future = pool.submit(_first_call_order)
        try:
            ids_first = future.result(timeout=60)
            ok("station_ids works before info is cached", len(ids_first) > 0)
        except FuturesTimeout:
            ok("station_ids works before info is cached (DEADLOCK)", False)
            return 1

    db = gf3d.Database(dbpath)
    ok("the database is open", not db.closed)
    print(f"       {db!r}")

    info = db.info
    ok("info is complete", info["nelem"] > 0 and info["nstations"] > 0 and info["dt"] > 0)
    ok("r_planet and rhoav are physical", info["r_planet"] > 1e6 and info["rhoav"] > 100)

    ids = db.station_ids
    on_disk = sorted(p.stem for p in (dbpath / "stations").glob("*.h5"))
    ok("station_ids are the station files, in order", ids == on_disk)
    print(f"       stations: {ids}")

    stations = db.stations
    ok("every station carries a network and a name",
       all(s.network and s.station and s.id == f"{s.network}.{s.station}" for s in stations))

    # ------------------------------------------------------------------
    print("\n 3. the source, and where it sits")

    cmt = gf3d.CMTSource.read(cmtpath)
    print(f"       {cmt.latitude}, {cmt.longitude} at {cmt.depth} km, "
          f"hdur {cmt.hdur} s, shift {cmt.time_shift} s, event {cmt.event_name}")
    ok("the moment tensor is six numbers in dyne-cm",
       len(cmt.tensor) == 6 and abs(cmt.Mrr) > 1e20)
    ok("the origin time was parsed", cmt.origin_time is not None)
    ok("the centroid time is the origin plus the shift",
       cmt.centroid_time is not None
       and abs((cmt.centroid_time - cmt.origin_time).total_seconds() - cmt.time_shift) < 1e-6)

    loc = db.locate(cmt.latitude, cmt.longitude, cmt.depth)
    print(f"       element {loc.ielem} ({loc.morton_hex}) "
          f"at {loc.xi:.6f} {loc.eta:.6f} {loc.gamma:.6f}")
    ok("the point was found where it was asked for", loc.distance_km < 1e-6)

    ok_raises("a NaN latitude", gf3d.GF_ERR_ARG, db.locate, float("nan"), 0.0, 10.0)
    ok_raises("the far side of the planet", gf3d.GF_ERR_NO_ELEMENT,
              db.locate, -cmt.latitude, cmt.longitude + 180.0, 10.0)

    plan = db.plan(cmt)
    print(f"       nt = {plan.nt}, dt_sub = {plan.dt_sub}, t_first = {plan.t_first}")
    ok("the plan is the stored length plus the padding", plan.nt == plan.nt_db + plan.npad)
    ok("a Heaviside conversion for a moment tensor", plan.kind_stf == 2)
    ok("the plan's own axis matches the arithmetic",
       np.allclose(plan.times, plan.t_first + np.arange(plan.nt) * plan.dt_sub, atol=0))

    # ------------------------------------------------------------------
    print("\n 4. seismograms and partials")

    r = db.partials(cmt)
    nsta = len(ids)
    ok("data is (nstations, 3, nt)", r.data.shape == (nsta, 3, plan.nt))
    ok("dp is (nstations, 10, 3, nt)", r.dp.shape == (nsta, gf3d.GF_NDP_LOC, 3, plan.nt))
    ok("t is (nt,)", r.t.shape == (plan.nt,))
    ok("everything is finite", np.isfinite(r.data).all() and np.isfinite(r.dp).all())
    ok("the traces are not all zero", np.abs(r.data).max() > 0)

    ok("the partial names are GF3DF's order",
       r.dp_names == ["Mrr", "Mtt", "Mpp", "Mrt", "Mrp", "Mtp", "lat", "lon", "dep", "tim"])
    ok("the units are per dyne-cm, degree, km and second",
       r.dp_units[0] == "m/dyne-cm" and r.dp_units[6] == "m/deg"
       and r.dp_units[8] == "m/km" and r.dp_units[9] == "m/s")

    # the identity of Stage 6, which no transposition of a four-dimensional
    # index survives
    lin = (r.dp[:, :6] * np.asarray(cmt.tensor)[None, :, None, None]).sum(axis=1)
    ok_err("sum(M_v dp_v) reproduces the seismogram",
           np.abs(lin - r.data).max() / np.abs(r.data).max())

    ok("trace() and partial() select the same data",
       np.array_equal(r.trace(ids[0], "Z"), r.data[0, 2])
       and np.array_equal(r.partial(ids[0], "N", "dep"), r.dp[0, 8, 0]))

    r1 = db.partials(cmt, kind=1)
    ok("kind=1 gives six partials", r1.dp.shape[1] == gf3d.GF_NDP_MT)
    ok_err("kind=1 and kind=2 agree on the six",
           np.abs(r1.dp - r.dp[:, :6]).max() / np.abs(r.dp[:, :6]).max())
    ok_err("and on the seismogram",
           np.abs(r1.data - r.data).max() / np.abs(r.data).max())

    try:
        db.partials(cmt, kind=3)
        ok("kind=3 refused", False)
    except ValueError:
        ok("kind=3 refused", True)

    # ------------------------------------------------------------------
    print("\n 5. against xgf3d --seis --format ascii")

    with tempfile.TemporaryDirectory() as tmp:
        cmd = [str(xgf3d), "--seis", str(dbpath), str(cmtpath), tmp, "--format", "ascii"]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        ok("xgf3d ran", proc.returncode == 0)
        if proc.returncode != 0:
            print(proc.stdout[-2000:])
            print(proc.stderr[-2000:])
        else:
            worst_t = 0.0
            worst_d = 0.0
            for ista, sid in enumerate(ids):
                t_ref, d_ref = read_ascii(Path(tmp) / f"{sid}.gf3d.txt")
                worst_t = max(worst_t, np.abs(t_ref - r.t).max())
                # the ASCII carries (nt, 3) as N, E, Z; the array is (3, nt)
                diff = np.abs(d_ref.T - r.data[ista])
                worst_d = max(worst_d, diff.max() / max(np.abs(d_ref).max(), 1e-300))
            ok_err("the time axis matches the executable", worst_t)
            ok_err("every station's trace matches the executable", worst_d)

    # ------------------------------------------------------------------
    print("\n 6. a force source")

    if forcepath is not None and forcepath.exists():
        force = gf3d.ForceSource.read(forcepath)
        print(f"       f0 {force.f0}, stf {force.stf}, factor {force.factor:g}, "
              f"direction {force.direction}")
        fr = db.seismograms(force)
        ok("force seismograms are (nstations, 3, nt)",
           fr.data.shape[0] == nsta and fr.data.shape[1] == 3)
        ok("they are finite and not all zero",
           np.isfinite(fr.data).all() and np.abs(fr.data).max() > 0)
        ok("a force source has no partials", fr.dp is None)
        ok_raises("partials of a force source", gf3d.GF_ERR_ARG, db.partials, force)
    else:
        print("       (no FORCESOLUTION in this example, skipped)")

    # ------------------------------------------------------------------
    print("\n 7. obspy, if it is installed")

    try:
        import obspy  # noqa: F401
    except ImportError:
        print("       (obspy not installed, skipped)")
    else:
        st = r.to_stream()
        ok("one trace per station and component", len(st) == 3 * nsta)
        tr = st[2]
        ok("the channel follows the solver's naming", tr.stats.channel == "BXZ")
        ok("the sample spacing is the stored one", abs(tr.stats.delta - plan.dt_sub) < 1e-6)
        ok("the data survived the trip", np.array_equal(tr.data, r.data[0, 2]))
        start = tr.stats.starttime
        expected = obspy.UTCDateTime(cmt.centroid_time) + float(r.t[0])
        ok("the start time is the centroid time plus t[0]", abs(start - expected) < 1e-6)

    # ------------------------------------------------------------------
    print("\n 8. closing, and a second database")

    if dbpath2 is not None and (dbpath2 / "mesh_info.h5").exists():
        # The two-database case: opening the second re-installs specfem's
        # process-wide globals, and the facade re-installs the first
        # handle's before every call. If it did not, this extraction would
        # come back subtly different -- the topography grid indexed with the
        # other database's dimensions.
        db2 = gf3d.Database(dbpath2)
        print(f"       second: {db2!r}")
        r2 = db2.seismograms(cmt)
        ok("the second database extracts", np.isfinite(r2.data).all())

        again = db.partials(cmt)
        ok("the first database is bit-for-bit unchanged",
           np.array_equal(again.data, r.data) and np.array_equal(again.dp, r.dp))
        db2.close()
    else:
        print("       (only one example database is built, skipped)")
        again = db.partials(cmt)
        ok("a repeated extraction is bit-for-bit identical",
           np.array_equal(again.data, r.data) and np.array_equal(again.dp, r.dp))

    db.close()
    ok("the database is closed", db.closed)
    db.close()
    ok("closing twice is harmless", db.closed)
    ok_raises("using a closed database", gf3d.GF_ERR_ARG, db.locate, 0.0, 0.0, 10.0)

    with gf3d.Database(dbpath) as ctx:
        ok("the context manager opens", not ctx.closed)
    ok("and closes on the way out", ctx.closed)

    print()
    if nfail == 0:
        print(" test_gf_python: all assertions passed\n")
        return 0
    print(f" test_gf_python: {nfail} assertion(s) FAILED\n")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
