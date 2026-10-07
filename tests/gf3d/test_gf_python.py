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

Usage: test_gf_python.py <xgf3d> <GFDB> <CMTSOLUTION> [FORCESOLUTION [block]]

With "block" as the fifth argument only sections 10 and 11, element_block and
the weights, are run. The runner does that on the fixture layout the full run did not use (the
contiguous fixture reads through the h5dread_f fallback, the chunked one
through the raw H5Dread_chunk route).
"""

from __future__ import annotations

import dataclasses
import math
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import TimeoutError as FuturesTimeout
from pathlib import Path

import numpy as np

import gf3d

TOL = 1.0e-12
# the weights against the library, relative to the sum of |terms| of a trace
TOLW = 1.0e-13

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


def ok_raises_py(name, exc_type, fn, *args, **kwargs):
    """The call must raise this Python exception, before reaching the library."""
    global nfail
    try:
        fn(*args, **kwargs)
    except exc_type as exc:
        print(f"  ok   {name:<40s} -> {type(exc).__name__}")
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


def element_block_section(dbpath, cmt):
    """Database.element_block: against the station files, itself, and the cache.

    The oracle for the values is h5py reading the station files directly:
    the same stored numbers through a second reader, with no arithmetic in
    between, so equality is bitwise and the reference is never recomputed.
    Everything else compares a selection or a prefix with the full read.
    """
    try:
        import h5py
    except ImportError:
        h5py = None

    # two elements: the extraction below leaves one in the cache
    db = gf3d.Database(dbpath, max_elements=2)
    info = db.info
    ids = db.station_ids
    nsta = len(ids)
    loc = db.locate(cmt.latitude, cmt.longitude, cmt.depth)
    ielem = loc.ielem
    nt_all = info["nt_subsampled"]

    # the whole record of every station is 3 GB on a real database: a prefix
    # of at most 256 MB there, everything on the fixture
    ntf = max(1, min(nt_all, int(256e6 // (nsta * 3 * 375 * 4))))
    if ntf < nt_all:
        print(f"       a large database: the full read is the first {ntf} of {nt_all} samples")
    nt7 = min(7, ntf)
    # the fixture's first stored samples are all zero, so a comparison over a
    # short prefix compares zeros: ntp ends past the onset, inside a chunk
    # (the chunked fixture's are 7 samples), and its last sample is checked
    # to be nonzero below
    ntp = max(1, ntf // 2 + 1)

    r = db.seismograms(cmt)
    ok("an extraction fills the cache", db.cache_stats["n_cached"] == 1)

    full = db.element_block(ielem, nt=ntf)
    ok("the block is float32, C-contiguous, (nstations, 3, nt, 375)",
       full.dtype == np.float32 and full.flags["C_CONTIGUOUS"]
       and full.shape == (nsta, 3, ntf, 375))
    ok("every value is finite and not all zero", np.isfinite(full).all() and np.any(full != 0))
    ok(f"sample {ntp} of {ntf}, the last of the test prefix, is not all zero",
       np.any(full[:, :, ntp - 1, :] != 0))
    ok("nt=None reads every stored sample",
       ntf < nt_all or db.element_block(ielem).shape == (nsta, 3, nt_all, 375))

    # ------------------------------------------------------------------
    # against the station files
    if h5py is None:
        print("       skipped: h5py not installed")
    else:
        # the library's own index is ascending Morton code; the directory
        # scan is an independent way to the same order
        dirs = sorted(p.name for p in (dbpath / "elements").iterdir() if p.is_dir())
        ok("the element directories are the index, in order",
           len(dirs) == info["nelem"] and dirs[ielem - 1] == loc.morton_hex)

        def reference(el, nt):
            """[s, a, t, m] from the files: h5py's (t, k, j, i, p, a)."""
            out = np.empty((nsta, 3, nt, 375), dtype=np.float32)
            for s, sid in enumerate(ids):
                with h5py.File(dbpath / "elements" / dirs[el - 1] / f"{sid}.h5", "r") as f:
                    d = f["displacement"][:nt]
                out[s] = np.transpose(d, (5, 0, 1, 2, 3, 4)).reshape(3, nt, 375)
            return out

        ok("the located element, every station, bitwise equal to the files",
           np.array_equal(full, reference(ielem, ntf)))

        # every element when the database is small, else the first few
        others = [e for e in range(1, info["nelem"] + 1) if e != ielem][:7]
        same = True
        for e in others:
            same = same and np.array_equal(db.element_block(e, nt=ntp), reference(e, ntp))
        ok(f"{len(others)} other elements, a prefix, bitwise equal to the files", same)

    # ------------------------------------------------------------------
    # prefixes and refusals
    print()
    same = True
    for n in sorted({1, nt7, ntp, ntf}):
        same = same and np.array_equal(db.element_block(ielem, nt=n), full[:, :, :n, :])
    ok(f"nt = 1, 7, {ntp} and the full length are prefixes of the full read, bitwise", same)

    ok_raises("nt = nt_subsampled + 1", gf3d.GF_ERR_ARG, db.element_block, ielem, nt=nt_all + 1)
    ok_raises("nt = 0", gf3d.GF_ERR_ARG, db.element_block, ielem, nt=0)
    ok_raises("nt = -1", gf3d.GF_ERR_ARG, db.element_block, ielem, nt=-1)
    ok_raises("ielem = 0", gf3d.GF_ERR_ARG, db.element_block, 0)
    ok_raises("ielem = nelem + 1", gf3d.GF_ERR_ARG, db.element_block, info["nelem"] + 1)
    ok_raises("ielem = 2**40", gf3d.GF_ERR_ARG, db.element_block, 2**40)

    # ------------------------------------------------------------------
    # stations
    if nsta >= 2:
        want = full[[1, 0]]
        ok("stations=[1, 0] is full[[1, 0]]",
           np.array_equal(db.element_block(ielem, stations=[1, 0], nt=ntf), want))
        ok("stations by id, reversed, is the same",
           np.array_equal(db.element_block(ielem, stations=[ids[1], ids[0]], nt=ntf), want))
        ok("stations mixed, id and index, is the same",
           np.array_equal(db.element_block(ielem, stations=[ids[1], 0], nt=ntf), want))
        ok("a repeated station is read twice",
           np.array_equal(db.element_block(ielem, stations=[0, 0], nt=ntp),
                          full[[0, 0], :, :ntp, :]))
        ok("a tuple and a numpy array select the same",
           np.array_equal(db.element_block(ielem, stations=(1, 0), nt=ntp),
                          db.element_block(ielem, stations=np.array([1, 0]), nt=ntp)))
    ok_raises_py("an unknown station id", ValueError,
                 db.element_block, ielem, stations=["XX.NOPE"])
    try:
        db.element_block(ielem, stations=["XX.NOPE"])
    except ValueError as exc:
        ok("the message names the id", "XX.NOPE" in str(exc))
    ok_raises(f"station index {nsta}", gf3d.GF_ERR_ARG,
              db.element_block, ielem, stations=[nsta])
    ok_raises("station index -1", gf3d.GF_ERR_ARG, db.element_block, ielem, stations=[-1])
    ok_raises("station index 2**40", gf3d.GF_ERR_ARG,
              db.element_block, ielem, stations=[2**40])
    ok_raises("no stations", gf3d.GF_ERR_ARG, db.element_block, ielem, stations=[])

    # ------------------------------------------------------------------
    # out=
    print()
    buf = np.empty((nsta, 3, ntp, 375), dtype=np.float32)
    got = db.element_block(ielem, nt=ntp, out=buf)
    ok("out= is returned itself", got is buf)
    ok("and filled", np.array_equal(buf, full[:, :, :ntp, :]))
    ok_raises_py("out= of the wrong shape", ValueError, db.element_block, ielem, nt=nt7,
                 out=np.empty((nsta, 3, nt7 + 1, 375), dtype=np.float32))
    ok_raises_py("out= of float64", ValueError, db.element_block, ielem, nt=nt7,
                 out=np.empty((nsta, 3, nt7, 375), dtype=np.float64))
    ok_raises_py("out= not contiguous (a slice)", ValueError, db.element_block, ielem, nt=nt7,
                 out=np.empty((nsta, 3, nt7, 376), dtype=np.float32)[..., :375])
    ok_raises_py("out= not contiguous (a transposed view)", ValueError,
                 db.element_block, ielem, nt=nt7,
                 out=np.empty((375, nt7, 3, nsta), dtype=np.float32).transpose(3, 2, 1, 0))
    ok_raises_py("out= not an ndarray", ValueError, db.element_block, ielem, nt=nt7,
                 out=buf.tolist())

    # ------------------------------------------------------------------
    # the cache is neither used nor filled, and files_read counts each file
    print()
    st0 = db.cache_stats
    db.element_block(ielem, stations=[0, 1] if nsta >= 2 else [0], nt=nt7)
    st1 = db.cache_stats
    nread = 2 if nsta >= 2 else 1
    ok("hits, misses, evictions and n_cached did not move",
       all(st1[k] == st0[k] for k in ("hits", "misses", "evictions", "n_cached")))
    ok(f"files_read grew by {nread}, one per station",
       st1["files_read"] - st0["files_read"] == nread)
    db.element_block(ielem, stations=[0], nt=nt7)
    st2 = db.cache_stats
    ok("one station: files_read grew by 1", st2["files_read"] - st1["files_read"] == 1)
    ok("and the cache still holds its element and counts",
       all(st2[k] == st0[k] for k in ("hits", "misses", "evictions", "n_cached")))
    try:
        db.element_block(0)
    except gf3d.GF3DError:
        pass
    ok("a refused call read no file", db.cache_stats["files_read"] == st2["files_read"])

    db.close()
    ok_raises("a closed database", gf3d.GF_ERR_ARG, db.element_block, ielem)


def stf_conv(x, k, kind, npad, dt_sub):
    """gf3d.h's conversion, vectorised, in float64: x (..., nt_db) -> (..., npad + nt_db)."""
    khalf = (len(k) - 1) // 2
    nd = x.shape[-1]
    n = npad + nd
    xpad = np.zeros(x.shape[:-1] + (n,))
    xpad[..., npad:] = x
    if kind == gf3d.GF_STF_NONE:
        return xpad
    win = np.zeros_like(xpad)
    for j in range(-khalf, khalf + 1):
        if abs(j) >= n:
            continue
        if j >= 0:
            win[..., j:] += k[j + khalf] * xpad[..., : n - j]
        else:
            win[..., : n + j] += k[j + khalf] * xpad[..., -j:]
    if kind == gf3d.GF_STF_GAUSS:
        return win
    cs = np.cumsum(xpad, axis=-1)
    shifted = np.zeros_like(xpad)
    if n > khalf + 1:
        shifted[..., khalf + 1:] = cs[..., : n - khalf - 1]
    return dt_sub * (shifted + win)


def trace_errors(y, yabs, ref):
    """Per trace: max|y - ref| relative to max of the sum of |terms|, and to the peak."""
    d = np.abs(y - ref).max(axis=-1)
    den = yabs.max(axis=-1)
    peak = np.abs(ref).max(axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        e_abs = np.where(den > 0, d / den, np.where(d > 0, np.inf, 0.0))
        e_peak = np.where(peak > 0, d / peak, np.where(d > 0, np.inf, 0.0))
    return e_abs.max(), e_peak.max(), den.max()


def weights_identity(db, src, kind, label):
    """STF(scale * block @ w) is the library's trace, for w and each dw column.

    Through public calls only, on a float64 contraction of the float32 block,
    so the error is rounding in the library's own contraction, not a layout.
    """
    nsta = db.info["nstations"]
    W = db.weights(src, kind)
    p = db.plan(src)
    r = db.seismograms(src) if kind == 0 else db.partials(src, kind=kind)
    nd = p.nt_db
    # the whole record of every station is 3 GB on a real database: the
    # stations that fit in 16 MB of float32, at least one
    nsel = max(1, min(nsta, int(16e6 // (3 * nd * 375 * 4))))
    sel = list(range(nsel))
    blk = db.element_block(W.location.ielem, stations=sel).astype(np.float64)
    k = gf3d.stf_kernel(p.kind_stf, p.hdur_corr, p.dt_sub, p.trunc)
    ok(f"{label}: the kernel has the plan's length", len(k) == 2 * p.khalf + 1)

    sc = W.scale[sel][:, None, None]
    cols = [(None, W.w, r.data[sel])]
    for c in range(W.dw.shape[0]):
        cols.append((c, W.dw[c], r.dp[sel, c]))
    for c, wv, ref in cols:
        x = sc * (blk @ wv)
        xa = np.abs(sc) * (np.abs(blk) @ np.abs(wv))
        y = stf_conv(x, k, p.kind_stf, p.npad, p.dt_sub)
        ya = stf_conv(xa, k, p.kind_stf, p.npad, p.dt_sub)
        e_abs, e_peak, den = trace_errors(y, ya, ref)
        name = "seismogram" if c is None else W.dw_names[c]
        ok(f"{label} {name}: the sum of |terms| is not zero", den > 0)
        ok_err(f"{label} {name}", e_abs, tol=TOLW)
        print(f"       {'':<10s} relative to the trace peak: {e_peak:12.5e}")
    return W, p, r, sel, blk, k


def weights_section(dbpath, cmt, force):
    """Database.weights and stf_kernel: the contraction identity, and the contract."""
    db = gf3d.Database(dbpath)
    info = db.info
    nsta = info["nstations"]

    # ------------------------------------------------------------------
    # moment tensor, kind 2: the seismogram and nine columns
    W, p, r, sel, blk, k = weights_identity(db, cmt, 2, "cmt")
    print("       the centroid time: a shift, not a weight; the library's column uses a second kernel")
    ok("w is (375,), dw (9, 375), scale (nstations,), float64",
       W.w.shape == (375,) and W.dw.shape == (9, 375) and W.scale.shape == (nsta,)
       and W.w.dtype == W.dw.dtype == W.scale.dtype == np.float64)
    ok("the partial names are the first nine",
       list(W.dw_names) == db.partial_names(10)[0][:9]
       and list(W.dw_units) == db.partial_names(10)[1][:9])
    ok("the names are tuples", isinstance(W.dw_names, tuple) and isinstance(W.dw_units, tuple))

    # the located element and the block are the library's own
    ok("the location is db.locate's, exactly", (W.location.ielem, W.location.xi,
       W.location.eta, W.location.gamma) == tuple(
           getattr(db.locate(cmt.latitude, cmt.longitude, cmt.depth), a)
           for a in ("ielem", "xi", "eta", "gamma")))

    # a time prefix gives the whole record's conversion on all but the last
    # khalf samples. The fixture's own kernel is wider than its record, so the
    # prefix is tried on the widest source of a shorter half duration whose
    # kernel leaves room for it.
    nd = p.nt_db
    ntp = nd // 2 + 1
    cs = cmt
    for h in cmt.hdur * np.geomspace(1.0, 0.02, 60):
        cs = dataclasses.replace(cmt, hdur=float(h))
        ps = db.plan(cs)
        if ps.npad + ntp - ps.khalf >= (ps.npad + ntp) // 2:
            break
    ps = db.plan(cs)
    nuse = ps.npad + ntp - ps.khalf
    print(f"       prefix: hdur {cs.hdur:.4g}, nt_db {nd}, npad {ps.npad}, khalf {ps.khalf}, "
          f"prefix {ntp}, {nuse} samples compared")
    ok("the prefix reaches past the kernel's half width, which is not zero",
       nuse >= (ps.npad + ntp) // 2 and ps.khalf >= 1)
    Ws = db.weights(cs, 0)
    rs = db.seismograms(cs)
    ks = gf3d.stf_kernel(ps.kind_stf, ps.hdur_corr, ps.dt_sub, ps.trunc)
    bp = db.element_block(Ws.location.ielem, stations=sel, nt=ntp).astype(np.float64)
    xp = Ws.scale[sel][:, None, None] * (bp @ Ws.w)
    xpa = np.abs(Ws.scale[sel])[:, None, None] * (np.abs(bp) @ np.abs(Ws.w))
    yp = stf_conv(xp, ks, ps.kind_stf, ps.npad, ps.dt_sub)[..., :nuse]
    ypa = stf_conv(xpa, ks, ps.kind_stf, ps.npad, ps.dt_sub)[..., :nuse]
    e_abs, e_peak, den = trace_errors(yp, ypa, rs.data[sel][..., :nuse])
    ok_err("prefix conversion, the whole record's samples", e_abs, tol=TOLW)
    print(f"       {'':<10s} relative to the trace peak: {e_peak:12.5e}")
    ok(f"sample {nuse - 1}, the prefix's last, is not all zero", np.any(yp[..., nuse - 1] != 0))

    # kind 1 and kind 0 are the leading part of kind 2
    W1 = db.weights(cmt, 1)
    W0 = db.weights(cmt, 0)
    ok("kind 1 is the first six columns of kind 2",
       W1.dw.shape == (6, 375) and np.array_equal(W1.dw, W.dw[:6])
       and list(W1.dw_names) == list(W.dw_names[:6]))
    ok("kind 0: w and scale the same, dw is (0, 375), no names",
       np.array_equal(W0.w, W.w) and np.array_equal(W0.scale, W.scale)
       and W0.dw.shape == (0, 375) and W0.dw_names == () and W0.dw_units == ())

    # frozen means frozen
    def assign(a):
        a[0] = 1.0
    ok_raises_py("w is read-only", ValueError, assign, W.w)
    ok_raises_py("dw is read-only", ValueError, assign, W.dw)
    ok_raises_py("scale is read-only", ValueError, assign, W.scale)
    ok_raises_py("the fields are frozen", dataclasses.FrozenInstanceError,
                 setattr, W, "w", W.w)

    ok_raises("kind 3", gf3d.GF_ERR_ARG, db.weights, cmt, 3)
    ok_raises("kind -1", gf3d.GF_ERR_ARG, db.weights, cmt, -1)

    # ------------------------------------------------------------------
    # a force source, kind 0
    if force is not None:
        print()
        Wf = weights_identity(db, force, 0, "force")[0]
        ok("force: dw is (0, 375)", Wf.dw.shape == (0, 375))
        ok("force: the location is db.locate's, exactly", (Wf.location.ielem, Wf.location.xi,
           Wf.location.eta, Wf.location.gamma) == tuple(
               getattr(db.locate(force.latitude, force.longitude, force.depth), a)
               for a in ("ielem", "xi", "eta", "gamma")))
        ok_raises("kind 1 for a force source", gf3d.GF_ERR_ARG, db.weights, force, 1)
        ok_raises("kind 2 for a force source", gf3d.GF_ERR_ARG, db.weights, force, 2)
    else:
        print("       (no FORCESOLUTION in this example, the force source skipped)")

    # ------------------------------------------------------------------
    # stf_kernel
    print()
    pc = p
    kp = gf3d.stf_kernel(pc.kind_stf, pc.hdur_corr, pc.dt_sub, pc.trunc)
    ok("at the plan's parameters: 2*khalf + 1 taps", kp.dtype == np.float64 and kp.shape == (2 * pc.khalf + 1,))
    hd = pc.hdur_corr if pc.hdur_corr > 0 else 8.0 * pc.dt_sub
    kh = gf3d.stf_kernel(gf3d.GF_STF_HEAVI, hd, pc.dt_sub, pc.trunc)
    m = (len(kh) - 1) // 2
    ok("a Heaviside kernel is wider than one tap", m > 0)
    ok("w(0) = 1/2", kh[m] == 0.5)
    ok("w(j) + w(-j) = 1, bitwise", all(kh[m + j] + kh[m - j] == 1.0 for j in range(1, m + 1)))
    kd = gf3d.stf_kernel(pc.kind_stf, 2.0 * pc.hdur_corr, pc.dt_sub, pc.trunc)
    ok("any parameters: twice the width, its own length",
       pc.hdur_corr <= 0
       or len(kd) == 2 * math.ceil(pc.trunc * (2.0 * pc.hdur_corr) / pc.dt_sub) + 1)
    ok("none is [1]", np.array_equal(gf3d.stf_kernel(gf3d.GF_STF_NONE, 3.0, 0.5, 4.0), [1.0]))
    ok("hdur 0, Heaviside: [1/2], the trapezoid",
       np.array_equal(gf3d.stf_kernel(gf3d.GF_STF_HEAVI, 0.0, 0.5, 4.0), [0.5]))
    ok("hdur 0, Gaussian: [1]",
       np.array_equal(gf3d.stf_kernel(gf3d.GF_STF_GAUSS, 0.0, 0.5, 4.0), [1.0]))
    ok_raises("dt 0", gf3d.GF_ERR_ARG, gf3d.stf_kernel, gf3d.GF_STF_HEAVI, 1.0, 0.0, 4.0)
    ok_raises("trunc 0", gf3d.GF_ERR_ARG, gf3d.stf_kernel, gf3d.GF_STF_HEAVI, 1.0, 0.5, 0.0)
    ok_raises("kind 5", gf3d.GF_ERR_ARG, gf3d.stf_kernel, 5, 1.0, 0.5, 4.0)

    db.close()
    ok_raises("a closed database", gf3d.GF_ERR_ARG, db.weights, cmt)


def main(argv):
    if len(argv) < 4:
        print(__doc__)
        return 2

    xgf3d = Path(argv[1]).resolve()
    dbpath = Path(argv[2]).resolve()
    cmtpath = Path(argv[3]).resolve()
    forcepath = Path(argv[4]).resolve() if len(argv) > 4 and argv[4] else None

    print()
    print(" ******************************")
    print(" test_gf_python")
    print(" ******************************")
    print()

    # the element block and the weights alone, on a database the full run did not use
    if len(argv) > 5 and argv[5] == "block":
        cmt = gf3d.CMTSource.read(cmtpath)
        force = (gf3d.ForceSource.read(forcepath)
                 if forcepath is not None and forcepath.exists() else None)
        print(f" 10. element_block, on {dbpath}")
        element_block_section(dbpath, cmt)
        print(f"\n 11. weights, on {dbpath}")
        weights_section(dbpath, cmt, force)
        print()
        if nfail == 0:
            print(" test_gf_python: all assertions passed\n")
            return 0
        print(f" test_gf_python: {nfail} assertion(s) FAILED\n")
        return 1

    # ------------------------------------------------------------------
    print(" 1. the package and the library it found")
    print(f"       library: {gf3d.library_path}")
    print(f"       version: {gf3d.library_version()}")
    ok("a version string came back", len(gf3d.library_version()) > 0)
    ok("the library's API version agrees with the package's",
       gf3d._lib.API_VERSION == gf3d._lib.lib.gf3d_api_version())
    # the struct layout check runs at import: a disagreement raises there

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
    # not Mrr in particular: a pure strike-slip source has Mrr = 0
    ok("the moment tensor is six numbers in dyne-cm",
       len(cmt.tensor) == 6 and max(abs(m) for m in cmt.tensor) > 0)
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
    # t0=None asks the library for specfem's own rule and gets the number back
    ok("t0_req is specfem's own start time",
       abs(plan.t0_req - 1.5 * cmt.hdur) <= 1e-12 * 1.5 * cmt.hdur)
    ok("an explicit t0 is reported as asked",
       abs(db.plan(cmt, t0=120.0).t0_req - 120.0) <= 1e-12 * 120.0)
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
    # index survives: that is an O(1) error. The bound is rounding: the
    # seismogram and each partial come from their own weight vector, and the
    # generated fixture's traces are up to ~5e5 times smaller than what is
    # summed into them -- a $GF3D_TEST_GFDB may be worse conditioned
    # (test_gf_partials_db, section 5, asserts it against that size),
    # so eps * 5e5 * a few is ~1e-10 of the peak.
    lin = (r.dp[:, :6] * np.asarray(cmt.tensor)[None, :, None, None]).sum(axis=1)
    ok_err("sum(M_v dp_v) reproduces the seismogram",
           np.abs(lin - r.data).max() / np.abs(r.data).max(), tol=1e-9)

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
        print("       skipped: obspy not installed")
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
    print("\n 8. keeping elements between extractions")

    # the eviction order is test_gf_cache's to pin; this is that the
    # package passes max_elements through, reports what the library counts,
    # and that a caching handle's arrays are the uncached one's, bitwise
    info = db.info
    ok("info reports bytes_per_element", info["bytes_per_element"] > 0)
    ok("an unnamed max_elements keeps nothing", db.max_elements == 0)

    ok_raises_py("max_elements = -1 refused", ValueError, gf3d.Database, dbpath, max_elements=-1)
    ok_raises_py("max_elements = 1.5 refused", TypeError, gf3d.Database, dbpath, max_elements=1.5)
    ok_raises_py("max_elements = True refused", TypeError, gf3d.Database, dbpath, max_elements=True)
    with gf3d.Database(dbpath, max_elements=2**40) as huge:
        ok("more than a C int opens, as 'all of them'", not huge.closed)

    # B: walk north from A until the element changes
    here = db.locate(cmt.latitude, cmt.longitude, cmt.depth)
    far = None
    for iw in range(1, 81):
        cand = dataclasses.replace(cmt, latitude=cmt.latitude + 0.25 * iw)
        if cand.latitude > 90.0:
            break
        try:
            if db.locate(cand.latitude, cand.longitude, cand.depth).ielem != here.ielem:
                far = cand
                break
        except gf3d.GF3DError:
            continue
    ok("a second element is reachable", far is not None)

    if far is not None and 3 * info["bytes_per_element"] < 2**30:
        dbc = gf3d.Database(dbpath, max_elements=2)
        ok("max_elements is kept", dbc.max_elements == 2)
        ok("a fresh cache has done nothing",
           dbc.cache_stats == dict(hits=0, misses=0, evictions=0, n_cached=0, files_read=0))

        same = True
        for k, src in enumerate((cmt, far, cmt, far)):
            if k == 2:
                st0 = dbc.cache_stats
            want, got = db.partials(src), dbc.partials(src)
            same = same and np.array_equal(want.data, got.data) and np.array_equal(want.dp, got.dp)
            if k == 2:
                st1 = dbc.cache_stats
                files_hit = st1["files_read"] - st0["files_read"]
        ok("A B A B: data and dp identical to the bit, cached or not", same)
        ok("returning to A was a hit", st1["hits"] - st0["hits"] == 1)
        # The handle also keeps the coordinates of the max(10, max_elements)
        # elements its locates used last (GF_NCOORD_MIN = 10). When that is
        # every element, as on the fixture, returning to A reads nothing at
        # all; on a larger database the locates at A and B may have tried
        # more candidates than it holds, and then only rereads of
        # coordinates, at most one locate's 10, are allowed.
        if info["nelem"] <= max(10, dbc.max_elements):
            ok("returning to A read no element file", files_hit == 0)
        else:
            ok("returning to A read no displacement, only coordinates", files_hit <= 10)

        st = dbc.cache_stats
        print(f"       cache_stats: {st}")
        ok("2 hits, 2 misses, 0 evictions, 2 held",
           (st["hits"], st["misses"], st["evictions"], st["n_cached"]) == (2, 2, 0, 2))

        with gf3d.Database(dbpath, max_elements=2) as other:
            ok("another handle has its own cache", other.cache_stats["misses"] == 0)

        dbc.close()
        ok_raises("cache_stats of a closed database", gf3d.GF_ERR_ARG, lambda: dbc.cache_stats)
    else:
        print("       no second element, or two elements exceed 1 GiB: not exercised")

    # ------------------------------------------------------------------
    print("\n 9. closing, and a second database")

    # A second handle on the same directory. Opening it re-installs specfem's
    # process-wide globals and takes the one kd-tree; the library re-installs
    # the first handle's state on every call, and if it did not, this
    # extraction would come back subtly different. The second open need not
    # be a *different* database for that -- it is a different open, which is
    # what the library keys on -- and requiring a second example meant this
    # never ran, both being gitignored.
    db2 = gf3d.Database(dbpath)
    print(f"       second handle: {db2!r}")
    r2 = db2.seismograms(cmt)
    ok("the second handle extracts", np.isfinite(r2.data).all())

    again = db.partials(cmt)
    ok("the first handle is bit-for-bit unchanged",
       np.array_equal(again.data, r.data) and np.array_equal(again.dp, r.dp))
    db2.close()

    db.close()
    ok("the database is closed", db.closed)
    db.close()
    ok("closing twice is harmless", db.closed)
    ok_raises("using a closed database", gf3d.GF_ERR_ARG, db.locate, 0.0, 0.0, 10.0)

    with gf3d.Database(dbpath) as ctx:
        ok("the context manager opens", not ctx.closed)
    ok("and closes on the way out", ctx.closed)

    # ------------------------------------------------------------------
    print("\n 10. element_block")

    element_block_section(dbpath, cmt)

    print("\n 11. weights")

    weights_section(dbpath, cmt, force if forcepath is not None and forcepath.exists() else None)

    print()
    if nfail == 0:
        print(" test_gf_python: all assertions passed\n")
        return 0
    print(f" test_gf_python: {nfail} assertion(s) FAILED\n")
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
