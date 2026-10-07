#!/usr/bin/env python
"""Using gf3d from a sampler: contract an element's raw data yourself.

`db.seismograms(src)` locates the source, reads its element, contracts it
with weights and applies the source time function. A sampler that moves a
source a little between evaluations repeats the expensive part (the read)
for an element it already had. The library therefore gives the three
pieces separately:

  W = db.weights(src, kind)       where the source is, and the weights
  u = db.element_block(ielem)     the element's raw data, u[s, a, t, m]
  k = gf3d.stf_kernel(...)        the source time function's taps

and the caller contracts, with `W.scale[s] * u[s, a, t, :] @ W.w` for the
seismogram and `... @ W.dw[c]` for partial c, then converts with `k`. The
conversion is linear and acts on time only, so it can equally be applied
to the element's data first, once per element, and every source in that
element is then a single contraction.

This script does it in numpy and checks it against the library, then runs
a small sampler loop that reads and converts each element once.

  1. the library's way, and the pieces' costs
  2. the contraction and the conversion in numpy, against the library
  3. the same for a force source, if the example has a FORCESOLUTION
  4. 200 moved sources, each element read and converted once

Usage
-----
  cd EXAMPLES/green_function_database/global/api_demos/python
  ../../../.venv/bin/python gf3d_sampler_demo.py [--db GFDB] [--png FILE] [--no-plot]

Not part of the Snakemake workflow; run it after the workflow has built
`../../GFDB`. Requires numpy, matplotlib (for the figure) and lib/libgf3d.so.
"""

from __future__ import annotations

import argparse
import copy
import sys
import time
from pathlib import Path

import numpy as np

try:
    import gf3d
except ImportError:
    _repo = Path(__file__).resolve().parents[5]
    sys.path.insert(0, str(_repo / "utils" / "green_function"))
    import gf3d

HERE = Path(__file__).resolve().parent
EXAMPLE = HERE.parents[1]          # .../green_function_database/global
DEFAULT_DB = EXAMPLE / "GFDB"
DATA = EXAMPLE / "validation_data"


def rule(title):
    print()
    print(title)
    print("-" * len(title))


def convert(x, k, p):
    """The library's source time function conversion of traces x[..., nt].

    The formula of the `gf3d.stf_kernel` docstring, over any leading axes.
    Returns the first p.nt samples, the library's output axis.
    """
    khalf = (len(k) - 1) // 2
    lead = x.shape[:-1]
    xpad = np.concatenate([np.zeros(lead + (p.npad,)), x], axis=-1)
    flat = xpad.reshape(-1, xpad.shape[-1])
    y = np.empty_like(flat)
    for i, tr in enumerate(flat):
        if p.kind_stf == gf3d.GF_STF_NONE:
            y[i] = tr
            continue
        win = np.convolve(tr, k)[khalf:khalf + len(tr)]       # zero outside
        if p.kind_stf == gf3d.GF_STF_GAUSS:
            y[i] = win
        else:                                                  # Heaviside
            P = np.concatenate([np.zeros(khalf + 1), np.cumsum(tr)])[:len(tr)]
            y[i] = p.dt_sub * (P + win)                        # P[i-khalf-1]
    return y.reshape(lead + (-1,))[..., :p.nt]


def contract(u, W, ndw=None):
    """Seismograms (s, a, t) and partials (s, c, a, t) on the stored grid."""
    u64 = u.astype(np.float64)
    x = W.scale[:, None, None] * (u64 @ W.w)
    if not len(W.dw):
        return x, None
    dx = W.scale[:, None, None, None] * (u64 @ W.dw.T)        # (s, a, t, c)
    return x, np.moveaxis(dx, -1, 1)


def table(rows):
    print(f"    {'column':<8s} {'unit':<12s} {'max rel diff':>12s}")
    for name, unit, d in rows:
        print(f"    {name:<8s} {unit:<12s} {d:12.2e}")


def rel(a, b):
    """Largest |a - b|, relative to the largest |b| over all traces."""
    return np.abs(a - b).max() / np.abs(b).max()


def check_source(db, src, kind, label):
    """Contract and convert in numpy, compare with the library. Returns
    the worst relative difference and what the plot needs."""
    rule(label)
    t0 = time.perf_counter()
    W = db.weights(src, kind=kind)
    t_w = (time.perf_counter() - t0) * 1e3
    t0 = time.perf_counter()
    u = db.element_block(W.location.ielem)
    t_u = (time.perf_counter() - t0) * 1e3
    print(f"  element {W.location.ielem}   block {u.shape} {u.dtype}, {u.nbytes / 2**20:.1f} MB")
    print(f"  db.weights        {t_w:8.2f} ms")
    print(f"  db.element_block  {t_u:8.2f} ms   (read from disk, every call)")

    p = db.plan(src)
    k = gf3d.stf_kernel(p.kind_stf, p.hdur_corr, p.dt_sub, p.trunc)
    print(f"  source time function: kind {p.kind_stf}, hdur {p.hdur_corr:.3g} s, "
          f"{len(k)} taps, npad {p.npad}")

    t0 = time.perf_counter()
    x, dx = contract(u, W)
    t_c = (time.perf_counter() - t0) * 1e3
    t0 = time.perf_counter()
    seis = convert(x, k, p)
    pars = convert(dx, k, p) if dx is not None else None
    t_k = (time.perf_counter() - t0) * 1e3
    print(f"  contraction       {t_c:8.2f} ms   conversion {t_k:.2f} ms")

    ref = db.seismograms(src).data
    rows = [("seis", "m", rel(seis, ref))]
    if pars is not None:
        # The time column (slot 9) has no weight: a centroid-time change
        # shifts the trace, so it is not a contraction of the element data.
        dp = db.partials(src, kind=kind).dp
        for c, (name, unit) in enumerate(zip(W.dw_names, W.dw_units)):
            rows.append((name, unit, rel(pars[:, c], dp[:, c])))
    print()
    table(rows)
    return max(d for *_, d in rows), seis, ref, p, W


def sampler(db, cmt, n=200, seed=1):
    """n sources moved around the CMT; each element read and converted once."""
    rule(f"4. {n} moved sources: each element read and converted once")
    rng = np.random.default_rng(seed)
    p = db.plan(cmt)             # hdur is the CMT's for every source, so is k
    k = gf3d.stf_kernel(p.kind_stf, p.hdur_corr, p.dt_sub, p.trunc)
    cache = {}                   # ielem -> U[s, a, t, m], the converted basis
    reads = 0
    t_fill = 0.0
    worst = 0.0
    t_total = 0.0
    done = 0
    for i in range(n):
        s = copy.deepcopy(cmt)
        s.latitude += rng.uniform(-0.05, 0.05)
        s.longitude += rng.uniform(-0.05, 0.05)
        s.depth += rng.uniform(-2.0, 2.0)
        t0 = time.perf_counter()
        try:
            W = db.weights(s, kind=0)
        except gf3d.GF3DError:
            continue                                    # left the database
        ielem = W.location.ielem
        if ielem not in cache:
            # once per element: read, convert every (s, a, m) trace along time
            t1 = time.perf_counter()
            u = db.element_block(ielem).astype(np.float64)
            cache[ielem] = np.moveaxis(convert(np.moveaxis(u, 2, 3), k, p), 3, 2)
            reads += 1
            t_fill += time.perf_counter() - t1
        y = W.scale[:, None, None] * (cache[ielem] @ W.w)     # the seismograms
        t_total += time.perf_counter() - t0
        done += 1
        if i < 3:          # the first few against the library, outside the timing
            worst = max(worst, rel(y, db.seismograms(s).data))
    print(f"  {done} sources evaluated, {len(cache)} distinct elements, {reads} reads")
    print(f"  reading and converting an element: {t_fill / reads * 1e3:.1f} ms each")
    print(f"  per evaluation otherwise (weights + one contraction): "
          f"{(t_total - t_fill) / done * 1e3:.2f} ms")
    print(f"  the first three against db.seismograms: {worst:.2e}")
    return worst


def plot(seis, ref, p, png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = p.t_first + np.arange(p.nt) * p.dt_sub
    fig, ax = plt.subplots(figsize=(10, 3.5))
    ax.plot(t, ref[0, 2], "k-", lw=1.5, label="db.seismograms")
    ax.plot(t, seis[0, 2], "r--", lw=1.0, label="element_block + weights + stf_kernel, numpy")
    ax.set_xlabel("time relative to the centroid time [s]")
    ax.set_ylabel("displacement [m]")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.savefig(png, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {png}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--db", default=str(DEFAULT_DB), help="the GFDB directory")
    ap.add_argument("--cmt", default=str(DATA / "CMTSOLUTION"), help="a CMTSOLUTION")
    ap.add_argument("--force", default=str(DATA / "FORCESOLUTION"), help="a FORCESOLUTION")
    ap.add_argument("--png", default=str(HERE / "gf3d_sampler_demo.png"), help="figure to write")
    ap.add_argument("--no-plot", action="store_true", help="numbers only")
    args = ap.parse_args(argv)

    if not (Path(args.db) / "mesh_info.h5").exists():
        raise SystemExit(f"no database at {args.db}; build it with snakemake in {EXAMPLE}")
    cmt = gf3d.CMTSource.read(args.cmt)

    worst = 0.0
    with gf3d.Database(args.db, max_elements=0) as db:
        d, seis, ref, p, _ = check_source(
            db, cmt, 2, "1-2. the CMT: weights, element_block, numpy, against the library")
        worst = max(worst, d)

        if Path(args.force).exists():
            force = gf3d.ForceSource.read(args.force)
            d, *_ = check_source(db, force, 0, "3. the force source (seismogram only)")
            worst = max(worst, d)

        worst = max(worst, sampler(db, cmt))

        if not args.no_plot:
            plot(seis, ref, p, args.png)

    print(f"\nlargest relative difference: {worst:.2e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
