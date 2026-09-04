#!/usr/bin/env python
"""Compare xgf3d --seis output against forward-modelled seismograms.

This is the physical acceptance gate of the Green function extraction
library (gf3df_integration_plan/stage_05_stf.md): the reconstructed trace
for a source file, against the trace specfem itself computes for that source
on the same mesh.

Truth sources, and nothing else:

  * the forward SAC traces in <forward>/OUTPUT_FILES/,
  * mesh_info.h5 and stations/NET.STA.h5 of the database,
  * the writer's own anti-alias filter, replicated line by line from
    src/specfem3D/green_function_stf.F90 (gf_butterworth_sos, gf_sosfiltfilt)
    and checked against the filtered STF the writer stored (--preflight).

Nothing is imported from gf_cross_validate.py: its comparison utilities taper
the two traces on different supports and detrend by endpoints, and inheriting
them would make those choices "the answer".

What is compared, and why on the coarse grid
--------------------------------------------
The reconstruction equals the forward trace lowpassed by the database's
Butterworth (the filter was applied to the reciprocal run's source time
function, so it band-limits the whole stored field). The forward trace is
therefore filtered with that same filter, at its native 0.1 s, with the same
zero initial conditions and no taper; its right-edge transient is excluded
from the window rather than hidden.

The comparison is made on the database's own sample times. The filtered
forward trace is oversampled a hundred-fold relative to its content (nothing
survives above f_cutoff, and the source spectrum is dead long before that),
so a cubic spline evaluated at the GF times is exact to ~1e-9 relative; the
GF trace is never interpolated. Both traces are band-limited below the coarse
Nyquist, so by Parseval the relative L2 on the coarse grid is the continuous
one, up to the window edges. Comparing at 0.1 s would only add the GF's own
interpolation error to the metric.

Nothing is detrended and no mean is removed: the static offset of a
Heaviside response is physics and part of the answer.

Output
------
Per station and component: relative L2 misfit ||g - f|| / ||f||, the best-fit
amplitude ratio <g,f>/<f,f> and the misfit after applying it, the peak ratio,
the cross-correlation lag (parabolic sub-sample refinement; positive means
the GF trace is *later* than the forward one), and the largest residual over
the forward peak. A rigid time shift and an amplitude error look identical in
relative L2; the lag and the amplitude ratio tell them apart, which is what
turns a number into a diagnosis.

Usage
-----
  gf_compare.py --gf <dir with NET.STA.gf3d.txt> --fwd <forward OUTPUT_FILES>
                --db <GFDB> [--json out.json] [--png <dir>] [--threshold X]
                [--stations NET.STA ...] [--trim-seconds S]

Exit status 1 if --threshold is given and the worst relative L2 exceeds it.
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import h5py
from scipy.interpolate import CubicSpline

COMPONENTS = ("N", "E", "Z")
CHANNEL_PREFIXES = ("BX", "BH", "MX")
FILTER_ORDER = 4


# ---------------------------------------------------------------------------
# The writer's filter, replicated from src/specfem3D/green_function_stf.F90
# ---------------------------------------------------------------------------

def gf_butterworth_sos(order, fc, fs):
    """gf_butterworth_sos (green_function_stf.F90:145-213), verbatim.

    Second-order sections of a Butterworth lowpass via the bilinear transform
    with pre-warped cutoff. Each row is [b0, b1, b2, 1, a1, a2].
    """
    nsections = order // 2
    wc = 2.0 * fs * np.tan(np.pi * fc / fs)
    wc2 = wc * wc
    sos = np.zeros((nsections, 6))
    for k in range(1, nsections + 1):
        pole_real = np.cos(np.pi * (2 * k + order - 1) / (2 * order)) * wc
        b0 = wc2
        b1 = 2.0 * wc2
        b2 = wc2
        a0 = 4.0 * fs * fs - 4.0 * fs * pole_real + wc2
        a1 = 2.0 * wc2 - 8.0 * fs * fs
        a2 = 4.0 * fs * fs + 4.0 * fs * pole_real + wc2
        sos[k - 1] = [b0 / a0, b1 / a0, b2 / a0, 1.0, a1 / a0, a2 / a0]
    return sos


def gf_sos_filter_forward(x, b0, b1, b2, a1, a2):
    """gf_sos_filter_forward (green_function_stf.F90:250-280): direct form II
    transposed, zero initial conditions."""
    y = np.empty_like(x)
    w1 = 0.0
    w2 = 0.0
    for i in range(len(x)):
        yi = b0 * x[i] + w1
        w1 = b1 * x[i] - a1 * yi + w2
        w2 = b2 * x[i] - a2 * yi
        y[i] = yi
    return y


def gf_sosfiltfilt(x, sos):
    """gf_sosfiltfilt (green_function_stf.F90:219-244): per section, forward
    pass, then the reversed array through the same forward pass. No padding,
    no taper, zero initial conditions on every pass."""
    y = np.array(x, dtype=np.float64)
    for sec in sos:
        b0, b1, b2, _, a1, a2 = sec
        y = gf_sos_filter_forward(y, b0, b1, b2, a1, a2)
        y = gf_sos_filter_forward(y[::-1], b0, b1, b2, a1, a2)[::-1]
    return y


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def read_gf3d(path):
    """Parse one NET.STA.gf3d.txt: the '#' header into a dict, the four
    columns into arrays."""
    header = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            body = line[1:].strip()
            if ":" in body:
                key, _, val = body.partition(":")
                header[key.strip()] = val.strip()
    data = np.loadtxt(path, comments="#")
    if data.ndim != 2 or data.shape[1] != 4:
        raise SystemExit(f"{path}: expected four columns, got shape {data.shape}")

    hd = [float(v) for v in header["stf hdur"].split()]
    kern = header["stf kernel"].split()
    ax = header["axis"].split()
    ax0 = [float(v) for v in header["axis t0"].split()]
    plan = {
        "kind": header["stf kind"],
        "hdur_src": hd[0], "hdur_target": hd[1], "hdur_db": hd[2], "hdur_corr": hd[3],
        "trunc": float(kern[0]), "khalf": int(kern[1]), "guard": kern[2].upper() == "T",
        "onset": float(header["stf onset"]),
        "note": header.get("stf note", ""),
        "dt": float(ax[0]), "subsample_step": int(ax[1]),
        "nt_db": int(ax[2]), "npad": int(ax[3]), "nt": int(ax[4]),
        "t0_db": ax0[0], "t0_req": ax0[1], "t0": ax0[2], "t_first": ax0[3],
        "element": header["element"],
    }
    if data.shape[0] != plan["nt"]:
        raise SystemExit(f"{path}: header says nt = {plan['nt']}, file has {data.shape[0]} rows")
    traces = {c: data[:, i + 1] for i, c in enumerate(COMPONENTS)}
    return header["station"], plan, data[:, 0], traces


def _scalar(v):
    """The writer stores scalar attributes as arrays of shape (1,)."""
    return np.asarray(v).reshape(-1)[0]


def read_db(db, station):
    """The database's own statement of dt, t0, subsample_step and the
    station's hdur, f_cutoff and stored (filtered) STF."""
    db = Path(db)
    with h5py.File(db / "mesh_info.h5", "r") as f:
        a = f.attrs
        mesh = {"dt": float(_scalar(a["dt"])), "t0": float(_scalar(a["t0"])),
                "nstep": int(_scalar(a["nstep"])), "subsample_step": int(_scalar(a["subsample_step"]))}
    with h5py.File(db / "stations" / f"{station}.h5", "r") as f:
        a = f.attrs
        sta = {"hdur": float(_scalar(a["hdur"])), "f_cutoff": float(_scalar(a["f_cutoff"])),
               "stf": np.array(f["stf"][:], dtype=np.float64)}
    return mesh, sta


def read_forward(fwd_dir, station, dt):
    """The forward SAC traces for a station, with a double-precision time axis.

    The SAC header's delta is float32; accumulating it over 19100 samples
    drifts by 3e-5 s, so the axis is built from the database's dt instead
    and delta is only checked against it.
    """
    from obspy import read as obspy_read

    fwd_dir = Path(fwd_dir)
    prefix = None
    for p in CHANNEL_PREFIXES:
        if (fwd_dir / f"{station}.{p}N.sem.sac").exists():
            prefix = p
            break
    if prefix is None:
        raise SystemExit(f"no forward SAC files for {station} in {fwd_dir}")

    traces = {}
    t_f = None
    for c in COMPONENTS:
        tr = obspy_read(str(fwd_dir / f"{station}.{prefix}{c}.sem.sac"))[0]
        delta = float(tr.stats.delta)
        if abs(delta - dt) > 1e-6 * dt:
            raise SystemExit(f"{station}.{prefix}{c}: SAC delta {delta} differs from the database dt {dt}")
        b = float(tr.stats.sac.b)
        t = b + np.arange(tr.stats.npts) * dt
        if t_f is None:
            t_f = t
        elif t_f.shape != t.shape or abs(t_f[0] - t[0]) > 1e-9:
            raise SystemExit(f"{station}: the three forward components have different axes")
        traces[c] = tr.data.astype(np.float64)
    return t_f, traces, prefix


# ---------------------------------------------------------------------------
# The checks
# ---------------------------------------------------------------------------

def preflight_filter(mesh, sta):
    """Rebuild the writer's stored STF from the analytic Gaussian and the
    replicated filter, and compare with what the writer stored (float32).

    This shares no code with the writer, so a wrong pre-warp or pole angle
    in the replication fails here, before any waveform is compared.
    """
    nstep, dt, t0, hdur = mesh["nstep"], mesh["dt"], mesh["t0"], sta["hdur"]
    t = (np.arange(1, nstep + 1) - 1) * dt - t0
    raw = np.exp(-(t / hdur) ** 2) / (np.sqrt(np.pi) * hdur)
    sos = gf_butterworth_sos(FILTER_ORDER, sta["f_cutoff"], 1.0 / dt)
    rebuilt = gf_sosfiltfilt(raw, sos)
    stored = sta["stf"]
    if stored.shape != rebuilt.shape:
        raise SystemExit(f"stored stf has {stored.shape[0]} samples, nstep is {nstep}")
    err = np.max(np.abs(rebuilt - stored)) / np.max(np.abs(stored))
    return err


def lag_by_xcorr(g, f, dt_sub):
    """Cross-correlation lag of g relative to f, in samples, refined with a
    parabola through the peak. Positive: g is later than f.

    The first differences are correlated, not the traces: a Heaviside
    response carries a static offset, and even after removing the mean a
    finite record keeps a pedestal whose self-correlation is a triangle
    peaking at zero lag, which pulls the peak inward (0.06 samples on the
    self-check below). Differencing removes DC exactly and leaves the
    impulse responses, whose correlation peak is where the delay is.
    """
    dg = np.diff(g)
    df = np.diff(f)
    n = len(dg)
    c = np.correlate(dg, df, mode="full")           # peaks at index (n - 1) + d for g delayed by d
    k = int(np.argmax(c))
    if 0 < k < len(c) - 1:
        cm, c0, cp = c[k - 1], c[k], c[k + 1]
        denom = cm - 2.0 * c0 + cp
        delta = 0.5 * (cm - cp) / denom if denom != 0.0 else 0.0
    else:
        delta = 0.0
    lag_samples = (k - (n - 1)) + delta
    return lag_samples, lag_samples * dt_sub


def lag_convention_selfcheck():
    """The sign convention above, asserted on a synthetic delay."""
    t = np.arange(400) * 0.5
    f = np.exp(-((t - 60.0) / 8.0) ** 2)
    g = np.exp(-((t - 61.5) / 8.0) ** 2)     # g is 3 samples (1.5 s) later
    lag_s, lag_sec = lag_by_xcorr(g, f, 0.5)
    if abs(lag_s - 3.0) > 0.05:
        raise SystemExit(f"lag convention self-check failed: {lag_s} samples for a +3 shift")


def compare_component(g, f):
    norm_f = np.sqrt(np.sum(f * f))
    if norm_f == 0.0:
        return None
    resid = g - f
    rel_l2 = np.sqrt(np.sum(resid * resid)) / norm_f
    amp = float(np.dot(g, f) / np.dot(f, f))
    resid_s = g - amp * f
    rel_l2_scaled = np.sqrt(np.sum(resid_s * resid_s)) / norm_f
    peak_f = np.max(np.abs(f))
    return {
        "rel_l2": float(rel_l2),
        "amp_ratio": amp,
        "rel_l2_after_scale": float(rel_l2_scaled),
        "peak_ratio": float(np.max(np.abs(g)) / peak_f),
        "max_resid_over_peak": float(np.max(np.abs(resid)) / peak_f),
    }


def compare_station(station, plan, t_gf, gf, t_f, fwd, sos, trim_seconds):
    dt_sub = plan["dt"] * plan["subsample_step"]
    khalf = plan["khalf"]

    # the forward trace, band-limited exactly as the database is
    fwd_filt = {c: gf_sosfiltfilt(fwd[c], sos) for c in COMPONENTS}

    # the window on the coarse grid: inside the forward record, clear of
    # the GF's trailing zero-extended samples and of the filter's right-edge
    # transient on the forward trace
    t_end = t_f[-1] - max(khalf * dt_sub, trim_seconds)
    idx = np.where((t_gf >= t_f[0]) & (t_gf <= t_end))[0]
    if khalf > 0:
        idx = idx[idx < plan["nt"] - khalf]
    if len(idx) < 10:
        raise SystemExit(f"{station}: fewer than 10 overlapping samples")

    result = {
        "plan": plan,
        "window": {"t_start": float(t_gf[idx[0]]), "t_end": float(t_gf[idx[-1]]),
                   "n": int(len(idx)), "dt_sub": dt_sub},
    }
    aligned = {}
    for c in COMPONENTS:
        spl = CubicSpline(t_f, fwd_filt[c])
        f_on_gf = spl(t_gf[idx])
        g = gf[c][idx]
        m = compare_component(g, f_on_gf)
        if m is None:
            result[c] = {"rel_l2": None}
            continue
        lag_s, lag_sec = lag_by_xcorr(g, f_on_gf, dt_sub)
        m["lag_samples"] = float(lag_s)
        m["lag_seconds"] = float(lag_sec)
        result[c] = m
        aligned[c] = (t_gf[idx], g, f_on_gf)
    return result, aligned


def plot_station(station, result, aligned, png_dir, source_label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    for ax, c in zip(axes, COMPONENTS):
        if c not in aligned:
            continue
        t, g, f = aligned[c]
        m = result[c]
        ax.plot(t, f, "k-", lw=0.9, label="forward (writer-filtered)")
        ax.plot(t, g, "r--", lw=0.9, label="gf3d")
        ax.plot(t, (g - f) * 10.0, "b-", lw=0.6, alpha=0.6, label="residual x10")
        ax.set_ylabel(f"{c} [m]")
        ax.set_title(f"{c}: rel L2 {m['rel_l2']:.3e}   amp {m['amp_ratio']:.5f}   "
                     f"lag {m['lag_seconds']:+.3f} s   max resid/peak {m['max_resid_over_peak']:.2e}",
                     fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right", fontsize=7)
    axes[-1].set_xlabel("time after origin [s]")
    p = result["plan"]
    fig.suptitle(f"{station}  {source_label}  {p['kind']}  hdur_corr {p['hdur_corr']:.3f} s  "
                 f"element {p['element']}", fontsize=10)
    fig.tight_layout()
    out = Path(png_dir) / f"{station}.gf_compare.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gf", required=True, help="directory holding NET.STA.gf3d.txt from xgf3d --seis")
    ap.add_argument("--fwd", required=True, help="the forward run's OUTPUT_FILES directory")
    ap.add_argument("--db", required=True, help="the Green function database directory")
    ap.add_argument("--json", help="write the results here")
    ap.add_argument("--png", help="write one overlay plot per station into this directory")
    ap.add_argument("--threshold", type=float, help="exit 1 if the worst relative L2 exceeds this")
    ap.add_argument("--stations", nargs="*", help="restrict to these NET.STA")
    ap.add_argument("--trim-seconds", type=float, default=60.0,
                    help="exclude at least this much of the forward record's end (filter transient)")
    ap.add_argument("--preflight-tol", type=float, default=1e-6,
                    help="allowed relative error when rebuilding the writer's stored STF")
    args = ap.parse_args()

    lag_convention_selfcheck()

    gf_files = sorted(Path(args.gf).glob("*.gf3d.txt"))
    if args.stations:
        gf_files = [p for p in gf_files if p.name[: -len(".gf3d.txt")] in set(args.stations)]
    if not gf_files:
        raise SystemExit(f"no *.gf3d.txt in {args.gf}")

    out = {"gf": str(args.gf), "fwd": str(args.fwd), "db": str(args.db),
           "filter": {"order": FILTER_ORDER}, "stations": {}}
    worst = 0.0
    sos = None

    for gf_path in gf_files:
        station, plan, t_gf, gf = read_gf3d(gf_path)
        mesh, sta = read_db(args.db, station)

        # the header and the database must agree on the axis
        for key, val in (("dt", mesh["dt"]), ("subsample_step", mesh["subsample_step"]),
                         ("t0_db", mesh["t0"]), ("hdur_db", sta["hdur"])):
            if abs(plan[key] - val) > 1e-9 * max(1.0, abs(val)):
                raise SystemExit(f"{station}: header {key} = {plan[key]} but the database says {val}")
        fc_expected = 1.0 / (2.0 * mesh["dt"] * mesh["subsample_step"])
        if abs(sta["f_cutoff"] - fc_expected) > 1e-9 * fc_expected:
            raise SystemExit(f"{station}: f_cutoff {sta['f_cutoff']} is not 1/(2 dt subsample_step) = {fc_expected}")

        if sos is None:
            err = preflight_filter(mesh, sta)
            out["filter"].update({"f_cutoff": sta["f_cutoff"], "fs": 1.0 / mesh["dt"],
                                  "preflight_stf_rebuild_rel_err": float(err)})
            print(f"preflight: writer's stored STF rebuilt with the replicated filter, "
                  f"max relative error {err:.3e} (float32 storage; tolerance {args.preflight_tol:.0e})")
            if err > args.preflight_tol:
                raise SystemExit("the filter replication does not reproduce the writer's STF; stopping")
            sos = gf_butterworth_sos(FILTER_ORDER, sta["f_cutoff"], 1.0 / mesh["dt"])

        t_f, fwd, prefix = read_forward(args.fwd, station, mesh["dt"])
        result, aligned = compare_station(station, plan, t_gf, gf, t_f, fwd, sos, args.trim_seconds)
        result["channel_prefix"] = prefix
        out["stations"][station] = result

        w = result["window"]
        print(f"\n{station}  ({plan['kind']}, hdur_corr {plan['hdur_corr']:.4f} s, khalf {plan['khalf']}, "
              f"onset {plan['onset']:.2e}, guard {plan['guard']})")
        print(f"  window {w['t_start']:.2f} .. {w['t_end']:.2f} s, {w['n']} samples at {w['dt_sub']} s")
        print(f"  {'comp':4s} {'rel L2':>11s} {'amp ratio':>11s} {'L2 scaled':>11s} {'peak ratio':>11s} "
              f"{'lag [smp]':>10s} {'lag [s]':>9s} {'maxres/pk':>10s}")
        for c in COMPONENTS:
            m = result[c]
            if m.get("rel_l2") is None:
                print(f"  {c:4s} (forward trace is identically zero)")
                continue
            print(f"  {c:4s} {m['rel_l2']:11.4e} {m['amp_ratio']:11.6f} {m['rel_l2_after_scale']:11.4e} "
                  f"{m['peak_ratio']:11.6f} {m['lag_samples']:+10.3f} {m['lag_seconds']:+9.3f} "
                  f"{m['max_resid_over_peak']:10.3e}")
            worst = max(worst, m["rel_l2"])

        if args.png:
            Path(args.png).mkdir(parents=True, exist_ok=True)
            label = "CMT" if plan["kind"] == "heaviside" else "force"
            p = plot_station(station, result, aligned, args.png, label)
            print(f"  plot: {p}")

    out["worst_rel_l2"] = float(worst)
    print(f"\nworst relative L2 over all stations and components: {worst:.4e}")

    if args.json:
        Path(args.json).parent.mkdir(parents=True, exist_ok=True)
        with open(args.json, "w") as f:
            json.dump(out, f, indent=2)
        print(f"wrote {args.json}")

    if args.threshold is not None and worst > args.threshold:
        print(f"FAIL: {worst:.4e} exceeds the threshold {args.threshold:.4e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
