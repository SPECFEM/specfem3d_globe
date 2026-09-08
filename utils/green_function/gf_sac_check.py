#!/usr/bin/env python
"""Check xgf3d's SAC output against its own ASCII output and against the
forward run's SAC headers.

Run xgf3d with ``--format all`` (and optionally ``--partials N``) into one
directory, then::

  gf_sac_check.py --dir <that directory> --fwd <forward_cmt/OUTPUT_FILES>

What is asserted, per station and component:

* obspy reads every file, seismograms and partials;
* ``delta`` is the stored spacing and ``npts`` the output length, both
  from the ASCII header (the forward run is on the solver grid, so neither
  is comparable with it), ``b`` is the ASCII header's ``t_first`` and ``o``
  is zero;
* the SAC data equal the ASCII columns rounded to single precision, which
  is what the writer stores, exactly;
* against the forward SAC file: the reference time ``nz*`` (the PDE time
  plus the CMT time shift, with the solver's rollover), ``evla/evlo/evdp``,
  ``stla/stlo``, ``cmpaz/cmpinc`` and ``o`` are equal, and ``b`` lies within
  one stored sample below the forward run's ``-t0`` (the planned axis
  starts at or before the requested ``t0`` by whole stored samples);
* a partial's file carries the same header as the seismogram apart from
  ``kuser1``/``kuser2``, which name the parameter and its unit, and its
  data equal the partials file's column.

The forward SAC files and the extraction's own ASCII are the oracles;
nothing here re-derives a header rule.
"""

import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gf_compare import read_gf3d  # noqa: E402

try:
    from obspy import read as obspy_read
except ImportError:  # pragma: no cover - the runner skips before this
    raise SystemExit("gf_sac_check.py needs obspy")

COMPONENTS = "NEZ"


class Checker:
    def __init__(self):
        self.nbad = 0
        self.nok = 0

    def check(self, cond, what):
        if cond:
            self.nok += 1
        else:
            self.nbad += 1
            print(f"  FAIL {what}")


def read_partials(path):
    """The '# partials' line and the columns of NET.STA.partials.txt."""
    names = None
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            body = line[1:].strip()
            if body.startswith("partials"):
                fields = body.partition(":")[2].split()
                names = fields[1:]
                assert int(fields[0]) == len(names)
    data = np.loadtxt(path, comments="#")
    ndp = len(names)
    assert data.shape[1] == 1 + 3 * ndp, (path, data.shape)
    cols = {}
    for ic, comp in enumerate(COMPONENTS):
        for ip, name in enumerate(names):
            cols[(comp, name)] = data[:, 1 + ndp * ic + ip]
    return names, data[:, 0], cols


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="xgf3d output written with --format all")
    ap.add_argument("--fwd", required=True, help="the forward run's OUTPUT_FILES directory")
    args = ap.parse_args()

    chk = Checker()
    txts = sorted(glob.glob(os.path.join(args.dir, "*.gf3d.txt")))
    if not txts:
        raise SystemExit(f"no *.gf3d.txt in {args.dir}; run xgf3d --seis ... --format all")

    for txt in txts:
        station, plan, t, traces = read_gf3d(txt)
        dt_sub = plan["dt"] * plan["subsample_step"]
        ptxt = txt.replace(".gf3d.txt", ".partials.txt")
        names, tp, pcols = ([], None, {})
        if os.path.exists(ptxt):
            names, tp, pcols = read_partials(ptxt)
            chk.check(np.array_equal(tp, t), f"{station}: partials time axis equals the seismogram's")
        print(f"{station}: nt = {plan['nt']}, dt_sub = {dt_sub}, t_first = {plan['t_first']}, "
              f"{len(names)} partials")

        for comp in COMPONENTS:
            chan = f"BX{comp}"
            sac = os.path.join(args.dir, f"{station}.{chan}.sem.sac")
            try:
                tr = obspy_read(sac)[0]
            except Exception as exc:  # noqa: BLE001
                chk.check(False, f"{station}.{chan}: obspy read failed: {exc}")
                continue
            h = tr.stats.sac

            # the axis, from the ASCII header
            chk.check(abs(h.delta - dt_sub) <= 1e-6 * dt_sub, f"{station}.{chan}: delta = dt_sub")
            chk.check(h.npts == plan["nt"] and tr.stats.npts == plan["nt"], f"{station}.{chan}: npts = nt")
            chk.check(abs(h.b - plan["t_first"]) <= 1e-5 * abs(plan["t_first"]),
                      f"{station}.{chan}: b = t_first ({h.b} vs {plan['t_first']})")
            chk.check(h.o == 0.0, f"{station}.{chan}: o = 0")

            # the data, single precision, exactly
            col = traces[comp].astype(np.float32)
            chk.check(tr.data.dtype == np.float32 and np.array_equal(tr.data, col),
                      f"{station}.{chan}: data equal the ASCII column in single precision")

            # the forward run's headers
            fsac = os.path.join(args.fwd, f"{station}.{chan}.sem.sac")
            if os.path.exists(fsac):
                fh = obspy_read(fsac)[0].stats.sac
                for key in ("nzyear", "nzjday", "nzhour", "nzmin", "nzsec", "nzmsec"):
                    chk.check(getattr(h, key) == getattr(fh, key),
                              f"{station}.{chan}: {key} = {getattr(h, key)} vs forward {getattr(fh, key)}")
                for key in ("evla", "evlo", "evdp", "stla", "stlo", "cmpaz", "cmpinc", "o"):
                    chk.check(np.float32(getattr(h, key)) == np.float32(getattr(fh, key)),
                              f"{station}.{chan}: {key} = {getattr(h, key)} vs forward {getattr(fh, key)}")
                # ours starts at or before the forward run's -t0, by less
                # than one stored sample
                chk.check(fh.b - h.delta < h.b <= fh.b + 1e-4,
                          f"{station}.{chan}: b in (fwd b - delta, fwd b] ({h.b} vs {fh.b}, delta {h.delta})")
                chk.check(str(h.kevnm).strip() == str(fh.kevnm).strip(),
                          f"{station}.{chan}: kevnm '{h.kevnm}' vs forward '{fh.kevnm}'")
                chk.check(h.kstnm.strip() == fh.kstnm.strip() and h.knetwk.strip() == fh.knetwk.strip(),
                          f"{station}.{chan}: kstnm/knetwk vs forward")
            else:
                print(f"  (no forward file {fsac}; header comparison skipped)")

            # the partials
            for name in names:
                psac = os.path.join(args.dir, f"{station}.{chan}.{name}.sem.sac")
                try:
                    ptr = obspy_read(psac)[0]
                except Exception as exc:  # noqa: BLE001
                    chk.check(False, f"{station}.{chan}.{name}: obspy read failed: {exc}")
                    continue
                ph = ptr.stats.sac
                chk.check(np.array_equal(ptr.data, pcols[(comp, name)].astype(np.float32)),
                          f"{station}.{chan}.{name}: data equal the partials column")
                chk.check(str(ph.kuser1).strip() == name, f"{station}.{chan}.{name}: kuser1 = {name}")
                same = all(getattr(ph, k) == getattr(h, k) for k in
                           ("b", "delta", "npts", "o", "nzyear", "nzjday", "nzhour", "nzmin", "nzsec",
                            "nzmsec", "evla", "evlo", "evdp", "stla", "stlo", "cmpaz", "cmpinc"))
                chk.check(same, f"{station}.{chan}.{name}: header equals the seismogram's")

    print(f"\n{chk.nok} checks passed, {chk.nbad} failed")
    return 1 if chk.nbad else 0


if __name__ == "__main__":
    sys.exit(main())
