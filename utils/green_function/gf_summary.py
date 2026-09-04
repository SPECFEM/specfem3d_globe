#!/usr/bin/env python
"""One overview figure from several gf_compare.py JSON files.

Each JSON is one case (an example and a source type). The figure has three
panels -- relative L2 (log scale), amplitude ratio, and cross-correlation
lag -- with the stations along the horizontal axis, a marker per component,
and a colour per case, so the whole gate is one picture: a station that is
wrong stands out in all three, an STF width error stands out in the
amplitude and spectrum numbers, a timing error in the lag.

Usage
-----
  gf_summary.py --json label=case.json [label=case.json ...] --png out.png
                [--threshold X] [--table out.md]

`label=` is optional; without it the JSON's own `source_label` is used.
"""

import argparse
import json
from pathlib import Path

import numpy as np

COMPONENTS = ("N", "E", "Z")
MARKERS = {"N": "o", "E": "s", "Z": "^"}


def load(spec):
    if "=" in spec:
        label, path = spec.split("=", 1)
    else:
        label, path = None, spec
    with open(path) as f:
        d = json.load(f)
    if label is None:
        label = d.get("source_label", Path(path).stem)
    return label, d


def rows_of(label, d):
    rows = []
    for station, r in d["stations"].items():
        for c in COMPONENTS:
            m = r.get(c, {})
            if m.get("rel_l2") is None:
                continue
            rows.append({
                "case": label, "station": station, "comp": c,
                "distance_deg": r["plan"].get("distance_deg", float("nan")),
                "rel_l2": m["rel_l2"], "amp_ratio": m["amp_ratio"],
                "lag_seconds": m["lag_seconds"],
                "spectrum_ratio_rms_log": m.get("spectrum_ratio_rms_log", float("nan")),
            })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", nargs="+", required=True, help="[label=]path.json, one per case")
    ap.add_argument("--png", required=True)
    ap.add_argument("--threshold", type=float, help="draw the relative-L2 gate")
    ap.add_argument("--table", help="also write the numbers as a Markdown table")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cases = [load(s) for s in args.json]
    rows = [r for label, d in cases for r in rows_of(label, d)]
    if not rows:
        raise SystemExit("no station results in the given JSON files")

    # stations by distance, cases in the order given
    stations = sorted({r["station"] for r in rows},
                      key=lambda s: min(r["distance_deg"] for r in rows if r["station"] == s))
    labels = [label for label, _ in cases]
    xpos = {s: i for i, s in enumerate(stations)}
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    ncase = len(labels)
    width = 0.8 / max(1, ncase * 3)

    fig, axes = plt.subplots(3, 1, figsize=(max(7, 1.6 * len(stations) + 4), 9.5), sharex=True)
    panels = (("rel_l2", "relative L2  ‖gf3d − fwd‖ / ‖fwd‖", True),
              ("amp_ratio", "amplitude ratio  ⟨g,f⟩/⟨f,f⟩", False),
              ("lag_seconds", "cross-correlation lag [s]  (+ = gf3d later)", False))
    for ax, (key, title, log) in zip(axes, panels):
        for k, label in enumerate(labels):
            for j, c in enumerate(COMPONENTS):
                pts = [(xpos[r["station"]] + (k * 3 + j - (ncase * 3 - 1) / 2) * width, r[key])
                       for r in rows if r["case"] == label and r["comp"] == c]
                if not pts:
                    continue
                xs, ys = zip(*pts)
                ax.scatter(xs, ys, marker=MARKERS[c], color=colours[k % len(colours)], s=34,
                           edgecolor="k", linewidth=0.4, label=f"{label} {c}" if ax is axes[0] else None)
        if log:
            ax.set_yscale("log")
            if args.threshold is not None:
                ax.axhline(args.threshold, color="r", lw=0.8, ls="--", label=f"gate {args.threshold:.0e}")
        elif key == "amp_ratio":
            ax.axhline(1.0, color="0.5", lw=0.6)
        else:
            ax.axhline(0.0, color="0.5", lw=0.6)
        ax.set_title(title, fontsize=9, loc="left")
        ax.grid(True, axis="y", alpha=0.3)
    axes[0].legend(fontsize=7, ncol=max(1, ncase), loc="upper right")
    axes[-1].set_xticks(range(len(stations)))
    axes[-1].set_xticklabels([f"{s}\n{min(r['distance_deg'] for r in rows if r['station'] == s):.1f}°"
                              for s in stations], fontsize=8)
    fig.suptitle("gf3d against the forward runs", fontsize=11)
    fig.tight_layout()
    Path(args.png).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.png, dpi=130)
    plt.close(fig)
    print(f"wrote {args.png}")

    if args.table:
        lines = ["| case | station | dist | comp | rel L2 | amp ratio | lag [s] | spec rms ln |",
                 "|---|---|---|---|---|---|---|---|"]
        for r in rows:
            lines.append(f"| {r['case']} | {r['station']} | {r['distance_deg']:.1f}° | {r['comp']} | "
                         f"{r['rel_l2']:.2e} | {r['amp_ratio']:.4f} | {r['lag_seconds']:+.3f} | "
                         f"{r['spectrum_ratio_rms_log']:.2e} |")
        Path(args.table).write_text("\n".join(lines) + "\n")
        print(f"wrote {args.table}")

    worst = max(r["rel_l2"] for r in rows)
    print(f"worst relative L2: {worst:.3e}")
    if args.threshold is not None and worst > args.threshold:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
