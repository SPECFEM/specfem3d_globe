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
  gf_summary.py --json label=case.json [label=case.json ...] [--png out.png]
                [--table out.md] [--threshold X] [--gate gate.json]

`label=` is optional; without it the JSON's own `source_label` is used.
`--threshold` draws the gate on the figure and reports against it; the
exit status is 1 only when `--gate` is also given and the worst relative L2
exceeds the threshold. That split is deliberate: in a workflow the figures
are produced by one rule and the pass/fail decision by a later one, so a
failed gate removes only its own gate file and never the figures you need
to see why it failed.
"""

import argparse
import json
from pathlib import Path

import numpy as np

COMPONENTS = ("N", "E", "Z")
VECTOR = "vec"                      # the three components together: the gate's measure
MARKERS = {"N": "o", "E": "s", "Z": "^", VECTOR: "*"}
SIZES = {"N": 34, "E": 34, "Z": 34, VECTOR: 110}


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
    """One row per component, plus one 'vec' row per station -- the three
    components together, which is what the gate is applied to. A component
    near a radiation node is small and its own relative error inflated; the
    vector row is not, and the component rows remain for diagnosis."""
    rows = []
    nan = float("nan")
    for station, r in d["stations"].items():
        dist = r["plan"].get("distance_deg", nan)
        for c in COMPONENTS:
            m = r.get(c, {})
            if m.get("rel_l2") is None:
                continue
            rows.append({
                "case": label, "station": station, "comp": c, "distance_deg": dist,
                "rel_l2": m["rel_l2"], "amp_ratio": m["amp_ratio"],
                "lag_seconds": m["lag_seconds"],
                "spectrum_ratio_rms_log": m.get("spectrum_ratio_rms_log", nan),
            })
        v = r.get("vector", {})
        if v.get("rel_l2") is not None:
            rows.append({
                "case": label, "station": station, "comp": VECTOR, "distance_deg": dist,
                "rel_l2": v["rel_l2"], "amp_ratio": v["amp_ratio"],
                "lag_seconds": nan, "spectrum_ratio_rms_log": nan,
            })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--json", nargs="+", required=True, help="[label=]path.json, one per case")
    ap.add_argument("--png", help="the overview figure")
    ap.add_argument("--threshold", type=float, help="the relative-L2 gate: drawn, and reported against")
    ap.add_argument("--table", help="also write the numbers as a Markdown table")
    ap.add_argument("--gate", help="write the pass/fail verdict here and exit 1 on failure")
    args = ap.parse_args()

    cases = [load(s) for s in args.json]
    rows = [r for label, d in cases for r in rows_of(label, d)]
    if not rows:
        raise SystemExit("no station results in the given JSON files")

    if args.png:
        draw(args, cases, rows)

    if args.table:
        write_table(args.table, rows)

    # the gate applies to the station vectors; components are reported
    vec = [r for r in rows if r["comp"] == VECTOR] or rows
    worst = max(r["rel_l2"] for r in vec)
    worst_comp = max(r["rel_l2"] for r in rows if r["comp"] != VECTOR) if any(r["comp"] != VECTOR for r in rows) else worst
    failing = [r for r in vec if args.threshold is not None and r["rel_l2"] > args.threshold]
    print(f"worst station-vector relative L2: {worst:.3e}  (worst single component {worst_comp:.3e})" +
          (f"  gate {args.threshold:.1e}: {'FAIL' if failing else 'pass'}" if args.threshold is not None else ""))
    for r in failing:
        print(f"  above the gate: {r['case']} {r['station']}  {r['rel_l2']:.3e}")

    if args.gate:
        verdict = {"threshold": args.threshold, "worst_rel_l2": worst, "pass": not failing,
                   "failing": [{k: r[k] for k in ("case", "station", "comp", "rel_l2", "amp_ratio", "lag_seconds")}
                               for r in failing],
                   "cases": [label for label, _ in cases]}
        Path(args.gate).parent.mkdir(parents=True, exist_ok=True)
        with open(args.gate, "w") as f:
            json.dump(verdict, f, indent=2)
        print(f"wrote {args.gate}")
        if failing:
            raise SystemExit(1)


def draw(args, cases, rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # stations by distance, cases in the order given
    stations = sorted({r["station"] for r in rows},
                      key=lambda s: min(r["distance_deg"] for r in rows if r["station"] == s))
    labels = [label for label, _ in cases]
    xpos = {s: i for i, s in enumerate(stations)}
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    ncase = len(labels)
    kinds = tuple(COMPONENTS) + (VECTOR,)
    width = 0.8 / max(1, ncase * len(kinds))

    fig, axes = plt.subplots(3, 1, figsize=(max(7, 1.6 * len(stations) + 4), 9.5), sharex=True)
    panels = (("rel_l2", "relative L2  ‖gf3d − fwd‖ / ‖fwd‖   (★ = three components together, the gate's measure)", True),
              ("amp_ratio", "amplitude ratio  ⟨g,f⟩/⟨f,f⟩", False),
              ("lag_seconds", "cross-correlation lag [s]  (+ = gf3d later)", False))
    for ax, (key, title, log) in zip(axes, panels):
        for k, label in enumerate(labels):
            for j, c in enumerate(kinds):
                pts = [(xpos[r["station"]] + (k * len(kinds) + j - (ncase * len(kinds) - 1) / 2) * width, r[key])
                       for r in rows if r["case"] == label and r["comp"] == c and np.isfinite(r[key])]
                if not pts:
                    continue
                xs, ys = zip(*pts)
                ax.scatter(xs, ys, marker=MARKERS[c], color=colours[k % len(colours)], s=SIZES[c],
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


def write_table(path, rows):
    lines = ["| case | station | dist | comp | rel L2 | amp ratio | lag [s] | spec rms ln |",
             "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['case']} | {r['station']} | {r['distance_deg']:.1f}° | {r['comp']} | "
                     f"{r['rel_l2']:.2e} | {r['amp_ratio']:.4f} | {r['lag_seconds']:+.3f} | "
                     f"{r['spectrum_ratio_rms_log']:.2e} |")
    Path(path).write_text("\n".join(lines) + "\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
