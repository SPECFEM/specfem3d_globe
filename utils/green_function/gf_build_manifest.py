#!/usr/bin/env python3
"""
Build manifest and centroids files from a completed Green function database.

Walks {GF_DATABASE_PATH}/elements/*/coordinates.h5 to extract centroids,
writes centroids.bin (binary) and manifest.csv (human-readable).
Also summarizes station completion status from {GF_DATABASE_PATH}/stations/*.h5.

Usage:
    python gf_build_manifest.py /path/to/gf_database/

Output files (written to GF_DATABASE_PATH):
    centroids.bin  — binary, 32 bytes/record: uint64 morton + 3x float64 (cx,cy,cz)
    manifest.csv   — text: morton_hex, cx, cy, cz
    completion.txt — text summary of per-station completion status
"""

import argparse
import struct
import sys
from pathlib import Path

try:
    import h5py
    import numpy as np
except ImportError:
    print("Error: h5py and numpy are required. Install with: pip install h5py numpy")
    sys.exit(1)


def morton_hex_to_int(hex_str):
    """Convert a 16-character hex string to a 64-bit integer."""
    return int(hex_str, 16)


def build_manifest(db_path):
    """Build centroids.bin and manifest.csv from coordinates.h5 files."""
    elements_dir = db_path / "elements"
    if not elements_dir.is_dir():
        print(f"Error: {elements_dir} does not exist")
        return False

    # find all coordinates.h5 files
    coord_files = sorted(elements_dir.glob("*/coordinates.h5"))
    if not coord_files:
        print(f"Warning: no coordinates.h5 files found in {elements_dir}")
        return False

    print(f"Found {len(coord_files)} element coordinate files")

    records = []
    for coord_file in coord_files:
        morton_hex = coord_file.parent.name

        with h5py.File(coord_file, "r") as f:
            # read centroid from attributes
            cx = float(np.asarray(f.attrs["cx"]).flat[0])
            cy = float(np.asarray(f.attrs["cy"]).flat[0])
            cz = float(np.asarray(f.attrs["cz"]).flat[0])

        morton_int = morton_hex_to_int(morton_hex)
        records.append((morton_int, morton_hex, cx, cy, cz))

    # sort by Morton code
    records.sort(key=lambda r: r[0])

    # write centroids.bin
    bin_path = db_path / "centroids.bin"
    with open(bin_path, "wb") as f:
        for morton_int, _, cx, cy, cz in records:
            f.write(struct.pack("<Qddd", morton_int, cx, cy, cz))
    print(f"Written {bin_path} ({len(records)} records, {len(records)*32} bytes)")

    # write manifest.csv
    csv_path = db_path / "manifest.csv"
    with open(csv_path, "w") as f:
        f.write("morton_hex,cx,cy,cz\n")
        for _, morton_hex, cx, cy, cz in records:
            # 17 significant digits, so that manifest.csv round-trips a float64
            # exactly. At .15e it does not, and a reader that falls back to the
            # text manifest gets centroids that differ from centroids.bin in the
            # last bit.
            f.write(f"{morton_hex},{cx:.16e},{cy:.16e},{cz:.16e}\n")
    print(f"Written {csv_path}")

    return True


def check_element_completion(db_path):
    """Summarize per-element completion status from element station files."""
    elements_dir = db_path / "elements"
    if not elements_dir.is_dir():
        print(f"No elements directory found at {elements_dir}")
        return

    # find all station h5 files (excluding coordinates.h5)
    station_files = sorted(
        f for f in elements_dir.glob("*/*.h5") if f.name != "coordinates.h5"
    )
    if not station_files:
        print("No element station files found")
        return

    # group by station name
    station_elements = {}  # station_name -> list of (morton_hex, comp_N, comp_E, comp_Z, comp_ALL)
    for sf in station_files:
        station_name = sf.stem  # e.g., "IU.SJG"
        morton_hex = sf.parent.name
        try:
            with h5py.File(sf, "r") as f:
                comp_n = int(np.asarray(f.attrs.get("computed_N", 0)).flat[0])
                comp_e = int(np.asarray(f.attrs.get("computed_E", 0)).flat[0])
                comp_z = int(np.asarray(f.attrs.get("computed_Z", 0)).flat[0])
                comp_all = int(np.asarray(f.attrs.get("computed_ALL", 0)).flat[0])
        except Exception as e:
            print(f"  {morton_hex}/{station_name}  Error reading: {e}")
            continue

        if station_name not in station_elements:
            station_elements[station_name] = []
        station_elements[station_name].append((morton_hex, comp_n, comp_e, comp_z, comp_all))

    # summarize per station
    completion_lines = []
    print(f"\nElement completion status ({len(station_files)} element-station files, "
          f"{len(station_elements)} stations):")
    print(f"{'Station':<20s} {'Elements':>8s} {'Complete':>8s} {'N':>5s} {'E':>5s} {'Z':>5s}")
    print("-" * 55)

    for station_name in sorted(station_elements):
        elems = station_elements[station_name]
        n_total = len(elems)
        n_all = sum(1 for _, _, _, _, a in elems if a == 1)
        n_n = sum(1 for _, n, _, _, _ in elems if n == 1)
        n_e = sum(1 for _, _, e, _, _ in elems if e == 1)
        n_z = sum(1 for _, _, _, z, _ in elems if z == 1)

        status_line = (f"{station_name:<20s} {n_total:>8d} {n_all:>8d} "
                       f"{n_n:>5d} {n_e:>5d} {n_z:>5d}")
        print(f"  {status_line}")
        completion_lines.append(status_line)

    # write completion summary
    completion_path = db_path / "completion.txt"
    with open(completion_path, "w") as f:
        f.write(f"Element completion status ({len(station_files)} element-station files, "
                f"{len(station_elements)} stations)\n")
        f.write(f"{'Station':<20s} {'Elements':>8s} {'Complete':>8s} {'N':>5s} {'E':>5s} {'Z':>5s}\n")
        f.write("-" * 55 + "\n")
        for line in completion_lines:
            f.write(line + "\n")
    print(f"Written {completion_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Build manifest and centroids files from a Green function database"
    )
    parser.add_argument(
        "db_path",
        type=Path,
        help="Path to the Green function database directory",
    )
    args = parser.parse_args()

    if not args.db_path.is_dir():
        print(f"Error: {args.db_path} is not a directory")
        sys.exit(1)

    build_manifest(args.db_path)
    check_element_completion(args.db_path)


if __name__ == "__main__":
    main()
