#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Vallés Puig, Ramon
"""Development-only independent CSPICE fixture generator; never used by Rust builds.

See README.md for the pinned environment and kernel set.
"""

# Independent reference generation: Python 3 + spiceypy 8.3.0 (CSPICE N0067).
import argparse
import csv
import json
import os
import shutil
import struct
from itertools import pairwise
from pathlib import Path
from tempfile import TemporaryDirectory

from provenance import check_kernels, check_toolkit

FIXTURES = Path(__file__).resolve().parents[2] / "tests/fixtures/lunar_de440"

ROTATED_FRAME = """KPL/FK
\\begindata
FRAME_TEST_LUNAR_ROTATED = 139999
FRAME_139999_NAME = 'TEST_LUNAR_ROTATED'
FRAME_139999_CLASS = 4
FRAME_139999_CLASS_ID = 139999
FRAME_139999_CENTER = 301
TKFRAME_139999_SPEC = 'ANGLES'
TKFRAME_139999_RELATIVE = 'MOON_ME'
TKFRAME_139999_ANGLES = ( 30 50 20 )
TKFRAME_139999_AXES = ( 1 2 3 )
TKFRAME_139999_UNITS = 'DEGREES'
"""


def generate(source, out, hashes, spice, np):
    files = list(hashes)
    for name in files:
        spice.furnsh(str(source / name))
    # Five original BPC records. Preserve exact DE440 coefficients with PCKW02.
    b = (source / "moon_pa_de440_200625.bpc").read_bytes()
    r = struct.unpack_from("<i", b, 76)[0]
    first, last, cls, _base, _typ, begin, finish = struct.unpack_from(
        "<2d5i", b, (r - 1) * 1024 + 24
    )
    init, interval, size, _count = struct.unpack_from("<4d", b, (finish - 4) * 8)
    k = int((820000000 - init) // interval)
    start = init + k * interval
    end = start + 4 * interval
    h = spice.pckopn(str(out / "excerpt.bpc"), "DE440 original coefficients excerpt", 0)
    records = np.array(
        struct.unpack_from(
            "<" + str(5 * int(size)) + "d", b, (begin - 1 + k * int(size)) * 8
        )
    ).reshape(5, int(size))
    cdata = records[:, 2:].flatten()
    spice.pckw02(
        h,
        cls,
        "J2000",
        start,
        end,
        "DE440 five original records",
        interval,
        5,
        int((size - 2) / 3) - 1,
        cdata,
        start,
    )
    spice.pckcls(h)
    h = spice.spkopn(str(out / "excerpt.bsp"), "DE440s original segment excerpt", 0)
    for target in [10, 3, 301, 399]:
        handle, desc, ident = spice.spksfs(target, (start + end) / 2, 80)
        spice.spksub(handle, desc, ident, start, end, h)
    spice.spkcls(h)
    for f in ["moon_de440_220930.tf", "naif0012.tls"]:
        shutil.copyfile(source / f, out / f)
    # A deliberately large three-axis fixed rotation exercises composition order.
    (out / "test_rotation.tf").write_text(ROTATED_FRAME, encoding="ascii", newline="\n")
    spice.furnsh(str(out / "test_rotation.tf"))
    # References always query the ORIGINAL kernels, not the excerpts or Siderust.
    sites = [
        (0, 0, 1737.4),
        (45, 30, 1738.0),
        (0, 90, 1737.4),
        (70, -90, 1737.4),
        (32.3, -89, 1738.0),
        (25.2, -80.13, 1736.0),
        (179.999, -45, 1737.4),
        (-179.999, 45, 1737.4),
    ]
    header = (
        ["et", "lon", "lat", "radius_km", "target", "frame"]
        + [f"m{i}" for i in range(9)]
        + [f"p{i}" for i in range(3)]
        + [f"v{i}" for i in range(3)]
        + ["alt", "az", "range", "above"]
    )

    def row(et, lon, lat, radius, target, frame):
        rot = spice.pxform("J2000", frame, et)
        pos = spice.spkpos(str(target), et, frame, "NONE", "301")[0]
        site = spice.latrec(radius, np.deg2rad(lon), np.deg2rad(lat))
        rel = pos - site
        longitude = np.deg2rad(lon)
        east = np.array([-np.sin(longitude), np.cos(longitude), 0])
        up = site / radius
        north = np.cross(up, east)
        e, n, u = np.dot(rel, east), np.dot(rel, north), np.dot(rel, up)
        alt = np.rad2deg(np.arctan2(u, np.hypot(n, e)))
        az = np.rad2deg(np.arctan2(e, n)) % 360
        return [
            et,
            lon,
            lat,
            radius,
            target,
            spice.namfrm(frame),
            *rot.flatten(),
            *pos,
            *rel,
            alt,
            az,
            np.linalg.norm(rel),
            int(alt > 0),
        ]

    with (out / "reference.csv").open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for et in np.linspace(start, end, 17):
            for lon, lat, radius in sites:
                for target in [10, 399]:
                    w.writerow(row(et, lon, lat, radius, target, "MOON_ME"))
        # Explicit Sun crossings of the horizon by bisection, just above/below.
        for lon, lat, radius in [sites[0], sites[4]]:
            ts = np.linspace(start, end, 257)
            for a, z in pairwise(ts):
                if (
                    row(a, lon, lat, radius, 10, "MOON_ME")[-4]
                    * row(z, lon, lat, radius, 10, "MOON_ME")[-4]
                    >= 0
                ):
                    continue
                for _ in range(45):
                    mid = (a + z) / 2
                    if (
                        row(a, lon, lat, radius, 10, "MOON_ME")[-4]
                        * row(mid, lon, lat, radius, 10, "MOON_ME")[-4]
                        <= 0
                    ):
                        z = mid
                    else:
                        a = mid
                for et in [(a + z) / 2 - 1, (a + z) / 2 + 1]:
                    w.writerow(row(et, lon, lat, radius, 10, "MOON_ME"))
        for et in [start, end]:
            for frame in ["MOON_PA_DE440", "MOON_ME_DE440_ME421", "TEST_LUNAR_ROTATED"]:
                w.writerow(row(et, 45, 30, 1737.4, 399, frame))
    # Broad external-kernel validation, plus actual BPC coverage edges for rotations.
    with (out / "external_reference.csv").open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for et in [-4733596800.0, 0, 717768000.0, 1700000000.0, 4735368000.0]:
            for target in [10, 399]:
                w.writerow(row(et, 32.3, -89, 1738.0, target, "MOON_ME"))
    with (out / "coverage_rotations.csv").open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        for et in [first, first + 1, last - 1, last]:
            w.writerow([et, *spice.pxform("J2000", "MOON_PA_DE440", et).flatten()])
    manifest = {
        "generator": "SpiceyPy " + spice.__version__,
        "toolkit": spice.tkvrsn("TOOLKIT"),
        "abcorr": "NONE",
        "excerpt_et": [start, end],
        "source_sha256": hashes,
    }
    (out / "provenance.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="ascii", newline="\n"
    )
    # Time references separately quantify tempoch vs LSK DELTET approximations.
    with (out / "time_reference.csv").open("w", encoding="ascii", newline="") as f:
        w = csv.writer(f)
        for utc in [
            "1972-07-01T00:00:00Z",
            "1999-01-01T00:00:00Z",
            "2016-12-31T23:59:59Z",
            "2016-12-31T23:59:60Z",
            "2017-01-01T00:00:00Z",
            "2026-10-29T12:00:00Z",
        ]:
            et = spice.str2et(utc)
            w.writerow([utc, spice.unitim(et, "TDB", "TDT"), et])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kernel-dir", type=Path, default=os.environ.get("KERNEL_DIR"))
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--verify", action="store_true", help="compare against checked-in fixtures"
    )
    mode.add_argument("--output-dir", type=Path, default=FIXTURES)
    args = parser.parse_args()
    if args.kernel_dir is None:
        parser.error("provide --kernel-dir or set KERNEL_DIR")
    try:
        import numpy as np
        import spiceypy as spice
    except ImportError as error:
        parser.exit(1, f"error: install scripts/cspice/requirements.txt ({error})\n")
    try:
        hashes = check_kernels(args.kernel_dir)
        check_toolkit(spice)
        destination = FIXTURES if args.verify else args.output_dir.resolve()
        # CSPICE writes newly created kernels. Stage the complete set before
        # touching the destination so a generation failure leaves it intact.
        parent = None if args.verify else destination.parent
        with TemporaryDirectory(prefix="lunar-de440-", dir=parent) as temp:
            generated = Path(temp)
            generate(args.kernel_dir, generated, hashes, spice, np)
            if args.verify:
                changed = [
                    path.name
                    for path in sorted(generated.iterdir())
                    if not (destination / path.name).is_file()
                    or path.read_bytes() != (destination / path.name).read_bytes()
                ]
                if changed:
                    raise ValueError("fixture differences: " + ", ".join(changed))
                print(f"Verified lunar DE440 fixtures in {destination}")
            else:
                destination.mkdir(exist_ok=True)
                for path in sorted(generated.iterdir()):
                    path.replace(destination / path.name)
                print(f"Generated lunar DE440 fixtures in {destination}")
    except (OSError, ValueError, spice.utils.exceptions.SpiceyError) as error:
        parser.exit(1, f"error: {error}\n")
    finally:
        spice.kclear()


if __name__ == "__main__":
    main()
