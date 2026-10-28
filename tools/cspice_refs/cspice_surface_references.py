#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 Vallés Puig, Ramon
"""Regenerate offline references for siderust surface Alt/Az tests.

Development-only helper. Not required to build or test siderust, and not a
runtime dependency.

## What CSPICE computes vs what this script computes independently

CSPICE (via SpiceyPy) is used only for:

* leap-second / timescale conversion:
  ``ET = unitim((jd_tt - 2451545)*86400, "TDT", "TDB")``
* lunar orientation:
  ``M = pxform("IAU_MOON", "J2000", ET)`` with ``v_j2000 = M @ v_body``
* Mars planetographic rectangular conversion for longitude interop checks:
  ``pgrrec("MARS", lon_west, lat, 0, a, f)``

Independent NumPy / trigonometry (not CSPICE) is used for:

* lunar sphere site Cartesian position
* local ENU basis on that sphere
* finite target offset along ENU
* geometric altitude, azimuth and slant range from the body-fixed LOS

Do **not** describe the stored Alt/Az/range literals as “full CSPICE end-to-end”
results. They are a composite: CSPICE orientation (+ optional pgrrec) with
independent ENU geometry.

Dependencies (optional):
  uv pip install spiceypy numpy

Kernels (NAIF generic; download once):
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc
  https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls

Usage:
  export CSPICE_PCK=/path/to/pck00011.tpc
  export CSPICE_LSK=/path/to/naif0012.tls
  python3 tools/cspice_refs/cspice_surface_references.py

Prints Rust-ready literals for tests/body_surface_horizontal.rs.
Validated against CSPICE N0067 / SpiceyPy 8.3.0.
"""

from __future__ import annotations

import math
import os
import sys


def require_spice():
    try:
        import numpy as np
        import spiceypy as spice
    except ImportError as exc:  # pragma: no cover
        raise SystemExit(
            "spiceypy and numpy are required to regenerate references:\n"
            "  uv pip install spiceypy numpy\n"
            f"({exc})"
        ) from exc
    return np, spice


def enu_basis_sphere(lon_rad: float, lat_rad: float, np):
    """Independent sphere ENU (east-positive longitude). Not a CSPICE call."""
    sin_lon, cos_lon = math.sin(lon_rad), math.cos(lon_rad)
    sin_lat, cos_lat = math.sin(lat_rad), math.cos(lat_rad)
    east = np.array([-sin_lon, cos_lon, 0.0])
    east /= np.linalg.norm(east)
    up = np.array([cos_lat * cos_lon, cos_lat * sin_lon, sin_lat])
    up /= np.linalg.norm(up)
    north = np.cross(up, east)
    north /= np.linalg.norm(north)
    return east, north, up


def altaz(los, east, north, up, np):
    """Independent geometric Alt/Az/range. Not CSPICE azlcpo."""
    e = float(np.dot(los, east))
    n = float(np.dot(los, north))
    u = float(np.dot(los, up))
    horiz = math.hypot(n, e)
    alt = math.degrees(math.atan2(u, horiz))
    # Match siderust's documented singularity: exact vertical → az = 0.
    az = 0.0 if horiz == 0.0 else math.degrees(math.atan2(e, n)) % 360.0
    return alt, az, float(np.linalg.norm(los))


def main() -> int:
    np, spice = require_spice()
    pck = os.environ.get("CSPICE_PCK", "pck00011.tpc")
    lsk = os.environ.get("CSPICE_LSK", "naif0012.tls")
    for path in (pck, lsk):
        if not os.path.isfile(path):
            print(f"missing kernel: {path}", file=sys.stderr)
            return 1

    spice.furnsh(pck)
    spice.furnsh(lsk)

    jd_tt = 2451545.0
    et = spice.unitim((jd_tt - 2451545.0) * 86400.0, "TDT", "TDB")
    # CSPICE: body-fixed → J2000. Columns of M are body axes in J2000.
    r_moon = spice.pxform("IAU_MOON", "J2000", et)

    print("// --- CSPICE orientation ---")
    print(f"// jd_tt={jd_tt}, ET via unitim(TDT→TDB)")
    print(f"// kernels={os.path.basename(pck)}+{os.path.basename(lsk)}")
    print("// M = pxform(IAU_MOON, J2000, ET); v_j2000 = M @ v_body")
    for i, axis in enumerate(["X", "Y", "Z"]):
        print(f"// body {axis} in J2000: {[float(r_moon[j, i]) for j in range(3)]}")
    print(f"// M row-major: {r_moon.tolist()}")

    r_km = 1737.4
    lon, lat = math.radians(20.0), math.radians(-3.5)
    obs = np.array(
        [
            r_km * math.cos(lat) * math.cos(lon),
            r_km * math.cos(lat) * math.sin(lon),
            r_km * math.sin(lat),
        ]
    )
    east, north, up = enu_basis_sphere(lon, lat, np)
    target = obs + 50.0 * east + 10.0 * north + 5.0 * up
    alt, az, rng = altaz(target - obs, east, north, up, np)
    print("\n// --- Independent sphere ENU / Alt/Az (not CSPICE) ---")
    print(f"// site lon=20°E lat=-3.5° h=0 R={r_km} km")
    print(f"// target_fixed_km = {target.tolist()}")
    print(f"// expected alt={alt}, az={az}, range_km={rng}")

    a, b = 3396.19, 3376.20
    f = (a - b) / a
    lon_west = math.radians(137.4)
    lat_pgr = math.radians(-4.6)
    # CSPICE planetographic → rectangular (default Mars west-positive lon).
    rec = spice.pgrrec("MARS", lon_west, lat_pgr, 0.0, a, f)
    print("\n// --- CSPICE pgrrec Mars longitude interop ---")
    print("// default BODY499 planetographic lon is west-positive")
    print("// (kernel-pool BODY*_PGR_POSITIVE_LON overrides are NOT applied here)")
    print(f"// pgrrec(west=137.4°, lat=-4.6°) km = {rec.tolist()}")
    print("// Siderust east-positive longitude = -137.4")

    spice.kclear()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
