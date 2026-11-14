"""Generate independent compact pointing fixtures using DE440s and ERFA.

Run: python compact_reference.py /path/to/de440s.bsp > compact_reference.csv
Requires numpy, pyerfa, and jplephem. Epochs and observer data are deterministic.
UT1 and polar motion below are synthetic inputs shared by both pipelines.
"""

import sys

import erfa
import numpy as np
from jplephem.spk import SPK


kernel = SPK.open(sys.argv[1])
au_km = 149597870.7
start, end = 2458849.5, 2470172.5
# Include the range boundaries, 2020 Jupiter/Saturn conjunction, 2022
# Mars opposition, and the April 2020 lunar perigee and apogee.
dates = [start, np.nextafter(end, start), 2459205.263, 2459921.75,
         2458947.26, 2458960.30, *np.linspace(start + 1, end - 1, 48)]
sites = [(0.0, 0.0, 0.0), (-0.1275, 51.507222, 20.0),
         (-70.403, -24.625, 2635.0), (18.96, 69.65, 100.0)]

print("tdb,lon,lat,height,body,gx,gy,gz,az_deg,alt_deg,vx,vy,vz,dut1,xp,yp,apparent_az,apparent_alt,tt_to_tdb_s,dpsi_2000a,deps_2000a")
for i, tdb in enumerate(dates):
    lon, lat, height = sites[i % len(sites)]
    lo, la = np.deg2rad([lon, lat])
    observer = erfa.gd2gc(1, lo, la, height) / 1000
    east = np.array([-np.sin(lo), np.cos(lo), 0])
    north = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    up = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    # Current atomic offset assumed for this test, including synthetic future epochs.
    tt_to_tdb = erfa.dtdb(tdb, 0.0, 0.0, 0.0, 0.0, 0.0)
    tt = erfa.tdbtt(tdb, 0.0, tt_to_tdb)
    nutation = erfa.nut06a(*tt)
    dut1, xp, yp = (0.4, 0.15, -0.25)
    ut1 = (tt[0], tt[1] - (69.184 - dut1) / 86400)
    rotation = erfa.c2t06a(*tt, *ut1, *np.deg2rad(np.array([xp, yp]) / 3600))
    emb, emb_v = kernel[0, 3].compute_and_differentiate(tdb)
    earth_offset, earth_offset_v = kernel[3, 399].compute_and_differentiate(tdb)
    earth = emb + earth_offset
    velocity = (emb_v + earth_offset_v) / au_km
    to_ecliptic = erfa.ecm06(2451545.0, 0.0)
    c_au_day = 299792.458 * 86400 / au_km

    def horizontal(vector):
        local = rotation @ vector - observer
        azimuth = np.arctan2(east @ local, north @ local) % (2 * np.pi)
        altitude = np.arctan2(up @ local, np.hypot(east @ local, north @ local))
        return np.rad2deg([azimuth, altitude])

    for body in (10, 301, 1, 2, 4, 5, 6, 7, 8):
        vector = (kernel[3, 301].compute(tdb) - earth_offset if body == 301
                  else kernel[0, body].compute(tdb) - earth)
        # Apply annual aberration to the geocentric direction, matching the
        # correction under test. Observer subtraction follows the correction.
        distance = np.linalg.norm(vector)
        beta = velocity / c_au_day
        direction = erfa.ab(vector / distance, beta, 1.0, np.sqrt(1 - beta @ beta))
        fields = [tdb, lon, lat, height, body, *(to_ecliptic @ vector),
                  *horizontal(vector), *(to_ecliptic @ velocity), dut1, xp, yp,
                  *horizontal(distance * direction), tt_to_tdb, *nutation]
        print(",".join(str(float(x)) for x in fields))
