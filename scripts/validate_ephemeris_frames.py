#!/usr/bin/env python3
"""Independent IMCCE + ERFA/CSPICE references (no Siderust matrix imports).

uv run --with pyerfa==2.0.1.5 --with spiceypy==8.3.0 python \
  scripts/validate_ephemeris_frames.py --kernel /path/to/de440s.bsp

--native-csv accepts the output of the pre-fix frame_consistency example and
prints before/frame-only/after residuals. --canonical-csv checks the corrected
example against the independently rotated pre-fix output. Kernels stay external.
"""
import argparse
import json
from pathlib import Path
import re
from urllib.request import urlopen

import erfa
import numpy as np
import spiceypy as spice

BASE = "https://ftp.imcce.fr/pub/ephem/planets/vsop87/"
# Published native-to-FK5 ties: vsop87.doc and ELP82B notice §8, p.12.
VSOP = np.array([[1, .000000440360, -.000000190919],
                 [-.000000479966, .917482137087, -.397776982902],
                 [0, .397776982902, .917482137087]])
ELP = np.array([[1, .000000437913, -.000000189859],
                [-.000000477299, .917482137607, -.397776981701],
                [0, .397776981701, .917482137607]])
ECM = erfa.ecm06(2451545., 0.)
R5H, _ = erfa.fk5hip()  # Fixed J2000 orientation; no catalogue spin.
RV = ECM @ R5H @ VSOP
RL = ECM @ R5H @ ELP


def et_from_tt(jd):
    """NAIF DELTET TT->ET convention, matching a standard generic LSK."""
    tt = (jd - 2451545.) * 86400.
    et = tt
    for _ in range(5):
        m = 6.239996 + 1.99096871e-7 * et
        et = tt + .001657 * np.sin(m + .01671 * np.sin(m))
    return et


def angle(a, b):
    return np.arctan2(np.linalg.norm(np.cross(a, b)), np.dot(a, b)) / erfa.DAS2R


def series_endpoints():
    """Direct scalar evaluation of original coefficients at nominal span edges."""
    rows = []
    for body, suffix, years in [('MERCURY', 'mer', 4000), ('EARTH', 'ear', 4000),
                                 ('MARS', 'mar', 4000), ('NEPTUNE', 'nep', 6000)]:
        for version in ['A', 'E']:
            text = urlopen(BASE + f'VSOP87{version}.{suffix}', timeout=30).read().decode('ascii')
            series = {}
            for line in text.splitlines():
                if 'VARIABLE' in line:
                    match = re.search(r'VARIABLE (\d).*\*T\*\*(\d)', line)
                    key = tuple(map(int, match.groups()))
                    series[key] = []
                elif line.strip():
                    series[key].append([float(line[79:97]), float(line[97:111]), float(line[111:131])])
            for sign in [-1, 1]:
                t = sign * years / 1000.
                p, v = np.zeros(3), np.zeros(3)
                for (axis, power), terms in series.items():
                    a, b, c = np.array(terms).T
                    arg = b + c*t
                    value = np.sum(a*np.cos(arg))
                    deriv = np.sum(-a*c*np.sin(arg))
                    p[axis-1] += t**power * value
                    v[axis-1] += (t**power * deriv + (power*t**(power-1)*value if power else 0)) / 365250.
                rows.append(dict(version=version, body=body, jd_tdb=2451545.+sign*years*365.25,
                                 position=(RV@p).tolist(), velocity=(RV@v).tolist()))
    return rows


def fixtures(kernel):
    checks = urlopen(BASE + 'vsop87.chk', timeout=30).read().decode('ascii')
    rows = []
    for version, body, jd, xyz, vel in re.findall(
        r'VSOP87([AE])\s+(\w+)\s+JD([\d.]+)[^\n]*\n([^\n]+)\n([^\n]+)', checks
    ):
        if body not in ['MERCURY', 'EARTH', 'MARS', 'NEPTUNE'] or float(jd) not in [2451545., 2415020., 2122820.]:
            continue
        p = [float(x) for x in re.findall(r"[xyz]\s+(-?\.?\d+\.\d+|-?\.\d+)", xyz)]
        v = [float(x) for x in re.findall(r"[xyz]'\s+(-?\.?\d+\.\d+|-?\.\d+)", vel)]
        assert len(p) == len(v) == 3
        rows.append(dict(version=version, body=body, jd_tdb=float(jd),
                         position=(RV @ p).tolist(), velocity=(RV @ v).tolist()))
    # ELP82B notice §9 Table H, PREC=0, native km (signs checked on rendered page).
    moon = [(2469000.5, [-361602.98536, 44996.99510, -30696.65316]),
            (2449000.5, [-363132.34248, 35863.65378, -33196.00409]),
            (2429000.5, [-371577.58161, 75271.14315, -32227.94618]),
            (2409000.5, [-373896.15893, 127406.79129, -30037.79225]),
            (2389000.5, [-346331.77361, 206365.40364, -28502.11732])]
    lunar = [dict(jd_tdb=jd, position=(RL @ p).tolist()) for jd, p in moon]
    spice.furnsh(str(kernel))
    jpl = []
    for jd in [2415020., 2451545., 2461344.5, 2488070.]:
        for body, center in [(10, 0), (399, 0), (399, 10), (301, 399)]:
            state, _ = spice.spkez(body, et_from_tt(jd), 'J2000', 'NONE', center)
            jpl.append(dict(jd_tt=jd, body=body, center=center,
                            position=(ECM @ state[:3]).tolist(),
                            velocity=(ECM @ state[3:] * 86400.).tolist()))
    return dict(provenance=dict(vsop=BASE+'vsop87.chk', elp='https://ftp.imcce.fr/pub/ephem/moon/elp82b/elp82b.ps',
                                erfa=erfa.__version__, cspice=spice.tkvrsn('TOOLKIT'), kernel=kernel.name),
                vsop=rows, vsop_endpoints=series_endpoints(), elp=lunar, jpl=jpl)


def compare(path, kernel, canonical):
    spice.furnsh(str(kernel))
    rows = np.loadtxt(path, delimiter=',')
    corrected = np.loadtxt(canonical, delimiter=',') if canonical else None
    print('JD,body,before_arcsec,before_km,frame_arcsec,frame_km,after_arcsec,after_km,before_m/s,frame_m/s,after_m/s,TT_TDB_effect_km')
    for i, row in enumerate(rows):
        jd, body, *state = row
        body = int(body)
        p, v = np.array(state[:3]), np.array(state[3:])
        r = RL if body == 301 else RV
        ref, _ = spice.spkez(body, et_from_tt(jd), 'J2000', 'NONE', 399 if body == 301 else 10)
        rp, rv = ECM @ ref[:3], ECM @ ref[3:] * 86400.
        cp, cv = r @ p, r @ v
        if body == 301:
            v[:] = np.nan
            cv[:] = np.nan
        if corrected is not None:
            # Lunar TT->TDB correction is separately allowed (at most 2 ms).
            tolerance = .003 if body == 301 else 1e-5
            assert np.linalg.norm(corrected[i, 2:5] - cp) < tolerance
            if body != 301:
                assert np.linalg.norm(corrected[i, 5:] - cv) / 86.4 < 1e-8
        dt = erfa.dtdb(2451545., jd-2451545., 0., 0., 0., 0.)
        values = [angle(p, rp), np.linalg.norm(p-rp), angle(p, cp), np.linalg.norm(p-cp),
                  angle(cp, rp), np.linalg.norm(cp-rp), np.linalg.norm(v-rv)/86.4,
                  np.linalg.norm(v-cv)/86.4, np.linalg.norm(cv-rv)/86.4,
                  np.linalg.norm(rv)*abs(dt)/86400.]
        print(f'{jd},{body},' + ','.join(f'{x:.9g}' for x in values))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--kernel', type=Path, required=True)
    parser.add_argument('--native-csv', type=Path)
    parser.add_argument('--canonical-csv', type=Path)
    parser.add_argument('--fixtures', type=Path)
    args = parser.parse_args()
    if args.fixtures:
        args.fixtures.write_text(json.dumps(fixtures(args.kernel), indent=2) + '\n')
    if args.native_csv:
        compare(args.native_csv, args.kernel, args.canonical_csv)
