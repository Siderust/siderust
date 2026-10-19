"""Regenerate independent IAU vectors using python-sgp4==2.25's C++ backend."""
from pathlib import Path

from sgp4 import __version__
from sgp4.api import Satrec, WGS72, WGS84, accelerated

assert __version__ == "2.25" and accelerated, "Requires python-sgp4 2.25's C++ extension"
root = Path(__file__).parent
lines = (root / "vallado.txt").read_text().splitlines()
cases = [(line, lines[i + 1]) for i, line in enumerate(lines) if line.startswith("1 ")]

with (root / "iau.txt").open("w") as output:
    output.write("# python-sgp4 2.25 Vallado C++ port: model minutes x y z vx vy vz\n")
    for index in [0, 1, 3, 16, 18, 23]:
        line1, line2 = cases[index]
        output.write(line1 + "\n" + line2 + "\n")
        for model, geopotential in [("72", WGS72), ("84", WGS84)]:
            satellite = Satrec.twoline2rv(line1, line2, geopotential)
            for minutes in [-1440.0, 0.0, 1440.0, 4320.0]:
                error, position, velocity = satellite.sgp4_tsince(minutes)
                assert not error, (index, model, minutes, error)
                output.write(" ".join(map(str, [model, minutes, *position, *velocity])) + "\n")

# Synthetic cases isolate piecewise eccentricity polynomials, polar/Lyddane
# boundaries, circular orbits, retrograde singularity guards and drag signs.
from math import pi

parameters = [
    (e, inclination, 2.0, 0.0)
    for e in [0.5, 0.64, 0.65, 0.65001, 0.69999, 0.7, 0.715, 0.71501, 0.8]
    for inclination in [2.9, 63.4, 179.0]
]
parameters += [
    (0.0, 0.0, 1.0, 0.0), (0.00005, 180.0, 14.0, 0.001),
    (0.00015, 90.0, 15.0, -0.001), (0.01, 180.0, 1.0, 0.0),
    (0.01, 3.1, 1.0, 0.0), (0.01, 11.46, 1.0, 0.0),
    (0.01, 90.0, 6.3999, 0.0), (0.01, 90.0, 6.4001, 0.0),
]
with (root / "synthetic.txt").open("w") as output:
    output.write("# Vallado C++: model e inclination rev/day bstar minutes x y z vx vy vz\n")
    for model, geopotential, mode in [("72a", WGS72, "a"), ("72i", WGS72, "i"), ("84i", WGS84, "i")]:
        for eccentricity, inclination, motion, bstar in parameters:
            satellite = Satrec()
            satellite.sgp4init(geopotential, mode, 99999, 25568.0, bstar, 0.0, 0.0,
                               eccentricity, 41.0 * pi / 180.0,
                               inclination * pi / 180.0, 17.0 * pi / 180.0,
                               motion * pi / 720.0, -14.0 * pi / 180.0)
            for minutes in [-1440.123456789, 0.0, 10080.0]:
                error, position, velocity = satellite.sgp4_tsince(minutes)
                assert error in (0, 6), (model, eccentricity, inclination, minutes, error)
                output.write(" ".join(map(str, [model, eccentricity, inclination, motion,
                                               bstar, minutes, *position, *velocity])) + "\n")
