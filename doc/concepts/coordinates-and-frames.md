# Coordinates, centers, and reference frames

A position gives the displacement from a chosen origin along chosen axes.
Siderust's Cartesian and spherical position types therefore carry a center,
a frame, and a length unit. The same physical point can have barycentric
(Solar System barycenter), heliocentric (Sun), and geocentric (Earth)
positions, all expressed along the same axes. Planetocentric positions use
the corresponding body's origin, subject to the selected ephemeris's
[target semantics](../guides/ephemeris-backends.md).

A direction carries a frame and has no distance. Rotating a direction changes
its components, but computing parallax requires a finite target position.
Subtracting two positions in the same coordinate system yields a displacement.
`qtty` represents physical units and supports explicit conversion. Prefer
the `siderust::qtty` re-export when passing quantities to Siderust.

The center type identifies an origin family. For static centers such as
`Geocentric`, the type identifies the origin completely. For `Topocentric`
the origin also depends on the observer's Earth site. `SurfaceObserver<B>`
carries site parameters for a specific body, and `Bodycentric` carries
runtime body parameters. The type alone cannot establish that two such
positions refer to the same origin.

## Celestial and orbital axes

| Frame | Meaning in Siderust |
| --- | --- |
| `ICRS`, `ICRF` | The International Celestial Reference System and its observational realization. The rotation providers treat their axes as aligned. |
| `GCRS` | Geocentric Celestial Reference System. The orientation provider aligns its axes with ICRS. A full relativistic barycentric-to-geocentric transformation requires additional physics. |
| `EquatorialMeanJ2000`, `EME2000` | Same mean equator and equinox axes at J2000. The fixed IAU 2006 frame bias connects these axes to ICRS. |
| `EquatorialMeanOfDate`, `EquatorialTrueOfDate` | Mean axes use precession to the epoch. True axes also include the selected nutation model. |
| `EclipticMeanJ2000` | IAU 2006 mean ecliptic and equinox at J2000, matching ERFA `ecm06` at that epoch. Native VSOP87 and ELP axes and SPICE `ECLIPJ2000` use different conventions. |
| `EclipticMeanOfDate`, `EclipticTrueOfDate` | Epoch-dependent ecliptic axes with mean or true equinox conventions. |
| `TEME` | True Equator, Mean Equinox, used by SGP4/SDP4. The provider connects TEME to the true equator and equinox of date with the approximate equation-of-equinoxes rotation Δψ cos εₐ. |
| `FK4B1950`, `Galactic` | Fixed catalog rotations. The FK4 path is a constant frame tie and omits E-term removal. The Galactic frame marker specifies axes independently of the coordinate center. |

The [ERFA ecliptic rotation](https://github.com/liberfa/erfa/blob/v2.0.1/src/ecm06.c)
and [IERS description of ICRF](https://www.iers.org/iers/en/dataproducts/icrf/icrf)
define the conventions used here. State-dependent orbit-relative frames such as
RTN, LVLH, and VNC use a concrete orbit state through
`astro::dynamics::frames`. The dynamics module also provides covariance
transport for these frames. A position rotation alone does not transform
velocity between rotating frames.

## Terrestrial and local axes

`CIRS` uses the Celestial Intermediate Pole and Celestial Intermediate Origin.
`TIRS` applies Earth rotation to those axes. The polar-motion step connects
TIRS to the terrestrial system (ITRS), whose practical realizations are ITRFs.
Siderust aligns its `ECEF` and `ITRF` Cartesian axes in the transform providers.
`Geodetic<ECEF>` uses the WGS84 ellipsoid and `Geodetic<ITRF>` uses GRS80.
The marker alone does not carry a station coordinate epoch, realization
identifier, or crustal velocity. See [time and Earth orientation](time-and-earth-orientation.md)
and the [IERS definition of ITRS](https://www.iers.org/iers/en/dataproducts/itrs/itrs).

A local horizontal frame uses east, north, and up at a site. Azimuth increases
from north through east. The site datum determines local up and the observer
position. A finite target's line of sight includes observer parallax.
Planetary body-fixed axes use the body's pole and prime meridian.
The [surface observation guide](../guides/lunar-observations.md) describes the
supported site models and lunar frame choices.
