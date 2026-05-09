# `coordinates`

## Scientific scope

This module defines the typed coordinate language used across `siderust`. A
coordinate is not just a tuple of numbers: it has a reference center, a
reference frame, and a magnitude semantics such as direction, position, or
velocity. The scientific purpose of the module is to preserve those semantics
in the type system so that invalid operations, such as mixing geocentric and
barycentric states or subtracting incompatible topocentric positions, are
rejected or explicitly checked.

The covered frames and centers are the ones needed for practical astronomy and
astrodynamics: inertial and Earth-fixed frames, horizontal coordinates,
geocentric / heliocentric / barycentric centers, and parameterized topocentric
or bodycentric centers. Observation-dependent effects are kept separate from
pure translations and rotations.

## Technical scope

What you will find here:

- `cartesian.rs` and `spherical.rs`: the core typed coordinate containers.
- `frames.rs`: frame marker types such as `ICRS`, `ECEF`, and equatorial /
  ecliptic variants.
- `centers.rs`: center marker types and parameterized site/body centers.
- `transform/`: frame rotations, center translations, and transformation
  context types.
- `observation/`: observer-state and observational-direction helpers.
- `horizontal.rs`: horizontal-coordinate convention helpers.
- `types.rs`: concise aliases for frequently used coordinate forms.

The module is built on top of `affn`, but it is the astronomy-specific layer
that binds geometric vectors to frames, centers, and observing semantics.

## References

- IERS Conventions (2010), IERS Technical Note 36.
  <https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36>
- International Astronomical Union. Resolutions adopted at the General
  Assemblies on celestial and terrestrial reference systems.
  <https://www.iau.org/Iau/Iau/Science/Resolutions.aspx>
- IAU SOFA Board. Standards of Fundamental Astronomy software collection.
  <https://www.iausofa.org/about-us>
