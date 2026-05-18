# `astro`

## Scientific scope

This module contains the astrometric reduction models that connect inertial or
catalogue quantities to physically observed directions and state vectors. It
covers the standard Earth-orientation chain used in modern positional
astronomy: precession, nutation, Celestial Intermediate Pole / Origin
quantities, Earth Rotation Angle, sidereal time, polar motion, stellar
aberration, gravitational light deflection, and proper-motion propagation.
Orbital and conic utilities live here as well because they provide the
time-dependent dynamical state used by the higher-level coordinate and event
pipelines.

The validity regime is the one implied by the IAU 2000/2006 resolutions, the
IERS terrestrial-to-celestial transformation conventions, and the SOFA-style
algorithm set from which these routines are derived. This is the canonical
place in `siderust` for standards-backed astrometry rather than numerical
ephemeris table evaluation.

## Technical scope

What you will find here:

- `aberration.rs`: stellar aberration using the full Lorentz transform.
- `cio.rs`: CIP `(X, Y)` and CIO locator `s`.
- `earth_rotation.rs`, `earth_rotation_provider.rs`, `era.rs`, `sidereal.rs`,
  `polar_motion.rs`: Earth-rotation and terrestrial/celestial orientation
  building blocks.
- `eop.rs`, `iers_data.rs`: typed Earth Orientation Parameter access.
- `precession.rs`, `nutation/`: IAU precession-nutation models and tables.
- `light_deflection.rs`: relativistic deflection by the Sun and planets.
- `proper_motion.rs`: catalogue-style proper-motion propagation.
- `conic.rs`, `orbit.rs`: Keplerian and conic-orbit abstractions.
- `orientation.rs`: IAU pole and prime-meridian rotation parameters for body
  frames.
- `units.rs`: astronomy-specific scalar units not provided by `qtty`.
- `dynamics/`: spacecraft dynamics — Cartesian inertial state, force models
  (two-body, J₂, drag, SRP, third-body), numerical integrators (RK4 / DOPRI5),
  state-transition matrices, covariance transport, gravity-field abstraction,
  and atmospheric density models.  All public APIs use `qtty` typed quantities.

All public entry points are expressed in typed `qtty` and `tempoch` values.
Where multiple standard models exist, selection is represented by marker types
at compile time rather than runtime enums.

## References

- IERS Conventions (2010), IERS Technical Note 36.
  <https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36>
- International Astronomical Union. Resolutions adopted at the 2000, 2006, and
  later General Assemblies.
  <https://www.iau.org/Iau/Iau/Science/Resolutions.aspx>
- IAU SOFA Board. Standards of Fundamental Astronomy software collection.
  <https://www.iausofa.org/about-us>
