# `bodies`

## Scientific scope

This module is the object catalogue of `siderust`: the Sun, planets, moons,
small bodies, and named stars that appear in the library's dynamical and
observational APIs. It gathers physical parameters and identity-level
metadata that are stable enough to be treated as canonical constants, while
leaving time-dependent dynamics to the ephemeris and orbit modules.

The intent is not to reproduce a complete astronomical database. Instead, the
module provides a curated set of frequently used bodies and catalog entries
with typed quantities so that mass, radius, albedo, and orbital-element data
can be consumed without unit ambiguity.

## Technical scope

What you will find here:

- `solar_system.rs`: top-level Solar System constants and aggregates.
- `planets.rs`: the `Planet` type and planetary constant records.
- `satellite.rs`: natural satellites and related constants.
- `comet.rs`: comet records and named examples such as Halley and Encke.
- `asteroid.rs`: asteroid records and example near-Earth objects.
- `stars.rs`: the `Star` type and stellar parameter helpers.
- `catalog.rs`: a curated bright-star catalogue.

The module re-exports the most commonly used body types and constants so that
callers can work from `siderust::bodies::*` without navigating the lower-level
files. All public scalar properties are expressed as typed `qtty` values.

## References

- NASA GSFC National Space Science Data Center. Planetary Fact Sheets.
  <https://nssdc.gsfc.nasa.gov/planetary/planetfact.html>
- NASA JPL Solar System Dynamics. Small-Body Database and related services.
  <https://ssd.jpl.nasa.gov/>
- International Astronomical Union (2006). Resolution B5, "Definition of a
  Planet in the Solar System."
  <https://www.iau.org/static/resolutions/Resolution_GA26-5-6.pdf>
