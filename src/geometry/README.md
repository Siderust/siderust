# `geometry`

## Scientific scope

This module contains geometric helpers that are astronomy-adjacent but not tied
to any particular ephemeris, body model, or time standard. The current focus is
sky tessellation: partitioning the visible hemisphere into a finite set of
sample directions for integrations, samplers, and sky-coverage calculations.

These tools belong below the level of astrometric conventions and above the
level of pure vector algebra. They are intended for reusable geometric work
such as sky maps, radiance integration, or Monte Carlo sampling in horizontal
coordinates.

## Technical scope

What you will find here:

- `sky_grid.rs`: the `SkyGrid` and `SkyGridCell` types for uniform and
  approximately equal-area altitude-azimuth tilings of the upper hemisphere.

The grid yields typed horizontal directions and typed per-cell solid angles, so
downstream code can remain dimensionally explicit even in purely geometric
sampling tasks.

## References

- Gorski, K. M., Hivon, E., Banday, A. J., et al. (2005). "HEALPix: A
  Framework for High-Resolution Discretization and Fast Analysis of Data
  Distributed on the Sphere." Astrophysical Journal 622, 759-771.
- Arvo, J. (1995). "Stratified sampling of spherical triangles." In SIGGRAPH
  '95 Proceedings, 437-438.
