# `observatories`

## Scientific scope

Ground-based astronomy is site-dependent. A topocentric computation requires a
geodetic position on a reference ellipsoid and, for many optical applications,
a representative local atmosphere. This module stores a small set of
citation-backed observatory constants so that the rest of the crate can refer to
canonical sites instead of ad hoc longitude / latitude / height triples.

The coordinates are geodetic WGS84 coordinates, not spherical Earth
approximations. That distinction matters for topocentric transforms, Earth
rotation, and any conversion to Earth-centered Earth-fixed coordinates.

## Technical scope

What you will find here:

- `Observatory`: named site record with geodetic position, reference pressure,
  and optional reference temperature / relative humidity.
- Named constants such as `ROQUE_DE_LOS_MUCHACHOS`, `EL_PARANAL`,
  `MAUNA_KEA`, and `LA_SILLA_OBSERVATORY`.

These values are the canonical site definitions used by higher-level
topocentric, atmospheric, and sky-brightness code in `siderust`.

## References

- National Imagery and Mapping Agency (2000). Department of Defense World
  Geodetic System 1984: Its Definition and Relationships with Local Geodetic
  Systems, TR8350.2.
- Bowring, B. R. (1976). "Transformation from spatial to geographical
  coordinates." Survey Review 23(181), 323-327.
- European Southern Observatory. Paranal Observatory site information.
  <https://www.hq.eso.org/sci/facilities/paranal/astroclimate/site.html>
- Instituto de Astrofisica de Canarias. Observatorios de Canarias.
  <https://www.iac.es/en/observatorios-de-canarias>
