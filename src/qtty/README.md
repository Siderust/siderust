# `qtty`

## Scientific scope

Physical astronomy code is only reliable when units are carried explicitly.
This module is the `siderust` compatibility facade over the standalone `qtty`
crate and exists so that public APIs can keep using typed lengths, angles,
times, pressures, optical depths, airmasses, and similar quantities without
unit ambiguity.

The scientific role of the module is therefore foundational rather than
algorithmic: it is the unit system that other modules build upon. In
particular, dimensionless astronomy-specific quantities such as the Celestial
Intermediate Pole coordinates are given their own type-level identity instead
of being treated as anonymous scalars.

## Technical scope

What you will find here:

- Re-exports of the canonical `qtty` unit markers, quantity aliases, and
  dimension families used throughout `siderust`.
- `dimensionless.rs`: compatibility exports for `Airmass`, `OpticalDepth`,
  `Albedo`, `IlluminationFraction`, `Refractivity`, and `Transmittance`.
- `CipCoordinate`: a local unit marker for IAU CIP `X` and `Y` coordinates.

This directory does not define the full quantity algebra itself; it stabilizes
the `siderust::qtty` surface while the canonical implementation lives in the
external `qtty` crate.

## References

- Bureau International des Poids et Mesures. The International System of Units
  (SI Brochure), 9th edition.
  <https://www.bipm.org/en/publications/si-brochure>
- IERS Conventions (2010), IERS Technical Note 36.
  <https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36>
