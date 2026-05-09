# `targets`

## Scientific scope

An astronomical target is a coordinate statement plus the epoch to which that
statement applies. For stars and many catalogue objects, the model also carries
proper motion so that the nominal direction can be propagated to nearby epochs.
This module packages that standard astrometric concept into the typed data
structures used by `siderust`.

The scope is intentionally limited to target identity and epoch propagation. It
does not itself implement parallax, aberration, or ephemeris generation; those
belong to `astro` and `calculus`.

## Technical scope

What you will find here:

- `target.rs`: `CoordinateWithPM<T>` and the backward-compatible `Target<T>`
  alias.
- `trackable.rs`: the `Trackable` trait for anything that can provide a
  position or direction at a requested epoch.
- `transform.rs`: support glue used by the target abstractions.

The key output of the module is a typed target record that couples a coordinate,
an epoch, and an optional proper-motion model without erasing frame, center, or
unit semantics.

## References

- International Astronomical Union (2006). Resolution B1.5 and associated
  proper-motion formalism for modern reference systems.
- Kovalevsky, J., and Seidelmann, P. K. (2004). Fundamentals of Astrometry.
  Cambridge University Press.
- ESA (1997). The Hipparcos and Tycho Catalogues, ESA SP-1200.
