# SGP4/SDP4 TLE Propagation

The `astro::sgp4` module implements the near-Earth SGP4 and deep-space SDP4
models described by [Vallado et al., AIAA 2006-6753](https://celestrak.org/publications/AIAA/2006-6753/).
The numerical implementation is written in Rust and uses the existing
`qtty`, `tempoch` and `affn` types. It supports `no_std` environments.

`TlePropagator` (also exported as `Sgp4Propagator`) accepts TLE **mean
elements**, not osculating Keplerian elements. The model accounts for secular
and periodic terms, BSTAR atmospheric drag, lunar and solar perturbations, and
deep-space synchronous and half-day resonances. TLE mean-motion derivatives
are metadata; they are not integrated separately.

## Gravity models and reference frame

| `GravityModel` | Geopotential | Sidereal/epoch convention | FFI code |
| --- | --- | --- | --- |
| `Wgs72` (default) | WGS-72 | AFSPC | 0 |
| `Wgs72Iau` | WGS-72 | IAU improved | 1 |
| `Wgs84` | WGS-84 | IAU improved | 2 |

The geopotential includes J2, J3 and J4. WGS-72 uses an Earth radius of
6378.135 km and gravitational parameter μ = 398600.8 km³/s²; WGS-84 uses
6378.137 km and μ = 398600.5 km³/s². Model selection affects both
initialization and propagation. The default is the Vallado/AFSPC-compatible
choice; outputs from different conventions are not interchangeable.

Propagation returns typed geocentric **TEME** Cartesian positions in
kilometres and velocities in kilometres per second. The propagator is
immutable and cloneable; evaluations do not depend on request order and
successful propagation does not allocate heap memory.

## Public API

Construct a propagator with `TlePropagator::from_tle` or select an explicit
gravity model using `TlePropagator::from_tle_with_model`. Both also have
context-aware counterparts,
`from_tle_with_context` and `from_tle_with_model_and_context`,
which retain an owned `tempoch::TimeContext` snapshot.

```rust
use siderust::astro::sgp4::{GravityModel, TlePropagator};
use siderust::formats::tle::parse_3le;
use siderust::qtty::Minutes;

let tle = parse_3le(
    "ISS (ZARYA)",
    "1 25544U 98067A   08264.51782528 -.00002182  00000-0 -11606-4 0  2927",
    "2 25544  51.6416 247.4627 0006703 130.5360 325.0288 15.72125391563537",
).unwrap();

let propagator = TlePropagator::from_tle_with_model(&tle, GravityModel::Wgs72).unwrap();
let state = propagator.propagate_minutes(Minutes::new(60.0)).unwrap();
assert!(state.position().distance().value() > 0.0);
```

Use `propagate_minutes(Minutes)` for an offset from the TLE epoch or
`propagate_at(JulianDate<UTC>)` for a UTC target epoch. Both return
`Result<TemeState, Sgp4Error>`.

## Time convention and precision

SGP4 uses **1,440 nominal minutes per civil UTC day**, not elapsed SI seconds.
Leap-second labels cannot be represented on this scalar UTC grid and are
rejected. Tempoch performs UTC/calendar validation using the propagator's
time-data context; calendar-based offset arithmetic retains nanosecond
resolution and avoids subtracting large floating-point Julian dates.

Calendar offset conversion is bounded by signed 64-bit nanoseconds
(approximately ±292 years). This is a representation limit, not an orbital
prediction accuracy guarantee. Scalar Julian dates have finite precision
(approximately 40 microseconds per ULP at modern epochs); use minute offsets
or precise calendar instants when sub-millisecond timing matters. Valid TLE
prediction horizons are typically much shorter and depend on the orbit and
the age of the elements.

## Errors and physical limitations

`Sgp4Error::InvalidElements` covers invalid or unsupported TLE mean
elements; `InvalidEpoch` and `TimeConversion` cover UTC/calendar failures;
`Propagation` covers numerical divergence or non-physical propagation
intermediates. As with the reference implementation, a propagated state may
lie below the Earth's surface. The propagator does not detect impact or
reentry.

A frame rotation may transform a position or an independent vector, but a
**coordinate velocity** in a time-dependent frame also requires the
`Rdot r` term and, when the center moves, the origin velocity. The available
frame-rotation providers do not themselves implement a complete six-dimensional
TEME state transformation.

## Numerical validation

Regression tests compare the typed API and FFI outputs against published
Vallado verification states and independently generated C++ reference
vectors. The cases exercise near-Earth and deep-space regimes, gravity model
selection, resonances, negative and fractional offsets, and error boundaries.
Successful reference states are checked against strict per-component absolute
tolerances of 1e-6 km (position) and 1e-9 km/s (velocity); bitwise agreement
across platforms is not guaranteed.

Fixture sources, licensing information and reproduction steps are documented
in [the SGP4 reference fixtures](../../../tests/fixtures/sgp4/README.md).
