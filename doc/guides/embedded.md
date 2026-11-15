# Use compact ephemerides offline

In your dependency's Cargo entry, set `default-features = false` and
`features = ["embedded-compact"]` to select the compact profile. The library
builds with `no_std` and `alloc`. The compact evaluator and the geometric
path exercised below require no heap allocation or I/O. Other library paths
may allocate, so check the application and its dependency features before
sizing memory. See [feature selection](../reference/features-and-data.md#cargo-features).

## Evaluate a supported epoch

The compact VSOP87E and ELP2000-82B series cover
**2458849.5 ≤ JD TDB < 2470172.5**, from 2020-01-01 up to 2051-01-01 on
the TDB calendar. Public queries take TT and convert to TDB before applying
the range check. Use the checked `Ephemeris::try_*` and `PlanetState::try_*`
methods for external dates. The infallible methods panic outside coverage.
If both ephemeris features are enabled, refer to `CompactEphemeris`
explicitly wherever the result must depend on the compact series.

The Sun and planets have barycentric positions in AU. The Moon has
geocentric positions in km. Both theories are rotated into the canonical
`EclipticMeanJ2000` axes using their respective native frame ties.
Earth barycentric velocity is in AU/day. The evaluator differentiates
positions using points 675 seconds on either side of the epoch, with a
second-order one-sided stencil at the interval edges.

## Compute lunar altitude, phase, and bright limb

From the checkout, run [the compact example](../../examples/25_embedded_compact.rs)
on a host:

```sh
cargo run --locked --example 25_embedded_compact --no-default-features --features embedded-compact
```

It evaluates checked geocentric Moon, heliocentric Earth, and heliocentric
Mars positions. Negating Earth's vector gives the geocentric Sun. Subtracting
Earth's vector from Mars's gives the geocentric Mars position. The example
forms topocentric true-of-date vectors once and reuses the Moon and Sun
vectors for phase and bright-limb orientation. It includes WGS84 parallax,
precession, explicit IAU 2000B nutation, and Earth rotation. It omits
light-time, aberration, gravitational deflection, and refraction. The example
uses `NullEop`, so it approximates UT1 by UTC and sets pole corrections
to zero. Use an appropriate [EOP provider](../concepts/time-and-earth-orientation.md)
if that approximation is too coarse.

[DE440s and ERFA comparisons](../../scripts/validation/README.md) check
sampled position, pointing, and Earth velocity against independent references.
The Archive reports document the coefficient selection and its error metrics for
[VSOP87](https://github.com/Siderust/archive/blob/main/reports/compact-vsop87.md)
and [ELP2000](https://github.com/Siderust/archive/blob/main/reports/compact-lunar-ephemeris.md).

The [ESP32 harness](../../scripts/compact-esp32) exercises the compact workload
with a null allocator. Application peak stack use depends on the target,
callers, and workload.
