# Full SGP4/SDP4 in Siderust 0.15.0

`astro::sgp4::{TlePropagator, Sgp4Propagator}` uses one native Rust numerical
implementation of [Vallado et al., AIAA 2006-6753](https://celestrak.org/publications/AIAA/2006-6753/).
The former secular J2/Kepler approximation is removed. No new dependency is
introduced: typed angles, time, lengths and TEME states use the existing
`qtty`, `tempoch` and `affn` stack, including their `no_std` mathematics support.

The model includes near-Earth SGP4, deep-space SDP4, atmospheric drag,
short/long-period corrections, lunar/solar effects, Lyddane corrections,
synchronous and half-day resonances. Internal coefficients use the published
normalized Earth-radius/minute system; dimensional inputs and outputs remain
typed. The propagator is immutable, cloneable and thread safe. Every request
starts from epoch constants, so negative and out-of-order offsets work without
retaining mutable resonance integrator state. Successful propagation allocates
no heap memory.

| Model | Geopotential | Epoch/sidereal convention | Lyddane propagation |
|---|---|---|---|
| `Wgs72` (default, FFI 0) | WGS-72 | AFSPC | AFSPC compatibility |
| `Wgs72Iau` (FFI 1) | WGS-72 | IAU | improved |
| `Wgs84` (FFI 2) | WGS-84 | IAU | improved |

Constants include the complete J2/J3/J4 geopotential. WGS-72 uses Earth radius
6378.135 km and μ = 398600.8 km³/s²; WGS-84 uses 6378.137 km and
μ = 398600.5 km³/s². Both the Kozai-to-Brouwer conversion and perturbations
use the selected geopotential. Model differences are physical/conventional
and must not be interpreted as rounding errors. The default reproduces
Vallado/AFSPC reference vectors at the tested tolerance; bitwise identity across
platforms and math implementations is not promised.

Input TLE mean elements must not be replaced with osculating elements.
SGP4 uses mean motion and BSTAR; TLE mean-motion derivatives are metadata
and are not integrated independently. Output is **geocentric TEME** position
in km and velocity in km/s, using the existing `affn`/`qtty` types. Existing
coordinate transforms provide other frames, with appropriate Earth orientation
data. This upgrade intentionally changes results from the simplified
propagator, including the epoch state.

## Time and errors

Use `propagate_minutes(qtty::Minutes::new(offset))` (the Siderust `qtty`
re-export) or `propagate_at(JulianDate<UTC>)`. Callers with a custom Tempoch
snapshot can use `TlePropagator::from_tle_with_context(&tle, context)` or
`TlePropagator::from_tle_with_model_and_context(&tle, model, context)`;
the propagator owns that snapshot for all epoch conversions. Both use 1440 nominal minutes per
civil UTC day; a leap insertion does not add an SGP4 model second. Tempoch
performs UTC/calendar validation and validates typed duration conversions.
Leap-second labels are rejected because they have no scalar civil JD
representation. Calendar differences use nanosecond resolution, avoiding
large-JD subtraction and SI-duration errors across leap seconds.

Calendar conversion limits offsets to signed 64-bit nanoseconds (approximately
±292 years). This is a representation limit, not an accuracy promise. Useful
TLE prediction spans are far shorter and depend strongly on the object and
element age. A caller-supplied scalar Julian date has roughly 40 µs per ULP
at modern epochs; prefer typed calendar instants or minute offsets when this
matters. Nanosecond calendar truncation and f64 arithmetic remain finite
precision boundaries. AFSPC's legacy scalar-JD epoch rounding is retained for
initializing its lunar/solar and sidereal terms, separately from precise offset
arithmetic.

`InvalidElements` covers invalid eccentricity, non-positive mean motion,
non-finite propagation elements and initialization errors. `InvalidEpoch` and
`TimeConversion` cover calendar/UTC failures. `Propagation` reports divergent
eccentricity, negative semilatus rectum, invalid motion or non-finite output.
The reference contract allows sub-Earth states for decaying trajectories;
this is not a reentry or impact detector. FFI errors remain `InvalidArgument`
and output buffers are untouched on failure. Unrecognized model codes retain
the established WGS-72 fallback.

## Validation and downstream migration

The typed API checks the complete 33-satellite Vallado corpus: 666 Cartesian
states and two expected divergence errors. It also checks 363 vectors from the
independent Vallado C++ port, covering all models and synthetic boundary cases. Every successful reference
vector uses absolute component tolerances below **1e-6 km** and **1e-9 km/s**.
Epoch/minute consistency, fractional offsets, leap boundaries, invalid values
and FFI parity have separate regression tests. Fixture provenance, licensing
and reproduction are in
[`tests/fixtures/sgp4/README.md`](../tests/fixtures/sgp4/README.md).
The Rust `sgp4` and python-sgp4 projects are validation references only; their
MIT notices accompany the test fixtures. Neither is a Siderust dependency.

[spod's migration](https://github.com/VPRamon/spod/issues/29) can replace local
`spod::sgp4` imports with `siderust::astro::sgp4`, update Siderust to 0.15.0,
wrap scalar minute arguments in `Minutes::new`, and use context-aware Tempoch
methods for civil UTC. It can remove its local wrapper and direct `sgp4`
dependency. The existing downstream vectors are included upstream.
`Wgs72Iau` intentionally fixes the downstream initializer's accidental WGS-84
selection; callers needing those previous results must select `Wgs84`
explicitly. No spod source changes or release publication are included.

## External backend status

The validated native kernel remains in place pending a corrected upstream
`sgp4` release. In `sgp4` 2.4.0, the half-day resonance `G520` polynomial
uses `e < 0.715`, whereas Vallado C++ uses `e <= 0.715` for that branch.
Exactly 15 of the 315 independent synthetic reference states fail at the
boundary, by up to 0.110733 km and 2.38195e-5 km/s. Correcting this operator
in an isolated experiment restores parity. No experimental backend, fork or
diagnostic executable is included as a production dependency. See
[PR #112](https://github.com/Siderust/siderust/pull/112) for the investigation.
