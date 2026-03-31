# FFI Bindings Matrix

This matrix is the checked-in parity checklist for the Rust-first binding stack.
Update it whenever a public Rust concept changes exposure in `qtty-ffi`,
`tempoch-ffi`, `siderust-ffi`, `siderust-cpp`, `siderust-py`, or `siderust-js`.

## Layering Rules

- `qtty-ffi` is the only ABI source for scalar quantity and derived-quantity carriers.
- `tempoch-ffi` is the only ABI source for JD/MJD/period carriers and time-scale IDs.
- `siderust-ffi` owns astronomy-domain APIs, coordinate semantics, orbit handles, subject
  dispatch, target handles, context handles, and runtime ephemerides.
- `affn` remains Rust-only. Its semantics are preserved in `siderust-ffi` POD structs and
  adapter type systems rather than a standalone `affn-ffi` crate.
- C++ consumes the C ABI. Python and JavaScript bind Rust directly but must preserve the
  same semantics and naming contracts as the canonical FFI layer.

## Concept Matrix

| Rust concept | Canonical crate | `qtty-ffi` | `tempoch-ffi` | `siderust-ffi` | C++ | Python | JS | Notes |
|---|---|---|---|---|---|---|---|---|
| Scalar quantity carrier | `qtty` | `qtty_quantity_t` | reused | reused | `qtty-cpp` wrappers | `qtty-py` `Quantity` | `qtty-js` typed quantity | Do not replace with raw floats once unit semantics matter. |
| Derived quantity carrier | `qtty` | `qtty_derived_quantity_t` | reused | reused | `qtty-cpp` wrappers | `qtty-py` derived quantity | `qtty-js` derived quantity | Used for velocities/frequencies/rates crossing ABI. |
| Time instant carrier | `tempoch` | reused | `tempoch_jd_t`, `tempoch_mjd_t` | re-exported | `tempoch-cpp` wrappers | `tempoch-py` time objects | `tempoch-js` time objects | `siderust-ffi` should not duplicate JD/MJD layouts. |
| Time interval / duration carrier | `tempoch` | `qtty` durations when scalarized | `tempoch_period_mjd_t` and scale conversions | reused | `tempoch-cpp` wrappers | `tempoch-py` periods | `tempoch-js` periods | Keep malformed interval validation in `tempoch-ffi`. |
| Reference frame IDs | `siderust` / `affn` | no | no | `SiderustFrame` | RAII enums | mirrored enums/classes | mirrored enums/classes | Stable discriminants across adapters. |
| Reference center IDs | `siderust` / `affn` | no | no | `SiderustCenter` | RAII enums | mirrored enums/classes | mirrored enums/classes | Preserve position vs displacement semantics. |
| Spherical direction | `siderust::coordinates::spherical::Direction` | qtty angles | time input only | `SiderustSphericalDir` | wrapper value type | typed object | typed object | Value type, no heap ownership needed. |
| Spherical position | `siderust::coordinates::spherical::Position` | qtty length | time input only | `SiderustSphericalPos` | wrapper value type | typed object | typed object | Distance unit must stay explicit. |
| Cartesian position | `siderust::coordinates::cartesian::Position` | qtty length | time input only | `SiderustCartesianPos` | wrapper value type | typed object | typed object | Keep frame and center tags attached. |
| Generic subject dispatch | `AltitudePeriodsProvider`, `AzimuthProvider`, `Trackable` | reused where needed | reused where needed | `SiderustSubject` | wrapper sum type | typed union/class | typed union/class | Avoid per-body API duplication. |
| Generic target with epoch and PM | `CoordinateWithPM<T>` | PM rates via qtty | epoch via `TempochJd` | `SiderustGenericTarget` + `SiderustGenericTargetData` | RAII handle | typed target class | typed target class | Preserve coord kind, epoch, PM, and RA convention. |
| Transform/model context | `AstroContext`, `EarthOrientationModel` | no | TT/UT1 carriers | `SiderustContext` / `SiderustEarthOrientationModel` | RAII handle | typed context object | typed context object | Use handle-based API instead of duplicating context-specific symbol families. |
| Runtime ephemeris | `RuntimeEphemeris` | no | time carriers reused | runtime ephemeris handle APIs | RAII handle | typed object | typed object | Ownership must be explicit with paired destroy functions. |
| Observatory/geodetic site | `Geodetic<ECEF>` | height may use quantity | time input only | `SiderustGeodetict` | wrapper value type | typed object | typed object | Preserve geodetic semantics, not `(lon, lat, h)` tuples alone. |
| Altitude/azimuth/twilight queries | `siderust::calculus` | thresholds via quantities when applicable | query times via `tempoch` | C ABI query/result structs + subject/context handles | thin wrappers | typed methods | typed methods | One parity fixture set should cover all adapters. |
| Body handles and catalogs | `Sun`, `Moon`, planets, `Star`, asteroid/comet/satellite types | no | reused | value enums or opaque handles as needed | thin wrappers | typed classes | typed classes | Prefer canonical Rust naming. |
| Orbit families | `KeplerianOrbit`, `MeanMotionOrbit`, `ConicOrbit`, `PreparedOrbit` | scalar orbital params can reuse quantities | epochs via `tempoch` | POD orbit structs + prepared-orbit handle | RAII wrappers | typed classes | typed classes | Preserve reference-center/frame metadata. |

## Current Consolidation Targets

- Keep `siderust-ffi` aligned to `tempoch-ffi 0.4.x` and `qtty-ffi 0.4.x`.
- Route all new dimensional/time payloads through `qtty-ffi` and `tempoch-ffi` carriers.
- Prefer a compact handle-based `siderust-ffi` surface for contexts, targets, prepared
  orbits, and runtime ephemerides instead of symbol-per-special-case growth.
- Keep adapter public APIs semantically aligned with the Rust concepts named above.
