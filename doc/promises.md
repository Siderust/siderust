# Project Promises (User-Facing Guarantees)

This file lists *behavioral guarantees* that we treat as part of the public
contract. If something here needs to change, it should be reflected as an API
change (and documented in the changelog).

## Safety & correctness

- `siderust` contains **no `unsafe` blocks** in `src/`.
- Coordinate semantics are enforced through types: you cannot combine different
  centers/frames/units without an explicit transform.
- Parameterized-center operations that require parameter equality either:
  - **panic** (default operators / convenience APIs), or
  - return a **typed error** via the `checked_*` / `try_*` APIs.

## Time axis conventions

- Time-dependent transforms require a **TT epoch** (`JulianDate` as used by the
  astronomy APIs). Where UT1 is required, APIs either compute it from a context
  or require it explicitly.
- The altitude API interprets `Period<MJD>` / `ModifiedJulianDate` inputs on the
  **TT axis** (see `CHANGELOG.md` for the rationale).

## Feature flags & build behavior

- JPL DE4xx ephemerides are loaded at **runtime** via [`RuntimeEphemeris`] and
  optional `runtime-data` (not compile-time `de440`/`de441` features).
- `SIDERUST_JPL_STUB` is retained for CI/local determinism with
  `cargo test --all-features`; it does not embed JPL data in the crate.

[`RuntimeEphemeris`]: https://docs.rs/siderust/latest/siderust/ephemeris/struct.RuntimeEphemeris.html

## Compatibility

- Public APIs follow SemVer (`Cargo.toml` versioning).
- The crate targets Rust 2021 edition. (MSRV is not currently pinned as a strict promise.)
