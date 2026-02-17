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

- `de440` / `de441` enable the corresponding JPL ephemeris backends.
- Enabling JPL features may trigger build-time downloads of NAIF BSP files.
- `SIDERUST_JPL_STUB` is an explicit opt-in for fast/offline builds:
  - `de441` is **mocked** to the analytical backend when stubbed (to keep tests runnable).
  - `de440` is **compile-only** when stubbed; runtime calls are expected to panic.

## Compatibility

- Public APIs follow SemVer (`Cargo.toml` versioning).
- The crate targets Rust 2021 edition. (MSRV is not currently pinned as a strict promise.)
