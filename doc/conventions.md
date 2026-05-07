# Authoring Conventions

This document is the canonical contributor reference for three cross-cutting
rules that every public surface of the `siderust` crate must follow:

1. **Academic-style module documentation.**
2. **Typed `qtty` and `tempoch` quantities at the API boundary.**
3. **Phantom-type model selection instead of runtime enums.**

These rules apply to all new code and to any module touched as part of an
existing change.

---

## 1. Module documentation template

Every `*.rs` file with public items begins with a module-level `//!` block
in the following shape:

```rust
//! # <Module title>
//!
//! ## Scientific scope
//!
//! 1–3 paragraphs framing the *physical / astronomical* problem this
//! module solves: what quantity is being modelled, the regime of
//! validity (e.g. zenith distances < 85°, optical wavelengths,
//! geocentric vs barycentric, equinox J2000.0), and the documented
//! limitations or known bias regimes.
//!
//! ## Technical scope
//!
//! What this module provides at the *API* level: the entry points, what
//! they take, what they return, and what they explicitly do *not* do
//! (i.e. delegations to other modules). State the units of inputs and
//! outputs in terms of their typed `qtty::*` / `tempoch::*` form.
//!
//! ## References
//!
//! - Author, A. B. (Year). "Title". Journal, Volume(Issue), pages.
//!   doi:DOI-or-URL
//! - Standards body (Year). Standard / convention name. Edition.
```

Every public item additionally carries `# Arguments`, `# Returns`,
`# Errors` (for `Result`-returning fns), `# Panics` (for fns that may
panic by design), and at least one runnable `# Examples` block for any
non-trivial entry point.

`# Examples` blocks must compile under `cargo test --doc`, so prefer
small, deterministic snippets over long pipelines.

The reference style above is mandatory; use the bibliographic block
even when there is only one citation. Cite IAU SOFA, NIST, IERS
conventions, and standard textbooks (Meeus, Seidelmann, Vallado) by
edition number where applicable.

## 2. Typed quantities at the API boundary

Public function signatures must express physical quantities as typed
`qtty::*` / `tempoch::*` values, never as bare `f64`:

```rust
// ❌ Discouraged
pub fn rayleigh_phase(cos_theta: f64) -> f64 { ... }

// ✅ Preferred
pub fn rayleigh_phase(cos_theta: Dimensionless) -> Dimensionless { ... }
```

Dimensionless physical ratios receive newtypes in `qtty`
(`OpticalDepth`, `Airmass`, `Albedo`, `IlluminationFraction`,
`Refractivity`, etc.) and are re-exported through `siderust::qtty` for
compatibility. Astronomy-specific dimensionless concepts such as
`CipCoordinate` remain in `siderust`. Use the existing newtypes when they
exist; add generic dimensionless concepts in the `qtty` crate first, then
re-export them from `siderust` if compatibility requires it.

`.value()` and `.raw()` are tools of the *math kernel*, not of the API.
They are permitted inside `pub(crate)` numerical inner loops (`vsop87`,
`elp2000`, `nut00a`, `jpl::eval`, `interp`, etc.) where peeling to `f64`
is required for performance, but the surface that calls those kernels
must be typed.

When you must round-trip through `f64` because an external source
(e.g. a built-in dataset table) is `f64`-valued, do the conversion in
*one* clearly-named helper local to the module rather than scattering
`.value()` through the public API.

## 3. Phantom-type model selection

When a function picks between alternative algorithms, conventions, or
models (nutation theories, airmass formulas, ephemeris backends,
extinction parameterizations…), the selection is expressed at compile
time via a marker trait + zero-sized type, never at runtime via an enum.

The canonical pattern lives in `atmosphere::airmass`:

```rust
pub trait AirmassFormula {
    const NAME: &'static str;
    fn airmass_from_radians(zenith: Radians) -> f64;
}

pub struct PlaneParallel;
impl AirmassFormula for PlaneParallel { /* ... */ }

pub struct Young1994;
impl AirmassFormula for Young1994 { /* ... */ }

pub fn airmass<F: AirmassFormula>(zenith: Radians) -> f64 {
    F::airmass_from_radians(zenith.value())
}
```

A second worked example is `astro::nutation::{NutationModel,
Iau2000A, Iau2000B, Iau2006, Iau2006A}` consumed by
`coordinates::transform::ModelContext`.

A small runtime-dispatch shim
(`fn dispatch_<name>(model: SomeEnum, …)`) may be retained at one
clearly marked seam if downstream code genuinely needs runtime
selection (e.g. a CLI flag → model). The enum must not appear inside
the canonical typed APIs, only inside that one shim.

## Co-authorship

Commits made under tasks that explicitly request "no Copilot
co-author" must omit the `Co-authored-by: Copilot` trailer. Otherwise,
the standard repository commit-trailer policy applies.
