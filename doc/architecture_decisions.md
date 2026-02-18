# Architectural Decisions

This document captures the deliberate design choices that define the crate’s public API and internal architecture. Its purpose is not merely descriptive, but normative: it establishes constraints and clarifies trade-offs so that future evolution remains coherent. Contributors should be able to understand not only what was done, but why it was done that way, and what consequences follow from those decisions.

## Coordinates as Type-Level Semantic Objects

Coordinates are modeled as semantic types whose parameters encode three distinct dimensions: center (`C`), reference frame (`F`), and unit (`U`). These parameters are part of the type itself for both Cartesian and spherical representations.

The motivation is correctness by construction. Physically invalid operations, such as mixing incompatible centers, frames, or units, are prevented at compile time rather than deferred to runtime checks or informal conventions. Instead of relying on documentation or discipline, the type system enforces invariants directly.

The result is an API that is more strongly typed and sometimes more verbose. However, this explicitness is intentional. The compiler becomes an active participant in maintaining physical consistency, significantly reducing the risk of subtle misuse in downstream code.

## Parameterized Centers and Runtime Equality

Certain centers, such as `Topocentric` (via `ObserverSite`) and `Bodycentric` (via `BodycentricParams`), cannot be fully expressed at the type level without causing an impractical explosion of distinct types. While the *kind* of center is encoded statically, the specific parameters associated with a given instance are validated at runtime.

This hybrid approach reflects a practical constraint. Encoding each possible observer site or body configuration as a unique type would be theoretically pure but operationally untenable. Instead, the design enforces structural correctness at compile time and parameter equality at runtime.

The consequence is that mismatches can trigger panics in default operations. For code that requires graceful failure semantics, fallible variants are provided. This keeps high-level usage ergonomic while preserving strict validation options for library or infrastructure code.

## Physically Meaningful Separation of Transforms

Transformations are separated according to their physical meaning rather than bundled under a single abstraction. Frame rotations are implemented as pure rotation matrices. Center shifts are modeled as affine translations. Observation effects, such as aberration or apparent place corrections, are treated as a distinct processing pipeline.

Each of these operations carries different invariants and requires different contextual inputs. For example, horizontal transformations depend on an observer site and Earth orientation parameters such as UT1 and polar motion. By keeping these domains distinct, the API makes such dependencies explicit rather than implicit.

This design avoids hidden assumptions, particularly regarding Earth orientation and observational context. The user is required to supply the appropriate information where necessary, reinforcing clarity over convenience.

## Explicit, Context-Driven Earth Orientation

High-precision topocentric and horizontal calculations depend on Earth orientation parameters. To accommodate this, the crate introduces `AstroContext`, which encapsulates EOP providers such as `IersEop` or `NullEop`.

The design principle is that precision should be an explicit choice. Some applications require high-fidelity UT1 and polar motion modeling; others may prioritize speed or operate in contexts where such corrections are unnecessary. By making the context explicit, callers can select the appropriate trade-off.

Convenience paths exist for common scenarios, but context-aware `_with_ctx` variants are available for full control. This dual-layer approach keeps typical workflows straightforward while preserving the ability to perform rigorous, high-precision computations.

## Pluggable Ephemeris Backends

Ephemeris computation is abstracted behind an `Ephemeris` trait. An analytical backend based on VSOP87 and ELP2000 is always available, ensuring a lightweight baseline. More precise JPL DE4xx backends, such as DE440 or DE441, are available behind feature flags.

This design reflects the diversity of use cases. Many applications do not require the mass and precision of DE4xx datasets, and including them by default would unnecessarily increase build size and complexity. Feature gating allows users to opt into heavier datasets only when needed.

At compile time, `DefaultEphemeris` selects the best available backend based on enabled features. This preserves a clean API surface while allowing build-time configurability.

## Build-Time Dataset Generation and Embedding

Datasets, including VSOP87, ELP2000, IERS EOP, and optionally DE4xx, are generated and embedded at build time. This eliminates runtime data dependencies and ensures deterministic behavior once the crate is compiled.

The advantage is twofold. First, runtime performance is improved because all required data is locally embedded. Second, operational environments remain simpler because no external data loading is required after build.

When JPL features are enabled, builds may download large datasets. To support fast or offline development workflows, the `SIDERUST_JPL_STUB` option provides a lightweight alternative for iteration.

## Safety-First Policy

The crate intentionally avoids `unsafe` blocks in its primary implementation. Given the numerical intensity of the domain, undefined behavior would be particularly difficult to detect and potentially catastrophic in its consequences.

Performance considerations are addressed through algorithmic improvements, careful use of safe Rust abstractions, and feature-gated dependencies rather than by relaxing safety guarantees. This policy reflects a prioritization of correctness, maintainability, and long-term reliability over marginal gains from low-level optimizations.
