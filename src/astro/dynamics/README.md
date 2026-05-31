# `astro/dynamics`

## Scientific scope

This module hosts astronomy-specific spacecraft dynamics on top of the generic
`principia` mechanics kernel. `siderust` keeps the runtime astronomy semantics
here: atmosphere and ephemeris context wiring, Earth gravity helpers,
spacecraft properties, and perturbation models such as drag, SRP, third-body,
geopotential, relativity, and empirical RTN accelerations.

## Module overview

```text
astro/dynamics
├── state.rs              — TT-fixed `OrbitState<C, F>` alias plus spacecraft helpers
├── errors.rs             — `DynamicsError`, `LocalFrameError`
├── context.rs            — `DynamicsContext` provider registry
├── units.rs              — typed aliases such as `GravitationalParameter`
├── atmosphere.rs         — density-provider models
├── forces/               — astronomy-specific perturbations + principia model re-exports
│   ├── mod.rs            — `AccelerationModel`, `CompositeModel`, `TwoBody`, `J2`
│   ├── drag.rs           — atmospheric drag
│   ├── srp.rs            — solar radiation pressure and eclipse models
│   ├── third_body.rs     — Sun and Moon perturbations
│   ├── geopotential.rs   — Earth spherical harmonics wrapper
│   ├── relativity.rs     — central-body 1PN relativity
│   └── empirical.rs      — constant RTN empirical acceleration
├── gravity/              — Earth gravity providers and low-degree models
│   ├── mod.rs            — provider/kernel re-exports from principia
│   └── egm_low_degree.rs — `TwoBodyEarth`, `LowDegreeEarth`
├── integrators/          — principia integrator re-exports
├── propagation/
│   ├── mod.rs            — propagation re-exports
│   └── propagator.rs     — astronomy-facing `Propagator` facade
└── variational/          — principia STM/variational re-exports
```

## Typical workflow

```rust
use siderust::astro::dynamics::{
    dopri5_propagate, CompositeModel, DynamicsContext, IntegratorTolerances,
    OrbitState, Position, Velocity,
};
use siderust::astro::dynamics::forces::{EARTH_J2, J2, TwoBody, GM_EARTH, R_EARTH};
use siderust::coordinates::frames::GCRS;
use siderust::qtty::Second;
use siderust::time::JulianDate;

let state = OrbitState::new(
    JulianDate::new(2_451_545.0).to_j2000s(),
    Position::<GCRS>::new(7000.0, 0.0, 0.0),
    Velocity::<GCRS>::new(0.0, 7.5, 0.0),
);
let model = CompositeModel::empty()
    .push(Box::new(TwoBody::new(GM_EARTH)))
    .push(Box::new(J2::new(GM_EARTH, R_EARTH, EARTH_J2)));
let ctx = DynamicsContext::empty();

let propagated = dopri5_propagate(
    &model,
    state,
    Second::new(60.0),
    IntegratorTolerances::uniform(1e-9, 1e-6, 1e-9),
    &ctx,
)?;

assert!(propagated.epoch.value() > state.epoch.value());
# Ok::<(), siderust::astro::dynamics::PrincipiaError>(())
```

## Technical scope

- Public propagated states are `principia::DynamicsState<TT, C, F>` through the
  `OrbitState<C, F>` alias.
- Generic propagators, integrators, covariance helpers, local-trajectory
  frames, and variational kernels are re-exported from `principia`.
- Astronomy-specific runtime dependencies stay in `DynamicsContext`; generic
  `principia` models do not own atmosphere, ephemeris, or EOP state.
- Earth low-degree gravity providers remain in `astro/dynamics/gravity`.
