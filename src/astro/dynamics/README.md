# `astro/dynamics`

## Scientific scope

This module implements Cartesian spacecraft dynamics: orbit propagation under
multiple perturbation force models, numerical state-transition matrices, and
ancillary science models (atmospheric density, gravity-field providers).

## Module overview

```
astro/dynamics
├── state.rs              — OrbitState<C, F>, StateDerivative, spacecraft state
├── errors.rs             — DynamicsError, LocalFrameError
├── context.rs            — DynamicsContext (provider registry)
├── units.rs              — Typed unit aliases (GravitationalParameter, etc.)
├── frames.rs             — Local orbital frames (RTN, VNC, LVLH)
├── covariance.rs         — 6×6 state covariance with frame transport
├── stm.rs                — Finite-difference state-transition matrices
├── atmosphere.rs         — Atmospheric density providers
├── forces/               — Force models (two-body, J2, drag, SRP, third-body, etc.)
│   ├── mod.rs            — ForceModel trait & re-exports
│   ├── traits.rs         — ForceModel, ForcePartials abstractions
│   ├── two_body.rs       — Newtonian central gravity
│   ├── j2.rs             — First zonal perturbation
│   ├── composite.rs      — Linear sum of force models
│   ├── drag.rs           — Atmospheric drag (cannonball)
│   ├── srp.rs            — Solar radiation pressure
│   ├── third_body.rs     — Sun & Moon perturbations (Battin's form)
│   ├── geopotential.rs   — Spherical-harmonic gravity
│   ├── relativity.rs     — Schwarzschild correction (1PN)
│   └── empirical.rs      — User-supplied RTN acceleration
├── gravity/              — Geopotential coefficient providers
│   ├── mod.rs            — Provider trait & kernel
│   ├── provider.rs       — GravityFieldProvider trait
│   ├── acceleration.rs   — Spherical-harmonic kernel
│   └── egm_low_degree.rs — TwoBodyEarth, LowDegreeEarth
├── integrators/          — Numerical ODE integrators
│   ├── mod.rs            — Trait abstractions
│   ├── rk4.rs            — Classical RK4 (fixed-step)
│   ├── dopri5.rs         — Adaptive DOPRI5
│   └── dop853.rs         — Adaptive DOP853
├── propagation/          — High-level propagation driver
│   ├── mod.rs            — Driver overview
│   ├── driver.rs         — propagate() entry point
│   ├── config.rs         — PropagationConfig
│   ├── result.rs         — PropagationResult
│   ├── events.rs         — Event detectors
│   └── error.rs          — PropagationError
└── variational/          — Analytic STM via variational equations
    ├── mod.rs            — Module overview
    ├── equations.rs      — A-matrix construction
    └── propagator.rs     — RK4 STM integration
```

## Typical workflow

```rust
use std::sync::Arc;
use siderust::astro::dynamics::{
    OrbitState, Position, Velocity,
    forces::{TwoBody, J2, CompositeForce, ForceModel},
    context::DynamicsContext,
    gravity::LowDegreeEarth,
    integrators::{AdaptiveStepper, Dopri5},
    propagation::{PropagationConfig, propagate},
    stm::finite_diff_stm,
};
use siderust::coordinates::frames::GCRS;
use siderust::time::JulianDate;
use siderust::qtty::{Second, IntegratorTolerances};

// 1. Build a typed orbit state in GCRS
let epoch = JulianDate::new(2_451_545.0).to_time();
let t_end = JulianDate::new(2_451_546.0).to_time();  // 1 day later
let pos = Position::new(7000.0, 0.0, 0.0);
let vel = Velocity::new(0.0, 7.5, 0.0);
let state = OrbitState::new(epoch, pos, vel);

// 2. Set up dynamics context with providers
let ctx = DynamicsContext::builder()
    .with_gravity(Arc::new(LowDegreeEarth))
    .build();

// 3. Compose a force model: two-body + J2
let two_body = TwoBody::earth();
let j2 = J2::earth();
let force = CompositeForce::empty()
    .push(Box::new(two_body))
    .push(Box::new(j2));

// 4. Propagate with adaptive integrator
let config = PropagationConfig::new(epoch, t_end)
    .with_initial_step(Second::new(30.0))
    .with_max_step(Second::new(600.0));

let integrator = Dopri5::new(IntegratorTolerances::uniform(1e-9, 1e-3, 1e-6));
let result = propagate(&integrator, &force, state, &config, &ctx)?;

// 5. Optionally compute a state-transition matrix (RK4 finite-difference)
let stm = finite_diff_stm(&force, state, Second::new(3600.0), 120, &ctx)?;

println!("Propagated {} steps", result.steps_taken);
println!("Final position: {:?}", result.samples.last().unwrap().position);
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Technical scope

All public APIs use [`qtty`] typed quantities — no raw `f64` values with
unit-encoded names escape the module boundary.  Coordinates are [`affn`]
vectors, giving compile-time frame and unit guarantees.

## Unit conventions

| Quantity | Type | Notes |
|---|---|---|
| Position | [`Kilometers`] | Geocentric, GCRS frame |
| Velocity | `km/s` (`Per<Kilometer, Second>`) | GCRS frame |
| Acceleration | `km/s²` | GCRS frame |
| Time step | [`Second`] | Passed to all integrators |
| GM | [`GravitationalParameter`] | km³/s² (EGM2008/WGS-84 convention) |
| Density | [`KilogramsPerCubicMeter`] | SI |
| Area/mass | `m²/kg` (`Per<SquareMeter, Kilogram>`) | SI |
| Pressure | [`Pascals`] | SI; used for SRP reference flux P₀ |

## Submodule reference

### Core state & context

- [`state`] — `OrbitState<C, F>`, `StateDerivative<F>` (frame-tagged typed state)
- [`errors`] — `DynamicsError`, `LocalFrameError` (comprehensive error types)
- [`context`] — `DynamicsContext` (provider registry thread through all calls)
- [`units`] — Typed unit aliases re-exported from `qtty`

### Frames & covariance

- [`frames`] — Local orbital frames (RTN, VNC, LVLH) with typed rotation matrices
- [`covariance`] — 6×6 block-diagonal state covariance with frame transport
- [`stm`] — Finite-difference state-transition matrices (legacy, alongside variational)

### Force models

- [`forces`] — Two-body, J2, drag, SRP, third-body, geopotential, relativity, empirical
- [`forces::composite`] — Linear sum of force models
- [`gravity`] — Geopotential field provider trait & low-degree implementations

### Integration & propagation

- [`integrators`] — RK4 (fixed), DOPRI5, DOP853 (adaptive)
- [`propagation`] — High-level driver with event detection
- [`variational`] — Analytic STM via variational equations (RK4)

### Atmosphere

- [`atmosphere`] — Density provider trait & built-in models (exponential, NRLMSISE-00 lite)

## Coverage

Every `.rs` file under `src/astro/dynamics/` is kept at ≥ 90 % line coverage.
Run the enforcement script from the `siderust` crate directory to verify:

```sh
bash scripts/check_dynamics_coverage.sh
```

The script runs `cargo llvm-cov`, parses the JSON report, and exits non-zero
(printing a table of offending files) if any dynamics file falls below the
threshold.
