# `astro/dynamics`

## Scientific scope

This module implements Cartesian spacecraft dynamics: orbit propagation under
multiple perturbation force models, numerical state-transition matrices, and
ancillary science models (atmospheric density, gravity-field providers).

Typical workflow:

1. Build an [`OrbitState`] with a typed `(r, v)` in GCRS.
2. Compose a [`CompositeForce`] from [`ForceModel`] implementations.
3. Propagate with [`rk4_propagate`] or [`dopri5_propagate`].
4. Optionally compute the state-transition matrix with [`finite_diff_stm`].

## Technical scope

All public APIs use [`qtty`] typed quantities — no raw `f64` values with
unit-encoded names escape the module boundary.  Coordinates are [`affn`]
vectors, giving compile-time frame and unit guarantees.

## Submodules

| Submodule | Description |
|---|---|
| `state` | [`OrbitState`], [`StateDerivative`], and [`SpacecraftState`] primitives. |
| `forces` | [`ForceModel`] trait, [`TwoBody`], [`J2`], [`DragForce`], [`ThirdBodySunMoon`], [`CannonballSrp`], [`CompositeForce`]. |
| `atmosphere` | [`DensityProvider`] trait, [`ExponentialAtmosphere`], [`ConstantDensity`]. |
| `gravity` | [`GravityFieldProvider`] trait, [`GravityConstants`], [`TwoBodyEarth`]. |
| `frames` | Local orbital frame (RSW / RTN) transforms. |
| `covariance` | Frame-tagged 6×6 state covariance with transport. |
| `stm` | Finite-difference state-transition matrix ([`finite_diff_stm`], [`finite_diff_stm_series`]). |
| `integrators` | [`rk4_step`] / [`rk4_propagate`] and adaptive [`dopri5_step`] / [`dopri5_propagate`]. |
| `units` | Local typed unit aliases ([`GravitationalParameter`]) not yet in `qtty`. |

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
