// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Cannonball solar radiation pressure — `CannonballSrp` force model.
//!
//! ## Scope
//!
//! Provides [`CannonballSrp`], a cannonball SRP model with selectable eclipse
//! shadow modelling via [`ShadowModel`].
//!
//! ## Equations
//!
//! ```text
//! a_srp = −ν · Cr · P₀ · (AU / |r_sun→sat|)² · (A/m) · r̂_sun→sat
//! ```
//!
//! where `ν ∈ [0, 1]` is the eclipse factor (1 = full sun, 0 = full shadow).
//! The acceleration pushes the spacecraft away from the Sun.
//!
//! ## Solar constant
//!
//! P₀ = 4.560 × 10⁻⁶ N/m² at 1 AU (IAU 2015 nominal solar irradiance
//! 1361 W/m² divided by speed of light 299792.458 km/s × 10⁻³).
//!
//! ## Units
//!
//! Position km, P₀ N/m², area-to-mass m²/kg → acceleration km/s² after
//! dividing by 1 000.
//!
//! ## Frame/center assumptions
//!
//! `<Geocentric, GCRS>`.
//!
//! ## Shadow models
//!
//! | [`ShadowModel`] | Description |
//! |-----------------|-------------|
//! | `None`          | No eclipse — ν = 1 always |
//! | `Cylindrical`   | Binary umbra/sunlight, Earth radius cylinder |
//! | `Conical`       | Continuous ν ∈ [0, 1], umbra + penumbra geometry (default) |
//!
//! ## References
//!
//! * Vallado, *Fundamentals of Astrodynamics and Applications* (4th ed.), §8, §3.5.
//! * Montenbruck & Gill, *Satellite Orbits*, §3.4, §3.4.2.

use affn::cartesian::Displacement;

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::GCRS;
use crate::qtty::{AreaToMass, Kilometer, SrpCoefficient};

use super::third_body::sun_geocentric;
use super::traits::{ForceModel, AU_IN_KM, P0, R_EARTH};

/// IAU 2015 nominal solar radius (km).
const R_SUN_KM: f64 = 695_700.0;

// =============================================================================
// ShadowModel
// =============================================================================

/// Earth shadow model used by [`CannonballSrp`].
///
/// Controls how the eclipse factor ν ∈ [0, 1] is computed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ShadowModel {
    /// No eclipse modelling — ν = 1.0 always (spacecraft is always fully sunlit).
    None,

    /// Cylindrical shadow: binary ν ∈ {0, 1}.
    ///
    /// Earth casts an infinite cylinder of radius R_earth in the anti-Sun
    /// direction.  Satellites inside the cylinder are fully in umbra (ν = 0);
    /// all others are fully lit (ν = 1).  Simple and fast, but ignores the
    /// penumbra transition zone.
    Cylindrical,

    /// Conical shadow: continuous ν ∈ [0, 1].
    ///
    /// Uses the apparent-angular-radius geometry (Montenbruck & Gill §3.4.2)
    /// to compute the fraction of the solar disk that is unobscured by Earth.
    /// Handles full umbra, penumbra, annular eclipse, and full sunlight.
    /// This is the default.
    #[default]
    Conical,
}

// =============================================================================
// CannonballSrp
// =============================================================================

/// Cannonball solar radiation pressure (SRP) force model.
///
/// ```text
/// a_srp = −ν · Cr · P₀ · (AU / |r_sun→sat|)² · (A/m) · r̂_sun→sat
/// ```
///
/// The acceleration pushes the spacecraft away from the Sun.  The eclipse
/// factor ν is computed by the selected [`ShadowModel`].
///
/// ## Builder
///
/// ```rust,ignore
/// use siderust::astro::dynamics::forces::{CannonballSrp, ShadowModel};
/// use siderust::qtty::{AreaToMass, SrpCoefficient};
///
/// // Default: conical shadow model.
/// let srp = CannonballSrp::new(SrpCoefficient::new(1.5), AreaToMass::new(0.01));
///
/// // Cylindrical shadow:
/// let srp = CannonballSrp::new(SrpCoefficient::new(1.5), AreaToMass::new(0.01))
///     .with_shadow(ShadowModel::Cylindrical);
///
/// // No eclipse:
/// let srp = CannonballSrp::new(SrpCoefficient::new(1.5), AreaToMass::new(0.01))
///     .with_shadow(ShadowModel::None);
/// ```
///
/// ## Sun position
///
/// The Sun geocentric position is fetched from the [`DynamicsContext`] ephemeris
/// at each call.
#[derive(Debug, Clone, Copy)]
pub struct CannonballSrp {
    /// Radiation-pressure coefficient (dimensionless).
    pub cr: SrpCoefficient,
    /// Area-to-mass ratio (m²/kg).
    pub area_to_mass: AreaToMass,
    /// Selected eclipse shadow model.
    pub shadow: ShadowModel,
}

impl CannonballSrp {
    /// Build a cannonball SRP model with the default [`ShadowModel::Conical`] eclipse model.
    pub fn new(cr: SrpCoefficient, area_to_mass: AreaToMass) -> Self {
        Self {
            cr,
            area_to_mass,
            shadow: ShadowModel::default(),
        }
    }

    /// Override the eclipse shadow model.
    pub fn with_shadow(mut self, shadow: ShadowModel) -> Self {
        self.shadow = shadow;
        self
    }
}

// =============================================================================
// Eclipse helpers (private)
// =============================================================================

/// Cylindrical shadow factor.
///
/// Returns `0.0` (full shadow) if the satellite is inside the Earth's shadow
/// cylinder, `1.0` (full sunlight) otherwise.
///
/// ## Algorithm
///
/// Earth casts an infinite cylinder of radius `r_earth_km` aligned with the
/// Sun–Earth direction extending in the anti-Sun half-space.  A satellite is
/// in the cylinder when both:
///
/// 1. Its projection onto the anti-Sun unit vector is **positive** (it is on
///    the anti-Sun side of the Earth centre).
/// 2. Its perpendicular distance from the cylinder axis is less than
///    `r_earth_km`.
///
/// ## Arguments
///
/// * `r_sat`      — satellite position in GCRS (km), geocentric.
/// * `r_sun`      — Sun position in GCRS (km), geocentric.
/// * `r_earth_km` — Earth equatorial radius (km).
///
/// ## Pathological inputs
///
/// If the Sun is at the origin (zero vector), shadow cannot be determined;
/// returns `1.0` (assume full sunlight).
fn cylindrical_shadow_factor(r_sat: [f64; 3], r_sun: [f64; 3], r_earth_km: f64) -> f64 {
    let r_sun_mag = (r_sun[0] * r_sun[0] + r_sun[1] * r_sun[1] + r_sun[2] * r_sun[2]).sqrt();
    if r_sun_mag == 0.0 {
        // Degenerate: Sun at Earth centre — assume full sunlight.
        return 1.0;
    }
    // Unit vector from Earth centre toward the Sun.
    let sun_hat = [
        r_sun[0] / r_sun_mag,
        r_sun[1] / r_sun_mag,
        r_sun[2] / r_sun_mag,
    ];

    // Projection of satellite onto the Sun direction.
    let proj = r_sat[0] * sun_hat[0] + r_sat[1] * sun_hat[1] + r_sat[2] * sun_hat[2];

    // Satellite must be on the anti-Sun side (negative projection).
    if proj >= 0.0 {
        return 1.0;
    }

    // Perpendicular component of r_sat with respect to the Sun–Earth axis.
    let perp_sq = (r_sat[0] - proj * sun_hat[0]).powi(2)
        + (r_sat[1] - proj * sun_hat[1]).powi(2)
        + (r_sat[2] - proj * sun_hat[2]).powi(2);

    if perp_sq < r_earth_km * r_earth_km {
        0.0 // inside the shadow cylinder
    } else {
        1.0
    }
}

/// Conical shadow factor (Montenbruck & Gill §3.4.2 / Vallado §3.5).
///
/// Returns a value in `[0, 1]` representing the fraction of the solar disk
/// that is unobscured by the Earth.
///
/// * `1.0` — full sunlight (no overlap between solar and Earth disks).
/// * `0.0` — full umbra (Earth disk entirely covers the solar disk).
/// * `(0, 1)` — penumbra or annular eclipse (partial overlap).
///
/// ## Algorithm
///
/// From the satellite, compute the apparent angular radii of the Sun
/// (`ρ_s = arcsin(R_sun / d_sun)`) and Earth (`ρ_e = arcsin(R_earth / d_earth)`),
/// and the angular separation `θ` between their centres.
///
/// Four cases:
///
/// 1. **Full sunlight**: `θ ≥ ρ_s + ρ_e` → return `1.0`.
/// 2. **Full umbra**: `ρ_e ≥ ρ_s` and `θ ≤ ρ_e − ρ_s` → return `0.0`.
/// 3. **Annular eclipse**: `ρ_s > ρ_e` and `θ ≤ ρ_s − ρ_e` → Earth disk fully
///    inside Sun disk; return `1 − (ρ_e / ρ_s)²` (unobscured solar fraction).
/// 4. **Partial penumbra**: otherwise → return `1 − A_overlap / (π ρ_s²)`
///    using the standard lens-area formula for two overlapping circles.
///
/// ## Arguments
///
/// * `r_sat`         — satellite position in GCRS (km), geocentric.
/// * `r_sun`         — Sun position in GCRS (km), geocentric.
/// * `r_earth_km`    — Earth equatorial radius (km).
/// * `r_sun_radius_km` — Solar radius (km).
///
/// ## Pathological inputs
///
/// If the satellite coincides with the Earth centre or Sun centre (zero
/// distance), shadow geometry is undefined; returns `1.0` (assume full
/// sunlight) to avoid a divide-by-zero.
fn conical_shadow_factor(
    r_sat: [f64; 3],
    r_sun: [f64; 3],
    r_earth_km: f64,
    r_sun_radius_km: f64,
) -> f64 {
    // Distance from satellite to Earth centre.
    let d_earth = (r_sat[0] * r_sat[0] + r_sat[1] * r_sat[1] + r_sat[2] * r_sat[2]).sqrt();
    if d_earth == 0.0 {
        return 1.0; // degenerate: satellite at Earth centre
    }

    // Vector and distance from satellite to Sun.
    let sat_to_sun = [
        r_sun[0] - r_sat[0],
        r_sun[1] - r_sat[1],
        r_sun[2] - r_sat[2],
    ];
    let d_sun = (sat_to_sun[0] * sat_to_sun[0]
        + sat_to_sun[1] * sat_to_sun[1]
        + sat_to_sun[2] * sat_to_sun[2])
        .sqrt();
    if d_sun == 0.0 {
        return 1.0; // degenerate: satellite at Sun centre
    }

    // Apparent angular radii (arcsin, clamped to [-1, 1] for robustness).
    let rho_e = (r_earth_km / d_earth).min(1.0).asin();
    let rho_s = (r_sun_radius_km / d_sun).min(1.0).asin();

    // Angular separation between Earth centre and Sun centre as seen from satellite.
    // Direction to Earth centre from satellite: -r_sat / d_earth
    // Direction to Sun from satellite: sat_to_sun / d_sun
    let cos_theta =
        (-r_sat[0] * sat_to_sun[0] - r_sat[1] * sat_to_sun[1] - r_sat[2] * sat_to_sun[2])
            / (d_earth * d_sun);
    let theta = cos_theta.clamp(-1.0, 1.0).acos();

    // Case 1: full sunlight.
    if theta >= rho_s + rho_e {
        return 1.0;
    }

    // Case 2: full umbra (Earth disk fully contains Sun disk).
    if rho_e >= rho_s && theta <= rho_e - rho_s {
        return 0.0;
    }

    // Case 3: annular eclipse (Sun disk fully contains Earth disk).
    if rho_s > rho_e && theta <= rho_s - rho_e {
        // Earth disk area fully inside Sun disk; fraction occluded = (rho_e/rho_s)^2.
        return 1.0 - (rho_e / rho_s).powi(2);
    }

    // Case 4: partial penumbra — two overlapping circles.
    //
    // Intersection area of circles with radii rho_s (Sun) and rho_e (Earth)
    // and centre separation theta (Montenbruck & Gill §3.4.2 eq. 3.88):
    //
    //   A = rho_s² · arccos(x / rho_s) + rho_e² · arccos(y / rho_e)
    //       − sqrt((rho_s + rho_e + theta)(rho_s + rho_e − theta)
    //              (theta + rho_s − rho_e)(theta − rho_s + rho_e)) / 2
    //
    // where x = (theta² + rho_s² − rho_e²) / (2·theta)
    //       y = (theta² + rho_e² − rho_s²) / (2·theta)
    let x = (theta * theta + rho_s * rho_s - rho_e * rho_e) / (2.0 * theta);
    let y = (theta * theta + rho_e * rho_e - rho_s * rho_s) / (2.0 * theta);

    let alpha = (x / rho_s).clamp(-1.0, 1.0).acos();
    let beta = (y / rho_e).clamp(-1.0, 1.0).acos();

    let semiperimeter_sq = (rho_s + rho_e + theta)
        * (rho_s + rho_e - theta)
        * (theta + rho_s - rho_e)
        * (theta - rho_s + rho_e);
    let triangle_term = semiperimeter_sq.max(0.0).sqrt() / 2.0;

    let a_overlap = rho_s * rho_s * alpha + rho_e * rho_e * beta - triangle_term;

    let a_sun = std::f64::consts::PI * rho_s * rho_s;
    (1.0 - a_overlap / a_sun).clamp(0.0, 1.0)
}

// =============================================================================
// ForceModel impl
// =============================================================================

impl ForceModel<Geocentric, GCRS> for CannonballSrp {
    #[inline]
    fn acceleration(
        &self,
        s: &OrbitState<Geocentric, GCRS>,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        let eph = ctx.require_ephemeris()?;
        let sun = sun_geocentric(eph, s.epoch_jd())?;

        // Vector from Sun to satellite (i.e. satellite minus Sun), GCRS km.
        let r_sun_sat = Displacement::<GCRS, Kilometer>::new(
            s.position.x().value() - sun.x().value(),
            s.position.y().value() - sun.y().value(),
            s.position.z().value() - sun.z().value(),
        );
        let r = r_sun_sat.magnitude().value();
        if r == 0.0 {
            // Satellite at Sun centre — return zero to avoid divide-by-zero.
            return Ok(Acceleration::<GCRS, AccelerationUnit>::new(0.0, 0.0, 0.0));
        }

        // Eclipse factor ν.
        let nu = match self.shadow {
            ShadowModel::None => 1.0,
            ShadowModel::Cylindrical => {
                let r_sat = [
                    s.position.x().value(),
                    s.position.y().value(),
                    s.position.z().value(),
                ];
                let r_sun_arr = [sun.x().value(), sun.y().value(), sun.z().value()];
                cylindrical_shadow_factor(r_sat, r_sun_arr, R_EARTH.value())
            }
            ShadowModel::Conical => {
                let r_sat = [
                    s.position.x().value(),
                    s.position.y().value(),
                    s.position.z().value(),
                ];
                let r_sun_arr = [sun.x().value(), sun.y().value(), sun.z().value()];
                conical_shadow_factor(r_sat, r_sun_arr, R_EARTH.value(), R_SUN_KM)
            }
        };

        let r2 = r * r;
        // N/m² · m²/kg = m/s²; divide by 1 000 to convert to km/s².
        let mag_km_s2 = nu
            * self.cr.value()
            * P0.value()
            * (AU_IN_KM * AU_IN_KM / r2)
            * self.area_to_mass.value()
            / 1_000.0;
        let inv_r = 1.0 / r;
        Ok(Acceleration::<GCRS, AccelerationUnit>::new(
            mag_km_s2 * r_sun_sat.x().value() * inv_r,
            mag_km_s2 * r_sun_sat.y().value() * inv_r,
            mag_km_s2 * r_sun_sat.z().value() * inv_r,
        ))
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::astro::dynamics::context::DynamicsContext;
    use crate::astro::dynamics::context::DynamicsContextBuilder;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::calculus::ephemeris::{AuPerDay, DynEphemeris, EphemerisError, Vsop87Ephemeris};
    use crate::coordinates::{
        cartesian::{Position as CartPosition, Velocity as CartVelocity},
        centers::{Barycentric, Geocentric as GeoCenter},
        frames::{EclipticMeanJ2000, GCRS},
    };
    use crate::qtty::{AstronomicalUnit, Kilometer};
    use crate::time::JulianDate;

    // ---- Stub ephemeris -------------------------------------------------------
    //
    // Returns a fixed Sun position at exactly 1 AU along the ecliptic +X axis
    // with Earth at the ecliptic origin.  After the ecliptic→GCRS rotation in
    // `sun_geocentric`, the Sun is ~1.496 × 10⁸ km from Earth, which yields a
    // physically correct SRP magnitude for tests that do not care about the
    // exact Sun direction.

    struct FixedSunEphemeris;

    impl DynEphemeris for FixedSunEphemeris {
        fn try_sun_barycentric(
            &self,
            _jd: JulianDate,
        ) -> Result<CartPosition<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError>
        {
            Ok(CartPosition::new(1.0_f64, 0.0_f64, 0.0_f64))
        }

        fn sun_barycentric(
            &self,
            _jd: JulianDate,
        ) -> CartPosition<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
            CartPosition::new(1.0_f64, 0.0_f64, 0.0_f64)
        }

        fn try_earth_barycentric(
            &self,
            _jd: JulianDate,
        ) -> Result<CartPosition<Barycentric, EclipticMeanJ2000, AstronomicalUnit>, EphemerisError>
        {
            Ok(CartPosition::new(0.0_f64, 0.0_f64, 0.0_f64))
        }

        fn earth_barycentric(
            &self,
            _jd: JulianDate,
        ) -> CartPosition<Barycentric, EclipticMeanJ2000, AstronomicalUnit> {
            CartPosition::new(0.0_f64, 0.0_f64, 0.0_f64)
        }

        fn try_earth_heliocentric(
            &self,
            _jd: JulianDate,
        ) -> Result<
            CartPosition<
                crate::coordinates::centers::Heliocentric,
                EclipticMeanJ2000,
                AstronomicalUnit,
            >,
            EphemerisError,
        > {
            Ok(CartPosition::new(0.0_f64, 0.0_f64, 0.0_f64))
        }

        fn earth_heliocentric(
            &self,
            _jd: JulianDate,
        ) -> CartPosition<
            crate::coordinates::centers::Heliocentric,
            EclipticMeanJ2000,
            AstronomicalUnit,
        > {
            CartPosition::new(0.0_f64, 0.0_f64, 0.0_f64)
        }

        fn try_earth_barycentric_velocity(
            &self,
            _jd: JulianDate,
        ) -> Result<CartVelocity<EclipticMeanJ2000, AuPerDay>, EphemerisError> {
            Ok(CartVelocity::new(0.0_f64, 0.0_f64, 0.0_f64))
        }

        fn earth_barycentric_velocity(
            &self,
            _jd: JulianDate,
        ) -> CartVelocity<EclipticMeanJ2000, AuPerDay> {
            CartVelocity::new(0.0_f64, 0.0_f64, 0.0_f64)
        }

        fn try_moon_geocentric(
            &self,
            _jd: JulianDate,
        ) -> Result<CartPosition<GeoCenter, EclipticMeanJ2000, Kilometer>, EphemerisError> {
            Ok(CartPosition::new(0.0_f64, 0.0_f64, 0.0_f64))
        }

        fn moon_geocentric(
            &self,
            _jd: JulianDate,
        ) -> CartPosition<GeoCenter, EclipticMeanJ2000, Kilometer> {
            CartPosition::new(0.0_f64, 0.0_f64, 0.0_f64)
        }
    }

    // ---- Helpers -------------------------------------------------------------

    fn ctx_stub() -> crate::astro::dynamics::context::DynamicsContext {
        DynamicsContextBuilder::new()
            .with_ephemeris(Arc::new(FixedSunEphemeris))
            .build()
    }

    fn ctx_vsop() -> crate::astro::dynamics::context::DynamicsContext {
        DynamicsContextBuilder::new()
            .with_ephemeris(Arc::new(Vsop87Ephemeris))
            .build()
    }

    fn leo() -> OrbitState {
        OrbitState::new_at_jd(
            JulianDate::JD_EPOCH_J2000_0,
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    // ---- ForceModel tests ----------------------------------------------------

    /// SRP at LEO, Cr = 1.5, A/m = 0.01 m²/kg, fully lit → ~7e-11 km/s² (1e-7 m/s²).
    #[test]
    fn srp_order_of_magnitude_at_leo() {
        let ctx = ctx_stub();
        let srp = CannonballSrp::new(SrpCoefficient::new(1.5), AreaToMass::new(0.01))
            .with_shadow(ShadowModel::None);
        let a = srp.acceleration(&leo(), &ctx).unwrap();
        let mag = a.magnitude().value();
        // Expected: ~6.8e-11 km/s²  (= 6.8e-8 m/s² ≈ 1e-7 m/s² order of magnitude).
        assert!(
            (1e-11..5e-10).contains(&mag),
            "SRP magnitude out of expected band: {mag:.3e} km/s²"
        );
    }

    /// A/m = 0 ⟹ zero acceleration regardless of ephemeris.
    #[test]
    fn srp_zero_when_area_is_zero() {
        let ctx = ctx_stub();
        let srp = CannonballSrp::new(SrpCoefficient::new(1.5), AreaToMass::new(0.0));
        let a = srp.acceleration(&leo(), &ctx).unwrap();
        assert_eq!(a.x().value(), 0.0);
        assert_eq!(a.y().value(), 0.0);
        assert_eq!(a.z().value(), 0.0);
    }

    /// `ShadowModel::None` must return non-zero acceleration even with the satellite
    /// on the anti-Sun side (no eclipse applied).
    #[test]
    fn shadow_model_none_ignores_eclipse() {
        let ctx = ctx_stub();
        // Put satellite in a position that would be shadowed under any eclipse model:
        // anti-Sun direction, well inside Earth radius cylinder.
        let anti_sun_state = OrbitState::new_at_jd(
            JulianDate::JD_EPOCH_J2000_0,
            Position::<GCRS>::new(-(R_EARTH.value() + 500.0), 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        let srp_no_eclipse = CannonballSrp::new(SrpCoefficient::new(1.5), AreaToMass::new(0.01))
            .with_shadow(ShadowModel::None);
        let a = srp_no_eclipse.acceleration(&anti_sun_state, &ctx).unwrap();
        // Must be non-zero (eclipse is ignored).
        let mag = a.magnitude().value();
        assert!(
            mag > 0.0,
            "ShadowModel::None should give non-zero SRP, got {mag}"
        );
    }

    /// Real ephemeris smoke-test: SRP is in the correct physical range at J2000.
    #[test]
    fn srp_real_ephemeris_magnitude() {
        let ctx = ctx_vsop();
        let srp = CannonballSrp::new(SrpCoefficient::new(1.5), AreaToMass::new(0.02))
            .with_shadow(ShadowModel::None);
        let s = OrbitState::new_at_jd(
            JulianDate::new(2_451_545.0),
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        );
        let a = srp.acceleration(&s, &ctx).unwrap();
        let mag = a.magnitude().value();
        assert!(
            (5e-11..5e-10).contains(&mag),
            "SRP magnitude out of expected band: {mag:.3e} km/s²"
        );
    }

    // ---- Cylindrical shadow tests --------------------------------------------

    /// Satellite on the anti-Sun axis, within the Earth shadow cylinder → ν = 0.
    #[test]
    fn cylindrical_shadow_on_axis_anti_sun() {
        // r_sun along +X (1 AU), satellite on −X at LEO altitude.
        let r_sun = [AU_IN_KM, 0.0, 0.0];
        let r_sat = [-(R_EARTH.value() + 100.0), 0.0, 0.0];
        let nu = cylindrical_shadow_factor(r_sat, r_sun, R_EARTH.value());
        assert_eq!(
            nu, 0.0,
            "satellite on anti-Sun axis should be fully shadowed"
        );
    }

    /// Satellite on the day side (+Sun direction) → ν = 1.
    #[test]
    fn cylindrical_shadow_day_side() {
        let r_sun = [AU_IN_KM, 0.0, 0.0];
        let r_sat = [R_EARTH.value() + 500.0, 0.0, 0.0];
        let nu = cylindrical_shadow_factor(r_sat, r_sun, R_EARTH.value());
        assert_eq!(nu, 1.0, "satellite on day side should be fully lit");
    }

    /// Satellite well outside the shadow cylinder (off-axis, anti-Sun) → ν = 1.
    #[test]
    fn cylindrical_shadow_outside_cylinder() {
        let r_sun = [AU_IN_KM, 0.0, 0.0];
        // Satellite perpendicular to the Sun-Earth axis; perp distance >> R_earth.
        let r_sat = [0.0, R_EARTH.value() * 2.0, 0.0];
        let nu = cylindrical_shadow_factor(r_sat, r_sun, R_EARTH.value());
        assert_eq!(nu, 1.0, "satellite outside cylinder should be fully lit");
    }

    // ---- Conical shadow tests ------------------------------------------------

    /// Full umbra: satellite on the anti-Sun axis at LEO → ν = 0.
    #[test]
    fn conical_shadow_full_umbra() {
        let r_sun_km = [AU_IN_KM, 0.0, 0.0];
        // Satellite directly behind Earth on the anti-Sun axis.
        let r_sat = [-(R_EARTH.value() + 500.0), 0.0, 0.0];
        let nu = conical_shadow_factor(r_sat, r_sun_km, R_EARTH.value(), R_SUN_KM);
        assert_eq!(
            nu, 0.0,
            "satellite in full umbra should have ν = 0; got {nu}"
        );
    }

    /// Full sunlight: satellite on the day side → ν = 1.
    #[test]
    fn conical_shadow_full_sun() {
        let r_sun_km = [AU_IN_KM, 0.0, 0.0];
        let r_sat = [R_EARTH.value() + 500.0, 0.0, 0.0];
        let nu = conical_shadow_factor(r_sat, r_sun_km, R_EARTH.value(), R_SUN_KM);
        assert_eq!(nu, 1.0, "satellite on day side should have ν = 1; got {nu}");
    }

    /// Partial penumbra: satellite off the umbra axis at high altitude → 0 < ν < 1.
    ///
    /// At 50,000 km altitude, ρ_e ≈ 0.128 rad.  By placing the satellite 6500 km
    /// off the anti-Sun axis the angular separation θ ≈ 0.130 rad falls inside
    /// the penumbra band [ρ_e − ρ_s, ρ_e + ρ_s] ≈ [0.123, 0.133] rad.
    #[test]
    fn conical_shadow_partial_penumbra() {
        let r_sun_km = [AU_IN_KM, 0.0, 0.0];
        // Satellite behind Earth (anti-Sun) at ~50,000 km altitude, offset 6500 km in Y.
        // θ ≈ 0.130 rad, rho_e ≈ 0.128 rad, rho_s ≈ 0.0047 rad → partial penumbra.
        let r_sat = [-50_000.0, 6_500.0, 0.0];
        let nu = conical_shadow_factor(r_sat, r_sun_km, R_EARTH.value(), R_SUN_KM);
        assert!(
            nu > 0.0 && nu < 1.0,
            "partial penumbra should give 0 < ν < 1; got {nu}"
        );
    }

    /// Annular eclipse: satellite behind Earth, Sun's apparent disk > Earth's.
    ///
    /// Placing the satellite 0.01 AU behind Earth (anti-Sun side, on axis) and the
    /// Sun at 0.1 AU from Earth: ρ_s ≈ 0.0423 rad > ρ_e ≈ 0.00426 rad, θ = 0 →
    /// annular eclipse → ν = 1 − (ρ_e / ρ_s)² ∈ (0, 1).
    #[test]
    fn conical_shadow_annular_eclipse() {
        let d_earth_km = 0.01 * AU_IN_KM; // satellite 0.01 AU from Earth, anti-Sun side
        let d_sun_km = 0.1 * AU_IN_KM; // Sun 0.1 AU from Earth in +X

        // Satellite on anti-Sun axis: r_sat = (-d_earth_km, 0, 0)
        // Sun in +X direction from Earth: r_sun = (d_sun_km, 0, 0)
        // → direction to Earth from satellite: (+X); direction to Sun from satellite: (+X)
        // → theta = 0 → annular eclipse condition.
        let r_sat = [-d_earth_km, 0.0, 0.0];
        let r_sun = [d_sun_km, 0.0, 0.0];

        let nu = conical_shadow_factor(r_sat, r_sun, R_EARTH.value(), R_SUN_KM);
        // nu = 1 − (rho_e / rho_s)^2 ≈ 1 − (0.00426 / 0.0423)^2 ≈ 0.990
        assert!(
            nu > 0.0 && nu < 1.0,
            "annular eclipse should give ν in (0, 1); got {nu}"
        );
        assert!(
            nu > 0.9,
            "annular eclipse: occulted fraction should be small; got nu={nu}"
        );
    }

    // ---- degenerate shadow edge cases ----------------------------------------

    #[test]
    fn cylindrical_shadow_zero_sun_returns_full_sun() {
        let nu = cylindrical_shadow_factor([7000.0, 0.0, 0.0], [0.0, 0.0, 0.0], 6378.0);
        assert_eq!(nu, 1.0, "degenerate sun at origin must return 1.0");
    }

    #[test]
    fn conical_shadow_satellite_at_earth_centre_returns_full_sun() {
        let nu = conical_shadow_factor(
            [0.0, 0.0, 0.0],
            [AU_IN_KM, 0.0, 0.0],
            R_EARTH.value(),
            R_SUN_KM,
        );
        assert_eq!(nu, 1.0, "satellite at Earth centre must return 1.0");
    }

    #[test]
    fn conical_shadow_satellite_at_sun_centre_returns_full_sun() {
        let r_sun = [AU_IN_KM, 0.0, 0.0];
        let nu = conical_shadow_factor(r_sun, r_sun, R_EARTH.value(), R_SUN_KM);
        assert_eq!(nu, 1.0, "satellite at Sun centre must return 1.0");
    }

    // ---- ForceModel::acceleration via Cylindrical/Conical models --------------

    #[test]
    fn srp_cylindrical_model_force_model() {
        let ctx = ctx_stub();
        let srp = CannonballSrp::new(SrpCoefficient::new(1.5), AreaToMass::new(0.01))
            .with_shadow(ShadowModel::Cylindrical);
        let a = srp.acceleration(&leo(), &ctx).unwrap();
        let mag = a.magnitude().value();
        assert!(
            mag > 0.0,
            "cylindrical-shadow SRP at day-side LEO must be non-zero"
        );
    }

    #[test]
    fn srp_conical_model_force_model() {
        let ctx = ctx_stub();
        let srp = CannonballSrp::new(SrpCoefficient::new(1.5), AreaToMass::new(0.01))
            .with_shadow(ShadowModel::Conical);
        let a = srp.acceleration(&leo(), &ctx).unwrap();
        let mag = a.magnitude().value();
        assert!(
            mag > 0.0,
            "conical-shadow SRP at day-side LEO must be non-zero"
        );
    }

    #[test]
    fn srp_no_ephemeris_returns_error() {
        let ctx = DynamicsContext::empty();
        let srp = CannonballSrp::new(SrpCoefficient::new(1.5), AreaToMass::new(0.01));
        let result = srp.acceleration(&leo(), &ctx);
        assert!(result.is_err(), "SRP without ephemeris must return error");
    }

    /// Tests the anti-Sun, outside-shadow-cylinder path (`1.0` branch in
    /// `cylindrical_shadow_factor`).  Satellite is behind Earth relative to
    /// the Sun but laterally displaced far enough to miss the shadow cylinder.
    #[test]
    fn cylindrical_shadow_antisun_outside_cylinder() {
        // Sun along +X.  Satellite at −8000 km in X (anti-Sun), 8000 km in Y
        // (outside the ~6378 km radius shadow cylinder).
        let nu = cylindrical_shadow_factor(
            [-8_000.0, 8_000.0, 0.0],
            [149_597_870.7, 0.0, 0.0],
            6_378.137,
        );
        assert_eq!(nu, 1.0, "anti-Sun but outside cylinder should be fully lit");
    }

    /// Exercises every non-primary method on `FixedSunEphemeris` so that the
    /// stub implementation lines are counted as covered.
    #[test]
    fn stub_ephemeris_all_methods_reachable() {
        let eph = FixedSunEphemeris;
        let jd = JulianDate::JD_EPOCH_J2000_0;

        let _ = eph.sun_barycentric(jd);
        let _ = eph.earth_barycentric(jd);
        let _ = eph.try_earth_heliocentric(jd).unwrap();
        let _ = eph.earth_heliocentric(jd);
        let _ = eph.try_earth_barycentric_velocity(jd).unwrap();
        let _ = eph.earth_barycentric_velocity(jd);
        let _ = eph.try_moon_geocentric(jd).unwrap();
        let _ = eph.moon_geocentric(jd);
    }
}
