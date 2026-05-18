// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Third-body perturbation force models — [`ThirdBody`] and [`ThirdBodyProvider`].
//!
//! ## Overview
//!
//! [`ThirdBody`] is a provider-backed, extensible third-body force model.
//! Any number of perturbing bodies can be registered via [`ThirdBodyProvider`]
//! implementations; built-in providers cover the Sun ([`SunPerturbation`])
//! and Moon ([`MoonPerturbation`]).
//!
//! ## Battin's Numerically-Stable Formula
//!
//! The naive direct form of the third-body acceleration,
//!
//! ```text
//! a = μ_b · [ (r_b − r_s) / |r_b − r_s|³  −  r_b / |r_b|³ ]
//! ```
//!
//! suffers **catastrophic cancellation** when the spacecraft is far from the
//! third body and the two terms nearly cancel (|r_b| ≫ |r_b − r_s|).  In
//! LEO the Sun is ~1.5 × 10⁸ km away while |r_s| ~ 7 × 10³ km; the two
//! subtracted vectors agree to ~5 significant figures, wasting ~5 decimal
//! digits of floating-point precision.
//!
//! **Battin's reformulation** (§10.2, also Vallado §8.7) avoids this by
//! factoring out |r_b|³ and using the identity |r_b−r_s|³ = |r_b|³·(1+q)^{3/2}:
//!
//! ```text
//! q    =  r_s · (r_s − 2 r_b)  /  |r_b|²
//! F(q) =  q · (3 + 3q + q²)  /  (1 + (1 + q)^{3/2})        [= (1+q)^{3/2} − 1]
//! a    =  −μ_b · (r_b · F(q) + r_s)  /  (|r_b|³ · (1 + F(q)))
//! ```
//!
//! The denominator `|r_b|³·(1+F(q)) = |r_b−r_s|³` is computed from the
//! scalar `F(q)` without subtracting large vectors.  The numerator
//! `r_b·F(q) + r_s` involves km-scale quantities of similar magnitude with
//! no catastrophic cancellation.
//!
//! ## Frame / center assumptions
//!
//! Frame: **GCRS** (geocentric celestial reference system, equatorial).
//! Center: **Geocentric**.
//!
//! Both the satellite position and the third-body positions must be expressed
//! in GCRS at the same epoch.  The built-in [`SunPerturbation`] and
//! [`MoonPerturbation`] providers rotate ephemeris output from ecliptic to
//! GCRS using [`ecliptic_of_date_to_mean_equatorial_matrix`] evaluated at
//! J2000 (the conventional sidereal reference epoch).
//!
//! ## References
//!
//! * Battin, R. H., *An Introduction to the Mathematics and Methods of
//!   Astrodynamics* (revised ed.), §10.2 — derivation of the F function.
//! * Vallado, D. A., *Fundamentals of Astrodynamics and Applications* (4th
//!   ed.), §8.7 — third-body perturbations.
//! * Montenbruck & Gill, *Satellite Orbits*, §3.2.

use std::sync::Arc;

use affn::cartesian::Displacement;

use crate::astro::dynamics::context::DynamicsContext;
use crate::astro::dynamics::errors::DynamicsError;
use crate::astro::dynamics::state::{Acceleration, AccelerationUnit, OrbitState};
use crate::astro::precession::ecliptic_of_date_to_mean_equatorial_matrix;
use crate::calculus::ephemeris::DynEphemeris;
use crate::coordinates::centers::Geocentric;
use crate::coordinates::frames::{EclipticMeanJ2000, GCRS};
use crate::qtty::{AstronomicalUnit, Kilometer};
use crate::time::JulianDate;

use super::traits::{ForceModel, GravitationalParameter, GM_MOON, GM_SUN};

// ---- Battin formula ---------------------------------------------------------

/// Numerically-stable third-body perturbation acceleration (Battin §10.2).
///
/// ## Inputs
///
/// * `r_s` — satellite geocentric position (km).
/// * `r_b` — third-body geocentric position (km).
/// * `mu`  — standard gravitational parameter μ = GM of the third body (km³/s²).
///
/// ## Derivation
///
/// Starting from the two-body perturbation:
///
/// ```text
/// a = μ · [ (r_b − r_s) / |r_b − r_s|³  −  r_b / |r_b|³ ]
/// ```
///
/// Note that |r_b − r_s|² = |r_b|² · (1 + q) where:
///
/// ```text
/// q = r_s · (r_s − 2 r_b) / |r_b|²
/// ```
///
/// so |r_b − r_s|³ = |r_b|³ · (1+q)^{3/2}.  Define the Battin F function:
///
/// ```text
/// F(q) = q · (3 + 3q + q²) / (1 + (1+q)^{3/2})
/// ```
///
/// which satisfies F(q) = (1+q)^{3/2} − 1 algebraically.  Substituting
/// and rearranging gives the stable form:
///
/// ```text
/// a = −μ / (|r_b|³ · (1 + F(q)))  ·  (r_b · F(q) + r_s)
/// ```
///
/// The denominator `|r_b|³ · (1 + F(q))` equals `|r_b − r_s|³` by
/// construction, so this is algebraically identical to the naive form.  The
/// numerator `r_b · F(q) + r_s` involves quantities of similar magnitude
/// (both ~km-scale) and avoids the km/s² cancellation of the naive form.
///
/// Returns `[0; 3]` if the third-body is at the origin (degenerate geometry).
fn battin_third_body_accel(r_s: [f64; 3], r_b: [f64; 3], mu: f64) -> [f64; 3] {
    let r_b_sq = r_b[0] * r_b[0] + r_b[1] * r_b[1] + r_b[2] * r_b[2];
    let r_b_norm = r_b_sq.sqrt();
    if r_b_norm == 0.0 {
        return [0.0; 3];
    }

    // q = r_s · (r_s − 2·r_b) / |r_b|²
    let q = (r_s[0] * (r_s[0] - 2.0 * r_b[0])
        + r_s[1] * (r_s[1] - 2.0 * r_b[1])
        + r_s[2] * (r_s[2] - 2.0 * r_b[2]))
        / r_b_sq;

    // F(q) = q · (3 + 3q + q²) / (1 + (1+q)^{3/2})
    //
    // Algebraically equal to (1+q)^{3/2} − 1, but computed via this ratio
    // to stay well-conditioned at small |q| (typical LEO case: q ≈ −1e-4).
    let fq = q * (3.0 + q * (3.0 + q)) / (1.0 + (1.0 + q).powf(1.5));

    // a = −μ / (|r_b|³ · (1 + F(q))) · (r_b · F(q) + r_s)
    //
    // Denominator = |r_b − r_s|³ (computed without subtracting large vectors).
    // Numerator terms r_b·F(q) and r_s are both km-scale; no large cancellation.
    let r_b3_times_1pfq = r_b_norm * r_b_sq * (1.0 + fq);
    if r_b3_times_1pfq == 0.0 {
        return [0.0; 3];
    }

    let scale = -mu / r_b3_times_1pfq;
    [
        scale * (r_b[0] * fq + r_s[0]),
        scale * (r_b[1] * fq + r_s[1]),
        scale * (r_b[2] * fq + r_s[2]),
    ]
}

// ---- Position helpers -------------------------------------------------------

/// Geocentric Sun displacement in GCRS (km), from an ephemeris provider.
pub(super) fn sun_geocentric(
    eph: &Arc<dyn DynEphemeris + Send + Sync>,
    jd: JulianDate,
) -> Result<Displacement<GCRS, Kilometer>, DynamicsError> {
    let sun_b = eph
        .try_sun_barycentric(jd)
        .map_err(|e| DynamicsError::EphemerisUnavailable {
            body: "Sun",
            source: Some(Box::new(e)),
        })?;
    let earth_b =
        eph.try_earth_barycentric(jd)
            .map_err(|e| DynamicsError::EphemerisUnavailable {
                body: "Earth",
                source: Some(Box::new(e)),
            })?;
    let d_ecl: Displacement<EclipticMeanJ2000, AstronomicalUnit> = sun_b - earth_b;
    let rot = ecliptic_of_date_to_mean_equatorial_matrix(JulianDate::JD_EPOCH_J2000_0);
    let d_gcrs: Displacement<GCRS, AstronomicalUnit> = rot.apply_vec(d_ecl);
    Ok(d_gcrs.to_unit::<Kilometer>())
}

/// Geocentric Moon displacement in GCRS (km), from an ephemeris provider.
pub(super) fn moon_geocentric(
    eph: &Arc<dyn DynEphemeris + Send + Sync>,
    jd: JulianDate,
) -> Result<Displacement<GCRS, Kilometer>, DynamicsError> {
    let m = eph
        .try_moon_geocentric(jd)
        .map_err(|e| DynamicsError::EphemerisUnavailable {
            body: "Moon",
            source: Some(Box::new(e)),
        })?;
    let m_ecl = Displacement::<EclipticMeanJ2000, Kilometer>::new(
        m.x().value(),
        m.y().value(),
        m.z().value(),
    );
    let rot = ecliptic_of_date_to_mean_equatorial_matrix(JulianDate::JD_EPOCH_J2000_0);
    Ok(rot.apply_vec(m_ecl))
}

// ---- ThirdBodyProvider trait ------------------------------------------------

/// Supplies the position and gravitational parameter of a perturbing body.
///
/// Implement this trait to add custom third-body perturbations beyond the
/// built-in [`SunPerturbation`] and [`MoonPerturbation`].  Inject instances
/// via [`ThirdBody::with`].
///
/// # Thread safety
///
/// Implementations must be `Send + Sync` so that a [`ThirdBody`] model can be
/// shared across threads without wrapping.
pub trait ThirdBodyProvider: Send + Sync {
    /// Human-readable body name for diagnostics and logging.
    fn name(&self) -> &'static str;

    /// Standard gravitational parameter μ = GM (km³/s²) of the body.
    fn gm(&self) -> GravitationalParameter;

    /// Geocentric GCRS position of the body at `epoch` (km), fetched from `eph`.
    ///
    /// The returned [`Displacement`] is measured from the Earth's geocenter to
    /// the body, expressed in the GCRS frame.
    fn position_relative_to_origin(
        &self,
        eph: &Arc<dyn DynEphemeris + Send + Sync>,
        epoch: JulianDate,
    ) -> Result<Displacement<GCRS, Kilometer>, DynamicsError>;
}

// ---- Built-in providers -----------------------------------------------------

/// Sun perturbation provider.
///
/// Fetches the Sun's geocentric GCRS position using barycentric ephemeris
/// subtraction (`r_sun_bary − r_earth_bary`) rotated from ecliptic to GCRS.
/// GM = [`GM_SUN`].
#[derive(Debug, Clone, Copy, Default)]
pub struct SunPerturbation;

impl ThirdBodyProvider for SunPerturbation {
    #[inline]
    fn name(&self) -> &'static str {
        "Sun"
    }

    #[inline]
    fn gm(&self) -> GravitationalParameter {
        GM_SUN
    }

    #[inline]
    fn position_relative_to_origin(
        &self,
        eph: &Arc<dyn DynEphemeris + Send + Sync>,
        epoch: JulianDate,
    ) -> Result<Displacement<GCRS, Kilometer>, DynamicsError> {
        sun_geocentric(eph, epoch)
    }
}

/// Moon perturbation provider.
///
/// Fetches the Moon's geocentric GCRS position from the ephemeris, rotated
/// from ecliptic to GCRS.
/// GM = [`GM_MOON`].
#[derive(Debug, Clone, Copy, Default)]
pub struct MoonPerturbation;

impl ThirdBodyProvider for MoonPerturbation {
    #[inline]
    fn name(&self) -> &'static str {
        "Moon"
    }

    #[inline]
    fn gm(&self) -> GravitationalParameter {
        GM_MOON
    }

    #[inline]
    fn position_relative_to_origin(
        &self,
        eph: &Arc<dyn DynEphemeris + Send + Sync>,
        epoch: JulianDate,
    ) -> Result<Displacement<GCRS, Kilometer>, DynamicsError> {
        moon_geocentric(eph, epoch)
    }
}

// ---- ThirdBody struct -------------------------------------------------------

/// Provider-backed, extensible third-body perturbation force model.
///
/// Accumulates accelerations from an arbitrary list of [`ThirdBodyProvider`]s.
/// Each provider contributes an acceleration computed with Battin's numerically
/// stable formula (see module docs).
///
/// ## Builder
///
/// ```rust,ignore
/// use siderust::astro::dynamics::forces::ThirdBody;
///
/// // Sun + Moon (most common case):
/// let model = ThirdBody::sun_and_moon();
///
/// // Sun only:
/// let model = ThirdBody::new().with_sun();
///
/// // Custom provider:
/// let model = ThirdBody::new().with(Box::new(my_provider));
/// ```
///
/// ## Empty list behaviour
///
/// If no bodies have been registered, [`ForceModel::acceleration`] returns
/// zero acceleration **without** requiring an ephemeris provider.  This
/// ensures that an empty `ThirdBody` model can be inserted into a composite
/// without causing spurious `EphemerisUnavailable` errors.
pub struct ThirdBody {
    bodies: Vec<Box<dyn ThirdBodyProvider>>,
}

impl ThirdBody {
    /// Create an empty model with no perturbing bodies.
    pub fn new() -> Self {
        Self { bodies: Vec::new() }
    }

    /// Add the Sun as a perturbing body.
    pub fn with_sun(mut self) -> Self {
        self.bodies.push(Box::new(SunPerturbation));
        self
    }

    /// Add the Moon as a perturbing body.
    pub fn with_moon(mut self) -> Self {
        self.bodies.push(Box::new(MoonPerturbation));
        self
    }

    /// Add a custom perturbing body.
    pub fn with(mut self, body: Box<dyn ThirdBodyProvider>) -> Self {
        self.bodies.push(body);
        self
    }

    /// Convenience constructor: Sun + Moon perturbations.
    ///
    /// This is the typical LEO/GEO configuration and replaces the old
    /// `ThirdBodySunMoon` struct.
    pub fn sun_and_moon() -> Self {
        Self::new().with_sun().with_moon()
    }
}

impl Default for ThirdBody {
    fn default() -> Self {
        Self::new()
    }
}

// ---- ForceModel impl --------------------------------------------------------

impl ForceModel<Geocentric, GCRS> for ThirdBody {
    fn acceleration(
        &self,
        s: &OrbitState<Geocentric, GCRS>,
        ctx: &DynamicsContext,
    ) -> Result<Acceleration<GCRS, AccelerationUnit>, DynamicsError> {
        // Short-circuit: no bodies registered → zero perturbation, no ephemeris needed.
        if self.bodies.is_empty() {
            return Ok(Acceleration::new(0.0, 0.0, 0.0));
        }

        let eph = ctx.require_ephemeris()?;
        let jd = s.epoch_jd();
        let r_s = [
            s.position.x().value(),
            s.position.y().value(),
            s.position.z().value(),
        ];

        let mut total = [0.0_f64; 3];
        for body in &self.bodies {
            let r_b_disp = body.position_relative_to_origin(eph, jd)?;
            let r_b = [
                r_b_disp.x().value(),
                r_b_disp.y().value(),
                r_b_disp.z().value(),
            ];
            let a = battin_third_body_accel(r_s, r_b, body.gm().value());
            total[0] += a[0];
            total[1] += a[1];
            total[2] += a[2];
        }

        Ok(Acceleration::new(total[0], total[1], total[2]))
    }
}

// ---- Backward-compatibility alias -------------------------------------------

/// Deprecated: use [`ThirdBody::sun_and_moon()`] or [`ThirdBody`] directly.
///
/// `ThirdBodySunMoon` was a zero-sized struct representing fixed Sun + Moon
/// perturbations.  The replacement [`ThirdBody`] is a provider-backed model
/// that supports arbitrary body lists.  For a drop-in replacement, call
/// `ThirdBody::sun_and_moon()`.
#[deprecated(
    since = "0.1.0",
    note = "use `ThirdBody::sun_and_moon()` or `ThirdBody` with explicit providers"
)]
pub type ThirdBodySunMoon = ThirdBody;

// ---- Tests ------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::astro::dynamics::context::DynamicsContextBuilder;
    use crate::astro::dynamics::state::OrbitState;
    use crate::astro::dynamics::{Position, Velocity};
    use crate::calculus::ephemeris::Vsop87Ephemeris;
    use crate::coordinates::frames::GCRS;
    use crate::time::JulianDate;

    // ---- Helpers ------------------------------------------------------------

    fn leo_at(epoch: JulianDate) -> OrbitState {
        OrbitState::new_at_jd(
            epoch,
            Position::<GCRS>::new(7000.0, 0.0, 0.0),
            Velocity::<GCRS>::new(0.0, 7.5, 0.0),
        )
    }

    fn ctx_with_eph() -> crate::astro::dynamics::context::DynamicsContext {
        DynamicsContextBuilder::new()
            .with_ephemeris(Arc::new(Vsop87Ephemeris))
            .build()
    }

    const J2000: JulianDate = JulianDate::JD_EPOCH_J2000_0;

    // ---- Battin formula vs naive (regression / sanity) ----------------------

    /// Naive third-body acceleration (direct formula, suffers catastrophic cancellation).
    fn naive_accel(r_s: [f64; 3], r_b: [f64; 3], mu: f64) -> [f64; 3] {
        let d = [r_b[0] - r_s[0], r_b[1] - r_s[1], r_b[2] - r_s[2]];
        let d_n = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let r_b_n = (r_b[0] * r_b[0] + r_b[1] * r_b[1] + r_b[2] * r_b[2]).sqrt();
        if d_n == 0.0 || r_b_n == 0.0 {
            return [0.0; 3];
        }
        let d3 = d_n * d_n * d_n;
        let r3 = r_b_n * r_b_n * r_b_n;
        [
            mu * (d[0] / d3 - r_b[0] / r3),
            mu * (d[1] / d3 - r_b[1] / r3),
            mu * (d[2] / d3 - r_b[2] / r3),
        ]
    }

    #[test]
    fn battin_vs_naive_typical_leo() {
        // Typical LEO position, Sun ~1 AU away – non-pathological geometry.
        // Both formulas should agree to high precision.
        let r_s = [7_000.0_f64, 0.0, 0.0]; // km, LEO
        let r_b = [1.496e8_f64, 0.0, 0.0]; // km, ~1 AU

        let mu = GM_SUN.value();

        let battin = battin_third_body_accel(r_s, r_b, mu);
        let naive = naive_accel(r_s, r_b, mu);

        let mag_b = (battin[0].powi(2) + battin[1].powi(2) + battin[2].powi(2)).sqrt();
        let mag_n = (naive[0].powi(2) + naive[1].powi(2) + naive[2].powi(2)).sqrt();

        // Both give finite, non-zero result.
        assert!(mag_b > 0.0 && mag_b.is_finite(), "Battin mag: {mag_b}");
        assert!(mag_n > 0.0 && mag_n.is_finite(), "Naive mag: {mag_n}");

        // Physical plausibility: Sun perturbation at LEO ~5e-10 km/s²
        assert!(
            mag_n > 1e-12 && mag_n < 1e-6,
            "naive result out of expected range: {mag_n:.4e} km/s²"
        );

        // Relative error between the two formulations should be very small.
        let rel_err = ((mag_b - mag_n) / mag_n).abs();
        println!("battin={mag_b:.6e}, naive={mag_n:.6e}, rel_err={rel_err:.2e}");
        assert!(
            rel_err < 1e-10,
            "Battin and naive disagree: rel_err={rel_err:.2e}"
        );
    }

    #[test]
    fn battin_near_collinear_pathological() {
        // Satellite almost collinear with the Sun (very close to the Sun-Earth
        // line), where the naive formula suffers catastrophic cancellation.
        // In this regime the two subtracted terms nearly cancel, losing ~5
        // significant digits; Battin's form retains full precision.
        let r_s = [7_000.0_f64, 1.0e-6, 0.0]; // nearly on the Sun-Earth axis
        let r_b = [1.496e8_f64, 0.0, 0.0];

        let mu = GM_SUN.value();

        let battin = battin_third_body_accel(r_s, r_b, mu);
        let mag = (battin[0].powi(2) + battin[1].powi(2) + battin[2].powi(2)).sqrt();

        // The result must be finite and physically plausible (~5e-10 km/s²).
        assert!(
            mag.is_finite() && mag > 0.0,
            "near-collinear: non-finite result {mag}"
        );
        assert!(
            mag < 1e-5,
            "near-collinear: unrealistically large {mag} km/s²"
        );
        println!("near-collinear |a_battin| = {mag:.4e} km/s²");
    }

    // ---- ForceModel on real ephemeris ----------------------------------------

    #[test]
    fn sun_and_moon_nonzero_acceleration() {
        let ctx = ctx_with_eph();
        let s = leo_at(J2000);
        let a = ThirdBody::sun_and_moon().acceleration(&s, &ctx).unwrap();
        let mag = a.magnitude().value();
        assert!(mag > 0.0, "third-body acceleration should be nonzero");
        assert!(
            mag < 1e-4,
            "third-body acceleration unrealistically large: {mag} km/s²"
        );
        println!("|a_sun_moon| = {mag:.4e} km/s²");
    }

    #[test]
    fn sun_only_acceleration_nonzero() {
        let ctx = ctx_with_eph();
        let s = leo_at(J2000);
        let a = ThirdBody::new().with_sun().acceleration(&s, &ctx).unwrap();
        let mag = a.magnitude().value();
        assert!(mag > 0.0, "Sun-only perturbation should be nonzero");
    }

    #[test]
    fn moon_only_acceleration_nonzero() {
        let ctx = ctx_with_eph();
        let s = leo_at(J2000);
        let a = ThirdBody::new().with_moon().acceleration(&s, &ctx).unwrap();
        let mag = a.magnitude().value();
        assert!(mag > 0.0, "Moon-only perturbation should be nonzero");
    }

    // ---- Empty body list ----------------------------------------------------

    #[test]
    fn empty_body_list_returns_zero_without_ephemeris() {
        // No ephemeris in context — empty ThirdBody must not request it.
        let ctx = crate::astro::dynamics::context::DynamicsContext::empty();
        let s = leo_at(J2000);
        let a = ThirdBody::new().acceleration(&s, &ctx).unwrap();
        assert_eq!(a.x().value(), 0.0);
        assert_eq!(a.y().value(), 0.0);
        assert_eq!(a.z().value(), 0.0);
    }

    // ---- Missing ephemeris error ---------------------------------------------

    #[test]
    fn missing_ephemeris_returns_error() {
        let ctx = crate::astro::dynamics::context::DynamicsContext::empty();
        let s = leo_at(J2000);
        let result = ThirdBody::sun_and_moon().acceleration(&s, &ctx);
        assert!(
            matches!(result, Err(DynamicsError::EphemerisUnavailable { .. })),
            "expected EphemerisUnavailable, got: {result:?}"
        );
    }

    // ---- Builder composition -------------------------------------------------

    #[test]
    fn builder_composition_sun_moon_consistent() {
        let ctx = ctx_with_eph();
        let s = leo_at(J2000);

        let a_combined = ThirdBody::sun_and_moon().acceleration(&s, &ctx).unwrap();
        let a_sun = ThirdBody::new().with_sun().acceleration(&s, &ctx).unwrap();
        let a_moon = ThirdBody::new().with_moon().acceleration(&s, &ctx).unwrap();

        // ThirdBody::sun_and_moon() == with_sun().with_moon() → linearity check.
        let sum_x = a_sun.x().value() + a_moon.x().value();
        let sum_y = a_sun.y().value() + a_moon.y().value();
        let sum_z = a_sun.z().value() + a_moon.z().value();

        let eps = 1e-20;
        assert!((a_combined.x().value() - sum_x).abs() < eps);
        assert!((a_combined.y().value() - sum_y).abs() < eps);
        assert!((a_combined.z().value() - sum_z).abs() < eps);
    }

    // ---- Custom provider via ThirdBodyProvider trait -------------------------

    struct StubBody {
        pos_km: [f64; 3],
        gm_val: GravitationalParameter,
    }

    impl ThirdBodyProvider for StubBody {
        fn name(&self) -> &'static str {
            "Stub"
        }
        fn gm(&self) -> GravitationalParameter {
            self.gm_val
        }
        fn position_relative_to_origin(
            &self,
            _eph: &Arc<dyn DynEphemeris + Send + Sync>,
            _epoch: JulianDate,
        ) -> Result<Displacement<GCRS, Kilometer>, DynamicsError> {
            Ok(Displacement::new(
                self.pos_km[0],
                self.pos_km[1],
                self.pos_km[2],
            ))
        }
    }

    // SAFETY: StubBody fields ([f64; 3] and GravitationalParameter) are both
    // Send + Sync; the unsafe impls are needed because Rust can't infer them
    // automatically for a struct defined inside a test module.
    unsafe impl Send for StubBody {}
    unsafe impl Sync for StubBody {}

    #[test]
    fn custom_provider_via_with() {
        let ctx_eph = ctx_with_eph();
        let s = leo_at(J2000);

        let model = ThirdBody::new().with(Box::new(StubBody {
            pos_km: [1.496e8, 0.0, 0.0],
            gm_val: GM_SUN,
        }));
        let a = model.acceleration(&s, &ctx_eph).unwrap();
        let mag = a.magnitude().value();
        assert!(
            mag > 0.0,
            "custom provider: expected nonzero accel, got {mag}"
        );
    }

    #[test]
    fn default_is_empty() {
        let ctx = crate::astro::dynamics::context::DynamicsContext::empty();
        let s = leo_at(J2000);
        let a = ThirdBody::default().acceleration(&s, &ctx).unwrap();
        assert_eq!(a.x().value(), 0.0);
    }
}
