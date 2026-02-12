// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Comprehensive tests for the [`calculus::ephemeris`] module.
//!
//! This test file aims to provide comprehensive coverage of the ephemeris module,
//! testing all backends (Vsop87Ephemeris, De440Ephemeris, De441Ephemeris) and
//! all trait methods across multiple epochs.

use qtty::*;
use siderust::calculus::ephemeris::{Ephemeris, Vsop87Ephemeris};
use siderust::time::JulianDate;

// ============================================================================
// Helper Functions
// ============================================================================

/// Build a [`JulianDate`] from a Julian Day Number (for clarity in tests).
fn jd_from_value(jd: f64) -> JulianDate {
    JulianDate::new(jd)
}

/// Standard test epoch: J2000.0
fn j2000() -> JulianDate {
    JulianDate::J2000
}

/// Test epoch: 2020-01-01 12:00 TT
fn epoch_2020() -> JulianDate {
    // JD 2458849.5 corresponds to 2020-01-01 00:00 UTC
    // Using 12:00 TT = JD 2458850.0
    jd_from_value(2458850.0)
}

/// Test epoch: 2026-02-15 12:00 TT
fn epoch_2026() -> JulianDate {
    // JD 2461079.0 approximately
    jd_from_value(2461079.0)
}

/// Epoch far in the past (for range testing)
fn epoch_1900() -> JulianDate {
    // JD 2415020.0 corresponds to 1900-01-01 12:00 TT
    jd_from_value(2415020.0)
}

/// Epoch far in the future (for range testing)
fn epoch_2100() -> JulianDate {
    // JD 2488069.5 corresponds to 2100-01-01 12:00 TT
    jd_from_value(2488069.5)
}

// ============================================================================
// VSOP87 Ephemeris Tests
// ============================================================================

mod vsop87_tests {
    use super::*;

    // -------------------------------------------------------------------------
    // sun_barycentric tests
    // -------------------------------------------------------------------------

    #[test]
    fn sun_barycentric_at_j2000() {
        let jd = j2000();
        let sun = Vsop87Ephemeris::sun_barycentric(jd);
        let pos = sun.get_position();

        // Sun should be very close to SSB at J2000 (within a few solar radii)
        // Sun-SSB offset is typically < 0.01 AU
        let dist = pos.distance();
        assert!(
            dist < AstronomicalUnits::new(0.02),
            "Sun should be close to SSB at J2000, got distance: {dist:?}"
        );

        // Verify components are reasonable (non-NaN, finite)
        assert!(pos.x().value().is_finite(), "x should be finite");
        assert!(pos.y().value().is_finite(), "y should be finite");
        assert!(pos.z().value().is_finite(), "z should be finite");
    }

    #[test]
    fn sun_barycentric_at_epoch_2020() {
        let jd = epoch_2020();
        let sun = Vsop87Ephemeris::sun_barycentric(jd);
        let pos = sun.get_position();

        // Sun-SSB offset varies but typically < 0.02 AU
        let dist = pos.distance();
        assert!(
            dist < AstronomicalUnits::new(0.03),
            "Sun should be reasonably close to SSB, got distance: {dist:?}"
        );
    }

    #[test]
    fn sun_barycentric_at_epoch_2026() {
        let jd = epoch_2026();
        let sun = Vsop87Ephemeris::sun_barycentric(jd);
        let pos = sun.get_position();

        let dist = pos.distance();
        assert!(
            dist < AstronomicalUnits::new(0.03),
            "Sun should be reasonably close to SSB at 2026"
        );
    }

    #[test]
    fn sun_barycentric_at_epoch_1900() {
        let jd = epoch_1900();
        let sun = Vsop87Ephemeris::sun_barycentric(jd);
        let pos = sun.get_position();

        // Even at extreme historical dates, Sun-SSB should be bounded
        let dist = pos.distance();
        assert!(
            dist < AstronomicalUnits::new(0.03),
            "Sun should be bounded even in 1900"
        );
    }

    #[test]
    fn sun_barycentric_at_epoch_2100() {
        let jd = epoch_2100();
        let sun = Vsop87Ephemeris::sun_barycentric(jd);
        let pos = sun.get_position();

        let dist = pos.distance();
        assert!(
            dist < AstronomicalUnits::new(0.03),
            "Sun should be bounded even in 2100"
        );
    }

    // -------------------------------------------------------------------------
    // earth_barycentric tests
    // -------------------------------------------------------------------------

    #[test]
    fn earth_barycentric_at_j2000() {
        let jd = j2000();
        let earth = Vsop87Ephemeris::earth_barycentric(jd);
        let pos = earth.get_position();

        // Earth should be ~1 AU from SSB
        let dist = pos.distance();
        assert!(
            dist > AstronomicalUnits::new(0.98),
            "Earth should be > 0.98 AU from SSB"
        );
        assert!(
            dist < AstronomicalUnits::new(1.02),
            "Earth should be < 1.02 AU from SSB"
        );
    }

    #[test]
    fn earth_barycentric_at_epoch_2020() {
        let jd = epoch_2020();
        let earth = Vsop87Ephemeris::earth_barycentric(jd);
        let pos = earth.get_position();

        let dist = pos.distance();
        assert!(
            dist > AstronomicalUnits::new(0.98) && dist < AstronomicalUnits::new(1.02),
            "Earth should be ~1 AU from SSB at 2020, got {dist:?}"
        );
    }

    #[test]
    fn earth_barycentric_at_epoch_2026() {
        let jd = epoch_2026();
        let earth = Vsop87Ephemeris::earth_barycentric(jd);
        let pos = earth.get_position();

        let dist = pos.distance();
        assert!(
            dist > AstronomicalUnits::new(0.98) && dist < AstronomicalUnits::new(1.02),
            "Earth should be ~1 AU from SSB at 2026"
        );
    }

    // -------------------------------------------------------------------------
    // earth_heliocentric tests
    // -------------------------------------------------------------------------

    #[test]
    fn earth_heliocentric_at_j2000() {
        let jd = j2000();
        let earth = Vsop87Ephemeris::earth_heliocentric(jd);
        let pos = earth.get_position();

        // Earth-Sun distance should be ~1 AU
        let dist = pos.distance();
        assert!(
            dist > AstronomicalUnits::new(0.98) && dist < AstronomicalUnits::new(1.02),
            "Earth heliocentric distance should be ~1 AU, got {dist:?}"
        );
    }

    #[test]
    fn earth_heliocentric_at_perihelion() {
        // Perihelion ~Jan 3 each year, Earth is closest to Sun
        // JD 2459946.0 ~ 2023-01-03
        let jd = jd_from_value(2459946.0);
        let earth = Vsop87Ephemeris::earth_heliocentric(jd);
        let pos = earth.get_position();

        let dist = pos.distance();
        // At perihelion, distance should be ~0.983 AU
        assert!(
            dist < AstronomicalUnits::new(0.99),
            "Earth at perihelion should be < 0.99 AU from Sun, got {dist:?}"
        );
    }

    #[test]
    fn earth_heliocentric_at_aphelion() {
        // Aphelion ~Jul 4 each year, Earth is farthest from Sun
        // JD 2460128.0 ~ 2023-07-04
        let jd = jd_from_value(2460128.0);
        let earth = Vsop87Ephemeris::earth_heliocentric(jd);
        let pos = earth.get_position();

        let dist = pos.distance();
        // At aphelion, distance should be ~1.017 AU
        assert!(
            dist > AstronomicalUnits::new(1.01),
            "Earth at aphelion should be > 1.01 AU from Sun, got {dist:?}"
        );
    }

    #[test]
    fn earth_heliocentric_at_epoch_2026() {
        let jd = epoch_2026();
        let earth = Vsop87Ephemeris::earth_heliocentric(jd);
        let pos = earth.get_position();

        let dist = pos.distance();
        assert!(
            dist > AstronomicalUnits::new(0.98) && dist < AstronomicalUnits::new(1.02),
            "Earth heliocentric at 2026 should be ~1 AU"
        );
    }

    // -------------------------------------------------------------------------
    // earth_barycentric_velocity tests
    // -------------------------------------------------------------------------

    #[test]
    fn earth_barycentric_velocity_at_j2000() {
        let jd = j2000();
        let vel = Vsop87Ephemeris::earth_barycentric_velocity(jd);

        // Earth's orbital velocity is ~30 km/s = ~0.0172 AU/day
        // Check that velocity magnitude is reasonable
        let vx = vel.x().value();
        let vy = vel.y().value();
        let vz = vel.z().value();
        let speed = (vx * vx + vy * vy + vz * vz).sqrt();

        // Expected: ~0.0172 AU/day (Earth mean orbital velocity)
        assert!(
            speed > 0.015 && speed < 0.02,
            "Earth orbital velocity should be ~0.017 AU/day, got {speed}"
        );
    }

    #[test]
    fn earth_barycentric_velocity_at_epoch_2020() {
        let jd = epoch_2020();
        let vel = Vsop87Ephemeris::earth_barycentric_velocity(jd);

        let vx = vel.x().value();
        let vy = vel.y().value();
        let vz = vel.z().value();
        let speed = (vx * vx + vy * vy + vz * vz).sqrt();

        assert!(
            speed > 0.015 && speed < 0.02,
            "Earth orbital velocity at 2020 should be ~0.017 AU/day"
        );
    }

    #[test]
    fn earth_barycentric_velocity_finite() {
        let jd = epoch_2026();
        let vel = Vsop87Ephemeris::earth_barycentric_velocity(jd);

        assert!(vel.x().value().is_finite(), "vx should be finite");
        assert!(vel.y().value().is_finite(), "vy should be finite");
        assert!(vel.z().value().is_finite(), "vz should be finite");
    }

    // -------------------------------------------------------------------------
    // moon_geocentric tests
    // -------------------------------------------------------------------------

    #[test]
    fn moon_geocentric_at_j2000() {
        let jd = j2000();
        let moon = Vsop87Ephemeris::moon_geocentric(jd);

        // Moon distance from Earth is ~384,400 km on average
        let dist = moon.distance();
        assert!(
            dist > Kilometers::new(356_000.0) && dist < Kilometers::new(407_000.0),
            "Moon distance should be 356k-407k km, got {dist:?}"
        );
    }

    #[test]
    fn moon_geocentric_at_perigee() {
        // Moon perigee example: varies each orbit
        // Average perigee ~356,500 km
        let jd = epoch_2020();
        let moon = Vsop87Ephemeris::moon_geocentric(jd);

        let dist = moon.distance();
        // Should still be within lunar distance bounds
        assert!(
            dist > Kilometers::new(350_000.0) && dist < Kilometers::new(410_000.0),
            "Moon distance at 2020 should be reasonable"
        );
    }

    #[test]
    fn moon_geocentric_at_epoch_2026() {
        let jd = epoch_2026();
        let moon = Vsop87Ephemeris::moon_geocentric(jd);

        let dist = moon.distance();
        assert!(
            dist > Kilometers::new(350_000.0) && dist < Kilometers::new(410_000.0),
            "Moon distance at 2026 should be 350k-410k km"
        );
    }

    #[test]
    fn moon_geocentric_components_finite() {
        let jd = epoch_2026();
        let moon = Vsop87Ephemeris::moon_geocentric(jd);

        assert!(moon.x().value().is_finite(), "Moon x should be finite");
        assert!(moon.y().value().is_finite(), "Moon y should be finite");
        assert!(moon.z().value().is_finite(), "Moon z should be finite");
    }

    // -------------------------------------------------------------------------
    // Consistency tests
    // -------------------------------------------------------------------------

    #[test]
    fn earth_bary_minus_sun_bary_equals_earth_helio() {
        // earth_heliocentric = earth_barycentric - sun_barycentric
        let jd = epoch_2020();
        let earth_bary = Vsop87Ephemeris::earth_barycentric(jd);
        let sun_bary = Vsop87Ephemeris::sun_barycentric(jd);
        let earth_helio = Vsop87Ephemeris::earth_heliocentric(jd);

        // Calculate expected heliocentric position
        let expected_x = earth_bary.get_position().x() - sun_bary.get_position().x();
        let expected_y = earth_bary.get_position().y() - sun_bary.get_position().y();
        let expected_z = earth_bary.get_position().z() - sun_bary.get_position().z();

        // Compare with actual heliocentric position (should be very close)
        let tol = AstronomicalUnits::new(1e-6);
        assert!(
            (earth_helio.get_position().x() - expected_x).abs() < tol,
            "x mismatch in heliocentric consistency check"
        );
        assert!(
            (earth_helio.get_position().y() - expected_y).abs() < tol,
            "y mismatch in heliocentric consistency check"
        );
        assert!(
            (earth_helio.get_position().z() - expected_z).abs() < tol,
            "z mismatch in heliocentric consistency check"
        );
    }

    #[test]
    fn position_continuous_over_time() {
        // Check that positions change smoothly (no jumps)
        let jd1 = epoch_2020();
        let jd2 = jd_from_value(epoch_2020().value() + 1.0);

        let earth1 = Vsop87Ephemeris::earth_barycentric(jd1);
        let earth2 = Vsop87Ephemeris::earth_barycentric(jd2);

        // Position should not change by more than ~0.02 AU in one day
        let dx = (earth2.get_position().x() - earth1.get_position().x()).abs();
        let dy = (earth2.get_position().y() - earth1.get_position().y()).abs();
        let dz = (earth2.get_position().z() - earth1.get_position().z()).abs();

        let daily_move = AstronomicalUnits::new(0.03); // Conservative bound
        assert!(
            dx < daily_move && dy < daily_move && dz < daily_move,
            "Earth should move smoothly between consecutive days"
        );
    }

    #[test]
    fn vsop87_ephemeris_is_zero_sized() {
        // Verify that Vsop87Ephemeris is a ZST (zero-sized type)
        assert_eq!(
            std::mem::size_of::<Vsop87Ephemeris>(),
            0,
            "Vsop87Ephemeris should be zero-sized"
        );
    }

    #[test]
    fn vsop87_ephemeris_default() {
        // Test that Vsop87Ephemeris implements Default
        fn assert_default<T: Default>() {}
        assert_default::<Vsop87Ephemeris>();
    }

    #[test]
    fn vsop87_ephemeris_clone_copy() {
        // Test Clone and Copy traits
        fn assert_clone<T: Clone>() {}
        fn assert_copy<T: Copy>() {}
        assert_clone::<Vsop87Ephemeris>();
        assert_copy::<Vsop87Ephemeris>();

        let eph1 = Vsop87Ephemeris;
        let eph2 = eph1;
        let eph3 = eph1;
        // All are the same ZST
        assert_eq!(std::mem::size_of_val(&eph2), std::mem::size_of_val(&eph3));
    }

    #[test]
    fn vsop87_ephemeris_debug() {
        // Test Debug trait
        let eph = Vsop87Ephemeris;
        let debug_str = format!("{:?}", eph);
        assert!(debug_str.contains("Vsop87Ephemeris"));
    }
}

// ============================================================================
// DE440 Ephemeris Tests (Feature-gated)
// ============================================================================

#[cfg(feature = "de440")]
mod de440_tests {
    use super::*;
    use siderust::calculus::ephemeris::De440Ephemeris;

    #[test]
    fn de440_sun_barycentric_at_j2000() {
        let jd = j2000();
        let sun = De440Ephemeris::sun_barycentric(jd);
        let pos = sun.get_position();

        // Sun should be very close to SSB
        let dist = pos.distance();
        assert!(
            dist < AstronomicalUnits::new(0.02),
            "DE440: Sun should be close to SSB at J2000"
        );
    }

    #[test]
    fn de440_earth_barycentric_at_j2000() {
        let jd = j2000();
        let earth = De440Ephemeris::earth_barycentric(jd);
        let pos = earth.get_position();

        let dist = pos.distance();
        assert!(
            dist > AstronomicalUnits::new(0.98) && dist < AstronomicalUnits::new(1.02),
            "DE440: Earth should be ~1 AU from SSB"
        );
    }

    #[test]
    fn de440_earth_heliocentric_at_j2000() {
        let jd = j2000();
        let earth = De440Ephemeris::earth_heliocentric(jd);
        let pos = earth.get_position();

        let dist = pos.distance();
        assert!(
            dist > AstronomicalUnits::new(0.98) && dist < AstronomicalUnits::new(1.02),
            "DE440: Earth heliocentric should be ~1 AU"
        );
    }

    #[test]
    fn de440_earth_barycentric_velocity_at_j2000() {
        let jd = j2000();
        let vel = De440Ephemeris::earth_barycentric_velocity(jd);

        let vx = vel.x().value();
        let vy = vel.y().value();
        let vz = vel.z().value();
        let speed = (vx * vx + vy * vy + vz * vz).sqrt();

        assert!(
            speed > 0.015 && speed < 0.02,
            "DE440: Earth velocity should be ~0.017 AU/day"
        );
    }

    #[test]
    fn de440_moon_geocentric_at_j2000() {
        let jd = j2000();
        let moon = De440Ephemeris::moon_geocentric(jd);

        let dist = moon.distance();
        assert!(
            dist > Kilometers::new(356_000.0) && dist < Kilometers::new(407_000.0),
            "DE440: Moon distance should be 356k-407k km"
        );
    }

    #[test]
    fn de440_consistency_sun_versus_vsop87() {
        let jd = j2000();
        let sun_de440 = De440Ephemeris::sun_barycentric(jd);
        let sun_vsop = Vsop87Ephemeris::sun_barycentric(jd);

        // DE440 and VSOP87 should give similar results (within ~1000 km = ~7e-6 AU)
        let tol = AstronomicalUnits::new(1e-4); // 1e-4 AU ~ 15,000 km tolerance
        let dx = (sun_de440.get_position().x() - sun_vsop.get_position().x()).abs();
        let dy = (sun_de440.get_position().y() - sun_vsop.get_position().y()).abs();
        let dz = (sun_de440.get_position().z() - sun_vsop.get_position().z()).abs();

        assert!(
            dx < tol && dy < tol && dz < tol,
            "DE440 and VSOP87 Sun positions should be similar at J2000"
        );
    }

    #[test]
    fn de440_consistency_earth_versus_vsop87() {
        let jd = j2000();
        let earth_de440 = De440Ephemeris::earth_barycentric(jd);
        let earth_vsop = Vsop87Ephemeris::earth_barycentric(jd);

        let tol = AstronomicalUnits::new(1e-4);
        let dx = (earth_de440.get_position().x() - earth_vsop.get_position().x()).abs();
        let dy = (earth_de440.get_position().y() - earth_vsop.get_position().y()).abs();
        let dz = (earth_de440.get_position().z() - earth_vsop.get_position().z()).abs();

        assert!(
            dx < tol && dy < tol && dz < tol,
            "DE440 and VSOP87 Earth positions should be similar at J2000"
        );
    }

    #[test]
    fn de440_at_epoch_2020() {
        let jd = epoch_2020();

        let sun = De440Ephemeris::sun_barycentric(jd);
        let earth_bary = De440Ephemeris::earth_barycentric(jd);
        let earth_helio = De440Ephemeris::earth_heliocentric(jd);
        let vel = De440Ephemeris::earth_barycentric_velocity(jd);
        let moon = De440Ephemeris::moon_geocentric(jd);

        // Verify all are finite and reasonable
        assert!(sun.get_position().distance() < AstronomicalUnits::new(0.03));
        assert!(earth_bary.get_position().distance() > AstronomicalUnits::new(0.9));
        assert!(earth_helio.get_position().distance() > AstronomicalUnits::new(0.9));
        assert!(vel.x().value().is_finite());
        assert!(moon.distance() > Kilometers::new(350_000.0));
    }

    #[test]
    fn de440_at_epoch_2026() {
        let jd = epoch_2026();

        let sun = De440Ephemeris::sun_barycentric(jd);
        let earth_bary = De440Ephemeris::earth_barycentric(jd);
        let earth_helio = De440Ephemeris::earth_heliocentric(jd);
        let vel = De440Ephemeris::earth_barycentric_velocity(jd);
        let moon = De440Ephemeris::moon_geocentric(jd);

        assert!(sun.get_position().x().value().is_finite());
        assert!(earth_bary.get_position().x().value().is_finite());
        assert!(earth_helio.get_position().x().value().is_finite());
        assert!(vel.x().value().is_finite());
        assert!(moon.x().value().is_finite());
    }
}

// ============================================================================
// DE441 Ephemeris Tests (Feature-gated)
// ============================================================================

#[cfg(feature = "de441")]
mod de441_tests {
    use super::*;
    use siderust::calculus::ephemeris::De441Ephemeris;

    #[test]
    fn de441_sun_barycentric_at_j2000() {
        let jd = j2000();
        let sun = De441Ephemeris::sun_barycentric(jd);
        let pos = sun.get_position();

        let dist = pos.distance();
        assert!(
            dist < AstronomicalUnits::new(0.02),
            "DE441: Sun should be close to SSB at J2000"
        );
    }

    #[test]
    fn de441_earth_barycentric_at_j2000() {
        let jd = j2000();
        let earth = De441Ephemeris::earth_barycentric(jd);
        let pos = earth.get_position();

        let dist = pos.distance();
        assert!(
            dist > AstronomicalUnits::new(0.98) && dist < AstronomicalUnits::new(1.02),
            "DE441: Earth should be ~1 AU from SSB"
        );
    }

    #[test]
    fn de441_earth_heliocentric_at_j2000() {
        let jd = j2000();
        let earth = De441Ephemeris::earth_heliocentric(jd);
        let pos = earth.get_position();

        let dist = pos.distance();
        assert!(
            dist > AstronomicalUnits::new(0.98) && dist < AstronomicalUnits::new(1.02),
            "DE441: Earth heliocentric should be ~1 AU"
        );
    }

    #[test]
    fn de441_earth_barycentric_velocity_at_j2000() {
        let jd = j2000();
        let vel = De441Ephemeris::earth_barycentric_velocity(jd);

        let vx = vel.x().value();
        let vy = vel.y().value();
        let vz = vel.z().value();
        let speed = (vx * vx + vy * vy + vz * vz).sqrt();

        assert!(
            speed > 0.015 && speed < 0.02,
            "DE441: Earth velocity should be ~0.017 AU/day"
        );
    }

    #[test]
    fn de441_moon_geocentric_at_j2000() {
        let jd = j2000();
        let moon = De441Ephemeris::moon_geocentric(jd);

        let dist = moon.distance();
        assert!(
            dist > Kilometers::new(356_000.0) && dist < Kilometers::new(407_000.0),
            "DE441: Moon distance should be 356k-407k km"
        );
    }

    #[test]
    fn de441_consistency_earth_versus_vsop87() {
        let jd = j2000();
        let earth_de441 = De441Ephemeris::earth_barycentric(jd);
        let earth_vsop = Vsop87Ephemeris::earth_barycentric(jd);

        let tol = AstronomicalUnits::new(1e-4);
        let dx = (earth_de441.get_position().x() - earth_vsop.get_position().x()).abs();
        let dy = (earth_de441.get_position().y() - earth_vsop.get_position().y()).abs();
        let dz = (earth_de441.get_position().z() - earth_vsop.get_position().z()).abs();

        assert!(
            dx < tol && dy < tol && dz < tol,
            "DE441 and VSOP87 Earth positions should be similar at J2000"
        );
    }

    #[test]
    fn de441_at_epoch_2020() {
        let jd = epoch_2020();

        let sun = De441Ephemeris::sun_barycentric(jd);
        let earth_bary = De441Ephemeris::earth_barycentric(jd);
        let earth_helio = De441Ephemeris::earth_heliocentric(jd);
        let vel = De441Ephemeris::earth_barycentric_velocity(jd);
        let moon = De441Ephemeris::moon_geocentric(jd);

        assert!(sun.get_position().distance() < AstronomicalUnits::new(0.03));
        assert!(earth_bary.get_position().distance() > AstronomicalUnits::new(0.9));
        assert!(earth_helio.get_position().distance() > AstronomicalUnits::new(0.9));
        assert!(vel.x().value().is_finite());
        assert!(moon.distance() > Kilometers::new(350_000.0));
    }

    #[test]
    fn de441_at_epoch_2026() {
        let jd = epoch_2026();

        let sun = De441Ephemeris::sun_barycentric(jd);
        let earth_bary = De441Ephemeris::earth_barycentric(jd);
        let earth_helio = De441Ephemeris::earth_heliocentric(jd);
        let vel = De441Ephemeris::earth_barycentric_velocity(jd);
        let moon = De441Ephemeris::moon_geocentric(jd);

        assert!(sun.get_position().x().value().is_finite());
        assert!(earth_bary.get_position().x().value().is_finite());
        assert!(earth_helio.get_position().x().value().is_finite());
        assert!(vel.x().value().is_finite());
        assert!(moon.x().value().is_finite());
    }

    #[test]
    fn de441_long_time_span() {
        // DE441 supports a long time span (-13200 to +17191 years)
        // Test a historical date
        let jd_historical = jd_from_value(2415020.0); // ~1900
        let earth_hist = De441Ephemeris::earth_barycentric(jd_historical);
        assert!(earth_hist.get_position().distance() > AstronomicalUnits::new(0.9));

        // Test a future date
        let jd_future = jd_from_value(2488069.5); // ~2100
        let earth_future = De441Ephemeris::earth_barycentric(jd_future);
        assert!(earth_future.get_position().distance() > AstronomicalUnits::new(0.9));
    }
}

// ============================================================================
// Generic Ephemeris Trait Tests
// ============================================================================

/// Tests that can be run with any Ephemeris implementation
mod generic_ephemeris_tests {
    use super::*;

    fn test_ephemeris_basic_properties<E: Ephemeris>() {
        let jd = j2000();

        // All positions should be finite
        let sun = E::sun_barycentric(jd);
        assert!(sun.get_position().x().value().is_finite());

        let earth_bary = E::earth_barycentric(jd);
        assert!(earth_bary.get_position().x().value().is_finite());

        let earth_helio = E::earth_heliocentric(jd);
        assert!(earth_helio.get_position().x().value().is_finite());

        let vel = E::earth_barycentric_velocity(jd);
        assert!(vel.x().value().is_finite());

        let moon = E::moon_geocentric(jd);
        assert!(moon.x().value().is_finite());
    }

    #[test]
    fn vsop87_basic_properties() {
        test_ephemeris_basic_properties::<Vsop87Ephemeris>();
    }

    #[cfg(feature = "de440")]
    #[test]
    fn de440_basic_properties() {
        use siderust::calculus::ephemeris::De440Ephemeris;
        test_ephemeris_basic_properties::<De440Ephemeris>();
    }

    #[cfg(feature = "de441")]
    #[test]
    fn de441_basic_properties() {
        use siderust::calculus::ephemeris::De441Ephemeris;
        test_ephemeris_basic_properties::<De441Ephemeris>();
    }

    /// Test that velocity is consistent with position change over time
    fn test_velocity_consistency<E: Ephemeris>() {
        let jd1 = epoch_2020();
        let jd2 = jd_from_value(epoch_2020().value() + 1.0); // 1 day later

        let earth1 = E::earth_barycentric(jd1);
        let earth2 = E::earth_barycentric(jd2);
        let vel = E::earth_barycentric_velocity(jd1);

        // Position change over 1 day should be approximately velocity * 1 day
        let dx = earth2.get_position().x() - earth1.get_position().x();
        let dy = earth2.get_position().y() - earth1.get_position().y();
        let dz = earth2.get_position().z() - earth1.get_position().z();

        // Velocity is in AU/day, so expected change ~ velocity * 1 day
        // Allow 10% tolerance for interpolation/nonlinearity
        let tol = AstronomicalUnits::new(0.002); // ~ 10% of daily motion

        let vx_expected = AstronomicalUnits::new(vel.x().value());
        let vy_expected = AstronomicalUnits::new(vel.y().value());
        let vz_expected = AstronomicalUnits::new(vel.z().value());

        assert!((dx - vx_expected).abs() < tol, "dx should match vx * 1 day");
        assert!((dy - vy_expected).abs() < tol, "dy should match vy * 1 day");
        assert!((dz - vz_expected).abs() < tol, "dz should match vz * 1 day");
    }

    #[test]
    fn vsop87_velocity_consistency() {
        test_velocity_consistency::<Vsop87Ephemeris>();
    }

    #[cfg(feature = "de440")]
    #[test]
    fn de440_velocity_consistency() {
        use siderust::calculus::ephemeris::De440Ephemeris;
        test_velocity_consistency::<De440Ephemeris>();
    }

    #[cfg(feature = "de441")]
    #[test]
    fn de441_velocity_consistency() {
        use siderust::calculus::ephemeris::De441Ephemeris;
        test_velocity_consistency::<De441Ephemeris>();
    }
}

// ============================================================================
// Time Edge Cases
// ============================================================================

mod edge_cases {
    use super::*;

    #[test]
    fn very_early_epoch() {
        // Test early epoch (1800)
        let jd = jd_from_value(2378497.0); // ~1800-01-01
        let earth = Vsop87Ephemeris::earth_barycentric(jd);
        assert!(earth.get_position().distance() > AstronomicalUnits::new(0.9));
    }

    #[test]
    fn very_late_epoch() {
        // Test late epoch (2200)
        let jd = jd_from_value(2524594.0); // ~2200-01-01
        let earth = Vsop87Ephemeris::earth_barycentric(jd);
        assert!(earth.get_position().distance() > AstronomicalUnits::new(0.9));
    }

    #[test]
    fn half_day_offset() {
        // Test at noon vs midnight
        let jd_noon = jd_from_value(2451545.0); // J2000.0 at 12:00
        let jd_midnight = jd_from_value(2451545.5); // J2000.0 at 00:00 next day

        let earth_noon = Vsop87Ephemeris::earth_barycentric(jd_noon);
        let earth_midnight = Vsop87Ephemeris::earth_barycentric(jd_midnight);

        // Position should differ by ~half a day's motion
        let dx = (earth_midnight.get_position().x() - earth_noon.get_position().x()).abs();
        assert!(
            dx > AstronomicalUnits::new(0.0001) && dx < AstronomicalUnits::new(0.02),
            "Earth should move noticeably in 12 hours"
        );
    }

    #[test]
    fn fractional_julian_day() {
        // Test with precise fractional JD
        let jd = jd_from_value(2451545.123456789);
        let earth = Vsop87Ephemeris::earth_barycentric(jd);
        assert!(earth.get_position().x().value().is_finite());
    }
}
