// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! # Test-support macros for typed coordinate equality
//!
//! ## Scientific scope
//!
//! Astronomical coordinates are 3-tuples of floating-point quantities
//! that are computed by long chains of trigonometric and matrix
//! operations (precession, nutation, frame rotations, light-time
//! corrections, …). At every link in those chains, round-off introduces
//! sub-ULP errors that accumulate in the last few bits of the mantissa
//! but are scientifically negligible. Bitwise `assert_eq!` is therefore
//! the wrong test primitive for coordinate values: the meaningful
//! question is whether two `Position`s agree to within a tolerance
//! comparable to the precision of the reference value, not whether their
//! `f64` bit patterns are identical.
//!
//! These macros provide that element-wise tolerance check at the typed
//! coordinate boundary, so that tests express intent ("agree to
//! 1e-6 AU") rather than re-implementing the comparison by hand.
//!
//! ## Technical scope
//!
//! This module provides:
//!
//! - [`assert_cartesian_eq!`](crate::assert_cartesian_eq) — assert that
//!   two [`cartesian::Position`](crate::coordinates::cartesian::Position)
//!   values agree component-wise within an `epsilon` tolerance, panicking
//!   with a coordinate-aware message on failure.
//! - [`assert_spherical_eq!`](crate::assert_spherical_eq) — same for
//!   [`spherical::Position`](crate::coordinates::spherical::Position),
//!   comparing distance in the typed `LengthUnit`, and polar/azimuth in
//!   degrees.
//! - `__assert_cartesian_eq` / `__assert_spherical_eq` — the
//!   `pub(crate)` implementations the macros expand to.
//!
//! `epsilon` is intentionally a raw `f64` (no typed tolerance newtype):
//! these are testing helpers, not part of the science API.
//!
//! ## References
//!
//! None — this is test infrastructure with no domain-specific
//! algorithm.

use crate::coordinates::{cartesian, centers::ReferenceCenter, frames::ReferenceFrame, spherical};
use crate::qtty::{Degrees, LengthUnit, Quantity};
use core::f64;

#[doc(hidden)]
pub(crate) fn __assert_cartesian_eq<C, F, U>(
    a: &cartesian::Position<C, F, U>,
    b: &cartesian::Position<C, F, U>,
    epsilon: f64,
    msg: Option<String>,
) where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
    Quantity<U>: std::cmp::PartialOrd + std::fmt::Display,
{
    let dx = (a.x() - b.x()).abs();
    let dy = (a.y() - b.y()).abs();
    let dz = (a.z() - b.z()).abs();
    if dx >= Quantity::<U>::new(epsilon)
        || dy >= Quantity::<U>::new(epsilon)
        || dz >= Quantity::<U>::new(epsilon)
    {
        if let Some(m) = msg {
            panic!(
                "{}. Cartesian coords differ: {} vs {} (ε = {})",
                m, a, b, epsilon
            );
        } else {
            panic!("Cartesian coords differ: {} vs {} (ε = {})", a, b, epsilon);
        }
    }
}

#[doc(hidden)]
pub(crate) fn __assert_spherical_eq<C, F, U>(
    a: &spherical::Position<C, F, U>,
    b: &spherical::Position<C, F, U>,
    epsilon: f64,
    msg: Option<String>,
) where
    C: ReferenceCenter,
    F: ReferenceFrame,
    U: LengthUnit,
    Quantity<U>: std::cmp::PartialOrd + std::fmt::Display,
{
    let d1 = a.distance;
    let d2 = b.distance;
    let dp = (a.polar - b.polar).abs();
    let da = (a.azimuth - b.azimuth).abs();
    if (d1 - d2).abs() >= Quantity::<U>::new(epsilon)
        || dp >= Degrees::new(epsilon)
        || da >= Degrees::new(epsilon)
    {
        if let Some(m) = msg {
            panic!(
                "{}. Spherical coords differ: {} vs {} (ε = {})",
                m, a, b, epsilon
            );
        } else {
            panic!("Spherical coords differ: {} vs {} (ε = {})", a, b, epsilon);
        }
    }
}

/// Assert that two [`cartesian::Position`](crate::coordinates::cartesian::Position) values agree
/// component-wise within `epsilon` (same type and units). Panics with a
/// coordinate-aware message on failure.
///
/// # Examples
/// ```rust,ignore
/// use siderust::coordinates::cartesian;
/// use siderust::qtty::{AstronomicalUnit, AU};
/// use siderust::assert_cartesian_eq;
///
/// let a = cartesian::position::ICRS::<AstronomicalUnit>::new(1.0 * AU, 0.0 * AU, 0.0 * AU);
/// let b = cartesian::position::ICRS::<AstronomicalUnit>::new(1.0 * AU, 0.0 * AU, 0.0 * AU);
/// assert_cartesian_eq!(a, b, 1e-12);
/// ```
#[macro_export]
macro_rules! assert_cartesian_eq {
    ($a:expr, $b:expr, $eps:expr $(,)?) => {{
        fn _check<T>(_: &T, _: &T) {}
        _check(&$a, &$b);
        $crate::macros::__assert_cartesian_eq(&$a, &$b, $eps, None);
    }};
    ($a:expr, $b:expr, $eps:expr, $($msg:tt)+) => {{
        fn _check<T>(_: &T, _: &T) {}
        _check(&$a, &$b);
        $crate::macros::__assert_cartesian_eq(
            &$a,
            &$b,
            $eps,
            Some(format!($($msg)+))
        );
    }};
}

/// Assert that two [`spherical::Position`](crate::coordinates::spherical::Position) values agree
/// within `epsilon` (distance in the typed length unit, angles in degrees).
/// Panics with a coordinate-aware message on failure.
///
/// # Examples
/// ```rust,ignore
/// use siderust::coordinates::spherical;
/// use siderust::qtty::{AstronomicalUnit, Degrees};
/// use siderust::assert_spherical_eq;
///
/// let a = spherical::position::EquatorialMeanJ2000::<AstronomicalUnit>::new(
///     Degrees::new(0.0), Degrees::new(0.0), 1.0,
/// );
/// let b = spherical::position::EquatorialMeanJ2000::<AstronomicalUnit>::new(
///     Degrees::new(0.0), Degrees::new(0.0), 1.0,
/// );
/// assert_spherical_eq!(a, b, 1e-12);
/// ```
#[macro_export]
macro_rules! assert_spherical_eq {
    ($a:expr, $b:expr, $eps:expr $(,)?) => {{
        fn _check<T>(_: &T, _: &T) {}
        _check(&$a, &$b);
        $crate::macros::__assert_spherical_eq(&$a, &$b, $eps, None);
    }};
    ($a:expr, $b:expr, $eps:expr, $($msg:tt)+) => {{
        fn _check<T>(_: &T, _: &T) {}
        _check(&$a, &$b);
        $crate::macros::__assert_spherical_eq(
            &$a,
            &$b,
            $eps,
            Some(format!($($msg)+))
        );
    }};
}

#[allow(unused_imports)]
pub(crate) use assert_cartesian_eq;
#[allow(unused_imports)]
pub(crate) use assert_spherical_eq;

#[cfg(test)]
mod tests {
    use crate::coordinates::{cartesian, spherical};
    use crate::qtty::{AstronomicalUnit, Degrees, AU};

    #[test]
    #[should_panic(expected = "Cartesian coords differ")]
    fn cartesian_macro_panics_on_mismatch() {
        let a = cartesian::position::ICRS::<AstronomicalUnit>::new(0.0 * AU, 0.0 * AU, 0.0 * AU);
        let b = cartesian::position::ICRS::<AstronomicalUnit>::new(1.0 * AU, 0.0 * AU, 0.0 * AU);
        assert_cartesian_eq!(a, b, 1e-6);
    }

    #[test]
    #[should_panic(expected = "custom cart message")]
    fn cartesian_macro_reports_custom_message() {
        let a = cartesian::position::ICRS::<AstronomicalUnit>::new(0.0 * AU, 0.0 * AU, 0.0 * AU);
        let b = cartesian::position::ICRS::<AstronomicalUnit>::new(0.0 * AU, 1.0 * AU, 0.0 * AU);
        assert_cartesian_eq!(a, b, 1e-8, "custom cart message");
    }

    #[test]
    #[should_panic(expected = "Spherical coords differ")]
    fn spherical_macro_panics_on_mismatch() {
        let a = spherical::position::EquatorialMeanJ2000::<AstronomicalUnit>::new(
            Degrees::new(0.0),
            Degrees::new(0.0),
            1.0,
        );
        let b = spherical::position::EquatorialMeanJ2000::<AstronomicalUnit>::new(
            Degrees::new(10.0),
            Degrees::new(0.0),
            1.0,
        );
        assert_spherical_eq!(a, b, 1e-6);
    }

    #[test]
    #[should_panic(expected = "custom spherical message")]
    fn spherical_macro_reports_custom_message() {
        let a = spherical::position::EquatorialMeanJ2000::<AstronomicalUnit>::new(
            Degrees::new(0.0),
            Degrees::new(0.0),
            1.0,
        );
        let b = spherical::position::EquatorialMeanJ2000::<AstronomicalUnit>::new(
            Degrees::new(0.0),
            Degrees::new(20.0),
            1.0,
        );
        assert_spherical_eq!(a, b, 1e-6, "custom spherical message");
    }
}
