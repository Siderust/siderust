// src/macros.rs
use crate::coordinates::{cartesian, centers::ReferenceCenter, frames::ReferenceFrame, spherical};
use crate::units::{LengthUnit, Quantity};
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
    let dp = (a.polar.value() - b.polar.value()).abs();
    let da = (a.azimuth.value() - b.azimuth.value()).abs();
    if (d1 - d2).abs() >= Quantity::<U>::new(epsilon) || dp >= epsilon || da >= epsilon {
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
