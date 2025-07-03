// src/macros.rs
use core::f64;
use crate::coordinates::{ cartesian, spherical, centers::ReferenceCenter, frames::ReferenceFrame };
use crate::units::Distance;

#[doc(hidden)]
pub(crate) fn __assert_cartesian_eq<C, F, U>(
    a: &cartesian::Position<C, F, U>,
    b: &cartesian::Position<C, F, U>,
    epsilon: f64,
    msg: Option<String>,
)
where C: ReferenceCenter, F: ReferenceFrame, U: Distance
{
    let dx = (a.x() - b.x()).abs();
    let dy = (a.y() - b.y()).abs();
    let dz = (a.z() - b.z()).abs();
    if dx >= epsilon.into() || dy >= epsilon.into() || dz >= epsilon.into() {
        if let Some(m) = msg {
            panic!("{}. Cartesian coords differ: {} vs {} (ε = {})", m, a, b, epsilon);
        } else {
            panic!("Cartesian coords differ: {} vs {} (ε = {})", a, b, epsilon);
        }
    }
}

#[doc(hidden)]
pub(crate) fn __assert_spherical_eq<C, F,  U>(
    a: &spherical::Position<C, F, U>,
    b: &spherical::Position<C, F, U>,
    epsilon: f64,
    msg: Option<String>,
)
where C: ReferenceCenter, F: ReferenceFrame, U: Distance
{
    let d1 = a.distance.unwrap_or(U::NAN);
    let d2 = b.distance.unwrap_or(U::NAN);
    let dp = (a.polar.as_f64()   - b.polar.as_f64()).abs();
    let da = (a.azimuth.as_f64() - b.azimuth.as_f64()).abs();
    if (d1 - d2).abs() >= epsilon.into() || dp >= epsilon || da >= epsilon {
        if let Some(m) = msg {
            panic!("{}. Spherical coords differ: {} vs {} (ε = {})", m, a, b, epsilon);
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
        // <-- note the extra `macros` path here:
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
