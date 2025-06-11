// src/macros.rs
use core::f64;
use crate::coordinates::{ cartesian, spherical, centers::ReferenceCenter, frames::ReferenceFrame };

#[doc(hidden)]
pub(crate) fn __assert_cartesian_eq<C, F>(
    a: &cartesian::Position<C, F>,
    b: &cartesian::Position<C, F>,
    epsilon: f64,
    msg: Option<String>,
)
where C: ReferenceCenter, F: ReferenceFrame
{
    let dx = (a.x() - b.x()).abs();
    let dy = (a.y() - b.y()).abs();
    let dz = (a.z() - b.z()).abs();
    if dx >= epsilon || dy >= epsilon || dz >= epsilon {
        if let Some(m) = msg {
            panic!("{}", m);
        } else {
            panic!("Cartesian coords differ: {} vs {} (ε = {})", a, b, epsilon);
        }
    }
}

#[doc(hidden)]
pub(crate) fn __assert_spherical_eq<C, F>(
    a: &spherical::Position<C, F>,
    b: &spherical::Position<C, F>,
    epsilon: f64,
    msg: Option<String>,
)
where C: ReferenceCenter, F: ReferenceFrame
{
    let d1 = a.distance.unwrap_or(f64::NAN);
    let d2 = b.distance.unwrap_or(f64::NAN);
    let dp = (a.polar.as_f64()   - b.polar.as_f64()).abs();
    let da = (a.azimuth.as_f64() - b.azimuth.as_f64()).abs();
    if (d1 - d2).abs() >= epsilon || dp >= epsilon || da >= epsilon {
        if let Some(m) = msg {
            panic!("{}", m);
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
