#[allow(unused_macros)]
macro_rules! assert_cartesian_eq {
    // ── standard form: a, b, ε ──────────────────────────────────────────────
    ($a:expr, $b:expr, $epsilon:expr $(,)?) => {{
        // compile-time type check (no run-time cost)
        fn _assert_same_type<T>(_: &T, _: &T) {}
        _assert_same_type(&$a, &$b);

        assert!(
            ($a.x() - $b.x()).abs() < $epsilon &&
            ($a.y() - $b.y()).abs() < $epsilon &&
            ($a.z() - $b.z()).abs() < $epsilon,
            "Cartesian coords differ: {} vs {} (ε = {})",
            $a, $b, $epsilon
        );
    }};

    // ── with user-supplied message: a, b, ε, "...", args… ───────────────────
    ($a:expr, $b:expr, $epsilon:expr, $($msg:tt)+) => {{
        fn _assert_same_type<T>(_: &T, _: &T) {}
        _assert_same_type(&$a, &$b);

        assert!(
            ($a.x() - $b.x()).abs() < $epsilon &&
            ($a.y() - $b.y()).abs() < $epsilon &&
            ($a.z() - $b.z()).abs() < $epsilon,
            $($msg)+               // forward the custom message & its args
        );
    }};
}

#[allow(unused_macros)]
macro_rules! assert_spherical_eq {
    ($a:expr, $b:expr, $epsilon:expr $(,)?) => {{
        fn _assert_same_type<T>(_: &T, _: &T) {}
        _assert_same_type(&$a, &$b);

        assert!(
            ($a.distance  - $b.distance).abs()  < $epsilon &&
            ($a.polar.as_f64()   - $b.polar.as_f64()).abs()   < $epsilon &&
            ($a.azimuth.as_f64() - $b.azimuth.as_f64()).abs() < $epsilon,
            "Spherical coords differ: {} vs {} (ε = {})",
            $a, $b, $epsilon
        );
    }};
    ($a:expr, $b:expr, $epsilon:expr, $($msg:tt)+) => {{
        fn _assert_same_type<T>(_: &T, _: &T) {}
        _assert_same_type(&$a, &$b);

        assert!(
            ($a.distance  - $b.distance).abs()  < $epsilon &&
            ($a.polar.as_f64()   - $b.polar.as_f64()).abs()   < $epsilon &&
            ($a.azimuth.as_f64() - $b.azimuth.as_f64()).abs() < $epsilon,
            $($msg)+
        );
    }};
}

#[allow(unused_imports)]
pub(crate) use assert_cartesian_eq;
#[allow(unused_imports)]
pub(crate) use assert_spherical_eq;
