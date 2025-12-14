//! # Reference Centers Module
//!
//! This module defines the concept of a *reference center* (origin) for astronomical and geodetic coordinate systems.
//! A reference center specifies the origin point from which positions are measured.
//!
//! ## Overview
//!
//! The [`ReferenceCenter`] trait provides a common interface for all reference center types. Each center is
//! represented as a zero-sized struct and implements the trait to provide its canonical name.
//!
//! The `new_center!` macro is used to conveniently declare new reference center types, ensuring consistency
//! and reducing boilerplate.
//!
//! ## Predefined Centers
//!
//! The following reference centers are provided out of the box:
//!
//! - `Barycentric`: Center of mass of the solar system.
//! - `Heliocentric`: Center of the Sun.
//! - `Geocentric`: Center of the Earth.
//! - `Topocentric`: Observer's location on the surface of the Earth.
//!
//! ## Extending
//!
//! To define a new reference center, use the `new_center!` macro.
//!
//! This creates a new zero-sized type `Lunarcentric` that implements [`ReferenceCenter`].
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::centers::{ReferenceCenter, Geocentric};
//!
//! let name = Geocentric::center_name();
//! assert_eq!(name, "Geocentric");
//! ```

/// A trait for defining a reference center (coordinate origin).
pub trait ReferenceCenter {
    const IS_GEOCENTRIC: bool;
    fn center_name() -> &'static str;
}

// Required for Transform specialization
pub trait NonGeocentric: ReferenceCenter {}

macro_rules! new_center {
    ($name:ident) => {
        #[derive(Debug, Copy, Clone)]
        pub struct $name;

        impl ReferenceCenter for $name {
            const IS_GEOCENTRIC: bool = false;

            fn center_name() -> &'static str {
                stringify!($name)
            }
        }
        impl NonGeocentric for $name {}
    };
}

new_center!(Barycentric);
new_center!(Heliocentric);
new_center!(Topocentric);

#[derive(Debug, Copy, Clone)]
pub struct Geocentric;
impl ReferenceCenter for Geocentric {
    const IS_GEOCENTRIC: bool = true;

    fn center_name() -> &'static str {
        "Geocentric"
    }
}

impl ReferenceCenter for () {
    const IS_GEOCENTRIC: bool = false;
    fn center_name() -> &'static str {
        ""
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to verify that a center implements `NonGeocentric`
    fn assert_non_geocentric<C: NonGeocentric>() {}

    #[test]
    fn center_names_are_correct() {
        assert_eq!(Barycentric::center_name(), "Barycentric");
        assert_eq!(Heliocentric::center_name(), "Heliocentric");
        assert_eq!(Topocentric::center_name(), "Topocentric");
        assert_eq!(Geocentric::center_name(), "Geocentric");
        assert_eq!(<() as ReferenceCenter>::center_name(), "");
    }

    #[test]
    fn geocentric_flags() {
        // These are compile-time constants, so we just verify they compile correctly
        // rather than asserting their values at runtime
        const _: () = {
            let _ = !<Barycentric as ReferenceCenter>::IS_GEOCENTRIC;
            let _ = !<Heliocentric as ReferenceCenter>::IS_GEOCENTRIC;
            let _ = !<Topocentric as ReferenceCenter>::IS_GEOCENTRIC;
            let _ = <Geocentric as ReferenceCenter>::IS_GEOCENTRIC;
            let _ = !<() as ReferenceCenter>::IS_GEOCENTRIC;
        };
    }

    #[test]
    fn non_geocentric_trait_implemented() {
        // These calls will fail to compile if the trait is not implemented
        assert_non_geocentric::<Barycentric>();
        assert_non_geocentric::<Heliocentric>();
        assert_non_geocentric::<Topocentric>();
    }
}
