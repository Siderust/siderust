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
    };
}

new_center!(Barycentric);
new_center!(Heliocentric);
new_center!(Topocentric);

// Required for Transform specialization
#[derive(Debug, Copy, Clone)]
pub struct Geocentric;
impl ReferenceCenter for Geocentric {
    const IS_GEOCENTRIC: bool = true;

    fn center_name() -> &'static str {
        "Geocentric"
    }
}

pub trait NonGeocentric: ReferenceCenter {}
impl NonGeocentric for Heliocentric {}
impl NonGeocentric for Barycentric {}
impl NonGeocentric for Topocentric {}
