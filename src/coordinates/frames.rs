//! # Reference Frames Module
//!
//! This module defines the concept of a *reference frame* for astronomical and geodetic coordinate systems.
//! A reference frame specifies the orientation of the axes used to describe positions in space.
//!
//! ## Overview
//!
//! The [`ReferenceFrame`] trait provides a common interface for all reference frame types. Each frame is
//! represented as a zero-sized struct and implements the trait to provide its canonical name.
//!
//! The `new_frame!` macro is used to conveniently declare new reference frame types, ensuring consistency
//! and reducing boilerplate.
//!
//! ## Predefined Frames
//!
//! The following reference frames are provided out of the box:
//!
//! - `ICRS`: International Celestial Reference System (quasi-inertial, used for most modern astronomy).
//! - `Horizontal`: Local horizon system (altitude-azimuth).
//! - `Equatorial`: Equatorial coordinate system (right ascension and declination).
//! - `Ecliptic`: Ecliptic coordinate system (based on the plane of Earth's orbit).
//! - `ITRF`: International Terrestrial Reference Frame (Earth-fixed).
//! - `ECEF`: Earth-Centered, Earth-Fixed (geocentric, rotating with the Earth).
//!
//! ## Extending
//!
//! To define a new reference frame, use the `new_frame!` macro.
//! This creates a new zero-sized type `Galactic` that implements [`ReferenceFrame`].
//!
//! ## Example
//!
//! ```rust
//! use siderust::coordinates::frames::{ReferenceFrame, ICRS};
//!
//! let name = ICRS::frame_name();
//! assert_eq!(name, "ICRS");
//! ```
/// A trait for defining a reference frame (orientation).
pub trait ReferenceFrame {
    fn frame_name() -> &'static str;
}

macro_rules! new_frame {
    ($name:ident) => {
        #[derive(Debug, Copy, Clone)]
        pub struct $name;

        impl ReferenceFrame for $name {
            fn frame_name() -> &'static str {
                stringify!($name)
            }
        }
    };
}

new_frame!(ICRS);
new_frame!(Horizontal);
new_frame!(Equatorial);
new_frame!(Ecliptic);
new_frame!(ITRF);
new_frame!(ECEF);


pub trait MutableFrame: ReferenceFrame {}
impl MutableFrame for ICRS {}
impl MutableFrame for Ecliptic {}
impl MutableFrame for Equatorial {}
