//! # Velocity Units Module
//!
//! Composite velocity units built from a length over a time. The underlying
//! implementation leverages the generic [`Per<N, D>`] type which composes two
//! phantom unit parameters. This allows the multiplication and division traits
//! to be implemented once for all unit pairs in [`Quantity`](super::Quantity).

use super::*;

/// Dimension alias for velocities (`Length / Time`).
pub type Velocity = DivDim<Length, Time>;

/// Marker trait for any unit with velocity dimension.
pub trait VelocityUnit: Unit<Dim = Velocity> {}
impl<T: Unit<Dim = Velocity>> VelocityUnit for T {}

pub type MeterPerSecond = Per<Meter, Second>;
pub type MetersPerSecond = Quantity<MeterPerSecond>;

pub type KilometerPerSecond = Per<Kilometer, Second>;
pub type KilometersPerSecond = Quantity<KilometerPerSecond>;

pub type KilometerPerHour = Per<Kilometer, Hour>;
pub type KilometersPerHour = Quantity<KilometerPerHour>;

pub type AuPerDay = Per<Au, Day>;
pub type AusPerDay = Quantity<AuPerDay>;
