//! # Frequency Units Module
//!
//! Frequency units are represented as a ratio between an angular measure and a
//! time unit using the generic [`Per<N, D>`] type. This enables automatic
//! multiplication and division without bespoke implementations for every pair.

use super::*;

/// Dimension alias for angular frequency (`Angular / Time`).
pub type Frequency = DivDim<Angular, Time>;

/// Marker trait for any unit with frequency dimension.
pub trait FrequencyUnit: Unit<Dim = Frequency> {}
impl<T: Unit<Dim = Frequency>> FrequencyUnit for T {}

pub type RadianPerDay = Per<Radian, Day>;
pub type RadiansPerDay = Quantity<RadianPerDay>;

pub type DegreePerDay = Per<Degree, Day>;
pub type DegreesPerDay = Quantity<DegreePerDay>;

pub type DegreePerYear = Per<Degree, Year>;
pub type DegreesPerYear = Quantity<DegreePerYear>;

pub type MilliArcsecondPerDay = Per<MilliArcsecond, Day>;
pub type MilliArcsecondsPerDay = Quantity<MilliArcsecondPerDay>;
