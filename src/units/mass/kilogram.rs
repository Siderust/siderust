//! Kilograms (kg) type and conversions.
//!
//! Provides a strongly-typed representation of a mass in Kilograms (kg).

pub const KG: Kilograms = Kilograms::new(1.0);

/// A strongly-typed representation of a mass in Kilograms (kg).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Kilograms(f64);

impl Kilograms {
    /// Creates a new `Kilograms` from a value in kg.
    ///
    /// # Arguments
    /// * `value` - The mass in kilograms.
    pub const fn new(value: f64) -> Self {
        Self(value)
    }

    /// Returns the inner value in kilograms.
    pub const fn value(&self) -> f64 {
        self.0
    }
}

impl std::fmt::Display for Kilograms {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} kg", self.0)
    }
}

// Arithmetic operations macro
crate::units::arithmetic_ops::impl_arithmetic_ops!(Kilograms);
