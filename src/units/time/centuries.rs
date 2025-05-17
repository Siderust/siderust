//! Time unit: Centuries
//!
//! This module provides the `Centuries` struct, representing a span of time in centuries (100-year periods).
//! It supports basic arithmetic operations and formatting.


/// Represents a span of time in centuries (100-year periods).
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct Centuries(pub f64);

impl Centuries {
        pub const fn new(centuries: f64) -> Self {
        Centuries(centuries)
    }

        pub fn value(&self) -> f64 {
        self.0
    }
}

impl std::fmt::Display for Centuries {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} Centuries", self.0)
    }
}

impl std::ops::Div<Centuries> for Centuries {
    type Output = f64;
    fn div(self, rhs: Centuries) -> Self::Output {
        self.value() / rhs.value()
    }
}

crate::units::arithmetic_ops::impl_arithmetic_ops!(Centuries);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_centuries_new_and_value() {
        let c = Centuries::new(2.5);
        assert_eq!(c.value(), 2.5);
    }

    #[test]
    fn test_centuries_display() {
        let c = Centuries::new(1.0);
        assert_eq!(format!("{}", c), "1 Centuries");
    }

    #[test]
    fn test_centuries_div() {
        let c1 = Centuries::new(4.0);
        let c2 = Centuries::new(2.0);
        assert_eq!(c1 / c2, 2.0);
    }
}
