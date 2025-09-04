use super::*;

pub enum Dimensionless {}
impl Dimension for Dimensionless {}
pub type Unitless = f64;

impl Unit for Unitless {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "";
}
impl std::fmt::Display for Quantity<Unitless> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.value())
    }
}

impl<U: LengthUnit> From<Quantity<U>> for Quantity<Unitless> {
    fn from(length: Quantity<U>) -> Self {
        Self::new(length.value())
    }
}
