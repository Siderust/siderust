use super::*;

pub enum Dimensionless {}
impl Dimension for Dimensionless {}

impl Unit for f64 {
    const RATIO: f64 = 1.0;
    type Dim = Dimensionless;
    const SYMBOL: &'static str = "";
}
impl std::fmt::Display for Quantity<f64> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {}", self.value(), "")
    }
}

impl From<Quantity<U>> for Quantity<f64>
where U: LengthUnit
{
    fn from(length: Quantity<U>) -> Self {
        Self::new(length.value())
    }
}
