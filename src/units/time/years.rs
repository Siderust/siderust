#[derive(Debug, Copy, Clone, PartialEq, PartialOrd)]
pub struct Years(pub f64);

impl Years {
    pub const fn new(years: f64) -> Self {
        Years(years)
    }

    pub fn value(&self) -> f64 {
        self.0
    }
}

impl std::fmt::Display for Years {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} Years", self.0)
    }
}

impl std::ops::Div<Years> for Years {
    type Output = f64;

    fn div(self, rhs: Years) -> Self::Output {
        self.value() / rhs.value()
    }
}