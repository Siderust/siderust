macro_rules! impl_arithmetic_ops {
    ($t:ty) => {

        impl std::ops::Add for $t {
            type Output = $t;
            fn add(self, rhs: $t) -> $t {
                <$t>::new(self.0 + rhs.0)
            }
        }

        impl std::ops::AddAssign for $t {
            fn add_assign(&mut self, other: $t) {
                self.0 += other.0;
            }
        }

        impl std::ops::Sub for $t {
            type Output = $t;
            fn sub(self, other: $t) -> $t {
                <$t>::new(self.0 - other.0)
            }
        }

        impl std::ops::SubAssign for $t {
            fn sub_assign(&mut self, other: $t) {
                self.0 -= other.0;
            }
        }

        impl std::ops::Mul<f64> for $t {
            type Output = $t;
            fn mul(self, scalar: f64) -> $t {
                <$t>::new(self.0 * scalar)
            }
        }


        impl  std::ops::Mul<$t> for f64 {
            type Output = $t;
        
            fn mul(self, rhs: $t) -> Self::Output {
                Self::Output::new(self * rhs.0)
            }
        }

        impl std::ops::Div<f64> for $t {
            type Output = $t;
            fn div(self, scalar: f64) -> $t {
                <$t>::new(self.0 / scalar)
            }
        }

        impl std::ops::Rem<f64> for $t {
            type Output = $t;
            fn rem(self, scalar: f64) -> $t {
                <$t>::new(self.0 % scalar)
            }
        }

        impl PartialEq<f64> for $t {
            fn eq(&self, other: &f64) -> bool {
                self.0 == *other
            }
        }

        impl std::ops::Neg for $t {
            type Output = $t;

            fn neg(self) -> $t {
                Self(-self.0)
            }
        }

    };
}


// Privately export it to be used by other modules within the crate
pub(crate) use impl_arithmetic_ops;
