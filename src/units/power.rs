//! # Power Units Module
//!
//! Commonly used power units in astronomy.
//!
//! ## Units
//! - **Watt (W)**: SI base unit of power.
//! - **Solar Luminosity (L☉)**: power radiated by the Sun.
//!
//! ## Example
//! ```rust
//! use siderust::units::*;
//!
//! // 2 kW
//! let p_w = Watts::new(2_000.0);
//!
//! // 3 L☉
//! let p_sol = SolarLuminosities::new(3.0);
//!
//! // Convenient conversion
//! let p_w_equiv = p_sol.to::<Watt>();
//! assert!((p_w_equiv.value() - 3.0 * 3.828e26).abs() < 1e15);
//! ```

use super::*;

/// Fundamental dimension – power.
pub enum Power {}
impl Dimension for Power {}

/// Marker trait for power units.
pub trait PowerUnit: Unit<Dim = Power> {}
impl<T: Unit<Dim = Power>> PowerUnit for T {}

define_unit!("W", Watt, Power, 1.0);
pub type W = Watt;
pub type Watts = Quantity<W>;
pub const WATT: Watts = Watts::new(1.0);

// 1 L☉ = 3.828 × 10²⁶ W (IAU 2015)
define_unit!("L☉", SolarLuminosity, Power, 3.828e26);
pub type SolarLuminosities = Quantity<SolarLuminosity>;
pub const L_SUN: SolarLuminosities = SolarLuminosities::new(1.0);
