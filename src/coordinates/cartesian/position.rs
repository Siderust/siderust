use crate::coordinates::{
    cartesian::Direction, centers, frames, kinds::{PositionKind, VelocityKind}
};
use super::Vector;
use crate::units::Unit;

pub type Position<C, F, U = f64>  = Vector<C, F, U, PositionKind>;
pub type Velocity<C, F, U = f64>  = Vector<C, F, U, VelocityKind>;

impl<C, F, U> Position<C, F, U>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
    U: Unit,
{
    /// Calculates the Euclidean distance with respect to the ReferenceCenter.
    ///
    /// # Returns
    /// The distance from the ReferenceCenter in units of U.
    pub fn distance(&self) -> U {
        (self.x().powi(2) + self.y().powi(2) + self.z().powi(2)).sqrt()
    }

    pub fn normalize(&self) -> Self {
        Self::from_vec3(
            self.as_vec3().normalize()
        )
    }

    pub fn direction(&self) -> super::Direction<C, F, U> {
        Direction::from_vec3(
            self.as_vec3().normalize()
        )
    }
}

// === ICRS-based Cartesian coordinate types ===
pub type ICRS<U = f64> = Position<centers::Barycentric,  frames::ICRS, U>;
pub type HCRS<U = f64> = Position<centers::Heliocentric, frames::ICRS, U>;
pub type GCRS<U = f64> = Position<centers::Geocentric,   frames::ICRS, U>;
pub type TCRS<U = f64> = Position<centers::Topocentric,  frames::ICRS, U>;

// === Ecliptic frame ===
pub type EclipticBarycentricCartesianCoord<U = f64>  = Position<centers::Barycentric,  frames::Ecliptic, U>;
pub type EclipticHeliocentricCartesianCoord<U = f64> = Position<centers::Heliocentric, frames::Ecliptic, U>;
pub type EclipticGeocentricCartesianCoord<U = f64>   = Position<centers::Geocentric,   frames::Ecliptic, U>;
pub type EclipticTopocentricCartesianCoord<U = f64>  = Position<centers::Topocentric,  frames::Ecliptic, U>;

// === Equatorial frame ===
pub type EquatorialBarycentricCartesianCoord<U = f64>  = Position<centers::Barycentric,  frames::Equatorial, U>;
pub type EquatorialHeliocentricCartesianCoord<U = f64> = Position<centers::Heliocentric, frames::Equatorial, U>;
pub type EquatorialGeocentricCartesianCoord<U = f64>   = Position<centers::Geocentric,   frames::Equatorial, U>;
pub type EquatorialTopocentricCartesianCoord<U = f64>  = Position<centers::Topocentric,  frames::Equatorial, U>;

// === Horizontal and Earth-fixed frame ===
pub type HorizontalTopocentricCartesianCoord<U = f64>  = Position<centers::Topocentric,  frames::Horizontal, U>;

// === Geographic and Earth-fixed frames ===
pub type GeographicCartesianCoord<U = f64> = Position<centers::Geocentric, frames::ECEF, U>;
