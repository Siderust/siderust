use crate::coordinates::{
    cartesian::Direction, centers, frames, kinds::PositionKind
};

use super::CartesianCoord;

pub type Position<C, F>  = CartesianCoord<C, F, PositionKind>;

impl<C, F> Position<C, F>
where
    C: centers::ReferenceCenter,
    F: frames::ReferenceFrame,
{
    /// Calculates the Euclidean distance with respect to the ReferenceCenter.
    ///
    /// # Returns
    /// The distance from the ReferenceCenter in AU.
    pub fn distance_from_origin(&self) -> f64 {
        (self.x().powi(2) + self.y().powi(2) + self.z().powi(2)).sqrt()
    }

    pub fn normalize(&self) -> Self {
        Self::from_vec3(self.as_vec3().normalize())
    }

    pub fn direction(&self) -> super::Direction<C, F> {
        Direction::from_vec3(self.as_vec3().normalize())
    }

}

// === ICRS-based Cartesian coordinate types ===
pub type ICRS = Position<centers::Barycentric,  frames::ICRS>;
pub type HCRS = Position<centers::Heliocentric, frames::ICRS>;
pub type GCRS = Position<centers::Geocentric,   frames::ICRS>;
pub type TCRS = Position<centers::Topocentric,  frames::ICRS>;

// === Ecliptic frame ===
pub type EclipticBarycentricCartesianCoord  = Position<centers::Barycentric,  frames::Ecliptic>;
pub type EclipticHeliocentricCartesianCoord = Position<centers::Heliocentric, frames::Ecliptic>;
pub type EclipticGeocentricCartesianCoord   = Position<centers::Geocentric,   frames::Ecliptic>;
pub type EclipticTopocentricCartesianCoord  = Position<centers::Topocentric,  frames::Ecliptic>;

// === Equatorial frame ===
pub type EquatorialBarycentricCartesianCoord  = Position<centers::Barycentric,  frames::Equatorial>;
pub type EquatorialHeliocentricCartesianCoord = Position<centers::Heliocentric, frames::Equatorial>;
pub type EquatorialGeocentricCartesianCoord   = Position<centers::Geocentric,   frames::Equatorial>;
pub type EquatorialTopocentricCartesianCoord  = Position<centers::Topocentric,  frames::Equatorial>;

// === Horizontal and Earth-fixed frame ===
pub type HorizontalTopocentricCartesianCoord  = Position<centers::Topocentric,  frames::Horizontal>;

// === Geographic and Earth-fixed frames ===
pub type GeographicCartesianCoord = Position<centers::Geocentric, frames::ECEF>;
