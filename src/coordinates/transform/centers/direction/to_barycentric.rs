use crate::astro::JulianDate;
use crate::coordinates::{
    centers::*,
    frames::MutableFrame,
    transform::TransformFrame,
    cartesian::direction::{Direction, Equatorial},
};
use crate::coordinates::transform::centers::TransformCenter;
use crate::astro::aberration::remove_aberration_from_direction;

// Heliocentric To Barycentric
impl<F: MutableFrame> TransformCenter<Direction<Barycentric, F>>
    for Direction<Heliocentric, F>
where
{
    fn to_center(&self, _jd: JulianDate) -> Direction<Barycentric, F> {
        Direction::from_vec3(self.as_vec3())
    }
}

// Geocentric To Barycentric
impl<F: MutableFrame> TransformCenter<Direction<Barycentric, F>>
    for Direction<Geocentric, F>
where
    Direction<Geocentric, F>: TransformFrame<Equatorial>, // ToEquatorial
    Equatorial: TransformFrame<Direction<Geocentric, F>>, // FromEquatorial
{
    #[inline]
    fn to_center(&self, jd: JulianDate) -> Direction<Barycentric, F> {
        // 1. Transform to Equatorial (already in Geocentric)
        let equatorial: Equatorial = self.to_frame();
        // 2. Remove aberration
        let deaberrated =  remove_aberration_from_direction(equatorial, jd);
        // 3. Recover target Frame
        let target_center: Direction::<Geocentric, F> = deaberrated.to_frame();
        // 4. Transform target Center
        Direction::<Barycentric, F>::from_vec3(target_center.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::aberration::apply_aberration_to_direction;
    use crate::coordinates::{cartesian::direction, centers};
    use crate::units::Degrees;

    const EPS: f64 = 1e-12;

    #[test]
    fn heliocentric_to_barycentric_identity() {
        let jd = JulianDate::J2000;
        let helio = direction::Equatorial::<centers::Heliocentric>::normalize(0.1, 0.2, 0.9);
        let bary: direction::Equatorial<centers::Barycentric> = helio.to_center(jd);

        assert!((helio.x() - bary.x()).abs().value() < EPS);
        assert!((helio.y() - bary.y()).abs().value() < EPS);
        assert!((helio.z() - bary.z()).abs().value() < EPS);
    }

    #[test]
    fn geocentric_apparent_to_barycentric_matches_mean() {
        let jd = JulianDate::J2000;
        // Start from a mean (aberration-free) direction defined via RA/Dec
        let mean_sph = crate::coordinates::spherical::direction::Equatorial::new(
            Degrees::new(120.0),
            Degrees::new(-30.0),
        );
        let mean = mean_sph.to_cartesian();
        // Add aberration to obtain an apparent geocentric direction
        let apparent = apply_aberration_to_direction(mean, jd);
        // Convert to barycentric, which should remove aberration
        let bary: direction::Equatorial<centers::Barycentric> = apparent.to_center(jd);

        assert!((bary.x() - mean.x()).abs().value() < 1e-8, "current {}, expected {}", bary.x(), mean.x());
        assert!((bary.y() - mean.y()).abs().value() < 1e-8, "current {}, expected {}", bary.x(), mean.x());
        assert!((bary.z() - mean.z()).abs().value() < 1e-8, "current {}, expected {}", bary.x(), mean.x());
    }
}
