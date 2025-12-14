use crate::astro::aberration::remove_aberration_from_direction;
use crate::astro::JulianDate;
use crate::coordinates::spherical::direction::DirectionUnit;
use crate::coordinates::transform::centers::TransformCenter;
use crate::coordinates::{
    cartesian::direction::Equatorial, cartesian::Vector, centers::*, frames::MutableFrame,
    transform::TransformFrame,
};

// Heliocentric To Barycentric (Direction only - uses DirectionUnit)
impl<F: MutableFrame> TransformCenter<Vector<Barycentric, F, DirectionUnit>>
    for Vector<Heliocentric, F, DirectionUnit>
{
    fn to_center(&self, _jd: JulianDate) -> Vector<Barycentric, F, DirectionUnit> {
        Vector::from_vec3(self.as_vec3())
    }
}

// Geocentric To Barycentric (Direction only - uses DirectionUnit)
impl<F: MutableFrame> TransformCenter<Vector<Barycentric, F, DirectionUnit>>
    for Vector<Geocentric, F, DirectionUnit>
where
    Vector<Geocentric, F, DirectionUnit>: TransformFrame<Equatorial>, // ToEquatorial
    Equatorial: TransformFrame<Vector<Geocentric, F, DirectionUnit>>, // FromEquatorial
{
    #[inline]
    fn to_center(&self, jd: JulianDate) -> Vector<Barycentric, F, DirectionUnit> {
        // 1. Transform to Equatorial (already in Geocentric)
        let equatorial: Equatorial = self.to_frame();
        // 2. Remove aberration
        let deaberrated = remove_aberration_from_direction(equatorial, jd);
        // 3. Recover target Frame
        let target_center: Vector<Geocentric, F, DirectionUnit> = deaberrated.to_frame();
        // 4. Transform target Center
        Vector::<Barycentric, F, DirectionUnit>::from_vec3(target_center.as_vec3())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::astro::aberration::apply_aberration_to_direction;
    use crate::coordinates::{cartesian::direction, centers};
    use qtty::Degrees;

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

        assert!((bary.x() - mean.x()).abs().value() < 1e-8);
        assert!((bary.y() - mean.y()).abs().value() < 1e-8);
        assert!((bary.z() - mean.z()).abs().value() < 1e-8);
    }
}
