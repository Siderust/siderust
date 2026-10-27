// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Body-fixed Cartesian conversion and local ENU basis.

use super::body::LatitudeConvention;
use super::error::SurfaceError;
use super::reference_surface::ReferenceSurface;
use crate::qtty::{Degrees, Meter, Meters, Quantity, Radian};
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar, Transcendental};

/// Orthonormal east–north–up basis in a body-fixed Cartesian frame.
///
/// Components are dimensionless direction cosines. The triad is right-handed:
/// `east × north = up`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EnuBasis {
    /// Unit east vector `(ex, ey, ez)`.
    pub east: [f64; 3],
    /// Unit north vector `(nx, ny, nz)`.
    pub north: [f64; 3],
    /// Unit up / local-vertical vector `(ux, uy, uz)`.
    pub up: [f64; 3],
}

impl EnuBasis {
    /// Validate orthonormality and right-handedness within `tol`.
    pub fn assert_orthonormal(&self, tol: f64) -> Result<(), SurfaceError> {
        let dots = [
            ("e·e", dot(self.east, self.east) - 1.0),
            ("n·n", dot(self.north, self.north) - 1.0),
            ("u·u", dot(self.up, self.up) - 1.0),
            ("e·n", dot(self.east, self.north)),
            ("e·u", dot(self.east, self.up)),
            ("n·u", dot(self.north, self.up)),
        ];
        for (label, v) in dots {
            if !v.is_finite() || v.abs() > tol {
                return Err(SurfaceError::DegenerateGeometry { reason: label });
            }
        }
        let cross = cross(self.east, self.north);
        let handed = dot(cross, self.up);
        if (handed - 1.0).abs() > tol {
            return Err(SurfaceError::DegenerateGeometry {
                reason: "ENU not right-handed",
            });
        }
        Ok(())
    }
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn norm(v: [f64; 3]) -> f64 {
    dot(v, v).sqrt()
}

#[inline]
fn normalize(v: [f64; 3]) -> Result<[f64; 3], SurfaceError> {
    let n = norm(v);
    if !n.is_finite() || n == 0.0 {
        return Err(SurfaceError::DegenerateGeometry {
            reason: "zero-length basis vector",
        });
    }
    Ok([v[0] / n, v[1] / n, v[2] / n])
}

/// Convert surface coordinates to body-fixed Cartesian metres.
///
/// Longitude is positive **east**. Latitude follows `convention`.
/// Height is measured along the reference-surface normal for both latitude
/// conventions (radially on a sphere). Planetocentric latitude selects the
/// ellipsoid intersection radially; it does not change the local vertical.
pub fn body_fixed_position(
    surface: &ReferenceSurface,
    lon: Degrees,
    lat: Degrees,
    height: Meters,
    convention: LatitudeConvention,
) -> Result<[Meters; 3], SurfaceError> {
    validate_angles(lon, lat, height)?;
    surface.validate()?;

    let lon_r = lon.to::<Radian>().value();
    let lat_r = lat.to::<Radian>().value();
    let (sin_lon, cos_lon) = lon_r.sin_cos();
    let (sin_lat, cos_lat) = lat_r.sin_cos();
    let h = height.value();
    // A subsurface observer is valid until the inward normal reaches the
    // ellipsoid's focal/evolute singularity or the centre-facing tangent plane.
    // This geometric domain is independent of the chosen latitude convention.
    let radial_limit = match *surface {
        ReferenceSurface::Sphere { radius } => radius.value(),
        ReferenceSurface::Ellipsoid { a, f } => {
            let up = local_enu_basis(surface, lon, lat, convention)?.up;
            let sin2 = up[2] * up[2];
            let e2 = 2.0 * f - f * f;
            let one_minus = 1.0 - e2 * sin2;
            let a = a.value();
            let meridional = a * (1.0 - e2) / (one_minus * one_minus.sqrt());
            let tangent_distance = a * one_minus.sqrt();
            meridional.min(tangent_distance)
        }
    };
    if !radial_limit.is_finite() || radial_limit <= 0.0 || h <= -radial_limit {
        return Err(SurfaceError::DegenerateGeometry {
            reason: "height reaches or crosses surface-normal singularity",
        });
    }

    let position: Result<[Meters; 3], SurfaceError> = match (*surface, convention) {
        (ReferenceSurface::Sphere { radius }, _) => {
            let r = radius.value() + h;
            Ok([
                Quantity::<Meter>::new(r * cos_lat * cos_lon),
                Quantity::<Meter>::new(r * cos_lat * sin_lon),
                Quantity::<Meter>::new(r * sin_lat),
            ])
        }
        (ReferenceSurface::Ellipsoid { a, f }, LatitudeConvention::Planetographic) => {
            let a = a.value();
            let e2 = 2.0 * f - f * f;
            let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
            Ok([
                Quantity::<Meter>::new((n + h) * cos_lat * cos_lon),
                Quantity::<Meter>::new((n + h) * cos_lat * sin_lon),
                Quantity::<Meter>::new((n * (1.0 - e2) + h) * sin_lat),
            ])
        }
        (ReferenceSurface::Ellipsoid { a, f }, LatitudeConvention::Planetocentric) => {
            // Planetocentric latitude φ′: radius of the ellipsoid surface
            // along the radius vector is
            //   ρ = ab / sqrt((b cos φ′)² + (a sin φ′)²)
            let a = a.value();
            let b = a * (1.0 - f);
            let denom = ((b * cos_lat).powi(2) + (a * sin_lat).powi(2)).sqrt();
            if denom == 0.0 {
                return Err(SurfaceError::DegenerateGeometry {
                    reason: "ellipsoid radius undefined",
                });
            }
            let rho = (a * b) / denom;
            // The planetocentric angle selects the intersection of the
            // radius vector with the ellipsoid. Height, however, is measured
            // along the same surface normal used by planetographic sites.
            // Thus both latitude conventions describe the same physical
            // location and horizon when their coordinates are converted.
            let normal = normalize([
                cos_lat * cos_lon / (a * a),
                cos_lat * sin_lon / (a * a),
                sin_lat / (b * b),
            ])?;
            Ok([
                Quantity::<Meter>::new(rho * cos_lat * cos_lon + h * normal[0]),
                Quantity::<Meter>::new(rho * cos_lat * sin_lon + h * normal[1]),
                Quantity::<Meter>::new(rho * sin_lat + h * normal[2]),
            ])
        }
    };
    let position = position?;
    if !position.iter().all(|value| value.value().is_finite()) {
        return Err(SurfaceError::NonFiniteInput {
            field: "body_fixed_position",
        });
    }
    Ok(position)
}

/// Local ENU basis at a surface site, expressed in the body-fixed frame.
///
/// At the geographic poles the east basis is formed from the longitude
/// meridian (`east = (−sin λ, cos λ, 0)`), which keeps azimuth well-defined
/// relative to the local meridian rather than returning an arbitrary value.
pub fn local_enu_basis(
    surface: &ReferenceSurface,
    lon: Degrees,
    lat: Degrees,
    convention: LatitudeConvention,
) -> Result<EnuBasis, SurfaceError> {
    validate_angles(lon, lat, Quantity::<Meter>::new(0.0))?;
    surface.validate()?;

    let lon_r = lon.to::<Radian>().value();
    let lat_r = lat.to::<Radian>().value();
    let (sin_lon, cos_lon) = lon_r.sin_cos();
    let (sin_lat, cos_lat) = lat_r.sin_cos();

    // East is always tangential to the parallel of latitude.
    let east = normalize([-sin_lon, cos_lon, 0.0])?;

    let up = match (*surface, convention) {
        (ReferenceSurface::Sphere { .. }, _) => {
            normalize([cos_lat * cos_lon, cos_lat * sin_lon, sin_lat])?
        }
        (ReferenceSurface::Ellipsoid { .. }, LatitudeConvention::Planetographic) => {
            // Geodetic latitude directly specifies the ellipsoid normal.
            normalize([cos_lat * cos_lon, cos_lat * sin_lon, sin_lat])?
        }
        (ReferenceSurface::Ellipsoid { a, f }, LatitudeConvention::Planetocentric) => {
            // Planetocentric latitude locates a point radially, but the
            // geometric local vertical remains the ellipsoid normal.
            let a = a.value();
            let b = a * (1.0 - f);
            normalize([
                cos_lat * cos_lon / (a * a),
                cos_lat * sin_lon / (a * a),
                sin_lat / (b * b),
            ])?
        }
    };

    // Up x East points north and preserves the longitude-defined meridian
    // at the poles, without assuming that planetocentric latitude equals
    // the surface-normal latitude.
    let north = normalize(cross(up, east))?;

    let basis = EnuBasis { east, north, up };
    basis.assert_orthonormal(1e-12)?;
    Ok(basis)
}

fn validate_angles(lon: Degrees, lat: Degrees, height: Meters) -> Result<(), SurfaceError> {
    if !lon.value().is_finite() {
        return Err(SurfaceError::NonFiniteInput { field: "longitude" });
    }
    if !lat.value().is_finite() {
        return Err(SurfaceError::NonFiniteInput { field: "latitude" });
    }
    if !height.value().is_finite() {
        return Err(SurfaceError::NonFiniteInput { field: "height" });
    }
    let lat_v = lat.value();
    if !(-90.0..=90.0).contains(&lat_v) {
        return Err(SurfaceError::InvalidLatitude { degrees: lat_v });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::qtty::{DEG, M};

    #[test]
    fn sphere_equator_position() {
        let s = ReferenceSurface::sphere(1_000.0 * M).unwrap();
        let p = body_fixed_position(
            &s,
            0.0 * DEG,
            0.0 * DEG,
            0.0 * M,
            LatitudeConvention::Planetocentric,
        )
        .unwrap();
        assert!((p[0].value() - 1000.0).abs() < 1e-12);
        assert!(p[1].value().abs() < 1e-12);
        assert!(p[2].value().abs() < 1e-12);
    }

    #[test]
    fn enu_orthonormal_mid_latitude() {
        let s = ReferenceSurface::sphere(1_737_400.0 * M).unwrap();
        let b = local_enu_basis(
            &s,
            45.0 * DEG,
            30.0 * DEG,
            LatitudeConvention::Planetocentric,
        )
        .unwrap();
        b.assert_orthonormal(1e-12).unwrap();
    }

    #[test]
    fn enu_at_north_pole_is_defined() {
        let s = ReferenceSurface::sphere(1_737_400.0 * M).unwrap();
        let b = local_enu_basis(
            &s,
            20.0 * DEG,
            90.0 * DEG,
            LatitudeConvention::Planetocentric,
        )
        .unwrap();
        b.assert_orthonormal(1e-12).unwrap();
        assert!((b.up[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ellipsoid_latitude_conventions_share_the_same_local_vertical() {
        let s = crate::coordinates::surface::MARS_IAU_ELLIPSOID;
        let lon = 137.0 * DEG;
        let graphic = 45.0_f64.to_radians();
        let e2 = s.e2();
        let centric = ((1.0 - e2) * graphic.tan()).atan().to_degrees();
        let g = local_enu_basis(&s, lon, 45.0 * DEG, LatitudeConvention::Planetographic).unwrap();
        let c =
            local_enu_basis(&s, lon, centric * DEG, LatitudeConvention::Planetocentric).unwrap();
        for i in 0..3 {
            assert!((c.up[i] - g.up[i]).abs() < 1e-14);
            assert!((c.north[i] - g.north[i]).abs() < 1e-14);
        }
        let p1 = body_fixed_position(
            &s,
            lon,
            45.0 * DEG,
            100.0 * M,
            LatitudeConvention::Planetographic,
        )
        .unwrap();
        let p2 = body_fixed_position(
            &s,
            lon,
            centric * DEG,
            100.0 * M,
            LatitudeConvention::Planetocentric,
        )
        .unwrap();
        for i in 0..3 {
            assert!((p1[i] - p2[i]).value().abs() < 1e-6);
        }
    }

    #[test]
    fn wgs84_matches_affn_order_of_magnitude() {
        // Greenwich equator on WGS84.
        let s = crate::coordinates::surface::EARTH_WGS84_SURFACE;
        let p = body_fixed_position(
            &s,
            0.0 * DEG,
            0.0 * DEG,
            0.0 * M,
            LatitudeConvention::Planetographic,
        )
        .unwrap();
        assert!((p[0].value() - 6_378_137.0).abs() < 1e-6);
    }

    #[test]
    fn direct_geometry_rejects_singular_interiors_and_bad_angles() {
        let sphere = ReferenceSurface::sphere(1_000.0 * M).unwrap();
        for height in [-1_000.0, -1_001.0] {
            let point = body_fixed_position(
                &sphere,
                0.0 * DEG,
                0.0 * DEG,
                height * M,
                LatitudeConvention::Planetocentric,
            );
            assert!(matches!(point, Err(SurfaceError::DegenerateGeometry { .. })));
        }
        let invalid_lat = local_enu_basis(
            &sphere,
            0.0 * DEG,
            100.0 * DEG,
            LatitudeConvention::Planetocentric,
        );
        assert!(matches!(invalid_lat, Err(SurfaceError::InvalidLatitude { .. })));
        let invalid_lon = local_enu_basis(
            &sphere,
            f64::NAN * DEG,
            0.0 * DEG,
            LatitudeConvention::Planetocentric,
        );
        assert!(matches!(invalid_lon, Err(SurfaceError::NonFiniteInput { .. })));

        let ellipsoid = ReferenceSurface::ellipsoid(4.0 * M, 0.25).unwrap();
        let meridional = 4.0 * (1.0 - ellipsoid.e2());
        let invalid_height = body_fixed_position(
            &ellipsoid,
            0.0 * DEG,
            0.0 * DEG,
            -meridional * M,
            LatitudeConvention::Planetographic,
        );
        assert!(matches!(
            invalid_height,
            Err(SurfaceError::DegenerateGeometry { .. })
        ));
    }

    #[test]
    fn enu_orthonormality_rejects_left_handed_or_nonunit_basis() {
        let sphere = ReferenceSurface::sphere(1_000.0 * M).unwrap();
        let basis = local_enu_basis(
            &sphere,
            0.0 * DEG,
            0.0 * DEG,
            LatitudeConvention::Planetocentric,
        )
        .unwrap();
        let left_handed = EnuBasis {
            east: basis.north,
            north: basis.east,
            up: basis.up,
        };
        assert!(left_handed.assert_orthonormal(1e-12).is_err());
        let nonunit = EnuBasis {
            east: [0.0, 0.0, 0.0],
            north: basis.north,
            up: basis.up,
        };
        assert!(nonunit.assert_orthonormal(1e-12).is_err());
    }
}
