// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Body-fixed Cartesian conversion and local ENU basis.

use super::body::LatitudeConvention;
use super::error::SurfaceError;
use super::reference_surface::ReferenceSurface;
use crate::qtty::{Degrees, Meter, Meters, Quantity, Radian};
#[cfg(not(feature = "std"))]
use crate::qtty::{Real, Scalar, Transcendental};
use core::f64::consts::FRAC_PI_2;

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
/// Height is measured along the local vertical defined by `convention`
/// (ellipsoid normal for planetographic; radius vector for planetocentric
/// on an ellipsoid; identical on a sphere).
pub fn body_fixed_position(
    surface: &ReferenceSurface,
    lon: Degrees,
    lat: Degrees,
    height: Meters,
    convention: LatitudeConvention,
) -> Result<[Meters; 3], SurfaceError> {
    validate_angles(lon, lat, height)?;

    let lon_r = lon.to::<Radian>().value();
    let lat_r = lat.to::<Radian>().value();
    let (sin_lon, cos_lon) = lon_r.sin_cos();
    let (sin_lat, cos_lat) = lat_r.sin_cos();
    let h = height.value();

    match (*surface, convention) {
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
            let r = rho + h;
            Ok([
                Quantity::<Meter>::new(r * cos_lat * cos_lon),
                Quantity::<Meter>::new(r * cos_lat * sin_lon),
                Quantity::<Meter>::new(r * sin_lat),
            ])
        }
    }
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

    let lon_r = lon.to::<Radian>().value();
    let lat_r = lat.to::<Radian>().value();
    let (sin_lon, cos_lon) = lon_r.sin_cos();
    let (sin_lat, cos_lat) = lat_r.sin_cos();

    // East is always tangential to the parallel of latitude.
    let east = normalize([-sin_lon, cos_lon, 0.0])?;

    let up = match (*surface, convention) {
        (ReferenceSurface::Sphere { .. }, _)
        | (ReferenceSurface::Ellipsoid { .. }, LatitudeConvention::Planetocentric) => {
            // Local vertical along the radius vector.
            normalize([cos_lat * cos_lon, cos_lat * sin_lon, sin_lat])?
        }
        (ReferenceSurface::Ellipsoid { a, f }, LatitudeConvention::Planetographic) => {
            // Ellipsoid normal (geodetic up).
            let a = a.value();
            let e2 = 2.0 * f - f * f;
            let n = a / (1.0 - e2 * sin_lat * sin_lat).sqrt();
            // Position on the ellipsoid (h=0), then normal from geodetic formula.
            let x = n * cos_lat * cos_lon;
            let y = n * cos_lat * sin_lon;
            let z = n * (1.0 - e2) * sin_lat;
            // For a rotational ellipsoid the unit normal equals the geodetic
            // (cos φ cos λ, cos φ sin λ, sin φ) direction — identical to the
            // spherical formula when expressed in geodetic latitude.
            let _ = (x, y, z);
            normalize([cos_lat * cos_lon, cos_lat * sin_lon, sin_lat])?
        }
    };

    // North = up × east keeps a right-handed ENU triad.
    // Near poles |cos lat|→0 so east×up from the spherical formula
    //   north = (−sin lat cos lon, −sin lat sin lon, cos lat)
    // is recovered; we compute it explicitly for numerical stability.
    let north = if lat_r.abs() > FRAC_PI_2 - 1e-14 {
        // At the pole: north points along −z for north pole? Convention:
        // at φ=+90°, north is undefined in the classical sense; we adopt
        // the meridian continuation: north = (−cos λ, −sin λ, 0) so that
        // azimuth 0° follows the longitude meridian toward decreasing z
        // when leaving the north pole along that meridian... Actually the
        // standard IAU/mission convention at the north pole uses:
        //   north = (−cos λ, −sin λ, 0)  (toward lon along decreasing z from +z)
        // Wait: at φ=90°, spherical north = (0,0,0) from the formula.
        // Longitude-dependent meridian: azimuth 0 points toward lon along
        // the great circle. Using north = (−cos λ, −sin λ, 0) makes
        // east×north = up for up=(0,0,1):
        // east=(−sin λ, cos λ, 0), north=(−cos λ, −sin λ, 0)
        // east×north = (0, 0, sin²λ + cos²λ) = (0,0,1). ✓
        if lat_r > 0.0 {
            normalize([-cos_lon, -sin_lon, 0.0])?
        } else {
            // South pole up=(0,0,−1): north = (cos λ, sin λ, 0)
            normalize([cos_lon, sin_lon, 0.0])?
        }
    } else {
        normalize([-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat])?
    };

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
}
