//! HEALPix RING-ordering numerical kernels.
//!
//! This private module implements the RING indexing formulae used by
//! [`HealpixGrid`](crate::healpix::HealpixGrid). The public API remains based on
//! typed Cartesian directions; this module converts unit vectors to the HEALPix
//! spherical variables `phi` and `z = cos(theta)` before applying the RING
//! indexing equations.
//!
//! # References
//!
//! - Gorski, K. M. et al. (2005), "HEALPix: A Framework for High-Resolution
//!   Discretization and Fast Analysis of Data Distributed on the Sphere",
//!   Astrophysical Journal, 622, 759.
//! - HEALPix Primer, section on RING and NESTED pixel numbering schemes.

use crate::coordinates::cartesian::Direction;
use crate::coordinates::frames::ReferenceFrame;
use crate::healpix::{HealpixError, HealpixGrid, HealpixIndex, HealpixOrdering, Result};
use std::f64::consts::{FRAC_PI_2, TAU};

pub(super) fn unit_vector_to_pixel(grid: &HealpixGrid, xyz: [f64; 3]) -> Result<HealpixIndex> {
    let [x, y, z] = xyz;
    let radius = x.hypot(y).hypot(z);
    if !radius.is_finite() || radius == 0.0 {
        return Err(HealpixError::InvalidAngles {
            reason: "direction vector must be finite and non-zero",
        });
    }
    let phi = y.atan2(x).rem_euclid(TAU);
    let z = (z / radius).clamp(-1.0, 1.0);
    phi_z_to_pixel(grid, phi, z)
}

fn phi_z_to_pixel(grid: &HealpixGrid, phi: f64, z: f64) -> Result<HealpixIndex> {
    if grid.ordering() != HealpixOrdering::Ring {
        return Err(HealpixError::UnsupportedOrdering(grid.ordering()));
    }
    if !phi.is_finite() || !z.is_finite() {
        return Err(HealpixError::InvalidAngles {
            reason: "HEALPix phi and z must be finite",
        });
    }
    if !(-1.0..=1.0).contains(&z) {
        return Err(HealpixError::InvalidAngles {
            reason: "HEALPix z must lie in [-1, 1]",
        });
    }

    let nside = u64::from(grid.nside().get());
    let nside_f = nside as f64;
    let nl4 = 4 * nside;
    let ncap = 2 * nside * (nside - 1);
    let phi = phi.rem_euclid(TAU);
    let za = z.abs();
    let tt = phi / FRAC_PI_2;

    let index = if za <= 2.0 / 3.0 {
        let jp = (nside_f * (0.5 + tt - 0.75 * z)).floor() as i64;
        let jm = (nside_f * (0.5 + tt + 0.75 * z)).floor() as i64;
        let ir = i64::try_from(nside).expect("nside fits i64") + 1 + jp - jm;
        let kshift = 1 - (ir & 1);
        let mut ip = (jp + jm - i64::try_from(nside).expect("nside fits i64") + kshift + 1) / 2 + 1;
        let nl4_i = i64::try_from(nl4).expect("nl4 fits i64");
        if ip > nl4_i {
            ip -= nl4_i;
        }
        if ip < 1 {
            ip += nl4_i;
        }
        ncap + u64::try_from(ir - 1).expect("equatorial ring is positive") * nl4
            + u64::try_from(ip - 1).expect("equatorial pixel is positive")
    } else {
        let tp = tt.fract();
        let tmp = nside_f * (3.0 * (1.0 - za)).sqrt();
        let jp = (tp * tmp).floor() as u64;
        let jm = ((1.0 - tp) * tmp).floor() as u64;
        let ir = jp + jm + 1;
        let mut ip = (tt * ir as f64).floor() as u64 + 1;
        let ring_pixels = 4 * ir;
        if ip > ring_pixels {
            ip -= ring_pixels;
        }
        if z >= 0.0 {
            2 * ir * (ir - 1) + ip - 1
        } else {
            grid.npix() - 2 * ir * (ir + 1) + ip - 1
        }
    };

    let pixel = HealpixIndex::new(index);
    grid.validate_index(pixel)?;
    Ok(pixel)
}

pub(super) fn pixel_to_theta_phi(grid: &HealpixGrid, index: HealpixIndex) -> Result<(f64, f64)> {
    if grid.ordering() != HealpixOrdering::Ring {
        return Err(HealpixError::UnsupportedOrdering(grid.ordering()));
    }
    grid.validate_index(index)?;

    let pix = index.get();
    let nside = u64::from(grid.nside().get());
    let nside_f = nside as f64;
    let nl4 = 4 * nside;
    let ncap = 2 * nside * (nside - 1);
    let equatorial_end = ncap + nl4 * (2 * nside + 1);
    let fact1 = 1.0 / (1.5 * nside_f);
    let fact2 = 1.0 / (3.0 * nside_f * nside_f);

    let (phi, z) = if pix < ncap {
        let iring = ((1.0 + (1.0 + 2.0 * pix as f64).sqrt()) * 0.5).floor() as u64;
        let iphi = pix + 1 - 2 * iring * (iring - 1);
        let z = 1.0 - iring as f64 * iring as f64 * fact2;
        let phi = (iphi as f64 - 0.5) * FRAC_PI_2 / iring as f64;
        (phi, z)
    } else if pix < equatorial_end {
        let ip = pix - ncap;
        let iring = ip / nl4 + nside;
        let iphi = ip % nl4 + 1;
        let fodd = if ((iring + nside) & 1) == 1 { 1.0 } else { 0.5 };
        let z = (2.0 * nside_f - iring as f64) * fact1;
        let phi = (iphi as f64 - fodd) * FRAC_PI_2 / nside_f;
        (phi, z)
    } else {
        let ip = grid.npix() - pix;
        let iring = ((1.0 + (2.0 * ip as f64 - 1.0).sqrt()) * 0.5).floor() as u64;
        let iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
        let z = -1.0 + iring as f64 * iring as f64 * fact2;
        let phi = (iphi as f64 - 0.5) * FRAC_PI_2 / iring as f64;
        (phi, z)
    };

    Ok((z.clamp(-1.0, 1.0).acos(), phi.rem_euclid(TAU)))
}

pub(crate) fn direction_from_theta_phi<F>(theta: f64, phi: f64) -> Direction<F>
where
    F: ReferenceFrame,
{
    let sin_theta = theta.sin();
    Direction::<F>::from_array([sin_theta * phi.cos(), sin_theta * phi.sin(), theta.cos()])
}
