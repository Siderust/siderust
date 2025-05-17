use rayon::join;
use crate::units::JulianDay;

/// A Vsop87 term (coefficient structure).
#[derive(Debug, Clone, Copy)]
pub struct Vsop87 {
    /// The amplitude coefficient.
    pub a: f64,
    /// The phase coefficient.
    pub b: f64,
    /// The frequency coefficient.
    pub c: f64,
}


/// Computes the sum of a Vsop87 series for a specific coordinate (X, Y, or Z).
///
/// `terms` is a slice of terms for a single power of T, like X0 or X1.
/// `t` is the Julian millennia from J2000
#[inline]
fn series_sum(terms: &[Vsop87], t: f64) -> f64 {
    terms.iter()
         .map(|term| term.a * (term.b + term.c * t).cos())
         .sum()
}

/// Computes the expansions for X, Y, or Z across multiple powers of T (X0..X5, etc.).
///
/// `expansions[i]` is the array of Vsop87 terms for T^i.
fn compute_coord_value(expansions: &[&[Vsop87]], t: f64) -> f64 {
    expansions.iter()
        .enumerate()
        .map(|(i, terms)| {
            // Each expansions[i] is multiplied by T^i
            let t_power = t.powi(i as i32);
            series_sum(terms, t) * t_power
        })
        .sum()
}


/// Computes heliocentric rectangular coordinates (X, Y, Z) in AU for a given JD.
///
/// # Arguments
/// - `jd`: Julian Date (TT)
/// - `x_expansions`: arrays for X (X0..X5)
/// - `y_expansions`: arrays for Y (Y0..Y5)
/// - `z_expansions`: arrays for Z (Z0..Z5)
///
/// # Returns
/// CartesianCoord<Heliocentric, Ecliptic>: X, Y, Z in AU
pub fn compute_vsop87(
    jd: JulianDay,
    x_expansions: &[&[Vsop87]],
    y_expansions: &[&[Vsop87]],
    z_expansions: &[&[Vsop87]],
) -> (f64, f64, f64) {

    let jd_tdb = JulianDay::tt_to_tdb(jd);
    let t = jd_tdb.julian_millennias().value();

    let (x, yz) = join(
        || compute_coord_value(x_expansions, t),
        || {
            let y = compute_coord_value(y_expansions, t);
            let z = compute_coord_value(z_expansions, t);
            (y, z)
        }
    );
    let (y, z) = yz;
    (x, y, z)
}
