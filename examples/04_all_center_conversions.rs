// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Example: all currently supported center conversions.
//!
//! This demonstrates all center-shift pairs implemented in `providers.rs`:
//! - Barycentric <-> Heliocentric
//! - Barycentric <-> Geocentric
//! - Heliocentric <-> Geocentric
//! - Identity shifts for each center
//!
//! It also demonstrates:
//! - **Bodycentric** conversions: from Barycentric, Heliocentric, and Geocentric into a
//!   body-centric frame (Mars-centric and ISS-centric) with round-trip verification.
//! - **Topocentric** conversions: observer-on-Earth parallax correction applied to
//!   positions originally expressed in each of the three standard centers.

use qtty::*;
use siderust::astro::orbit::Orbit;
use siderust::coordinates::cartesian::Position;
use siderust::coordinates::centers::{
    Barycentric, Bodycentric, BodycentricParams, Geocentric, Geodetic, Heliocentric,
    ReferenceCenter, Topocentric,
};
use siderust::coordinates::frames::{EclipticMeanJ2000, ECEF};
use siderust::coordinates::transform::{CenterShiftProvider, TransformCenter};
use siderust::time::JulianDate;

type F = EclipticMeanJ2000;
type U = AstronomicalUnit;

// ─── Standard center shifts ──────────────────────────────────────────────────

fn show_center_conversion<C1, C2>(jd: &JulianDate, src: &Position<C1, F, U>)
where
    C1: ReferenceCenter<Params = ()>,
    C2: ReferenceCenter<Params = ()>,
    (): CenterShiftProvider<C1, C2, F>,
    (): CenterShiftProvider<C2, C1, F>,
{
    let out: Position<C2, F, U> = src.to_center(*jd);
    let back: Position<C1, F, U> = out.to_center(*jd);
    let err = (*src - back).magnitude();

    println!(
        "{:<12} -> {:<12} out=({:+.9})  roundtrip={:.3e}",
        C1::center_name(),
        C2::center_name(),
        out,
        err
    );
}

// ─── Bodycentric ─────────────────────────────────────────────────────────────

/// Transform `src` (in center `C`) into body-centric coordinates and back.
///
/// Round-trip: C → Bodycentric → Geocentric → C.
fn show_bodycentric_conversion<C>(
    label: &str,
    jd: &JulianDate,
    src: &Position<C, F, U>,
    params: BodycentricParams,
) where
    C: ReferenceCenter<Params = ()>,
    Position<C, F, U>: TransformCenter<Bodycentric, F, U>,
    Position<Bodycentric, F, U>: TransformCenter<Geocentric, F, U>,
    (): CenterShiftProvider<Geocentric, C, F>,
{
    let bary: Position<Bodycentric, F, U> = src.to_center((params, *jd));
    let recovered_geo: Position<Geocentric, F, U> = bary.to_center(*jd);
    let recovered: Position<C, F, U> = recovered_geo.to_center(*jd);
    let err = (*src - recovered).magnitude();

    println!(
        "{:<12} -> {:<12} dist={:.6}  roundtrip={:.3e}",
        label,
        Bodycentric::center_name(),
        bary.distance(),
        err
    );
}

// ─── Topocentric ─────────────────────────────────────────────────────────────

/// Transform a geocentric position to topocentric and back.
///
/// The `label` argument names the original center (before it was shifted to
/// geocentric), so the output shows the full chain (e.g. Barycentric → Geo →
/// Topocentric).
fn show_topocentric_conversion(
    label: &str,
    jd: &JulianDate,
    geo: &Position<Geocentric, F, U>,
    site: Geodetic<ECEF>,
) {
    let topo: Position<Topocentric, F, U> = geo.to_center((site, *jd));
    let recovered: Position<Geocentric, F, U> = topo.to_center(*jd);
    let err = (*geo - recovered).magnitude();

    println!(
        "{:<12} -> {:<12} out=({:+.6})  roundtrip={:.3e}",
        label, "Topocentric", topo, err
    );
}

// ─── main ─────────────────────────────────────────────────────────────────────

fn main() {
    let jd = JulianDate::new(2_460_000.5);
    println!("Center conversion demo at JD(TT) = {:.1}\n", jd);

    let p_bary = Position::<Barycentric, F, U>::new(0.40, -0.10, 1.20);
    let p_helio: Position<Heliocentric, F, U> = p_bary.to_center(jd);
    let p_geo: Position<Geocentric, F, U> = p_bary.to_center(jd);

    // ── Standard center shifts via CenterShiftProvider ────────────────────────
    println!("── Standard center shifts ─────────────────────────────────────────────");

    // Barycentric source
    show_center_conversion::<Barycentric, Barycentric>(&jd, &p_bary);
    show_center_conversion::<Barycentric, Heliocentric>(&jd, &p_bary);
    show_center_conversion::<Barycentric, Geocentric>(&jd, &p_bary);

    // Heliocentric source
    show_center_conversion::<Heliocentric, Heliocentric>(&jd, &p_helio);
    show_center_conversion::<Heliocentric, Barycentric>(&jd, &p_helio);
    show_center_conversion::<Heliocentric, Geocentric>(&jd, &p_helio);

    // Geocentric source
    show_center_conversion::<Geocentric, Geocentric>(&jd, &p_geo);
    show_center_conversion::<Geocentric, Barycentric>(&jd, &p_geo);
    show_center_conversion::<Geocentric, Heliocentric>(&jd, &p_geo);

    // ── Bodycentric: Mars-like orbit (heliocentric reference) ──────────────────
    println!("\n── Bodycentric – Mars-like orbit (heliocentric ref) ───────────────────");
    let mars_orbit = Orbit::new(
        1.524 * AU,
        0.0934,
        Degrees::new(1.85),
        Degrees::new(49.56),
        Degrees::new(286.5),
        Degrees::new(19.41),
        jd,
    );
    let mars_params = BodycentricParams::heliocentric(mars_orbit);

    show_bodycentric_conversion("Heliocentric", &jd, &p_helio, mars_params);
    show_bodycentric_conversion("Barycentric", &jd, &p_bary, mars_params);
    show_bodycentric_conversion("Geocentric", &jd, &p_geo, mars_params);

    // ── Bodycentric: ISS-like orbit (geocentric reference) ────────────────────
    println!("\n── Bodycentric – ISS-like orbit (geocentric ref) ──────────────────────");
    let iss_orbit = Orbit::new(
        0.0000426 * AU, // ~6 378 km in AU
        0.001,
        Degrees::new(51.6),
        Degrees::new(0.0),
        Degrees::new(0.0),
        Degrees::new(0.0),
        jd,
    );
    let iss_params = BodycentricParams::geocentric(iss_orbit);

    show_bodycentric_conversion("Heliocentric", &jd, &p_helio, iss_params);
    show_bodycentric_conversion("Barycentric", &jd, &p_bary, iss_params);
    show_bodycentric_conversion("Geocentric", &jd, &p_geo, iss_params);

    // ── Topocentric: observer at Barcelona ────────────────────────────────────
    println!("\n── Topocentric – Barcelona (lon=2.17°, lat=41.39°, h=12 m) ───────────");
    // Topocentric is geocentric-relative: shift Bary/Helio to Geocentric first,
    // then apply the parallax correction.
    let site = Geodetic::<ECEF>::new(Degrees::new(2.17), Degrees::new(41.39), 12.0 * M);

    let p_geo_from_bary: Position<Geocentric, F, U> = p_bary.to_center(jd);
    let p_geo_from_helio: Position<Geocentric, F, U> = p_helio.to_center(jd);

    show_topocentric_conversion("Barycentric", &jd, &p_geo_from_bary, site);
    show_topocentric_conversion("Heliocentric", &jd, &p_geo_from_helio, site);
    show_topocentric_conversion("Geocentric", &jd, &p_geo, site);
}
