// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

#![no_std]
#![no_main]

use core::{
    alloc::{GlobalAlloc, Layout},
    hint::black_box,
    panic::PanicInfo,
};

use esp_hal::{clock::CpuClock, time::Instant};
use esp_println::println;
use siderust::astro::{eop::NullEop, nutation::Iau2000B};
use siderust::coordinates::{
    cartesian::Position,
    centers::{Geocentric, Geodetic},
    frames::{EclipticMeanJ2000, ECEF},
    transform::{AstroContext, TransformFrame},
};
use siderust::ephemeris::{CompactEphemeris, Ephemeris, EphemerisError, MajorPlanet};
use siderust::event::{
    horizontal::{
        equatorial_to_horizontal_true_of_date_with_ctx,
        geocentric_j2000_to_apparent_topocentric_with_ctx,
    },
    lunar::{
        phase::moon_phase_from_topocentric, position_angle::moon_bright_limb_from_topocentric,
    },
};
use siderust::qtty::{AstronomicalUnit, Degrees, Meters};
use siderust::JulianDate;

esp_bootloader_esp_idf::esp_app_desc!();

struct NoHeap;

// SAFETY: a null allocation is a valid failure. This allocator never hands
// out any pointers and needs no heap reservation.
unsafe impl GlobalAlloc for NoHeap {
    unsafe fn alloc(&self, _: Layout) -> *mut u8 {
        core::ptr::null_mut()
    }
    unsafe fn dealloc(&self, _: *mut u8, _: Layout) {}
}

#[global_allocator]
static ALLOCATOR: NoHeap = NoHeap;

#[panic_handler]
fn panic(info: &PanicInfo<'_>) -> ! {
    println!("{info}");
    halt()
}

fn halt() -> ! {
    loop {
        core::hint::spin_loop();
    }
}

fn observe(jd: JulianDate) -> Result<(), EphemerisError> {
    let site = Geodetic::<ECEF>::new(
        Degrees::new(-0.1275),
        Degrees::new(51.507222),
        Meters::new(20.0),
    );
    let ctx = AstroContext::<CompactEphemeris, NullEop>::default();
    let ctx = ctx.with_model::<Iau2000B>();
    let earth = CompactEphemeris::try_earth_barycentric(jd)?;
    let sun = CompactEphemeris::try_sun_barycentric(jd)?;
    let geocentric = |body: siderust::coordinates::cartesian::Position<_, _, AstronomicalUnit>| {
        Position::<Geocentric, EclipticMeanJ2000, AstronomicalUnit>::new(
            body.x() - earth.x(),
            body.y() - earth.y(),
            body.z() - earth.z(),
        )
    };
    let moon = CompactEphemeris::try_moon_geocentric(jd)?;
    let moon = geocentric_j2000_to_apparent_topocentric_with_ctx(&moon.to_frame(), site, jd, &ctx);
    let sun = geocentric_j2000_to_apparent_topocentric_with_ctx(
        &geocentric(sun).to_frame(),
        site,
        jd,
        &ctx,
    );
    black_box(equatorial_to_horizontal_true_of_date_with_ctx(
        &moon, site, jd, &ctx,
    ));
    black_box(equatorial_to_horizontal_true_of_date_with_ctx(
        &sun, site, jd, &ctx,
    ));
    black_box(moon_phase_from_topocentric(&moon, &sun));
    black_box(moon_bright_limb_from_topocentric(&moon, &sun));
    for planet in [
        MajorPlanet::Mercury,
        MajorPlanet::Venus,
        MajorPlanet::Mars,
        MajorPlanet::Jupiter,
        MajorPlanet::Saturn,
        MajorPlanet::Uranus,
        MajorPlanet::Neptune,
    ] {
        let body = match CompactEphemeris::try_planet_barycentric(planet, jd) {
            Ok(body) => body,
            Err(error) => panic!("{error}"),
        };
        let topo = geocentric_j2000_to_apparent_topocentric_with_ctx(
            &geocentric(body).to_frame(),
            site,
            jd,
            &ctx,
        );
        black_box(equatorial_to_horizontal_true_of_date_with_ctx(
            &topo, site, jd, &ctx,
        ));
    }
    Ok(())
}

#[esp_hal::main]
fn main() -> ! {
    let _peripherals = esp_hal::init(esp_hal::Config::default().with_cpu_clock(CpuClock::_80MHz));
    let mut samples = [0; 7];
    for sample in &mut samples {
        let start = Instant::now();
        for index in 0..128 {
            let jd = JulianDate::new(black_box(2_458_850.0 + 88.0 * f64::from(index)));
            if let Err(error) = observe(jd) {
                println!("{error}");
                halt();
            }
        }
        *sample = start.elapsed().as_micros();
    }
    samples.sort_unstable();
    println!("80 MHz, microseconds per 128 epochs: {samples:?}");
    halt()
}
