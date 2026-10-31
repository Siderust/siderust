//! Print geometric ephemeris states for independent frame validation.
#![allow(clippy::print_stdout)]
use siderust::bodies::solar_system::{Earth, Mars, Mercury, Moon, Neptune};
use siderust::qtty::{AstronomicalUnit, Kilometer};
use siderust::time::JulianDate;
fn main() {
    for jd in [
        2396759.5, 2415020.0, 2451545.0, 2461344.5, 2488070.0, 2506000.0,
    ] {
        let epoch = JulianDate::new(jd);
        macro_rules! planet {
            ($body:ident, $id:expr) => {{
                let (p, v) = $body::vsop87a_pos_vel(epoch);
                let p = p.to_unit::<Kilometer>();
                let v = v.to_unit::<siderust::qtty::Per<Kilometer, siderust::qtty::Day>>();
                println!(
                    "{jd},{},{},{},{},{},{},{}",
                    $id,
                    p.x().value(),
                    p.y().value(),
                    p.z().value(),
                    v.x().value(),
                    v.y().value(),
                    v.z().value()
                );
            }};
        }
        planet!(Mercury, 1);
        planet!(Earth, 399);
        planet!(Mars, 4);
        planet!(Neptune, 8);
        let p = Moon::get_geo_position::<AstronomicalUnit>(epoch).to_unit::<Kilometer>();
        // ELP exposes positions only. Do not manufacture a high precision
        // velocity from finite differences of a one-part Julian Date.
        println!(
            "{jd},301,{},{},{},NaN,NaN,NaN",
            p.x().value(),
            p.y().value(),
            p.z().value()
        );
    }
}
