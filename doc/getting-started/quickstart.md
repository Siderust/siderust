# Compute the Earth–Mars separation

In a Cargo application, run `cargo add siderust` and leave the default
features enabled. To use this checkout, pass its path to
`cargo add siderust --path`. Put the following in `src/main.rs` and run
`cargo run`.

```rust
use siderust::bodies::solar_system::{Earth, Mars};
use siderust::time::J2000;

fn main() {
    let earth = Earth::vsop87a(J2000);
    let mars = Mars::vsop87a(J2000);
    println!("Earth–Mars separation: {}", earth.distance_to(&mars));
}
```

Both positions use the Sun as their origin, the mean ecliptic and equinox of
J2000 as their axes, and astronomical units as their length unit. The result
is the straight-line separation in AU at J2000.0 (JD 2451545.0 TT). This
simultaneous geometric calculation uses bundled VSOP87A coefficients and
requires no JPL kernel. A signal travel-time calculation would need to
account for the different emission and reception epochs.

For a tour of the solar-system accessors, run
`cargo run --example 08_solar_system` from the Siderust checkout. To switch
models, check the [ephemeris guide](../guides/ephemeris-backends.md).
