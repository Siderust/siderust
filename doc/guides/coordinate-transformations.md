# Transform a position between coordinate systems

Start by identifying the input origin, axes, unit, and epoch. A frame rotation
preserves the center. A center shift preserves the frame and needs the relative
position of the two origins. Cartesian-to-spherical conversion preserves both.

With default features, this program takes a barycentric position in the mean
ecliptic and equinox of J2000, then expresses it relative to Mars at J2000 TT.
Save it as `src/main.rs` in the [quickstart application](../getting-started/quickstart.md).

```rust
use siderust::coordinates::{
    cartesian::Position,
    centers::{Barycentric, Marscentric},
    frames::EclipticMeanJ2000,
    transform::{try_to_center_with, AstroContext},
};
use siderust::qtty::AstronomicalUnit;
use siderust::J2000;

fn main() {
    let context = AstroContext::new();
    let point = Position::<Barycentric, EclipticMeanJ2000, AstronomicalUnit>::new(1., 0., 0.);
    let shifted: Position<Marscentric, EclipticMeanJ2000, AstronomicalUnit> =
        try_to_center_with(&point, J2000, &context).expect("Mars shift at J2000");
    println!("{shifted}");
}
```

The Mars translation uses the context's ephemeris. In application code,
propagate `CenterShiftError` when the provider, target, or epoch may vary.
`try_to_center_with` is the checked entry point for Cartesian positions with
unit center parameters (`Params = ()`). The `to_center*` conveniences can
panic if a required shift is unavailable.

For a rotation, use `TransformFrame` and let the destination type select the
axes. Supply a TT epoch and a context for time-dependent frames. The
[coordinate transformation example](../../examples/02_coordinate_transformations.rs)
shows both operations with `cargo run --example 02_coordinate_transformations`.
For a site, construct `Geodetic<ECEF>` with east-positive longitude,
geodetic latitude, and height above the WGS84 ellipsoid, in that order.
Its Cartesian conversion gives an Earth-fixed, geocentric position. The
[horizontal example](../../examples/19_horizontal_coordinates.rs) shows
the Earth-specific mean and true equatorial paths. Check
[Earth orientation coverage](../concepts/time-and-earth-orientation.md#earth-rotation-and-pole-motion)
before using them for precision pointing. For a named site, use the
built-in observatory catalog, generated from [observatories.toml](../../data/observatories.toml).
With `serde`, `ObservatoryCatalog::from_toml` accepts that format. Loading
from a path additionally requires `std`.

For a catalog direction with no distance, rotate the direction directly.
For a moving-frame state, use an API that includes the frame-rate and origin
velocity terms needed by the calculation. See [coordinate contracts](../reference/guarantees.md).
