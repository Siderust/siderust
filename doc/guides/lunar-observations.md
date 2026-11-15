# Observe from the lunar surface

Use the analytical surface pipeline for a geometric line of sight on the
IAU lunar reference sphere without external kernels. From the checkout, run:

```sh
cargo run --example 22_moon_surface_quickstart
cargo run --example 23_lunar_landing_site_navigator -- -89.5 0 2461358.5
```

The [first example](../../examples/22_moon_surface_quickstart.rs) rotates an
ICRS direction into the Moon's body-fixed frame at J2000 TT and projects it
onto the site's east, north, and up axes. The direction has no distance, so
the calculation includes no parallax. The
[second example](../../examples/23_lunar_landing_site_navigator.rs) samples
finite Sun and Earth positions from the supplied latitude, east longitude,
and starting JD TT. A finite target is translated to the Moon's center,
rotated into lunar axes, and differenced with the site's position before
projection. This includes observer parallax. Visibility means that the
target center lies above the geometric horizon. Model terrain, disc
extent, aberration, light-time, and refraction separately where applicable.

## Select the site datum and orientation

`SurfaceSite::<Moon>::new` takes longitude, latitude, and height in that
order. The default site uses east-positive longitude, planetocentric
latitude, and height in metres above a sphere of radius 1737.4 km. On an
ellipsoid, planetocentric latitude uses the radius vector, while
planetographic latitude and local up use the surface normal. Preserve this
distinction when importing site coordinates.

The analytical `MoonIauMeanEarth` provider uses the NAIF `IAU_MOON`
orientation, including the periodic terms and the quadratic prime
meridian term. It converts TT to TDB before evaluating the lunar arguments.
`Moon::ROTATION` stores the secular elements only. `MoonPrincipalAxes`
is a distinct marker without an analytical provider. For landing
coordinates or maps defined in a DE-specific principal-axes or Mean-Earth
realization, use [matching SPICE kernels](spice-kernels.md).

## Apply the same geometry to other bodies

The [surface example](../../examples/21_body_surface_horizontal.rs) also
uses Mars and Earth sites. Earth uses WGS84 and its Earth orientation chain.
Other bodies use their own pole and prime meridian models. The nonlunar
planetary providers use linear elements with TT as an approximation to the
orientation time argument, so they omit periodic terms found in some NAIF
models.

Mars defaults to a planetographic site on the IAU ellipsoid (equatorial
radius 3396.19 km, polar radius 3376.20 km). Mercury, Venus, and Pluto use
mean spheres. A body-fixed frame alone does not imply a default surface
model. In particular, the gas giants lack a default `BodySurfaceModel`.

Siderust uses east-positive longitude, matching SPICE planetocentric
longitude. SPICE planetographic conventions can differ by body, for example
west-positive on Mars by default. Use
`SurfaceSite::from_planetographic_longitude` for the body's default external
sense or specify an explicit sense for a product that uses a different one.
The API does not read SPICE kernel-pool longitude overrides. See the
[NAIF planetary constants and orientation models](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00011.tpc).
