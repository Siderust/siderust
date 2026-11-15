# Compact ephemeris validation

Measurements from the compact integration at Siderust commit `c15ccf1`.
These are sample maxima over the fixture and its synthetic Earth
orientation inputs.

The [siderust-py reference generator](https://github.com/Siderust/siderust-py/tree/fb8069a2d5710d450d40f1abcbbc3df3960a8f76/scripts/validation#compact-ephemeris)
uses DE440s and ERFA via jplephem and PyERFA. Generate or verify the fixture
by following that guide. Run the Rust regression test with:

```sh
cargo test --locked -p siderust --test compact_ephemeris --no-default-features --features embedded-compact -- --nocapture
```

Across 486 body/site/epoch cases, including both ends of the supported
interval and targets within 0.051° of the horizon:

| Maximum discrepancy | Measured |
| --- | ---: |
| Geocentric angular separation | 3.89″ |
| Topocentric angular separation | 4.16″ |
| Raw azimuth difference | 23.57″ |
| Altitude difference | 3.99″ |
| Topocentric separation with annual aberration | 4.16″ |
| Lunar radial distance | 0.91 km |
| Earth barycentric velocity | 0.182 m/s |
| First-order aberration bound from velocity difference | 0.000126″ |
| TT–TDB offset difference against ERFA | 8.02 µs |
| IAU 2000B vs ERFA IAU 2006A longitude / obliquity | 1.73 / 0.65 mas |
| Changing synthetic EOP to UT1=UTC and zero polar motion | 6.07″ |

The geometric reference uses shared synthetic DUT1 = 0.4 seconds, xp = 0.15″
and yp = −0.25″, including for future dates. This checks transformations
independently of the availability of measured EOP.


## SGP4 vectors

The [siderust-py SGP4 generator](https://github.com/Siderust/siderust-py/tree/fb8069a2d5710d450d40f1abcbbc3df3960a8f76/scripts/validation#sgp4)
calls Vallado's compiled C++ implementation through the pinned python-sgp4
package. See the [SGP4 fixture notes](../../tests/fixtures/sgp4/README.md)
for provenance and tolerances.
