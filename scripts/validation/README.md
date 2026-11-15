# Compact ephemeris validation

The [siderust-py reference generator](https://github.com/Siderust/siderust-py/tree/fb8069a2d5710d450d40f1abcbbc3df3960a8f76/scripts/validation#compact-ephemeris)
uses DE440s and ERFA via jplephem and PyERFA. Generate or verify the fixture
by following that guide. Run the Rust regression test with:

```sh
cargo test --locked -p siderust --test compact_ephemeris --no-default-features --features embedded-compact
```

The regression requires geometric and annually aberrated pointing within
one arcminute of the reference for each sampled target, epoch and site.
It uses great-circle separation, including near the horizon, and checks
lunar range and Earth barycentric velocity separately. The samples include
both supported interval edges.

The geometric reference uses shared synthetic DUT1 = 0.4 seconds, xp = 0.15″
and yp = −0.25″, including for future dates. This checks transformations
independently of the availability of measured EOP.


## SGP4 vectors

The [siderust-py SGP4 generator](https://github.com/Siderust/siderust-py/tree/fb8069a2d5710d450d40f1abcbbc3df3960a8f76/scripts/validation#sgp4)
calls Vallado's compiled C++ implementation through the pinned python-sgp4
package. See the [SGP4 fixture notes](../../tests/fixtures/sgp4/README.md)
for provenance and tolerances.
