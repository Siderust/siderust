# Compact ephemeris validation

Measurements from the compact integration at Siderust commit `c15ccf1`.
These are sample maxima over the fixture and its synthetic Earth
orientation inputs.

The [generator](compact_reference.py) uses DE440s and
ERFA. Python keeps the generator on the existing ERFA and jplephem reference
implementations and preserves independence from Siderust.
The SHA-256 of `de440s.bsp` is
`c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2`.

```sh
uv run --python 3.14 --with numpy==2.5.4 --with pyerfa==2.0.1.5 --with jplephem==2.24 \
  scripts/validation/compact_reference.py /path/to/de440s.bsp > /tmp/compact_reference.csv
cmp /tmp/compact_reference.csv tests/fixtures/compact_reference.csv
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

`generate_iau.py` uses the pinned python-sgp4 package to call Vallado's
compiled C++ implementation independently of the Rust SGP4 implementation.
See the [SGP4 fixture notes](../../tests/fixtures/sgp4/README.md) for
the regeneration command and tolerances.
