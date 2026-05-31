# siderust-tle test fixtures

All fixtures here are **synthesised** from public-domain reference TLEs
that have been republished by Celestrak for decades (the canonical
"ISS (ZARYA) at epoch 2008-264" sample and the Vanguard 1 sample). They
are not redistributions of any rate-limited or licensed dataset.

Layout:

* `tle/iss_zarya.3le` — classic 3-line TLE.
* `omm/iss_zarya.kvn` — OMM keyword-value notation (CCSDS 502.0).
* `omm/iss_zarya.xml` — OMM XML (CCSDS 502.0 Annex E).
* `omm/iss_zarya.json` — Celestrak-style flat JSON object.
* `omm/celestrak_sample.json` — JSON array combining ISS and Vanguard 1.

Tests in `src/tests.rs::fixtures_round_trip` parse each of these and
assert structural fields. To extend coverage, add new files under
`test-data/` and wire them into a fresh `#[test]`.
