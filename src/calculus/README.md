# `calculus`

## Scientific scope

This module is the numerical backbone of `siderust`. It contains the
ephemeris evaluators, orbital solvers, event-search kernels, and reusable
root-finding / interval machinery that higher-level modules depend on. If
`astro` expresses standards and reduction models, `calculus` expresses the
actual numerical procedures used to evaluate them over time.

The scientific range is broad: analytical planetary and lunar theories,
optional JPL DE4xx backends, altitude and azimuth event searches, lunar phase
and photometry helpers, and generic mathematical primitives for solving and
bracketing one-dimensional astronomical problems.

## Technical scope

What you will find here:

- `vsop87/`: analytical planetary series evaluation.
- `elp2000/`: analytical lunar theory evaluation.
- `jpl/` and `ephemeris/`: JPL DE440/DE441 and backend abstraction.
- `solar/`, `lunar/`, `stellar/`: body-specific event and geometry helpers.
- `altitude/`, `azimuth/`: unified event-search APIs and provider traits.
- `conic_equations.rs`: heliocentric orbit wrappers around reusable
  `keplerian` anomaly solvers.
- `math_core/`: reusable Brent, bracketing, extrema, and interval logic.
- `horizontal.rs`, `pluto.rs`: focused numerical helpers used elsewhere.

The public APIs remain typed, but the internal kernels are the sanctioned
location where raw `f64` math is performed for speed and parity with published
series coefficients.

## References

- Bretagnon, P., and Francou, G. (1988). "Planetary theories in rectangular
  and spherical variables: VSOP87 solution." Astronomy and Astrophysics 202,
  309-315.
- Chapront-Touze, M., and Chapront, J. (1988). "ELP 2000-85: A semi-analytical
  lunar ephemeris adequate for historical times." Astronomy and Astrophysics
  190, 342-352.
- Park, R. S., Folkner, W. M., Williams, J. G., and Boggs, D. H. (2021). "The
  JPL Planetary and Lunar Ephemerides DE440 and DE441." Astronomical Journal
  161, 105.
- NASA NAIF. SPK Required Reading.
  <https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html>
