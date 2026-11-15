# Lunar position and nutation interpolation

A long lunar search evaluates a topocentric altitude signal many times.
[MoonAltitudeContext](../../src/event/lunar/moon_cache.rs) precomputes geocentric
lunar Cartesian position from the configured default ephemeris (full ELP2000-82B
in the standard benchmark build). Each four-day segment uses a degree-eight
Chebyshev fit. The nine cosine-spaced nodes are mapped from `[-1, 1]` to the
segment's time interval. A cosine transform yields coefficients for each
axis and Clenshaw recurrence evaluates the series in time linear in its degree.
The same module caches IAU 2000B nutation longitude, obliquity and mean obliquity
on a two-hour grid and linearly interpolates those parameters.

At each query the context still applies observer-dependent WGS84 parallax,
Earth rotation, precession and horizontal projection. Cached values outside
their usable domain fall back to the direct models. The context is built per
search. The two thresholds of an altitude-band search share one context.

The position cache and crossing search approximate different functions.
The former approximates the three Cartesian components. The latter fits the
cached evaluator's altitude residual, checks candidate crossings against that
evaluator and can revert to scan and Brent for difficult segments. See
[lunar event entry points](../../src/event/lunar/altitude_periods.rs) and the
[shared crossing search](../../src/event/search/crossings.rs).

The [cache tests](../../src/event/lunar/moon_cache.rs) assert maximum-axis
position error below 1 km and nutation parameter error below `5e-10` radians
at sampled times over 30 days near J2000. A seven-day comparison at Roque
asserts altitude differences below 0.02° against the direct evaluator.
These are sampled limits against the underlying library models. Near a
crossing, the timing effect of an altitude error depends on the local slope.
A shallow or grazing crossing requires separate validation.

The [lunar benchmark report](moon_altitude_perf_report.md) describes the current
search-engine experiment.
