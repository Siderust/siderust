# Lunar position interpolation

[MoonAltitudeContext](../../src/event/lunar/moon_cache.rs) precomputes geocentric
lunar Cartesian position from the configured default ephemeris (full ELP2000-82B
in the standard benchmark build). Each four-day segment uses a degree-eight
Chebyshev fit. The nine cosine-spaced nodes are mapped from `[-1, 1]` to the
segment's time interval. A cosine transform yields coefficients for each
axis and Clenshaw recurrence evaluates the series in time linear in its degree.

Each query rotates the cached position into the Earth-fixed frame using the
configured Earth orientation, subtracts the WGS84 site, and projects onto
the local horizon. The default context supplies the Earth-orientation and
nutation models at each TT query. Queries outside the cache use the direct
ephemeris. The two thresholds of an altitude-band search share one context.

The position cache and crossing search approximate different functions.
The former approximates the three Cartesian components. The latter fits the
cached evaluator's altitude residual, checks candidate crossings against that
evaluator and can revert to scan and Brent for difficult segments. See
[lunar event entry points](../../src/event/lunar/altitude_periods.rs) and the
[shared crossing search](../../src/event/search/crossings.rs).

The [cache tests](../../src/event/lunar/moon_cache.rs) compare sampled
positions and altitudes with direct evaluations of the same ephemeris.
Interpolation error and the local altitude slope affect crossing times.
The root solver's time tolerance applies to the interpolated signal.

The [lunar benchmark workloads](moon_altitude_perf_report.md) exercise cache
construction, repeated altitude evaluations, and event searches.
