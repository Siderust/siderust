# Time and Earth orientation

The time **scale** determines what a time coordinate means. Julian Date (JD)
and Modified Julian Date (MJD) are formats, so converting JD to MJD leaves
the scale unchanged. Siderust's `JulianDate` and `ModifiedJulianDate` aliases
default to Terrestrial Time (TT). Most astronomical transforms and ephemeris
queries take TT. The lunar SPICE observer takes `Time<TDB>` in its default
J2000-seconds format.

| Scale | Role |
| --- | --- |
| UTC | Coordinated Universal Time, interpreted with leap-second and historical offset data. |
| TAI | International Atomic Time, a continuous scale. TT is offset from TAI by 32.184 seconds. |
| TT | Time argument for terrestrial astrometric models and most Siderust entry points. |
| UT1 | Universal Time tied to Earth rotation. DUT1 means UT1 minus UTC. |
| TDB | Barycentric Dynamical Time, used in planetary and lunar ephemeris evaluation, analytical lunar orientation, and SPICE ephemeris time. |

Use [Tempoch's typed conversions](https://docs.rs/tempoch) to move between
scales. Use a `TimeContext` for civil JD/MJD encoding and context-dependent
UT1 conversions. Civil UTC JD/MJD use nominal 86,400-second days and cannot
encode a leap-second instant, even though Tempoch can represent that instant
internally. For example, MJD 57754.0 UTC is 2017-01-01 midnight, 86,401 SI
seconds after the previous UTC midnight. TLE propagation uses nominal civil
minutes (1,440 per day). Event-search `Period<MJD>` inputs on the Siderust
altitude API use TT.

## Earth rotation and pole motion

The celestial-to-terrestrial chain moves from GCRS through CIRS and TIRS to
terrestrial axes. Precession and nutation determine the celestial pole and
origin, with Earth orientation parameters (EOP) providing corrections
`dX` and `dY`. The Earth Rotation Angle uses UT1. Polar motion uses
`xp`, `yp` and the modeled terrestrial origin locator. An equinox-based
path uses sidereal time. Match mean or apparent sidereal time to the
coordinate frame expected by the method. The explicit UT1/TT
horizontal traits use equinox-based sidereal rotation. The context-aware
horizontal projection and its inverse share the full Earth-orientation
chain with observer parallax. Observed `dX` and `dY` correct the modeled
celestial pole in the CIO-based chain, while the equinox-of-date axes
follow the selected precession-nutation model. See
[IERS Conventions, Chapter 5](https://iers-conventions.obspm.fr/chapter5.php).

`AstroContext` supplies the ephemeris and EOP providers. `with_model`
selects a nutation model for APIs that accept a `ModelContext`. The default
`IersEop` reads Tempoch's active time-data bundle and queries daily EOP by
UTC MJD. To load daily tables through Tempoch, enable
`tempoch/runtime-data-fetch` in the application and call
`tempoch::update_runtime_time_data` before evaluating. This may fetch a bundle
if a valid cache is unavailable. Siderust's `runtime-data` feature handles
ephemeris acquisition. With the bundled time data alone, the checked EOP
lookup reports missing data. When daily values are unavailable, infallible lookup
derives DUT1 from modeled ΔT (TT minus UT1) and sets polar motion and
celestial pole offsets to zero. Infallible lookup can panic if the time model cannot evaluate
the epoch. For measured or predicted daily EOP, check
`AstroContext::try_eop_at_tt`, table provenance, and field completeness.
Missing optional pole fields can become zero even on successful lookup.
`NullEop` explicitly uses UT1 ≈ UTC and zero pole corrections. DUT1 alone
can matter at arcsecond scales because Earth rotates approximately 15 arcseconds
per second of UT1.

`IersEop` reads active data at lookup time. For a reproducible sequence,
keep the active data fixed during the run or supply an EOP provider backed
by your own pinned data. Record whether daily values, predictions, a ΔT
model, or the zero-EOP approximation were used.
