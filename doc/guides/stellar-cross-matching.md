# Cross-match stellar catalogs across epochs

Reference epochs matter when comparing stellar catalogs. A star with a proper
motion of 1.8 arcseconds per year shifts by about 29 arcseconds between J2000.0
and J2016.0. A sub-arcsecond search at the stored positions can miss that source.

Run the [Gaia-style example](../../examples/26_stellar_cross_match.rs) with the
default build. It uses synthetic CSV data embedded in the source.

```sh
cargo run --example 26_stellar_cross_match
```

The example searches around the Gaia-style directions at J2000.0. One source
crosses 0° right ascension during the backward propagation and has two nearby
historical candidates. A second has one candidate, and a third emits no rows.

## Choose an epoch, radius and policy

Call [`LargeStarCatalog::cross_match`][cross-match] with the other catalog and
[`CrossMatchOptions`][options]. Supply a common reference epoch as
`qtty::JulianYears` and an inclusive radius as `qtty::Radians`. A quantity such
as `Arcseconds::new(0.5).to()` converts directly to the expected radius type.
The epochs must use a consistent time scale. Gaia DR3 publishes its reference
epoch in TCB. Prepare the input directions in ICRS and express RA and
declination in radians, even if the input CSV used degrees.

At the comparison epoch, a row can be used directly, including a Gaia
position-only solution. Away from that epoch, a row with either proper-motion
component missing follows the selected [`MissingProperMotion`][missing] policy.

- `Skip` leaves the row out of the candidate search.
- `AssumeStationary` uses the recorded position and ignores any available
  partial component. The candidate endpoint carries `assumed_stationary = true`.
- `Error` returns a [`CrossMatchError`][error] identifying the catalog side and
  row index.

The matcher rejects invalid coordinates, non-finite epochs or proper motion,
and invalid radii, even when the opposing catalog is empty. The valid radius
range is zero through π radians.

## Interpret the candidates

Each [`CrossMatchCandidate`][candidate] contains both input row indices, their
catalog-supplied IDs, and the great-circle separation in radians. The row
indices identify the records within their respective input catalogs, including
when IDs repeat.

`AllCandidates` reports every pair within the radius. `Nearest` returns one
pair per matched left row, resolving exact numerical ties by lowest right row
index. Both modes order by left row index, separation, then right row index.
Each emitted pair carries the number of candidates for its left and right
endpoints in the full search. Counts above one indicate geometric ambiguity,
including choices hidden by `Nearest`. Several left rows may select the same
right row. Assess physical identity separately using evidence such as
astrometric uncertainties, photometry, or independent observations.

## Propagation model and search cost

[`CatalogRecord::propagate_to`][propagate] uses the Gaia convention
`pm_ra_cosdec = μᵅ⋆ = (dα/dt) cos δ` and `pm_dec = μᵟ` in milliarcseconds per
Julian year. It adds the linear east and north displacements to the reference
unit vector and takes the direction of the result. This is the transverse
case of the rectilinear model described in the
[Gaia DR3 astrometry documentation][gaia-model]. The stored RA sets the
orientation of the tangent basis at an exact pole. Angular separation uses
the spherical Vincenty formula supplied by `affn`.

The model omits radial velocity, parallax, light-time effects, and uncertainty
propagation. Perspective acceleration can matter for nearby high-velocity
stars and long epoch baselines. Set the radius using the uncertainty and
resolution of the input catalogs and the suitability of the motion model for
the epoch interval.

Rows are evaluated once per call. The search sorts the right catalog by
declination, then tests full spherical separations only in the declination
window of each left row. With L left rows, R right rows, and P pairs in the
windows, this costs O(R log(R+1) + L log(R+2) + P) time, plus sorting each
left row's candidate list in all-candidates mode. Storage is O(R) plus the
output. Broad windows or sources concentrated in the same declination strip
can lead to O(LR) pair tests.

[cross-match]: https://docs.rs/siderust/latest/siderust/catalogs/struct.LargeStarCatalog.html#method.cross_match
[options]: https://docs.rs/siderust/latest/siderust/catalogs/struct.CrossMatchOptions.html
[candidate]: https://docs.rs/siderust/latest/siderust/catalogs/struct.CrossMatchCandidate.html
[missing]: https://docs.rs/siderust/latest/siderust/catalogs/enum.MissingProperMotion.html
[error]: https://docs.rs/siderust/latest/siderust/catalogs/enum.CrossMatchError.html
[propagate]: https://docs.rs/siderust/latest/siderust/catalogs/struct.CatalogRecord.html#method.propagate_to
[gaia-model]: https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu3ast/sec_cu3ast_intro/ssec_cu3ast_intro_tansforms.html
