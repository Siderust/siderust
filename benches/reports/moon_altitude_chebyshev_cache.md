# Moon Altitude Speedup (365-day): Chebyshev Ephemeris Cache + Nutation Interpolation

> Note: This report lives under `benches/reports/` because it documents
> benchmark-driven performance work.

This document explains the method implemented in:

- `../../src/calculus/lunar/moon_cache.rs`
- `../../src/calculus/lunar/altitude_periods.rs`

The goal is to keep long-horizon Moon altitude window finding fast by replacing
“recompute a heavy ephemeris at every sample time” with “precompute and
interpolate smoothly within short segments”, while preserving the physical
meaning of “topocentric Moon altitude”.

---

## 1) Problem: why `find_moon_above_horizon_365day` was slow

The “Moon above horizon for 365 days” workflow is an **interval finding** problem:

1. Scan the period at a coarse step (2 hours) and detect sign changes of
   \(f(t) - h_{thr}\), where \(f(t)\) is topocentric altitude and \(h_{thr}\) is the threshold (typically 0°).
2. Refine each bracketed root with Brent’s method.
3. Assemble [enter, exit] time intervals.

The expensive part is **evaluating \(f(t)\)**. In the uncached implementation,
each altitude evaluation runs the full chain: ELP2000-82B lunar theory (the
dominant cost), then a sequence of transforms (ecliptic→equatorial J2000, an
observer-dependent topocentric parallax translation using Earth-rotation terms,
IAU 2006 precession, IAU 2000B nutation) and finally equatorial→horizontal to
extract altitude.

Even with the improved 2-hour scan step, a 365-day search still requires **thousands** of altitude evaluations (scan points + Brent refinement + classification). That dominated runtime.

---

## 2) New idea: treat the Moon like an ephemeris table (interpolate, don’t recompute)

The key observation:

- Over hours-to-days, the Moon’s geocentric position varies **smoothly**.
- The interval finder needs altitude at many closely spaced times.

So we build a fast evaluator:

1. **Precompute** Moon geocentric ecliptic \(\vec r(t) = (x,y,z)\) at carefully chosen nodes.
2. **Interpolate** \(\vec r(t)\) at arbitrary times during the search.
3. Run the same physical transforms (parallax / precession / nutation / horizontal) using the interpolated \(\vec r(t)\).

This is conceptually similar to how high-end ephemerides (e.g., JPL DE series) are served: **Chebyshev polynomials over time segments**.

---

## 3) Chebyshev interpolation: math background

### 3.1 Why Chebyshev?

Polynomial interpolation on equally spaced points can be numerically unstable (Runge phenomenon). Chebyshev nodes minimize the maximum interpolation error for many smooth functions.

For \(n\) nodes, Chebyshev nodes on \([-1,1]\) are:

$$
\xi_k = \cos\left(\frac{\pi(2k+1)}{2n}\right),\quad k = 0,\dots,n-1
$$

We map a time \(t\) in a segment \([t_0, t_1]\) to \(x\in[-1,1]\) by:

$$
\text{mid} = \frac{t_0+t_1}{2},\quad \text{half} = \frac{t_1-t_0}{2},\quad x = \frac{t-\text{mid}}{\text{half}}
$$

Then each coordinate is approximated by a degree-\(d\) Chebyshev series:

$$
X(t) \approx \sum_{j=0}^{d} c_j\,T_j(x)
$$

where \(T_j\) are Chebyshev polynomials of the first kind.

### 3.2 Computing coefficients (discrete cosine transform form)

Given values \(f(\xi_k)\) at the Chebyshev nodes, the coefficients can be computed via:

$$
 c_0 = \frac{1}{n}\sum_{k=0}^{n-1} f(\xi_k),\qquad
 c_j = \frac{2}{n}\sum_{k=0}^{n-1} f(\xi_k)\cos\left(\frac{\pi j(2k+1)}{2n}\right)\; (j\ge 1)
$$

We do this independently for \(x(t)\), \(y(t)\), and \(z(t)\).

### 3.3 Fast evaluation (Clenshaw recurrence)

Evaluating \(\sum c_jT_j(x)\) is done efficiently and stably by Clenshaw recurrence:

- It is \(O(d)\)
- It avoids explicitly constructing high-degree polynomials

Implementation: `clenshaw_eval()` in `../../src/calculus/lunar/moon_cache.rs`.

---

## 4) What we interpolate (and what we do not)

### 4.1 Interpolated: geocentric ecliptic Cartesian position

We cache the output of:

- `Moon::get_geo_position::<Kilometer>(jd)` (ELP2000) from `../../src/calculus/elp2000/elp_series.rs`

But instead of returning a `Position` object, the cache stores three scalar series:

- \(x_{ecl}(t)\), \(y_{ecl}(t)\), \(z_{ecl}(t)\)

This is what ELP2000 computes most expensively.

### 4.2 Not interpolated: topocentric corrections and Earth rotation

Some parts depend strongly on the observer and Earth rotation:
site-vector rotation (UT1 + sidereal time terms), parallax translation, and
local-sidereal-time / hour-angle computations.

Those are computed per query. They are comparatively cheap.

### 4.3 Also interpolated: nutation parameters

IAU 2000B nutation (77 terms) is not huge compared to ELP2000, but it is still
a meaningful cost at the scale of thousands of evaluations. We cache and
linearly interpolate the nutation triplet \((\Delta\psi, \Delta\epsilon, \epsilon_0)\)
to keep those series evaluations out of the per-sample hot path.

The nutation rotation matrix is built from these values:

$$
R = R_1(\epsilon_0 + \Delta\epsilon)\;R_3(\Delta\psi)\;R_1(-\epsilon_0)
$$

See `NutationCache::nutation_rotation()` in `../../src/calculus/lunar/moon_cache.rs`.

---

## 5) Astro chain (cached evaluator reproduces the original physics)

`MoonAltitudeContext::altitude_rad(jd)` reconstructs the same transform pipeline as the original `Moon::get_horizontal()`:

1. **Geocentric ecliptic** \(\vec r_{ecl}(t)\)  
   - Before: full ELP2000 series each call  
   - Now: Chebyshev interpolation

2. **Ecliptic → EquatorialMeanJ2000**  
   Rotation about +X by J2000 mean obliquity \(\epsilon\):

$$
\begin{bmatrix}x_{eq}\\y_{eq}\\z_{eq}\end{bmatrix} =
\begin{bmatrix}
1 & 0 & 0\\
0 & \cos\epsilon & -\sin\epsilon\\
0 & \sin\epsilon & \cos\epsilon
\end{bmatrix}
\begin{bmatrix}x_{ecl}\\y_{ecl}\\z_{ecl}\end{bmatrix}
$$

3. **Topocentric parallax**  
   Compute observer ECEF vector from WGS84, rotate it into the relevant equatorial
   basis using Earth-orientation terms (UT1/EOP), then:

$$
\vec r_{topo} = \vec r_{geo} - \vec r_{site}
$$

4. **Precession (IAU 2006)** to move from J2000 to date-dependent equatorial axes

5. **Nutation (IAU 2000B)** to obtain true-of-date axes  
   - Before: compute nutation series every call  
   - Now: interpolated nutation parameters

6. **Equatorial → Horizontal**  
   Compute HA from LST, then:

$$
\sin(alt)=\sin\delta\sin\phi + \cos\delta\cos\phi\cos(HA)
$$

This keeps the scientific meaning of “Moon altitude” unchanged; we only accelerate how the geocentric ephemeris is obtained.

---

## 6) Interval finding improvement: label crossing direction without probes

The generic interval engine in `../../src/calculus/math_core/intervals.rs` labels
crossings by probing \(f(t\pm\varepsilon)\), which costs **two extra altitude
evaluations per crossing**.

In the Moon case, we can infer direction from the scan bracket:

- If \(g(t_a)=f(t_a)-h_{thr} < 0\) and \(g(t_b) > 0\), then the function entered “above threshold” → direction \(+1\).
- If \(g(t_a) > 0\) and \(g(t_b) < 0\), it exited → direction \(-1\).

The function `find_and_label_crossings()` (in `../../src/calculus/lunar/moon_cache.rs`)
does scan, Brent refinement, and direction labeling in one pass.

---

## 7) Benchmarked comparison (2026-01-01, Roque de los Muchachos)

From `cargo bench --bench moon_altitude "365day"`:

| Method | Time (approx.) |
|---|---:|
| 10-minute scan (validation baseline) | ~17.33 s |
| 2-hour scan, uncached (previous) | ~2.81 s |
| 2-hour scan + Chebyshev cache (new) | **~0.257 s** |

That’s a **~10.9×** speedup vs the previous optimized approach and **~67×** vs the 10-minute scan baseline.

---

## 8) Parameter choices and tuning knobs

Current defaults (see `../../src/calculus/lunar/moon_cache.rs`):

| Parameter | Default |
|---|---|
| Segment length | 4 days |
| Chebyshev degree | 8 (9 nodes) |
| Nutation step | 2 hours |

Shorter segments or higher degree increase precompute work (and usually improve
accuracy). Longer segments or lower degree reduce precompute work but can
degrade interpolation accuracy, which then destabilizes root times.

For horizon crossing detection, you generally want altitude errors well below ~0.01–0.05° to keep root times stable to within ~seconds to tens of seconds.

---

## 9) How this compares to other potential approaches

### A) “Just scan finer”

- Pros: simple, robust
- Cons: scales linearly with evaluations; expensive because ELP2000 dominates

Result here: 10-minute scan is ~17 s (365-day).

### B) Keep 2-hour scan + Brent, but no interpolation

- Pros: already much faster than 10-minute scan
- Cons: still thousands of altitude calls; each call runs ELP2000

Result here: ~2.8 s.

### C) Two-tier ephemeris (cheap model for scan, full ELP2000 for refine)

- Pros: likely good win, simpler than Chebyshev
- Cons: needs careful error control; still calls ELP2000 many times in refine

Chebyshev cache is effectively a high-quality “cheap model” that stays accurate enough for *both* scan and refinement.

### D) Parallelize the year

- Pros: can reduce wall time on multi-core machines
- Cons: adds complexity; makes sense only once single-thread cost is low

With the cache, total time is already ~0.26 s; parallelism is usually unnecessary.

### E) “SDIM” / SIMD

- SIMD already exists in ELP2000 and VSOP87 series evaluation.
- This optimization is different: it reduces the **number of ELP2000 evaluations** dramatically.

Chebyshev caching composes well with SIMD: the precompute stage still benefits from ELP2000’s SIMD summation.

---

## 10) Notes / limitations

Today the cache is built per `find_moon_*` call. If an application issues many
queries over the same span/site, exposing a reusable context would avoid
rebuilding. The cache interpolates *geocentric* position (not topocentric),
which keeps it site-agnostic and preserves correctness; observer-dependent
parallax is applied after interpolation. Nutation is interpolated linearly; on
a 2-hour grid this is conservative for the smooth IAU 2000B corrections.

---

## 11) Where to look in code

- Chebyshev coefficients + evaluation:
  - `compute_cheb_coeffs()` and `clenshaw_eval()` in `../../src/calculus/lunar/moon_cache.rs`
- Position cache builder:
  - `MoonPositionCache::new()` in `../../src/calculus/lunar/moon_cache.rs`
- Nutation interpolation:
  - `NutationCache` in `../../src/calculus/lunar/moon_cache.rs`
- Fast altitude evaluator:
  - `MoonAltitudeContext::altitude_rad()` in `../../src/calculus/lunar/moon_cache.rs`
- Wiring into API:
  - `find_moon_above_horizon()` in `../../src/calculus/lunar/altitude_periods.rs`

---

## 12) Future work ideas

- Cache precession matrices similarly (they vary even more slowly than nutation).
- Add an adaptive degree/segment auto-tuner to guarantee a target error bound.
- Expose a public reusable `MoonAltitudeContext` builder for repeated queries.
