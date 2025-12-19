Below is a **structured technical report** of the issues identified in the current `coordinates` module design.
Issues are **grouped by domain**, each with **context**, **problem statement**, and **why it matters**.
No solutions are prescribed beyond minimal clarification, so you can address them one by one deliberately.

---

# Coordinates Module – Design Issues Report

## Domain C — Reference Frames and Time Modeling

### C1. ICRS and Equatorial treated as identical

**Context**
ICRS ↔ Equatorial is modeled as identity; ecliptic/equatorial conversion uses constant J2000 obliquity.

**Problem**
ICRS, mean equator J2000, and true equator of date are distinct reference systems. Treating them as identical while naming them distinctly is misleading.

**Why it matters**
The type names promise physical rigor that the implementation does not uphold, leading to false confidence in results.

---

### C2. Epoch assumptions hidden inside `From` conversions

**Context**
Some `From`/`Into` conversions implicitly assume `JulianDate::J2000`.

**Problem**
`From` implies a semantics-preserving conversion. Introducing time-dependent behavior silently violates that expectation.

**Why it matters**
This makes it easy to accidentally mix epochs without noticing, producing subtly wrong results.

---

## Domain D — Type-System and API Correctness

### D1. Runtime parameter equality checked only with `debug_assert!`

**Context**
Operations like vector addition rely on `debug_assert!(center_params == other.center_params)`.

**Problem**
In release builds, incompatible coordinates can be combined silently.

**Why it matters**
This is not a recoverable numerical error; it is a logical violation. Silent failure here is dangerous in scientific code.

---

### D2. Parameterized centers weaken compile-time guarantees

**Context**
Centers like `Topocentric` rely on runtime parameters (`ObserverSite`).

**Problem**
The system transitions from compile-time guarantees to runtime checks without clear API boundaries indicating the loss of safety.

**Why it matters**
Users may assume type-level safety still applies when it no longer does.

---

## Domain E — API Structure and Ergonomics

### E1. Cartesian/Spherical × Position/Direction module explosion

**Context**
Parallel modules exist for each combination, leading to repeated names:
`cartesian::position::Ecliptic`, `spherical::direction::Ecliptic`, etc.

**Problem**
This creates a large surface area with duplicated concepts and long import paths.

**Why it matters**
Ergonomic friction encourages incorrect aliasing or misuse and increases maintenance cost.

---

### E2. Frame types duplicated across representations

**Context**
`Ecliptic`, `Equatorial`, etc. appear separately under cartesian and spherical modules.

**Problem**
Frames are representation-independent mathematical objects. Duplicating them ties semantics to layout.

**Why it matters**
This complicates transformations and makes frame identity harder to reason about.

---

## Domain F — Coordinate Conventions

### F1. Horizontal frame axis and azimuth conventions are implicit

**Context**
Horizontal Cartesian axes (`x=North, y=West, z=Zenith`) and azimuth sign conventions are implemented but not encoded in types.

**Problem**
Multiple conventions exist in astronomy and aerospace. The current design relies on documentation and internal knowledge.

**Why it matters**
Users comparing results with external software may get sign-flipped or rotated results without understanding why.

---

## Domain G — Architectural Boundaries

### G1. Geometry, ephemerides, and observation physics are mixed

**Context**
Center transforms combine:

* geometric translation
* ephemeris evaluation
* observational effects

**Problem**
These concerns evolve independently and obey different physical rules, yet are tightly coupled.

**Why it matters**
This reduces extensibility and makes correctness reasoning harder as the system grows.

---

## Summary

The dominant pattern across all domains is **semantic overloading**:

* centers vs observers
* directions vs lines-of-sight
* geometry vs observation physics
* compile-time vs runtime guarantees

None of these are implementation bugs; they are **modeling boundary issues**, which is good news: addressing them systematically will significantly strengthen the design without rewriting core math.

When you are ready, we can:

* pick **one domain** (e.g. Directions semantics or Observation modeling),
* restate its invariants formally,
* and derive the minimal type/API changes needed before moving on to the next.
