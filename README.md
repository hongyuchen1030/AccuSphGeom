# AccuSphGeom

`AccuSphGeom` is a C++ library for **robust spherical geometry on the unit sphere**.

The library is organized around two fundamental problem classes:

1. **Construction problems**
   Computing geometric quantities such as intersection points.

2. **Predicate problems**
   Robust classification and decision logic, such as point-in-polygon and arc membership.

---

# 1. Installation & Quick Start

## Installation

```bash
git clone https://github.com/hongyuchen1030/AccuSphGeom.git
```

Compile with:

```bash
g++ program.cpp -I/path/to/AccuSphGeom/include -std=c++17
```

Include the library:

```cpp
#include <accusphgeom/accusphgeom.hpp>
```

---

# 2. Construction Problems (Intersection Computation)

The construction layer provides **accurate and vectorizable intersection point computation on the sphere**.

This work is based on:

* Chen, H. *Accurate and Robust Great Circle Arc Intersection and Great Circle Arc Constant Latitude Intersection on the Sphere.*
  SIAM Journal on Scientific Computing
  https://epubs.siam.org/doi/full/10.1137/25M1737614

---

## 2.1 Supported Intersection Types

The library provides two intersection types:

### Great-circle arc × Great-circle arc (GCA–GCA)

```cpp
#include <accusphgeom/constructions/gca_gca_intersection.hpp>
```

### Great-circle arc × Constant-latitude line (GCA–ConstLat)

```cpp
#include <accusphgeom/constructions/gca_constlat_intersection.hpp>
```

---

## 2.2 API Design

Each intersection routine provides two interfaces.

### Scalar API (convenience interface)

```cpp
auto p = accusphgeom::constructions::gca_gca_intersection(a0, a1, b0, b1);
```

* Returns the unique intersection point.
* Throws `std::domain_error` if:

  * no valid intersection exists, or
  * both antipodal candidates satisfy the minor-arc constraints.

---

### Vectorized / Batch API (performance interface)

```cpp
auto r = accusphgeom::constructions::try_gca_gca_intersection(a0, a1, b0, b1);
```

Returns:

* `r.point` — mask-selected intersection point
* `r.status` — classification flag

Example:

```cpp
if (r.status == 0) {
    use(r.point);
} else {
    // handle ambiguity or failure
}
```

This interface is designed for:

* SIMD vectorization
* batched execution
* branch-free construction kernels

---

### GCA–ConstLat follows the same pattern

```cpp
auto r = accusphgeom::constructions::try_gca_constlat_intersection(a0, a1, z0);
```

---

## 2.3 Design Principle

* **Scalar API**: convenience interface with exception-based handling
* **try_ API**: performance interface returning explicit status for downstream control

---

# 3. Predicate Problems

The predicate layer provides **robust geometric classification and decision logic**.

This work is based on:

* Chen, H. (2026). *Accurate and Robust Algorithms for Spherical Polygon Operations.*
  EGUsphere preprint
  https://egusphere.copernicus.org/preprints/2026/egusphere-2026-636/
  PDF: https://egusphere.copernicus.org/preprints/2026/egusphere-2026-636/egusphere-2026-636.pdf

---

## 3.1 Minor Arc Predicate

```cpp
#include <accusphgeom/predicates/on_minor_arc.hpp>

bool on_arc = accusphgeom::predicates::on_minor_arc(q, a, b);
```

---

## 3.2 Spherical Point-in-Polygon (SPIP)

```cpp
#include <accusphgeom/algorithms/point_in_polygon_sphere.hpp>

auto loc = accusphgeom::algorithms::point_in_polygon_sphere(q, poly);
```

Returns:

* `Inside`
* `Outside`
* `OnVertex`
* `OnEdge`

---

## 3.3 Core Algorithm Notes for Spherical Point in Polygon (SPIP)

The ray endpoint `R` is normally chosen as a perturbed antipode of `q`, so the
ray is geometrically well separated from the query.

For each polygon edge `AB`, the crossing logic uses four orientation signs:

* `s_qR_A = orient(q, R, A)`
* `s_qR_B = orient(q, R, B)`
* `s_AB_q = orient(A, B, q)`
* `s_AB_R = orient(A, B, R)`

In the nondegenerate case, the implementation uses the strict 4-sign crossing
theorem: the arcs `qR` and `AB` cross if and only if the endpoint orientations
are strictly separated on both supporting great circles.

Boundary handling is performed before parity counting:

* if `q` exactly matches a polygon vertex, return `OnVertex`
* if `q` lies on a polygon edge minor arc, return `OnEdge`

The ray endpoint `R` is constructed as a perturbed antipode of `q`:

* start from `-q`
* perturb the least dominant coordinate
* renormalize to the sphere

If any polygon edge yields `s_AB_R == 0`, the algorithm throws an error. This
typically indicates a polygon that is too large, close to hemispherical, or
otherwise poorly separated from the ray construction.

---

## 3.4 Robustness Tiers (With Global IDs)

When global IDs are available, **Simulation of Simplicity (SoS)** resolves
ray-vertex degeneracies.

| Tier       | Name           | Description                                                      |
| :--------- | :------------- | :--------------------------------------------------------------- |
| **Tier 1** | Full Global    | User provides `q`, `R`, and vertex IDs. Strongest mode.          |
| **Tier 2** | Semi-Specified | User provides `q` and vertex IDs; library infers `R` and its ID. |
| **Tier 3** | Local/Internal | User provides vertex IDs; library infers IDs for `q` and `R`.    |

---

# 4. Architecture Overview

The library is organized into three layers:

1. **Predicates**
   Robust classification and sign evaluation

2. **Constructions**
   High-accuracy geometric computations

3. **Algorithms**
   Built on top of predicates and constructions

---

# 5. Usage Guidance for Performance-Critical Workflows

For high-performance applications:

* Prefer `try_...` APIs
* Batch inputs where possible
* Use returned `status` values as masks
* Avoid scalar APIs in performance-critical loops

---

# 6. Third-Party Foundation

`AccuSphGeom` is implemented on top of bundled third-party multi-precision
support in `third_party/geogram_psm`, and Jonathan Shewchuk's adaptive
predicates for core orientation-sign evaluation, including the vendored
`predicates.c`.

* The vendored `predicates.c` is public domain
* `third_party/predicates.c` remains part of the repository's third-party foundation
* These components provide the exact and adaptive arithmetic used by the predicate layer

---

# 7. License

See `LICENSES/` and `REUSE.toml`
