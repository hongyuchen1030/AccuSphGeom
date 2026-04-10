# AccuSphGeom

`AccuSphGeom` is a C++ library for robust spherical geometry on the unit sphere.
The current public API is organized into two categories:

1. Predicate problems: spherical point-in-polygon (SPIP) together with
   supporting predicate helpers such as `on_minor_arc` and `quadruple3D`.
2. Construction problems: intersection-point computations for great-circle arc
   with great-circle arc, and great-circle arc with constant latitude. For
   convenience, the library also provides an accurate cross-product API,
   `accucross`.

## SPIP, SPherical Point in Polygon Polygon
This repository contains the implementation. The main SPIP method is described
in:

- Chen, H. (2026). Accurate and Robust Algorithms for Spherical Polygon
  Operations. EGUsphere preprint.
  https://egusphere.copernicus.org/preprints/2026/egusphere-2026-636/
  PDF: https://egusphere.copernicus.org/preprints/2026/egusphere-2026-636/egusphere-2026-636.pdf

Which has following properties 
- **Robustness:** Handles common geometric degeneracies and numerical instabilities.
- **Tiered Reliability:** Supports three robustness tiers when global IDs are available.
- **Adaptive SPIP:** Spherical point-in-polygon uses adaptive predicates with exact fallback.
- **Layered Architecture:** Predicate programs and construction programs are parallel first-class layers; algorithms sit above them.

The library classifies a query point as:

- `Outside`
- `Inside`
- `OnVertex`
- `OnEdge`

The implementation also builds on top of the bundled third-party 

The robust pipeline uses Jonathan Shewchuk's adaptive predicates for core
orientation-sign evaluation, including the vendored `predicates.c`, which is in public domain


### Simulation of Simplicity (SoS)

Simulation of Simplicity is used for global-ID degeneracy resolution in the
Tier 1, Tier 2, and Tier 3 global-ID paths.
  Transactions on Graphics, 9(1), 66-104.
  https://doi.org/10.1145/77635.77639

### Half-Open Rule

The non-global-ID branch uses a deterministic half-open convention for
ray-vertex tie-breaking.
 Important note:
- the rule is used as deterministic tie-breaking in non-global-ID mode
- for full robustness, users should use Tier 1 global-ID mode

### Geogram Multiprecision

Geogram multiprecision components are used as the exact fallback for predicate
evaluation in the robust pipeline.
- Geogram PCK-related components included in the vendored code layout used by
  this project

## Intersection Point Constructions
This is based on 
- Chen, H. Accurate and Robust Great Circle Arc Intersection and Great
  Circle Arc Constant Latitude Intersection on the Sphere. SIAM Journal on
  Scientific Computing.
  https://epubs.siam.org/doi/full/10.1137/25M1737614

We provide following API:
- accurate cross product: `include/accusphgeom/constructions/accucross.hpp`
- Intersection between GCA and ConstLat: `include/accusphgeom/constructions/gca_constlat_intersection.hpp`
- Intersection between GCA and GCA:`include/accusphgeom/constructions/gca_gca_intersection.hpp`
It can provide fast and accurate intersection point calculation: it achieves as calculating in doubling the working prevision and round it back to the working precision, and without performance overhead when vectorized and parallized properly.

The `AccuCross`, in `include/accusphgeom/constructions/accucross.hpp`; The internal calculation function `AccuX_GCA` in `include/accusphgeom/constructions/gca_gca_intersection.hpp` and smilarly the `AccuX_constlat` in `include/accusphgeom/constructions/gca_constlat_intersection.hpp` are already SIMD vectorized and support batch processing, further opimization will be implemented for the final final API call for the intersection calculation.
## 1. Installation & Quick Start

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/hongyuchen1030/AccuSphGeom.git
   ```

2. Include the umbrella header or the specific API header you need:
   ```cpp
   #include <accusphgeom/accusphgeom.hpp>
   ```

3. Add the `include/` directory to your include path (`C++17` required):
   ```bash
   g++ program.cpp -I/path/to/AccuSphGeom/include -std=c++17
   ```

### Quick Start Examples

```cpp
#include <array>
#include <accusphgeom/algorithms/point_in_polygon_sphere.hpp>

const std::array<double, 3> q = {
    0.5773502691896257, 0.5773502691896257, 0.5773502691896257};
const std::array<std::array<double, 3>, 3> poly = {{
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
}};

const auto loc = accusphgeom::algorithms::point_in_polygon_sphere(q, poly);
```

```cpp
#include <array>
#include <cstdint>
#include <accusphgeom/algorithms/point_in_polygon_sphere.hpp>

const std::array<double, 3> q = {
    0.5773502691896257, 0.5773502691896257, 0.5773502691896257};
const std::array<double, 3> r = {
    -0.5773502691896257, -0.5773502691896257, -0.5773502691896257};
const std::array<std::array<double, 3>, 3> poly = {{
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
}};
const std::array<std::int64_t, 3> vertex_ids = {10, 20, 30};

const auto loc = accusphgeom::algorithms::point_in_polygon_sphere(
    q, 40, r, 50, poly, vertex_ids);
```

This library provides API coverage across robustness modes and container
interfaces. In total, there are **12 API combinations**, formed by:

- **4 robustness modes**:
  - no global ID
  - Tier 1 (full global robustness)
  - Tier 2 (semi-specified global robustness)
  - Tier 3 (local/internal robustness)
- **3 API overload families**:
  - raw pointer
  - `std::array`
  - `std::vector`

Total combinations:
`4 (robustness) x 3 (overloads) = 12 API cases`

All of these combinations are explicitly tested in
`tests/test_library_usage.cpp`.

## 2. Using The Public API

### SPIP algorithm

Header:

```cpp
#include <accusphgeom/algorithms/point_in_polygon_sphere.hpp>
```

Primary call pattern:

```cpp
const auto loc = accusphgeom::algorithms::point_in_polygon_sphere(q, poly);
```

Use this when you want to classify a query point on the unit sphere against a
spherical polygon. The result is one of `Outside`, `Inside`, `OnVertex`, or
`OnEdge`.

If you use the SPIP algorithm in research software or publications, cite the
EGUsphere preprint above as the algorithmic source.

### `on_minor_arc` predicate API

Header:

```cpp
#include <accusphgeom/predicates/on_minor_arc.hpp>
```

Supported overload families:

- raw pointer: `on_minor_arc(q_ptr, a_ptr, b_ptr, tolerance)`
- `std::array`: `on_minor_arc(q, a, b, tolerance)`
- `std::vector`: `on_minor_arc(q, a, b, tolerance)`

Example:

```cpp
#include <array>
#include <accusphgeom/predicates/on_minor_arc.hpp>

const std::array<double, 3> q = {0.0, 1.0, 0.0};
const std::array<double, 3> a = {1.0, 0.0, 0.0};
const std::array<double, 3> b = {0.0, 1.0, 0.0};

const bool on_arc = accusphgeom::predicates::on_minor_arc(q, a, b);
```

This predicate answers whether `q` lies on the minor arc from `a` to `b`.
Passing `tolerance = 0` uses the exact predicate path for coplanarity testing.
Passing a nonzero tolerance switches the coplanarity gate to a floating-point
filter before the endpoint-order test.

If you use this predicate API as part of SPIP or related robust spherical
predicate work, cite the EGUsphere preprint above.

### Construction API: GCA-GCA intersection

Header:

```cpp
#include <accusphgeom/constructions/gca_gca_intersection.hpp>
```

Primary call pattern:

```cpp
#include <array>
#include <accusphgeom/constructions/gca_gca_intersection.hpp>

const std::array<double, 3> a0 = {1.0, 0.0, 0.0};
const std::array<double, 3> a1 = {0.0, 1.0, 0.0};
const std::array<double, 3> b0 = {0.0, 0.0, 1.0};
const std::array<double, 3> b1 = {0.0, 1.0, 0.0};

const auto p =
    accusphgeom::constructions::gca_gca_intersection(a0, a1, b0, b1);
```

Use `gca_gca_intersection` when you want the unique intersection point that
lies on both minor arcs. If you need both antipodal great-circle candidates,
call `accux_gca(...)` and inspect `point_pos` and `point_neg`.

This API throws `std::domain_error` when the supporting great circles are
degenerate or coincident, when both antipodal candidates lie on both minor
arcs, or when neither candidate lies on both minor arcs.

If you use this construction API, cite the SIAM paper above.

### Construction API: `AccuCross`

Header:

```cpp
#include <accusphgeom/constructions/accucross.hpp>
```

Primary call pattern:

```cpp
#include <array>
#include <accusphgeom/constructions/accucross.hpp>

const std::array<double, 3> a = {1.0, 0.0, 0.0};
const std::array<double, 3> b = {0.0, 1.0, 0.0};

const auto c = accusphgeom::constructions::accucross(a, b);
```

For convenience, the package also provides an accurate cross product
construction API, `AccuCross`, in
`include/accusphgeom/constructions/accucross.hpp`. The return value stores the
cross product as a high/low expansion pair through `c.hi` and `c.lo`.

If you already have split high/low inputs, the API also provides:

```cpp
const auto c =
    accusphgeom::constructions::accucross(a_hi, a_lo, b_hi, b_lo);
```

`AccuCross` belongs to the construction problem. If you use it in research
software or publications, cite the SIAM paper above.

### Construction API: GCA-constLat intersection

Header:

```cpp
#include <accusphgeom/constructions/gca_constlat_intersection.hpp>
```

Primary call pattern:

```cpp
#include <array>
#include <accusphgeom/constructions/gca_constlat_intersection.hpp>

const std::array<double, 3> a0 = {1.0, 0.0, 0.0};
const std::array<double, 3> a1 = {0.0, 1.0, 0.0};
const double z0 = 0.5;

const auto p =
    accusphgeom::constructions::gca_constlat_intersection(a0, a1, z0);
```

Use `gca_constlat_intersection` when you want the unique point on the minor
arc from `a0` to `a1` whose `z` coordinate equals `z0`. If you need both
constant-latitude candidates together with the high/low split normal, call
`accux_constlat(...)`.

This API throws `std::domain_error` when the great-circle normal is degenerate,
when both antipodal candidates lie on the same minor arc, or when neither
candidate lies on the minor arc.

If you use this construction API, cite the SIAM paper above.

## 3. Core Algorithm Notes for Spherical Point in Polygon (SPIP)

The ray endpoint `R` is normally chosen as a perturbed antipode of `q`, so the
ray is geometrically well separated from the query.

For each polygon edge `AB`, the crossing logic uses four orientation signs:

- `s_qR_A = orient(q, R, A)`
- `s_qR_B = orient(q, R, B)`
- `s_AB_q = orient(A, B, q)`
- `s_AB_R = orient(A, B, R)`

In the nondegenerate case, the implementation uses the strict 4-sign crossing
theorem: the arcs `qR` and `AB` cross if and only if the endpoint orientations
are strictly separated on both supporting great circles.

Boundary handling is performed before parity counting:

- if `q` exactly matches a polygon vertex, return `OnVertex`
- if `q` lies on a polygon edge minor arc, return `OnEdge`

The ray endpoint `R` is constructed as a perturbed antipode of `q`:

- start from `-q`
- perturb the least dominant coordinate
- renormalize to the sphere

If any polygon edge yields `s_AB_R == 0`, the algorithm throws an error. This
typically indicates a polygon that is too large, close to hemispherical, or
otherwise poorly separated from the ray construction.

## 4. Robustness Tiers (With Global IDs)

When global IDs are available, **Simulation of Simplicity (SoS)** resolves
ray-vertex degeneracies.

| Tier | Name | Description |
| :--- | :--- | :--- |
| **Tier 1** | Full Global | User provides `q`, `R`, and vertex IDs. Strongest mode. |
| **Tier 2** | Semi-Specified | User provides `q` and vertex IDs; library infers `R` and its ID. |
| **Tier 3** | Local/Internal | User provides vertex IDs; library infers IDs for `q` and `R`. |

## 5. Package Architecture

The package is organized into three layers:

1. **Predicate programs**
   Robust sign, incidence, and classification routines based on adaptive
   precision. These do not perform geometric constructions.
2. **Construction programs**
   High-accuracy geometric quantity computations built on reusable EFT numeric
   support. The low-level EFT arithmetic building blocks remain part of the
   package.
3. **Algorithms**
   Higher-level geometry algorithms built on predicates and/or constructions.
   Spherical point-in-polygon currently uses the adaptive predicate layer.

Current public include layout:

```text
include/accusphgeom/
  accusphgeom.hpp

  numeric/
    simd_fma.hpp
    eft.hpp

  predicates/
    on_minor_arc.hpp
    orient3d.hpp
    quadruple3d.hpp

  constructions/
    accucross.hpp
    gca_constlat_intersection.hpp
    gca_gca_intersection.hpp

  algorithms/
    point_in_polygon_sphere.hpp
```

## 6. Third-Party Foundation

`AccuSphGeom` is implemented on top of bundled third-party multi-precision
support in `third_party/geogram_psm`. And Jonathan Shewchuk's adaptive predicates for core
orientation-sign evaluation, including the vendored `predicates.c`
- the vendored `predicates.c` is public domain
- `third_party/predicates.c` remains part of the repository's third-party
  foundation
package foundation for the exact and adaptive arithmetic used by the predicate
layer. Keep the repository license metadata and third-party notices intact.

## 7. License

This repository is REUSE-compliant. See `REUSE.toml` and `LICENSES/` for the
package metadata and license texts.
