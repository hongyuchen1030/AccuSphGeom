# AccuSphGeom

`AccuSphGeom` is a C++ library for robust spherical geometry on the unit sphere.
The current public algorithm layer provides **spherical point-in-polygon (PIP)**
classification together with the supporting robust predicates and low-level EFT
numeric building blocks used by the construction layer.

The library classifies a query point as:

- `Outside`
- `Inside`
- `OnVertex`
- `OnEdge`

### Key Properties

- **Robustness:** Handles common geometric degeneracies and numerical instabilities.
- **Tiered Reliability:** Supports three robustness tiers when global IDs are available.
- **Adaptive PIP:** Spherical point-in-polygon uses adaptive predicates with exact fallback.
- **Layered Architecture:** Predicate programs and construction programs are parallel first-class layers; algorithms sit above them.

This repository contains the implementation. The associated paper discusses the
algorithms in more detail and provides broader algorithmic context:

- Chen, H. (2026). Accurate and Robust Algorithms for Spherical Polygon
  Operations. EGUsphere preprint.
  https://egusphere.copernicus.org/preprints/2026/egusphere-2026-636/

## 1. Installation & Quick Start

### Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/hongyuchen1030/AccuSphGeom.git
   ```

2. Include the header:
   ```cpp
   #include <accusphgeom/algorithms/point_in_polygon_sphere.hpp>
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

## 2. Core Algorithm (Spherical PIP)

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

## 3. Robustness Tiers (With Global IDs)

When global IDs are available, **Simulation of Simplicity (SoS)** resolves
ray-vertex degeneracies.

| Tier | Name | Description |
| :--- | :--- | :--- |
| **Tier 1** | Full Global | User provides `q`, `R`, and vertex IDs. Strongest mode. |
| **Tier 2** | Semi-Specified | User provides `q` and vertex IDs; library infers `R` and its ID. |
| **Tier 3** | Local/Internal | User provides vertex IDs; library infers IDs for `q` and `R`. |

## 4. Package Architecture

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
    orient3d.hpp
    quadruple3d.hpp

  constructions/
    accucross.hpp
    gca_constlat_intersection.hpp
    gca_gca_intersection.hpp

  algorithms/
    point_in_polygon_sphere.hpp
```

## 5. License

This repository is REUSE-compliant. See `REUSE.toml` and `LICENSES/` for the
package metadata and license texts.
