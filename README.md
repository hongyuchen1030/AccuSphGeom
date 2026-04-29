

# AccuSphGeom

`AccuSphGeom` is a C++ library for **robust spherical geometry on the unit sphere**.

The algorithms implemented in this repository are described in:

Chen, H. (2026)
*Accurate and Robust Algorithms for Spherical Polygon Operations*
[https://egusphere.copernicus.org/preprints/2026/egusphere-2026-636/](https://egusphere.copernicus.org/preprints/2026/egusphere-2026-636/)

This repository provides a **reference implementation of all related algorithms** in that work.
If you use the algorithm APIs (e.g., spherical point-in-polygon or latitude–longitude bounds), please cite the above paper.

---

In addition, the construction layer is based on:

Chen, H.
*Accurate and Robust Great Circle Arc Intersection and Great Circle Arc Constant Latitude Intersection on the Sphere.*
SIAM Journal on Scientific Computing
[https://epubs.siam.org/doi/full/10.1137/25M1737614](https://epubs.siam.org/doi/full/10.1137/25M1737614)

If you use the **construction APIs** (e.g., GCA–GCA or GCA–ConstLat intersection), please also cite this paper.

---

The library is organized in **two conceptual levels**:

## Level 1: Fundamental Building Blocks

* **Construction programs**
  Accurate computation of geometric quantities (e.g., intersection points)

* **Predicate programs**
  Robust classification and decision logic (e.g., orientation, arc membership)

## Level 2: Algorithms

* **Spherical Point-in-Polygon (SPIP)**
* **Latitude–Longitude Bounds (LatLon Bounds)**


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

# 2. Algorithms

This section introduces the **algorithm-level APIs**, which combine predicates and constructions.

---

## 2.1 Latitude–Longitude Bounds

The lat-lon bounds algorithm computes the bounding box of a spherical face.

### API Overview

* `get_face_location_info(...)`

  ```cpp
  template <typename T, std::size_t N>
  FaceLocationInfo<T> get_face_location_info(
      const std::array<numeric::Vec3<T>, N>& face_vertices,
      T polar_cap_lat_deg = static_cast<T>(default_polar_cap_lat_deg));
  ```

  Classifies one face using the user-provided `polar_cap_lat_deg`.
  This cap separates local faces from north/south pole-candidate faces.
  Returns both the classification label and latitude extrema (`face_z_max`, `face_z_min`).

* `generate_lat_lon_bounds_local(...)`

  ```cpp
  template <typename T, std::size_t N>
  LatLonBounds<T> generate_lat_lon_bounds_local(
      const std::array<numeric::Vec3<T>, N>& face_vertices,
      const FaceLocationInfo<T>& info,
      T endpoint_lat_snap_tol_deg =
          static_cast<T>(default_endpoint_lat_snap_tol_deg));
  ```

  Handles non-pole faces.
  This is the dominant fast path and is fully SIMD-friendly.

* `generate_lat_lon_bounds_pole(...)`

  ```cpp
  template <std::size_t N>
  LatLonBounds<double> generate_lat_lon_bounds_pole(
      const std::array<numeric::Vec3<double>, N>& face_vertices,
      const FaceLocationInfo<double>& info,
      const std::int64_t* global_vertex_ids = nullptr,
      const std::int64_t* pole_query_id = nullptr,
      const std::int64_t* waypoint_id = nullptr,
      double endpoint_lat_snap_tol_deg = default_endpoint_lat_snap_tol_deg);
  ```

  Used only for pole-candidate faces.
  It uses robust point-in-polygon logic and is not fully vectorized.

* `generate_lat_lon_bounds(...)`

  ```cpp
  template <std::size_t N>
  LatLonBounds<double> generate_lat_lon_bounds(
      const std::array<numeric::Vec3<double>, N>& face_vertices,
      double polar_cap_lat_deg = default_polar_cap_lat_deg,
      const std::int64_t* global_vertex_ids = nullptr,
      const std::int64_t* pole_query_id = nullptr,
      const std::int64_t* waypoint_id = nullptr,
      double endpoint_lat_snap_tol_deg = default_endpoint_lat_snap_tol_deg);
  ```

  Convenience dispatcher.
  It calls `get_face_location_info(...)`, then routes to the local or pole path.
  For performance-critical workflows, manually classify and dispatch.

### Batch Dispatch Example

```cpp
#include <accusphgeom/algorithms/lat_lon_bounds.hpp>

#include <array>
#include <vector>

using accusphgeom::algorithms::FaceLocationLabel;
using accusphgeom::algorithms::LatLonBounds;
using accusphgeom::algorithms::generate_lat_lon_bounds_local;
using accusphgeom::algorithms::generate_lat_lon_bounds_pole;
using accusphgeom::algorithms::get_face_location_info;
using accusphgeom::numeric::Vec3;

std::vector<LatLonBounds<double>> compute_bounds(
    const std::vector<std::array<Vec3<double>, 4>>& faces) {
  constexpr double polar_cap_lat_deg = 80.0;

  std::vector<LatLonBounds<double>> bounds;
  bounds.reserve(faces.size());

  for (const auto& face : faces) {
    const auto info = get_face_location_info(face, polar_cap_lat_deg);

    if (info.label == FaceLocationLabel::Local) {
      bounds.push_back(generate_lat_lon_bounds_local(face, info));
    } else {
      bounds.push_back(generate_lat_lon_bounds_pole(face, info));
    }
  }

  return bounds;
}
```

Recommended high-performance workflow:

1. Call `get_face_location_info(...)`.
2. Route local faces to `generate_lat_lon_bounds_local(...)`.
3. Route pole-candidate faces to `generate_lat_lon_bounds_pole(...)`.

---

### Design Notes

* The **local path dominates performance** and is designed for SIMD/vectorized execution
* The **pole path is rare** and focuses on robustness rather than SIMD efficiency
* Only the pole path incurs non-vectorizable overhead due to robust predicate usage

---

## 2.2 Spherical Point-in-Polygon (SPIP)

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

### Core Algorithm Notes

The ray endpoint `R` is normally chosen as a perturbed antipode of `q`.

For each polygon edge `AB`, the crossing logic uses four orientation signs:

* `s_qR_A = orient(q, R, A)`
* `s_qR_B = orient(q, R, B)`
* `s_AB_q = orient(A, B, q)`
* `s_AB_R = orient(A, B, R)`

Crossing occurs if and only if the endpoints are strictly separated on both great circles.

Boundary handling:

* Exact vertex match → `OnVertex`
* On minor arc → `OnEdge`

Ray construction:

* Start from `-q`
* Perturb the least dominant coordinate
* Renormalize

Degenerate cases (`s_AB_R == 0`) trigger an error.

---

### Robustness Tiers (With Global IDs)

Simulation of Simplicity (SoS) resolves degeneracies.

| Tier           | Description                            |
| -------------- | -------------------------------------- |
| Full Global    | User provides `q`, `R`, and vertex IDs |
| Semi-Specified | User provides `q` and vertex IDs       |
| Local/Internal | Library assigns IDs                    |

---

# 3. Construction Programs

The construction layer provides **accurate and vectorizable intersection computations**.

Based on:

Chen, H.
*Accurate and Robust Great Circle Arc Intersection and Great Circle Arc Constant Latitude Intersection on the Sphere.*
SIAM Journal on Scientific Computing
[https://epubs.siam.org/doi/full/10.1137/25M1737614](https://epubs.siam.org/doi/full/10.1137/25M1737614)

---

## Supported Intersection Types

### GCA–GCA

```cpp
#include <accusphgeom/constructions/gca_gca_intersection.hpp>
```

### GCA–ConstLat

```cpp
#include <accusphgeom/constructions/gca_constlat_intersection.hpp>
```

---

## API Design

### Scalar API

```cpp
auto p = accusphgeom::constructions::gca_gca_intersection(a0, a1, b0, b1);
```

* Returns intersection point
* Throws on invalid or ambiguous cases

---

### `try_` API for Packed / SIMD-Style Execution

```cpp
auto r = accusphgeom::constructions::try_gca_gca_intersection(a0, a1, b0, b1);
```

Returns:

* `r.point`
* `r.status`

Designed for:

* packed / SIMD-style execution
* batched workflows
* branch-free kernels

The `try_...` construction APIs are the vectorizable interface.
They support scalar `double`, and they also support packed numeric backends that
behave like lane-wise arithmetic types.

Current packed support includes Eigen arrays through the Eigen adapter header:

```cpp
#include <accusphgeom/adapters/eigen/numeric.hpp>
```

This enables packed usage such as:

```cpp
using Pack = EigenPack<4>;
auto r = accusphgeom::constructions::try_gca_gca_intersection(a0, a1, b0, b1);
```

See the working examples in:

* `tests/test_gca_constlat_eigen_pack.cpp`
* `tests/test_gca_gca_eigen_pack.cpp`

These examples show how to:

* define packed `EigenPack<N>` coordinates
* build `numeric::Vec3<Pack>` inputs
* call the `try_...` APIs for batched lane-wise evaluation
* compare packed outputs against scalar reference results

---

### GCA–ConstLat

```cpp
auto r = accusphgeom::constructions::try_gca_constlat_intersection(a0, a1, z0);
```

This follows the same packed / SIMD-style model and also works with the Eigen
adapter shown above.

---

## Design Principle

* Scalar API: convenience
* `try_` API: performance, vectorization, and control

Today, Eigen is the supported packed backend for the construction `try_...`
APIs. We plan to add `std::simd` support in the future using the same packed
interface style.

---

# 4. Predicate Programs

The predicate layer provides **robust geometric classification and decision logic**.

---

## Implemented Predicates

### 1. Minor Arc Membership

```cpp
#include <accusphgeom/predicates/on_minor_arc.hpp>

bool on_arc = accusphgeom::predicates::on_minor_arc(q, a, b);
```

---

### 2. Orientation in 3D

```cpp
orient(a, b, c)
```

---

### 3. Quadruple Product Predicate

```cpp
quadruple(a, b, c, d)
```

---

These predicates provide the **robust foundation** for all higher-level algorithms.

---

# 5. Usage Guidance

For performance-critical workflows:

* Prefer `try_...` APIs for intersection calculations
* For lat–lon bounds, first call `get_face_location_info(...)` (vectorized classification), then:

  * use `generate_lat_lon_bounds_local(...)` for the dominant non-pole faces (SIMD-friendly)
  * use `generate_lat_lon_bounds_pole(...)` only for pole candidates (robust, not fully vectorized)
* Batch inputs where possible
* Use status flags as masks
* Avoid scalar APIs in tight loops

For GPU workflows:

* Perform classification first
* Separate local and pole workloads
* Avoid mixed branching within kernels


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
