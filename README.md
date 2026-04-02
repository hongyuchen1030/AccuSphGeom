# AccuSphGeom

`AccuSphGeom` is a C++ geometry package for robust spherical predicates, EFT-based constructions, and algorithms built on top of those two parallel layers.

## Architecture

The package is organized as:

- `numeric/`: reusable EFT and FMA kernels
- `predicates/`: adaptive sign and incidence programs
- `constructions/`: accurate geometric constructions using EFT
- `algorithms/`: higher-level routines built from predicates and constructions

The public umbrella header is:

```cpp
#include <accusphgeom/accusphgeom.hpp>
```

## Spherical Point-In-Polygon

`accusphgeom::algorithms::point_in_polygon_sphere()` uses only the adaptive predicate path. The old EFT-based PIP route has been removed.

The classifier returns:

- `Location::Outside`
- `Location::Inside`
- `Location::OnVertex`
- `Location::OnEdge`

The API supports:

- raw pointers
- `std::array`
- `std::vector`
- no global IDs
- Tier 1, Tier 2, and Tier 3 symbolic robustness modes

Example:

```cpp
#include <array>
#include <accusphgeom/algorithms/point_in_polygon_sphere.hpp>

const std::array<double, 3> q = {
    0.5773502691896257, 0.5773502691896257, 0.5773502691896257};
const std::array<std::array<double, 3>, 3> polygon = {{
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
}};

const auto loc = accusphgeom::algorithms::point_in_polygon_sphere(q, polygon);
```

## Construction Modules

The construction layer currently exposes:

- `accusphgeom::constructions::accucross`
- `accusphgeom::constructions::gca_constlat_intersection`
- `accusphgeom::constructions::gca_gca_intersection`

These are implemented in terms of reusable kernels in `accusphgeom/numeric/eft.hpp`, including:

- `two_prod_fma`
- `two_sum`
- `accurate_difference_of_products`
- `compensated_dot_product`
- `sum_of_squares_c`
- `acc_sqrt_re`
* maintains robustness near degeneracies

### Design Summary

* Robust pipeline: correctness-first, always safe
* EFT pipeline: performance-oriented, with controlled fallback

Both pipelines:

* share identical geometric logic
* differ only in numerical evaluation strategy

References for the EFT design are given in Chen et al. (2025) and Ogita et al. (2005) later in the References and Acknowledgements section


-----

## 5\. Non-Global-ID Strategy

When IDs are absent, the library falls back to a **deterministic local tie-breaking rule** (Half-Open Convention) rather than symbolic perturbation.

  * **Crossing Rule:** If $s_{qR,A}$ or $s_{qR,B}$ is zero, a crossing is counted only if `(below_a XOR below_b)`.
  * **Important:** This provides deterministic behavior but is not a substitute for full global symbolic robustness. Use Tier 1 if global SoS is required.

-----

## 6\. Internal ID Assignment

If the library must assign IDs internally (Tiers 2/3):

1.  It calculates the maximum vertex ID.
2.  It assigns `q_id` and `r_id` as `max_id + 1` and `max_id + 2` respectively.
3.  If IDs would overflow, it uses a **MEX** (Minimum Excluded value) over nonnegative integers.

-----

## 7\. Tests and Visualization

Detailed tests are located in [`tests/`](tests/).

  * **Unit Tests:** Covers both pipelines and all robustness tiers.
  * **Visualization:** `tests/test_pip_complicated_visualization.nb` (Mathematica) visualizes the polygon, query points, and great-circle arcs for complex cases.

-----

## 8. References and Acknowledgements

### 8.1 Associated Paper

The algorithms implemented in this repository are discussed in detail in:

- Chen, H. (2026). Accurate and Robust Algorithms for Spherical Polygon
  Operations. EGUsphere preprint.
  https://egusphere.copernicus.org/preprints/2026/egusphere-2026-636/

This repository contains the implementation. The paper provides the full
algorithmic context and derivations.

### 8.2 Robust Geometric Predicates (Shewchuk)

The robust pipeline uses Jonathan Shewchuk's adaptive predicates for core
orientation-sign evaluation, including the vendored `predicates.c`
implementation.

- Shewchuk, J. R. (1997). Adaptive Precision Floating-Point Arithmetic and Fast
  Robust Geometric Predicates. Discrete & Computational Geometry, 18(3),
  305-363.

Notes:

- the vendored `predicates.c` is public domain
- see [third_party/LICENSE_shewchuk.txt](third_party/LICENSE_shewchuk.txt)

### 8.3 Simulation of Simplicity (SoS)

Simulation of Simplicity is used for global-ID degeneracy resolution in the
Tier 1, Tier 2, and Tier 3 global-ID paths.

- Edelsbrunner, H., & Mucke, E. P. (1990). Simulation of Simplicity: A
  Technique to Cope with Degenerate Cases in Geometric Algorithms. ACM
  Transactions on Graphics, 9(1), 66-104.
  https://doi.org/10.1145/77635.77639

### 8.4 Half-Open Rule

The non-global-ID branch uses a deterministic half-open convention for
ray-vertex tie-breaking.

- Hormann, K., & Agathos, A. (2001). The Point in Polygon Problem for
  Arbitrary Polygons. Computational Geometry, 20(3), 131-144.
  https://doi.org/10.1016/S0925-7721(01)00012-8

Important note:

- the original work is planar
- this repository applies the rule on the sphere
- no global theoretical proof is claimed here for the spherical case
- the rule is used as deterministic tie-breaking in non-global-ID mode
- for full robustness, users should use Tier 1 global-ID mode

### 8.5 Geogram Multiprecision

Geogram multiprecision components are used as the exact fallback for predicate
evaluation in the robust pipeline.

- Geogram PSM exact arithmetic
- Geogram PCK-related components included in the vendored code layout used by
  this project

### 8.6 Error-Free Transformations (EFT)

The EFT pipeline is informed by the compensated-arithmetic literature used for
accurate determinant and sign evaluation.

- Chen, H., Ullrich, P. A., and Panetta, J. (2025):
Fast and Accurate Intersections on a Sphere,
arXiv:2510.09892
https://arxiv.org/abs/2510.09892
- Ogita, T., Rump, S. M., and Oishi, S. (2005):
Accurate sum and dot product,
SIAM J. Sci. Comput., 26(6), 1955–1988
https://doi.org/10.1137/030601818
- Also based on: Kahan, W. (1996). Lecture Notes on the Status of IEEE 754.
  http://www.cs.berkeley.edu/~wkahan/ieee754status/IEEE754.PDF

### 8.7 Linear Algebra Backend

The project depends on the Eigen library for fixed-size vectors and the
templated EFT implementation.

- Eigen library

### 8.8 Summary Mapping

- robust pipeline -> Shewchuk + Geogram
- SoS -> Edelsbrunner & Mucke
- non-global-ID -> Hormann & Agathos (adapted)
- EFT -> Ogita-Rump-Oishi + Kahan
- spherical algorithm -> Chen (2026)
