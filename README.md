# Spherical-Point-In-Polygon

## 1. Project Overview

`Spherical-Point-In-Polygon` is a C++ library for:

- spherical point-in-polygon (PIP)
- robust geometric predicates on the sphere

The main target problem is classifying a query point on the unit sphere as:

- `Outside`
- `Inside`
- `OnVertex`
- `OnEdge`

Key properties:

- robust to common degeneracies
- three robustness tiers when global IDs are available
- two precision pipelines:
  - a pure robust pipeline
  - an EFT pipeline based on compensated arithmetic
- designed for correctness first, with performance-aware implementation choices

This repository contains the implementation. The associated paper discusses the
algorithms in more detail and provides broader algorithmic context:


- Chen, H. (2026). Accurate and Robust Algorithms for Spherical Polygon
  Operations. EGUsphere preprint.
  https://egusphere.copernicus.org/preprints/2026/egusphere-2026-636/


## 2. Core Algorithm (Spherical PIP)

The point-in-polygon test is parity-based. Given:

- `q`: the query point on the sphere
- `R`: a ray endpoint on the sphere
- `AB`: a polygon edge, interpreted as the minor great-circle arc from vertex `A` to vertex `B`

the algorithm counts how many polygon edges intersect the minor arc `qR`. An odd number of crossings means `Inside`; an even number means `Outside`.

The ray endpoint `R` is normally chosen as a perturbed antipode of `q`, so the ray is geometrically well separated from the query.

For each polygon edge `AB`, the crossing logic uses four orientation signs:

- `s_qR_A = orient(q, R, A)`
- `s_qR_B = orient(q, R, B)`
- `s_AB_q = orient(A, B, q)`
- `s_AB_R = orient(A, B, R)`

In the nondegenerate case, the implementation uses the strict 4-sign crossing theorem: the arcs `qR` and `AB` cross if and only if the endpoint orientations are strictly separated on both supporting great circles. In code, this is the branch where all four signs are nonzero and the sign pattern is checked directly.

When `s_qR_A` or `s_qR_B` is zero in the non-global-ID mode, the algorithm falls back to a half-open convention rather than symbolic perturbation. This is a deterministic local tie-breaking rule for degenerate ray/vertex configurations and avoids double counting.

Boundary handling is performed before parity counting:

- if `q` exactly matches a polygon vertex, return `OnVertex`
- if `q` lies on a polygon edge minor arc, return `OnEdge`


## 3. Ray Endpoint (R) Policy

The ray endpoint `R` is constructed as a perturbed antipode of `q`:

- start from `-q`
- perturb the least dominant coordinate
- renormalize to the sphere

The perturbation is needed to avoid the most obvious antipodal degeneracies.

Design decision:

If any polygon edge yields

`s_AB_R == 0`

then the algorithm throws an error.

Interpretation:

- this typically indicates a polygon that is too large, close to hemispherical, or otherwise poorly separated from the ray construction
- retrying with another `R` is intentionally not used, because that would make the classification rule path-dependent and harder to reason about

What users should do instead:

- split the polygon into smaller pieces, or
- use Tier 1 global-ID mode so endpoint degeneracies can be resolved symbolically where applicable


## 4. Precision Pipelines

The repository provides two precision pipelines.

### (1) Pure robust pipeline

- Shewchuk adaptive orientation predicates
- filtered quadruple product
- exact fallback via Geogram multiprecision

### (2) EFT pipeline

- compensated arithmetic based on error-free transforms
- same high-level algorithmic structure
- tolerance-based zero detection

**Important design rule:**

- The algorithm is intended to be identical across pipelines.
- Only the numerical evaluation of predicate signs differs.

The EFT implementation is largely informed by:

- Chen, H., Ullrich, P. A., and Panetta, J. (2025):  
  *Fast and Accurate Intersections on a Sphere*,  
  arXiv:2510.09892, https://arxiv.org/abs/2510.09892.

- Ogita, T., Rump, S. M., and Oishi, S. (2005):  
  *Accurate sum and dot product*,  
  SIAM Journal on Scientific Computing, 26(6), 1955–1988,  
  https://doi.org/10.1137/030601818.

Both are related to:

- W. Kahan, "Lecture Notes on the Status of IEEE-754", 1996,
  http://www.cs.berkeley.edu/~wkahan/ieee754status/IEEE754.PDF

in particular for the 2x2 determinant and FMA idea.


## 5. Non-global-ID Strategy

When global IDs are not provided:

- no Simulation of Simplicity (SoS) is used
- degeneracy is handled by a deterministic geometric rule
- this is local deterministic tie-breaking, not full global symbolic robustness

The non-global-ID crossing policy is:

- `s_AB_q == 0` -> treated as not crossing
  - boundary cases have already been handled by the explicit `OnVertex` / `OnEdge` checks
- `s_AB_R == 0` -> error
  - the ray construction is considered invalid for this polygon/query configuration
- if `s_qR_A` and `s_qR_B` are both nonzero
  - use the strict 4-sign theorem
- otherwise
  - use the half-open rule

Half-open rule:

```text
below_a = (s_qR_A < 0)
below_b = (s_qR_B < 0)

crossing only if:
(below_a XOR below_b)
```

This is a standard half-open convention in the style of Hormann and Agathos
(2001).

Important limitation:

- Hormann and Agathos discuss the half-open rule in the planar setting
- this repository works on the spherical surface
- this repository does not currently claim a global theoretical proof for the
  spherical no-global-ID tie-breaking rule
- in this repository, the half-open rule is used as a deterministic tie-breaking
  convention when no global ID is present and a degenerate ray/vertex case
  occurs
- behavior is deterministic and consistent
- but it is not equivalent to full global symbolic robustness
- if full global robustness is required, use Tier 1 global-ID mode


## 6. Robustness Tiers With Global IDs

When global IDs are available, the library supports three robustness tiers.

### Tier 1: full global robustness

User provides:

- polygon vertex coordinates and vertex IDs
- query point `q` and `q_id`
- outside point `R` and `r_id`

Properties:

- all participating points belong to one explicit symbolic ordering
- this is the strongest and most predictable mode
- degeneracies are resolved relative to caller-controlled global IDs

### Tier 2: semi-specified robustness

User provides:

- polygon vertex coordinates and vertex IDs
- query point `q` and `q_id`

Library does:

- infer `R`
- assign an internal ID to `R`

Properties:

- still globally consistent with respect to the polygon and `q`
- only `R` is synthetic

### Tier 3: local/internal robustness

User provides:

- polygon vertex coordinates and vertex IDs

Library does:

- assign an internal ID to `q`
- infer `R`
- assign an internal ID to `R`

Properties:

- deterministic
- consistent for a given call
- ordering involving `q` and `R` is library-defined
- this is local robustness only, not full global SoS


## 7. SoS (Simulation of Simplicity)

SoS is used only in global-ID paths.

It is applied only to the endpoint degeneracies:

- `s_qR_A == 0`
- `s_qR_B == 0`

It is not applied to:

- `s_AB_q`
- `s_AB_R`

Those remain exact geometric tests:

- `s_AB_q == 0` -> not crossing
- `s_AB_R == 0` -> error

Purpose of SoS:

- resolve ray-vertex degeneracy using symbolic ordering
- preserve deterministic crossing decisions when the ray passes through a vertex-level degenerate configuration


## 8. Internal ID Assignment

When the library assigns IDs internally:

1. compute the maximum polygon vertex ID
2. if `max_id + 2` is safe:
   - assign `q_id = max_id + 1` when needed
   - assign `r_id = max_id + 2`
3. otherwise:
   - compute a MEX over nonnegative integers
   - use the first unused ID for `q_id` when needed
   - use the second unused ID for `r_id`

Constraints:

- no negative IDs are used
- polygon vertex IDs are never modified
- internal IDs are always nonnegative


## 9. API Reference and Usage

### 9.1 Precision Pipelines

The library provides two precision pipelines:

- robust pipeline: the default double-based implementation
- EFT pipeline: templated evaluation using `spip::V3_T<T>` and compensated arithmetic

The algorithm is identical in both pipelines. Only the numerical evaluation of
predicate signs differs.

Overload families:

- robust pipeline: raw pointers, `std::array<std::array<double, 3>, N>`,
  `std::vector<std::array<double, 3>>`, and `std::vector<std::vector<double>>`
- EFT pipeline: `V3_T<T>*`, `std::vector<V3_T<T>>`, and
  `std::array<V3_T<T>, N>`

### 9.2 Return Type

`spip::pip::Location` has four values:

- `Location::Inside`
- `Location::Outside`
- `Location::OnEdge`
- `Location::OnVertex`

For complete API usage examples, including all overloads and robustness tiers,
see `tests/test_pip_basic.cpp`.

### 9.3 Representative Examples

#### 9.3.1 Robust Pipeline (No Global ID)

A representative robust no-global-ID example:

```cpp
#include <array>

#include "spip/algorithms/point_in_polygon_sphere.hpp"

const std::array<double, 3> q = {0.5773502691896257, 0.5773502691896257,
                                 0.5773502691896257};
const std::array<std::array<double, 3>, 3> poly = {{
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
}};

const auto loc = spip::pip::point_in_polygon_sphere(q, poly);
```

#### 9.3.2 Robust Pipeline (Tier 1)

A representative robust Tier 1 example:

```cpp
#include <cstdint>

#include "spip/algorithms/point_in_polygon_sphere.hpp"

double q[3] = {0.5773502691896257, 0.5773502691896257, 0.5773502691896257};
double R[3] = {-0.5773502691896257, -0.5773502691896257, -0.5773502691896257};
double poly_storage[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
const double* poly[3] = {poly_storage[0], poly_storage[1], poly_storage[2]};
const std::int64_t vertex_ids[3] = {10, 20, 30};

const auto loc =
    spip::pip::point_in_polygon_sphere(q, 40, R, 50, poly, vertex_ids, 3);
```

#### 9.3.3 EFT Pipeline

A representative EFT example:

```cpp
#include <array>

#include "spip/algorithms/point_in_polygon_sphere.hpp"

using spip::V3_T;

const V3_T<double> q(0.5773502691896257, 0.5773502691896257, 0.5773502691896257);
const std::array<V3_T<double>, 3> poly = {
    V3_T<double>(1.0, 0.0, 0.0),
    V3_T<double>(0.0, 1.0, 0.0),
    V3_T<double>(0.0, 0.0, 1.0),
};

const auto loc = spip::pip::point_in_polygon_sphere(q, poly);
```

### 9.6 Important Notes

- global IDs must be unique
- IDs must be non-negative integers
- vertex IDs are never modified
- `q` and `R` IDs may be internally assigned
- for full robustness, use Tier 1
- `s_AB_R == 0` triggers an error for an invalid face/ray configuration

### 9.7 Link to Tests

See `tests/test_pip_complicated.cpp` for a more involved case, and
`tests/test_pip_complicated_visualization.nb` for the corresponding
visualization notebook.


## 10. Tests and Visualization

Tests live in [tests/](tests/) and in the repository test directory on GitHub:

- https://github.com/hongyuchen1030/Spherical-Point-In-Polygon/tree/main/tests

Notable coverage includes:

- robust pipeline unit tests
- EFT pipeline unit tests
- robustness tiers with global IDs
- a more complicated polygon/query configuration

The Mathematica notebook
[tests/test_pip_complicated_visualization.nb](tests/test_pip_complicated_visualization.nb)
visualizes:

- the polygon
- the query points
- the query-to-ray great-circle arcs

That notebook is the visualization companion for the complicated PIP test case.


## 11. References and Acknowledgements

### 11.1 Associated Paper

The algorithms implemented in this repository are discussed in detail in:

- Chen, H. (2026). Accurate and Robust Algorithms for Spherical Polygon
  Operations. EGUsphere preprint.
  https://egusphere.copernicus.org/preprints/2026/egusphere-2026-636/

This repository contains the implementation. The paper provides the full
algorithmic context and derivations.

### 11.2 Robust Geometric Predicates (Shewchuk)

The robust pipeline uses Jonathan Shewchuk's adaptive predicates for core
orientation-sign evaluation, including the vendored `predicates.c`
implementation.

- Shewchuk, J. R. (1997). Adaptive Precision Floating-Point Arithmetic and Fast
  Robust Geometric Predicates. Discrete & Computational Geometry, 18(3),
  305-363.

Notes:

- the vendored `predicates.c` is public domain
- see [third_party/LICENSE_shewchuk.txt](third_party/LICENSE_shewchuk.txt)

### 11.3 Simulation of Simplicity (SoS)

Simulation of Simplicity is used for global-ID degeneracy resolution in the
Tier 1, Tier 2, and Tier 3 global-ID paths.

- Edelsbrunner, H., & Mucke, E. P. (1990). Simulation of Simplicity: A
  Technique to Cope with Degenerate Cases in Geometric Algorithms. ACM
  Transactions on Graphics, 9(1), 66-104.
  https://doi.org/10.1145/77635.77639

### 11.4 Half-Open Rule

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

### 11.5 Geogram Multiprecision

Geogram multiprecision components are used as the exact fallback for predicate
evaluation in the robust pipeline.

- Geogram PSM exact arithmetic
- Geogram PCK-related components included in the vendored code layout used by
  this project

### 11.6 Error-Free Transformations (EFT)

The EFT pipeline is informed by the compensated-arithmetic literature used for
accurate determinant and sign evaluation.

- Ogita, T., Rump, S. M., & Oishi, S. (2005). Accurate Sum and Dot Product.
  SIAM Journal on Scientific Computing, 26(6), 1955-1988.
  https://doi.org/10.1137/030601818
- Additional reference: https://arxiv.org/abs/2510.09892
- Based on: Kahan, W. (1996). Lecture Notes on the Status of IEEE 754.
  http://www.cs.berkeley.edu/~wkahan/ieee754status/IEEE754.PDF

### 11.7 Linear Algebra Backend

The project depends on the Eigen library for fixed-size vectors and the
templated EFT implementation.

- Eigen library

### 11.8 Summary Mapping

- robust pipeline -> Shewchuk + Geogram
- SoS -> Edelsbrunner & Mucke
- non-global-ID -> Hormann & Agathos (adapted)
- EFT -> Ogita-Rump-Oishi + Kahan
- spherical algorithm -> Chen (2026)
