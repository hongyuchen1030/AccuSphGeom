
# Spherical-Point-In-Polygon

`Spherical-Point-In-Polygon` is a C++ library for robust **spherical point-in-polygon (PIP)** classification and geometric predicates on the unit sphere.

The library classifies a query point as:

  * **`Outside`**
  * **`Inside`**
  * **`OnVertex`**
  * **`OnEdge`**

### Key Properties

  * **Robustness:** Handles common geometric degeneracies and numerical instabilities.
  * **Tiered Reliability:** Supports three robustness tiers when global IDs are available.
  * **Precision Pipelines:**
      * **Pure Robust:** Adaptive predicates (Shewchuk) + exact fallback (Geogram).
      * **EFT:** Compensated arithmetic based on Error-Free Transformations.
  * **Header-Only:** Designed for "correctness first" with easy integration.

This repository contains the implementation. The associated paper discusses the
algorithms in more detail and provides broader algorithmic context:

- Chen, H. (2026). Accurate and Robust Algorithms for Spherical Polygon
  Operations. EGUsphere preprint.
  https://egusphere.copernicus.org/preprints/2026/egusphere-2026-636/

-----

## 1\. Installation & Quick Start

### Installation

This library is header-only. No installation step or package manager is required.

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/hongyuchen1030/Spherical-Point-In-Polygon.git](https://github.com/hongyuchen1030/Spherical-Point-In-Polygon.git)
    ```

3.  **Include the header:**
    ```cpp
    include <spip/algorithms/point_in_polygon_sphere.hpp>
    ```

5.  **Compiler Path:** Add the `include/` directory to your include path (C++17 required).
    ```bash
    g++ program.cpp -I/path/to/Spherical-Point-In-Polygon/include -std=c++17
    ```

### Quick Start Examples

**Standard Usage (No Global IDs)**

```cpp
#include <array>
#include <spip/algorithms/point_in_polygon_sphere.hpp>

const std::array<double, 3> q = {0.5773502691896257, 0.5773502691896257, 0.5773502691896257};
const std::array<std::array<double, 3>, 3> poly = {{
    {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}
}};

const auto loc = spip::pip::point_in_polygon_sphere(q, poly);
```

**Tier 1 Global-ID Usage**

```cpp
#include <cstdint>
#include <spip/algorithms/point_in_polygon_sphere.hpp>

double q[3] = {0.5773502691896257, 0.5773502691896257, 0.5773502691896257};
double R[3] = {-0.5773502691896257, -0.5773502691896257, -0.5773502691896257};
double poly_storage[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
const double* poly[3] = {poly_storage[0], poly_storage[1], poly_storage[2]};
const std::int64_t vertex_ids[3] = {10, 20, 30};

const auto loc = spip::pip::point_in_polygon_sphere(q, 40, R, 50, poly, vertex_ids, 3);
```


This library provides a complete API coverage across precision models, robustness modes, and container interfaces. In total, there are **24 API combinations**, formed by:

* **2 precision models**:

  * robust (adaptive exact predicates)
  * EFT (compensated arithmetic)

* **4 robustness modes**:

  * no global ID
  * Tier 1 (full global robustness)
  * Tier 2 (semi-specified global robustness)
  * Tier 3 (local/internal robustness)

* **3 API overload families**:

  * raw pointer
  * std::array
  * std::vector

These combinations ensure that users can choose the appropriate trade-off between robustness, performance, and interface convenience.

| Precision Model | Robustness Mode | Raw Pointer API | std::array API | std::vector API |
| --------------- | --------------- | --------------- | -------------- | --------------- |
| Robust          | No Global ID    | ✓               | ✓              | ✓               |
| Robust          | Tier 1          | ✓               | ✓              | ✓               |
| Robust          | Tier 2          | ✓               | ✓              | ✓               |
| Robust          | Tier 3          | ✓               | ✓              | ✓               |
| EFT             | No Global ID    | ✓               | ✓              | ✓               |
| EFT             | Tier 1          | ✓               | ✓              | ✓               |
| EFT             | Tier 2          | ✓               | ✓              | ✓               |
| EFT             | Tier 3          | ✓               | ✓              | ✓               |

Total combinations:
2 (precision) × 4 (robustness) × 3 (overloads) = **24 API cases**

All of these combinations are explicitly tested in:

[tests/test_library_usage.cpp](tests/test_library_usage.cpp)

-----

## 2\. Core Algorithm (Spherical PIP)

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

-----

## 3\. Robustness Tiers (With Global IDs)

When global IDs are available, **Simulation of Simplicity (SoS)** resolves ray-vertex degeneracies.

| Tier | Name | Description |
| :--- | :--- | :--- |
| **Tier 1** | Full Global | User provides $q$, $R$, and vertex IDs. Strongest mode. |
| **Tier 2** | Semi-Specified | User provides $q$ and vertex IDs; Library infers $R$ and its ID. |
| **Tier 3** | Local/Internal | User provides vertex IDs; Library infers IDs for $q$ and $R$. |

-----

## 4. Precision Pipelines

| Pipeline        | Technology                | Description                                                                                                                    |
| :-------------- | :------------------------ | :----------------------------------------------------------------------------------------------------------------------------- |
| **Pure Robust** | Shewchuk + Geogram        | Uses adaptive orientation predicates with filtered evaluation and exact multiprecision fallback.                               |
| **EFT**         | EFT + Shewchuk (fallback) | Uses compensated arithmetic (error-free transformations) for most evaluations, with robust predicate fallback near degeneracy. |

The two pipelines implement **the same algorithmic logic**.
They differ only in how predicate signs are evaluated.

### Pure Robust Pipeline

* Uses Shewchuk’s adaptive `orient3d`
* Uses filtered quadruple product
* Falls back to exact multiprecision (Geogram) when needed
* Fully robust by construction

Got it — here is the **corrected raw markdown text** with:

* **no bullet references**
* **only in-text citations**
* clean scientific tone
* no re-rendering issues

---

## 4. Precision Pipelines

| Pipeline        | Technology                | Description                                                                                                                    |
| :-------------- | :------------------------ | :----------------------------------------------------------------------------------------------------------------------------- |
| **Pure Robust** | Shewchuk + Geogram        | Uses adaptive orientation predicates with filtered evaluation and exact multiprecision fallback.                               |
| **EFT**         | EFT + Shewchuk (fallback) | Uses compensated arithmetic (error-free transformations) for most evaluations, with robust predicate fallback near degeneracy. |

The two pipelines implement **the same algorithmic logic**.
They differ only in how predicate signs are evaluated.

### Pure Robust Pipeline

* Uses Shewchuk’s adaptive `orient3d`
* Uses filtered quadruple product
* Falls back to exact multiprecision (Geogram) when needed
* Fully robust by construction

### EFT Pipeline

The EFT pipeline uses compensated arithmetic to simulate higher working precision and reduce reliance on branching-based filtering. The implementation follows the general approach of error-free transformations and compensated evaluation as described in Chen et al. (2025) and Ogita et al. (2005).

Core idea:

* Use **error-free transformations (EFT)** to approximate double working precision
* Keep computation largely **branch-free and SIMD-friendly**
* Avoid frequent fallback to exact arithmetic in non-degenerate cases

### Zero Detection and Fallback

In the EFT pipeline:

* A value is treated as zero if:

  abs(value) < tol

* The default tolerance is:

  tol = 1e-8

When a value falls within this tolerance:

* The computation falls back to the robust predicate:

  * Shewchuk adaptive `orient3d`
  * with exact fallback when necessary

This hybrid design:

* preserves high performance in regular cases
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

