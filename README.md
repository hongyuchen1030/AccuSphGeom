
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

The test is parity-based. Given a query point `q`, a ray endpoint `R`, and an edge `AB` (minor great-circle arc), the algorithm counts intersections with the arc `qR`. An odd count indicates `Inside`.

### Orientation Logic

The implementation uses four orientation signs:

  - $s_{qR,A} = \text{orient}(q, R, A)$
  - $s_{qR,B} = \text{orient}(q, R, B)$
  - $s_{AB,q} = \text{orient}(A, B, q)$
  - $s_{AB,R} = \text{orient}(A, B, R)$

In nondegenerate cases, the **4-sign crossing theorem** is applied. Boundary cases (`OnVertex`/`OnEdge`) are handled explicitly before parity counting.

### Ray Endpoint ($R$) Policy

`R` is a perturbed antipode of `q`.

  * **Design Decision:** If any edge yields $s_{AB,R} == 0$, the algorithm throws an error rather than retrying, ensuring path-independent classification.
  * **Recommended Fix:** Split the polygon into smaller segments or use **Tier 1 Global-ID mode**.

-----

## 3\. Robustness Tiers (With Global IDs)

When global IDs are available, **Simulation of Simplicity (SoS)** resolves ray-vertex degeneracies.

| Tier | Name | Description |
| :--- | :--- | :--- |
| **Tier 1** | Full Global | User provides $q$, $R$, and vertex IDs. Strongest mode. |
| **Tier 2** | Semi-Specified | User provides $q$ and vertex IDs; Library infers $R$ and its ID. |
| **Tier 3** | Local/Internal | User provides vertex IDs; Library infers IDs for $q$ and $R$. |

-----

## 4\. Precision Pipelines

| Pipeline | Technology | Description |
| :--- | :--- | :--- |
| **Pure Robust** | Shewchuk + Geogram | Uses adaptive predicates and exact multiprecision fallback. |
| **EFT** | Ogita-Rump-Oishi | Uses Error-Free Transformations and compensated arithmetic. |

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

Detailed tests are located in [`tests/`](https://www.google.com/search?q=%5Bhttps://github.com/hongyuchen1030/Spherical-Point-In-Polygon/tree/main/tests%5D\(https://github.com/hongyuchen1030/Spherical-Point-In-Polygon/tree/main/tests\)).

  * **Unit Tests:** Covers both pipelines and all robustness tiers.
  * **Visualization:** `tests/test_pip_complicated_visualization.nb` (Mathematica) visualizes the polygon, query points, and great-circle arcs for complex cases.

-----

## 8. References and Acknowledgements

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

