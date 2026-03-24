# spherical-pip

## 1. Project Overview

`spherical-pip` is a C++ library for:

- spherical point-in-polygon (PIP)
- robust geometric predicates on the sphere

The main target problem is classifying a query point on the unit sphere as:

- `Outside`
- `Inside`
- `OnVertex`
- `OnEdge`

Key properties:

- robust to common degeneracies
- two precision pipelines:
  - a pure robust pipeline
  - an EFT pipeline based on compensated arithmetic
- designed for correctness first, with performance-aware implementation choices


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

When `s_qR_A` or `s_qR_B` is zero in the non-global-ID mode, the algorithm falls back to a half-open convention rather than symbolic perturbation. This avoids double counting at ray/vertex configurations.

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
- use a global-ID SoS mode so endpoint degeneracies can be resolved symbolically where applicable


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

Important design rule:

- the algorithm is intended to be identical across pipelines
- only the numerical evaluation of predicate signs differs


## 5. Non-global-ID Strategy

When global IDs are not provided:

- no Simulation of Simplicity (SoS) is used
- degeneracy is handled by a deterministic geometric rule

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

This is a standard half-open convention in the style of Hormann and Agathos (2001).

Important limitation:

- no symbolic perturbation is used in this mode
- behavior is deterministic and consistent
- but it is not equivalent to full global symbolic robustness


## 6. Robustness Tiers With Using Global ID

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


## 9. Usage (Minimal)

Basic point-in-polygon call without global IDs:

```cpp
#include <array>
#include <vector>

#include "spip/algorithms/point_in_polygon_sphere.hpp"

int main() {
  const std::vector<std::array<double, 3>> poly = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };
  const std::array<double, 3> q = {
      0.5773502691896257,
      0.5773502691896257,
      0.5773502691896257,
  };

  const spip::pip::Location loc = spip::pip::point_in_polygon_sphere(q, poly);
  return static_cast<int>(loc);
}
```

Global-ID call, Tier 2:

```cpp
#include <array>
#include <cstdint>
#include <vector>

#include "spip/algorithms/point_in_polygon_sphere.hpp"

int main() {
  const std::vector<std::array<double, 3>> poly = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };
  const std::vector<std::int64_t> vertex_ids = {10, 20, 30};
  const std::array<double, 3> q = {
      0.5773502691896257,
      0.5773502691896257,
      0.5773502691896257,
  };

  const spip::pip::Location loc =
      spip::pip::point_in_polygon_sphere(q, 40, poly, vertex_ids);
  return static_cast<int>(loc);
}
```


## 10. Tests and Visualization

Tests live in [tests/](/global/u1/h/hyvchen/spherical-pip/tests).

Notable coverage includes:

- basic spherical PIP cases
- robustness tiers with global IDs
- EFT path compilation and runtime behavior
- a more complicated polygon/query configuration

The Mathematica notebook
[tests/test_pip_complicated_visualization.nb](/global/u1/h/hyvchen/spherical-pip/tests/test_pip_complicated_visualization.nb)
visualizes:

- the polygon
- the query points
- the query-to-ray great-circle arcs

That notebook is the visualization companion for the complicated PIP test case.


## 11. Acknowledgements

This repository uses Jonathan Shewchuk’s robust geometric predicates:

- Jonathan Richard Shewchuk,
  "Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates",
  1997
- the vendored `predicates.c` is in the public domain
- see [third_party/LICENSE_shewchuk.txt](/global/u1/h/hyvchen/spherical-pip/third_party/LICENSE_shewchuk.txt)

This repository also uses Geogram multiprecision components:

- Geogram PSM exact arithmetic
- Geogram PCK-related support included in the vendored code layout used by this project

The project depends on Eigen for fixed-size vector and packet-related types used in the C++ implementation.
