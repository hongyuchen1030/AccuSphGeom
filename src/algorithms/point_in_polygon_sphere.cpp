#include <spip/algorithms/point_in_polygon_sphere.hpp>

#include <array>
#include <cstddef>
#include <cmath>      // std::sqrt
#include <stdexcept>  // invalid_argument, domain_error
#include <string>
#include <vector>

#include <spip/predicates/orient3d.hpp>

namespace spip::pip {

namespace {

using spip::predicates::orient3d_on_sphere;
using spip::predicates::Sign;

inline void require_nonnull3(const double* p, const char* name) {
  if (!p) {
    throw std::invalid_argument(std::string("point_in_polygon_sphere: null pointer for ") + name);
  }
}

inline bool equal3(const double* a, const double* b) {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

// -----------------------------------------------------------------------------
// TEMP helper for "on minor arc".
// Placeholder only. Replace dot computations with your robust dot sign predicate.
// -----------------------------------------------------------------------------
inline double dot3_naive(const double* u, const double* v) {
  return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

// q lies on the *minor* great-circle arc AB if:
//   (1) q is on the supporting great circle: orient(A,B,q) == 0
//   (2) q is between A and B on the shorter arc (simple dot-based test for now)
inline bool on_minor_arc_simple(const double* q, const double* A, const double* B) {
  if (orient3d_on_sphere(A, B, q) != Sign::Zero) return false;

  // Between-ness on minor arc for unit vectors:
  // q is between A and B on the shorter arc iff:
  //   A·q >= A·B and B·q >= A·B
  //
  // NOTE: placeholder, swap with robust dot predicate later.
  const double ab = dot3_naive(A, B);
  const double aq = dot3_naive(A, q);
  const double bq = dot3_naive(B, q);
  return (aq >= ab) && (bq >= ab);
}

// -----------------------------------------------------------------------------
// Ray endpoint: antipode with deterministic perturbation.
// Placeholder strategy. Replace with your preferred perturb policy.
// -----------------------------------------------------------------------------
inline void make_perturbed_antipode_simple(const double* q, double* r_out) {
  // start from antipode
  r_out[0] = -q[0];
  r_out[1] = -q[1];
  r_out[2] = -q[2];

  // deterministic perturbation (placeholder)
  const double ax = (r_out[0] < 0.0) ? -r_out[0] : r_out[0];
  const double ay = (r_out[1] < 0.0) ? -r_out[1] : r_out[1];
  const double az = (r_out[2] < 0.0) ? -r_out[2] : r_out[2];

  constexpr double eps = 1e-15;
  if (ax <= ay && ax <= az) r_out[0] += eps;
  else if (ay <= ax && ay <= az) r_out[1] += eps;
  else r_out[2] += eps;

  // renormalize to stay on S^2 (placeholder)
  const double n2 = r_out[0] * r_out[0] + r_out[1] * r_out[1] + r_out[2] * r_out[2];
  if (n2 > 0.0) {
    const double inv = 1.0 / std::sqrt(n2);
    r_out[0] *= inv;
    r_out[1] *= inv;
    r_out[2] *= inv;
  }
}

// -----------------------------------------------------------------------------
// We consider the ray minor arc from query point q to ray endpoint R (perturbed
// antipode). For a polygon edge minor arc AB, this predicate returns true iff
// AB contributes one crossing to the ray crossing count, using the half-open
// convention of Hormann and Agathos (2001) to avoid double counting at vertices.
//
// The underlying arc–arc crossing characterization can be written using 4
// orient3D signs (Chen et al., 2026 preprint):
//   orient(C,D,A) = orient(A,B,D) = -orient(A,B,C) = -orient(C,D,B),
// where we set C=q and D=R.
//
// References:
//   Chen, H., Ullrich, P. A., Panetta, J., Marsico, D., Hanke, M., Jain, R.,
//   Zhang, C., and Jacob, R. L. Accurate and Robust Geometric Algorithms for
//   Regridding on the Sphere. EGUsphere [preprint], 2026.
//   DOI: https://doi.org/10.5194/egusphere-2026-636
//
//   Hormann, K., and Agathos, A. The point in polygon problem for arbitrary
//   polygons. Computational Geometry, 20(3), 2001.
//   DOI: https://doi.org/10.1016/S0925-7721(01)00012-8
// -----------------------------------------------------------------------------
inline bool counts_as_ray_crossing_half_open(const double* A,
                                            const double* B,
                                            const double* q,
                                            const double* R) {
  // Variable mapping from Chen et al to Hormann and Agathos
  //   s1 = orient(C,D,A) -> orient(q,R,A)
  //   s4 = orient(C,D,B) -> orient(q,R,B)
  //   s3 = orient(A,B,C) -> orient(A,B,q)
  //   s2 = orient(A,B,D) -> orient(A,B,R)
  const Sign s_qR_A = orient3d_on_sphere(q, R, A);
  const Sign s_qR_B = orient3d_on_sphere(q, R, B);
  const Sign s_AB_q = orient3d_on_sphere(A, B, q);
  const Sign s_AB_R = orient3d_on_sphere(A, B, R);

  // Step 1: half-open straddle test about the ray great circle.
  // below(x) := (x < 0). Zero is treated as non-negative.
  const bool below_a = (s_qR_A == Sign::Negative);
  const bool below_b = (s_qR_B == Sign::Negative);
  if (!(below_a ^ below_b)) return false;

  // Step 2: in-front test about the edge great circle.
  // After vertex/edge prechecks, these should not be Zero. If they are, the
  // caller's boundary handling is inconsistent with the prechecks.
  if (s_AB_q == Sign::Zero) {
    throw std::domain_error(
        "counts_as_ray_crossing_half_open: q lies on AB great circle (boundary case)");
  }
  if (s_AB_R == Sign::Zero) {
    throw std::domain_error(
        "counts_as_ray_crossing_half_open: R lies on AB great circle (ray degeneracy)");
  }

  const bool q_side = (s_AB_q == Sign::Negative);
  const bool r_side = (s_AB_R == Sign::Negative);
  return (q_side ^ r_side);
}

}  // namespace

// -----------------------------------------------------------------------------
// Core: pointer-based polygon as array of pointers to 3 doubles
// poly[i] points to 3 doubles (x,y,z). poly is CCW on S^2.
// -----------------------------------------------------------------------------
Location point_in_polygon_sphere(const double* q,
                                 const double* const* poly,
                                 std::size_t n) {
  require_nonnull3(q, "q");
  if (!poly) {
    throw std::invalid_argument("point_in_polygon_sphere: poly is null");
  }
  if (n < 3) {
    throw std::invalid_argument("point_in_polygon_sphere: polygon must have >= 3 vertices");
  }

  // Validate polygon vertex pointers once
  for (std::size_t i = 0; i < n; ++i) {
    require_nonnull3(poly[i], "poly[i]");
  }

  // 1) Vertex test (exact equality)
  for (std::size_t i = 0; i < n; ++i) {
    if (equal3(q, poly[i])) {
      return Location::OnVertex;
    }
  }

  // 2) Edge test (simple minor-arc test for now)
  for (std::size_t i = 0; i < n; ++i) {
    const double* A = poly[i];
    const double* B = poly[(i + 1) % n];
    if (on_minor_arc_simple(q, A, B)) {
      return Location::OnEdge;
    }
  }

  // 3) Ray endpoint: perturbed antipode
  double r_arr[3];
  make_perturbed_antipode_simple(q, r_arr);
  const double* R = r_arr;

  // 4) Parity counting (half-open rule)
  bool inside = false;
  for (std::size_t i = 0; i < n; ++i) {
    const double* A = poly[i];
    const double* B = poly[(i + 1) % n];

    if (counts_as_ray_crossing_half_open(A, B, q, R)) {
      inside = !inside;
    }
  }

  return inside ? Location::Inside : Location::Outside;
}

// -----------------------------------------------------------------------------
// Adapter overloads
// -----------------------------------------------------------------------------
Location point_in_polygon_sphere(const std::array<double, 3>& q,
                                 const std::vector<std::array<double, 3>>& poly) {
  if (poly.size() < 3) {
    throw std::invalid_argument("point_in_polygon_sphere: polygon must have >= 3 vertices");
  }

  std::vector<const double*> ptrs;
  ptrs.reserve(poly.size());
  for (const auto& v : poly) ptrs.push_back(v.data());

  return point_in_polygon_sphere(q.data(), ptrs.data(), ptrs.size());
}

Location point_in_polygon_sphere(const std::vector<double>& q,
                                 const std::vector<std::vector<double>>& poly) {
  if (q.size() != 3) {
    throw std::invalid_argument("point_in_polygon_sphere: q must have size 3");
  }
  if (poly.size() < 3) {
    throw std::invalid_argument("point_in_polygon_sphere: polygon must have >= 3 vertices");
  }

  std::vector<const double*> ptrs;
  ptrs.reserve(poly.size());
  for (std::size_t i = 0; i < poly.size(); ++i) {
    if (poly[i].size() != 3) {
      throw std::invalid_argument("point_in_polygon_sphere: poly[i] must have size 3");
    }
    ptrs.push_back(poly[i].data());
  }

  return point_in_polygon_sphere(q.data(), ptrs.data(), ptrs.size());
}

}  // namespace spip::pip