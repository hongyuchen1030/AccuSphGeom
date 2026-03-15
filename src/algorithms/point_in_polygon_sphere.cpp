#include <spip/algorithms/point_in_polygon_sphere.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <spip/kernels/pip_kernel_adaptive.hpp>


// Debug only

#include <iomanip>
#include <limits>

namespace spip::pip {

namespace {

using Kernel = spip::kernels::PIPKernelAdaptive;
using Sign = Kernel::Sign;

// Require a valid pointer to three coordinates.
inline void require_nonnull3(const double* p, const char* name) {
  if (!p) {
    throw std::invalid_argument(
        std::string("point_in_polygon_sphere: null pointer for ") + name);
  }
}

// Exact coordinatewise equality for 3-vectors.
inline bool equal3(const double* a, const double* b) {
  return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

// Negate a predicate sign.
inline Sign flip_sign(Sign s) {
  if (s == Sign::Positive) return Sign::Negative;
  if (s == Sign::Negative) return Sign::Positive;
  return Sign::Zero;
}

// Exact sign of a stored scalar.
inline Sign exact_sign_coord(double x) {
  if (x > 0.0) return Sign::Positive;
  if (x < 0.0) return Sign::Negative;
  return Sign::Zero;
}

// Evaluate the sign of the 2 x 2 determinant
//
//   | a00 a01 |
//   | a10 a11 |
//
// by embedding it into the 3 x 3 determinant
//
//   | a00 a01 0 |
//   | a10 a11 0 |
//   |  0   0  1 | .
//
// This reuses the adaptive kernel orient3d predicate.
inline Sign exact_sign_det2(double a00, double a01,
                            double a10, double a11) {
  const double r0[3] = {a00, a01, 0.0};
  const double r1[3] = {a10, a11, 0.0};
  const double r2[3] = {0.0, 0.0, 1.0};
  return Kernel::orient3d_on_sphere(r0, r1, r2);
}

// Return true iff q lies on the closed non-antipodal minor arc AB.
//
// The test consists of two conditions:
//
// (1) q lies on the supporting great circle of AB:
//       orient(A, B, q) = 0
//
// (2) q lies between A and B on the minor arc. Equivalently, the normals
//     A x q and q x B do not point in opposite directions:
//       (A x q) . (q x B) >= 0
//
// The second condition is evaluated as the scalar quadruple product
// quadruple3d(A, q, q, B).
//
// Preconditions:
// - q, A, and B are valid pointers to 3 coordinates
// - vertex coincidence q == A or q == B has already been excluded
inline bool on_minor_arc(const double* q, const double* A, const double* B) {
  if (equal3(A, B)) {
    return false;
  }
  if (Kernel::orient3d_on_sphere(A, B, q) != Sign::Zero) {
    return false;
  }
  return Kernel::quadruple3d(A, q, q, B) != Sign::Negative;
}

// Validate the raw-pointer polygon input and test whether q matches a polygon
// vertex exactly. The function returns true iff q coincides with some poly[i].
inline bool validate_polygon_and_check_vertex(const double* q,
                                              const double* const* poly,
                                              std::size_t n) {
  require_nonnull3(q, "q");
  if (!poly) {
    throw std::invalid_argument("point_in_polygon_sphere: poly is null");
  }
  if (n < 3) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: polygon must have >= 3 vertices");
  }

  for (std::size_t i = 0; i < n; ++i) {
    require_nonnull3(poly[i], "poly[i]");
    if (equal3(q, poly[i])) {
      return true;
    }
  }
  return false;
}

// Exact boundary test over all polygon edges.
//
// Preconditions:
// - polygon input has already been validated
// - exact vertex coincidence has already been excluded
inline bool point_on_polygon_edge_exact(const double* q,
                                        const double* const* poly,
                                        std::size_t n) {
  const double* A = poly[n - 1];
  for (std::size_t i = 0; i < n; ++i) {
    const double* B = poly[i];
    if (on_minor_arc(q, A, B)) {
      return true;
    }
    A = B;
  }
  return false;
}

// Construct a deterministic perturbation of the antipode of q.
//
// This endpoint is used only by the non-SoS overload. A small perturbation is
// added to the antipode and the result is renormalized to the unit sphere.
inline void make_perturbed_antipode_simple(const double* q, double* r_out) {
  r_out[0] = -q[0];
  r_out[1] = -q[1];
  r_out[2] = -q[2];

  const double ax = std::fabs(r_out[0]);
  const double ay = std::fabs(r_out[1]);
  const double az = std::fabs(r_out[2]);

  constexpr double eps = 1e-8;
  if (ax <= ay && ax <= az) {
    r_out[0] += eps;
  } else if (ay <= ax && ay <= az) {
    r_out[1] += eps;
  } else {
    r_out[2] += eps;
  }

  const double n2 =
      r_out[0] * r_out[0] +
      r_out[1] * r_out[1] +
      r_out[2] * r_out[2];

  if (n2 > 0.0) {
    const double inv_n = 1.0 / std::sqrt(n2);
    r_out[0] *= inv_n;
    r_out[1] *= inv_n;
    r_out[2] *= inv_n;
  }
}

// Return true iff the polygon edge AB contributes one crossing to the ray
// crossing count for the minor arc qR, where R is a perturbed ray endpoint.
//
// This routine uses two different decision rules.
//
// (1) Nondegenerate case: theorem check
//     If both ray-plane endpoint signs are nonzero,
//       s_qR_A = orient(q, R, A),
//       s_qR_B = orient(q, R, B),
//     then the crossing decision is made from the 4-sign consistency rule for
//     minor-arc intersection:
//
//       orient(q, R, A) = orient(A, B, R)
//                       = -orient(A, B, q)
//                       = -orient(q, R, B).
//
//     Equivalently, with
//       s_qR_A = orient(q, R, A),
//       s_qR_B = orient(q, R, B),
//       s_AB_q = orient(A, B, q),
//       s_AB_R = orient(A, B, R),
//
//     the edge AB intersects the ray minor arc qR iff
//
//       s_qR_A == s_AB_R &&
//       s_qR_A == flip_sign(s_AB_q) &&
//       s_qR_A == flip_sign(s_qR_B).
//
//     This is the correct strict minor-arc crossing test in the absence of
//     degeneracy.
//
// (2) Degenerate ray/vertex case: half-open rule
//     If s_qR_A == 0 or s_qR_B == 0, then a polygon vertex lies on the ray
//     great circle. In this case the strict 4-sign theorem does not apply.
//     Instead, we use the half-open convention of Hormann and Agathos (2001):
//
//       - an endpoint is treated as "below" iff orient(q, R, endpoint) < 0
//       - zero is treated as non-negative
//
//     so the ray-plane straddle test becomes
//
//       (s_qR_A < 0) xor (s_qR_B < 0).
//
//     If this straddle test passes, the second-stage edge-plane test is
//
//       (s_AB_q < 0) xor (s_AB_R < 0).
//
//     Since we assume here that s_AB_q and s_AB_R are already nonzero,
//     this yields the half-open crossing decision.
//
// Preconditions:
// - exact boundary prechecks have already excluded q on a polygon vertex or edge
// - s_AB_q and s_AB_R are therefore expected to be nonzero
// - the only degeneracies handled here are s_qR_A == 0 and/or s_qR_B == 0
inline bool counts_as_ray_crossing_half_open(const double* A,
                                             const double* B,
                                             const double* q,
                                             const double* R) {
  const Sign s_qR_A = Kernel::orient3d_on_sphere(q, R, A);
  const Sign s_qR_B = Kernel::orient3d_on_sphere(q, R, B);
  const Sign s_AB_q = Kernel::orient3d_on_sphere(A, B, q);
  const Sign s_AB_R = Kernel::orient3d_on_sphere(A, B, R);

  std::cout << "  orient(q,R,A) = " << static_cast<int>(s_qR_A) << '\n';
  std::cout << "  orient(q,R,B) = " << static_cast<int>(s_qR_B) << '\n';
  std::cout << "  orient(A,B,q) = " << static_cast<int>(s_AB_q) << '\n';
  std::cout << "  orient(A,B,R) = " << static_cast<int>(s_AB_R) << '\n';

  if (s_AB_q == Sign::Zero) {
    throw std::domain_error(
        "counts_as_ray_crossing_half_open: q lies on AB great circle "
        "(boundary case)");
  }
  if (s_AB_R == Sign::Zero) {
    throw std::domain_error(
        "counts_as_ray_crossing_half_open: R lies on AB great circle "
        "(ray degeneracy)");
  }

  const bool has_ray_plane_zero =
      (s_qR_A == Sign::Zero) || (s_qR_B == Sign::Zero);

  std::cout << "  has_ray_plane_zero = " << has_ray_plane_zero << '\n';

  if (!has_ray_plane_zero) {
    const bool theorem_check =
        (s_qR_A == s_AB_R) &&
        (s_qR_A == flip_sign(s_AB_q)) &&
        (s_qR_A == flip_sign(s_qR_B));

    std::cout << "  nondegenerate case: use 4-sign theorem\n";
    std::cout << "  theorem_check = " << theorem_check << '\n';
    return theorem_check;
  }

  const bool below_a = (s_qR_A == Sign::Negative);
  const bool below_b = (s_qR_B == Sign::Negative);

  std::cout << "  degenerate ray/vertex case: use half-open rule\n";
  std::cout << "  below_a = " << below_a << '\n';
  std::cout << "  below_b = " << below_b << '\n';

  if (!(below_a ^ below_b)) {
    std::cout << "  half-open straddle test failed\n";
    return false;
  }

  std::cout << "  half-open straddle test passed\n";

  const bool q_side = (s_AB_q == Sign::Negative);
  const bool r_side = (s_AB_R == Sign::Negative);

  std::cout << "  q_side = " << q_side << '\n';
  std::cout << "  r_side = " << r_side << '\n';

  const bool crossing = (q_side ^ r_side);
  std::cout << "  half-open crossing decision = " << crossing << '\n';

  return crossing;
}

// Build a strict symbolic rank ordering from global vertex IDs.
//
// Smaller rank means earlier symbolic perturbation priority. Ranks 0 and 1 are
// reserved for q and R, respectively, so polygon vertices begin at rank 2.
//
// Preconditions:
// - global_vertex_ids is non-null
// - IDs are unique over the polygon vertices
inline std::vector<int> build_vertex_ranks(const std::int64_t* global_vertex_ids,
                                           std::size_t n) {
  if (!global_vertex_ids) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global_vertex_ids is null");
  }

  std::vector<std::pair<std::int64_t, std::size_t>> ids;
  ids.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    ids.emplace_back(global_vertex_ids[i], i);
  }

  std::sort(ids.begin(), ids.end(),
            [](const auto& a, const auto& b) {
              if (a.first != b.first) return a.first < b.first;
              return a.second < b.second;
            });

  for (std::size_t i = 1; i < ids.size(); ++i) {
    if (ids[i - 1].first == ids[i].first) {
      throw std::invalid_argument(
          "point_in_polygon_sphere: global_vertex_ids must be unique");
    }
  }

  std::vector<int> ranks(n, -1);
  for (std::size_t rank = 0; rank < n; ++rank) {
    ranks[ids[rank].second] = static_cast<int>(rank + 2);
  }
  return ranks;
}

// Evaluate the symbolic sign prescribed by Simulation of Simplicity for a
// degenerate 3 x 3 determinant whose rows are already ordered by increasing
// symbolic rank:
//
//   rank(a) < rank(b) < rank(c).
//
// This implementation follows Table 14-II of:
//
//   Herbert Edelsbrunner and Ernst Peter Mücke. 1990. Simulation of
//   simplicity: a technique to cope with degenerate cases in geometric
//   algorithms. ACM Trans. Graph. 9, 1 (Jan. 1990), 66–104.
//   https://doi.org/10.1145/77635.77639
//
// The coefficient sequence below is evaluated in decreasing perturbation
// significance. Each 2 x 2 coefficient is reduced to an exact determinant
// sign test, and each 1 x 1 coefficient is reduced to the exact sign of the
// stored scalar.
//
// Precondition:
// - Kernel::orient3d_on_sphere(a, b, c) == 0
inline Sign symbolically_perturbed_sign_sorted(const double* a,
                                               const double* b,
                                               const double* c) {
  Sign s = exact_sign_det2(b[0], b[1], c[0], c[1]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2(b[2], b[0], c[2], c[0]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2(b[1], b[2], c[1], c[2]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2(c[0], c[1], a[0], a[1]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord(c[0]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord(-c[1]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2(c[2], c[0], a[2], a[0]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord(c[2]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2(a[0], a[1], b[0], b[1]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord(-b[0]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord(b[1]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord(a[0]);
  if (s != Sign::Zero) return s;

  return Sign::Positive;
}

// Exact orient3d sign with Simulation of Simplicity fallback.
//
// The exact predicate is evaluated first. If it is nonzero, that sign is
// returned directly. Otherwise, the degenerate case is resolved by symbolic
// perturbation using the supplied distinct symbolic ranks.
inline Sign orient3d_on_sphere_sos(const double* A, int rankA,
                                   const double* B, int rankB,
                                   const double* C, int rankC) {
  const Sign exact = Kernel::orient3d_on_sphere(A, B, C);
  if (exact != Sign::Zero) {
    return exact;
  }

  if (rankA == rankB || rankA == rankC || rankB == rankC) {
    throw std::invalid_argument(
        "orient3d_on_sphere_sos: symbolic ranks must be distinct");
  }

  struct Row {
    const double* p;
    int rank;
    int original_index;
  };

  Row rows[3] = {
      {A, rankA, 0},
      {B, rankB, 1},
      {C, rankC, 2},
  };

  std::sort(std::begin(rows), std::end(rows),
            [](const Row& x, const Row& y) { return x.rank < y.rank; });

  int inversions = 0;
  if (rows[0].original_index > rows[1].original_index) ++inversions;
  if (rows[0].original_index > rows[2].original_index) ++inversions;
  if (rows[1].original_index > rows[2].original_index) ++inversions;

  Sign s = symbolically_perturbed_sign_sorted(rows[0].p, rows[1].p, rows[2].p);
  if (inversions & 1) {
    s = flip_sign(s);
  }
  return s;
}

// SoS version of the ray crossing predicate.
//
// All four orient tests are evaluated with symbolic perturbation. This removes
// the need for an explicitly perturbed antipode and yields a deterministic
// crossing decision from the exact combinatorial ordering of q, R, and the
// polygon vertices.
inline bool counts_as_ray_crossing_sos(const double* A, int rankA,
                                       const double* B, int rankB,
                                       const double* q, int rankQ,
                                       const double* R, int rankR) {
  const Sign s_qR_A = orient3d_on_sphere_sos(q, rankQ, R, rankR, A, rankA);
  const Sign s_qR_B = orient3d_on_sphere_sos(q, rankQ, R, rankR, B, rankB);
  const Sign s_AB_q = orient3d_on_sphere_sos(A, rankA, B, rankB, q, rankQ);
  const Sign s_AB_R = orient3d_on_sphere_sos(A, rankA, B, rankB, R, rankR);

  const bool below_a = (s_qR_A == Sign::Negative);
  const bool below_b = (s_qR_B == Sign::Negative);
  if (!(below_a ^ below_b)) {
    return false;
  }

  const bool q_side = (s_AB_q == Sign::Negative);
  const bool r_side = (s_AB_R == Sign::Negative);
  return (q_side ^ r_side);
}

}  // namespace

Location point_in_polygon_sphere(const double* q,
                                 const double* const* poly,
                                 std::size_t n) {
  std::cout << "\n==== point_in_polygon_sphere debug ====\n";
  std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
std::cout << "q = (" << q[0] << ", " << q[1] << ", " << q[2] << ")\n";

  if (validate_polygon_and_check_vertex(q, poly, n)) {
    std::cout << "result: OnVertex\n";
    return Location::OnVertex;
  }
  std::cout << "vertex check passed\n";

  if (point_on_polygon_edge_exact(q, poly, n)) {
    std::cout << "result: OnEdge\n";
    return Location::OnEdge;
  }
  std::cout << "edge check passed\n";

  double r_arr[3];
  make_perturbed_antipode_simple(q, r_arr);
  const double* R = r_arr;

  std::cout << "R = (" << R[0] << ", " << R[1] << ", " << R[2] << ")\n";

  bool inside = false;
  const double* A = poly[n - 1];
  for (std::size_t i = 0; i < n; ++i) {
    const double* B = poly[i];

    std::cout << "\nedge " << i << '\n';
    std::cout << "A = (" << A[0] << ", " << A[1] << ", " << A[2] << ")\n";
    std::cout << "B = (" << B[0] << ", " << B[1] << ", " << B[2] << ")\n";

    const bool crossing = counts_as_ray_crossing_half_open(A, B, q, R);
    std::cout << "crossing = " << crossing << '\n';

    if (crossing) {
      inside = !inside;
    }

    std::cout << "inside = " << inside << '\n';
    A = B;
  }

  std::cout << "final result = "
            << (inside ? "Inside" : "Outside") << '\n';

  return inside ? Location::Inside : Location::Outside;
}


Location point_in_polygon_sphere(const std::array<double, 3>& q,
                                 const std::vector<std::array<double, 3>>& poly) {
  if (poly.size() < 3) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: polygon must have >= 3 vertices");
  }

  std::vector<const double*> ptrs;
  ptrs.reserve(poly.size());
  for (const auto& v : poly) {
    ptrs.push_back(v.data());
  }

  return point_in_polygon_sphere(q.data(), ptrs.data(), ptrs.size());
}

Location point_in_polygon_sphere(const std::vector<double>& q,
                                 const std::vector<std::vector<double>>& poly) {
  if (q.size() != 3) {
    throw std::invalid_argument("point_in_polygon_sphere: q must have size 3");
  }
  if (poly.size() < 3) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: polygon must have >= 3 vertices");
  }

  std::vector<const double*> ptrs;
  ptrs.reserve(poly.size());
  for (std::size_t i = 0; i < poly.size(); ++i) {
    if (poly[i].size() != 3) {
      throw std::invalid_argument(
          "point_in_polygon_sphere: poly[i] must have size 3");
    }
    ptrs.push_back(poly[i].data());
  }

  return point_in_polygon_sphere(q.data(), ptrs.data(), ptrs.size());
}

Location point_in_polygon_sphere(const double* q,
                                 const double* const* poly,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n) {
  if (validate_polygon_and_check_vertex(q, poly, n)) {
    return Location::OnVertex;
  }

  if (point_on_polygon_edge_exact(q, poly, n)) {
    return Location::OnEdge;
  }

  const std::vector<int> vertex_ranks = build_vertex_ranks(global_vertex_ids, n);

  const double r_arr[3] = {-q[0], -q[1], -q[2]};
  const double* R = r_arr;

  constexpr int rankQ = 0;
  constexpr int rankR = 1;

  bool inside = false;
  const double* A = poly[n - 1];
  int rankA = vertex_ranks[n - 1];

  for (std::size_t i = 0; i < n; ++i) {
    const double* B = poly[i];
    const int rankB = vertex_ranks[i];

    if (counts_as_ray_crossing_sos(A, rankA, B, rankB, q, rankQ, R, rankR)) {
      inside = !inside;
    }

    A = B;
    rankA = rankB;
  }

  return inside ? Location::Inside : Location::Outside;
}

Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    const std::vector<std::array<double, 3>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids) {
  if (poly.size() < 3) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: polygon must have >= 3 vertices");
  }
  if (global_vertex_ids.size() != poly.size()) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global_vertex_ids size must match poly size");
  }

  std::vector<const double*> ptrs;
  ptrs.reserve(poly.size());
  for (const auto& v : poly) {
    ptrs.push_back(v.data());
  }

  return point_in_polygon_sphere(
      q.data(), ptrs.data(), global_vertex_ids.data(), ptrs.size());
}

Location point_in_polygon_sphere(
    const std::vector<double>& q,
    const std::vector<std::vector<double>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids) {
  if (q.size() != 3) {
    throw std::invalid_argument("point_in_polygon_sphere: q must have size 3");
  }
  if (poly.size() < 3) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: polygon must have >= 3 vertices");
  }
  if (global_vertex_ids.size() != poly.size()) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global_vertex_ids size must match poly size");
  }

  std::vector<const double*> ptrs;
  ptrs.reserve(poly.size());
  for (std::size_t i = 0; i < poly.size(); ++i) {
    if (poly[i].size() != 3) {
      throw std::invalid_argument(
          "point_in_polygon_sphere: poly[i] must have size 3");
    }
    ptrs.push_back(poly[i].data());
  }

  return point_in_polygon_sphere(
      q.data(), ptrs.data(), global_vertex_ids.data(), ptrs.size());
}

}  // namespace spip::pip