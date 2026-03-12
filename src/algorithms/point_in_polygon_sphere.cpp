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

#include <spip/predicates/orient3d.hpp>
#include <spip/predicates/quadruple3d.hpp>

namespace spip::pip {

namespace {

using spip::predicates::orient3d_on_sphere;
using spip::predicates::quadruple3d;
using spip::predicates::Sign;

inline void require_nonnull3(const double* p, const char* name) {
  if (!p) {
    throw std::invalid_argument(
        std::string("point_in_polygon_sphere: null pointer for ") + name);
  }
}

inline bool equal3(const double* a, const double* b) {
  return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

inline Sign flip_sign(Sign s) {
  if (s == Sign::Positive) return Sign::Negative;
  if (s == Sign::Negative) return Sign::Positive;
  return Sign::Zero;
}

inline Sign exact_sign_coord(double x) {
  if (x > 0.0) return Sign::Positive;
  if (x < 0.0) return Sign::Negative;
  return Sign::Zero;
}

// Reuse the existing exact orient3d predicate to evaluate the sign of a 2x2
// determinant:
//
//   | a00 a01 |
//   | a10 a11 |
//
// by embedding it as a 3x3 determinant:
//
//   | a00 a01 0 |
//   | a10 a11 0 |
//   |  0   0  1 |
inline Sign exact_sign_det2(double a00, double a01,
                            double a10, double a11) {
  const double r0[3] = {a00, a01, 0.0};
  const double r1[3] = {a10, a11, 0.0};
  const double r2[3] = {0.0, 0.0, 1.0};
  return orient3d_on_sphere(r0, r1, r2);
}

// q lies on the interior/boundary of the non-antipodal minor arc AB iff:
//   (1) q lies on the supporting great circle of AB:
//         orient(A,B,0,q) == 0
//   (2) the sub-arc normals A×q and q×B do not point in opposite directions:
//         (A×q)·(q×B) >= 0
//
// The second condition is the scalar quadruple product:
//   quadruple3d(A, q, q, B) = (A×q)·(q×B)
//
// Preconditions:
// - q, A, B are valid non-null pointers to 3 doubles
// - vertex coincidence q==A or q==B has already been handled by caller
inline bool on_minor_arc(const double* q, const double* A, const double* B) {
  if (equal3(A, B)) {
    return false;
  }
  if (orient3d_on_sphere(A, B, q) != Sign::Zero) {
    return false;
  }
  return quadruple3d(A, q, q, B) != Sign::Negative;
}

// Shared validation for raw-pointer polygon input and exact vertex test.
// Returns true iff q exactly matches a polygon vertex.
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

// Shared exact boundary test over all edges. Assumes exact OnVertex has already
// been excluded.
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

// Ray endpoint: antipode with deterministic perturbation.
// This is used only by the non-SoS overload.
inline void make_perturbed_antipode_simple(const double* q, double* r_out) {
  r_out[0] = -q[0];
  r_out[1] = -q[1];
  r_out[2] = -q[2];

  const double ax = std::fabs(r_out[0]);
  const double ay = std::fabs(r_out[1]);
  const double az = std::fabs(r_out[2]);

  constexpr double eps = 1e-15;
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

// We consider the ray minor arc from query point q to ray endpoint R
// (perturbed antipode). For a polygon edge minor arc AB, this predicate
// returns true iff AB contributes one crossing to the ray crossing count,
// using the half-open convention of Hormann and Agathos (2001).
inline bool counts_as_ray_crossing_half_open(const double* A,
                                             const double* B,
                                             const double* q,
                                             const double* R) {
  const Sign s_qR_A = orient3d_on_sphere(q, R, A);
  const Sign s_qR_B = orient3d_on_sphere(q, R, B);
  const Sign s_AB_q = orient3d_on_sphere(A, B, q);
  const Sign s_AB_R = orient3d_on_sphere(A, B, R);

  // Step 1: half-open straddle test about the ray great circle.
  // below(x) := (x < 0). Zero is treated as non-negative.
  const bool below_a = (s_qR_A == Sign::Negative);
  const bool below_b = (s_qR_B == Sign::Negative);
  if (!(below_a ^ below_b)) {
    return false;
  }

  // Step 2: in-front test about the edge great circle.
  // After exact boundary prechecks, these should not be Zero.
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

  const bool q_side = (s_AB_q == Sign::Negative);
  const bool r_side = (s_AB_R == Sign::Negative);
  return (q_side ^ r_side);
}

// Build a strict symbolic order from global vertex IDs.
// Smaller rank means larger symbolic perturbation priority.
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

  // Reserve ranks 0 and 1 for q and R.
  std::vector<int> ranks(n, -1);
  for (std::size_t rank = 0; rank < n; ++rank) {
    ranks[ids[rank].second] = static_cast<int>(rank + 2);
  }
  return ranks;
}

// S2-style symbolic perturbation sign for a 3x3 determinant, with rows
// already sorted by increasing symbolic rank: rank(a) < rank(b) < rank(c).
// Precondition: orient3d_on_sphere(a,b,c) == 0.
inline Sign symbolically_perturbed_sign_sorted(const double* a,
                                               const double* b,
                                               const double* c) {
  // Coefficients are tested in decreasing perturbation significance.
  // 2x2 coefficients are evaluated through the existing exact orient3d
  // predicate, while 1x1 coefficients use the exact sign of the stored value.

  // da[2], da[1], da[0] from b x c
  Sign s = exact_sign_det2(b[0], b[1], c[0], c[1]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2(b[2], b[0], c[2], c[0]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2(b[1], b[2], c[1], c[2]);
  if (s != Sign::Zero) return s;

  // db[2], db[2]*da[1], db[2]*da[0], db[1], db[1]*da[0]
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

  // db[0] is redundant here, same as in S2.

  // dc[2], dc[2]*da[1], dc[2]*da[0], dc[2]*db[1]
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

inline Sign orient3d_on_sphere_sos(const double* A, int rankA,
                                   const double* B, int rankB,
                                   const double* C, int rankC) {
  const Sign exact = orient3d_on_sphere(A, B, C);
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
  if (validate_polygon_and_check_vertex(q, poly, n)) {
    return Location::OnVertex;
  }

  if (point_on_polygon_edge_exact(q, poly, n)) {
    return Location::OnEdge;
  }

  double r_arr[3];
  make_perturbed_antipode_simple(q, r_arr);
  const double* R = r_arr;

  bool inside = false;
  const double* A = poly[n - 1];
  for (std::size_t i = 0; i < n; ++i) {
    const double* B = poly[i];

    if (counts_as_ray_crossing_half_open(A, B, q, R)) {
      inside = !inside;
    }

    A = B;
  }

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