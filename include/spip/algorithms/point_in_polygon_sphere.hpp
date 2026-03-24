#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "spip/core/types.hpp"
#include "spip/kernels/pip_kernel_eft.hpp"

namespace spip::pip {

enum class Location : std::uint8_t { Outside, Inside, OnVertex, OnEdge };

// No global IDs: use the perturbed-antipode ray and exact boundary checks.
Location point_in_polygon_sphere(const double* q,
                                 const double* const* poly,
                                 std::size_t n);

// Tier 3: vertex IDs are user-provided; q and the inferred perturbed-antipode
// endpoint R receive internal IDs. The resulting SoS ordering is deterministic
// within the call, but it is only a local symbolic hierarchy.
Location point_in_polygon_sphere(const double* q,
                                 const double* const* poly,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

// Tier 2: vertex IDs and q_id are user-provided; the inferred perturbed-
// antipode endpoint R receives an internal ID.
Location point_in_polygon_sphere(const double* q,
                                 std::int64_t q_id,
                                 const double* const* poly,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

// Tier 1: vertex IDs, q_id, the outside point R, and r_id are all
// user-provided, so the symbolic ordering is fully controlled by the caller.
Location point_in_polygon_sphere(const double* q,
                                 std::int64_t q_id,
                                 const double* R,
                                 std::int64_t r_id,
                                 const double* const* poly,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

Location point_in_polygon_sphere(const std::array<double, 3>& q,
                                 const std::vector<std::array<double, 3>>& poly);

Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    const std::vector<std::array<double, 3>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    std::int64_t q_id,
    const std::vector<std::array<double, 3>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    std::int64_t q_id,
    const std::array<double, 3>& R,
    std::int64_t r_id,
    const std::vector<std::array<double, 3>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(const std::vector<double>& q,
                                 const std::vector<std::vector<double>>& poly);

Location point_in_polygon_sphere(
    const std::vector<double>& q,
    const std::vector<std::vector<double>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(
    const std::vector<double>& q,
    std::int64_t q_id,
    const std::vector<std::vector<double>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(
    const std::vector<double>& q,
    std::int64_t q_id,
    const std::vector<double>& R,
    std::int64_t r_id,
    const std::vector<std::vector<double>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids);

// EFT overloads. These are templated and therefore implemented in this header.
template <typename T>
Location point_in_polygon_sphere(const V3_T<T>& q,
                                 const V3_T<T>* poly,
                                 std::size_t n);

template <typename T>
Location point_in_polygon_sphere(const V3_T<T>& q,
                                 const V3_T<T>* poly,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

template <typename T>
Location point_in_polygon_sphere(const V3_T<T>& q,
                                 std::int64_t q_id,
                                 const V3_T<T>* poly,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

template <typename T>
Location point_in_polygon_sphere(const V3_T<T>& q,
                                 std::int64_t q_id,
                                 const V3_T<T>& R,
                                 std::int64_t r_id,
                                 const V3_T<T>* poly,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

template <typename T>
inline Location point_in_polygon_sphere(const V3_T<T>& q,
                                        const std::vector<V3_T<T>>& poly) {
  return point_in_polygon_sphere(q, poly.data(), poly.size());
}

template <typename T>
inline Location point_in_polygon_sphere(
    const V3_T<T>& q,
    const std::vector<V3_T<T>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids) {
  if (global_vertex_ids.size() != poly.size()) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global_vertex_ids size must match poly size");
  }
  return point_in_polygon_sphere(
      q, poly.data(), global_vertex_ids.data(), poly.size());
}

template <typename T>
inline Location point_in_polygon_sphere(
    const V3_T<T>& q,
    std::int64_t q_id,
    const std::vector<V3_T<T>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids) {
  if (global_vertex_ids.size() != poly.size()) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global_vertex_ids size must match poly size");
  }
  return point_in_polygon_sphere(
      q, q_id, poly.data(), global_vertex_ids.data(), poly.size());
}

template <typename T>
inline Location point_in_polygon_sphere(
    const V3_T<T>& q,
    std::int64_t q_id,
    const V3_T<T>& R,
    std::int64_t r_id,
    const std::vector<V3_T<T>>& poly,
    const std::vector<std::int64_t>& global_vertex_ids) {
  if (global_vertex_ids.size() != poly.size()) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global_vertex_ids size must match poly size");
  }
  return point_in_polygon_sphere(
      q, q_id, R, r_id, poly.data(), global_vertex_ids.data(), poly.size());
}

template <typename T, std::size_t N>
inline Location point_in_polygon_sphere(const V3_T<T>& q,
                                        const std::array<V3_T<T>, N>& poly) {
  return point_in_polygon_sphere(q, poly.data(), poly.size());
}

template <typename T, std::size_t N>
inline Location point_in_polygon_sphere(
    const V3_T<T>& q,
    const std::array<V3_T<T>, N>& poly,
    const std::array<std::int64_t, N>& global_vertex_ids) {
  return point_in_polygon_sphere(
      q, poly.data(), global_vertex_ids.data(), poly.size());
}

template <typename T, std::size_t N>
inline Location point_in_polygon_sphere(
    const V3_T<T>& q,
    std::int64_t q_id,
    const std::array<V3_T<T>, N>& poly,
    const std::array<std::int64_t, N>& global_vertex_ids) {
  return point_in_polygon_sphere(
      q, q_id, poly.data(), global_vertex_ids.data(), poly.size());
}

template <typename T, std::size_t N>
inline Location point_in_polygon_sphere(
    const V3_T<T>& q,
    std::int64_t q_id,
    const V3_T<T>& R,
    std::int64_t r_id,
    const std::array<V3_T<T>, N>& poly,
    const std::array<std::int64_t, N>& global_vertex_ids) {
  return point_in_polygon_sphere(
      q, q_id, R, r_id, poly.data(), global_vertex_ids.data(), poly.size());
}

namespace detail {

enum class CrossingResult : std::uint8_t { NoCrossing, Crossing };
enum class FixedRayResult : std::uint8_t { Outside, Inside };

struct SymbolicRanks {
  std::vector<std::int64_t> vertex_ranks;
  std::int64_t q_rank = -1;
  std::int64_t r_rank = -1;
};

struct InternalSymbolicIds {
  std::int64_t q_id = -1;
  std::int64_t r_id = -1;
};

template <typename T>
inline void require_valid_polygon(const V3_T<T>* poly, std::size_t n) {
  if (poly == nullptr) {
    throw std::invalid_argument("point_in_polygon_sphere: poly is null");
  }
  if (n < 3) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: polygon must have >= 3 vertices");
  }
}

template <typename T>
inline bool equal3(const V3_T<T>& a, const V3_T<T>& b) {
  return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

inline spip::predicates::Sign flip_sign(spip::predicates::Sign s) {
  using Sign = spip::predicates::Sign;
  if (s == Sign::Positive) return Sign::Negative;
  if (s == Sign::Negative) return Sign::Positive;
  return Sign::Zero;
}

inline spip::predicates::Sign exact_sign_coord_sos(double x) {
  using Sign = spip::predicates::Sign;
  if (x > 0.0) return Sign::Positive;
  if (x < 0.0) return Sign::Negative;
  return Sign::Zero;
}

inline spip::predicates::Sign exact_sign_det2_sos(double a00, double a01,
                                                  double a10, double a11) {
  const double r0[3] = {a00, a01, 0.0};
  const double r1[3] = {a10, a11, 0.0};
  const double r2[3] = {0.0, 0.0, 1.0};
  return spip::predicates::orient3d_on_sphere(r0, r1, r2);
}

inline SymbolicRanks build_symbolic_ranks(const std::int64_t* global_vertex_ids,
                                          std::size_t n,
                                          std::int64_t q_id,
                                          std::int64_t r_id) {
  if (!global_vertex_ids) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global_vertex_ids is null");
  }
  if (q_id < 0 || r_id < 0) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global IDs must be nonnegative");
  }

  std::vector<std::pair<std::int64_t, std::size_t>> ids;
  ids.reserve(n + 2);
  for (std::size_t i = 0; i < n; ++i) {
    if (global_vertex_ids[i] < 0) {
      throw std::invalid_argument(
          "point_in_polygon_sphere: global IDs must be nonnegative");
    }
    ids.emplace_back(global_vertex_ids[i], i);
  }
  const std::size_t q_slot = n;
  const std::size_t r_slot = n + 1;
  ids.emplace_back(q_id, q_slot);
  ids.emplace_back(r_id, r_slot);

  std::sort(ids.begin(), ids.end(),
            [](const auto& a, const auto& b) {
              if (a.first != b.first) return a.first < b.first;
              return a.second < b.second;
            });

  for (std::size_t i = 1; i < ids.size(); ++i) {
    if (ids[i - 1].first == ids[i].first) {
      throw std::invalid_argument(
          "point_in_polygon_sphere: global IDs must be unique");
    }
  }

  SymbolicRanks ranks;
  ranks.vertex_ranks.assign(n, -1);
  for (std::size_t rank = 0; rank < ids.size(); ++rank) {
    const std::size_t slot = ids[rank].second;
    const std::int64_t compact_rank = static_cast<std::int64_t>(rank);
    if (slot < n) {
      ranks.vertex_ranks[slot] = compact_rank;
    } else if (slot == q_slot) {
      ranks.q_rank = compact_rank;
    } else {
      ranks.r_rank = compact_rank;
    }
  }
  return ranks;
}

// Assign internal IDs for Tier 2 or Tier 3 without modifying vertex IDs.
//
// First try max_id + {1,2}. If that would overflow, fall back to a MEX scan
// over the nonnegative integers induced by the vertex IDs alone.
inline InternalSymbolicIds assign_internal_symbolic_ids(
    const std::int64_t* global_vertex_ids,
    std::size_t n,
    bool assign_q_id) {
  if (!global_vertex_ids) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global_vertex_ids is null");
  }

  std::int64_t max_id = -1;
  std::vector<std::int64_t> ids;
  ids.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::int64_t id = global_vertex_ids[i];
    if (id < 0) {
      throw std::invalid_argument(
          "point_in_polygon_sphere: global IDs must be nonnegative");
    }
    if (id > max_id) {
      max_id = id;
    }
    ids.push_back(id);
  }

  InternalSymbolicIds assigned;
  if (max_id <= std::numeric_limits<std::int64_t>::max() - 2) {
    if (assign_q_id) {
      assigned.q_id = max_id + 1;
    }
    assigned.r_id = max_id + 2;
    return assigned;
  }

  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

  std::int64_t next = 0;
  std::int64_t first_unused = -1;
  std::int64_t second_unused = -1;
  std::size_t i = 0;
  while (second_unused < 0) {
    while (i < ids.size() && ids[i] < next) {
      ++i;
    }
    if (i < ids.size() && ids[i] == next) {
      ++next;
      ++i;
      continue;
    }
    if (first_unused < 0) {
      first_unused = next;
    } else {
      second_unused = next;
    }
    ++next;
  }

  if (assign_q_id) {
    assigned.q_id = first_unused;
  }
  assigned.r_id = second_unused;
  return assigned;
}

inline spip::predicates::Sign symbolically_perturbed_sign_sorted(
    const double* a, const double* b, const double* c) {
  using Sign = spip::predicates::Sign;

  Sign s = exact_sign_det2_sos(b[0], b[1], c[0], c[1]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2_sos(b[2], b[0], c[2], c[0]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2_sos(b[1], b[2], c[1], c[2]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2_sos(c[0], c[1], a[0], a[1]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord_sos(c[0]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord_sos(-c[1]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2_sos(c[2], c[0], a[2], a[0]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord_sos(c[2]);
  if (s != Sign::Zero) return s;

  s = exact_sign_det2_sos(a[0], a[1], b[0], b[1]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord_sos(-b[0]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord_sos(b[1]);
  if (s != Sign::Zero) return s;

  s = exact_sign_coord_sos(a[0]);
  if (s != Sign::Zero) return s;

  return Sign::Positive;
}

inline spip::predicates::Sign orient3d_on_sphere_sos_from_doubles(
    const double* A, std::int64_t rankA,
    const double* B, std::int64_t rankB,
    const double* C, std::int64_t rankC) {
  using Sign = spip::predicates::Sign;

  const Sign exact = spip::predicates::orient3d_on_sphere(A, B, C);
  if (exact != Sign::Zero) {
    return exact;
  }

  if (rankA == rankB || rankA == rankC || rankB == rankC) {
    throw std::invalid_argument(
        "orient3d_on_sphere_sos: symbolic ranks must be distinct");
  }

  struct Row {
    const double* p;
    std::int64_t rank;
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

template <typename T>
inline std::array<double, 3> to_double_array(const V3_T<T>& v) {
  return {static_cast<double>(v[0]),
          static_cast<double>(v[1]),
          static_cast<double>(v[2])};
}

template <typename T>
inline spip::predicates::Sign orient3d_on_sphere_sos_eft(
    const V3_T<T>& A, std::int64_t rankA,
    const V3_T<T>& B, std::int64_t rankB,
    const V3_T<T>& C, std::int64_t rankC) {
  using Sign = spip::predicates::Sign;

  const Sign approx = spip::kernels::PIPKernelEFT::orient3d_on_sphere(A, B, C);
  if (approx != Sign::Zero) {
    return approx;
  }

  const std::array<double, 3> Ad = to_double_array(A);
  const std::array<double, 3> Bd = to_double_array(B);
  const std::array<double, 3> Cd = to_double_array(C);
  return orient3d_on_sphere_sos_from_doubles(
      Ad.data(), rankA, Bd.data(), rankB, Cd.data(), rankC);
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
template <typename T>
inline bool on_minor_arc_eft(const V3_T<T>& q,
                             const V3_T<T>& A,
                             const V3_T<T>& B) {
  using Kernel = spip::kernels::PIPKernelEFT;
  using Sign = spip::predicates::Sign;

  if (equal3(A, B)) {
    return false;
  }
  if (Kernel::orient3d_on_sphere(A, B, q) != Sign::Zero) {
    return false;
  }
  return Kernel::quadruple3d(A, q, q, B) != Sign::Negative;
}

template <typename T>
inline bool point_on_polygon_edge_eft(const V3_T<T>& q,
                                      const V3_T<T>* poly,
                                      std::size_t n) {
  const V3_T<T>* A = &poly[n - 1];
  for (std::size_t i = 0; i < n; ++i) {
    const V3_T<T>* B = &poly[i];
    if (on_minor_arc_eft(q, *A, *B)) {
      return true;
    }
    A = B;
  }
  return false;
}

// Construct a deterministic perturbation of the antipode of q.
//
// This endpoint is the unique ray endpoint used by this algorithm. If some
// polygon edge satisfies orient(A, B, R) = 0, the ray construction is
// degenerate and the query is rejected.
template <typename T>
inline V3_T<T> make_perturbed_antipode_eft(const V3_T<T>& q) {
  V3_T<T> r = -q;

  const T ax = std::abs(r[0]);
  const T ay = std::abs(r[1]);
  const T az = std::abs(r[2]);

  const T eps = T(1e-15);
  if (ax <= ay && ax <= az) {
    r[0] += eps;
  } else if (ay <= ax && ay <= az) {
    r[1] += eps;
  } else {
    r[2] += eps;
  }

  const T n2 = r.dot(r);
  if (n2 > T(0)) {
    r /= std::sqrt(n2);
  }
  return r;
}

[[noreturn]] inline void throw_ray_endpoint_degeneracy() {
  throw std::domain_error(
      "point_in_polygon_sphere: degenerate ray endpoint because orient(A, B, R) "
      "== 0 for some polygon edge; the polygon is too large or close to "
      "hemispherical, the perturbed-antipode ray construction becomes "
      "degenerate, and you should split the face into smaller pieces or use "
      "the global-ID API so Simulation of Simplicity can resolve degeneracies");
}

// Return true iff the polygon edge AB contributes one crossing to the ray
// crossing count for the minor arc qR.
//
// The decision uses four predicate signs:
//
//   s_qR_A = orient(q, R, A)
//   s_qR_B = orient(q, R, B)
//   s_AB_q = orient(A, B, q)
//   s_AB_R = orient(A, B, R)
//
// For the non-global-ID path:
// - s_AB_q == 0 means not crossing
// - s_AB_R == 0 means the perturbed-antipode ray is invalid and an exception is
//   thrown
// - if s_qR_A and s_qR_B are both nonzero, use the strict 4-sign theorem
// - otherwise use the half-open branch structure
template <typename T>
inline CrossingResult counts_as_ray_crossing_eft(const V3_T<T>& A,
                                                 const V3_T<T>& B,
                                                 const V3_T<T>& q,
                                                 const V3_T<T>& R) {
  using Kernel = spip::kernels::PIPKernelEFT;
  using Sign = spip::predicates::Sign;

  const Sign s_qR_A = Kernel::orient3d_on_sphere(q, R, A);
  const Sign s_qR_B = Kernel::orient3d_on_sphere(q, R, B);
  const Sign s_AB_q = Kernel::orient3d_on_sphere(A, B, q);
  const Sign s_AB_R = Kernel::orient3d_on_sphere(A, B, R);

  if (s_AB_q == Sign::Zero) {
    return CrossingResult::NoCrossing;
  }
  if (s_AB_R == Sign::Zero) {
    throw_ray_endpoint_degeneracy();
  }

  if (s_qR_A != Sign::Zero && s_qR_B != Sign::Zero) {
    const bool crossing =
        (s_qR_A == s_AB_R) &&
        (s_qR_A == flip_sign(s_AB_q)) &&
        (s_qR_A == flip_sign(s_qR_B));
    return crossing ? CrossingResult::Crossing : CrossingResult::NoCrossing;
  }

  const bool below_a = (s_qR_A == Sign::Negative);
  const bool below_b = (s_qR_B == Sign::Negative);
  if (!(below_a ^ below_b)) {
    return CrossingResult::NoCrossing;
  }

  const bool q_side = (s_AB_q == Sign::Negative);
  const bool r_side = (s_AB_R == Sign::Negative);
  return (q_side ^ r_side) ? CrossingResult::Crossing
                           : CrossingResult::NoCrossing;
}

// Return true iff the polygon edge AB contributes one crossing to the ray
// crossing count for the minor arc qR in the global-ID path.
//
// The symbolic ordering is derived from the participating global IDs:
// vertices plus user-provided or internally assigned IDs for q and R.
//
// The exact orient signs s_AB_q and s_AB_R are never perturbed:
// - s_AB_q == 0 means not crossing
// - s_AB_R == 0 means the perturbed-antipode ray is invalid and an exception is
//   thrown
//
// The only SoS-resolved degeneracies are the endpoint tests s_qR_A and s_qR_B.
// After resolving those signs, the decision uses the same strict 4-sign theorem
// as the nondegenerate non-global-ID path.
template <typename T>
inline CrossingResult counts_as_ray_crossing_sos_eft(
    const V3_T<T>& A, std::int64_t rankA,
    const V3_T<T>& B, std::int64_t rankB,
    const V3_T<T>& q, std::int64_t rankQ,
    const V3_T<T>& R, std::int64_t rankR) {
  using Kernel = spip::kernels::PIPKernelEFT;
  using Sign = spip::predicates::Sign;

  const Sign s_AB_q = Kernel::orient3d_on_sphere(A, B, q);
  if (s_AB_q == Sign::Zero) {
    return CrossingResult::NoCrossing;
  }

  const Sign s_AB_R = Kernel::orient3d_on_sphere(A, B, R);
  if (s_AB_R == Sign::Zero) {
    throw_ray_endpoint_degeneracy();
  }

  Sign s_qR_A = Kernel::orient3d_on_sphere(q, R, A);
  if (s_qR_A == Sign::Zero) {
    s_qR_A = orient3d_on_sphere_sos_eft(q, rankQ, R, rankR, A, rankA);
  }

  Sign s_qR_B = Kernel::orient3d_on_sphere(q, R, B);
  if (s_qR_B == Sign::Zero) {
    s_qR_B = orient3d_on_sphere_sos_eft(q, rankQ, R, rankR, B, rankB);
  }

  const bool crossing =
      (s_qR_A == s_AB_R) &&
      (s_qR_A == flip_sign(s_AB_q)) &&
      (s_qR_A == flip_sign(s_qR_B));
  return crossing ? CrossingResult::Crossing : CrossingResult::NoCrossing;
}

template <typename T>
inline FixedRayResult classify_with_fixed_ray_eft(const V3_T<T>& q,
                                                  const V3_T<T>* poly,
                                                  std::size_t n,
                                                  const V3_T<T>& R) {
  bool inside = false;
  const V3_T<T>* A = &poly[n - 1];
  for (std::size_t i = 0; i < n; ++i) {
    const V3_T<T>* B = &poly[i];
    const CrossingResult crossing = counts_as_ray_crossing_eft(*A, *B, q, R);
    if (crossing == CrossingResult::Crossing) {
      inside = !inside;
    }
    A = B;
  }

  return inside ? FixedRayResult::Inside : FixedRayResult::Outside;
}

template <typename T>
inline FixedRayResult classify_with_fixed_ray_sos_eft(
    const V3_T<T>& q,
    const V3_T<T>* poly,
    const SymbolicRanks& ranks,
    std::size_t n,
    const V3_T<T>& R) {
  bool inside = false;
  const V3_T<T>* A = &poly[n - 1];
  std::int64_t rankA = ranks.vertex_ranks[n - 1];

  for (std::size_t i = 0; i < n; ++i) {
    const V3_T<T>* B = &poly[i];
    const std::int64_t rankB = ranks.vertex_ranks[i];

    const CrossingResult crossing =
        counts_as_ray_crossing_sos_eft(
            *A, rankA, *B, rankB, q, ranks.q_rank, R, ranks.r_rank);
    if (crossing == CrossingResult::Crossing) {
      inside = !inside;
    }

    A = B;
    rankA = rankB;
  }

  return inside ? FixedRayResult::Inside : FixedRayResult::Outside;
}

}  // namespace detail

template <typename T>
inline Location point_in_polygon_sphere(const V3_T<T>& q,
                                        const V3_T<T>* poly,
                                        std::size_t n) {
  detail::require_valid_polygon(poly, n);

  for (std::size_t i = 0; i < n; ++i) {
    if (detail::equal3(q, poly[i])) {
      return Location::OnVertex;
    }
  }

  if (detail::point_on_polygon_edge_eft(q, poly, n)) {
    return Location::OnEdge;
  }

  const detail::FixedRayResult primary =
      detail::classify_with_fixed_ray_eft(
          q, poly, n, detail::make_perturbed_antipode_eft(q));
  return primary == detail::FixedRayResult::Inside
             ? Location::Inside
             : Location::Outside;
}

template <typename T>
inline Location point_in_polygon_sphere(const V3_T<T>& q,
                                        const V3_T<T>* poly,
                                        const std::int64_t* global_vertex_ids,
                                        std::size_t n) {
  const detail::InternalSymbolicIds assigned =
      detail::assign_internal_symbolic_ids(global_vertex_ids, n, true);
  return point_in_polygon_sphere(
      q, assigned.q_id, poly, global_vertex_ids, n);
}

template <typename T>
inline Location point_in_polygon_sphere(const V3_T<T>& q,
                                        std::int64_t q_id,
                                        const V3_T<T>* poly,
                                        const std::int64_t* global_vertex_ids,
                                        std::size_t n) {
  detail::require_valid_polygon(poly, n);

  for (std::size_t i = 0; i < n; ++i) {
    if (detail::equal3(q, poly[i])) {
      return Location::OnVertex;
    }
  }

  if (detail::point_on_polygon_edge_eft(q, poly, n)) {
    return Location::OnEdge;
  }

  const detail::InternalSymbolicIds assigned =
      detail::assign_internal_symbolic_ids(global_vertex_ids, n, false);
  const detail::SymbolicRanks ranks =
      detail::build_symbolic_ranks(
          global_vertex_ids, n, q_id, assigned.r_id);

  const detail::FixedRayResult primary =
      detail::classify_with_fixed_ray_sos_eft(
          q, poly, ranks, n, detail::make_perturbed_antipode_eft(q));
  return primary == detail::FixedRayResult::Inside
             ? Location::Inside
             : Location::Outside;
}

template <typename T>
inline Location point_in_polygon_sphere(const V3_T<T>& q,
                                        std::int64_t q_id,
                                        const V3_T<T>& R,
                                        std::int64_t r_id,
                                        const V3_T<T>* poly,
                                        const std::int64_t* global_vertex_ids,
                                        std::size_t n) {
  detail::require_valid_polygon(poly, n);

  for (std::size_t i = 0; i < n; ++i) {
    if (detail::equal3(q, poly[i])) {
      return Location::OnVertex;
    }
  }

  if (detail::point_on_polygon_edge_eft(q, poly, n)) {
    return Location::OnEdge;
  }

  const detail::SymbolicRanks ranks =
      detail::build_symbolic_ranks(global_vertex_ids, n, q_id, r_id);

  const detail::FixedRayResult primary =
      detail::classify_with_fixed_ray_sos_eft(q, poly, ranks, n, R);
  return primary == detail::FixedRayResult::Inside
             ? Location::Inside
             : Location::Outside;
}

}  // namespace spip::pip
