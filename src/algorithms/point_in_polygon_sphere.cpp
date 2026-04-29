#include <accusphgeom/algorithms/point_in_polygon_sphere.hpp>
#include <accusphgeom/predicates/on_minor_arc.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace accusphgeom::algorithms {

namespace {

  // Helper: Convert a container of 3D points to std::vector<const double*>
  template <typename Container>
  std::vector<const double*> to_ptr_vector(const Container& poly) {
    std::vector<const double*> ptrs;
    ptrs.reserve(poly.size());
    for (const auto& v : poly) {
      ptrs.push_back(v.data());
    }
    return ptrs;
  }

  // Overload for std::vector<std::vector<double>>: check size==3
  inline std::vector<const double*> to_ptr_vector(
  const std::vector<std::vector<double>>& poly) {
    std::vector<const double*> ptrs;
    ptrs.reserve(poly.size());
    for (const auto& v : poly) {
      if (v.size() != 3) {
        throw std::invalid_argument(
          "point_in_polygon_sphere: poly[i] must have size 3");
      }
      ptrs.push_back(v.data());
    }
    return ptrs;
  }

// Helper: Validate polygon, ID, and query/endpoint sizes
//         for all container types
inline void validate_poly_and_query_sizes(
    std::size_t poly_size,
    std::size_t id_size = SIZE_MAX,
    std::size_t q_size = SIZE_MAX,
    std::size_t r_size = SIZE_MAX) {
  if (poly_size < 3) {
    throw std::invalid_argument(
      "point_in_polygon_sphere: polygon must have >= 3 vertices");
  }
  if (id_size != SIZE_MAX && poly_size != id_size) {
    throw std::invalid_argument(
      "point_in_polygon_sphere: global_vertex_ids size must match poly size");
  }
  if (q_size != SIZE_MAX && q_size != 3) {
    throw std::invalid_argument(
      "point_in_polygon_sphere: q must have size 3");
  }
  if (r_size != SIZE_MAX && r_size != 3) {
    throw std::invalid_argument(
      "point_in_polygon_sphere: R must have size 3");
  }
}

using Sign = accusphgeom::predicates::Sign;

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
    if (accusphgeom::predicates::on_minor_arc(q, A, B)) {
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
inline std::array<double, 3> make_perturbed_antipode(const double* q) {
  std::array<double, 3> r = {-q[0], -q[1], -q[2]};

  const double ax = std::fabs(r[0]);
  const double ay = std::fabs(r[1]);
  const double az = std::fabs(r[2]);

  constexpr double eps = 1e-8;
  if (ax <= ay && ax <= az) {
    r[0] += eps;
  } else if (ay <= ax && ay <= az) {
    r[1] += eps;
  } else {
    r[2] += eps;
  }

  const double n2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
  if (n2 > 0.0) {
    const double inv_n = 1.0 / std::sqrt(n2);
    r[0] *= inv_n;
    r[1] *= inv_n;
    r[2] *= inv_n;
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
inline detail::CrossingResult counts_as_ray_crossing(
    const double* A, const double* B, const double* q, const double* R) {
  const Sign s_qR_A = accusphgeom::predicates::orient3d_on_sphere(q, R, A);
  const Sign s_qR_B = accusphgeom::predicates::orient3d_on_sphere(q, R, B);
  const Sign s_AB_q = accusphgeom::predicates::orient3d_on_sphere(A, B, q);
  const Sign s_AB_R = accusphgeom::predicates::orient3d_on_sphere(A, B, R);

  if (s_AB_q == Sign::Zero) {
    return detail::CrossingResult::NoCrossing;
  }
  if (s_AB_R == Sign::Zero) {
    throw_ray_endpoint_degeneracy();
  }

  if (s_qR_A != Sign::Zero && s_qR_B != Sign::Zero) {
    const bool crossing =
        (s_qR_A == s_AB_R) &&
        (s_qR_A == detail::flip_sign(s_AB_q)) &&
        (s_qR_A == detail::flip_sign(s_qR_B));
    return crossing ? detail::CrossingResult::Crossing
                    : detail::CrossingResult::NoCrossing;
  }

  const bool below_a = (s_qR_A == Sign::Negative);
  const bool below_b = (s_qR_B == Sign::Negative);
  if (!(below_a ^ below_b)) {
    return detail::CrossingResult::NoCrossing;
  }

  const bool q_side = (s_AB_q == Sign::Negative);
  const bool r_side = (s_AB_R == Sign::Negative);
  return (q_side ^ r_side) ? detail::CrossingResult::Crossing
                           : detail::CrossingResult::NoCrossing;
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
inline detail::CrossingResult counts_as_ray_crossing_sos(
    const double* A, std::int64_t rankA,
    const double* B, std::int64_t rankB,
    const double* q, std::int64_t rankQ,
    const double* R, std::int64_t rankR) {
  const Sign s_AB_q = accusphgeom::predicates::orient3d_on_sphere(A, B, q);
  if (s_AB_q == Sign::Zero) {
    return detail::CrossingResult::NoCrossing;
  }

  const Sign s_AB_R = accusphgeom::predicates::orient3d_on_sphere(A, B, R);
  if (s_AB_R == Sign::Zero) {
    throw_ray_endpoint_degeneracy();
  }

  Sign s_qR_A = accusphgeom::predicates::orient3d_on_sphere(q, R, A);
  if (s_qR_A == Sign::Zero) {
    s_qR_A = detail::orient3d_on_sphere_sos_from_doubles(
        q, rankQ, R, rankR, A, rankA);
  }

  Sign s_qR_B = accusphgeom::predicates::orient3d_on_sphere(q, R, B);
  if (s_qR_B == Sign::Zero) {
    s_qR_B = detail::orient3d_on_sphere_sos_from_doubles(
        q, rankQ, R, rankR, B, rankB);
  }

  const bool crossing =
      (s_qR_A == s_AB_R) &&
      (s_qR_A == detail::flip_sign(s_AB_q)) &&
      (s_qR_A == detail::flip_sign(s_qR_B));
  return crossing ? detail::CrossingResult::Crossing
                  : detail::CrossingResult::NoCrossing;
}

inline detail::FixedRayResult classify_with_fixed_ray(
    const double* q, const double* const* poly, std::size_t n, const double* R) {
  bool inside = false;
  const double* A = poly[n - 1];
  for (std::size_t i = 0; i < n; ++i) {
    const double* B = poly[i];
    const detail::CrossingResult crossing = counts_as_ray_crossing(A, B, q, R);
    if (crossing == detail::CrossingResult::Crossing) {
      inside = !inside;
    }
    A = B;
  }
  return inside ? detail::FixedRayResult::Inside
                : detail::FixedRayResult::Outside;
}

inline detail::FixedRayResult classify_with_fixed_ray_sos(
    const double* q,
    const double* const* poly,
    const detail::SymbolicRanks& ranks,
    std::size_t n,
    const double* R) {
  bool inside = false;
  const double* A = poly[n - 1];
  std::int64_t rankA = ranks.vertex_ranks[n - 1];

  for (std::size_t i = 0; i < n; ++i) {
    const double* B = poly[i];
    const std::int64_t rankB = ranks.vertex_ranks[i];

    const detail::CrossingResult crossing =
        counts_as_ray_crossing_sos(
            A, rankA, B, rankB, q, ranks.q_rank, R, ranks.r_rank);
    if (crossing == detail::CrossingResult::Crossing) {
      inside = !inside;
    }

    A = B;
    rankA = rankB;
  }

  return inside ? detail::FixedRayResult::Inside
                : detail::FixedRayResult::Outside;
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

  const std::array<double, 3> primary_ray = make_perturbed_antipode(q);
  const detail::FixedRayResult primary =
      classify_with_fixed_ray(q, poly, n, primary_ray.data());
  return primary == detail::FixedRayResult::Inside
             ? Location::Inside
             : Location::Outside;
}

Location point_in_polygon_sphere(const std::array<double, 3>& q,
                                 const std::vector<std::array<double, 3>>& poly) {
  validate_poly_and_query_sizes(poly.size(), SIZE_MAX, 3);
  auto ptrs = to_ptr_vector(poly);
  return point_in_polygon_sphere(q.data(), ptrs.data(), ptrs.size());
}

Location point_in_polygon_sphere(const std::vector<double>& q,
                                 const std::vector<std::vector<double>>& poly) {
  validate_poly_and_query_sizes(poly.size(), SIZE_MAX, q.size());
  auto ptrs = to_ptr_vector(poly);
  return point_in_polygon_sphere(q.data(), ptrs.data(), ptrs.size());
}

Location detail::point_in_polygon_sphere_impl(
    const double* q,
    const double* const* poly,
    const std::int64_t* compact_vertex_ids,
    std::size_t n) {
  // Tier 3: vertices have user IDs; q and the inferred perturbed-antipode
  // endpoint R receive internal IDs. The induced SoS ordering is therefore
  // deterministic per call, but not a fully user-controlled global ordering.
  const detail::InternalSymbolicIds assigned =
      detail::assign_internal_symbolic_ids(compact_vertex_ids, n, true);
  return detail::point_in_polygon_sphere_impl(q, assigned.q_id, poly,
                                              compact_vertex_ids, n);
}

Location detail::point_in_polygon_sphere_impl(
    const double* q,
    std::int64_t q_id,
    const double* const* poly,
    const std::int64_t* compact_vertex_ids,
    std::size_t n) {
  // Tier 2: vertices and q have user IDs; only the inferred ray endpoint R is
  // assigned internally.
  if (validate_polygon_and_check_vertex(q, poly, n)) {
    return Location::OnVertex;
  }

  if (point_on_polygon_edge_exact(q, poly, n)) {
    return Location::OnEdge;
  }

  const detail::InternalSymbolicIds assigned =
      detail::assign_internal_symbolic_ids(compact_vertex_ids, n, false);
  const detail::SymbolicRanks ranks =
      detail::build_symbolic_ranks(
          compact_vertex_ids, n, q_id, assigned.r_id);

  const std::array<double, 3> primary_ray = make_perturbed_antipode(q);
  const detail::FixedRayResult primary =
      classify_with_fixed_ray_sos(q, poly, ranks, n, primary_ray.data());
  return primary == detail::FixedRayResult::Inside
             ? Location::Inside
             : Location::Outside;
}

Location detail::point_in_polygon_sphere_impl(
    const double* q,
    std::int64_t q_id,
    const double* R,
    std::int64_t r_id,
    const double* const* poly,
    const std::int64_t* compact_vertex_ids,
    std::size_t n) {
  // Tier 1: q, R, and every vertex participate in one caller-defined symbolic
  // ordering.
  require_nonnull3(R, "R");

  if (validate_polygon_and_check_vertex(q, poly, n)) {
    return Location::OnVertex;
  }

  if (point_on_polygon_edge_exact(q, poly, n)) {
    return Location::OnEdge;
  }

  const detail::SymbolicRanks ranks =
      detail::build_symbolic_ranks(compact_vertex_ids, n, q_id, r_id);

  const detail::FixedRayResult primary =
      classify_with_fixed_ray_sos(q, poly, ranks, n, R);
  return primary == detail::FixedRayResult::Inside
             ? Location::Inside
             : Location::Outside;
}

}  // namespace accusphgeom::algorithms
