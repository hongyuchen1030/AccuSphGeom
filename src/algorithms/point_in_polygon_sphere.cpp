#include "accusphgeom/algorithms/point_in_polygon_sphere.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace accusphgeom::algorithms {

namespace {

using predicates::Sign;

void require_nonnull3(const double* p, const char* name) {
  if (!p) {
    throw std::invalid_argument(
        std::string("point_in_polygon_sphere: null pointer for ") + name);
  }
}

bool equal3(const double* a, const double* b) {
  return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
}

bool on_minor_arc(const double* q, const double* a, const double* b) {
  if (equal3(a, b)) {
    return false;
  }
  if (predicates::orient3d_on_sphere(a, b, q) != Sign::Zero) {
    return false;
  }
  return predicates::quadruple3d(a, q, q, b) != Sign::Negative;
}

bool validate_polygon_and_check_vertex(const double* q,
                                       const double* const* polygon,
                                       std::size_t n) {
  require_nonnull3(q, "q");
  if (!polygon) {
    throw std::invalid_argument("point_in_polygon_sphere: polygon is null");
  }
  if (n < 3) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: polygon must have at least 3 vertices");
  }
  for (std::size_t i = 0; i < n; ++i) {
    require_nonnull3(polygon[i], "polygon[i]");
    if (equal3(q, polygon[i])) {
      return true;
    }
  }
  return false;
}

bool point_on_polygon_edge_exact(const double* q, const double* const* polygon,
                                 std::size_t n) {
  const double* a = polygon[n - 1];
  for (std::size_t i = 0; i < n; ++i) {
    const double* b = polygon[i];
    if (on_minor_arc(q, a, b)) {
      return true;
    }
    a = b;
  }
  return false;
}

std::array<double, 3> make_perturbed_antipode(const double* q) {
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
  const double inv_n = 1.0 / std::sqrt(n2);
  r[0] *= inv_n;
  r[1] *= inv_n;
  r[2] *= inv_n;
  return r;
}

[[noreturn]] void throw_ray_endpoint_degeneracy() {
  throw std::domain_error(
      "point_in_polygon_sphere: degenerate ray endpoint because orient(A, B, R) "
      "== 0 for some polygon edge; split the polygon or use the global-ID API");
}

detail::CrossingResult counts_as_ray_crossing(const double* a, const double* b,
                                              const double* q,
                                              const double* r) {
  const Sign s_qr_a = predicates::orient3d_on_sphere(q, r, a);
  const Sign s_qr_b = predicates::orient3d_on_sphere(q, r, b);
  const Sign s_ab_q = predicates::orient3d_on_sphere(a, b, q);
  const Sign s_ab_r = predicates::orient3d_on_sphere(a, b, r);

  if (s_ab_q == Sign::Zero) {
    return detail::CrossingResult::NoCrossing;
  }
  if (s_ab_r == Sign::Zero) {
    throw_ray_endpoint_degeneracy();
  }
  if (s_qr_a != Sign::Zero && s_qr_b != Sign::Zero) {
    const bool crossing = (s_qr_a == s_ab_r) &&
                          (s_qr_a == detail::flip_sign(s_ab_q)) &&
                          (s_qr_a == detail::flip_sign(s_qr_b));
    return crossing ? detail::CrossingResult::Crossing
                    : detail::CrossingResult::NoCrossing;
  }

  const bool below_a = s_qr_a == Sign::Negative;
  const bool below_b = s_qr_b == Sign::Negative;
  if (!(below_a ^ below_b)) {
    return detail::CrossingResult::NoCrossing;
  }

  const bool q_side = s_ab_q == Sign::Negative;
  const bool r_side = s_ab_r == Sign::Negative;
  return (q_side ^ r_side) ? detail::CrossingResult::Crossing
                           : detail::CrossingResult::NoCrossing;
}

detail::CrossingResult counts_as_ray_crossing_sos(
    const double* a, std::int64_t rank_a, const double* b, std::int64_t rank_b,
    const double* q, std::int64_t rank_q, const double* r, std::int64_t rank_r) {
  const Sign s_ab_q = predicates::orient3d_on_sphere(a, b, q);
  if (s_ab_q == Sign::Zero) {
    return detail::CrossingResult::NoCrossing;
  }

  const Sign s_ab_r = predicates::orient3d_on_sphere(a, b, r);
  if (s_ab_r == Sign::Zero) {
    throw_ray_endpoint_degeneracy();
  }

  Sign s_qr_a = predicates::orient3d_on_sphere(q, r, a);
  if (s_qr_a == Sign::Zero) {
    s_qr_a = detail::orient3d_on_sphere_sos_from_doubles(q, rank_q, r, rank_r,
                                                         a, rank_a);
  }

  Sign s_qr_b = predicates::orient3d_on_sphere(q, r, b);
  if (s_qr_b == Sign::Zero) {
    s_qr_b = detail::orient3d_on_sphere_sos_from_doubles(q, rank_q, r, rank_r,
                                                         b, rank_b);
  }

  const bool crossing = (s_qr_a == s_ab_r) &&
                        (s_qr_a == detail::flip_sign(s_ab_q)) &&
                        (s_qr_a == detail::flip_sign(s_qr_b));
  return crossing ? detail::CrossingResult::Crossing
                  : detail::CrossingResult::NoCrossing;
}

detail::FixedRayResult classify_with_fixed_ray(const double* q,
                                               const double* const* polygon,
                                               std::size_t n,
                                               const double* r) {
  bool inside = false;
  const double* a = polygon[n - 1];
  for (std::size_t i = 0; i < n; ++i) {
    const double* b = polygon[i];
    if (counts_as_ray_crossing(a, b, q, r) == detail::CrossingResult::Crossing) {
      inside = !inside;
    }
    a = b;
  }
  return inside ? detail::FixedRayResult::Inside
                : detail::FixedRayResult::Outside;
}

detail::FixedRayResult classify_with_fixed_ray_sos(
    const double* q, const double* const* polygon,
    const detail::SymbolicRanks& ranks, std::size_t n, const double* r) {
  bool inside = false;
  const double* a = polygon[n - 1];
  std::int64_t rank_a = ranks.vertex_ranks[n - 1];
  for (std::size_t i = 0; i < n; ++i) {
    const double* b = polygon[i];
    const std::int64_t rank_b = ranks.vertex_ranks[i];
    if (counts_as_ray_crossing_sos(a, rank_a, b, rank_b, q, ranks.q_rank, r,
                                   ranks.r_rank) ==
        detail::CrossingResult::Crossing) {
      inside = !inside;
    }
    a = b;
    rank_a = rank_b;
  }
  return inside ? detail::FixedRayResult::Inside
                : detail::FixedRayResult::Outside;
}

void require_matching_id_count(std::size_t polygon_size, std::size_t ids_size) {
  if (polygon_size != ids_size) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global_vertex_ids size must match polygon size");
  }
}

}  // namespace

namespace detail {

SymbolicRanks build_symbolic_ranks(const std::int64_t* global_vertex_ids,
                                   std::size_t n, std::int64_t q_id,
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
  ids.emplace_back(q_id, n);
  ids.emplace_back(r_id, n + 1);

  std::sort(ids.begin(), ids.end(),
            [](const auto& lhs, const auto& rhs) {
              if (lhs.first != rhs.first) {
                return lhs.first < rhs.first;
              }
              return lhs.second < rhs.second;
            });

  for (std::size_t i = 1; i < ids.size(); ++i) {
    if (ids[i - 1].first == ids[i].first) {
      throw std::invalid_argument(
          "point_in_polygon_sphere: global IDs must be unique");
    }
  }

  SymbolicRanks out;
  out.vertex_ranks.assign(n, -1);
  for (std::size_t rank = 0; rank < ids.size(); ++rank) {
    const auto slot = ids[rank].second;
    if (slot < n) {
      out.vertex_ranks[slot] = static_cast<std::int64_t>(rank);
    } else if (slot == n) {
      out.q_rank = static_cast<std::int64_t>(rank);
    } else {
      out.r_rank = static_cast<std::int64_t>(rank);
    }
  }
  return out;
}

InternalSymbolicIds assign_internal_symbolic_ids(
    const std::int64_t* global_vertex_ids, std::size_t n, bool assign_q_id) {
  if (!global_vertex_ids) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global_vertex_ids is null");
  }

  std::int64_t max_id = -1;
  std::vector<std::int64_t> ids;
  ids.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto id = global_vertex_ids[i];
    if (id < 0) {
      throw std::invalid_argument(
          "point_in_polygon_sphere: global IDs must be nonnegative");
    }
    max_id = std::max(max_id, id);
    ids.push_back(id);
  }

  InternalSymbolicIds out;
  if (max_id <= std::numeric_limits<std::int64_t>::max() - 2) {
    if (assign_q_id) {
      out.q_id = max_id + 1;
    }
    out.r_id = max_id + 2;
    return out;
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
    out.q_id = first_unused;
  }
  out.r_id = second_unused;
  return out;
}

predicates::Sign symbolically_perturbed_sign_sorted(const double* a,
                                                    const double* b,
                                                    const double* c) {
  predicates::Sign s = exact_sign_det2(b[0], b[1], c[0], c[1]);
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

predicates::Sign orient3d_on_sphere_sos_from_doubles(
    const double* a, std::int64_t rank_a, const double* b, std::int64_t rank_b,
    const double* c, std::int64_t rank_c) {
  const Sign exact = predicates::orient3d_on_sphere(a, b, c);
  if (exact != Sign::Zero) {
    return exact;
  }
  if (rank_a == rank_b || rank_a == rank_c || rank_b == rank_c) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: symbolic ranks must be distinct");
  }

  struct Row {
    const double* p;
    std::int64_t rank;
    int original_index;
  } rows[3] = {{a, rank_a, 0}, {b, rank_b, 1}, {c, rank_c, 2}};

  std::sort(std::begin(rows), std::end(rows),
            [](const Row& lhs, const Row& rhs) { return lhs.rank < rhs.rank; });

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

}  // namespace detail

Location point_in_polygon_sphere(const double* q, const double* const* polygon,
                                 std::size_t n) {
  if (validate_polygon_and_check_vertex(q, polygon, n)) {
    return Location::OnVertex;
  }
  if (point_on_polygon_edge_exact(q, polygon, n)) {
    return Location::OnEdge;
  }
  const auto ray = make_perturbed_antipode(q);
  return classify_with_fixed_ray(q, polygon, n, ray.data()) ==
                 detail::FixedRayResult::Inside
             ? Location::Inside
             : Location::Outside;
}

Location point_in_polygon_sphere(const double* q, const double* const* polygon,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n) {
  const detail::InternalSymbolicIds ids =
      detail::assign_internal_symbolic_ids(global_vertex_ids, n, true);
  return point_in_polygon_sphere(q, ids.q_id, polygon, global_vertex_ids, n);
}

Location point_in_polygon_sphere(const double* q, std::int64_t q_id,
                                 const double* const* polygon,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n) {
  if (validate_polygon_and_check_vertex(q, polygon, n)) {
    return Location::OnVertex;
  }
  if (point_on_polygon_edge_exact(q, polygon, n)) {
    return Location::OnEdge;
  }
  const detail::InternalSymbolicIds ids =
      detail::assign_internal_symbolic_ids(global_vertex_ids, n, false);
  const detail::SymbolicRanks ranks =
      detail::build_symbolic_ranks(global_vertex_ids, n, q_id, ids.r_id);
  const auto ray = make_perturbed_antipode(q);
  return classify_with_fixed_ray_sos(q, polygon, ranks, n, ray.data()) ==
                 detail::FixedRayResult::Inside
             ? Location::Inside
             : Location::Outside;
}

Location point_in_polygon_sphere(const double* q, std::int64_t q_id,
                                 const double* r, std::int64_t r_id,
                                 const double* const* polygon,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n) {
  require_nonnull3(r, "r");
  if (validate_polygon_and_check_vertex(q, polygon, n)) {
    return Location::OnVertex;
  }
  if (point_on_polygon_edge_exact(q, polygon, n)) {
    return Location::OnEdge;
  }
  const detail::SymbolicRanks ranks =
      detail::build_symbolic_ranks(global_vertex_ids, n, q_id, r_id);
  return classify_with_fixed_ray_sos(q, polygon, ranks, n, r) ==
                 detail::FixedRayResult::Inside
             ? Location::Inside
             : Location::Outside;
}

Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    const std::vector<std::array<double, 3>>& polygon) {
  std::vector<const double*> ptrs;
  ptrs.reserve(polygon.size());
  for (const auto& v : polygon) {
    ptrs.push_back(v.data());
  }
  return point_in_polygon_sphere(q.data(), ptrs.data(), ptrs.size());
}

Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    const std::vector<std::array<double, 3>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids) {
  require_matching_id_count(polygon.size(), global_vertex_ids.size());
  std::vector<const double*> ptrs;
  ptrs.reserve(polygon.size());
  for (const auto& v : polygon) {
    ptrs.push_back(v.data());
  }
  return point_in_polygon_sphere(q.data(), ptrs.data(), global_vertex_ids.data(),
                                 ptrs.size());
}

Location point_in_polygon_sphere(
    const std::array<double, 3>& q, std::int64_t q_id,
    const std::vector<std::array<double, 3>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids) {
  require_matching_id_count(polygon.size(), global_vertex_ids.size());
  std::vector<const double*> ptrs;
  ptrs.reserve(polygon.size());
  for (const auto& v : polygon) {
    ptrs.push_back(v.data());
  }
  return point_in_polygon_sphere(q.data(), q_id, ptrs.data(),
                                 global_vertex_ids.data(), ptrs.size());
}

Location point_in_polygon_sphere(
    const std::array<double, 3>& q, std::int64_t q_id,
    const std::array<double, 3>& r, std::int64_t r_id,
    const std::vector<std::array<double, 3>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids) {
  require_matching_id_count(polygon.size(), global_vertex_ids.size());
  std::vector<const double*> ptrs;
  ptrs.reserve(polygon.size());
  for (const auto& v : polygon) {
    ptrs.push_back(v.data());
  }
  return point_in_polygon_sphere(q.data(), q_id, r.data(), r_id, ptrs.data(),
                                 global_vertex_ids.data(), ptrs.size());
}

Location point_in_polygon_sphere(const std::vector<double>& q,
                                 const std::vector<std::vector<double>>& polygon) {
  predicates::detail::require_size3(q.size(), "q");
  std::vector<const double*> ptrs;
  ptrs.reserve(polygon.size());
  for (const auto& v : polygon) {
    predicates::detail::require_size3(v.size(), "polygon[i]");
    ptrs.push_back(v.data());
  }
  return point_in_polygon_sphere(q.data(), ptrs.data(), ptrs.size());
}

Location point_in_polygon_sphere(
    const std::vector<double>& q, const std::vector<std::vector<double>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids) {
  predicates::detail::require_size3(q.size(), "q");
  require_matching_id_count(polygon.size(), global_vertex_ids.size());
  std::vector<const double*> ptrs;
  ptrs.reserve(polygon.size());
  for (const auto& v : polygon) {
    predicates::detail::require_size3(v.size(), "polygon[i]");
    ptrs.push_back(v.data());
  }
  return point_in_polygon_sphere(q.data(), ptrs.data(), global_vertex_ids.data(),
                                 ptrs.size());
}

Location point_in_polygon_sphere(
    const std::vector<double>& q, std::int64_t q_id,
    const std::vector<std::vector<double>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids) {
  predicates::detail::require_size3(q.size(), "q");
  require_matching_id_count(polygon.size(), global_vertex_ids.size());
  std::vector<const double*> ptrs;
  ptrs.reserve(polygon.size());
  for (const auto& v : polygon) {
    predicates::detail::require_size3(v.size(), "polygon[i]");
    ptrs.push_back(v.data());
  }
  return point_in_polygon_sphere(q.data(), q_id, ptrs.data(),
                                 global_vertex_ids.data(), ptrs.size());
}

Location point_in_polygon_sphere(
    const std::vector<double>& q, std::int64_t q_id, const std::vector<double>& r,
    std::int64_t r_id, const std::vector<std::vector<double>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids) {
  predicates::detail::require_size3(q.size(), "q");
  predicates::detail::require_size3(r.size(), "r");
  require_matching_id_count(polygon.size(), global_vertex_ids.size());
  std::vector<const double*> ptrs;
  ptrs.reserve(polygon.size());
  for (const auto& v : polygon) {
    predicates::detail::require_size3(v.size(), "polygon[i]");
    ptrs.push_back(v.data());
  }
  return point_in_polygon_sphere(q.data(), q_id, r.data(), r_id, ptrs.data(),
                                 global_vertex_ids.data(), ptrs.size());
}

}  // namespace accusphgeom::algorithms
