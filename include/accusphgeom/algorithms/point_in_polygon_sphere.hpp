#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "accusphgeom/predicates/orient3d.hpp"
#include "accusphgeom/predicates/quadruple3d.hpp"

namespace accusphgeom::algorithms {

enum class Location : std::uint8_t { Outside, Inside, OnVertex, OnEdge };

Location point_in_polygon_sphere(const double* q,
                                 const double* const* polygon,
                                 std::size_t n);

Location point_in_polygon_sphere(const double* q,
                                 const double* const* polygon,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

Location point_in_polygon_sphere(const double* q,
                                 std::int64_t q_id,
                                 const double* const* polygon,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

Location point_in_polygon_sphere(const double* q,
                                 std::int64_t q_id,
                                 const double* r,
                                 std::int64_t r_id,
                                 const double* const* polygon,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

template <std::size_t N>
inline Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    const std::array<std::array<double, 3>, N>& polygon) {
  std::array<const double*, N> ptrs{};
  for (std::size_t i = 0; i < N; ++i) {
    ptrs[i] = polygon[i].data();
  }
  return point_in_polygon_sphere(q.data(), ptrs.data(), N);
}

Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    const std::vector<std::array<double, 3>>& polygon);

template <std::size_t N>
inline Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    const std::array<std::array<double, 3>, N>& polygon,
    const std::array<std::int64_t, N>& global_vertex_ids) {
  std::array<const double*, N> ptrs{};
  for (std::size_t i = 0; i < N; ++i) {
    ptrs[i] = polygon[i].data();
  }
  return point_in_polygon_sphere(q.data(), ptrs.data(), global_vertex_ids.data(),
                                 N);
}

Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    const std::vector<std::array<double, 3>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids);

template <std::size_t N>
inline Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    std::int64_t q_id,
    const std::array<std::array<double, 3>, N>& polygon,
    const std::array<std::int64_t, N>& global_vertex_ids) {
  std::array<const double*, N> ptrs{};
  for (std::size_t i = 0; i < N; ++i) {
    ptrs[i] = polygon[i].data();
  }
  return point_in_polygon_sphere(q.data(), q_id, ptrs.data(),
                                 global_vertex_ids.data(), N);
}

Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    std::int64_t q_id,
    const std::vector<std::array<double, 3>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids);

template <std::size_t N>
inline Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    std::int64_t q_id,
    const std::array<double, 3>& r,
    std::int64_t r_id,
    const std::array<std::array<double, 3>, N>& polygon,
    const std::array<std::int64_t, N>& global_vertex_ids) {
  std::array<const double*, N> ptrs{};
  for (std::size_t i = 0; i < N; ++i) {
    ptrs[i] = polygon[i].data();
  }
  return point_in_polygon_sphere(q.data(), q_id, r.data(), r_id, ptrs.data(),
                                 global_vertex_ids.data(), N);
}

Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    std::int64_t q_id,
    const std::array<double, 3>& r,
    std::int64_t r_id,
    const std::vector<std::array<double, 3>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(
    const std::vector<double>& q,
    const std::vector<std::vector<double>>& polygon);

Location point_in_polygon_sphere(
    const std::vector<double>& q,
    const std::vector<std::vector<double>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(
    const std::vector<double>& q,
    std::int64_t q_id,
    const std::vector<std::vector<double>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(
    const std::vector<double>& q,
    std::int64_t q_id,
    const std::vector<double>& r,
    std::int64_t r_id,
    const std::vector<std::vector<double>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids);

namespace detail {

template <typename GlobalId>
using EnableIfSupportedGlobalId = std::enable_if_t<
    std::is_integral_v<GlobalId> &&
    !std::is_same_v<std::remove_cv_t<GlobalId>, std::int64_t>,
    int>;

template <typename GlobalId, EnableIfSupportedGlobalId<GlobalId> = 0>
inline std::int64_t convert_global_id(GlobalId id) {
  if constexpr (std::numeric_limits<GlobalId>::digits >
                std::numeric_limits<std::int64_t>::digits) {
    if (id < static_cast<GlobalId>(std::numeric_limits<std::int64_t>::min()) ||
        id > static_cast<GlobalId>(std::numeric_limits<std::int64_t>::max())) {
      throw std::invalid_argument(
          "point_in_polygon_sphere: global IDs must fit in int64_t");
    }
  }
  return static_cast<std::int64_t>(id);
}

template <typename GlobalId,
          std::size_t N,
          EnableIfSupportedGlobalId<GlobalId> = 0>
inline std::array<std::int64_t, N> convert_global_vertex_ids(
    const std::array<GlobalId, N>& global_vertex_ids) {
  std::array<std::int64_t, N> converted{};
  for (std::size_t i = 0; i < N; ++i) {
    converted[i] = convert_global_id(global_vertex_ids[i]);
  }
  return converted;
}

template <typename GlobalId, EnableIfSupportedGlobalId<GlobalId> = 0>
inline std::vector<std::int64_t> convert_global_vertex_ids(
    const std::vector<GlobalId>& global_vertex_ids) {
  std::vector<std::int64_t> converted;
  converted.reserve(global_vertex_ids.size());
  for (const GlobalId id : global_vertex_ids) {
    converted.push_back(convert_global_id(id));
  }
  return converted;
}

template <typename GlobalId, EnableIfSupportedGlobalId<GlobalId> = 0>
inline std::vector<std::int64_t> convert_global_vertex_ids(
    const GlobalId* global_vertex_ids,
    std::size_t n) {
  if (!global_vertex_ids) {
    throw std::invalid_argument(
        "point_in_polygon_sphere: global_vertex_ids is null");
  }
  std::vector<std::int64_t> converted;
  converted.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    converted.push_back(convert_global_id(global_vertex_ids[i]));
  }
  return converted;
}

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

inline predicates::Sign flip_sign(predicates::Sign sign) {
  if (sign == predicates::Sign::Positive) {
    return predicates::Sign::Negative;
  }
  if (sign == predicates::Sign::Negative) {
    return predicates::Sign::Positive;
  }
  return predicates::Sign::Zero;
}

inline predicates::Sign exact_sign_coord_sos(double x) {
  if (x > 0.0) {
    return predicates::Sign::Positive;
  }
  if (x < 0.0) {
    return predicates::Sign::Negative;
  }
  return predicates::Sign::Zero;
}

inline predicates::Sign exact_sign_det2_sos(double a00, double a01, double a10,
                                            double a11) {
  const double r0[3] = {a00, a01, 0.0};
  const double r1[3] = {a10, a11, 0.0};
  const double r2[3] = {0.0, 0.0, 1.0};
  return predicates::orient3d_on_sphere(r0, r1, r2);
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

inline predicates::Sign symbolically_perturbed_sign_sorted(
    const double* a,
    const double* b,
    const double* c) {
  using Sign = predicates::Sign;

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

inline predicates::Sign orient3d_on_sphere_sos_from_doubles(
    const double* a,
    std::int64_t rank_a,
    const double* b,
    std::int64_t rank_b,
    const double* c,
    std::int64_t rank_c) {
  using Sign = predicates::Sign;

  const Sign exact = predicates::orient3d_on_sphere(a, b, c);
  if (exact != Sign::Zero) {
    return exact;
  }

  if (rank_a == rank_b || rank_a == rank_c || rank_b == rank_c) {
    throw std::invalid_argument(
        "orient3d_on_sphere_sos: symbolic ranks must be distinct");
  }

  struct Row {
    const double* point;
    std::int64_t rank;
    int original_index;
  };

  Row rows[3] = {
      {a, rank_a, 0},
      {b, rank_b, 1},
      {c, rank_c, 2},
  };

  std::sort(std::begin(rows), std::end(rows),
            [](const Row& lhs, const Row& rhs) {
              return lhs.rank < rhs.rank;
            });

  int inversions = 0;
  if (rows[0].original_index > rows[1].original_index) ++inversions;
  if (rows[0].original_index > rows[2].original_index) ++inversions;
  if (rows[1].original_index > rows[2].original_index) ++inversions;

  Sign sos = symbolically_perturbed_sign_sorted(
      rows[0].point, rows[1].point, rows[2].point);
  if (inversions & 1) {
    sos = flip_sign(sos);
  }
  return sos;
}

}  // namespace detail

template <typename GlobalId, detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(const double* q,
                                        const double* const* polygon,
                                        const GlobalId* global_vertex_ids,
                                        std::size_t n) {
  const std::vector<std::int64_t> converted =
      detail::convert_global_vertex_ids(global_vertex_ids, n);
  return point_in_polygon_sphere(q, polygon, converted.data(), n);
}

template <typename GlobalId, detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(const double* q,
                                        std::int64_t q_id,
                                        const double* const* polygon,
                                        const GlobalId* global_vertex_ids,
                                        std::size_t n) {
  const std::vector<std::int64_t> converted =
      detail::convert_global_vertex_ids(global_vertex_ids, n);
  return point_in_polygon_sphere(q, q_id, polygon, converted.data(), n);
}

template <typename GlobalId, detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(const double* q,
                                        std::int64_t q_id,
                                        const double* r,
                                        std::int64_t r_id,
                                        const double* const* polygon,
                                        const GlobalId* global_vertex_ids,
                                        std::size_t n) {
  const std::vector<std::int64_t> converted =
      detail::convert_global_vertex_ids(global_vertex_ids, n);
  return point_in_polygon_sphere(
      q, q_id, r, r_id, polygon, converted.data(), n);
}

template <std::size_t N,
          typename GlobalId,
          detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    const std::array<std::array<double, 3>, N>& polygon,
    const std::array<GlobalId, N>& global_vertex_ids) {
  const std::array<std::int64_t, N> converted =
      detail::convert_global_vertex_ids(global_vertex_ids);
  return point_in_polygon_sphere(q, polygon, converted);
}

template <typename GlobalId, detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    const std::vector<std::array<double, 3>>& polygon,
    const std::vector<GlobalId>& global_vertex_ids) {
  const std::vector<std::int64_t> converted =
      detail::convert_global_vertex_ids(global_vertex_ids);
  return point_in_polygon_sphere(q, polygon, converted);
}

template <std::size_t N,
          typename GlobalId,
          detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    std::int64_t q_id,
    const std::array<std::array<double, 3>, N>& polygon,
    const std::array<GlobalId, N>& global_vertex_ids) {
  const std::array<std::int64_t, N> converted =
      detail::convert_global_vertex_ids(global_vertex_ids);
  return point_in_polygon_sphere(q, q_id, polygon, converted);
}

template <typename GlobalId, detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    std::int64_t q_id,
    const std::vector<std::array<double, 3>>& polygon,
    const std::vector<GlobalId>& global_vertex_ids) {
  const std::vector<std::int64_t> converted =
      detail::convert_global_vertex_ids(global_vertex_ids);
  return point_in_polygon_sphere(q, q_id, polygon, converted);
}

template <std::size_t N,
          typename GlobalId,
          detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    std::int64_t q_id,
    const std::array<double, 3>& r,
    std::int64_t r_id,
    const std::array<std::array<double, 3>, N>& polygon,
    const std::array<GlobalId, N>& global_vertex_ids) {
  const std::array<std::int64_t, N> converted =
      detail::convert_global_vertex_ids(global_vertex_ids);
  return point_in_polygon_sphere(q, q_id, r, r_id, polygon, converted);
}

template <typename GlobalId, detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(
    const std::array<double, 3>& q,
    std::int64_t q_id,
    const std::array<double, 3>& r,
    std::int64_t r_id,
    const std::vector<std::array<double, 3>>& polygon,
    const std::vector<GlobalId>& global_vertex_ids) {
  const std::vector<std::int64_t> converted =
      detail::convert_global_vertex_ids(global_vertex_ids);
  return point_in_polygon_sphere(q, q_id, r, r_id, polygon, converted);
}

template <typename GlobalId, detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(
    const std::vector<double>& q,
    const std::vector<std::vector<double>>& polygon,
    const std::vector<GlobalId>& global_vertex_ids) {
  const std::vector<std::int64_t> converted =
      detail::convert_global_vertex_ids(global_vertex_ids);
  return point_in_polygon_sphere(q, polygon, converted);
}

template <typename GlobalId, detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(
    const std::vector<double>& q,
    std::int64_t q_id,
    const std::vector<std::vector<double>>& polygon,
    const std::vector<GlobalId>& global_vertex_ids) {
  const std::vector<std::int64_t> converted =
      detail::convert_global_vertex_ids(global_vertex_ids);
  return point_in_polygon_sphere(q, q_id, polygon, converted);
}

template <typename GlobalId, detail::EnableIfSupportedGlobalId<GlobalId> = 0>
inline Location point_in_polygon_sphere(
    const std::vector<double>& q,
    std::int64_t q_id,
    const std::vector<double>& r,
    std::int64_t r_id,
    const std::vector<std::vector<double>>& polygon,
    const std::vector<GlobalId>& global_vertex_ids) {
  const std::vector<std::int64_t> converted =
      detail::convert_global_vertex_ids(global_vertex_ids);
  return point_in_polygon_sphere(q, q_id, r, r_id, polygon, converted);
}

}  // namespace accusphgeom::algorithms
