#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
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

Location point_in_polygon_sphere(const double* q, std::int64_t q_id,
                                 const double* const* polygon,
                                 const std::int64_t* global_vertex_ids,
                                 std::size_t n);

Location point_in_polygon_sphere(const double* q, std::int64_t q_id,
                                 const double* r, std::int64_t r_id,
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
    const std::array<double, 3>& q, std::int64_t q_id,
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
    const std::array<double, 3>& q, std::int64_t q_id,
    const std::vector<std::array<double, 3>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids);

template <std::size_t N>
inline Location point_in_polygon_sphere(
    const std::array<double, 3>& q, std::int64_t q_id,
    const std::array<double, 3>& r, std::int64_t r_id,
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
    const std::array<double, 3>& q, std::int64_t q_id,
    const std::array<double, 3>& r, std::int64_t r_id,
    const std::vector<std::array<double, 3>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(const std::vector<double>& q,
                                 const std::vector<std::vector<double>>& polygon);

Location point_in_polygon_sphere(
    const std::vector<double>& q, const std::vector<std::vector<double>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(
    const std::vector<double>& q, std::int64_t q_id,
    const std::vector<std::vector<double>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids);

Location point_in_polygon_sphere(
    const std::vector<double>& q, std::int64_t q_id, const std::vector<double>& r,
    std::int64_t r_id, const std::vector<std::vector<double>>& polygon,
    const std::vector<std::int64_t>& global_vertex_ids);

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

inline predicates::Sign flip_sign(predicates::Sign sign) {
  if (sign == predicates::Sign::Positive) {
    return predicates::Sign::Negative;
  }
  if (sign == predicates::Sign::Negative) {
    return predicates::Sign::Positive;
  }
  return predicates::Sign::Zero;
}

inline predicates::Sign exact_sign_coord(double x) {
  if (x > 0.0) {
    return predicates::Sign::Positive;
  }
  if (x < 0.0) {
    return predicates::Sign::Negative;
  }
  return predicates::Sign::Zero;
}

inline predicates::Sign exact_sign_det2(double a00, double a01, double a10,
                                        double a11) {
  const double r0[3] = {a00, a01, 0.0};
  const double r1[3] = {a10, a11, 0.0};
  const double r2[3] = {0.0, 0.0, 1.0};
  return predicates::orient3d_on_sphere(r0, r1, r2);
}

SymbolicRanks build_symbolic_ranks(const std::int64_t* global_vertex_ids,
                                   std::size_t n, std::int64_t q_id,
                                   std::int64_t r_id);

InternalSymbolicIds assign_internal_symbolic_ids(
    const std::int64_t* global_vertex_ids, std::size_t n, bool assign_q_id);

predicates::Sign orient3d_on_sphere_sos_from_doubles(
    const double* a, std::int64_t rank_a, const double* b, std::int64_t rank_b,
    const double* c, std::int64_t rank_c);

}  // namespace detail

}  // namespace accusphgeom::algorithms
