#pragma once

#include <array>
#include <cmath>
#include <stdexcept>

#include "accusphgeom/constructions/accucross.hpp"
#include "accusphgeom/predicates/on_minor_arc.hpp"

namespace accusphgeom::constructions {

namespace internal {

template <typename T>
inline constexpr T gca_gca_minor_arc_tol =
    static_cast<T>(1e-8);

}  // namespace internal

template <typename T>
struct GcaGcaIntersections {
  numeric::Vec3<T> point_pos{};
  numeric::Vec3<T> point_neg{};
};

template <typename T>
struct GcaGcaTryResult {
  numeric::Vec3<T> point{};
  int status{};  // 0 ok, 1 both valid, 2 none valid
};

template <typename T>
inline GcaGcaIntersections<T> accux_gca(const numeric::Vec3<T>& a0,
                                        const numeric::Vec3<T>& a1,
                                        const numeric::Vec3<T>& b0,
                                        const numeric::Vec3<T>& b1) {
  const auto n1 = accucross(a0, a1);
  const auto n2 = accucross(b0, b1);
  const auto v = accucross(n1.hi, n1.lo, n2.hi, n2.lo);
  const auto sum = numeric::sum_of_squares_c<T, 3>(v.hi, v.lo);
  const auto norm = numeric::acc_sqrt_re(sum.hi, sum.lo);
  const T n = norm.hi + norm.lo;

  const numeric::Vec3<T> point_pos = {(v.hi[0] + v.lo[0]) / n,
                                      (v.hi[1] + v.lo[1]) / n,
                                      (v.hi[2] + v.lo[2]) / n};
  const numeric::Vec3<T> point_neg = {-point_pos[0], -point_pos[1],
                                      -point_pos[2]};
  return {point_pos, point_neg};
}

template <typename T>
inline GcaGcaTryResult<T> try_gca_gca_intersection(
    const numeric::Vec3<T>& a0,
    const numeric::Vec3<T>& a1,
    const numeric::Vec3<T>& b0,
    const numeric::Vec3<T>& b1) {
  const auto candidates = accux_gca(a0, a1, b0, b1);
  constexpr T minor_arc_tolerance =
      internal::gca_gca_minor_arc_tol<T>;
  const bool pos_finite =
      std::isfinite(candidates.point_pos[0]) &
      std::isfinite(candidates.point_pos[1]) &
      std::isfinite(candidates.point_pos[2]);
  const bool neg_finite =
      std::isfinite(candidates.point_neg[0]) &
      std::isfinite(candidates.point_neg[1]) &
      std::isfinite(candidates.point_neg[2]);
  const bool pos_on_arc_a =
      predicates::on_minor_arc_tol_ptr(candidates.point_pos, a0, a1,
                                       minor_arc_tolerance);
  const bool pos_on_arc_b =
      predicates::on_minor_arc_tol_ptr(candidates.point_pos, b0, b1,
                                       minor_arc_tolerance);
  const bool neg_on_arc_a =
      predicates::on_minor_arc_tol_ptr(candidates.point_neg, a0, a1,
                                       minor_arc_tolerance);
  const bool neg_on_arc_b =
      predicates::on_minor_arc_tol_ptr(candidates.point_neg, b0, b1,
                                       minor_arc_tolerance);
  const int pos_valid = pos_finite & pos_on_arc_a & pos_on_arc_b;
  const int neg_valid = neg_finite & neg_on_arc_a & neg_on_arc_b;

  const T pos_mask = static_cast<T>(pos_valid & !neg_valid);
  const T neg_mask = static_cast<T>(neg_valid & !pos_valid);

  numeric::Vec3<T> out{};
  out[0] = pos_mask * candidates.point_pos[0] +
           neg_mask * candidates.point_neg[0];
  out[1] = pos_mask * candidates.point_pos[1] +
           neg_mask * candidates.point_neg[1];
  out[2] = pos_mask * candidates.point_pos[2] +
           neg_mask * candidates.point_neg[2];

  const int both = pos_valid & neg_valid;
  const int none = (!pos_valid) & (!neg_valid);
  const int status = both + 2 * none;

  return {out, status};
}

template <typename T>
inline numeric::Vec3<T> gca_gca_intersection(const numeric::Vec3<T>& a0,
                                             const numeric::Vec3<T>& a1,
                                             const numeric::Vec3<T>& b0,
                                             const numeric::Vec3<T>& b1) {
  const auto result = try_gca_gca_intersection(a0, a1, b0, b1);

  if (result.status == 0) return result.point;

  if (result.status == 1) {
    throw std::domain_error(
        "gca_gca_intersection: both intersections lie on both minor arcs");
  }
  throw std::domain_error(
      "gca_gca_intersection: no intersection lies on both minor arcs");
}

}  // namespace accusphgeom::constructions
