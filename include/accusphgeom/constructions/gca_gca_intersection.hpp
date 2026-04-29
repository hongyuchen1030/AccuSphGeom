#pragma once
#include <type_traits>
#include <array>
#include <cmath>
#include <stdexcept>

#include "accusphgeom/constructions/accucross.hpp"
#include "accusphgeom/predicates/on_minor_arc.hpp"

namespace accusphgeom::constructions {

namespace internal {

inline constexpr double gca_gca_minor_arc_tol = 1e-8;

}  // namespace internal

template <typename T>
struct GcaGcaIntersections {
  numeric::Vec3<T> point_pos{};
  numeric::Vec3<T> point_neg{};
};

template <typename T, typename StatusT>
struct GcaGcaTryResult {
  numeric::Vec3<T> point{};
  StatusT status{};  // 0 ok, 1 both valid, 2 none valid
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
inline auto try_gca_gca_intersection(const numeric::Vec3<T>& a0,
                                     const numeric::Vec3<T>& a1,
                                     const numeric::Vec3<T>& b0,
                                     const numeric::Vec3<T>& b1) {
  const auto candidates = accux_gca(a0, a1, b0, b1);
  const T tol = T(internal::gca_gca_minor_arc_tol);
  const T one = T(1);

  const auto pos_finite =
      numeric::isfinite_mask(candidates.point_pos[0]) *
      numeric::isfinite_mask(candidates.point_pos[1]) *
      numeric::isfinite_mask(candidates.point_pos[2]);

  const auto neg_finite =
      numeric::isfinite_mask(candidates.point_neg[0]) *
      numeric::isfinite_mask(candidates.point_neg[1]) *
      numeric::isfinite_mask(candidates.point_neg[2]);

  const auto pos_on_arc_a =
      pos_finite *
      predicates::internal::on_minor_arc_tol_ptr(candidates.point_pos.data(),
                                                 a0.data(), a1.data(), tol);
  const auto pos_on_arc_b =
      pos_finite *
      predicates::internal::on_minor_arc_tol_ptr(candidates.point_pos.data(),
                                                 b0.data(), b1.data(), tol);

  const auto neg_on_arc_a =
      neg_finite *
      predicates::internal::on_minor_arc_tol_ptr(candidates.point_neg.data(),
                                                 a0.data(), a1.data(), tol);
  const auto neg_on_arc_b =
      neg_finite *
      predicates::internal::on_minor_arc_tol_ptr(candidates.point_neg.data(),
                                                 b0.data(), b1.data(), tol);

  const auto pos_valid = pos_finite * pos_on_arc_a * pos_on_arc_b;
  const auto neg_valid = neg_finite * neg_on_arc_a * neg_on_arc_b;

  const auto pos_mask = pos_valid * (one - neg_valid);
  const auto neg_mask = neg_valid * (one - pos_valid);

  numeric::Vec3<T> out{};
  out[0] = pos_mask * candidates.point_pos[0] +
           neg_mask * candidates.point_neg[0];
  out[1] = pos_mask * candidates.point_pos[1] +
           neg_mask * candidates.point_neg[1];
  out[2] = pos_mask * candidates.point_pos[2] +
           neg_mask * candidates.point_neg[2];

  const auto both = pos_valid * neg_valid;
  const auto none = (one - pos_valid) * (one - neg_valid);
  const auto status = both + none * T(2);

  using StatusT = decltype(status);
  return GcaGcaTryResult<T, StatusT>{out, status};
}

template <typename T>
inline numeric::Vec3<T> gca_gca_intersection(const numeric::Vec3<T>& a0,
                                             const numeric::Vec3<T>& a1,
                                             const numeric::Vec3<T>& b0,
                                             const numeric::Vec3<T>& b1) {
  static_assert(std::is_arithmetic_v<T>,
                "gca_gca_intersection is scalar-only; use "
                "try_gca_gca_intersection for packed/SIMD types.");

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
