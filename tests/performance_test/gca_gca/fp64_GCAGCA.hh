#pragma once

#include <cmath>

#include "accusphgeom/adapters/eigen/numeric.hpp"
#include "accusphgeom/constructions/gca_gca_intersection.hpp"

namespace accusphgeom::performance_test {

template <typename T>
inline accusphgeom::numeric::Vec3<T> fp64_cross(
    const accusphgeom::numeric::Vec3<T>& a,
    const accusphgeom::numeric::Vec3<T>& b) {
  return {
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0],
  };
}

inline accusphgeom::numeric::Vec3<double> fp64_normalize(
    const accusphgeom::numeric::Vec3<double>& v) {
  const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  return {v[0] / n, v[1] / n, v[2] / n};
}

template <int N>
inline accusphgeom::numeric::Vec3<EigenPack<N>> fp64_normalize(
    const accusphgeom::numeric::Vec3<EigenPack<N>>& v) {
  const EigenPack<N> n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
  return {v[0] / n, v[1] / n, v[2] / n};
}

inline accusphgeom::constructions::GcaGcaIntersections<double> fp64_gca_gca(
    const accusphgeom::numeric::Vec3<double>& a0,
    const accusphgeom::numeric::Vec3<double>& a1,
    const accusphgeom::numeric::Vec3<double>& b0,
    const accusphgeom::numeric::Vec3<double>& b1) {
  const auto n1 = fp64_cross(a0, a1);
  const auto n2 = fp64_cross(b0, b1);
  const auto v = fp64_normalize(fp64_cross(n1, n2));

  const accusphgeom::numeric::Vec3<double> point_pos = v;
  const accusphgeom::numeric::Vec3<double> point_neg = {-v[0], -v[1], -v[2]};
  return {point_pos, point_neg};
}

template <int N>
inline accusphgeom::constructions::GcaGcaIntersections<EigenPack<N>>
fp64_gca_gca(const accusphgeom::numeric::Vec3<EigenPack<N>>& a0,
             const accusphgeom::numeric::Vec3<EigenPack<N>>& a1,
             const accusphgeom::numeric::Vec3<EigenPack<N>>& b0,
             const accusphgeom::numeric::Vec3<EigenPack<N>>& b1) {
  const auto n1 = fp64_cross(a0, a1);
  const auto n2 = fp64_cross(b0, b1);
  const auto v = fp64_normalize(fp64_cross(n1, n2));

  const accusphgeom::numeric::Vec3<EigenPack<N>> point_pos = v;
  const accusphgeom::numeric::Vec3<EigenPack<N>> point_neg = {
      -v[0], -v[1], -v[2]};
  return {point_pos, point_neg};
}

template <typename T>
inline auto fp64_try_gca_gca_intersection(
    const accusphgeom::numeric::Vec3<T>& a0,
    const accusphgeom::numeric::Vec3<T>& a1,
    const accusphgeom::numeric::Vec3<T>& b0,
    const accusphgeom::numeric::Vec3<T>& b1) {
  namespace numeric = accusphgeom::numeric;
  namespace predicates = accusphgeom::predicates;

  const auto candidates = fp64_gca_gca(a0, a1, b0, b1);
  const T tol = T(accusphgeom::constructions::internal::gca_gca_minor_arc_tol);
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

  accusphgeom::numeric::Vec3<T> out{};
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
  return accusphgeom::constructions::GcaGcaTryResult<T, StatusT>{out, status};
}

}  // namespace accusphgeom::performance_test
