#pragma once

#include <cmath>

#include "accusphgeom/adapters/eigen/numeric.hpp"
#include "accusphgeom/constructions/gca_constlat_intersection.hpp"

namespace accusphgeom::performance_test {

inline accusphgeom::constructions::GcaConstLatIntersections<double>
fp64_gca_constlat(
    const accusphgeom::numeric::Vec3<double>& pointA,
    const accusphgeom::numeric::Vec3<double>& pointB,
    double constZ) {
  using accusphgeom::constructions::GcaConstLatIntersections;

  const double nx = pointA[1] * pointB[2] - pointA[2] * pointB[1];
  const double ny = pointA[2] * pointB[0] - pointA[0] * pointB[2];
  const double nz = pointA[0] * pointB[1] - pointA[1] * pointB[0];

  const double nx2 = nx * nx;
  const double ny2 = ny * ny;
  const double nz2 = nz * nz;

  const double denom = nx2 + ny2;
  const double norm_n2 = denom + nz2;
  const double s = std::sqrt(denom - norm_n2 * constZ * constZ);

  GcaConstLatIntersections<double> out{};

  out.point_pos = {
      -(constZ * nx * nz - s * ny) / denom,
      -(constZ * ny * nz + s * nx) / denom,
      constZ};

  out.point_neg = {
      -(constZ * nx * nz + s * ny) / denom,
      -(constZ * ny * nz - s * nx) / denom,
      constZ};

  return out;
}

template <int N>
inline accusphgeom::constructions::GcaConstLatIntersections<EigenPack<N>>
fp64_gca_constlat(
    const accusphgeom::numeric::Vec3<EigenPack<N>>& pointA,
    const accusphgeom::numeric::Vec3<EigenPack<N>>& pointB,
    const EigenPack<N>& constZ) {
  using accusphgeom::constructions::GcaConstLatIntersections;
  using T = EigenPack<N>;
  using std::sqrt;

  const T nx = pointA[1] * pointB[2] - pointA[2] * pointB[1];
  const T ny = pointA[2] * pointB[0] - pointA[0] * pointB[2];
  const T nz = pointA[0] * pointB[1] - pointA[1] * pointB[0];

  const T nx2 = nx * nx;
  const T ny2 = ny * ny;
  const T nz2 = nz * nz;

  const T denom = nx2 + ny2;
  const T norm_n2 = denom + nz2;
  const T s = sqrt(denom - norm_n2 * constZ * constZ);

  GcaConstLatIntersections<T> out{};

  out.point_pos = {
      -(constZ * nx * nz - s * ny) / denom,
      -(constZ * ny * nz + s * nx) / denom,
      constZ};

  out.point_neg = {
      -(constZ * nx * nz + s * ny) / denom,
      -(constZ * ny * nz - s * nx) / denom,
      constZ};

  return out;
}


template <typename T>
inline accusphgeom::constructions::GcaConstLatTryResult<T>
fp64_try_gca_constlat_intersection(
    const accusphgeom::numeric::Vec3<T>& a,
    const accusphgeom::numeric::Vec3<T>& b,
    T z0) {
  namespace numeric = accusphgeom::numeric;
  namespace predicates = accusphgeom::predicates;

  const auto c = fp64_gca_constlat(a, b, z0);

  const T tol =
      T(accusphgeom::constructions::internal::gca_constlat_minor_arc_tol);
  const T one = T(1);

  const auto pos_finite =
      numeric::isfinite_mask(c.point_pos[0]) *
      numeric::isfinite_mask(c.point_pos[1]);

  const auto neg_finite =
      numeric::isfinite_mask(c.point_neg[0]) *
      numeric::isfinite_mask(c.point_neg[1]);

  const auto pos_on_arc =
      pos_finite *
      predicates::internal::on_minor_arc_tol_ptr(c.point_pos.data(), a.data(),
                                                 b.data(), tol);

  const auto neg_on_arc =
      neg_finite *
      predicates::internal::on_minor_arc_tol_ptr(c.point_neg.data(), a.data(),
                                                 b.data(), tol);

  const auto pos_valid = pos_finite * pos_on_arc;
  const auto neg_valid = neg_finite * neg_on_arc;

  const auto pos_mask = pos_valid * (one - neg_valid);
  const auto neg_mask = neg_valid * (one - pos_valid);

  accusphgeom::numeric::Vec3<T> out{};
  out[0] = pos_mask * c.point_pos[0] + neg_mask * c.point_neg[0];
  out[1] = pos_mask * c.point_pos[1] + neg_mask * c.point_neg[1];
  out[2] = pos_mask * c.point_pos[2] + neg_mask * c.point_neg[2];

  const auto both = pos_valid * neg_valid;
  const auto none = (one - pos_valid) * (one - neg_valid);
  const auto status = both + none * T(2);

  return accusphgeom::constructions::GcaConstLatTryResult<T>{out, status};
}


}  // namespace accusphgeom::performance_test
