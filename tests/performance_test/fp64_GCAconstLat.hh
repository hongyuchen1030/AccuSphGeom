#pragma once

#include <cmath>

#include "accusphgeom/adapters/eigen/numeric.hpp"
#include "accusphgeom/constructions/gca_constlat_intersection.hpp"

namespace accusphgeom::performance_test {

inline accusphgeom::constructions::GcaConstLatIntersections<double>
fp64_gca_constlat_pure_fp(
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
fp64_gca_constlat_pure_fp(
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

}  // namespace accusphgeom::performance_test
