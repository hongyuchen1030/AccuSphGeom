#pragma once

#include <Eigen/Dense>

#include <tuple>

#include "accusphgeom/numeric/eft.hpp"

namespace accusphgeom::performance_test {

template <typename T>
std::tuple<T, T> fp64_gca_constlat_pure_fp_xy(
    const accusphgeom::numeric::Vec3<T>& pointA,
    const accusphgeom::numeric::Vec3<T>& pointB,
    const T& constZ);

inline std::tuple<double, double>
fp64_gca_constlat_pure_fp_xy(const Eigen::Vector3d& pointA,
                             const Eigen::Vector3d& pointB,
                             const double constZ) {
  return fp64_gca_constlat_pure_fp_xy(
      accusphgeom::numeric::Vec3<double>{pointA[0], pointA[1], pointA[2]},
      accusphgeom::numeric::Vec3<double>{pointB[0], pointB[1], pointB[2]},
      constZ);
}

}  // namespace accusphgeom::performance_test
