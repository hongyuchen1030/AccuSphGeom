#pragma once

#include <Eigen/Dense>

#include "accusphgeom/numeric/mask.hpp"

template <int N>
using EigenPack = Eigen::Array<double, N, 1>;

namespace accusphgeom::numeric {

template <int N>
inline Eigen::Array<double, N, 1>
isfinite_mask(const Eigen::Array<double, N, 1>& x) {
  return x.isFinite().template cast<double>();
}

template <int N>
inline Eigen::Array<double, N, 1>
numeric_abs(const Eigen::Array<double, N, 1>& x) {
  return x.abs();
}

template <int N>
inline Eigen::Array<double, N, 1>
mask_equal(const Eigen::Array<double, N, 1>& a,
           const Eigen::Array<double, N, 1>& b) {
  return (a == b).template cast<double>();
}

template <int N>
inline Eigen::Array<double, N, 1>
mask_le(const Eigen::Array<double, N, 1>& a,
        const Eigen::Array<double, N, 1>& b) {
  return (a <= b).template cast<double>();
}

template <int N>
inline Eigen::Array<double, N, 1>
mask_ge(const Eigen::Array<double, N, 1>& a,
        const Eigen::Array<double, N, 1>& b) {
  return (a >= b).template cast<double>();
}

}  // namespace accusphgeom::numeric
