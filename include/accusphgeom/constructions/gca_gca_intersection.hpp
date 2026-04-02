#pragma once

#include <stdexcept>

#include "accusphgeom/constructions/accucross.hpp"

namespace accusphgeom::constructions {

template <typename T>
inline numeric::Vec3<T> gca_gca_intersection(const numeric::Vec3<T>& a0,
                                             const numeric::Vec3<T>& a1,
                                             const numeric::Vec3<T>& b0,
                                             const numeric::Vec3<T>& b1) {
  const auto n1 = accucross(a0, a1);
  const auto n2 = accucross(b0, b1);
  const auto v = accucross(n1.hi, n1.lo, n2.hi, n2.lo);
  const auto sum = numeric::sum_of_squares_c<T, 3>(v.hi, v.lo);
  const auto norm = numeric::acc_sqrt_re(sum.hi, sum.lo);
  const T n = norm.hi + norm.lo;
  if (n == T(0)) {
    throw std::domain_error(
        "gca_gca_intersection: great circles are degenerate or coincident");
  }
  return {(v.hi[0] + v.lo[0]) / n, (v.hi[1] + v.lo[1]) / n,
          (v.hi[2] + v.lo[2]) / n};
}

}  // namespace accusphgeom::constructions
