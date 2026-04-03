#pragma once

#include <array>
#include <stdexcept>

#include "accusphgeom/constructions/accucross.hpp"

namespace accusphgeom::constructions {

template <typename T>
struct GcaConstLatIntersection {
  numeric::Vec3<T> point{};
  numeric::Vec3<T> normal_high{};
  numeric::Vec3<T> normal_low{};
};

template <typename T>
inline GcaConstLatIntersection<T> gca_constlat_intersection(
    const numeric::Vec3<T>& a, const numeric::Vec3<T>& b, T z0,
    int branch_sign = 1) {
  if (branch_sign != 1 && branch_sign != -1) {
    throw std::invalid_argument(
        "gca_constlat_intersection: branch_sign must be +1 or -1");
  }

  const auto normal = accucross(a, b);
  const auto s2 = numeric::sum_of_squares_c<T, 2>(
      {normal.hi[0], normal.hi[1]}, {normal.lo[0], normal.lo[1]});
  const auto s3 = numeric::sum_of_squares_c<T, 3>(normal.hi, normal.lo);
  const auto zsq = numeric::two_prod_fma(z0, z0);
  const auto d = numeric::compensated_dot_product(
      std::array<T, 4>{s3.hi, s3.hi, s3.lo, s3.lo},
      std::array<T, 4>{zsq.hi, zsq.lo, zsq.hi, zsq.lo});
  const auto e = numeric::two_sum(s2.hi, -d.hi);
  const T planar_sq = e.hi + (e.lo + s2.lo - d.lo);
  const auto s = numeric::acc_sqrt_re(planar_sq);

  const T nx = normal.hi[0] + normal.lo[0];
  const T ny = normal.hi[1] + normal.lo[1];
  const T nz = normal.hi[2] + normal.lo[2];
  const T planar = s.hi + s.lo;
  const T denom = s2.hi + s2.lo;
  if (denom == T(0)) {
    throw std::domain_error(
        "gca_constlat_intersection: degenerate great-circle normal");
  }

  const auto x_num = numeric::compensated_dot_product(
      std::array<T, 2>{nx * nz, -static_cast<T>(branch_sign) * ny},
      std::array<T, 2>{z0, planar});
  const auto y_num = numeric::compensated_dot_product(
      std::array<T, 2>{ny * nz, static_cast<T>(branch_sign) * nx},
      std::array<T, 2>{z0, planar});

  GcaConstLatIntersection<T> out{};
  out.normal_high = normal.hi;
  out.normal_low = normal.lo;
  out.point = {-(x_num.hi + x_num.lo) / denom, -(y_num.hi + y_num.lo) / denom,
               z0};
  return out;
}

}  // namespace accusphgeom::constructions
