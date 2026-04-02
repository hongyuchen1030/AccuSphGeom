#pragma once

#include <array>

#include "accusphgeom/numeric/eft.hpp"

namespace accusphgeom::constructions {

template <typename T>
using Vec3 = numeric::Vec3<T>;

template <typename T>
using Vec3Expansion2 = numeric::Vec3Expansion2<T>;

template <typename T>
inline Vec3Expansion2<T> accucross(const Vec3<T>& a, const Vec3<T>& b) {
  Vec3Expansion2<T> out{};
  const numeric::Expansion2<T> x =
      numeric::accurate_difference_of_products(a[1], b[2], a[2], b[1]);
  const numeric::Expansion2<T> y =
      numeric::accurate_difference_of_products(a[2], b[0], a[0], b[2]);
  const numeric::Expansion2<T> z =
      numeric::accurate_difference_of_products(a[0], b[1], a[1], b[0]);
  out.hi = {x.hi, y.hi, z.hi};
  out.lo = {x.lo, y.lo, z.lo};
  return out;
}

template <typename T>
inline Vec3Expansion2<T> accucross(const Vec3<T>& a_hi, const Vec3<T>& a_lo,
                                   const Vec3<T>& b_hi, const Vec3<T>& b_lo) {
  const bool a_exact =
      a_lo[0] == T(0) && a_lo[1] == T(0) && a_lo[2] == T(0);
  const bool b_exact =
      b_lo[0] == T(0) && b_lo[1] == T(0) && b_lo[2] == T(0);
  if (a_exact && b_exact) {
    return accucross(a_hi, b_hi);
  }

  const std::array<T, 8> x_lhs = {a_hi[1], a_hi[1], a_lo[1], a_lo[1],
                                  -a_hi[2], -a_hi[2], -a_lo[2], -a_lo[2]};
  const std::array<T, 8> x_rhs = {b_hi[2], b_lo[2], b_hi[2], b_lo[2],
                                  b_hi[1], b_lo[1], b_hi[1], b_lo[1]};
  const std::array<T, 8> y_lhs = {a_hi[2], a_hi[2], a_lo[2], a_lo[2],
                                  -a_hi[0], -a_hi[0], -a_lo[0], -a_lo[0]};
  const std::array<T, 8> y_rhs = {b_hi[0], b_lo[0], b_hi[0], b_lo[0],
                                  b_hi[2], b_lo[2], b_hi[2], b_lo[2]};
  const std::array<T, 8> z_lhs = {a_hi[0], a_hi[0], a_lo[0], a_lo[0],
                                  -a_hi[1], -a_hi[1], -a_lo[1], -a_lo[1]};
  const std::array<T, 8> z_rhs = {b_hi[1], b_lo[1], b_hi[1], b_lo[1],
                                  b_hi[0], b_lo[0], b_hi[0], b_lo[0]};

  const numeric::Expansion2<T> x = numeric::compensated_dot_product(x_lhs, x_rhs);
  const numeric::Expansion2<T> y = numeric::compensated_dot_product(y_lhs, y_rhs);
  const numeric::Expansion2<T> z = numeric::compensated_dot_product(z_lhs, z_rhs);

  return {{x.hi, y.hi, z.hi}, {x.lo, y.lo, z.lo}};
}

}  // namespace accusphgeom::constructions
