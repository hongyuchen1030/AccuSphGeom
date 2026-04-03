#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "accusphgeom/numeric/simd_fma.hpp"

namespace accusphgeom::numeric {

template <typename T>
using Vec3 = std::array<T, 3>;

template <typename T>
struct Expansion2 {
  T hi{};
  T lo{};
};

template <typename T>
struct Vec3Expansion2 {
  Vec3<T> hi{};
  Vec3<T> lo{};
};

template <typename T>
inline Expansion2<T> two_prod_fma(T a, T b) {
  const T x = a * b;
  const T y = simd_fma(a, b, -x);
  return {x, y};
}

template <typename T>
inline Expansion2<T> two_square_fma(T a) {
  return two_prod_fma(a, a);
}

template <typename T>
inline Expansion2<T> two_sum(T a, T b) {
  const T x = a + b;
  const T z = x - a;
  const T y = (a - (x - z)) + (b - z);
  return {x, y};
}

template <typename T>
inline Expansion2<T> accurate_difference_of_products(T a, T b, T c, T d) {
  const Expansion2<T> prod1 = two_prod_fma(a, b);
  const Expansion2<T> prod2 = two_prod_fma(c, -d);
  const Expansion2<T> sum = two_sum(prod1.hi, prod2.hi);
  return {sum.hi, prod1.lo + (sum.lo + prod2.lo)};
}

template <typename T, std::size_t N>
inline Expansion2<T> compensated_dot_product(const std::array<T, N>& lhs,
                                             const std::array<T, N>& rhs) {
  Expansion2<T> acc = two_prod_fma(lhs[0], rhs[0]);
  for (std::size_t i = 1; i < N; ++i) {
    const Expansion2<T> prod = two_prod_fma(lhs[i], rhs[i]);
    const Expansion2<T> sum = two_sum(acc.hi, prod.hi);
    acc.hi = sum.hi;
    acc.lo += prod.lo + sum.lo;
  }
  return acc;
}

template <typename T, std::size_t N>
inline Expansion2<T> sum_of_squares_c(const std::array<T, N>& high,
                                      const std::array<T, N>& low) {
  std::array<T, 2 * N> lhs{};
  std::array<T, 2 * N> rhs{};
  for (std::size_t i = 0; i < N; ++i) {
    lhs[2 * i] = high[i];
    rhs[2 * i] = high[i];
    lhs[2 * i + 1] = low[i];
    rhs[2 * i + 1] = low[i];
  }
  return compensated_dot_product(lhs, rhs);
}

template <typename T>
inline Expansion2<T> acc_sqrt_re(T value, T error = T(0)) {
  if (value < T(0)) {
    throw std::domain_error("acc_sqrt_re: negative radicand");
  }
  if (value == T(0)) {
    return {T(0), T(0)};
  }
  const T root = std::sqrt(value);
  const Expansion2<T> square = two_square_fma(root);
  const T residual = (value - square.hi) - square.lo + error;
  const T correction = residual / (T(2) * root);
  return {root, correction};
}

template <typename T>
inline T expansion_value(const Expansion2<T>& x) {
  return x.hi + x.lo;
}

template <typename T>
inline T dot3(const Vec3<T>& a, const Vec3<T>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

template <typename T>
inline T norm3(const Vec3<T>& v) {
  return std::sqrt(dot3(v, v));
}

template <typename T>
inline Vec3<T> normalize(const Vec3<T>& v) {
  const T n = norm3(v);
  if (n == T(0)) {
    throw std::invalid_argument("normalize: zero vector");
  }
  return {v[0] / n, v[1] / n, v[2] / n};
}

}  // namespace accusphgeom::numeric
