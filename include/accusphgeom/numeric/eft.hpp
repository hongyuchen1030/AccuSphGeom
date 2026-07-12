#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "accusphgeom/numeric/simd_fma.hpp"
#include "accusphgeom/numeric/sqrt.hpp"
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
  const T neg_x = -x;
  const T y = simd_fma(a, b, neg_x);
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
inline Expansion2<T> fast_two_sum(T a, T b) {
  T x = a + b;
  T y = (a - x) + b;
  return {x, y};
}

template <typename T>
inline Expansion2<T> accurate_difference_of_products(T a, T b, T c, T d) {
  const Expansion2<T> prod1 = two_prod_fma(a, b);
  const T neg_d = -d;
  const Expansion2<T> prod2 = two_prod_fma(c, neg_d);
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
inline T accurate_dot_product_fma(const std::array<T, N>& lhs,
                                             const std::array<T, N>& rhs) {
  Expansion2<T> acc = two_prod_fma(lhs[0], rhs[0]);
  for (std::size_t i = 1; i < N; ++i) {
    const Expansion2<T> prod = two_prod_fma(lhs[i], rhs[i]);
    const Expansion2<T> sum = two_sum(acc.hi, prod.hi);
    acc.hi = sum.hi;
    acc.lo += prod.lo + sum.lo;
  }
  return acc.hi + acc.lo;
}

template <typename T>
inline Expansion2<T> sum_non_neg(const Expansion2<T> & A,
                                    const Expansion2<T> & B) {
  const T Ah = A.hi;
  const T Al = A.lo;
  const T Bh = B.hi;
  const T Bl = B.lo;
  auto [H, h] = two_sum(Ah, Bh);
  T c = Al + Bl;
  T d = h + c;
  return fast_two_sum(H, d);
}

template <typename T, std::size_t N>
inline Expansion2<T> sum_of_squares_c(const std::array<T, N>& vals,
                                      const std::array<T, N>& errs) {
  T S = T(0), s = T(0);
  for (std::size_t i = 0; i < N; ++i) {
    auto [P, p] = two_prod_fma(vals[i], vals[i]);
    const Expansion2<T> sum =
        sum_non_neg(Expansion2<T>{S, s}, Expansion2<T>{P, p});
    S = sum.hi;
    s = sum.lo;
  }
  T R = accurate_dot_product_fma(vals, errs);
  T err = T(2) * R + s;
  return fast_two_sum(S, err);
}

template <typename T>
inline Expansion2<T> acc_sqrt_re(T value, T error = T(0)) {
  const T root = numeric_sqrt(value);
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
  return numeric_sqrt(dot3(v, v));
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
