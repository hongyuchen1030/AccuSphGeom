#pragma once

#include <array>
#include <cstddef>
#include <tuple>

#include "spip/core/types.hpp"
#include "spip/predicates/eft/simd_fma.hh"

namespace spip {
namespace predicates {
namespace eft {

// -----------------------------------------------------------------------------
// Error-free product decomposition.
// Returns (hi, lo) such that
//
//   a b = hi + lo,
//
// where hi is the rounded product a * b in the working precision, and lo is the
// exact residual recovered by one fused multiply-add. Thus hi + lo represents
// the exact product in a two-term expansion.
// -----------------------------------------------------------------------------
template <typename T>
inline std::tuple<T, T> two_prod_fma(T a, T b) {
  const T hi = a * b;
  const T lo = ::simd_fma(a, b, -hi);
  return std::make_tuple(hi, lo);
}

// -----------------------------------------------------------------------------
// Error-free sum decomposition.
// Returns (hi, lo) such that
//
//   a + b = hi + lo,
//
// where hi is the rounded sum in the working precision and lo is the exact
// rounding residual. Thus hi + lo represents the exact sum in a two-term
// expansion.
// -----------------------------------------------------------------------------
template <typename T>
inline std::tuple<T, T> two_sum(T a, T b) {
  const T hi = a + b;
  const T z = hi - a;
  const T lo = (a - (hi - z)) + (b - z);
  return std::make_tuple(hi, lo);
}

// -----------------------------------------------------------------------------
// Accurate difference of two products.
// Returns (hi, lo) such that
//
//   a b - c d = hi + lo.
//
// This is the compensated evaluation of a 2x2 determinant. The routine first
// decomposes the two products exactly,
//
//   a b = p1 + e1,
//   c (-d) = p2_neg + e2,
//
// then combines their leading parts by an error-free summation, and finally
// accumulates all residual terms into lo.
// -----------------------------------------------------------------------------
template <typename T>
inline std::tuple<T, T> accu_dop(T a, T b, T c, T d) {
  T p1, e1;
  std::tie(p1, e1) = two_prod_fma(a, b);

  T p2_neg, e2;
  std::tie(p2_neg, e2) = two_prod_fma(c, -d);

  T hi, carry;
  std::tie(hi, carry) = two_sum(p1, p2_neg);

  const T lo = e1 + (carry + e2);
  return std::make_tuple(hi, lo);
}

// -----------------------------------------------------------------------------
// Compensated dot product.
// Returns (hi, lo) such that
//
//   sum_{i=0}^{N-1} a[i] b[i] = hi + lo.
//
// Each product a[i] b[i] is decomposed by two_prod_fma, and the running sum of
// leading terms is updated by two_sum. The residual product terms and summation
// residuals are accumulated into lo. The pair (hi, lo) is therefore a two-term
// compensated representation of the dot product.
// -----------------------------------------------------------------------------
template <typename T, std::size_t N>
inline std::tuple<T, T> compensated_dot_product(const std::array<T, N>& a,
                                                const std::array<T, N>& b) {
  T sum_hi, sum_lo;
  std::tie(sum_hi, sum_lo) = two_prod_fma(a[0], b[0]);

  for (std::size_t i = 1; i < N; ++i) {
    T prod_hi, prod_lo;
    std::tie(prod_hi, prod_lo) = two_prod_fma(a[i], b[i]);

    T new_hi, sigma;
    std::tie(new_hi, sigma) = two_sum(sum_hi, prod_hi);

    sum_hi = new_hi;
    sum_lo += prod_lo + sigma;
  }

  return std::make_tuple(sum_hi, sum_lo);
}

// -----------------------------------------------------------------------------
// Compensated cross product of two perturbed vectors.
//
// Inputs are interpreted as
//
//   u = v1 + ev1,
//   v = v2 + ev2,
//
// and the routine returns (hi, lo) such that
//
//   u x v = hi + lo
//
// componentwise, where hi and lo are 3-vectors.
//
// In the exact-input fast path, each component is evaluated as a compensated
// difference of two products:
//
//   (u x v)_x = u_y v_z - u_z v_y,
//   (u x v)_y = u_z v_x - u_x v_z,
//   (u x v)_z = u_x v_y - u_y v_x.
//
// In the general path, each component is expanded into the 8-term bilinear form
// induced by
//
//   (v1 + ev1) x (v2 + ev2),
//
// and the resulting component sum is evaluated by a compensated dot product.
// -----------------------------------------------------------------------------
template <typename T>
inline std::pair<V3_T<T>, V3_T<T>> compensated_cross_product(
    const V3_T<T>& v1, const V3_T<T>& ev1,
    const V3_T<T>& v2, const V3_T<T>& ev2) {
  V3_T<T> hi;
  V3_T<T> lo;

  const bool zero_ev1 = (ev1[0] == T(0) && ev1[1] == T(0) && ev1[2] == T(0));
  const bool zero_ev2 = (ev2[0] == T(0) && ev2[1] == T(0) && ev2[2] == T(0));

  if (zero_ev1 && zero_ev2) {
    T h, l;

    std::tie(h, l) = accu_dop(v1[1], v2[2], v1[2], v2[1]);
    hi[0] = h;
    lo[0] = l;

    std::tie(h, l) = accu_dop(v1[2], v2[0], v1[0], v2[2]);
    hi[1] = h;
    lo[1] = l;

    std::tie(h, l) = accu_dop(v1[0], v2[1], v1[1], v2[0]);
    hi[2] = h;
    lo[2] = l;

    return std::make_pair(hi, lo);
  }

  {
    const std::array<T, 8> a = {
        v1[1],  v1[1],  ev1[1],  ev1[1],
       -v1[2], -v1[2], -ev1[2], -ev1[2]};
    const std::array<T, 8> b = {
        v2[2],  ev2[2], v2[2],  ev2[2],
        v2[1],  ev2[1], v2[1],  ev2[1]};

    T h, l;
    std::tie(h, l) = compensated_dot_product<T, 8>(a, b);
    hi[0] = h;
    lo[0] = l;
  }

  {
    const std::array<T, 8> a = {
        v1[2],  v1[2],  ev1[2],  ev1[2],
       -v1[0], -v1[0], -ev1[0], -ev1[0]};
    const std::array<T, 8> b = {
        v2[0],  ev2[0], v2[0],  ev2[0],
        v2[2],  ev2[2], v2[2],  ev2[2]};

    T h, l;
    std::tie(h, l) = compensated_dot_product<T, 8>(a, b);
    hi[1] = h;
    lo[1] = l;
  }

  {
    const std::array<T, 8> a = {
        v1[0],  v1[0],  ev1[0],  ev1[0],
       -v1[1], -v1[1], -ev1[1], -ev1[1]};
    const std::array<T, 8> b = {
        v2[1],  ev2[1], v2[1],  ev2[1],
        v2[0],  ev2[0], v2[0],  ev2[0]};

    T h, l;
    std::tie(h, l) = compensated_dot_product<T, 8>(a, b);
    hi[2] = h;
    lo[2] = l;
  }

  return std::make_pair(hi, lo);
}

// -----------------------------------------------------------------------------
// Convenience overload for exact inputs.
// Returns (hi, lo) such that
//
//   v1 x v2 = hi + lo
//
// componentwise, with zero input perturbation vectors.
// -----------------------------------------------------------------------------
template <typename T>
inline std::pair<V3_T<T>, V3_T<T>> compensated_cross_product(
    const V3_T<T>& v1, const V3_T<T>& v2) {
  const V3_T<T> zero(T(0), T(0), T(0));
  return compensated_cross_product(v1, zero, v2, zero);
}

// -----------------------------------------------------------------------------
// Compensated scalar triple product.
// Returns an accurate approximation to
//
//   a · (b x c).
//
// The cross product b x c is first represented componentwise in two-term form,
//
//   (b x c)_k = h_k + l_k,   k = 0,1,2,
//
// using three compensated differences of products. The final triple product is
// then evaluated as the compensated dot product
//
//   a_0 (h_0 + l_0) + a_1 (h_1 + l_1) + a_2 (h_2 + l_2).
//
// This quantity is the determinant det[a, b, c], i.e., the orient3d sign on the
// sphere up to the usual geometric interpretation.
// -----------------------------------------------------------------------------
template <typename T>
inline T compensated_triple_product(const V3_T<T>& a,
                                    const V3_T<T>& b,
                                    const V3_T<T>& c) {
  T h0, l0;
  T h1, l1;
  T h2, l2;

  std::tie(h0, l0) = accu_dop(b[1], c[2], b[2], c[1]);
  std::tie(h1, l1) = accu_dop(b[2], c[0], b[0], c[2]);
  std::tie(h2, l2) = accu_dop(b[0], c[1], b[1], c[0]);

  const std::array<T, 6> lhs = {a[0], a[0], a[1], a[1], a[2], a[2]};
  const std::array<T, 6> rhs = {h0, l0, h1, l1, h2, l2};

  T hi, lo;
  std::tie(hi, lo) = compensated_dot_product<T, 6>(lhs, rhs);
  return hi + lo;
}

}  // namespace eft
}  // namespace predicates
}  // namespace spip