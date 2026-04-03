#pragma once

#include <cstddef>

#include "spip/core/types.hpp"
#include "spip/predicates/eft/simd_fma.hh"

namespace spip {
namespace predicates {
namespace eft {

template <typename T>
struct TwoTerm {
  T hi;
  T lo;
};

template <typename T>
struct TwoTermVec3 {
  V3_T<T> hi;
  V3_T<T> lo;
};

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
inline TwoTerm<T> two_prod_fma(T a, T b) {
  const T hi = a * b;
  const T lo = ::simd_fma(a, b, -hi);
  return {hi, lo};
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
inline TwoTerm<T> two_sum(T a, T b) {
  const T hi = a + b;
  const T z = hi - a;
  const T lo = (a - (hi - z)) + (b - z);
  return {hi, lo};
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
inline TwoTerm<T> accu_dop(T a, T b, T c, T d) {
  const TwoTerm<T> prod1 = two_prod_fma(a, b);
  const TwoTerm<T> prod2_neg = two_prod_fma(c, -d);
  const TwoTerm<T> sum = two_sum(prod1.hi, prod2_neg.hi);
  return {sum.hi, prod1.lo + (sum.lo + prod2_neg.lo)};
}

template <typename T>
inline TwoTerm<T> compensated_dot_product_6(
    T a0, T b0, T a1, T b1, T a2, T b2,
    T a3, T b3, T a4, T b4, T a5, T b5) {
  TwoTerm<T> sum = two_prod_fma(a0, b0);

  const TwoTerm<T> prod1 = two_prod_fma(a1, b1);
  const TwoTerm<T> add1 = two_sum(sum.hi, prod1.hi);
  sum.hi = add1.hi;
  sum.lo += prod1.lo + add1.lo;

  const TwoTerm<T> prod2 = two_prod_fma(a2, b2);
  const TwoTerm<T> add2 = two_sum(sum.hi, prod2.hi);
  sum.hi = add2.hi;
  sum.lo += prod2.lo + add2.lo;

  const TwoTerm<T> prod3 = two_prod_fma(a3, b3);
  const TwoTerm<T> add3 = two_sum(sum.hi, prod3.hi);
  sum.hi = add3.hi;
  sum.lo += prod3.lo + add3.lo;

  const TwoTerm<T> prod4 = two_prod_fma(a4, b4);
  const TwoTerm<T> add4 = two_sum(sum.hi, prod4.hi);
  sum.hi = add4.hi;
  sum.lo += prod4.lo + add4.lo;

  const TwoTerm<T> prod5 = two_prod_fma(a5, b5);
  const TwoTerm<T> add5 = two_sum(sum.hi, prod5.hi);
  sum.hi = add5.hi;
  sum.lo += prod5.lo + add5.lo;

  return sum;
}

template <typename T>
inline TwoTerm<T> compensated_dot_product_8(
    T a0, T b0, T a1, T b1, T a2, T b2, T a3, T b3,
    T a4, T b4, T a5, T b5, T a6, T b6, T a7, T b7) {
  TwoTerm<T> sum = two_prod_fma(a0, b0);

  const TwoTerm<T> prod1 = two_prod_fma(a1, b1);
  const TwoTerm<T> add1 = two_sum(sum.hi, prod1.hi);
  sum.hi = add1.hi;
  sum.lo += prod1.lo + add1.lo;

  const TwoTerm<T> prod2 = two_prod_fma(a2, b2);
  const TwoTerm<T> add2 = two_sum(sum.hi, prod2.hi);
  sum.hi = add2.hi;
  sum.lo += prod2.lo + add2.lo;

  const TwoTerm<T> prod3 = two_prod_fma(a3, b3);
  const TwoTerm<T> add3 = two_sum(sum.hi, prod3.hi);
  sum.hi = add3.hi;
  sum.lo += prod3.lo + add3.lo;

  const TwoTerm<T> prod4 = two_prod_fma(a4, b4);
  const TwoTerm<T> add4 = two_sum(sum.hi, prod4.hi);
  sum.hi = add4.hi;
  sum.lo += prod4.lo + add4.lo;

  const TwoTerm<T> prod5 = two_prod_fma(a5, b5);
  const TwoTerm<T> add5 = two_sum(sum.hi, prod5.hi);
  sum.hi = add5.hi;
  sum.lo += prod5.lo + add5.lo;

  const TwoTerm<T> prod6 = two_prod_fma(a6, b6);
  const TwoTerm<T> add6 = two_sum(sum.hi, prod6.hi);
  sum.hi = add6.hi;
  sum.lo += prod6.lo + add6.lo;

  const TwoTerm<T> prod7 = two_prod_fma(a7, b7);
  const TwoTerm<T> add7 = two_sum(sum.hi, prod7.hi);
  sum.hi = add7.hi;
  sum.lo += prod7.lo + add7.lo;

  return sum;
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
inline TwoTerm<T> compensated_dot_product(const T (&a)[N], const T (&b)[N]) {
  if constexpr (N == 6) {
    return compensated_dot_product_6(
        a[0], b[0], a[1], b[1], a[2], b[2],
        a[3], b[3], a[4], b[4], a[5], b[5]);
  } else if constexpr (N == 8) {
    return compensated_dot_product_8(
        a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3],
        a[4], b[4], a[5], b[5], a[6], b[6], a[7], b[7]);
  }

  TwoTerm<T> sum = two_prod_fma(a[0], b[0]);
  for (std::size_t i = 1; i < N; ++i) {
    const TwoTerm<T> prod = two_prod_fma(a[i], b[i]);
    const TwoTerm<T> add = two_sum(sum.hi, prod.hi);
    sum.hi = add.hi;
    sum.lo += prod.lo + add.lo;
  }
  return sum;
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
inline TwoTermVec3<T> compensated_cross_product(
    const V3_T<T>& v1, const V3_T<T>& ev1,
    const V3_T<T>& v2, const V3_T<T>& ev2) {
  TwoTermVec3<T> result;

  const bool zero_ev1 = (ev1[0] == T(0) && ev1[1] == T(0) && ev1[2] == T(0));
  const bool zero_ev2 = (ev2[0] == T(0) && ev2[1] == T(0) && ev2[2] == T(0));

  if (zero_ev1 && zero_ev2) {
    const TwoTerm<T> x = accu_dop(v1[1], v2[2], v1[2], v2[1]);
    const TwoTerm<T> y = accu_dop(v1[2], v2[0], v1[0], v2[2]);
    const TwoTerm<T> z = accu_dop(v1[0], v2[1], v1[1], v2[0]);
    result.hi[0] = x.hi;
    result.lo[0] = x.lo;
    result.hi[1] = y.hi;
    result.lo[1] = y.lo;
    result.hi[2] = z.hi;
    result.lo[2] = z.lo;
    return result;
  }

  result.hi.setZero();
  result.lo.setZero();

  const TwoTerm<T> x = compensated_dot_product_8(
      v1[1], v2[2], v1[1], ev2[2], ev1[1], v2[2], ev1[1], ev2[2],
      -v1[2], v2[1], -v1[2], ev2[1], -ev1[2], v2[1], -ev1[2], ev2[1]);
  result.hi[0] = x.hi;
  result.lo[0] = x.lo;

  const TwoTerm<T> y = compensated_dot_product_8(
      v1[2], v2[0], v1[2], ev2[0], ev1[2], v2[0], ev1[2], ev2[0],
      -v1[0], v2[2], -v1[0], ev2[2], -ev1[0], v2[2], -ev1[0], ev2[2]);
  result.hi[1] = y.hi;
  result.lo[1] = y.lo;

  const TwoTerm<T> z = compensated_dot_product_8(
      v1[0], v2[1], v1[0], ev2[1], ev1[0], v2[1], ev1[0], ev2[1],
      -v1[1], v2[0], -v1[1], ev2[0], -ev1[1], v2[0], -ev1[1], ev2[0]);
  result.hi[2] = z.hi;
  result.lo[2] = z.lo;

  return result;
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
inline TwoTermVec3<T> compensated_cross_product(
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
  const TwoTerm<T> x = accu_dop(b[1], c[2], b[2], c[1]);
  const TwoTerm<T> y = accu_dop(b[2], c[0], b[0], c[2]);
  const TwoTerm<T> z = accu_dop(b[0], c[1], b[1], c[0]);
  const TwoTerm<T> dot = compensated_dot_product_6(
      a[0], x.hi, a[0], x.lo,
      a[1], y.hi, a[1], y.lo,
      a[2], z.hi, a[2], z.lo);
  return dot.hi + dot.lo;
}

}  // namespace eft
}  // namespace predicates
}  // namespace spip
