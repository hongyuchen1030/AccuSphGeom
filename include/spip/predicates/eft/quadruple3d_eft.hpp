#pragma once

#include <cmath>
#include <limits>

#include "spip/core/types.hpp"
#include "spip/predicates/quadruple3d.hpp"
#include "spip/predicates/eft/basic.hpp"
#include "spip/predicates/eft/orient3d_eft.hpp"

namespace spip {
namespace predicates {
namespace eft {

// -----------------------------------------------------------------------------
// EFT-based quadruple product
//
//   (a x b) · (c x d).
//
// Let
//
//   a x b = u_h + u_l,
//   c x d = v_h + v_l,
//
// where each component of u_h, u_l, v_h, and v_l is produced by
// compensated_cross_product(). The final dot product is then evaluated as
//
//   (u_h + u_l) · (v_h + v_l)
//
// by one compensated dot product over the 6-term vectors
//
//   [u_hx, u_lx, u_hy, u_ly, u_hz, u_lz],
//   [v_hx, v_lx, v_hy, v_ly, v_hz, v_lz].
//
// The return value is the compensated pair (hi, lo) such that
//
//   (a x b) · (c x d) = hi + lo.
// -----------------------------------------------------------------------------
template <typename T>
inline TwoTerm<T> compensated_quadruple_product(const V3_T<T>& a,
                                                const V3_T<T>& b,
                                                const V3_T<T>& c,
                                                const V3_T<T>& d) {
  const TwoTermVec3<T> ab = compensated_cross_product(a, b);
  const TwoTermVec3<T> cd = compensated_cross_product(c, d);
  return compensated_dot_product_6(
      ab.hi[0], cd.hi[0], ab.lo[0], cd.lo[0],
      ab.hi[1], cd.hi[1], ab.lo[1], cd.lo[1],
      ab.hi[2], cd.hi[2], ab.lo[2], cd.lo[2]);
}

// -----------------------------------------------------------------------------
// EFT-based sign of the quadruple product
//
//   (a x b) · (c x d).
//
// The quantity is evaluated by compensated_quadruple_product() and classified
// by the sign of hi + lo.
// -----------------------------------------------------------------------------
template <typename T>
inline Sign quadruple3d(const V3_T<T>& a,
                        const V3_T<T>& b,
                        const V3_T<T>& c,
                        const V3_T<T>& d) {
  const TwoTerm<T> value_terms = compensated_quadruple_product(a, b, c, d);
  const T value = value_terms.hi + value_terms.lo;
  if (std::abs(value) < zero_tolerance<T>) {
    return Sign::Zero;
  }
  if (value > T(0)) {
    return Sign::Positive;
  }
  if (value < T(0)) {
    return Sign::Negative;
  }
  return Sign::Zero;
}

}  // namespace eft
}  // namespace predicates
}  // namespace spip
