#pragma once

#include <array>
#include <tuple>

#include "spip/core/types.hpp"
#include "spip/predicates/quadruple3d.hpp"
#include "spip/predicates/eft/basic.hpp"

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
inline std::tuple<T, T> compensated_quadruple_product(const V3_T<T>& a,
                                                      const V3_T<T>& b,
                                                      const V3_T<T>& c,
                                                      const V3_T<T>& d) {
  V3_T<T> ab_hi, ab_lo;
  V3_T<T> cd_hi, cd_lo;

  std::tie(ab_hi, ab_lo) = compensated_cross_product(a, b);
  std::tie(cd_hi, cd_lo) = compensated_cross_product(c, d);

  const std::array<T, 6> lhs = {
      ab_hi[0], ab_lo[0],
      ab_hi[1], ab_lo[1],
      ab_hi[2], ab_lo[2]};

  const std::array<T, 6> rhs = {
      cd_hi[0], cd_lo[0],
      cd_hi[1], cd_lo[1],
      cd_hi[2], cd_lo[2]};

  return compensated_dot_product<T, 6>(lhs, rhs);
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
  T hi, lo;
  std::tie(hi, lo) = compensated_quadruple_product(a, b, c, d);

  const T value = hi + lo;
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