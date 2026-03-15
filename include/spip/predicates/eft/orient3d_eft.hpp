#pragma once

#include "spip/core/types.hpp"
#include "spip/predicates/orient3d.hpp"
#include "spip/predicates/eft/basic.hpp"

namespace spip {
namespace predicates {
namespace eft {

// -----------------------------------------------------------------------------
// EFT-based sign of the scalar triple product
//
//   det[a, b, c] = a · (b x c).
//
// For spherical predicates, this is the sign of orient3d_on_sphere(a, b, c).
// The triple product is evaluated by compensated arithmetic in basic.hpp.
// -----------------------------------------------------------------------------
template <typename T>
inline Sign orient3d_on_sphere(const V3_T<T>& a,
                               const V3_T<T>& b,
                               const V3_T<T>& c) {
  const T det = compensated_triple_product(a, b, c);

  if (det > T(0)) {
    return Sign::Positive;
  }
  if (det < T(0)) {
    return Sign::Negative;
  }
  return Sign::Zero;
}

}  // namespace eft
}  // namespace predicates
}  // namespace spip