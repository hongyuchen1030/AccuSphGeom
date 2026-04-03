#pragma once

#include "spip/core/types.hpp"
#include "spip/predicates/orient3d.hpp"
#include "spip/predicates/eft/orient3d_eft.hpp"
#include "spip/predicates/eft/quadruple3d_eft.hpp"

namespace spip {
namespace kernels {

struct PIPKernelEFT {
  using Sign = spip::predicates::Sign;

  template <typename T>
  static Sign orient3d_on_sphere(const V3_T<T>& a,
                                 const V3_T<T>& b,
                                 const V3_T<T>& c) {
    return spip::predicates::eft::orient3d_on_sphere(a, b, c);
  }

  template <typename T>
  static Sign quadruple3d(const V3_T<T>& a,
                          const V3_T<T>& b,
                          const V3_T<T>& c,
                          const V3_T<T>& d) {
    return spip::predicates::eft::quadruple3d(a, b, c, d);
  }
};

}  // namespace kernels
}  // namespace spip