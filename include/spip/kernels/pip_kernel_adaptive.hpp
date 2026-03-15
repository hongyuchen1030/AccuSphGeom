#pragma once

#include "spip/predicates/orient3d.hpp"
#include "spip/predicates/quadruple3d.hpp"

namespace spip {
namespace kernels {

struct PIPKernelAdaptive {
  using Sign = spip::predicates::Sign;

  static Sign orient3d_on_sphere(const double* a,
                                 const double* b,
                                 const double* c) {
    return spip::predicates::orient3d_on_sphere(a, b, c);
  }

  static Sign quadruple3d(const double* a,
                          const double* b,
                          const double* c,
                          const double* d) {
    return spip::predicates::quadruple3d(a, b, c, d);
  }
};

}  // namespace kernels
}  // namespace spip