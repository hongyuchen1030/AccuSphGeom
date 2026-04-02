#include "accusphgeom/predicates/quadruple3d.hpp"

#include "MultiPrecision_psm.h"

namespace accusphgeom::predicates::detail {

namespace {

Sign sign_from_geo(GEO::Sign sign) {
  switch (sign) {
    case GEO::POSITIVE:
      return Sign::Positive;
    case GEO::NEGATIVE:
      return Sign::Negative;
    default:
      return Sign::Zero;
  }
}

}  // namespace

Sign quadruple3d_exact_fallback(const double* a, const double* b,
                                const double* c, const double* d) {
  using GEO::expansion_nt;

  const auto dot3 = [](const double* u, const double* v) {
    return expansion_nt(u[0]) * expansion_nt(v[0]) +
           expansion_nt(u[1]) * expansion_nt(v[1]) +
           expansion_nt(u[2]) * expansion_nt(v[2]);
  };

  const expansion_nt ac = dot3(a, c);
  const expansion_nt bd = dot3(b, d);
  const expansion_nt ad = dot3(a, d);
  const expansion_nt bc = dot3(b, c);
  return sign_from_geo((ac * bd - ad * bc).sign());
}

}  // namespace accusphgeom::predicates::detail
