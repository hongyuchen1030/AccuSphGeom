#include "accusphgeom/predicates/quadruple3d.hpp"

// Geogram PSM exact arithmetic
#include "MultiPrecision_psm.h"

namespace accusphgeom::predicates::internal {

static inline Sign sign_from_geo(GEO::Sign s) {
  switch (s) {
    case GEO::POSITIVE: return Sign::Positive;
    case GEO::NEGATIVE: return Sign::Negative;
    default:            return Sign::Zero;
  }
}

Sign quadruple3d_exact_fallback(const double* a, const double* b,
                                const double* c, const double* d) {
  using GEO::expansion_nt;

  auto dot3 = [](const double* u, const double* v) -> expansion_nt {
    return expansion_nt(u[0]) * expansion_nt(v[0]) +
           expansion_nt(u[1]) * expansion_nt(v[1]) +
           expansion_nt(u[2]) * expansion_nt(v[2]);
  };

  const expansion_nt ac = dot3(a, c);
  const expansion_nt bd = dot3(b, d);
  const expansion_nt ad = dot3(a, d);
  const expansion_nt bc = dot3(b, c);

  const expansion_nt Delta = ac * bd - ad * bc;
  return sign_from_geo(Delta.sign());
}

}  // namespace accusphgeom::predicates::internal
