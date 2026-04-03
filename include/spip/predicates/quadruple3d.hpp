#pragma once

#include <array>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <vector>

// Reuse the shared Sign enum and common internal helpers.
#include "spip/predicates/orient3d.hpp"

// Shewchuk wrapper header used by orient3d.hpp for SPIP_PREDICATES_REAL +
// exactinit().
#include "predicates.h"

// Config: defines FPG_UNCERTAIN_VALUE and pulls in fabs().
#include "spip/predicates/pck_filter_config.hpp"

// Vendored PCK filter. Your build already exposes third_party/ as an include
// dir.
#include "quadruple3d.h"

namespace spip::predicates {

namespace internal {

using Real = SPIP_PREDICATES_REAL;

// Exact fallback implemented in src/predicates/quadruple3d_exact.cpp.
Sign quadruple3d_exact_fallback(const double* a, const double* b,
                                const double* c, const double* d);

inline Sign sign_from_int(int v) {
  if (v > 0) return Sign::Positive;
  if (v < 0) return Sign::Negative;
  return Sign::Zero;
}

inline Sign quadruple3d_raw_ptr(const Real* a,
                                const Real* b,
                                const Real* c,
                                const Real* d) {
  require_nonnull(a, "a");
  require_nonnull(b, "b");
  require_nonnull(c, "c");
  require_nonnull(d, "d");

  // Keep the same initialization pattern as orient3d.hpp.
  ensure_exactinit();

#ifdef SPIP_PREDICATES_USE_FLOAT
  static_assert(!std::is_same_v<Real, float>,
                "spip::predicates::quadruple3d: float mode not supported yet "
                "(PCK filter is double).");
#endif

  const int s = quadruple_3d_filter(
      reinterpret_cast<const double*>(a),
      reinterpret_cast<const double*>(b),
      reinterpret_cast<const double*>(c),
      reinterpret_cast<const double*>(d));

  if (s != FPG_UNCERTAIN_VALUE) {
    return sign_from_int(s);
  }

  return quadruple3d_exact_fallback(
      reinterpret_cast<const double*>(a),
      reinterpret_cast<const double*>(b),
      reinterpret_cast<const double*>(c),
      reinterpret_cast<const double*>(d));
}

}  // namespace internal

// Pointer API.
template <class T>
inline Sign quadruple3d(const T* a, const T* b, const T* c, const T* d) {
  internal::require_supported_and_matches_real<T>();
  return internal::quadruple3d_raw_ptr(a, b, c, d);
}

// std::array API.
template <class T>
inline Sign quadruple3d(const std::array<T, 3>& a,
                        const std::array<T, 3>& b,
                        const std::array<T, 3>& c,
                        const std::array<T, 3>& d) {
  internal::require_supported_and_matches_real<T>();
  return internal::quadruple3d_raw_ptr(a.data(), b.data(), c.data(), d.data());
}

// std::vector API.
template <class T>
inline Sign quadruple3d(const std::vector<T>& a,
                        const std::vector<T>& b,
                        const std::vector<T>& c,
                        const std::vector<T>& d) {
  internal::require_supported_and_matches_real<T>();
  internal::require_size3(a.size(), "a");
  internal::require_size3(b.size(), "b");
  internal::require_size3(c.size(), "c");
  internal::require_size3(d.size(), "d");
  return internal::quadruple3d_raw_ptr(a.data(), b.data(), c.data(), d.data());
}

}  // namespace spip::predicates