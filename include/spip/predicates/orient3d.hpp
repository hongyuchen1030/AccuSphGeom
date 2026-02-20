// include/spip/predicates/orient3d.hpp
#pragma once

#include <array>
#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <type_traits>
#include <vector>

// third_party/predicates.h (Shewchuk wrapper header)
// Your build must add `third_party/` to include dirs.
#include "predicates.h"

namespace spip::predicates {

enum class Sign : int8_t { Negative = -1, Zero = 0, Positive = 1 };

namespace internal {

using Real = SPIP_PREDICATES_REAL;

inline void ensure_exactinit() {
  static std::once_flag flag;
  std::call_once(flag, []() { exactinit(); });
}

inline Sign sign_from_value(Real v) {
  if (v > Real(0)) return Sign::Positive;
  if (v < Real(0)) return Sign::Negative;
  return Sign::Zero;
}

inline void require_nonnull(const void* p, const char* name) {
  if (!p) {
    throw std::invalid_argument(std::string("spip::predicates::orient3d: null pointer for ") + name);
  }
}

inline void require_size3(std::size_t n, const char* name) {
  if (n != 3) {
    throw std::invalid_argument(std::string("spip::predicates::orient3d: ") + name + " must have size 3");
  }
}

template <class T>
inline constexpr bool supported_scalar_v =
    std::is_same_v<T, float> || std::is_same_v<T, double>;

template <class T>
inline void require_supported_and_matches_real() {
  static_assert(supported_scalar_v<T>, "spip::predicates::orient3d: only float/double are supported");
  static_assert(std::is_same_v<T, Real>,
                "spip::predicates::orient3d: scalar type must match SPIP_PREDICATES_REAL "
                "(define/un-define SPIP_PREDICATES_USE_FLOAT consistently)");
}

inline Sign orient3d_raw_ptr(const Real* a,
                             const Real* b,
                             const Real* o,
                             const Real* q) {
  require_nonnull(a, "a");
  require_nonnull(b, "b");
  require_nonnull(o, "o");
  require_nonnull(q, "q");

  ensure_exactinit();

  // Shewchuk expects non-const arrays; it does not mutate them.
  Real pa[3] = {a[0], a[1], a[2]};
  Real pb[3] = {b[0], b[1], b[2]};
  Real pc[3] = {o[0], o[1], o[2]};
  Real pd[3] = {q[0], q[1], q[2]};

  const Real v = orient3d(pa, pb, pc, pd);
  return sign_from_value(v);
}

inline Sign orient3d_on_sphere_raw_ptr(const Real* a,
                                       const Real* b,
                                       const Real* q) {
  const Real o[3] = {Real(0), Real(0), Real(0)};
  return orient3d_raw_ptr(a, b, o, q);
}

}  // namespace internal

// -----------------------------------------------------------------------------
// Public API: pointer core (best for broad compatibility)
// -----------------------------------------------------------------------------

template <class T>
inline Sign orient3d(const T* a, const T* b, const T* o, const T* q) {
  internal::require_supported_and_matches_real<T>();
  return internal::orient3d_raw_ptr(a, b, o, q);
}

/// Unit-sphere specialization (O is origin): sign( (A x B) · Q )
template <class T>
inline Sign orient3d_on_sphere(const T* a, const T* b, const T* q) {
  internal::require_supported_and_matches_real<T>();
  return internal::orient3d_on_sphere_raw_ptr(a, b, q);
}

// -----------------------------------------------------------------------------
// Public API: std::array overloads
// -----------------------------------------------------------------------------

template <class T>
inline Sign orient3d(const std::array<T, 3>& a,
                     const std::array<T, 3>& b,
                     const std::array<T, 3>& o,
                     const std::array<T, 3>& q) {
  internal::require_supported_and_matches_real<T>();
  return internal::orient3d_raw_ptr(a.data(), b.data(), o.data(), q.data());
}

template <class T>
inline Sign orient3d_on_sphere(const std::array<T, 3>& a,
                              const std::array<T, 3>& b,
                              const std::array<T, 3>& q) {
  internal::require_supported_and_matches_real<T>();
  return internal::orient3d_on_sphere_raw_ptr(a.data(), b.data(), q.data());
}

// -----------------------------------------------------------------------------
// Public API: std::vector overloads (expects size == 3)
// -----------------------------------------------------------------------------

template <class T>
inline Sign orient3d(const std::vector<T>& a,
                     const std::vector<T>& b,
                     const std::vector<T>& o,
                     const std::vector<T>& q) {
  internal::require_supported_and_matches_real<T>();
  internal::require_size3(a.size(), "a");
  internal::require_size3(b.size(), "b");
  internal::require_size3(o.size(), "o");
  internal::require_size3(q.size(), "q");
  return internal::orient3d_raw_ptr(a.data(), b.data(), o.data(), q.data());
}

template <class T>
inline Sign orient3d_on_sphere(const std::vector<T>& a,
                              const std::vector<T>& b,
                              const std::vector<T>& q) {
  internal::require_supported_and_matches_real<T>();
  internal::require_size3(a.size(), "a");
  internal::require_size3(b.size(), "b");
  internal::require_size3(q.size(), "q");
  return internal::orient3d_on_sphere_raw_ptr(a.data(), b.data(), q.data());
}

}  // namespace spip::predicates
