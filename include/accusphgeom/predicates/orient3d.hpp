#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#include "predicates.h"

namespace accusphgeom::predicates {

enum class Sign : std::int8_t { Negative = -1, Zero = 0, Positive = 1 };

namespace detail {

inline void ensure_exactinit() {
  static std::once_flag flag;
  std::call_once(flag, []() { exactinit(); });
}

inline void require_nonnull(const void* p, const char* name) {
  if (!p) {
    throw std::invalid_argument(std::string("accusphgeom::orient3d: null ") +
                                name);
  }
}

inline void require_size3(std::size_t n, const char* name) {
  if (n != 3) {
    throw std::invalid_argument(std::string("accusphgeom::orient3d: ") + name +
                                " must have size 3");
  }
}

inline Sign sign_from_value(double value) {
  if (value > 0.0) {
    return Sign::Positive;
  }
  if (value < 0.0) {
    return Sign::Negative;
  }
  return Sign::Zero;
}

template <typename T>
inline std::array<double, 3> cast3(const T* p) {
  require_nonnull(p, "pointer");
  return {static_cast<double>(p[0]), static_cast<double>(p[1]),
          static_cast<double>(p[2])};
}

template <typename A, typename B, typename O, typename Q>
inline Sign orient3d_ptr(const A* a, const B* b, const O* o, const Q* q) {
  ensure_exactinit();
  const std::array<double, 3> pa = cast3(a);
  const std::array<double, 3> pb = cast3(b);
  const std::array<double, 3> po = cast3(o);
  const std::array<double, 3> pq = cast3(q);
  double aa[3] = {pa[0], pa[1], pa[2]};
  double bb[3] = {pb[0], pb[1], pb[2]};
  double oo[3] = {po[0], po[1], po[2]};
  double qq[3] = {pq[0], pq[1], pq[2]};
  return sign_from_value(::orient3d(aa, bb, oo, qq));
}

}  // namespace detail

template <typename A, typename B, typename O, typename Q>
inline Sign orient3d(const A* a, const B* b, const O* o, const Q* q) {
  return detail::orient3d_ptr(a, b, o, q);
}

template <typename A, typename B, typename Q>
inline Sign orient3d_on_sphere(const A* a, const B* b, const Q* q) {
  const double origin[3] = {0.0, 0.0, 0.0};
  return detail::orient3d_ptr(a, b, origin, q);
}

template <typename A, typename B, typename O, typename Q>
inline Sign orient3d(const std::array<A, 3>& a, const std::array<B, 3>& b,
                     const std::array<O, 3>& o, const std::array<Q, 3>& q) {
  return orient3d(a.data(), b.data(), o.data(), q.data());
}

template <typename A, typename B, typename Q>
inline Sign orient3d_on_sphere(const std::array<A, 3>& a,
                               const std::array<B, 3>& b,
                               const std::array<Q, 3>& q) {
  return orient3d_on_sphere(a.data(), b.data(), q.data());
}

template <typename A, typename B, typename O, typename Q>
inline Sign orient3d(const std::vector<A>& a, const std::vector<B>& b,
                     const std::vector<O>& o, const std::vector<Q>& q) {
  detail::require_size3(a.size(), "a");
  detail::require_size3(b.size(), "b");
  detail::require_size3(o.size(), "o");
  detail::require_size3(q.size(), "q");
  return orient3d(a.data(), b.data(), o.data(), q.data());
}

template <typename A, typename B, typename Q>
inline Sign orient3d_on_sphere(const std::vector<A>& a,
                               const std::vector<B>& b,
                               const std::vector<Q>& q) {
  detail::require_size3(a.size(), "a");
  detail::require_size3(b.size(), "b");
  detail::require_size3(q.size(), "q");
  return orient3d_on_sphere(a.data(), b.data(), q.data());
}

}  // namespace accusphgeom::predicates
