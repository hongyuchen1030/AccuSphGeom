#pragma once

#include <array>
#include <cstddef>
#include <cmath>
#include <vector>

#include "MultiPrecision_psm.h"
#include "accusphgeom/predicates/orient3d.hpp"
#include "quadruple3d.h"

namespace accusphgeom::predicates {

namespace detail {

Sign quadruple3d_exact_fallback(const double* a, const double* b,
                                const double* c, const double* d);

template <typename A, typename B, typename C, typename D>
inline Sign quadruple3d_ptr(const A* a, const B* b, const C* c, const D* d) {
  const std::array<double, 3> aa = cast3(a);
  const std::array<double, 3> bb = cast3(b);
  const std::array<double, 3> cc = cast3(c);
  const std::array<double, 3> dd = cast3(d);

  const int sign = quadruple_3d_filter(aa.data(), bb.data(), cc.data(), dd.data());
  if (sign > 0) {
    return Sign::Positive;
  }
  if (sign < 0) {
    return Sign::Negative;
  }
  if (sign != FPG_UNCERTAIN_VALUE) {
    return Sign::Zero;
  }
  return quadruple3d_exact_fallback(aa.data(), bb.data(), cc.data(), dd.data());
}

}  // namespace detail

template <typename A, typename B, typename C, typename D>
inline Sign quadruple3d(const A* a, const B* b, const C* c, const D* d) {
  return detail::quadruple3d_ptr(a, b, c, d);
}

template <typename A, typename B, typename C, typename D>
inline Sign quadruple3d(const std::array<A, 3>& a, const std::array<B, 3>& b,
                        const std::array<C, 3>& c, const std::array<D, 3>& d) {
  return quadruple3d(a.data(), b.data(), c.data(), d.data());
}

template <typename A, typename B, typename C, typename D>
inline Sign quadruple3d(const std::vector<A>& a, const std::vector<B>& b,
                        const std::vector<C>& c, const std::vector<D>& d) {
  detail::require_size3(a.size(), "a");
  detail::require_size3(b.size(), "b");
  detail::require_size3(c.size(), "c");
  detail::require_size3(d.size(), "d");
  return quadruple3d(a.data(), b.data(), c.data(), d.data());
}

}  // namespace accusphgeom::predicates
