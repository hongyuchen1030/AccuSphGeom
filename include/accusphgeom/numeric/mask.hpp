#pragma once

#include <cmath>
#include <type_traits>

namespace accusphgeom::numeric {

template <typename T>
inline std::enable_if_t<std::is_arithmetic_v<T>, T>
isfinite_mask(const T& x) {
  return std::isfinite(x) ? T(1) : T(0);
}

template <typename T>
inline std::enable_if_t<std::is_arithmetic_v<T>, T>
mask_equal(const T& a, const T& b) {
  return a == b ? T(1) : T(0);
}

template <typename T>
inline std::enable_if_t<std::is_arithmetic_v<T>, T>
mask_le(const T& a, const T& b) {
  return a <= b ? T(1) : T(0);
}

template <typename T>
inline std::enable_if_t<std::is_arithmetic_v<T>, T>
mask_ge(const T& a, const T& b) {
  return a >= b ? T(1) : T(0);
}

}  // namespace accusphgeom::numeric
