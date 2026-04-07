#pragma once

#include <cmath>
#include <type_traits>

namespace accusphgeom::numeric {

template <typename T>
inline std::enable_if_t<std::is_arithmetic_v<T>, T> simd_fma(T a, T b, T c) {
  return std::fma(a, b, c);
}

}  // namespace accusphgeom::numeric
