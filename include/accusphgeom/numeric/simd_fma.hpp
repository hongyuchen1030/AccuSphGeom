#pragma once

#include <cmath>
#include <type_traits>

#if __has_include(<Eigen/Core>)
#include <Eigen/Core>
#define ACCUSPHGEOM_HAS_EIGEN 1
#else
#define ACCUSPHGEOM_HAS_EIGEN 0
#endif

namespace accusphgeom::numeric {

template <typename T>
inline std::enable_if_t<std::is_arithmetic_v<T>, T> simd_fma(T a, T b, T c) {
  return std::fma(a, b, c);
}

#if ACCUSPHGEOM_HAS_EIGEN
template <typename Derived>
inline std::enable_if_t<
    std::is_base_of_v<Eigen::DenseBase<Derived>, Derived>, Derived>
simd_fma(const Derived& a, const Derived& b, const Derived& c) {
  Derived out;
  for (int i = 0; i < a.size(); ++i) {
    out.derived().data()[i] =
        std::fma(a.derived().data()[i], b.derived().data()[i],
                 c.derived().data()[i]);
  }
  return out;
}
#endif

}  // namespace accusphgeom::numeric
