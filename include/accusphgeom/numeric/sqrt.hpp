#pragma once

#include <cmath>

namespace accusphgeom::numeric {

template <typename T>
inline T numeric_sqrt(const T& x) {
  using std::sqrt;
  return sqrt(x);
}

}  // namespace accusphgeom::numeric