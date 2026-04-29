#pragma once

#include <cmath>

namespace accusphgeom::numeric {

template <typename T>
inline T numeric_abs(const T& x) {
  using std::abs;
  return abs(x);
}

}  // namespace accusphgeom::numeric
