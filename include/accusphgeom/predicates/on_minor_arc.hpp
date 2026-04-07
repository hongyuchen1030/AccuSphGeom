#pragma once

#include <array>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include "accusphgeom/predicates/orient3d.hpp"
#include "accusphgeom/predicates/quadruple3d.hpp"

namespace accusphgeom::predicates {

namespace internal {

template <class T>
inline void require_nonnull_on_minor_arc(const T* p, const char* name) {
  if (!p) {
    throw std::invalid_argument(
        std::string("accusphgeom::predicates::on_minor_arc: null pointer for ") +
        name);
  }
}

template <class T>
inline bool equal3_exact(const T* a, const T* b) {
  require_nonnull_on_minor_arc(a, "a");
  require_nonnull_on_minor_arc(b, "b");
  return (a[0] == b[0]) && (a[1] == b[1]) && (a[2] == b[2]);
}

template <class T>
inline bool on_minor_arc_raw_ptr(const T* q, const T* a, const T* b) {
  require_nonnull_on_minor_arc(q, "q");
  require_nonnull_on_minor_arc(a, "a");
  require_nonnull_on_minor_arc(b, "b");

  if (equal3_exact(a, b)) {
    return false;
  }
  if (orient3d_on_sphere(a, b, q) != Sign::Zero) {
    return false;
  }
  return quadruple3d(a, q, q, b) != Sign::Negative;
}

}  // namespace internal

template <class T>
inline bool on_minor_arc(const T* q, const T* a, const T* b) {
  internal::require_supported_and_matches_real<T>();
  return internal::on_minor_arc_raw_ptr(q, a, b);
}

template <class T>
inline bool on_minor_arc(const std::array<T, 3>& q,
                         const std::array<T, 3>& a,
                         const std::array<T, 3>& b) {
  internal::require_supported_and_matches_real<T>();
  return internal::on_minor_arc_raw_ptr(q.data(), a.data(), b.data());
}

template <class T>
inline bool on_minor_arc(const std::vector<T>& q,
                         const std::vector<T>& a,
                         const std::vector<T>& b) {
  internal::require_supported_and_matches_real<T>();
  internal::require_size3(q.size(), "q");
  internal::require_size3(a.size(), "a");
  internal::require_size3(b.size(), "b");
  return internal::on_minor_arc_raw_ptr(q.data(), a.data(), b.data());
}

}  // namespace accusphgeom::predicates
