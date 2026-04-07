#pragma once

#include <array>
#include <stdexcept>

#include "accusphgeom/constructions/accucross.hpp"
#include "accusphgeom/predicates/on_minor_arc.hpp"

namespace accusphgeom::constructions {

template <typename T>
struct GcaGcaIntersections {
  numeric::Vec3<T> point_pos{};
  numeric::Vec3<T> point_neg{};
};

template <typename T>
inline GcaGcaIntersections<T> accux_gca(const numeric::Vec3<T>& a0,
                                        const numeric::Vec3<T>& a1,
                                        const numeric::Vec3<T>& b0,
                                        const numeric::Vec3<T>& b1) {
  const auto n1 = accucross(a0, a1);
  const auto n2 = accucross(b0, b1);
  const auto v = accucross(n1.hi, n1.lo, n2.hi, n2.lo);
  const auto sum = numeric::sum_of_squares_c<T, 3>(v.hi, v.lo);
  const auto norm = numeric::acc_sqrt_re(sum.hi, sum.lo);
  const T n = norm.hi + norm.lo;
  if (n == T(0)) {
    throw std::domain_error(
        "gca_gca_intersection: great circles are degenerate or coincident");
  }

  const numeric::Vec3<T> point_pos = {(v.hi[0] + v.lo[0]) / n,
                                      (v.hi[1] + v.lo[1]) / n,
                                      (v.hi[2] + v.lo[2]) / n};
  const numeric::Vec3<T> point_neg = {-point_pos[0], -point_pos[1],
                                      -point_pos[2]};
  return {point_pos, point_neg};
}

template <typename T>
inline numeric::Vec3<T> gca_gca_intersection(const numeric::Vec3<T>& a0,
                                             const numeric::Vec3<T>& a1,
                                             const numeric::Vec3<T>& b0,
                                             const numeric::Vec3<T>& b1) {
  const auto candidates = accux_gca(a0, a1, b0, b1);
  const bool pos_on_minor_arcs =
      predicates::on_minor_arc(candidates.point_pos, a0, a1) &&
      predicates::on_minor_arc(candidates.point_pos, b0, b1);
  const bool neg_on_minor_arcs =
      predicates::on_minor_arc(candidates.point_neg, a0, a1) &&
      predicates::on_minor_arc(candidates.point_neg, b0, b1);

  if (pos_on_minor_arcs && !neg_on_minor_arcs) {
    return candidates.point_pos;
  }
  if (neg_on_minor_arcs && !pos_on_minor_arcs) {
    return candidates.point_neg;
  }
  if (pos_on_minor_arcs && neg_on_minor_arcs) {
    throw std::domain_error(
        "gca_gca_intersection: both intersections lie on both minor arcs");
  }
  throw std::domain_error(
      "gca_gca_intersection: no intersection lies on both minor arcs");
}

}  // namespace accusphgeom::constructions
