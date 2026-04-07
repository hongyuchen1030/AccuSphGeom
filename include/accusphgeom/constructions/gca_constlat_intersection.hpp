#pragma once

#include <array>
#include <stdexcept>

#include "accusphgeom/constructions/accucross.hpp"
#include "accusphgeom/predicates/on_minor_arc.hpp"

namespace accusphgeom::constructions {

template <typename T>
struct GcaConstLatIntersections {
  numeric::Vec3<T> point_pos{};
  numeric::Vec3<T> point_neg{};
  numeric::Vec3<T> normal_high{};
  numeric::Vec3<T> normal_low{};
};

template <typename T>
inline GcaConstLatIntersections<T> accux_constlat(const numeric::Vec3<T>& a,
                                                  const numeric::Vec3<T>& b,
                                                  T z0) {
  const auto normal = accucross(a, b);
  const auto s2 = numeric::sum_of_squares_c<T, 2>(
      {normal.hi[0], normal.hi[1]}, {normal.lo[0], normal.lo[1]});
  const auto s3 = numeric::sum_of_squares_c<T, 3>(normal.hi, normal.lo);
  const auto zsq = numeric::two_prod_fma(z0, z0);
  const auto d = numeric::compensated_dot_product(
      std::array<T, 4>{s3.hi, s3.hi, s3.lo, s3.lo},
      std::array<T, 4>{zsq.hi, zsq.lo, zsq.hi, zsq.lo});
  const auto e = numeric::two_sum(s2.hi, -d.hi);
  const T planar_sq = e.hi + (e.lo + s2.lo - d.lo);
  const auto s = numeric::acc_sqrt_re(planar_sq);

  const T nx = normal.hi[0] + normal.lo[0];
  const T ny = normal.hi[1] + normal.lo[1];
  const T nz = normal.hi[2] + normal.lo[2];
  const T planar = s.hi + s.lo;
  const T denom = s2.hi + s2.lo;
  if (denom == T(0)) {
    throw std::domain_error(
        "gca_constlat_intersection: degenerate great-circle normal");
  }

  const auto x_num_pos = numeric::compensated_dot_product(
      std::array<T, 2>{nx * nz, -ny},
      std::array<T, 2>{z0, planar});
  const auto y_num_pos = numeric::compensated_dot_product(
      std::array<T, 2>{ny * nz, nx},
      std::array<T, 2>{z0, planar});
  const auto x_num_neg = numeric::compensated_dot_product(
      std::array<T, 2>{nx * nz, ny},
      std::array<T, 2>{z0, planar});
  const auto y_num_neg = numeric::compensated_dot_product(
      std::array<T, 2>{ny * nz, -nx},
      std::array<T, 2>{z0, planar});

  GcaConstLatIntersections<T> out{};
  out.normal_high = normal.hi;
  out.normal_low = normal.lo;
  out.point_pos = {-(x_num_pos.hi + x_num_pos.lo) / denom,
                   -(y_num_pos.hi + y_num_pos.lo) / denom, z0};
  out.point_neg = {-(x_num_neg.hi + x_num_neg.lo) / denom,
                   -(y_num_neg.hi + y_num_neg.lo) / denom, z0};
  return out;
}

template <typename T>
inline numeric::Vec3<T> gca_constlat_intersection(const numeric::Vec3<T>& a,
                                                  const numeric::Vec3<T>& b,
                                                  T z0) {
  const auto candidates = accux_constlat(a, b, z0);
  const bool pos_on_arc = predicates::on_minor_arc(candidates.point_pos, a, b);
  const bool neg_on_arc = predicates::on_minor_arc(candidates.point_neg, a, b);

  if (pos_on_arc && !neg_on_arc) {
    return candidates.point_pos;
  }
  if (neg_on_arc && !pos_on_arc) {
    return candidates.point_neg;
  }
  if (pos_on_arc && neg_on_arc) {
    throw std::domain_error(
        "gca_constlat_intersection: both intersections lie on the minor arc");
  }
  throw std::domain_error(
      "gca_constlat_intersection: no intersection lies on the minor arc");
}

}  // namespace accusphgeom::constructions
