#pragma once
#include <type_traits>
#include <array>
#include <stdexcept>

#include "accusphgeom/constructions/accucross.hpp"
#include "accusphgeom/numeric/mask.hpp"
#include "accusphgeom/predicates/on_minor_arc.hpp"

namespace accusphgeom::constructions {

namespace internal {

inline constexpr double gca_constlat_minor_arc_tol = 1e-8;

}  // namespace internal


template <typename T>
struct GcaConstLatIntersections {
  numeric::Vec3<T> point_pos{};
  numeric::Vec3<T> point_neg{};
};


template <typename T>
struct GcaConstLatTryResult {
  numeric::Vec3<T> point{};
  T status{};
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
  const T neg_d_hi = -d.hi;
  const auto e = numeric::two_sum(s2.hi, neg_d_hi);
  const T planar_sq = e.hi + (e.lo + s2.lo - d.lo);
  const auto s = numeric::acc_sqrt_re(planar_sq);

  const T nx = normal.hi[0] + normal.lo[0];
  const T ny = normal.hi[1] + normal.lo[1];
  const T nz = normal.hi[2] + normal.lo[2];
  const T planar = s.hi + s.lo;
  const T denom = s2.hi + s2.lo;

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
  out.point_pos = {-(x_num_pos.hi + x_num_pos.lo) / denom,
                   -(y_num_pos.hi + y_num_pos.lo) / denom, z0};
  out.point_neg = {-(x_num_neg.hi + x_num_neg.lo) / denom,
                   -(y_num_neg.hi + y_num_neg.lo) / denom, z0};
  return out;
}

template <typename T>
inline auto
try_gca_constlat_intersection(const numeric::Vec3<T>& a,
                              const numeric::Vec3<T>& b,
                              T z0) {
  const auto c = accux_constlat(a, b, z0);
  const T tol = T(internal::gca_constlat_minor_arc_tol);
  const T one = T(1);

  const auto pos_finite =
      numeric::isfinite_mask(c.point_pos[0]) *
      numeric::isfinite_mask(c.point_pos[1]);

  const auto neg_finite =
      numeric::isfinite_mask(c.point_neg[0]) *
      numeric::isfinite_mask(c.point_neg[1]);

  const auto pos_on_arc =
      pos_finite *
      predicates::internal::on_minor_arc_tol_ptr(c.point_pos.data(), a.data(),
                                                 b.data(), tol);

  const auto neg_on_arc =
      neg_finite *
      predicates::internal::on_minor_arc_tol_ptr(c.point_neg.data(), a.data(),
                                                 b.data(), tol);

  const auto pos_valid = pos_finite * pos_on_arc;
  const auto neg_valid = neg_finite * neg_on_arc;

  const auto pos_mask = pos_valid * (one - neg_valid);
  const auto neg_mask = neg_valid * (one - pos_valid);

  numeric::Vec3<T> out{};
  out[0] = pos_mask * c.point_pos[0] + neg_mask * c.point_neg[0];
  out[1] = pos_mask * c.point_pos[1] + neg_mask * c.point_neg[1];
  out[2] = pos_mask * c.point_pos[2] + neg_mask * c.point_neg[2];

  const auto both = pos_valid * neg_valid;
  const auto none = (one - pos_valid) * (one - neg_valid);
  const auto status = both + none * T(2);

  return GcaConstLatTryResult<T>{out, status};
}

template <typename T>
inline numeric::Vec3<T>
gca_constlat_intersection(const numeric::Vec3<T>& a,
                          const numeric::Vec3<T>& b,
                          T z0) {
  static_assert(std::is_arithmetic_v<T>,
                "gca_constlat_intersection is scalar-only; use "
                "try_gca_constlat_intersection for packed/SIMD types.");

  const auto r = try_gca_constlat_intersection(a, b, z0);

  if (r.status == T(0)) return r.point;

  if (r.status == T(1)) {
    throw std::domain_error("both intersections lie on the minor arc");
  }

  throw std::domain_error("no intersection lies on the minor arc");
}

}  // namespace accusphgeom::constructions
