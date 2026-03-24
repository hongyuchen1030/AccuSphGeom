#include "MultiPrecision_psm.h"
#include "spip/algorithms/point_in_polygon_sphere.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <vector>

int main() {
  GEO::expansion_nt a(1.0), b(2.0);
  GEO::expansion_nt c = a*b - a;

  using spip::V3_T;
  using spip::pip::Location;

  const V3_T<double> q(1.0, 0.0, 0.0);
  const V3_T<double> r(-1.0, -1.0, -1.0);
  const std::array<V3_T<double>, 3> tri = {
      V3_T<double>(1.0, 0.0, 0.0),
      V3_T<double>(0.0, 1.0, 0.0),
      V3_T<double>(0.0, 0.0, 1.0),
  };
  const std::array<std::int64_t, 3> ids = {10, 20, 30};
  const std::vector<V3_T<double>> tri_vec = {tri[0], tri[1], tri[2]};
  const std::vector<std::int64_t> ids_vec = {10, 20, 30};
  const std::array<std::int64_t, 3> overflow_ids = {
      std::numeric_limits<std::int64_t>::max() - 2,
      std::numeric_limits<std::int64_t>::max() - 1,
      std::numeric_limits<std::int64_t>::max(),
  };

  const Location l0 = spip::pip::point_in_polygon_sphere(q, tri);
  const Location l1 = spip::pip::point_in_polygon_sphere(q, tri, ids);
  const Location l1_t2 = spip::pip::point_in_polygon_sphere(q, 40, tri, ids);
  const Location l1_t1 =
      spip::pip::point_in_polygon_sphere(q, 40, r, 50, tri, ids);
  const Location l2 = spip::pip::point_in_polygon_sphere(q, tri_vec);
  const Location l3 = spip::pip::point_in_polygon_sphere(q, tri_vec, ids_vec);
  const Location l3_t2 =
      spip::pip::point_in_polygon_sphere(q, 40, tri_vec, ids_vec);
  const Location l3_t1 =
      spip::pip::point_in_polygon_sphere(q, 40, r, 50, tri_vec, ids_vec);
  const Location l4 = spip::pip::point_in_polygon_sphere(
      q, tri_vec.data(), tri_vec.size());
  const Location l5 = spip::pip::point_in_polygon_sphere(
      q, tri_vec.data(), ids_vec.data(), tri_vec.size());
  const Location l6 = spip::pip::point_in_polygon_sphere(
      q, 40, tri_vec.data(), ids_vec.data(), tri_vec.size());
  const Location l7 = spip::pip::point_in_polygon_sphere(
      q, 40, r, 50, tri_vec.data(), ids_vec.data(), tri_vec.size());
  const Location l8 = spip::pip::point_in_polygon_sphere(q, tri, overflow_ids);

  (void)l0;
  (void)l1;
  (void)l1_t2;
  (void)l1_t1;
  (void)l2;
  (void)l3;
  (void)l3_t2;
  (void)l3_t1;
  (void)l4;
  (void)l5;
  (void)l6;
  (void)l7;
  (void)l8;

  return (c.sign() == GEO::POSITIVE) ? 0 : 1;
}
