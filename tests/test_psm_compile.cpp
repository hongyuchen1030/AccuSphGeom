#include "MultiPrecision_psm.h"
#include "accusphgeom/algorithms/point_in_polygon_sphere.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <vector>

int main() {
  GEO::expansion_nt a(1.0), b(2.0);
  GEO::expansion_nt c = a * b - a;

  using accusphgeom::algorithms::Location;

  const std::array<double, 3> q = {1.0, 0.0, 0.0};
  const std::array<double, 3> r = {-1.0, -1.0, -1.0};
  const std::array<std::array<double, 3>, 3> tri = {{
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  }};
  using GlobalId = std::int64_t;
  const std::array<GlobalId, 3> ids = {10, 20, 30};
  const std::vector<std::array<double, 3>> tri_vec = {tri[0], tri[1], tri[2]};
  const std::vector<GlobalId> ids_vec = {10, 20, 30};
  const std::array<GlobalId, 3> overflow_ids = {
      std::numeric_limits<GlobalId>::max() - 2,
      std::numeric_limits<GlobalId>::max() - 1,
      std::numeric_limits<GlobalId>::max(),
  };
  const GlobalId q_id = 40;
  const GlobalId r_id = 50;

  const Location l0 = accusphgeom::algorithms::point_in_polygon_sphere(q, tri);
  const Location l1 =
      accusphgeom::algorithms::point_in_polygon_sphere(q, tri, ids);
  const Location l1_t2 =
      accusphgeom::algorithms::point_in_polygon_sphere(q, q_id, tri, ids);
  const Location l1_t1 =
      accusphgeom::algorithms::point_in_polygon_sphere(q, q_id, r, r_id, tri, ids);
  const Location l2 =
      accusphgeom::algorithms::point_in_polygon_sphere(q, tri_vec);
  const Location l3 =
      accusphgeom::algorithms::point_in_polygon_sphere(q, tri_vec, ids_vec);
  const Location l3_t2 =
      accusphgeom::algorithms::point_in_polygon_sphere(q, q_id, tri_vec, ids_vec);
  const Location l3_t1 = accusphgeom::algorithms::point_in_polygon_sphere(
      q, q_id, r, r_id, tri_vec, ids_vec);
  const double* tri_ptrs[3] = {tri_vec[0].data(), tri_vec[1].data(), tri_vec[2].data()};
  const Location l4 = accusphgeom::algorithms::point_in_polygon_sphere(
      q.data(), tri_ptrs, tri_vec.size());
  const Location l5 = accusphgeom::algorithms::point_in_polygon_sphere(
      q.data(), tri_ptrs, ids_vec.data(), tri_vec.size());
  const Location l6 = accusphgeom::algorithms::point_in_polygon_sphere(
      q.data(), q_id, tri_ptrs, ids_vec.data(), tri_vec.size());
  const Location l7 = accusphgeom::algorithms::point_in_polygon_sphere(
      q.data(), q_id, r.data(), r_id, tri_ptrs, ids_vec.data(), tri_vec.size());
  const Location l8 =
      accusphgeom::algorithms::point_in_polygon_sphere(q, tri, overflow_ids);

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
