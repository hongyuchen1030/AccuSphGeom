#include "accusphgeom/algorithms/point_in_polygon_sphere.hpp"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

using accusphgeom::algorithms::Location;
using accusphgeom::algorithms::point_in_polygon_sphere;

namespace {

void expect(Location got, Location expected, const char* name) {
  if (got != expected) {
    std::cerr << "[FAIL] " << name << '\n';
    std::exit(EXIT_FAILURE);
  }
  std::cout << "[PASS] " << name << '\n';
}

}  // namespace

int main() {
  const std::vector<std::array<double, 3>> polygon = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };
  const std::array<double, 3> q_vertex = polygon[0];
  const std::array<double, 3> q_edge = {0.70710678118654752,
                                        0.70710678118654752, 0.0};
  const std::array<double, 3> q_inside = {0.5773502691896257,
                                          0.5773502691896257,
                                          0.5773502691896257};
  const std::array<double, 3> r_inside = {-0.5773502691896257,
                                          -0.5773502691896257,
                                          -0.5773502691896257};
  const std::vector<std::int64_t> ids = {10, 20, 30};
  const std::vector<std::int64_t> overflow_ids = {
      std::numeric_limits<std::int64_t>::max() - 2,
      std::numeric_limits<std::int64_t>::max() - 1,
      std::numeric_limits<std::int64_t>::max(),
  };

  expect(point_in_polygon_sphere(q_vertex, polygon), Location::OnVertex,
         "vertex classification");
  expect(point_in_polygon_sphere(q_edge, polygon), Location::OnEdge,
         "edge classification");
  expect(point_in_polygon_sphere(q_inside, polygon), Location::Inside,
         "inside classification");
  expect(point_in_polygon_sphere(q_inside, polygon, ids), Location::Inside,
         "tier 3 classification");
  expect(point_in_polygon_sphere(q_inside, 40, polygon, ids), Location::Inside,
         "tier 2 classification");
  expect(point_in_polygon_sphere(q_inside, 40, r_inside, 50, polygon, ids),
         Location::Inside, "tier 1 classification");
  expect(point_in_polygon_sphere(q_inside, polygon, overflow_ids),
         Location::Inside, "tier 3 overflow ids");
  expect(point_in_polygon_sphere(q_inside, 7, polygon, overflow_ids),
         Location::Inside, "tier 2 overflow ids");
  return EXIT_SUCCESS;
}
