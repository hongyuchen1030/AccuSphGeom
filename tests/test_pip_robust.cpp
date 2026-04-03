#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
#include <array>

#include "accusphgeom/algorithms/point_in_polygon_sphere.hpp"

namespace {

using accusphgeom::algorithms::Location;
using accusphgeom::algorithms::point_in_polygon_sphere;

std::array<double, 3> normalize(double x, double y, double z) {
  const double n = std::sqrt(x * x + y * y + z * z);
  if (n == 0.0) {
    throw std::invalid_argument("normalize: zero vector");
  }
  return {x / n, y / n, z / n};
}

void expect_equal(Location got,
                  Location expected,
                  const char* test_name) {
  if (got != expected) {
    std::cerr << "[FAIL] " << test_name
              << ": expected " << static_cast<int>(expected)
              << ", got " << static_cast<int>(got) << '\n';
    std::exit(EXIT_FAILURE);
  }
  std::cout << "[PASS] " << test_name << '\n';
}

}  // namespace

int main() {
  // Spherical triangle:
  //   A = (1,0,0)
  //   B = (0,1,0)
  //   C = (0,0,1)
  //
  // Edge AB is the minor arc on the equator (z = 0).
  const std::vector<std::array<double, 3>> poly = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };

  // Scenario 1: query point exactly on a vertex.
  const std::array<double, 3> q_vertex = poly[0];

    // Scenario 2: query point exactly on the equatorial edge AB.
    // Computed in Wolfram Mathematica with adaptive precision and
    // rounded to 17 decimal digits.
    //
    // q = Normalize[{1,1,0}]
    //
    const std::array<double, 3> q_edge = {
    0.70710678118654752,
    0.70710678118654752,
    0.0
    };

  // Scenario 3: query point strictly inside the spherical triangle.
  // This point is not on any edge or vertex. For the standard ray-crossing
  // construction, it yields a genuine crossing configuration.
  const std::array<double, 3> q_inside = normalize(1.0, 1.0, 1.0);

  const Location loc_vertex = point_in_polygon_sphere(q_vertex, poly);
  const Location loc_edge = point_in_polygon_sphere(q_edge, poly);
  const Location loc_inside = point_in_polygon_sphere(q_inside, poly);
  const std::vector<std::int64_t> vertex_ids = {10, 20, 30};
  const std::array<double, 3> r_outside = normalize(-1.0, -1.0, -1.0);
  const Location loc_inside_tier3 =
      point_in_polygon_sphere(q_inside, poly, vertex_ids);
  const Location loc_inside_tier2 =
      point_in_polygon_sphere(q_inside, 40, poly, vertex_ids);
  const Location loc_inside_tier1 =
      point_in_polygon_sphere(q_inside, 40, r_outside, 50, poly, vertex_ids);
  const std::vector<std::int64_t> overflow_vertex_ids = {
      std::numeric_limits<std::int64_t>::max() - 2,
      std::numeric_limits<std::int64_t>::max() - 1,
      std::numeric_limits<std::int64_t>::max(),
  };
  const Location loc_inside_overflow_tier3 =
      point_in_polygon_sphere(q_inside, poly, overflow_vertex_ids);
  const Location loc_inside_overflow_tier2 =
      point_in_polygon_sphere(q_inside, 7, poly, overflow_vertex_ids);

  expect_equal(loc_vertex, Location::OnVertex,
               "adaptive/no-global-id: query on vertex");
  expect_equal(loc_edge, Location::OnEdge,
               "adaptive/no-global-id: query on edge");
  expect_equal(loc_inside, Location::Inside,
               "adaptive/no-global-id: query strictly inside");
  expect_equal(loc_inside_tier3, Location::Inside,
               "adaptive/tier3: query strictly inside");
  expect_equal(loc_inside_tier2, Location::Inside,
               "adaptive/tier2: query strictly inside");
  expect_equal(loc_inside_tier1, Location::Inside,
               "adaptive/tier1: query strictly inside");
  expect_equal(loc_inside_overflow_tier3, Location::Inside,
               "adaptive/tier3 overflow->mex: query strictly inside");
  expect_equal(loc_inside_overflow_tier2, Location::Inside,
               "adaptive/tier2 overflow->mex: query strictly inside");

  return EXIT_SUCCESS;
}
