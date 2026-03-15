#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <array>

#include "spip/algorithms/point_in_polygon_sphere.hpp"

namespace {

using spip::pip::Location;
using spip::pip::point_in_polygon_sphere;

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

  expect_equal(loc_vertex, Location::OnVertex,
               "adaptive/no-global-id: query on vertex");
  expect_equal(loc_edge, Location::OnEdge,
               "adaptive/no-global-id: query on edge");
  expect_equal(loc_inside, Location::Inside,
               "adaptive/no-global-id: query strictly inside");

  return EXIT_SUCCESS;
}