#include "accusphgeom/accusphgeom.hpp"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

using accusphgeom::algorithms::Location;
using accusphgeom::algorithms::point_in_polygon_sphere;

namespace {

const char* to_string(Location x) {
  switch (x) {
    case Location::Outside:
      return "Outside";
    case Location::Inside:
      return "Inside";
    case Location::OnVertex:
      return "OnVertex";
    case Location::OnEdge:
      return "OnEdge";
  }
  return "Unknown";
}

void expect(Location got, Location expected, const char* name) {
  if (got != expected) {
    std::cerr << "[FAIL] " << name << ": expected " << to_string(expected)
              << ", got " << to_string(got) << '\n';
    std::exit(EXIT_FAILURE);
  }
  std::cout << "[PASS] " << name << '\n';
}

}  // namespace

int main() {
  const std::array<double, 3> q = {
      0.5773502691896257, 0.5773502691896257, 0.5773502691896257};
  const std::array<double, 3> r = {
      -0.5773502691896257, -0.5773502691896257, -0.5773502691896257};
  const std::array<std::array<double, 3>, 3> polygon = {{
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  }};
  const std::array<std::int64_t, 3> ids = {10, 20, 30};

  double q_raw[3] = {q[0], q[1], q[2]};
  double r_raw[3] = {r[0], r[1], r[2]};
  double polygon_storage[3][3] = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };
  const double* polygon_raw[3] = {
      polygon_storage[0], polygon_storage[1], polygon_storage[2]};
  const std::int64_t ids_raw[3] = {10, 20, 30};

  const std::vector<double> q_vec = {q[0], q[1], q[2]};
  const std::vector<double> r_vec = {r[0], r[1], r[2]};
  const std::vector<std::vector<double>> polygon_vec = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };
  const std::vector<std::int64_t> ids_vec = {10, 20, 30};

  expect(point_in_polygon_sphere(q_raw, polygon_raw, 3), Location::Inside,
         "raw / no global ids");
  expect(point_in_polygon_sphere(q, polygon), Location::Inside,
         "array / no global ids");
  expect(point_in_polygon_sphere(q_vec, polygon_vec), Location::Inside,
         "vector / no global ids");

  expect(point_in_polygon_sphere(q_raw, polygon_raw, ids_raw, 3),
         Location::Inside, "raw / tier 3");
  expect(point_in_polygon_sphere(q, polygon, ids), Location::Inside,
         "array / tier 3");
  expect(point_in_polygon_sphere(q_vec, polygon_vec, ids_vec), Location::Inside,
         "vector / tier 3");

  expect(point_in_polygon_sphere(q_raw, 40, polygon_raw, ids_raw, 3),
         Location::Inside, "raw / tier 2");
  expect(point_in_polygon_sphere(q, 40, polygon, ids), Location::Inside,
         "array / tier 2");
  expect(point_in_polygon_sphere(q_vec, 40, polygon_vec, ids_vec),
         Location::Inside, "vector / tier 2");

  expect(point_in_polygon_sphere(q_raw, 40, r_raw, 50, polygon_raw, ids_raw, 3),
         Location::Inside, "raw / tier 1");
  expect(point_in_polygon_sphere(q, 40, r, 50, polygon, ids), Location::Inside,
         "array / tier 1");
  expect(point_in_polygon_sphere(q_vec, 40, r_vec, 50, polygon_vec, ids_vec),
         Location::Inside, "vector / tier 1");

  return EXIT_SUCCESS;
}
