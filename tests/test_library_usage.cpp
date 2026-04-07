#include <accusphgeom/algorithms/point_in_polygon_sphere.hpp>

#include <array>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <type_traits>
#include <vector>

using accusphgeom::algorithms::Location;
using accusphgeom::algorithms::point_in_polygon_sphere;

namespace {

const char* to_str(Location x) {
  switch (x) {
    case Location::Outside: return "Outside";
    case Location::Inside: return "Inside";
    case Location::OnVertex: return "OnVertex";
    case Location::OnEdge: return "OnEdge";
  }
  return "Unknown";
}

void expect_equal(Location got, Location expected, const char* test_name) {
  if (got != expected) {
    std::cerr << "[FAIL] " << test_name << ": expected " << to_str(expected)
              << ", got " << to_str(got) << '\n';
    std::exit(EXIT_FAILURE);
  }
  std::cout << "[PASS] " << test_name << '\n';
}

template <typename GlobalId>
void expect_inside_with_global_id_type(
    const std::array<double, 3>& q,
    const std::array<std::array<double, 3>, 3>& poly,
    const char* test_name) {
  static_assert(std::is_integral_v<GlobalId>,
                "GlobalId must be an integral type");
  const std::array<GlobalId, 3> vertex_ids = {
      static_cast<GlobalId>(10),
      static_cast<GlobalId>(20),
      static_cast<GlobalId>(30),
  };
  expect_equal(point_in_polygon_sphere(q, poly, vertex_ids),
               Location::Inside,
               test_name);
}

}  // namespace

int main() {
  const std::int64_t q_id = 40;
  const std::int64_t r_id = 50;

  const std::array<double, 3> q_array = {
      0.5773502691896257, 0.5773502691896257, 0.5773502691896257};
  const std::array<double, 3> r_array = {
      -0.5773502691896257, -0.5773502691896257, -0.5773502691896257};
  const std::array<std::array<double, 3>, 3> poly_array = {{
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  }};
  const std::array<std::int64_t, 3> vertex_ids_array = {10, 20, 30};

  const std::vector<double> q_vec = {
      0.5773502691896257, 0.5773502691896257, 0.5773502691896257};
  const std::vector<double> r_vec = {
      -0.5773502691896257, -0.5773502691896257, -0.5773502691896257};
  const std::vector<std::vector<double>> poly_vec = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };
  const std::vector<std::int64_t> vertex_ids_vec = {10, 20, 30};

  double q_raw[3] = {0.5773502691896257, 0.5773502691896257, 0.5773502691896257};
  double r_raw[3] = {-0.5773502691896257, -0.5773502691896257,
                     -0.5773502691896257};
  double poly_storage[3][3] = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };
  const double* poly_raw[3] = {
      poly_storage[0], poly_storage[1], poly_storage[2]};
  const std::int64_t vertex_ids_raw[3] = {10, 20, 30};

  // robust / no-global-ID / raw pointer overload
  expect_equal(point_in_polygon_sphere(q_raw, poly_raw, 3),
               Location::Inside,
               "robust / no-global-ID / raw pointer overload");

  // robust / no-global-ID / std::array overload
  expect_equal(point_in_polygon_sphere(q_array, poly_array),
               Location::Inside,
               "robust / no-global-ID / std::array overload");

  // robust / no-global-ID / std::vector overload
  expect_equal(point_in_polygon_sphere(q_vec, poly_vec),
               Location::Inside,
               "robust / no-global-ID / std::vector overload");

  // robust / Tier 1 / raw pointer overload
  expect_equal(point_in_polygon_sphere(
                   q_raw, q_id, r_raw, r_id, poly_raw, vertex_ids_raw, 3),
               Location::Inside,
               "robust / Tier 1 / raw pointer overload");

  // robust / Tier 1 / std::array overload
  expect_equal(point_in_polygon_sphere(
                   q_array, q_id, r_array, r_id, poly_array, vertex_ids_array),
               Location::Inside,
               "robust / Tier 1 / std::array overload");

  // robust / Tier 1 / std::vector overload
  expect_equal(point_in_polygon_sphere(
                   q_vec, q_id, r_vec, r_id, poly_vec, vertex_ids_vec),
               Location::Inside,
               "robust / Tier 1 / std::vector overload");

  // robust / Tier 2 / raw pointer overload
  expect_equal(point_in_polygon_sphere(q_raw, q_id, poly_raw, vertex_ids_raw, 3),
               Location::Inside,
               "robust / Tier 2 / raw pointer overload");

  // robust / Tier 2 / std::array overload
  expect_equal(point_in_polygon_sphere(q_array, q_id, poly_array, vertex_ids_array),
               Location::Inside,
               "robust / Tier 2 / std::array overload");

  // robust / Tier 2 / std::vector overload
  expect_equal(point_in_polygon_sphere(q_vec, q_id, poly_vec, vertex_ids_vec),
               Location::Inside,
               "robust / Tier 2 / std::vector overload");

  // robust / Tier 3 / raw pointer overload
  expect_equal(point_in_polygon_sphere(q_raw, poly_raw, vertex_ids_raw, 3),
               Location::Inside,
               "robust / Tier 3 / raw pointer overload");

  // robust / Tier 3 / std::array overload
  expect_equal(point_in_polygon_sphere(q_array, poly_array, vertex_ids_array),
               Location::Inside,
               "robust / Tier 3 / std::array overload");

  // robust / Tier 3 / std::vector overload
  expect_equal(point_in_polygon_sphere(q_vec, poly_vec, vertex_ids_vec),
               Location::Inside,
               "robust / Tier 3 / std::vector overload");

  // robust / Tier 3 / std::array overload / short global IDs
  expect_inside_with_global_id_type<short>(
      q_array,
      poly_array,
      "robust / Tier 3 / std::array overload / short global IDs");

  // robust / Tier 3 / std::array overload / int global IDs
  expect_inside_with_global_id_type<int>(
      q_array,
      poly_array,
      "robust / Tier 3 / std::array overload / int global IDs");

  // robust / Tier 3 / std::array overload / long global IDs
  expect_inside_with_global_id_type<long>(
      q_array,
      poly_array,
      "robust / Tier 3 / std::array overload / long global IDs");

  // robust / Tier 3 / std::array overload / long long global IDs
  expect_inside_with_global_id_type<long long>(
      q_array,
      poly_array,
      "robust / Tier 3 / std::array overload / long long global IDs");

  std::cout << "OK\n";
  return 0;
}
