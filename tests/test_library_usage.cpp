#include <spip/algorithms/point_in_polygon_sphere.hpp>

#include <array>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <vector>

using spip::V3_T;
using spip::pip::Location;
using spip::pip::point_in_polygon_sphere;

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

  const V3_T<double> q_eft(
      0.5773502691896257, 0.5773502691896257, 0.5773502691896257);
  const V3_T<double> r_eft(
      -0.5773502691896257, -0.5773502691896257, -0.5773502691896257);
  const std::array<V3_T<double>, 3> poly_eft_array = {
      V3_T<double>(1.0, 0.0, 0.0),
      V3_T<double>(0.0, 1.0, 0.0),
      V3_T<double>(0.0, 0.0, 1.0),
  };
  const std::vector<V3_T<double>> poly_eft_vec = {
      V3_T<double>(1.0, 0.0, 0.0),
      V3_T<double>(0.0, 1.0, 0.0),
      V3_T<double>(0.0, 0.0, 1.0),
  };
  const V3_T<double> poly_eft_raw[3] = {
      V3_T<double>(1.0, 0.0, 0.0),
      V3_T<double>(0.0, 1.0, 0.0),
      V3_T<double>(0.0, 0.0, 1.0),
  };

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

  // EFT / no-global-ID / raw pointer overload
  expect_equal(point_in_polygon_sphere(q_eft, poly_eft_raw, 3),
               Location::Inside,
               "EFT / no-global-ID / raw pointer overload");

  // EFT / no-global-ID / std::array overload
  expect_equal(point_in_polygon_sphere(q_eft, poly_eft_array),
               Location::Inside,
               "EFT / no-global-ID / std::array overload");

  // EFT / no-global-ID / std::vector overload
  expect_equal(point_in_polygon_sphere(q_eft, poly_eft_vec),
               Location::Inside,
               "EFT / no-global-ID / std::vector overload");

  // EFT / Tier 1 / raw pointer overload
  expect_equal(point_in_polygon_sphere(
                   q_eft, q_id, r_eft, r_id, poly_eft_raw, vertex_ids_raw, 3),
               Location::Inside,
               "EFT / Tier 1 / raw pointer overload");

  // EFT / Tier 1 / std::array overload
  expect_equal(point_in_polygon_sphere(
                   q_eft, q_id, r_eft, r_id, poly_eft_array, vertex_ids_array),
               Location::Inside,
               "EFT / Tier 1 / std::array overload");

  // EFT / Tier 1 / std::vector overload
  expect_equal(point_in_polygon_sphere(
                   q_eft, q_id, r_eft, r_id, poly_eft_vec, vertex_ids_vec),
               Location::Inside,
               "EFT / Tier 1 / std::vector overload");

  // EFT / Tier 2 / raw pointer overload
  expect_equal(point_in_polygon_sphere(q_eft, q_id, poly_eft_raw, vertex_ids_raw, 3),
               Location::Inside,
               "EFT / Tier 2 / raw pointer overload");

  // EFT / Tier 2 / std::array overload
  expect_equal(point_in_polygon_sphere(q_eft, q_id, poly_eft_array, vertex_ids_array),
               Location::Inside,
               "EFT / Tier 2 / std::array overload");

  // EFT / Tier 2 / std::vector overload
  expect_equal(point_in_polygon_sphere(q_eft, q_id, poly_eft_vec, vertex_ids_vec),
               Location::Inside,
               "EFT / Tier 2 / std::vector overload");

  // EFT / Tier 3 / raw pointer overload
  expect_equal(point_in_polygon_sphere(q_eft, poly_eft_raw, vertex_ids_raw, 3),
               Location::Inside,
               "EFT / Tier 3 / raw pointer overload");

  // EFT / Tier 3 / std::array overload
  expect_equal(point_in_polygon_sphere(q_eft, poly_eft_array, vertex_ids_array),
               Location::Inside,
               "EFT / Tier 3 / std::array overload");

  // EFT / Tier 3 / std::vector overload
  expect_equal(point_in_polygon_sphere(q_eft, poly_eft_vec, vertex_ids_vec),
               Location::Inside,
               "EFT / Tier 3 / std::vector overload");

  std::cout << "OK\n";
  return 0;
}
