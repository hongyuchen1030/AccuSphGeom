#include <spip/algorithms/point_in_polygon_sphere.hpp>

#include <array>
#include <cassert>
#include <iostream>
#include <vector>

using spip::pip::Location;
using spip::pip::point_in_polygon_sphere;

static const char* to_str(Location x) {
  switch (x) {
    case Location::Outside:  return "Outside";
    case Location::Inside:   return "Inside";
    case Location::OnVertex: return "OnVertex";
    case Location::OnEdge:   return "OnEdge";
  }
  return "Unknown";
}

int main() {
  // A simple CCW spherical triangle on the unit sphere (not spanning > hemisphere).
  const std::vector<std::array<double, 3>> tri = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0},
  };

  // Query: "roughly inside" (normalize is not done here; assume inputs are unit as your code expects).
  const std::array<double, 3> q_inside = {0.5773502691896257,
                                          0.5773502691896257,
                                          0.5773502691896257};

  // Query: exactly a vertex
  const std::array<double, 3> q_vertex = {1.0, 0.0, 0.0};

  Location r1 = point_in_polygon_sphere(q_inside, tri);
  Location r2 = point_in_polygon_sphere(q_vertex, tri);

  std::cout << "q_inside: " << to_str(r1) << "\n";
  std::cout << "q_vertex: " << to_str(r2) << "\n";

  // These asserts assume your current implementation classifies these correctly.
  // If the triangle orientation or assumptions differ, adjust expected values.
  assert(r2 == Location::OnVertex);
  assert(r1 == Location::Inside || r1 == Location::OnEdge);

  std::cout << "OK\n";
  return 0;
}