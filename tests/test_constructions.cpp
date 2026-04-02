#include "accusphgeom/accusphgeom.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>

using accusphgeom::constructions::accucross;
using accusphgeom::constructions::gca_constlat_intersection;
using accusphgeom::constructions::gca_gca_intersection;
using accusphgeom::numeric::Vec3;

namespace {

bool approx(double a, double b, double tol = 1e-12) {
  return std::fabs(a - b) <= tol;
}

void expect(bool ok, const char* name) {
  if (!ok) {
    std::cerr << "[FAIL] " << name << '\n';
    std::exit(EXIT_FAILURE);
  }
  std::cout << "[PASS] " << name << '\n';
}

}  // namespace

int main() {
  const Vec3<double> ex = {1.0, 0.0, 0.0};
  const Vec3<double> ey = {0.0, 1.0, 0.0};
  const Vec3<double> ez = {0.0, 0.0, 1.0};

  const auto cross = accucross(ex, ey);
  expect(approx(cross.hi[2] + cross.lo[2], 1.0), "accucross basic");

  const auto gca_lat = gca_constlat_intersection(ey, ez, 0.5, 1);
  expect(approx(gca_lat.point[2], 0.5), "gca const-lat z");

  const auto inter = gca_gca_intersection(ex, ey, ey, ez);
  expect(approx(std::fabs(inter[1]), 1.0), "gca-gca intersection");
  return EXIT_SUCCESS;
}
