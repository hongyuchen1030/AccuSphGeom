#include "accusphgeom/predicates/quadruple3d.hpp"

#include <array>
#include <cstdlib>
#include <iostream>

using accusphgeom::predicates::Sign;

namespace {

void expect(Sign got, Sign expected, const char* name) {
  if (got != expected) {
    std::cerr << "[FAIL] " << name << '\n';
    std::exit(EXIT_FAILURE);
  }
  std::cout << "[PASS] " << name << '\n';
}

}  // namespace

int main() {
  const std::array<double, 3> e1 = {1.0, 0.0, 0.0};
  const std::array<double, 3> e2 = {0.0, 1.0, 0.0};

  expect(accusphgeom::predicates::quadruple3d(e1, e2, e1, e2), Sign::Positive,
         "positive sign");
  expect(accusphgeom::predicates::quadruple3d(e1, e2, e2, e1), Sign::Negative,
         "negative sign");
  return EXIT_SUCCESS;
}
