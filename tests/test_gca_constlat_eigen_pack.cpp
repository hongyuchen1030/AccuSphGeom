#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "accusphgeom/adapters/eigen/numeric.hpp"
#include "accusphgeom/constructions/gca_constlat_intersection.hpp"

namespace {

using accusphgeom::constructions::accux_constlat;
using accusphgeom::numeric::Vec3;
using accusphgeom::constructions::try_gca_constlat_intersection;

constexpr int N = 4;
using Pack = EigenPack<N>;
constexpr double kTol = 0.0;

bool close(double a, double b) {
  return std::abs(a - b) <= kTol;
}

}  // namespace

int main() {
  Eigen::Matrix<double, N, 3> ptsA;
  Eigen::Matrix<double, N, 3> ptsB;
  Eigen::Array<double, N, 1> z0;

  ptsA <<
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.7071067811865475, 0.7071067811865475, 0.0,
      0.0, 0.7071067811865475, 0.7071067811865475;

  ptsB <<
      0.0, 0.0, 1.0,
      0.0, 0.0, 1.0,
      0.0, 0.0, 1.0,
      1.0, 0.0, 0.0;

  z0 <<
      0.5,
      -0.25,
      0.125,
      -0.5;

  Pack ax = ptsA.template block<N, 1>(0, 0).array();
  Pack ay = ptsA.template block<N, 1>(0, 1).array();
  Pack az = ptsA.template block<N, 1>(0, 2).array();

  Pack bx = ptsB.template block<N, 1>(0, 0).array();
  Pack by = ptsB.template block<N, 1>(0, 1).array();
  Pack bz = ptsB.template block<N, 1>(0, 2).array();

  Pack zz = z0;

  Vec3<Pack> a{ax, ay, az};
  Vec3<Pack> b{bx, by, bz};

  const auto packed = accux_constlat(a, b, zz);

  for (int i = 0; i < N; ++i) {
    Vec3<double> ai{ptsA(i, 0), ptsA(i, 1), ptsA(i, 2)};
    Vec3<double> bi{ptsB(i, 0), ptsB(i, 1), ptsB(i, 2)};
    const auto scalar = accux_constlat(ai, bi, z0(i));

    if (!close(packed.point_pos[0](i), scalar.point_pos[0]) ||
        !close(packed.point_pos[1](i), scalar.point_pos[1]) ||
        !close(packed.point_pos[2](i), scalar.point_pos[2]) ||
        !close(packed.point_neg[0](i), scalar.point_neg[0]) ||
        !close(packed.point_neg[1](i), scalar.point_neg[1]) ||
        !close(packed.point_neg[2](i), scalar.point_neg[2])) {
      std::cerr << "[FAIL] lane " << i << " mismatch\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "[PASS] accux_constlat Eigen-pack test passed.\n";


  const auto packed_try = try_gca_constlat_intersection(a, b, zz);

  for (int i = 0; i < N; ++i) {
    Vec3<double> ai{ptsA(i, 0), ptsA(i, 1), ptsA(i, 2)};
    Vec3<double> bi{ptsB(i, 0), ptsB(i, 1), ptsB(i, 2)};
    const auto scalar_try = try_gca_constlat_intersection(ai, bi, z0(i));

    if (!close(packed_try.point[0](i), scalar_try.point[0]) ||
        !close(packed_try.point[1](i), scalar_try.point[1]) ||
        !close(packed_try.point[2](i), scalar_try.point[2]) ||
        packed_try.status(i) != scalar_try.status) {
      std::cerr << "[FAIL] lane " << i
                << " try_gca_constlat_intersection mismatch\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "[PASS] gca_constlat Eigen-pack tests passed.\n";
  return EXIT_SUCCESS;
}
