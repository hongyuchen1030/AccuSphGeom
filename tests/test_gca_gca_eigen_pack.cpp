#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "accusphgeom/adapters/eigen/numeric.hpp"
#include "accusphgeom/constructions/gca_gca_intersection.hpp"

namespace {

using accusphgeom::constructions::accux_gca;
using accusphgeom::constructions::try_gca_gca_intersection;
using accusphgeom::numeric::Vec3;

constexpr int N = 4;
using Pack = EigenPack<N>;
constexpr double kTol = 0.0;

bool close(double a, double b) {
  return std::abs(a - b) <= kTol;
}

}  // namespace

int main() {
  Eigen::Matrix<double, N, 3> ptsA0;
  Eigen::Matrix<double, N, 3> ptsA1;
  Eigen::Matrix<double, N, 3> ptsB0;
  Eigen::Matrix<double, N, 3> ptsB1;

  ptsA0 <<
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.7071067811865475, 0.7071067811865475, 0.0,
      0.0, 0.0, 1.0;

  ptsA1 <<
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0,
      0.0, 0.0, 1.0,
      1.0, 0.0, 0.0;

  ptsB0 <<
      1.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      0.7071067811865475, 0.0, 0.7071067811865475;

  ptsB1 <<
      0.0, 0.0, 1.0,
      0.0, 1.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 1.0, 0.0;

  Pack a0x = ptsA0.template block<N, 1>(0, 0).array();
  Pack a0y = ptsA0.template block<N, 1>(0, 1).array();
  Pack a0z = ptsA0.template block<N, 1>(0, 2).array();

  Pack a1x = ptsA1.template block<N, 1>(0, 0).array();
  Pack a1y = ptsA1.template block<N, 1>(0, 1).array();
  Pack a1z = ptsA1.template block<N, 1>(0, 2).array();

  Pack b0x = ptsB0.template block<N, 1>(0, 0).array();
  Pack b0y = ptsB0.template block<N, 1>(0, 1).array();
  Pack b0z = ptsB0.template block<N, 1>(0, 2).array();

  Pack b1x = ptsB1.template block<N, 1>(0, 0).array();
  Pack b1y = ptsB1.template block<N, 1>(0, 1).array();
  Pack b1z = ptsB1.template block<N, 1>(0, 2).array();

  Vec3<Pack> a0{a0x, a0y, a0z};
  Vec3<Pack> a1{a1x, a1y, a1z};
  Vec3<Pack> b0{b0x, b0y, b0z};
  Vec3<Pack> b1{b1x, b1y, b1z};

  const auto packed = accux_gca(a0, a1, b0, b1);

  for (int i = 0; i < N; ++i) {
    Vec3<double> a0i{ptsA0(i, 0), ptsA0(i, 1), ptsA0(i, 2)};
    Vec3<double> a1i{ptsA1(i, 0), ptsA1(i, 1), ptsA1(i, 2)};
    Vec3<double> b0i{ptsB0(i, 0), ptsB0(i, 1), ptsB0(i, 2)};
    Vec3<double> b1i{ptsB1(i, 0), ptsB1(i, 1), ptsB1(i, 2)};
    const auto scalar = accux_gca(a0i, a1i, b0i, b1i);

    if (!close(packed.point_pos[0](i), scalar.point_pos[0]) ||
        !close(packed.point_pos[1](i), scalar.point_pos[1]) ||
        !close(packed.point_pos[2](i), scalar.point_pos[2]) ||
        !close(packed.point_neg[0](i), scalar.point_neg[0]) ||
        !close(packed.point_neg[1](i), scalar.point_neg[1]) ||
        !close(packed.point_neg[2](i), scalar.point_neg[2])) {
      std::cerr << "[FAIL] lane " << i << " accux_gca mismatch\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "[PASS] accux_gca Eigen-pack test passed.\n";

  const auto packed_try = try_gca_gca_intersection(a0, a1, b0, b1);

  for (int i = 0; i < N; ++i) {
    Vec3<double> a0i{ptsA0(i, 0), ptsA0(i, 1), ptsA0(i, 2)};
    Vec3<double> a1i{ptsA1(i, 0), ptsA1(i, 1), ptsA1(i, 2)};
    Vec3<double> b0i{ptsB0(i, 0), ptsB0(i, 1), ptsB0(i, 2)};
    Vec3<double> b1i{ptsB1(i, 0), ptsB1(i, 1), ptsB1(i, 2)};
    const auto scalar_try = try_gca_gca_intersection(a0i, a1i, b0i, b1i);

    if (!close(packed_try.point[0](i), scalar_try.point[0]) ||
        !close(packed_try.point[1](i), scalar_try.point[1]) ||
        !close(packed_try.point[2](i), scalar_try.point[2]) ||
        packed_try.status(i) != scalar_try.status) {
      std::cerr << "[FAIL] lane " << i
                << " try_gca_gca_intersection mismatch\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "[PASS] gca_gca Eigen-pack tests passed.\n";
  return EXIT_SUCCESS;
}
