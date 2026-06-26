#include <Eigen/Dense>

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <type_traits>

#include "accusphgeom/adapters/eigen/numeric.hpp"
#include "fp64_GCAconstLat.hh"

namespace {

using accusphgeom::numeric::Vec3;
using accusphgeom::performance_test::fp64_gca_constlat_pure_fp_xy;

constexpr int kMaxVecWidth = 8;
using Pack2 = EigenPack<2>;
using Pack4 = EigenPack<4>;
using Pack8 = EigenPack<8>;

constexpr std::size_t kDefaultDataSize = 100000;
constexpr std::size_t kDefaultNumTests = 100;

template <typename T>
constexpr int vec_width() {
  if constexpr (std::is_arithmetic_v<T>) {
    return 1;
  } else {
    static_assert(T::SizeAtCompileTime != Eigen::Dynamic,
                  "Pack type must have fixed compile-time size");
    return T::SizeAtCompileTime;
  }
}

static_assert(vec_width<double>() == 1);
static_assert(vec_width<Pack2>() == 2);
static_assert(vec_width<Pack4>() == 4);
static_assert(vec_width<Pack8>() == 8);

void create_output_directory(const std::string& directory) {
  struct stat info;
  if (stat(directory.c_str(), &info) != 0) {
    if (mkdir(directory.c_str(), 0777) == -1) {
      std::cerr << "Error creating directory: " << directory << '\n';
    }
  }
}

std::size_t parse_arg(const char* arg, const char* name) {
  const long long value = std::atoll(arg);
  if (value <= 0) {
    throw std::invalid_argument(std::string(name) + " must be positive");
  }
  return static_cast<std::size_t>(value);
}

template <typename T>
void do_not_optimize(const T& value) {
#if defined(__GNUC__) || defined(__clang__)
  asm volatile("" : : "g"(&value) : "memory");
#else
  (void)value;
#endif
}

template <typename T>
void compute_values_pure_fp(const Eigen::MatrixXd& ptsA,
                            const Eigen::MatrixXd& ptsB,
                            const Eigen::VectorXd& latitudes,
                            Eigen::VectorXd& x,
                            Eigen::VectorXd& y) {
  constexpr int stride = vec_width<T>();
  static_assert(stride > 0, "stride must be positive");

  const std::size_t size = static_cast<std::size_t>(ptsA.rows());
  assert(size % stride == 0);

  if ((ptsB.rows() != ptsA.rows()) || (latitudes.rows() != ptsA.rows()) ||
      (x.rows() != ptsA.rows()) || (y.rows() != ptsA.rows())) {
    throw std::runtime_error("All passed arrays must have the same size");
  }

  for (std::size_t i = 0; i < size; i += stride) {
    Vec3<T> a{};
    Vec3<T> b{};
    T z0{};

    if constexpr (std::is_arithmetic_v<T>) {
      for (int c = 0; c < 3; ++c) {
        a[c] = ptsA(static_cast<Eigen::Index>(i), c);
        b[c] = ptsB(static_cast<Eigen::Index>(i), c);
      }
      z0 = latitudes(static_cast<Eigen::Index>(i));
    } else {
      for (int c = 0; c < 3; ++c) {
        a[c] =
            ptsA.template block<stride, 1>(static_cast<Eigen::Index>(i), c)
                .array();
        b[c] =
            ptsB.template block<stride, 1>(static_cast<Eigen::Index>(i), c)
                .array();
      }
      z0 = latitudes.template segment<stride>(static_cast<Eigen::Index>(i))
               .array();
    }

    const auto [px, py] = fp64_gca_constlat_pure_fp_xy(a, b, z0);

    if constexpr (std::is_arithmetic_v<T>) {
      x(static_cast<Eigen::Index>(i)) = px;
      y(static_cast<Eigen::Index>(i)) = py;
    } else {
      x.template segment<stride>(static_cast<Eigen::Index>(i)).array() = px;
      y.template segment<stride>(static_cast<Eigen::Index>(i)).array() = py;
    }
  }
}

template <typename T>
double run_benchmark_pure_fp(std::size_t num_tests,
                             const Eigen::MatrixXd& ptsA,
                             const Eigen::MatrixXd& ptsB,
                             const Eigen::VectorXd& latitudes,
                             Eigen::VectorXd& x,
                             Eigen::VectorXd& y) {
  // One warm-up run outside timing.
  compute_values_pure_fp<T>(ptsA, ptsB, latitudes, x, y);

  const auto start = std::chrono::high_resolution_clock::now();

  for (std::size_t t = 0; t < num_tests; ++t) {
    compute_values_pure_fp<T>(ptsA, ptsB, latitudes, x, y);
    do_not_optimize(x.data());
    do_not_optimize(y.data());
  }

  const std::chrono::duration<double> elapsed =
      std::chrono::high_resolution_clock::now() - start;

  return elapsed.count();
}

double l2_error(const Eigen::VectorXd& lhs, const Eigen::VectorXd& rhs) {
  return (lhs - rhs).norm();
}

}  // namespace

int main(int argc, char** argv) {
#ifndef NDEBUG
  std::cerr << "Warning: benchmark is not compiled with NDEBUG. "
               "Use Release mode for meaningful timings.\n";
#endif

  const std::size_t data_size =
      (argc > 1) ? parse_arg(argv[1], "data_size") : kDefaultDataSize;
  const std::size_t num_tests =
      (argc > 2) ? parse_arg(argv[2], "num_tests") : kDefaultNumTests;

  if ((data_size % kMaxVecWidth) != 0) {
    throw std::invalid_argument("data_size must be divisible by 8");
  }

  std::srand(12345);

  Eigen::MatrixXd ptsA(data_size, 3);
  Eigen::MatrixXd ptsB(data_size, 3);
  Eigen::VectorXd latitudes(data_size);

  Eigen::VectorXd x_scalar(data_size);
  Eigen::VectorXd y_scalar(data_size);
  Eigen::VectorXd x_pack2(data_size);
  Eigen::VectorXd y_pack2(data_size);
  Eigen::VectorXd x_pack4(data_size);
  Eigen::VectorXd y_pack4(data_size);
  Eigen::VectorXd x_pack8(data_size);
  Eigen::VectorXd y_pack8(data_size);

  double pointA_x_significand = 5028374390644146;
  int pointA_x_exponent = -53;
  double pointA_y_significand = -7472957205960756;
  int pointA_y_exponent = -53;
  double pointA_z_significand = 6254432438282003;
  int pointA_z_exponent = -79;

  double pointB_x_significand = 5167685454902838;
  int pointB_x_exponent = -53;
  double pointB_y_significand = -7377307466399399;
  int pointB_y_exponent = -53;
  double pointB_z_significand = 4525606513452550;
  int pointB_z_exponent = -78;

  double constZ_significand = 7998403280412384;
  int constZ_exponent = -79;

  Eigen::Vector3d pointA(
      pointA_x_significand * std::pow(2.0, pointA_x_exponent),
      pointA_y_significand * std::pow(2.0, pointA_y_exponent),
      pointA_z_significand * std::pow(2.0, pointA_z_exponent));

  Eigen::Vector3d pointB(
      pointB_x_significand * std::pow(2.0, pointB_x_exponent),
      pointB_y_significand * std::pow(2.0, pointB_y_exponent),
      pointB_z_significand * std::pow(2.0, pointB_z_exponent));

  const double constZ = constZ_significand * std::pow(2.0, constZ_exponent);

  for (std::size_t i = 0; i < data_size; ++i) {
    ptsA.row(static_cast<Eigen::Index>(i)) = pointA;
    ptsB.row(static_cast<Eigen::Index>(i)) = pointB;
  }

  latitudes.setConstant(constZ);

  ptsA += 1e-12 * Eigen::MatrixXd::Random(ptsA.rows(), ptsA.cols());
  ptsB += 1e-12 * Eigen::MatrixXd::Random(ptsB.rows(), ptsB.cols());
  latitudes += 1e-12 * Eigen::VectorXd::Random(latitudes.size());

  x_scalar.setZero();
  y_scalar.setZero();
  const double scalar_time =
      run_benchmark_pure_fp<double>(num_tests, ptsA, ptsB, latitudes,
                                    x_scalar, y_scalar);

  x_pack2.setZero();
  y_pack2.setZero();
  const double pack2_time =
      run_benchmark_pure_fp<Pack2>(num_tests, ptsA, ptsB, latitudes,
                                   x_pack2, y_pack2);

  x_pack4.setZero();
  y_pack4.setZero();
  const double pack4_time =
      run_benchmark_pure_fp<Pack4>(num_tests, ptsA, ptsB, latitudes,
                                  x_pack4, y_pack4);

  x_pack8.setZero();
  y_pack8.setZero();
  const double pack8_time =
      run_benchmark_pure_fp<Pack8>(num_tests, ptsA, ptsB, latitudes,
                                   x_pack8, y_pack8);

  const double scalar_time_per_point =
      scalar_time / static_cast<double>(num_tests * data_size);
  const double pack2_time_per_point =
      pack2_time / static_cast<double>(num_tests * data_size);
  const double pack4_time_per_point =
      pack4_time / static_cast<double>(num_tests * data_size);
  const double pack8_time_per_point =
      pack8_time / static_cast<double>(num_tests * data_size);

  const double x_l2_pack2 = l2_error(x_pack2, x_scalar);
  const double y_l2_pack2 = l2_error(y_pack2, y_scalar);
  const double x_l2 = l2_error(x_pack4, x_scalar);
  const double y_l2 = l2_error(y_pack4, y_scalar);
  const double x_l2_pack8 = l2_error(x_pack8, x_scalar);
  const double y_l2_pack8 = l2_error(y_pack8, y_scalar);

  create_output_directory("tests/performance_test/output");
  std::ofstream csv("tests/performance_test/output/"
                    "gca_constlat_pure_fp_SIMDPack_timing.csv");
  csv << std::scientific << std::setprecision(16);
  csv << "method,threadsNum,vec_width,time\n";
  csv << "pure_fp,1,1," << scalar_time_per_point << "\n";
  csv << "pure_fp,1,2," << pack2_time_per_point << "\n";
  csv << "pure_fp,1,4," << pack4_time_per_point << "\n";
  csv << "pure_fp,1,8," << pack8_time_per_point << "\n";

  std::cout << std::fixed << std::setprecision(3);
  std::cout << "gca_constlat pure-FP benchmark complete\n";
  std::cout << "  data_size=" << data_size << ", num_tests=" << num_tests
            << "\n";
  std::cout << "  pure_fp width=1 : " << scalar_time_per_point
            << " s/point\n";
  std::cout << "  pure_fp width=2 : " << pack2_time_per_point
            << " s/point\n";
  std::cout << "  pure_fp width=4 : " << pack4_time_per_point
            << " s/point\n";
  std::cout << "  pure_fp width=8 : " << pack8_time_per_point
            << " s/point\n";
  std::cout << "  width=2 vs width=1 x_l2=" << std::scientific
            << x_l2_pack2 << ", y_l2=" << y_l2_pack2 << "\n";
  std::cout << "  width=4 vs width=1 x_l2=" << x_l2
            << ", y_l2=" << y_l2 << "\n";
  std::cout << "  width=8 vs width=1 x_l2=" << x_l2_pack8
            << ", y_l2=" << y_l2_pack8 << "\n";
  std::cout << "  timing CSV: tests/performance_test/output/"
               "gca_constlat_pure_fp_SIMDPack_timing.csv\n";

  return EXIT_SUCCESS;
}
